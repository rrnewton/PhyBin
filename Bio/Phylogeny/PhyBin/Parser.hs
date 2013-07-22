{-# LANGUAGE ScopedTypeVariables, BangPatterns #-}
{-# OPTIONS_GHC -fwarn-incomplete-patterns #-}
{-# OPTIONS_GHC -fwarn-unused-imports #-}

module Bio.Phylogeny.PhyBin.Parser
       (newick_parser, parseNewick, parseNewicks, unitTests)
       where
import           Control.Exception  (evaluate, handle, SomeException)
import qualified Data.ByteString.Lazy.Char8 as B
import           Data.Char          (isSpace)
import           Data.Map as M
import           Data.Set as S
import           Data.List as L
import           Text.Parsec
import           Text.Parsec.ByteString.Lazy
import           Test.HUnit          ((~:),(~=?),Test,test,assertFailure)
import           Bio.Phylogeny.PhyBin.CoreTypes (NewickTree(..), DefDecor, treeSize, LabelTable)
import           Prelude as P

type NameHack = (String->String)

-- | Parse a bytestring into a NewickTree with branch lengths.  The
--   first argument is file from which the data came and is just for
--   better error messages.
parseNewick :: LabelTable -> NameHack -> String -> B.ByteString -> (LabelTable, NewickTree DefDecor)
parseNewick tbl0 name_hack file input =
  extractLabelTable tbl0 $ 
  runB file (newick_parser name_hack) $
  B.filter (not . isSpace) input

-- | Parse a list of trees, starting with an empty map of labels and accumulating a final map.
parseNewicks :: NameHack -> [(String,B.ByteString)] -> (LabelTable, [NewickTree DefDecor])
parseNewicks name_hack pairs =
  P.foldr fn (M.empty,[]) pairs
 where
   fn (file,bstr) (!acc,!ls) =
     let (acc',tr) = parseNewick acc name_hack file bstr
     in (acc', tr:ls)

runB :: Show a => String -> Parser a -> B.ByteString -> a
runB file p input = case (parse p "" input) of
	         Left err -> error ("parse error in file "++ show file ++" at "++ show err)
		 Right x  -> x

extractLabelTable :: LabelTable -> TempTree -> (LabelTable, NewickTree DefDecor)
extractLabelTable tbl0 tr = (finMap,finTree)
 where
   flipped = M.fromList $ L.map (\(x,y)->(y,x)) $ M.toList tbl0
   -- (_,finMap,finTree) = loop (S.fromList (M.elems tbl0)) tbl0 tr
   (_,finMap,finTree) = loop flipped tbl0 tr
   
   loop seen acc (NTLeaf (d,Just nm) _)
     | M.member nm seen = (seen, acc, NTLeaf d (seen M.! nm))
     | otherwise = let nxt = M.size acc in
                   (M.insert nm nxt seen,
                    M.insert nxt nm acc,  NTLeaf d nxt)
   loop seen1 acc1 (NTInterior (d,Nothing) chlds) =
     let (seen',acc',ls') = 
          P.foldr (\ x (seen2,acc2,ls) ->
                   let (seen3,acc3,x') = loop seen2 acc2 x in
                   (seen3, acc3, x':ls))
                  (seen1,acc1,[])
                  chlds
     in (seen',acc', NTInterior d ls')

----------------------------------------------------------------------------------------------------
-- Newick file format parser definitions:
----------------------------------------------------------------------------------------------------

-- | Hack: we store the names in the leaves.
type TempTree = NewickTree (DefDecor,Maybe String)

tag :: DefDecor -> TempTree -> TempTree
tag l s =
  case s of 
    NTLeaf (_,m) n            -> NTLeaf (l,m) n
    NTInterior (_,Nothing) ls -> NTInterior (l,Nothing) ls

-- | This parser ASSUMES that whitespace has been prefiltered from the input.
newick_parser :: NameHack -> Parser TempTree
newick_parser name_hack = 
   do x <- subtree name_hack
      -- Get the top-level metadata:
      l <- branchMetadat name_hack
      _ <- char ';'
      return$ tag l x

subtree :: NameHack -> Parser TempTree
subtree name_hack = internal name_hack <|> leaf name_hack

defaultMeta :: (Maybe Int, Double)
defaultMeta = (Nothing,0.0)

leaf :: NameHack -> Parser TempTree
leaf name_hack =
  do nm <- name
     let nm' = name_hack nm
     return$ NTLeaf (defaultMeta,Just nm') 0 

internal :: NameHack -> Parser TempTree
internal name_hack =
   do _   <- char '('       
      bs  <- branchset name_hack
      _   <- char ')'       
      _nm <- name -- IGNORED, internal names.
      return$ NTInterior (defaultMeta,Nothing) bs

branchset :: NameHack -> Parser [TempTree]
branchset name_hack =
    do b <- branch name_hack <?> "at least one branch"
       rest <- option [] $ try$ do char ','; branchset name_hack
       return (b:rest)

branch :: NameHack -> Parser TempTree
branch name_hack =
         do s <- subtree name_hack
            l <- branchMetadat name_hack 
            return$ tag l s

-- If the length is omitted, it is implicitly zero.
branchMetadat :: NameHack -> Parser DefDecor 
branchMetadat name_hack = option defaultMeta $ do
    char ':'
    n <- (try sciNotation <|> number)
    -- IF the branch length succeeds then we go for the bracketed bootstrap value also:
    bootstrap <- option Nothing $ do
      char '['
      s <- many1 digit
      char ']'
      return (Just (read s))
    return (bootstrap,n)

-- | Parse a normal, decimal number.
number :: Parser Double
number = 
  do sign   <- option "" $ string "-"
     first  <- many1 digit
     second <- option "0" $ try$ do char '.'; many1 digit
     return (read (sign ++ first++"."++second) :: Double)

-- | Parse a number in scientific notation.
sciNotation :: Parser Double
sciNotation =
  do coeff <- do first <- many1 digit
                 second <- option "0" $ try$ do char '.'; many1 digit
                 return $ first++"."++second
     char 'e'
     sign  <- option "" $ string "-"
     expon <- many1 digit
     return (read (coeff++"e"++sign++expon))

-- | Names are a mess... they contain all kinds of garbage sometimes it seems.
--   Thus we are very permissive.  We allow anything that is not something we specifically need to reserve.
name :: Parser String
name = option "" $ many1 (noneOf "()[]:;, \t\n\r\f\v")
-- name = option "" $ many1 (letter <|> digit <|> oneOf "_.-")


--------------------------------------------------------------------------------
-- Unit Tests
--------------------------------------------------------------------------------

tre1 :: NewickTree DefDecor
tre1 = snd $ parseNewick M.empty id "" $ B.pack "(A:0.1,B:0.2,(C:0.3,D:0.4):0.5);"

unitTests :: Test
unitTests = test
   [ "test name"   ~: "foo" ~=?  run name "foo"
   , "test number" ~:  3.3  ~=?  run number "3.3"
   , "test number" ~:  3.0  ~=?  run number "3"
   , "test number" ~:  -3.0 ~=?  run number "-3"

   , "leaf"     ~: ntl "A" ~=?  run (leaf id)    "A"
   , "subtree"  ~: ntl "A" ~=?  run (subtree id) "A"

     -- These are not allowed:
   , "null branchset" ~: errortest$ run (branchset id) ""
     
   , "internal" ~: NTInterior nada [ntl "A"] ~=?  run (internal id) "(A);"
   , "example: no nodes are named"  ~: NTInterior nada
                                         [ntl "", ntl "", NTInterior nada [ntl "", ntl ""]]
        			   ~=? run (newick_parser id) "(,,(,));"
   , "example: leaf nodes are named" ~: 6 ~=?  treeSize (run (newick_parser id) "(A,B,(C,D));")
   , "example: all nodes are named"  ~: 6 ~=?  treeSize (run (newick_parser id) "(A,B,(C,D)E)F;")

   , "example: all but root node have a distance to parent"  ~: 6 ~=? treeSize (run (newick_parser id) "(:0.1,:0.2,(:0.3,:0.4):0.5);")
   , "example: all have a distance to parent"              ~: 6 ~=? treeSize (run (newick_parser id) "(:0.1,:0.2,(:0.3,:0.4):0.5):0.6;")
   , "example: distances and leaf names (popular)"         ~: 6 ~=? treeSize tre1
   , "example: distances and all names"                    ~: 6 ~=? treeSize (run (newick_parser id) "(A:0.1,B:0.2,(C:0.3,D:0.4)E:0.5)F;")
   , "example: a tree rooted on a leaf node (rare)"        ~: 6 ~=? treeSize (run (newick_parser id) "((B:0.2,(C:0.3,D:0.4)E:0.5)F:0.1)A;")
   ]
 where
  ntl s = NTLeaf ((Nothing,0.0),Just s) 0
  nada = ((Nothing,0),Nothing)

run :: Show a => Parser a -> String -> a
run p input = runB "<unknown>" p (B.pack input)

errortest :: t -> IO ()
errortest x = 
   --() ~=?
    handle (\ (e::SomeException) -> return ()) $ 
      do evaluate x
         assertFailure "test was expected to throw an error"
              
