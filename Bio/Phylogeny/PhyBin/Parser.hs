{-# LANGUAGE ScopedTypeVariables #-}
{-# OPTIONS_GHC -fwarn-incomplete-patterns #-}
{-# OPTIONS_GHC -fwarn-unused-imports #-}

module Bio.Phylogeny.PhyBin.Parser
       (newick_parser, parseNewick, unitTests)
       where
import           Control.Exception  (evaluate, handle, SomeException)
import qualified Data.ByteString.Lazy.Char8 as B
import           Data.Char          (isSpace)
import           Data.Map as M
import           Text.Parsec
import           Text.Parsec.ByteString.Lazy
import           Test.HUnit          ((~:),(~=?),Test,test,assertFailure)
import           Bio.Phylogeny.PhyBin.CoreTypes (NewickTree(..), DefDecor, treeSize, LabelTable)
import           Prelude as P


-- | Parse a bytestring into a NewickTree with branch lengths.  The
--   first argument is file from which the data came and is just for
--   better error messages.
parseNewick :: String -> B.ByteString -> (LabelTable, NewickTree DefDecor)
parseNewick file input =
  extractLabelTable $ 
  runB file newick_parser $
  B.filter (not . isSpace) input

runB :: Show a => String -> Parser a -> B.ByteString -> a
runB file p input = case (parse p "" input) of
	         Left err -> error ("parse error in file "++ show file ++" at "++ show err)
		 Right x  -> x


extractLabelTable :: TempTree -> (LabelTable, NewickTree DefDecor)
extractLabelTable tr = loop M.empty tr
 where
   loop acc (NTLeaf (d,Just nm) _) =
     let nxt = M.size acc in
     (M.insert nxt nm acc,  NTLeaf d nxt)
   loop acc1 (NTInterior (d,Nothing) chlds) =
     let (acc',ls') = 
          P.foldr (\ x (acc2,ls) ->
                   let (acc3,x') = loop acc2 x in
                   (acc3, x':ls))
                 (acc1,[]) chlds
     in (acc', NTInterior d ls')
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
newick_parser :: Parser TempTree
newick_parser = 
   do x <- subtree
      -- Get the top-level metadata:
      l <- branchMetadat
      _ <- char ';'
      return$ tag l x

subtree :: Parser TempTree
subtree = internal <|> leaf

defaultMeta :: (Maybe Int, Double)
defaultMeta = (Nothing,0.0)

leaf :: Parser TempTree
leaf = do nm<-name; return$ NTLeaf (defaultMeta,Just nm) 0 

internal :: Parser TempTree
internal = do _   <- char '('       
	      bs  <- branchset
	      _   <- char ')'       
              _nm <- name -- IGNORED, internal names.
              return$ NTInterior (defaultMeta,Nothing) bs

branchset :: Parser [TempTree]
branchset =
    do b <- branch <?> "at least one branch"
       rest <- option [] $ try$ do char ','; branchset
       return (b:rest)

branch :: Parser TempTree
branch = do s<-subtree; l<-branchMetadat; 
	    return$ tag l s

-- If the length is omitted, it is implicitly zero.
branchMetadat :: Parser DefDecor 
branchMetadat = option defaultMeta $ do
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
tre1 = snd $ parseNewick "" $ B.pack "(A:0.1,B:0.2,(C:0.3,D:0.4):0.5);"

unitTests :: Test
unitTests = test
   [ "test name"   ~: "foo" ~=?  run name "foo"
   , "test number" ~:  3.3  ~=?  run number "3.3"
   , "test number" ~:  3.0  ~=?  run number "3"
   , "test number" ~:  -3.0 ~=?  run number "-3"

   , "leaf"     ~: ntl "A" ~=?  run leaf    "A"
   , "subtree"  ~: ntl "A" ~=?  run subtree "A"

     -- These are not allowed:
   , "null branchset" ~: errortest$ run branchset ""
     
   , "internal" ~: NTInterior nada [ntl "A"] ~=?  run internal "(A);"
   , "example: no nodes are named"  ~: NTInterior nada
                                         [ntl "", ntl "", NTInterior nada [ntl "", ntl ""]]
        			   ~=? run newick_parser "(,,(,));"
   , "example: leaf nodes are named" ~: 6 ~=?  treeSize (run newick_parser "(A,B,(C,D));")
   , "example: all nodes are named"  ~: 6 ~=?  treeSize (run newick_parser "(A,B,(C,D)E)F;")

   , "example: all but root node have a distance to parent"  ~: 6 ~=? treeSize (run newick_parser "(:0.1,:0.2,(:0.3,:0.4):0.5);")
   , "example: all have a distance to parent"              ~: 6 ~=? treeSize (run newick_parser "(:0.1,:0.2,(:0.3,:0.4):0.5):0.6;")
   , "example: distances and leaf names (popular)"         ~: 6 ~=? treeSize tre1
   , "example: distances and all names"                    ~: 6 ~=? treeSize (run newick_parser "(A:0.1,B:0.2,(C:0.3,D:0.4)E:0.5)F;")
   , "example: a tree rooted on a leaf node (rare)"        ~: 6 ~=? treeSize (run newick_parser "((B:0.2,(C:0.3,D:0.4)E:0.5)F:0.1)A;")
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
              
