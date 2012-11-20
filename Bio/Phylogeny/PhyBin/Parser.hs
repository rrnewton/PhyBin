{-# OPTIONS_GHC -fwarn-incomplete-patterns #-}
{-# OPTIONS_GHC -fwarn-unused-imports #-}

module Bio.Phylogeny.PhyBin.Parser
       (newick_parser, parseNewick)
       where
import qualified Data.ByteString.Lazy.Char8 as B
import           Data.Char          (isSpace)
import           Text.Parsec
import           Text.Parsec.ByteString.Lazy
import           Test.HUnit         ()
import           Bio.Phylogeny.PhyBin.CoreTypes (NewickTree(..), DefDecor, toLabel)


-- | Parse a bytestring into a NewickTree with branch lengths.  The
--   first argument is file from which the data came and is just for
--   error error messages.
parseNewick :: String -> B.ByteString -> NewickTree DefDecor
parseNewick file input = 
  runB file newick_parser $
  B.filter (not . isSpace) input

runB :: Show a => String -> Parser a -> B.ByteString -> a
runB file p input = case (parse p "" input) of
	         Left err -> error ("parse error in file "++ show file ++" at "++ show err)
		 Right x  -> x

----------------------------------------------------------------------------------------------------
-- Newick file format parser definitions:
----------------------------------------------------------------------------------------------------

tag :: a -> NewickTree a -> NewickTree a
tag l s =
  case s of 
    NTLeaf _ n      -> NTLeaf l n
    NTInterior _ ls -> NTInterior l ls

-- | This parser ASSUMES that whitespace has been prefiltered from the input.
newick_parser :: Parser (NewickTree DefDecor)
newick_parser = 
   do x <- subtree
      -- Get the top-level metadata:
      l <- branchMetadat
      _ <- char ';'
      return$ tag l x

subtree :: Parser (NewickTree DefDecor)
subtree = internal <|> leaf

defaultMeta :: (Maybe Int, Double)
defaultMeta = (Nothing,0.0)

leaf :: Parser (NewickTree DefDecor)
leaf = do n<-name; return$ NTLeaf defaultMeta (toLabel n)

internal :: Parser (NewickTree DefDecor)
internal = do _   <- char '('       
	      bs  <- branchset
	      _   <- char ')'       
              _nm <- name -- IGNORED
              return$ NTInterior defaultMeta bs

branchset :: Parser [NewickTree DefDecor]
branchset =
    do b <- branch <?> "at least one branch"
       rest <- option [] $ try$ do char ','; branchset
       return (b:rest)

branch :: Parser (NewickTree DefDecor)
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

name :: Parser String
name = option "" $ many1 (letter <|> digit <|> oneOf "_.-")



--------------------------------------------------------------------------------
-- Unit Tests
--------------------------------------------------------------------------------

-- tests = test 
--    [ "test name"   ~: "foo" ~=?  run name "foo"
--    , "test number" ~:  3.3  ~=?  run number "3.3"
--    , "test number" ~:  3.0  ~=?  run number "3"
--    , "test number" ~:  -3.0 ~=?  run number "-3"

--    , "leaf"     ~: ntl "A" ~=?  run leaf    "A"
--    , "subtree"  ~: ntl "A" ~=?  run subtree "A"

--    -- These are not allowed:
--    , "null branchset" ~: errortest$ run branchset ""

--    , "internal" ~: NTInterior (Nothing,0.0) [ntl "A"] ~=?  run internal "(A);"
--    ]
--  where
--    ntl s = NTLeaf (Nothing,0.0) (toLabel s)


   --   "example: no nodes are named"  ~: NTInterior (Nothing,0)
   --                                       [ntl "", ntl "", NTInterior (Nothing,0) [ntl "", ntl ""]]
   --      			   ~=? run newick_parser "(,,(,));"
   -- , "example: leaf nodes are named" ~: 6 ~=?  cnt (run newick_parser "(A,B,(C,D));")
   -- , "example: all nodes are named"  ~: 6 ~=?  cnt (run newick_parser "(A,B,(C,D)E)F;")

   -- , "example: all but root node have a distance to parent"  ~: 6 ~=? cnt (run newick_parser "(:0.1,:0.2,(:0.3,:0.4):0.5);")
   -- , "example: all have a distance to parent"              ~: 6 ~=? cnt (run newick_parser "(:0.1,:0.2,(:0.3,:0.4):0.5):0.6;")
   -- , "example: distances and leaf names (popular)"         ~: 6 ~=? cnt tre1
   -- , "example: distances and all names"                    ~: 6 ~=? cnt (run newick_parser "(A:0.1,B:0.2,(C:0.3,D:0.4)E:0.5)F;")
   -- , "example: a tree rooted on a leaf node (rare)"        ~: 6 ~=? cnt (run newick_parser "((B:0.2,(C:0.3,D:0.4)E:0.5)F:0.1)A;")


run :: Show a => Parser a -> String -> a
run p input = runB "<unknown>" p (B.pack input)
              