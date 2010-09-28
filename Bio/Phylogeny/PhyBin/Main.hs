{-# LANGUAGE ScopedTypeVariables, RecordWildCards, TypeSynonymInstances, CPP #-}
{-# OPTIONS_GHC -fwarn-incomplete-patterns #-}

--module Bio.Phylogeny.PhyBin.Main where
module Main where

import Text.Printf
import Text.Parsec
import Text.Parsec.Char
import Text.Parsec.Combinator
import Text.Parsec.ByteString.Lazy 

import Numeric
import Data.Function
import Data.List
import Data.Maybe
import Data.Char
import qualified Data.ByteString.Lazy.Char8 as B
import qualified Data.Map as M
import qualified Data.Set as S

import Control.Monad
import Control.Exception hiding (try)
--import Control.DeepSeq
import Control.Concurrent

-- Wow, I actually couldn't figure out how to open a file (and get a
-- HANDLE) so that I could then use getFileAttributes under
-- System.Win32.  Giving up because I think I can just use the
-- OS-independent System.Directory

import System.Directory
import System.FilePath

import System.Environment
import System.Directory
import System.Console.GetOpt
import System.Exit
import System.IO

import Test.HUnit
import HSH hiding (run)
import qualified HSH 

-- For vizualization:
import Data.Graph.Inductive as G  hiding (run) 
import Data.GraphViz        as Gv hiding (Label) 
import qualified Data.Graph.Inductive as G 
import qualified Data.GraphViz        as Gv
--import Data.Graph.Inductive.Query.DFS

import Text.PrettyPrint.HughesPJClass hiding (char, Style)

import Debug.Trace

phybin_version = "0.1.2" -- NOTE: Remember to keep me in sync with the .cabal file

----------------------------------------------------------------------------------------------------
-- Type definitions
----------------------------------------------------------------------------------------------------

type BranchLen = Double

-- Even though the Newick format allows it, ignoring interior node
-- labels. (They are not commonly used.)
data NewickTree a = 
   NTLeaf     a Label
 | NTInterior a [NewickTree a]
 deriving (Show, Eq, Ord)

{-
-- [2010.09.22] Disabling:
instance NFData Atom where
  rnf a = rnf (fromAtom a :: Int)

instance NFData a => NFData (NewickTree a) where
  rnf (NTLeaf l n)      = rnf (l,n)
  rnf (NTInterior l ls) = rnf (l,ls)
-}

instance Pretty (NewickTree dec) where 
 pPrint (NTLeaf dec name)   = text (fromLabel name)
 pPrint (NTInterior dec ls) = 
     --parens$ commasep ls
     (parens$ sep$ map_but_last (<>text",") $ map pPrint ls)


-- Experimental: toggle this to change the representation of labels:
----------------------------------------
--type Label = Atom; (toLabel, fromLabel) = (toAtom, fromAtom)
----------------------------------------
type Label = String; (toLabel, fromLabel) = (id, id)
----------------------------------------
fromLabel :: Label -> String

----------------------------------------------------------------------------------------------------
-- OS specific bits:
----------------------------------------------------------------------------------------------------
-- #ifdef WIN32
-- is_regular_file = undefined
-- is_directory path = 
--   getFileAttributes
-- --getFileInformationByHandle
-- --    bhfiFileAttributes
-- file_exists = undefined
-- #else
-- is_regular_file :: FilePath -> IO Bool
-- is_regular_file file = 
--   do stat <- getFileStatus file; 
--      -- Hmm, this is probably bad practice... hard to know its exhaustive:
--      return$ isRegularFile stat || isNamedPipe stat || isSymbolicLink stat
-- is_directory :: FilePath -> IO Bool
-- is_directory path = 
--   do stat <- getFileStatus path
--      return (isDirectory stat)
-- file_exists = fileExist
-- #endif

-- Here we ASSUME it exists, then these functions are good enough:
is_regular_file = doesFileExist
is_directory = doesDirectoryExist 
file_exists path = 
  do f <- doesFileExist path
     d <- doesDirectoryExist path
     return (f || d)

----------------------------------------------------------------------------------------------------
-- General helper/utility functions:
----------------------------------------------------------------------------------------------------

--commacat ls = hcat (intersperse (text ", ") $ map pPrint ls)
commasep ls = sep (intersperse (text ", ") $ map pPrint ls)

map_but_last fn [] = []
map_but_last fn [h] = [h]
map_but_last fn (h:t) = fn h : map_but_last fn t

fst3 (a,_,_) = a
snd3 (_,b,_) = b
thd3 (_,_,c) = c

merge [] ls = ls
merge ls [] = ls
merge l@(a:b) r@(x:y) = 
  if a < x
  then a : merge b r
  else x : merge y l 

-- Set subtraction for sorted lists:
demerge ls [] = ls
demerge [] ls = error$ "demerge: first list did not contain all of second, remaining: " ++ show ls
demerge l@(a:b) r@(x:y) = 
  case a `compare` x of
   EQ -> demerge b y
   LT -> a : demerge b r 
   GT -> error$ "demerge: element was missing from first list: "++ show x

maybeCons Nothing  ls = ls
maybeCons (Just x) ls = x : ls

maybeInsert fn Nothing  ls = ls
maybeInsert fn (Just x) ls = insertBy fn x ls
--maybeInsert  Nothing  ls = ls
--maybeInsert  (Just x) ls = insert (x) ls -- Relies on ORD

----------------------------------------------------------------------------------------------------
-- Newick file format parser definitions:
----------------------------------------------------------------------------------------------------

tag l s =
  case s of 
    NTLeaf _ n      -> NTLeaf l n
    NTInterior _ ls -> NTInterior l ls

-- This parser ASSUMES that whitespace has been prefiltered from the input.
newick_parser :: Parser (NewickTree BranchLen)
newick_parser = 
   do x <- subtree
      l<-len
      char ';'
      return$ tag l x

subtree :: Parser (NewickTree BranchLen)
subtree = internal <|> leaf

leaf :: Parser (NewickTree BranchLen)
leaf = do n<-name; return$ NTLeaf 0.0 (toLabel n)

internal :: Parser (NewickTree BranchLen)
internal = do char '('       
	      bs <- branchset
	      char ')'       
              nm <- name -- IGNORED
              return$ NTInterior 0.0 bs

branchset :: Parser [NewickTree BranchLen]
branchset =
    do b <- branch <?> "at least one branch"
       rest <- option [] $ try$ do char ','; branchset
       return (b:rest)

branch :: Parser (NewickTree BranchLen)
branch = do s<-subtree; l<-len; 
	    return$ tag l s

len :: Parser Double
len = option 0.0 $ do char ':'; number

number :: Parser Double
number = 
  do sign <- option "" $ string "-"
     fst <- many1 digit
     snd <- option "0" $ try$ do char '.'; many1 digit
     return (read (sign ++ fst++"."++snd) :: Double)

name :: Parser String
name = option "" $ many1 (letter <|> digit <|> oneOf "_.-")


----------------------------------------------------------------------------------------------------
-- Normal form for unordered, unrooted trees
----------------------------------------------------------------------------------------------------

-- The basic idea is that what we *want* is the following, 
--   ROOT: most balanced point
--   ORDER: sorted in increasing subtree weight

-- But that's not quite good enough.  There are ties to break.  To do
-- that we fall back on the (totally ordered) leaf labels.

--------------------------------------------------------------------------------

-- A common type of tree is "AnnotatedTree", which contains three things:
--  (1) branch length from parent to "this" node
--  (2) subtree weights for future use
--      (defined as number of LEAVES, not counting intermediate nodes)
--  (3) sorted lists of labels for symmetry breaking
type AnnotatedTree = NewickTree (BranchLen, Int, [Label])


annotateWLabLists :: NewickTree BranchLen -> AnnotatedTree
annotateWLabLists tr = case tr of 
  NTLeaf bl n      -> NTLeaf (bl,1,[n]) n
  NTInterior bl ls -> 
      let children = map annotateWLabLists ls in 
      NTInterior (bl,      sum $ map get_weight children,
		  foldl1 merge $ map get_label_list   children)
		 children

instance Functor NewickTree where 
   fmap fn (NTLeaf dec x)      = NTLeaf (fn dec) x 
   fmap fn (NTInterior dec ls) = NTInterior (fn dec) (map (fmap fn) ls)

map_labels fn (NTLeaf     dec lab) = NTLeaf dec $ fn lab
map_labels fn (NTInterior dec ls)  = NTInterior dec$ map (map_labels fn) ls

all_labels (NTLeaf     _ lab) = [lab]
all_labels (NTInterior _ ls)  = concat$ map all_labels ls

get_dec (NTLeaf     dec _) = dec
get_dec (NTInterior dec _) = dec

-- Set all the decorations to a constant:
set_dec d = fmap (const d)
--set_dec d (NTLeaf _ x) = NTLeaf d x
--set_dec d (NTInterior _ ls) = NTInterior d $ map (set_dec d) ls

get_children (NTLeaf _ _) = []
get_children (NTInterior _ ls) = ls

-- Number of LEAVES contained in subtree:
get_weight = snd3 . get_dec

-- Sorted list of leaf labels contained in subtree
get_label_list   = thd3 . get_dec

add_weight (l1,w1,sorted1) node  = 
  let (_,w2,sorted2) = get_dec node in 
  (l1, w1+w2, merge sorted1 sorted2)

-- Remove the influence of one subtree from the metadata of another.
subtract_weight (l1,w1,sorted1) node =  
  let (_,w2,sorted2) = get_dec node in 
  (l1, w1-w2, demerge sorted1 sorted2)

-- Turn on for extra invariant checking:
debug = False
	
--------------------------------------------------------------------------------

-- I ran into a nasty bug as a result of "deriving Ord".  But I didn't end up getting rid of it.
--instance Ord AnnotatedTree where 
--  compare (NTLeaf _ _) (NTInterior _ _) = LT
--  compare (NTLeaf _ _) (NTLeaf _ _)     = EQ
--compare_nodes :: AnnotatedTree -> AnnotatedTree -> Ordering
-- Our sorting criteria for the children of interior nodes:

compare_childtrees node1 node2 = 
    case get_weight node1 `compare` get_weight node2 of 
     -- Comparisons on atoms cause problems WRT to determinism between runs if parallelism is introduced.
     -- Can consider it an optimization for the serial case perhaps:
--     EQ -> case map deAtom (get_label_list node1) `compare` 
--	        map deAtom (get_label_list node2) of
     EQ -> case get_label_list node1 `compare` get_label_list node2 of
            --EQ -> error "FIXME EQ"
            EQ -> error$ "Internal invariant broken.  These two children have equal ordering priority:\n" 
		  ++ "Pretty printing:\n  "
		  ++ show (pPrint node1) ++ "\n  " ++ show (pPrint node2)
		  ++ "\nFull data structures:\n  "
		  ++ show (node1) ++ "\n  " ++ show (node2)

	    x  -> x
     x -> x


-- This is it, here's the routine that transforms a tree into normal form.
-- This relies HEAVILY on lazy evaluation.
normalize :: AnnotatedTree -> AnnotatedTree
normalize tree = snd$ loop tree Nothing
 where 

  add_context dec Nothing  = dec
  add_context dec (Just c) = add_weight dec c

  -- loop: Walk over the tree.
  -- Inputs: 
  --    1. node: the NewickTree node to process ("us")
  --    2. context: all nodes connected through the parent, "flipped" as though *we* were root
  --                The "flipped" part has ALREADY been normalized.
  -- Outputs: 
  --    1. new node
  --    3. the best candidate root anywhere under this subtree
  loop :: AnnotatedTree -> Maybe AnnotatedTree -> (AnnotatedTree, AnnotatedTree)
  loop node context  = case node of
    NTLeaf dec@(l,w,sorted) name -> 
	(node, 
	 -- If the leaf becomes the root... we could introduce another node:
	 NTInterior (add_context (0,w,sorted) context) $
	            (verify_sorted "1" id$ maybeInsert compare_childtrees context [node])

	 -- It may be reasonable to not support leaves becoming root.. that changes the number of nodes!
	            --error "normalize: leaf becoming root not currently supported."
	)
    
    NTInterior dec@(l,w,_) ls -> 
     let 
         -- If this node becomes the root, the parent becomes one of our children:
         inverted = NTInterior inverted_dec inverted_children
	 inverted_dec      = add_context dec context
         inverted_children = verify_sorted "2" id$ maybeInsert compare_childtrees context newchildren

	 newchildren = --trace ("SORTED "++ show (map (get_label_list . fst) sorted)) $
		       map fst sorted
         sorted = sortBy (compare_childtrees `on` fst) possibs

         possibs = 
	  flip map ls $ \ child -> 
	   let 

	       -- Will this diverge???  Probably depends on how equality (for delete) is defined... 

	       -- Reconstruct the current node missing one child (because it became a parent):
	       -- Update its metadata appropriately:
	       newinverted = NTInterior (subtract_weight inverted_dec child) 
			                (verify_sorted "3" id$ delete newnode inverted_children)
	       (newnode, _) = result

  	       result = loop child (Just newinverted) 
	   in
	       result
	 
         -- Either us or a candidate suggested by one of the children:
         rootcandidates = inverted : map snd sorted

         -- Who wins?  The "most balanced".  Minimize max subtree weight.
	 -- The compare operator is NOT allowed to return EQ here.  Therefore there will be a unique minima.
	 winner = --trace ("Candidates: \n"++ show (nest 6$ vcat (map pPrint (zip (map max_subtree_weight rootcandidates) rootcandidates )))) $ 
		  minimumBy cmpr_subtree_weight rootcandidates

	 max_subtree_weight = maximum . map get_weight . get_children 
	 fat_id = map get_label_list . get_children 

         cmpr_subtree_weight tr1 tr2 = 
           case max_subtree_weight  tr1 `compare` max_subtree_weight tr2 of
	     EQ -> -- As a fallback we compare the alphabetic order of the "bignames" of the children:
                   case fat_id tr1 `compare` fat_id tr2 of 
		     EQ -> error$ "\nInternal invariant broken.  These two were equally good roots:\n" 
			          ++ show (pPrint tr1) ++ "\n" ++ show (pPrint tr2)
		     x -> x
	     x -> x

     in (NTInterior dec newchildren, winner)


-- Verify that our invariants are met:
verify_sorted msg = 
 if debug 
 then \ project nodes ->
  let weights = map (get_weight . project) nodes in 
    if sort weights == weights
    then nodes
--    else error$ "Child list failed verification: "++ show (pPrint nodes)
    else error$ msg ++ ": Child list failed verification, not sorted: "++ show (weights)
	        ++"\n  "++ show (sep $ map pPrint nodes) ++ 
                "\n\nFull output:\n  " ++ (concat$ intersperse "\n  " $ map show nodes)
 else \ _ nodes -> nodes


-- TODO: Salvage any of these tests that are worthwhile and get them into the unit tests:	        	
tt = normalize $ annotateWLabLists $ run newick_parser "(A,(C,D,E),B);"

tt0 = drawNewickTree "tt0" $ annotateWLabLists $ run newick_parser "(A,(C,D,E),B);"

tt2 = toGraph tt
tt3 = drawNewickTree "tt3" tt

norm4 = norm "((C,D,E),B,A);"
tt4 = drawNewickTree "tt4"$ trace ("FINAL: "++ show (pPrint norm4)) $ norm4

norm5 = normalize$ annotateWLabLists$ run newick_parser "(D,E,C,(B,A));"
tt5 = drawNewickTree "tt5"$ norm5

tt5' = prettyPrint' $ dotNewickTree "norm5" 1.0 norm5

ttall = do tt3; tt4; tt5

----------------------------------------------------------------------------------------------------
-- Equivalence classes:
----------------------------------------------------------------------------------------------------

data BinEntry = BE {
   members :: [String], 
   trees   :: [AnnotatedTree]
}
  deriving Show 

-- We index the results of binning by topology-only trees that have their decorations removed.
-- (But we leave the weights on and leave the type as AnnotatedTree so as to acces Ord.)
type BinResults = M.Map AnnotatedTree BinEntry

-- Takes labeled trees, classifies labels into equivalence classes.
--binthem :: [(String, NewickTree BranchLen)] -> M.Map (NewickTree ()) BinEntry
binthem :: [(String, NewickTree BranchLen)] -> BinResults
binthem ls = binthem_normed normalized
 where
  normalized = map (\ (lab,tree) -> (lab, normalize $ annotateWLabLists tree)) ls

-- This version accepts trees that are already normalized:
binthem_normed :: [(String, AnnotatedTree)] -> BinResults
binthem_normed normalized = 
--   foldl (\ acc (lab,tree) -> M.insertWith update tree (BE{ members=[lab] }) acc)
   foldl (\ acc (lab,tree) -> M.insertWith update (anonymize_annotated tree) (BE [lab] [tree]) acc)
	 M.empty normalized
	 --(map (mapSnd$ fmap (const ())) normalized) -- still need to STRIP them
 where 
 --(++)
-- update new old = BE{ members= (members new ++ members old) }
 update new old = BE (members new ++ members old) (trees new ++ trees old)
 --strip = fmap (const ())

-- Remove branch lengths and labels but leave weights
anonymize_annotated :: AnnotatedTree -> AnnotatedTree
anonymize_annotated = fmap (\ (bl, w, labs) -> (0, w, []))


----------------------------------------------------------------------------------------------------
-- Other tools and algorithms.
----------------------------------------------------------------------------------------------------

-- Extract all edges connected to a particular node in every tree.  Return branch lengths.
all_edge_weights lab trees = 
     concat$ map (loop []) trees
  where 
 loop acc (NTLeaf len name) | lab == name = len:acc
 loop acc (NTLeaf _ _)                    = acc
 loop acc (NTInterior _ ls) = foldl loop acc ls


----------------------------------------------------------------------------------------------------
-- Bitvector based normalization.
----------------------------------------------------------------------------------------------------

-- TODO: This approach is probably faster. Give it a try.

{-
int NumberOfSetBits(int i)
{
    i = i - ((i >> 1) & 0x55555555);
    i = (i & 0x33333333) + ((i >> 2) & 0x33333333);
    return ((i + (i >> 4) & 0xF0F0F0F) * 0x1010101) >> 24;
}

int __builtin_popcount (unsigned int x);
-}


----------------------------------------------------------------------------------------------------
-- Visualization with GraphViz and FGL:
----------------------------------------------------------------------------------------------------

-- First we need to be able to convert to FGL graphs:
toGraph :: AnnotatedTree -> Gr String Double
toGraph tree = run_ G.empty $ loop tree
  where
 loop (NTLeaf (l,_,_) name) = 
    do let str = fromLabel name
       G.insMapNodeM str
       return str
 loop (NTInterior (l,_,sorted) ls) =
    do let bigname = concat$ map fromLabel sorted
       names <- mapM loop ls
       G.insMapNodeM bigname
       mapM_ (\x -> insMapEdgeM (bigname, x, 0.0)) names
       return bigname

-- This version uses the tree nodes themselves as labels.
toGraph2 :: AnnotatedTree -> Gr AnnotatedTree Double
toGraph2 tree = run_ G.empty $ loop tree
  where
 loop node@(NTLeaf _  _) =  
    do G.insMapNodeM node 
       return ()
 loop node@(NTInterior (l,_,sorted) ls) =
    do mapM_ loop ls
       G.insMapNodeM node 
       mapM_ (\x -> insMapEdgeM (node, x, fst3$ get_dec x)) ls
       return ()


-- The channel retuned will carry a single message to signal
-- completion of the subprocess.
drawNewickTree :: String -> AnnotatedTree -> IO (Chan (), AnnotatedTree)
drawNewickTree title tree =
  do chan <- newChan
     let dot = dotNewickTree title (1.0 / avg_branchlen [tree])
	                     tree
	 runit = do runGraphvizCanvas default_cmd dot Xlib
		    writeChan chan ()
     --str <- prettyPrint d
     --putStrLn$ "Generating the following graphviz tree:\n " ++ str
     forkIO runit
       --do runGraphvizCanvas Dot dot Xlib; return ()
       
     return (chan, tree)

--default_cmd = TwoPi -- Totally ignores edge lengths.
default_cmd = Neato

-- Show a float without scientific notation:
myShowFloat :: Double -> String
-- showFloat weight = showEFloat (Just 2) weight ""
myShowFloat fl = printf "%.4f" fl


dotNewickTree :: String -> Double -> AnnotatedTree -> DotGraph G.Node
dotNewickTree title edge_scale tree = 
    --trace ("EDGE SCALE: " ++ show edge_scale) $
    graphToDot myparams graph
 where 
  graph = toGraph2 tree
  myparams :: GraphvizParams AnnotatedTree Double () AnnotatedTree
  myparams = defaultParams { globalAttributes= [GraphAttrs [Gv.Label (StrLabel title)]],
			     fmtNode= nodeAttrs, fmtEdge= edgeAttrs }
  nodeAttrs :: (Int,AnnotatedTree) -> [Attribute]
  nodeAttrs (num,node) =
    let children = get_children node in 
    [ Gv.Label (StrLabel$ concat$ map fromLabel$ thd3$ get_dec node)
    , Shape (if null children then {-PlainText-} Ellipse else PointShape)
    , Style [SItem Filled []]
    ]

  -- TOGGLE:
  --  edgeAttrs (_,_,weight) = [ArrowHead noArrow, Len (weight * edge_scale + bump), Gv.Label (StrLabel$ show (weight))]
  edgeAttrs (_,_,weight) = 
                           let draw_weight = compute_draw_weight weight edge_scale in
                           --trace ("EDGE WEIGHT "++ show weight ++ " drawn at "++ show draw_weight) $
			   [ArrowHead noArrow, Gv.Label (StrLabel$ myShowFloat weight)] ++ -- TEMPTOGGLE
			   --[ArrowHead noArrow, Gv.Label (StrLabel$ show draw_weight)] ++ -- TEMPTOGGLE
			    if weight == 0.0
			    then [Color [X11Color Red], Len minlen]
			    else [Len draw_weight]
  minlen = 0.7
  maxlen = 3.0
  compute_draw_weight w scale = 
    let scaled = (abs w) * scale + minlen in 
    -- Don't draw them too big or it gets annoying:
    (min scaled maxlen)


-- This version shows the ordered/rooted structure of the normalized tree.
dotNewickTree_debug :: String -> AnnotatedTree -> DotGraph G.Node
dotNewickTree_debug title tree = graphToDot myparams graph
 where 
  graph = toGraph2 tree
  myparams :: GraphvizParams AnnotatedTree Double () AnnotatedTree
  myparams = defaultParams { globalAttributes= [GraphAttrs [Gv.Label (StrLabel title)]],
			     fmtNode= nodeAttrs, fmtEdge= edgeAttrs }
  nodeAttrs :: (Int,AnnotatedTree) -> [Attribute]
  nodeAttrs (num,node) =
    let children = get_children node in 
    [ Gv.Label (if null children 
  	        then StrLabel$ concat$ map fromLabel$ thd3$ get_dec node
	        else RecordLabel$ take (length children) $ 
                                  -- This will leave interior nodes unlabeled:
	                          map (PortName . PN) $ map show [1..]
		                  -- This version gives some kind of name to interior nodes:
--	                          map (\ (i,ls) -> LabelledTarget (PN$ show i) (fromLabel$ head ls)) $ 
--                                       zip [1..] (map (thd3 . get_dec) children)
               )
    , Shape Record
    , Style [SItem Filled []]
    ]

  edgeAttrs (num1,num2,weight) = 
    let node1 = fromJust$ lab graph num1 
	node2 = fromJust$ lab graph num2 	
	ind = fromJust$ elemIndex node2 (get_children node1)
    in [TailPort$ LabelledPort (PN$ show$ 1+ind) (Just South)]



----------------------------------------------------------------------------------------------------
-- Utilities and UNIT TESTING
----------------------------------------------------------------------------------------------------

parseNewick :: String -> B.ByteString -> NewickTree BranchLen
parseNewick file input = 
  runB file newick_parser $
  B.filter (not . isSpace) input

runB :: Show a => String -> Parser a -> B.ByteString -> a
runB file p input = case (parse p "" input) of
	         Left err -> error ("parse error in file "++ show file ++" at "++ show err)
		 Right x  -> x

runPr prs str = print (run prs str)
run p input = runB "<unknown>" p (B.pack input)

errortest :: t -> IO ()
errortest x = 
   --() ~=?
    handle (\ (e::SomeException) -> return ()) $ 
      do evaluate x
         assertFailure "test was expected to throw an error"

cnt :: NewickTree a -> Int
cnt (NTLeaf _ _) = 1
cnt (NTInterior _ ls) = 1 + sum (map cnt ls)

tr1 = run newick_parser "(A:0.1,B:0.2,(C:0.3,D:0.4):0.5);"
tr1draw = drawNewickTree "tr1"$ annotateWLabLists tr1
tr1dot = putStrLn$ prettyPrint' $ dotNewickTree "" 1.0 $ annotateWLabLists tr1


norm = normalize . annotateWLabLists . run newick_parser
norm2 = normalize . annotateWLabLists . parseNewick "test"
tests = 
  let 
      ntl s = NTLeaf 0.0 (toLabel s)
  in 
  test 
   [ "test name"   ~: "foo" ~=?  run name "foo"
   , "test number" ~:  3.3  ~=?  run number "3.3"
   , "test number" ~:  3.0  ~=?  run number "3"
   , "test number" ~:  -3.0 ~=?  run number "-3"

   , "leaf"     ~: ntl "A" ~=?  run leaf    "A"
   , "subtree"  ~: ntl "A" ~=?  run subtree "A"

   -- These are not allowed:
   , "null branchset" ~: errortest$ run branchset ""

   , "internal" ~: NTInterior 0.0 [ntl "A"] ~=?  run internal "(A);"

   , "example: no nodes are named"  ~: NTInterior 0 [ntl "", ntl "",NTInterior 0 [ntl "", ntl ""]]
				   ~=? run newick_parser "(,,(,));"
   , "example: leaf nodes are named" ~: 6 ~=?  cnt (run newick_parser "(A,B,(C,D));")
   , "example: all nodes are named"  ~: 6 ~=?  cnt (run newick_parser "(A,B,(C,D)E)F;")

   , "example: all but root node have a distance to parent"  ~: 6 ~=? cnt (run newick_parser "(:0.1,:0.2,(:0.3,:0.4):0.5);")
   , "example: all have a distance to parent"              ~: 6 ~=? cnt (run newick_parser "(:0.1,:0.2,(:0.3,:0.4):0.5):0.6;")
   , "example: distances and leaf names (popular)"         ~: 6 ~=? cnt tr1
   , "example: distances and all names"                    ~: 6 ~=? cnt (run newick_parser "(A:0.1,B:0.2,(C:0.3,D:0.4)E:0.5)F;")
   , "example: a tree rooted on a leaf node (rare)"        ~: 6 ~=? cnt (run newick_parser "((B:0.2,(C:0.3,D:0.4)E:0.5)F:0.1)A;")

   , "merge" ~: [1,2,3,4,5,6] ~=? merge [1,3,5] [2,4,6]

   , "demerge" ~: [2,4,6] ~=? demerge [1,2,3,4,5,6] [1,3,5]
   , "demerge" ~: [1,3,5] ~=? demerge [1,2,3,4,5,6] [2,4,6]

   , "annotateWLabLists" ~: 
     --NTInterior (0.0,[A,B,C,D]) [NTLeaf (0.1,[A]) A,NTLeaf (0.2,[B]) B,NTInterior (0.5,[C,D]) [NTLeaf (0.3,[C]) C,NTLeaf (0.4,[D]) D]]
        map toLabel ["A","B","C","D"] -- ORD on atoms is expensive... it must use the whole string.
    ~=? get_label_list (annotateWLabLists tr1)

   -- Make sure that all of these normalize to the same thing.
   , "normalize1" ~: "(C, D, E, (A, B))" ~=?  show (pPrint$ norm "(A,(C,D,E),B);")
   , "normalize2" ~: "(C, D, E, (A, B))" ~=?  show (pPrint$ norm "((C,D,E),B,A);")
   , "normalize2" ~: "(C, D, E, (A, B))" ~=?  show (pPrint$ norm "(D,E,C,(B,A));")

   -- Here's an example from the rhizobia datasetsthat that caused my branch-sorting to fail.
   , "normalize3" ~:  "(((BB, BJ)), (MB, ML), (RE, (SD, SM)))" 
		 ~=? show (pPrint$ norm2 (B.pack "(((ML,MB),(RE,(SD,SM))),(BB,BJ));"))

-- "((BB: 2.691831, BJ: 1.179707): 0.000000, ((ML: 0.952401, MB: 1.020319): 0.000000, (RE: 2.031345, (SD: 0.180786, SM: 0.059988): 0.861187): 0.717913): 0.000000);"


   , "dotConversion" ~: True ~=? 100 < length (prettyPrint' $ dotNewickTree "" 1.0$ norm "(D,E,C,(B,A));") -- 444

   
   , "phbin: these 3 trees should fall in the same category" ~: 
      1 ~=? (length $ M.toList $
             binthem [("one",   run newick_parser "(A,(C,D,E),B);"),
 		      ("two",   run newick_parser "((C,D,E),B,A);"),
		      ("three", run newick_parser "(D,E,C,(B,A));")])

   ]

run_tests = runTestTT tests
t = run_tests

   
----------------------------------------------------------------------------------------------------
-- Driver to put the pieces together (parse, normalize, bin)
----------------------------------------------------------------------------------------------------

-- Due to the number of configuration options for the driver, we pack them into a record:
data PhyBinConfig = 
  PBC { verbose :: Bool
      , num_taxa :: Int
      , name_hack :: Label -> Label
      , output_dir :: String
      , inputs :: [String]
      , do_graph :: Bool
      , do_draw :: Bool
      }

default_phybin_config = 
 PBC { verbose = False
      , num_taxa = error "must be able to determine the number of taxa expected in the dataset"
      , name_hack = id -- Default, no transformation of leaf-labels
      , output_dir = "./"
      , inputs = []
      , do_graph = False
      , do_draw = False
     }


driver :: PhyBinConfig -> IO ()
driver PBC{..} =
 do 
    --------------------------------------------------------------------------------
    -- First, find out where we are and open the files:
    --------------------------------------------------------------------------------
    cd <- getCurrentDirectory 
    --putStrLn$ "PHYBIN RUNNING IN DIRECTORY: "++ cd

    all :: [[String]] <- forM inputs $ \ path -> do
      exists <- file_exists path 

      --stat   <- if exists then getFileStatus path else return (error "internal invariant")
      -- [2010.09.23] This is no longer really necessary:
      if not exists then do 
	 putStr$ "Input not a file/directory, assuming wildcard, using 'find' for expansion"
	 entries <- HSH.run$ "find " ++ path	 
	 putStrLn$ "("++show (length entries)++" files found):  "++ show path
	 return entries
       else do
	 isdir <- is_directory path
	 reg  <- is_regular_file path
	 if isdir then do 
	    putStr$ "Input is a directory, reading all regular files contained "
	    children <- getDirectoryContents path
	    filtered <- filterM is_regular_file $ map (combine path) children
	    putStrLn$ "("++show (length filtered)++" regular files found):  "++ show path
	    return$ filtered
          else if reg then do 
	    return [path]
	  else error$ "phybin: Unhandled input path: " ++ path

    let files = concat all -- take 10 $ concat all
	num_files = length files

    putStrLn$ "Parsing "++show num_files++" Newick tree files."
    --putStrLn$ "\nFirst ten \n"++ concat (map (++"\n") $ map show $ take 10 files)

    --------------------------------------------------------------------------------
    -- Next, parse the files and do error checking and annotation.
    --------------------------------------------------------------------------------
    --
    -- results contains: num-nodes, parsed, warning-files
    results <- forM files $ \ file -> 
      do --stat <- getFileStatus file		 
	 reg <- is_regular_file file
	 if not reg then return (0,[],[(-1, file)]) else do

           h <- openFile file ReadMode 
	   bstr <- B.hGetContents h

           -- Clip off the first three characters:
           let 
	       parsed = map_labels name_hack $ parseNewick file bstr
	       annot  = annotateWLabLists parsed
	       normal = normalize annot
	       weight = get_weight annot

           -- TEMPTOGGLE
	   when False $ do putStrLn$ "DRAWING TREE";  drawNewickTree "Annotated" annot;  drawNewickTree "Normalized" normal
			   putStrLn$ "WEIGHTS OF NORMALIZED' CHILDREN: "++ show (map get_weight$ get_children normal)

           if not$ weight == num_taxa
	    then do --putStrLn$ "\n WARNING: file contained an empty or single-node tree: "++ show file
 		    when verbose$ putStrLn$ "\n WARNING: file contained unexpected number of leaves ("
					    ++ show weight ++"): "++ show file
		    return (0,[], [(weight, file)])
	    else do 
	     when verbose$ putStr "."

	     --evaluate$ deepseq$ runB newick_parser bstr
	     --evaluate$ cnt$ runB newick_parser bstr
	     num <- evaluate$ cnt parsed
	     --num <- evaluate$ cnt normal

	     hClose h
	     --return$ (num, [normal])
	     return$ (num, [parsed], [])

    putStrLn$ "\nNumber of input trees: " ++ show num_files
    putStrLn$ "Number of VALID trees (correct # of leaves/taxa): " ++ show (length$ concat$ map snd3 results)
    putStrLn$ "Total tree nodes contained in valid trees: "++ show (sum$ map fst3 results)

    --------------------------------------------------------------------------------
    -- Do the actual binning:
    --------------------------------------------------------------------------------

    putStrLn$ "Creating equivalence classes (bins), bin sizes:"

    let classes = --binthem_normed$ zip files $ concat$ map snd3 results
	          binthem$  zip files $ concat$ map snd3 results
	binlist = reverse $ sortBy (compare `on` fst3) $
		  map (\ (tr,ls) -> (length (members ls), tr,ls)) $ M.toList classes
	numbins = length binlist
	taxa = S.unions$ map (S.fromList . all_labels . snd3) binlist
	warnings = concat $ map thd3 results
	base i size = combine output_dir ("bin" ++ show i ++"_"++ show size)

    ----------------------------------------
    -- TEST, TEMPTOGGLE: print out edge weights :
    -- forM_ (map snd3 results) $ \parsed -> do 
    --    let weights = all_edge_weights (head$ S.toList taxa) parsed
    --    trace ("weights of "++ show parsed ++" "++ show weights) $
    --      return ()
    -- exitSuccess
    ----------------------------------------

    --------------------------------------------------------------------------------
    -- Finally, produce all the required outputs.
    --------------------------------------------------------------------------------

    forM_ binlist $ \ (len, tr, _) -> do
       when (len > 1) $ -- Omit that long tail of single element classes...
          putStrLn$ "  "++ show (pPrint tr) ++" members: "++ show len

    putStrLn$ "\nTotal unique taxa ("++ show (S.size taxa) ++"):  "++ 
	      show (sep $ map (text . fromLabel) $ S.toList taxa)

    putStrLn$ "Final number of tree bins: "++ show (M.size classes)
    forM_ (zip [1..] binlist) $ \ (i, (size, tr, bentry)) -> do
       --putStrLn$ ("  WRITING " ++ combine output_dir ("bin" ++ show i ++"_"++ show size ++".txt"))
       writeFile (base i size ++".txt") (concat$ map (++"\n") (members bentry))
       writeFile (base i size ++".tr")  (show (pPrint tr) ++ ";\n")
--       writeFile (base i size ++".rawtree")  (show tr ++ ";\n") -- TempToggle

    when (not$ null warnings) $
	writeFile (combine output_dir "bin_WARNINGS.txt")
		  ("This file was generated to record all of the files which WERE NOT incorporated successfully into the results.\n" ++
		   "Each of these files had some kind of problem, likely one of the following:\n"++
		   "  (1) a mismatched number of taxa (leaves) in the tree relative to the rest of the dataset\n"++
		   "  (2) a file that could not be read.\n"++
		   "  (3) a file that could not be parsed.\n\n"++
		   concat (map (\ (n,file) -> 
				(if n == -1 
				 then "Not a regular/readable file: "++ file 
				 else "Wrong number of taxa ("++ show n ++"): "++ file)
				++"\n") 
		           warnings))
    putStrLn$ "Wrote contents of each bin to bin<N>_<binsize>.txt"
    putStrLn$ "Wrote representative trees to bin<N>_<binsize>.tr"

    when (do_graph) $ do
      forM_ (zip [1..] binlist) $ \ (i, (size, tr, bentry)) -> do
         let 
             dot = dotNewickTree ("bin #"++ show i) (1.0 / avg_branchlen (trees bentry))
		                 --(annotateWLabLists$ fmap (const 0) tr)
		                 -- TEMP FIXME -- using just ONE representative tree:
		                 (--trace ("WEIGHTED: "++ show (head$ trees bentry)) $ 
		                  --(head$ trees bentry) )
				  (avg_trees$ trees bentry))
	 when (size > 1 || numbins < 100) $ do 
	   runGraphvizCommand default_cmd dot Pdf (base i size ++ ".pdf")
	   return ()
      putStrLn$ "Wrote visual representations of trees to bin<N>_<binsize>.pdf"

    --putStrLn$ "Wrote representative tree to bin<N>_<binsize>.tr"
    putStrLn$ "Finished."
    --------------------------------------------------------------------------------
    -- End driver
    --------------------------------------------------------------------------------


-- Average branch length across several trees.
avg_branchlen :: [AnnotatedTree] -> Double
avg_branchlen ls = fst total / snd total
  where
   total = sum_ls $ map sum_tree ls
   sum_ls ls = (sum$ map fst ls, sum$ map snd ls)
{-
   sum_tree (NTLeaf (l,_,_) _) | l < 0 = 
       trace ("!!! GOT NEGATIVE BRANCH LENGTH: "++ show l) $
       (0,0)
-}
   sum_tree (NTLeaf (0,_,_) _)         = (0,0)
   sum_tree (NTLeaf (l,_,_) _)         = (abs l,1)
   sum_tree (NTInterior (l,_,_) ls) = 
       let (x,y) = sum_ls$ map sum_tree ls in
       if l == 0 then (x, y) else ((abs l) + x, 1+y)

{-
nonzero_blens :: AnnotatedTree -> Int
nonzero_blens node = 
    let children = sum $ map nonzero_blens $ get_children node in
    if (fst3 $ get_dec node) == 0 
    then children
    else children + 1
-}

-- Come up with an average tree from a list of isomorphic trees.
-- This comes up with some blending of edge lengths.
avg_trees :: [AnnotatedTree] -> AnnotatedTree
avg_trees ls = --summed -- TEMPTOGGLE
    fmap (\ (blen,w,ls) -> (blen / count, w,ls)) summed
 where
  summed = foldl1 sum_2trees ls
  count = fromIntegral$ length ls

  sum_2trees a b = case (a,b) of
    (NTLeaf (l1,w,ls) nm, NTLeaf (l2,_,_) _) -> NTLeaf (l1+l2,w,ls) nm
    (NTInterior (l1,w,ls) ls1, 
     NTInterior (l2,_,_)  ls2) -> NTInterior (l1+l2,w,ls) $ 
                                             map (uncurry sum_2trees) $ zip ls1 ls2
    _ -> error "avg_trees: applied to non-isomorphic trees"


bump = 0.00001 -- for DIRTY HACKS

{- 
 ----------------------------------------
 PARSING TIMING TEST:
 ----------------------------------------

 Compiling this with GHC 6.12 on my laptop -O2...
 It takes 0.043s startup to parse ten files.
 And 0.316 seconds to parse 2648.. so we can say that's almost all time spent parsing/building/traversing.
 (All nodes summed to 14966)
  (The tested version uses Strings for labels... not Atoms)

 Comparing against the original mzscheme version (with Racket 5.0)
 with default optimization (there's no obvious -O2), well the
 generated .exe has a ~0.5 second startup time overhead...
   0.881 seconds total to do the parsing, or about 380ms just for parsing.
   But that doesn't do the counting!
   Ok, this mzscheme version is in a messed up state at this point, but hacking
   it to do a count (and it gets a different one.. 12319), I get 0.882 seconds real time, 
   that is neglibly more.
   
 If anything parsec should be at a disadvantage because of the lack of
 a preprocessing phase to generate the FSM...

 Btw, switching node labels over to Atoms made no difference. (But
 didn't slow down at least.)  We wouldn't expect this to save anything
 on the construction side... parsec still allocates/builds the strings
 before we intern them.

 -}


----------------------------------------------------------------------------------------------------
-- MAIN script: Read command line options and call the program.
----------------------------------------------------------------------------------------------------

-- Note: ORDER is important here, we process options in this order:
data Flag 
    = Verbose  
    | Version 
    | Output String
    | NumTaxa Int
    | NullOpt
    | Graph | Draw
    | Force 
    | View
    | TabDelimited Int Int

    | SelfTest

    | NameCutoff String
    | NamePrefix Int
    | NameTable String  -- Must come after Cutoff/Prefix

  deriving (Show, Eq, Ord) -- ORD is really used.

parseTabDelim str = 
  TabDelimited 8 9
    
options :: [OptDescr Flag]
options =
     [ Option ['v']     ["verbose"] (NoArg Verbose)    "print WARNINGS and other information (recommended at first)"
     , Option ['V']     ["version"] (NoArg Version)    "show version number"

     , Option ['o']     ["output"]  (ReqArg Output "DIR")  "set directory to contain all output files (default \"./\")"

     , Option []     ["selftest"]   (NoArg SelfTest)   "run internal unit tests"


{- -- TODO: FIXME: IMPLEMENT THIS:
     , Option []        []          (NoArg NullOpt)  ""
     , Option ['t']     ["tabbed"]  (ReqArg parseTabDelim "NUM1:NUM2")$  "assume the input is a tab-delimited file with gene names \n"++
		                                                        "in column NUM1 and Newick trees in NUM2"
-}

     , Option []        []          (NoArg NullOpt)  ""
     , Option []        []  (NoArg$ error "internal problem")  "----------------------------- Visualization --------------------------------"
     , Option ['g']     ["graphbins"] (NoArg Graph)  "use graphviz to produce .dot and .pdf output files named bin1.*, bin2.*, etc"
     , Option ['d']     ["drawbins"]  (NoArg Draw)   "like -g, but open GUI windows to show a tree for each bin"

     , Option ['w']     ["view"]    (NoArg View)$  "for convenience, \"view mode\" simply displays input Newick files without binning" 

     , Option []        []          (NoArg NullOpt)  ""
     , Option []        []  (NoArg$ error "internal problem")  "--------------------------- Handling taxa names ----------------------------"
--     , Option ['n']     ["numtaxa"] (ReqArg (NumTaxa . read) "NUM") "expect NUM taxa for this dataset (otherwise it will guess)"
-- ^^ TODO, FIXME: The "guessing" part doesn't actually work yet, implement it!!
     , Option ['n']     ["numtaxa"] (ReqArg (NumTaxa . read) "NUM") "expect NUM taxa for this dataset"

{- -- TODO: FIXME: IMPLEMENT THIS:
     , Option ['f']     ["force"]   (NoArg Force)    "force phybin to consume and bin trees with different numbers of taxa"
-}
     , Option []        []          (NoArg NullOpt)  ""
     , Option ['p']     ["nameprefix"]  (ReqArg (NamePrefix . read) "NUM") $ 
		  "Leaf names in the input Newick trees are usually gene names, not taxa.\n"++
    	    	  "It is typical to extract taxa names from genes.  This option extracts a\n"++
                  "prefix of NUM characters to serve as the taxa name."

     , Option []        []          (NoArg NullOpt)  ""
     , Option ['s']     ["namesep"]   (ReqArg NameCutoff "STR") $  --"\n"++
		  "An alternative to --nameprefix, STR provides a set of delimeter characters,\n"++
                  "for example '-' or '0123456789'.  The taxa name is then a variable-length\n"++
		  "prefix of each gene name up to but not including any character in STR."

     , Option []        []          (NoArg NullOpt)  ""
     , Option ['m']     ["namemap"]  (ReqArg NameTable "FILE") $ 
		  "Even once prefixes are extracted it may be necessary to use a lookup table\n"++
		  "to compute taxa names, e.g. if multiple genes/plasmids map onto one taxa.\n"++
		  "This option specifies a text file with find/replace entries of the form\n"++
		  "\"<string> <taxaname>\", which are applied AFTER -s and -p."
     ]

usage = "\nUsage: phybin [OPTION...] files or directories...\n\n"++

        "PhyBin takes Newick tree files as input and produces, at minimum, files of the form\n"++
        "binXX_YY.tr, each containing a list of input file paths that fall into that bin.\n\n"++

	"USAGE NOTES: Currently phybin ignores input trees with the wrong number of taxa.\n"++
	"If given a directory as input phybin will assume all contained files are Newick trees.\n\n"++

	"\nOptions include:\n"
defaultErr errs = error $ "ERROR!\n" ++ (concat errs ++ usageInfo usage options)


main = 
  do argv <- getArgs 

     (opts,files) <- 
       case getOpt Permute options argv of
	 (o,n,[]  ) -> return (o,n)
         (_,_,errs) -> defaultErr errs

     let process_opt cfg opt = case opt of 
	   NullOpt -> return cfg
	   Verbose -> return cfg { verbose= True } 
	   Version -> do putStrLn$ "phybin version "++phybin_version; exitSuccess

	   SelfTest -> do run_tests; exitSuccess

	   Output s -> return cfg { output_dir= s }

	   NumTaxa n -> return cfg { num_taxa= n }
	   Graph     -> return cfg { do_graph= True } 
	   Draw	     -> return cfg { do_draw = True } 
	   View      -> return cfg -- Handled below

	   TabDelimited _ _ -> error "tabbed option not handled yet"
	   Force            -> error "force option not handled yet"


	   NameCutoff str -> let set = S.fromList str 
				 new = toLabel . takeWhile (not . flip S.member set) . fromLabel
			     in return cfg { name_hack = new . name_hack cfg }
	   NamePrefix n   -> let new = toLabel . (take n) . fromLabel 
			     in return cfg { name_hack = new . name_hack cfg }

           -- This should always be after cutoff/prefix:
	   NameTable file -> do reader <- name_table_reader file
				return cfg { name_hack = reader . name_hack cfg }


     config <- foldM process_opt default_phybin_config{ inputs=files } 
	             (sort opts) -- NOTE: options processed in sorted order.

     when (null files) $ do
	defaultErr ["No file arguments!\n"]

     if View `elem` opts 
      then view_graphs config
      --else driver config{ name_hack= name_hack_legionella }
      else driver config

view_graphs :: PhyBinConfig -> IO ()
view_graphs PBC{..} = 
           do chans <- forM inputs $ \ file -> do 
                putStrLn$ "Drawing "++ file ++"...\n"
		str <- B.readFile file
		putStrLn$ "Parsed: " ++ (B.unpack str)
 	        (chan, tr) <- drawNewickTree file $ 
			      annotateWLabLists$ 
			      map_labels name_hack $ 
			      parseNewick file str
	        return chan
	      forM_ chans readChan 
	      return ()


--------------------------------------------------------------------------------
-- Every dataset it seems needs a new hack on the names.
name_table_reader file = 
  do contents <- readFile file
     let mp = M.fromList $ 
	      map (\ls -> case ls of 
		           [a,b] -> (a,b)
		           _ -> error$ "Each line of "++file++"must contain two whitespace free strings: "++ unwords ls) $ 
	      filter (not . null) $
	      map words $ 
	      lines contents
	    
     return$ 
       \ name_to_hack -> 
	   case M.lookup name_to_hack mp of -- Could use a trie
	     Just x -> x
	     Nothing -> name_to_hack

temp = driver default_phybin_config{ num_taxa=7, inputs=["../datasets/test.tr"] }

----------------------------------------------------------------------------------------------------
-- TODO: expose a command line argument for testing.
-- The below test exposed my normalization bug relating to "deriving Ord".
-- I need to transform it into one or more proper unit tests.

main_test = 
 withArgs ["-w","~/newton_and_newton_local/datasets/yersinia/yersinia_trees/111.dnd","-m","../datasets/yersinia/name_table_hack_yersinia.txt"]
	  main 

--a :: AnnotatedTree
-- annotateWLabLists$ 
a :: NewickTree Double
a = set_dec 1 $ 
    NTInterior () [NTInterior () [NTLeaf () "RE",NTInterior () [NTLeaf () "SD",NTLeaf () "SM"]],NTInterior () [NTLeaf () "BB",NTLeaf () "BJ"],NTInterior () [NTLeaf () "MB",NTLeaf () "ML"]]

b :: NewickTree Double
b = set_dec 1 $ 
    NTInterior () [NTInterior () [NTLeaf () "BB",NTLeaf () "BJ"],NTInterior () [NTLeaf () "MB",NTLeaf () "ML"],NTInterior () [NTLeaf () "RE",NTInterior () [NTLeaf () "SD",NTLeaf () "SM"]]]

ls = [("a",a),("b",b)]

-- This is one:
num_binned = M.size $ binthem ls

a_ =  ("980.dnd",NTInterior 0.0 [NTInterior 5.697e-2 [NTLeaf 3.95e-2 "SM",NTLeaf 5.977e-2 "SD"],NTLeaf 0.13143 "RE",NTInterior 0.13899 [NTInterior 9.019e-2 [NTLeaf 0.11856 "BB",NTLeaf 0.13592 "BJ"],NTInterior 0.13194 [NTLeaf 0.19456 "MB",NTLeaf 0.16603 "ML"]]])

b_ = ("999.dnd",NTInterior 0.0 [NTInterior 6.527e-2 [NTInterior 0.13734 [NTLeaf 2.975e-2 "SM",NTLeaf 3.002e-2 "SD"],NTLeaf 0.18443 "RE"],NTInterior 6.621e-2 [NTLeaf 0.16184 "MB",NTLeaf 0.15233 "ML"],NTInterior 0.23143 [NTLeaf 9.192e-2 "BB",NTLeaf 0.10125 "BJ"]])

-- But THIS is two:  ack!
num2 = M.size $ binthem [a_,b_]

-- Here's the test that breaks things:
a_norm = normalize (annotateWLabLists$ snd a_)

--a_norm = NTInterior (0.13899,7,["BB","BJ","MB","ML","RE","SD","SM"]) [NTInterior (0.0,3,["RE","SD","SM"]) [NTLeaf (0.13143,1,["RE"]) "RE",NTInterior (5.697e-2,2,["SD","SM"]) [NTLeaf (5.977e-2,1,["SD"]) "SD",NTLeaf (3.95e-2,1,["SM"]) "SM"]],NTInterior (9.019e-2,2,["BB","BJ"]) [NTLeaf (0.11856,1,["BB"]) "BB",NTLeaf (0.13592,1,["BJ"]) "BJ"],NTInterior (0.13194,2,["MB","ML"]) [NTLeaf (0.19456,1,["MB"]) "MB",NTLeaf (0.16603,1,["ML"]) "ML"]]

b_norm = NTInterior (0.0,7,["BB","BJ","MB","ML","RE","SD","SM"]) [NTInterior (0.23143,2,["BB","BJ"]) [NTLeaf (9.192e-2,1,["BB"]) "BB",NTLeaf (0.10125,1,["BJ"]) "BJ"],NTInterior (6.621e-2,2,["MB","ML"]) [NTLeaf (0.16184,1,["MB"]) "MB",NTLeaf (0.15233,1,["ML"]) "ML"],NTInterior (6.527e-2,3,["RE","SD","SM"]) [NTLeaf (0.18443,1,["RE"]) "RE",NTInterior (0.13734,2,["SD","SM"]) [NTLeaf (3.002e-2,1,["SD"]) "SD",NTLeaf (2.975e-2,1,["SM"]) "SM"]]]

d1 = drawNewickTree "" a_norm
d2 = drawNewickTree "" b_norm

d1_ = forkIO $ do runGraphvizCanvas Dot (dotNewickTree_debug "" a_norm) Xlib; return ()
d2_ = forkIO $ do runGraphvizCanvas Dot (dotNewickTree_debug "" b_norm) Xlib; return ()
