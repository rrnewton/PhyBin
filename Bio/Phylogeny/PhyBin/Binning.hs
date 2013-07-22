{-# LANGUAGE ScopedTypeVariables, RecordWildCards, TypeSynonymInstances, CPP #-}
{-# LANGUAGE NamedFieldPuns #-}
{-# LANGUAGE OverloadedStrings #-}
{-# OPTIONS_GHC -fwarn-incomplete-patterns #-}
{-# OPTIONS_GHC -fwarn-unused-imports #-}

-- | This module contains the code that does the tree normalization and binning.

module Bio.Phylogeny.PhyBin.Binning
       ( -- * Binning and normalization
         binthem, normalize, annotateWLabLists, 
         deAnnotate, OneCluster(..), BinResults, StrippedTree,
         -- * Utilities and unit tests
         get_weight, unitTests
       )
       where

import qualified Data.Foldable as F
import           Data.Function       (on)
import           Data.List           (delete, minimumBy, sortBy, insertBy, intersperse, sort)
import           Data.Maybe          (fromMaybe, catMaybes)
import qualified Data.ByteString.Lazy.Char8 as B
import qualified Data.Map                   as M
import qualified Data.Set                   as S
import           Control.Monad       (forM, forM_, filterM, when, unless)
import           Control.Exception   (evaluate)
import           Control.Applicative ((<$>),(<*>))
import           Control.Concurrent  (Chan)
import           System.FilePath     (combine)
import           System.Directory    (doesFileExist, doesDirectoryExist,
                                      getDirectoryContents, getCurrentDirectory)
import           System.IO           (openFile, hClose, IOMode(ReadMode))
import           System.Process      (system)
import           System.Exit         (ExitCode(..))
import           Test.HUnit          ((~:),(~=?),Test,test)
import qualified HSH 

-- For vizualization:
import           Text.PrettyPrint.HughesPJClass hiding (char, Style)
import           Bio.Phylogeny.PhyBin.CoreTypes
import           Bio.Phylogeny.PhyBin.Parser (parseNewicks, parseNewick)
import           Bio.Phylogeny.PhyBin.PreProcessor (collapseBranches)
import           Bio.Phylogeny.PhyBin.Visualize (dotToPDF, dotNewickTree, viewNewickTree)
import           Bio.Phylogeny.PhyBin.RFDistance
import           Bio.Phylogeny.PhyBin.Util

-- Turn on for extra invariant checking:
debug :: Bool
debug = True

----------------------------------------------------------------------------------------------------
-- Normal form for unordered, unrooted trees
----------------------------------------------------------------------------------------------------

-- The basic idea is that what we *want* is the following, 
--   ROOT: most balanced point
--   ORDER: sorted in increasing subtree weight

-- But that's not quite good enough.  There are ties to break.  To do
-- that we fall back on the (totally ordered) leaf labels.

--------------------------------------------------------------------------------


-- Our sorting criteria for the children of interior nodes:
compare_childtrees :: AnnotatedTree -> AnnotatedTree -> Ordering
compare_childtrees node1 node2 = 
    case (subtreeWeight $ get_dec node1) `compare` (subtreeWeight $ get_dec node2) of 
     -- Comparisons on atoms cause problems WRT to determinism between runs if parallelism is introduced.
     -- Can consider it an optimization for the serial case perhaps:
--     EQ -> case map deAtom (get_label_list node1) `compare` 
--	        map deAtom (get_label_list node2) of
     EQ -> case (sortedLabels$ get_dec node1) `compare` (sortedLabels$ get_dec node2) of
            EQ -> error$ "Internal invariant broken.  These two children have equal ordering priority:\n" 
		  ++ "Pretty printing:\n  "
		  ++ show (pPrint node1) ++ "\n  " ++ show (pPrint node2)
		  ++ "\nFull data structures:\n  "
		  ++ show (node1) ++ "\n  " ++ show (node2)
	    x  -> x
     x -> x


-- | This is it, here's the routine that transforms a tree into normal form.
--   This relies HEAVILY on lazy evaluation.
normalize :: AnnotatedTree -> AnnotatedTree
normalize tree = snd$ loop tree Nothing
 where 

  add_context dec Nothing  = dec
  add_context dec (Just c) = add_weight dec c

  -- loop: Walk over the tree, turning it inside-out in the process.
  -- Inputs: 
  --    1. node: the NewickTree node to process ("us")
  --    2. context: all nodes connected through the parent, "flipped" as though *we* were root
  --                The "flipped" part has ALREADY been normalized.
  -- Outputs: 
  --    1. new node
  --    3. the best candidate root anywhere under this subtree
  loop :: AnnotatedTree -> Maybe (AnnotatedTree) -> (AnnotatedTree, AnnotatedTree)
  loop node ctxt  = case node of
    NTLeaf (StandardDecor _ _ w sorted) _name -> 
	(node, 
	 -- If the leaf becomes the root... we could introduce another node:
	 NTInterior (add_context (StandardDecor 0 Nothing w sorted) ctxt) $
	            (verify_sorted "1" id$ maybeInsert compare_childtrees ctxt [node])

	 -- It may be reasonable to not support leaves becoming root.. that changes the number of nodes!
	            --error "normalize: leaf becoming root not currently supported."
	)
    
    NTInterior dec ls -> 
     let 
         -- If this node becomes the root, the parent becomes one of our children:
         inverted = NTInterior inverted_dec inverted_children
	 inverted_dec      = add_context dec ctxt
         inverted_children = verify_sorted "2" id$ maybeInsert compare_childtrees ctxt newchildren

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
verify_sorted :: (Show a, Pretty a) => String -> (a -> AnnotatedTree) -> [a] -> [a]
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
tt :: AnnotatedTree
tt = normalize $ annotateWLabLists $ snd$ parseNewick M.empty id "" "(A,(C,D,E),B);"

norm4 :: FullTree StandardDecor
norm4 = norm "((C,D,E),B,A);"

norm5 :: AnnotatedTree
norm5 = normalize$ annotateWLabLists$ snd$ parseNewick M.empty id "" "(D,E,C,(B,A));"


----------------------------------------------------------------------------------------------------
-- Equivalence classes on Trees:
----------------------------------------------------------------------------------------------------

-- | When binning, the members of a OneCluster are isomorphic trees.  When clustering
-- based on robinson-foulds distance they are merely similar trees.
data OneCluster a = OneCluster {
   members  :: [TreeName], 
   bintrees :: [NewickTree a]
   -- TODO : Get rid of this in favor of a simple list of FullTree..
}
  deriving Show 

-- | Ignore metadata (but keep weights) for the purpose of binning
type StrippedTree = NewickTree Int

-- | Index the results of binning by topology-only stripped trees
--   that have their decorations removed.
type BinResults a = M.Map StrippedTree (OneCluster a)

-- | The binning function.
--   Takes labeled trees, classifies labels into equivalence classes.
binthem :: [FullTree DefDecor] -> BinResults StandardDecor
binthem ls = binthem_normed normalized
 where
  normalized = map (\ (FullTree n lab tree) ->
                     FullTree n lab (normalize $ annotateWLabLists tree)) ls

-- | This version accepts trees that are already normalized:
binthem_normed :: [FullTree StandardDecor] -> BinResults StandardDecor
binthem_normed normalized = 
--   foldl (\ acc (lab,tree) -> M.insertWith update tree (OneCluster{ members=[lab] }) acc)
   foldl (\ acc (FullTree treename _ tree) ->
           M.insertWith update (anonymize_annotated tree) (OneCluster [treename] [tree]) acc)
	 M.empty normalized
	 --(map (mapSnd$ fmap (const ())) normalized) -- still need to STRIP them
 where 
-- update new old = OneCluster{ members= (members new ++ members old) }
 update new old = OneCluster (members new ++ members old) (bintrees new ++ bintrees old)
 --strip = fmap (const ())

-- | For binning. Remove branch lengths and labels but leave weights.
anonymize_annotated :: AnnotatedTree -> StrippedTree
anonymize_annotated = fmap (\ (StandardDecor bl bs w labs) -> w)


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
-- General helper/utility functions:
----------------------------------------------------------------------------------------------------

merge :: Ord a => [a] -> [a] -> [a]
merge [] ls = ls
merge ls [] = ls
merge l@(a:b) r@(x:y) = 
  if a < x
  then a : merge b r
  else x : merge y l 

-- Set subtraction for sorted lists:
demerge :: (Ord a, Show a) => [a] -> [a] -> [a]
demerge ls [] = ls
demerge [] ls = error$ "demerge: first list did not contain all of second, remaining: " ++ show ls
demerge (a:b) r@(x:y) = 
  case a `compare` x of
   EQ -> demerge b y
   LT -> a : demerge b r 
   GT -> error$ "demerge: element was missing from first list: "++ show x

-- maybeCons :: Maybe a -> [a] -> [a]
-- maybeCons Nothing  ls = ls
-- maybeCons (Just x) ls = x : ls

maybeInsert :: (a -> a -> Ordering) -> Maybe a -> [a] -> [a]
maybeInsert _  Nothing  ls = ls
maybeInsert fn (Just x) ls = insertBy fn x ls

-- | Add the metadata that is used for binning
annotateWLabLists :: NewickTree DefDecor -> AnnotatedTree
annotateWLabLists tr = case tr of 
  NTLeaf (bs,bl) n      -> NTLeaf (StandardDecor bl bs 1 [n]) n
  NTInterior (bs,bl) ls -> 
      let children = map annotateWLabLists ls in 
      NTInterior (StandardDecor bl bs
                  (sum $ map (subtreeWeight . get_dec) children)
		  (foldl1 merge $ map (sortedLabels . get_dec) children))
		 children

----------------------------------------------------------------------------------------------------
-- Simple Helper Functions
----------------------------------------------------------------------------------------------------

-- | Take the extra annotations away.  Inverse of `annotateWLabLists`.
deAnnotate :: FullTree StandardDecor -> FullTree DefDecor
deAnnotate (FullTree a b tr) = FullTree a b (fmap (\ (StandardDecor bl bs _ _) -> (bs,bl)) tr)

-- Number of LEAVES contained in subtree:
get_weight :: AnnotatedTree -> Int
get_weight = subtreeWeight . get_dec

-- Sorted list of leaf labels contained in subtree
get_label_list :: AnnotatedTree -> [Label]
get_label_list   = sortedLabels . get_dec

add_weight :: StandardDecor -> AnnotatedTree -> StandardDecor
add_weight (StandardDecor l1 bs1 w1 sorted1) node  = 
  let (StandardDecor _ bs2 w2 sorted2) = get_dec node in 
  StandardDecor l1 ((+) <$> bs1 <*> bs2) (w1+w2) (merge sorted1 sorted2)

-- Remove the influence of one subtree from the metadata of another.
subtract_weight :: StandardDecor -> AnnotatedTree -> StandardDecor
subtract_weight (StandardDecor l1 bs1 w1 sorted1) node =  
  let (StandardDecor _ bs2 w2 sorted2) = get_dec node in 
  StandardDecor l1 ((-) <$> bs1 <*> bs2) (w1-w2) (demerge sorted1 sorted2)
	

----------------------------------------------------------------------------------------------------
-- UNIT TESTING
----------------------------------------------------------------------------------------------------

tre1 :: (LabelTable, NewickTree DefDecor)
tre1 = parseNewick M.empty id "" "(A:0.1,B:0.2,(C:0.3,D:0.4):0.5);"

tre1draw :: IO (Chan (), FullTree StandardDecor)
tre1draw = viewNewickTree "tre1"$ (FullTree "" (fst tre1) (annotateWLabLists (snd tre1)))

tre1dot :: IO ()
tre1dot = print $ dotNewickTree "" 1.0 $ (FullTree "" (fst tre1) (annotateWLabLists$ snd tre1))

norm :: String -> FullTree StandardDecor
norm = norm2 . B.pack

norm2 :: B.ByteString -> FullTree StandardDecor
norm2 bstr = FullTree "" tbl (normalize $ annotateWLabLists tr)
  where
    (tbl,tr) = parseNewick M.empty id "test" bstr

unitTests :: Test
unitTests = 
  let ntl s = NTLeaf (Nothing,0.0) s
      (tbl,tre1_) = tre1
  in 
  test 
   [ 
     "merge" ~: [1,2,3,4,5,6] ~=? merge [1,3,5] [2,4,6::Int]
   , "demerge" ~: [2,4,6] ~=? demerge [1,2,3,4,5,6] [1,3,5::Int]
   , "demerge" ~: [1,3,5] ~=? demerge [1,2,3,4,5,6] [2,4,6::Int]
   , "annotateWLabLists" ~: 
     ["A","B","C","D"] -- ORD on atoms is expensive... it must use the whole string.
      ~=? map (tbl M.!) 
          (sortedLabels (get_dec (annotateWLabLists tre1_)))

   -- Make sure that all of these normalize to the same thing.
   , "normalize1" ~: "(C, D, E, (A, B))" ~=?  show (displayDefaultTree $ deAnnotate$ norm "(A,(C,D,E),B);")
   , "normalize2" ~: "(C, D, E, (A, B))" ~=?  show (pPrint$ norm "((C,D,E),B,A);")
   , "normalize2" ~: "(C, D, E, (A, B))" ~=?  show (pPrint$ norm "(D,E,C,(B,A));")

   -- Here's an example from the rhizobia datasetsthat that caused my branch-sorting to fail.
   , "normalize3" ~:  "(((BB, BJ)), (MB, ML), (RE, (SD, SM)))" 
		 ~=? show (pPrint$ norm2 (B.pack "(((ML,MB),(RE,(SD,SM))),(BB,BJ));"))

-- "((BB: 2.691831, BJ: 1.179707): 0.000000, ((ML: 0.952401, MB: 1.020319): 0.000000, (RE: 2.031345, (SD: 0.180786, SM: 0.059988): 0.861187): 0.717913): 0.000000);"

   , "phbin: these 3 trees should fall in the same category" ~: 
      1 ~=? (length $ M.toList $
             binthem $ snd $ 
              parseNewicks id [("one",  "(A,(C,D,E),B);"), 
                               ("two",  "((C,D,E),B,A);"),
                               ("three","(D,E,C,(B,A));")]
             -- [("one",   snd$parseNewick M.empty id "" "(A,(C,D,E),B);"),
             --  ("two",   snd$parseNewick M.empty id "" "((C,D,E),B,A);"),
             --  ("three", snd$parseNewick M.empty id "" "(D,E,C,(B,A));")]
            )
      
   , "dotConversion" ~: True ~=? 100 < length (show $ dotNewickTree "" 1.0$ norm "(D,E,C,(B,A));") -- 444
   ]

--------------------------------------------------------------------------------
