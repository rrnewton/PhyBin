{-# LANGUAGE ScopedTypeVariables, RecordWildCards, TypeSynonymInstances, CPP #-}
{-# LANGUAGE NamedFieldPuns #-}
{-# LANGUAGE OverloadedStrings #-}
{-# OPTIONS_GHC -fwarn-incomplete-patterns #-}
{-# OPTIONS_GHC -fwarn-unused-imports #-}

-- | This module contains the code that does the tree normalization and binning.
--   It's the heart of the prgoram.

module Bio.Phylogeny.PhyBin
       ( driver, binthem, normalize, annotateWLabLists, unitTests )
       where

import           Data.Function      (on)
import           Data.List          (delete, minimumBy, sortBy, insertBy, intersperse, sort)
import           Data.Maybe         (fromMaybe)
import qualified Data.ByteString.Lazy.Char8 as B
import qualified Data.Map                   as M
import qualified Data.Set                   as S
import           Control.Monad       (forM, forM_, filterM, when)
import           Control.Exception   (evaluate)
import           Control.Applicative ((<$>),(<*>))
import           Control.Concurrent  (Chan)
import           System.FilePath     (combine)
import           System.Directory    (doesFileExist, doesDirectoryExist,
                                      getDirectoryContents, getCurrentDirectory)
import           System.IO           (openFile, hClose, IOMode(ReadMode))
import           Test.HUnit          ((~:),(~=?),Test,test)
import qualified HSH 

-- For vizualization:
import           Text.PrettyPrint.HughesPJClass hiding (char, Style)
import           Bio.Phylogeny.PhyBin.CoreTypes
import           Bio.Phylogeny.PhyBin.Parser (parseNewick)
import           Bio.Phylogeny.PhyBin.PreProcessor (collapseBranches)
import           Bio.Phylogeny.PhyBin.Visualize (dotToPDF, dotNewickTree, viewNewickTree)

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

-- A common type of tree is "AnnotatedTree", which contains the standard decorator.
type AnnotatedTree = NewickTree StandardDecor

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
  loop :: AnnotatedTree -> Maybe AnnotatedTree -> (AnnotatedTree, AnnotatedTree)
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
verify_sorted :: (Show a, Pretty a) => String -> (a -> NewickTree StandardDecor) -> [a] -> [a]
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
tt = normalize $ annotateWLabLists $ parseNewick "" "(A,(C,D,E),B);"

norm4 :: AnnotatedTree
norm4 = norm "((C,D,E),B,A);"

norm5 :: AnnotatedTree
norm5 = normalize$ annotateWLabLists$ parseNewick "" "(D,E,C,(B,A));"


----------------------------------------------------------------------------------------------------
-- Equivalence classes on Trees:
----------------------------------------------------------------------------------------------------

-- | All the members of a BinEntry are isomorphic trees.
data BinEntry = BE {
   members :: [String], 
   trees   :: [AnnotatedTree]
}
  deriving Show 

-- | Ignore metadata (but keep weights) for the purpose of binning
type StrippedTree = NewickTree Int

-- | Index the results of binning by topology-only stripped trees
--   that have their decorations removed.
type BinResults = M.Map StrippedTree BinEntry

-- | The binning function.
--   Takes labeled trees, classifies labels into equivalence classes.
binthem :: [(String, NewickTree DefDecor)] -> BinResults
binthem ls = binthem_normed normalized
 where
  normalized = map (\ (lab,tree) -> (lab, normalize $ annotateWLabLists tree)) ls


-- | This version accepts trees that are already normalized:
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

-- | Driver to put all the pieces together (parse, normalize, bin)
driver :: PhyBinConfig -> IO ()
driver PBC{ verbose, num_taxa, name_hack, output_dir, inputs, do_graph, branch_collapse_thresh } =
   -- Unused: do_draw
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
    results :: [(Int, [NewickTree DefDecor], [(Int, String)])] <- forM files $ \ file -> 
      do --stat <- getFileStatus file		 
	 reg <- is_regular_file file
	 if not reg then return (0,[],[(-1, file)]) else do

           h <- openFile file ReadMode 
	   bstr <- B.hGetContents h

           -- Clip off the first three characters:
           let 
	       parsed = map_labels name_hack $ parseNewick file bstr
               pruned = case branch_collapse_thresh of 
                          Nothing  -> parsed
                          Just thr -> collapseBranches ((< thr) . snd) 
                                         (error "FINISH THIS COLLAPSER") parsed
	       annot  = annotateWLabLists pruned
	       normal = normalize annot
	       weight = get_weight annot

           -- TEMPTOGGLE
	   when False $ do putStrLn$ "DRAWING TREE";  viewNewickTree "Annotated" annot;  viewNewickTree "Normalized" normal
			   putStrLn$ "WEIGHTS OF NORMALIZED' CHILDREN: "++ show (map get_weight$ get_children normal)

           if not$ weight == num_taxa
	    then do --putStrLn$ "\n WARNING: file contained an empty or single-node tree: "++ show file
 		    when verbose$ putStrLn$ "\n WARNING: file contained unexpected number of leaves ("
					    ++ show weight ++"): "++ show file
		    return (0,[], [(weight, file)])
	    else do 
	     when verbose$ putStr "."

	     num <- evaluate$ treeSize pruned
	     --num <- evaluate$ cnt normal

	     hClose h
	     --return$ (num, [normal])
	     return$ (num, [pruned], [])

    putStrLn$ "\nNumber of input trees: " ++ show num_files
    putStrLn$ "Number of VALID trees (correct # of leaves/taxa): " ++ show (length$ concat$ map snd3 results)
    putStrLn$ "Total tree nodes contained in valid trees: "++ show (sum$ map fst3 results)

    --------------------------------------------------------------------------------
    -- Do the actual binning:
    --------------------------------------------------------------------------------

    putStrLn$ "Creating equivalence classes (bins)..."

    let classes = --binthem_normed$ zip files $ concat$ map snd3 results
	          binthem$  zip files $ concat$ map snd3 results
	binlist = reverse $ sortBy (compare `on` fst3) $
		  map (\ (tr,ls) -> (length (members ls), tr, ls)) $ M.toList classes
	numbins = length binlist
	taxa = S.unions$ map (S.fromList . all_labels . snd3) binlist
	warnings = concat $ map thd3 results
	base i size = combine output_dir ("bin" ++ show i ++"_"++ show size)

    putStrLn$ "  "++show numbins++" bins found.  Bin sizes, excluding singletons:"

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

    forM_ binlist $ \ (len, _, _) -> do
       when (len > 1) $ -- Omit that long tail of single element classes...
          -- putStrLn$ "  "++ show (pPrint tr) ++" members: "++ show len
          putStrLn$ "  * members: "++ show len

    putStrLn$ "\nTotal unique taxa ("++ show (S.size taxa) ++"):\n"++ 
	      show (sep $ map (text . fromLabel) $ S.toList taxa)

    putStrLn$ "Final number of tree bins: "++ show (M.size classes)
    forM_ (zip [1::Int ..] binlist) $ \ (i, (size, _tr, bentry)) -> do
       --putStrLn$ ("  WRITING " ++ combine output_dir ("bin" ++ show i ++"_"++ show size ++".txt"))
       writeFile (base i size ++".txt") (concat$ map (++"\n") (members bentry))
       -- writeFile (base i size ++".tr")  (show (pPrint tr) ++ ";\n")
       -- Printing the average tree instead of the stripped cannonical one:
       let avg = avg_trees$ trees bentry 
       when debug$ do
         writeFile (base i size ++".dbg") (show (pPrint (avg_trees$ trees bentry)) ++ "\n")
       writeFile   (base i size ++".tr")  (show (displayDefaultTree$ deAnnotate$ avg) ++ ";\n")

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
    putStrLn$ "[finished] Wrote contents of each bin to bin<N>_<binsize>.txt"
    putStrLn$ "           Wrote representative trees to bin<N>_<binsize>.tr"  
    when (do_graph) $ do
      putStrLn$ "Next do the time consuming operation of writing out graphviz visualizations:"
      forM_ (zip [1::Int ..] binlist) $ \ (i, (size, _tr, bentry)) -> do
--	 when (size > 1 || numbins < 100) $ do 
           let dot = dotNewickTree ("bin #"++ show i) (1.0 / avg_branchlen (trees bentry))
                                   --(annotateWLabLists$ fmap (const 0) tr)
                                   -- TEMP FIXME -- using just ONE representative tree:
                                   ( --trace ("WEIGHTED: "++ show (head$ trees bentry)) $ 
                                     --(head$ trees bentry) 
 				    (avg_trees$ trees bentry) )
	   _ <- dotToPDF dot (base i size ++ ".pdf")
	   return ()
      putStrLn$ "[finished] Wrote visual representations of trees to bin<N>_<binsize>.pdf"

    --putStrLn$ "Wrote representative tree to bin<N>_<binsize>.tr"
    putStrLn$ "Finished."
    --------------------------------------------------------------------------------
    -- End driver
    --------------------------------------------------------------------------------

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
avg_trees origls = 
     fmap unlift $ 
     foldIsomorphicTrees (foldl1 sumthem) $ 
     map (fmap lift) origls     
  where
    totalCount = fromIntegral$ length origls

    -- Here we do the actual averaging:
    unlift :: TempDecor -> StandardDecor
    unlift (bl, (bs, bootcount), wt, ls) =
      (StandardDecor (bl / totalCount)
                     (case bootcount of
                        0 -> Nothing
                        _ -> Just$ round (fromIntegral bs / fromIntegral bootcount))
                     wt ls)
      
    lift :: StandardDecor -> TempDecor
    lift (StandardDecor bl bs wt ls) = (bl, (fromMaybe 0 bs, countMayb bs), wt, ls)

    sumthem :: TempDecor -> TempDecor -> TempDecor
    sumthem (bl1, (bs1, cnt1), wt, ls)
            (bl2, (bs2, cnt2), _,  _) =
            (bl1+bl2, (bs1 + bs2, cnt1+cnt2), wt, ls)
    countMayb Nothing  = 0
    countMayb (Just _) = 1

-- Used only by avg_trees above...
type TempDecor = (Double, (Int, Int), Int, [Label])

    
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
is_regular_file :: FilePath -> IO Bool
is_regular_file = doesFileExist

is_directory :: FilePath -> IO Bool
is_directory = doesDirectoryExist 

file_exists :: FilePath -> IO Bool
file_exists path = 
  do f <- doesFileExist path
     d <- doesDirectoryExist path
     return (f || d)

----------------------------------------------------------------------------------------------------
-- General helper/utility functions:
----------------------------------------------------------------------------------------------------


fst3 :: (t, t1, t2) -> t
snd3 :: (t, t1, t2) -> t1
thd3 :: (t, t1, t2) -> t2
fst3 (a,_,_) = a
snd3 (_,b,_) = b
thd3 (_,_,c) = c


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

-- annotateWLabLists :: NewickTree BranchLen -> AnnotatedTree
annotateWLabLists :: NewickTree DefDecor -> AnnotatedTree
annotateWLabLists tr = case tr of 
  NTLeaf (bs,bl) n      -> NTLeaf (StandardDecor bl bs 1 [n]) n
  NTInterior (bs,bl) ls -> 
      let children = map annotateWLabLists ls in 
      NTInterior (StandardDecor bl bs
                  (sum $ map (subtreeWeight . get_dec) children)
		  (foldl1 merge $ map (sortedLabels . get_dec) children))
		 children

-- | Take the extra annotations away.  Inverse of `annotateWLabLists`.
deAnnotate :: AnnotatedTree -> NewickTree DefDecor 
deAnnotate = fmap (\ (StandardDecor bl bs _ _) -> (bs,bl))


-- Number of LEAVES contained in subtree:
get_weight :: NewickTree StandardDecor -> Int
get_weight = subtreeWeight . get_dec

-- Sorted list of leaf labels contained in subtree
get_label_list :: NewickTree StandardDecor -> [Label]
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

tre1 :: NewickTree DefDecor
tre1 = parseNewick "" "(A:0.1,B:0.2,(C:0.3,D:0.4):0.5);"

tre1draw :: IO (Chan (), AnnotatedTree)
tre1draw = viewNewickTree "tre1"$ annotateWLabLists tre1

tre1dot :: IO ()
tre1dot = print $ dotNewickTree "" 1.0 $ annotateWLabLists tre1


norm :: String -> AnnotatedTree
norm = norm2 . B.pack

norm2 :: B.ByteString -> AnnotatedTree
norm2 = normalize . annotateWLabLists . parseNewick "test"

unitTests :: Test
unitTests = 
  let 
      ntl s = NTLeaf (Nothing,0.0) (toLabel s)
  in 
  test 
   [ 
     "merge" ~: [1,2,3,4,5,6] ~=? merge [1,3,5] [2,4,6::Int]
   , "demerge" ~: [2,4,6] ~=? demerge [1,2,3,4,5,6] [1,3,5::Int]
   , "demerge" ~: [1,3,5] ~=? demerge [1,2,3,4,5,6] [2,4,6::Int]
   , "annotateWLabLists" ~: 
     --NTInterior (0.0,[A,B,C,D]) [NTLeaf (0.1,[A]) A,NTLeaf (0.2,[B]) B,NTInterior (0.5,[C,D]) [NTLeaf (0.3,[C]) C,NTLeaf (0.4,[D]) D]]
        map toLabel ["A","B","C","D"] -- ORD on atoms is expensive... it must use the whole string.
    ~=? sortedLabels (get_dec (annotateWLabLists tre1))

   -- Make sure that all of these normalize to the same thing.
   , "normalize1" ~: "(C, D, E, (A, B))" ~=?  show (pPrint$ norm "(A,(C,D,E),B);")
   , "normalize2" ~: "(C, D, E, (A, B))" ~=?  show (pPrint$ norm "((C,D,E),B,A);")
   , "normalize2" ~: "(C, D, E, (A, B))" ~=?  show (pPrint$ norm "(D,E,C,(B,A));")

   -- Here's an example from the rhizobia datasetsthat that caused my branch-sorting to fail.
   , "normalize3" ~:  "(((BB, BJ)), (MB, ML), (RE, (SD, SM)))" 
		 ~=? show (pPrint$ norm2 (B.pack "(((ML,MB),(RE,(SD,SM))),(BB,BJ));"))

-- "((BB: 2.691831, BJ: 1.179707): 0.000000, ((ML: 0.952401, MB: 1.020319): 0.000000, (RE: 2.031345, (SD: 0.180786, SM: 0.059988): 0.861187): 0.717913): 0.000000);"

   , "phbin: these 3 trees should fall in the same category" ~: 
      1 ~=? (length $ M.toList $
             binthem [("one",   parseNewick "" "(A,(C,D,E),B);"),
 		      ("two",   parseNewick "" "((C,D,E),B,A);"),
		      ("three", parseNewick "" "(D,E,C,(B,A));")])
      
   , "dotConversion" ~: True ~=? 100 < length (show $ dotNewickTree "" 1.0$ norm "(D,E,C,(B,A));") -- 444
   ]

--------------------------------------------------------------------------------
-- No longer used:
--------------------------------------------------------------------------------

-- runPr :: Show a => Parser a -> String -> IO ()
-- runPr prs str = print (run prs str)


-- commasep :: Pretty a => [a] -> Doc
-- commasep ls = sep (intersperse (text ", ") $ map pPrint ls)

-- bump :: Double
-- bump = 0.00001 -- for DIRTY HACKS
