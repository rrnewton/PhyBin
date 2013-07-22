{-# LANGUAGE ScopedTypeVariables, RecordWildCards, TypeSynonymInstances, CPP #-}
{-# LANGUAGE NamedFieldPuns #-}
{-# LANGUAGE OverloadedStrings #-}
{-# OPTIONS_GHC -fwarn-incomplete-patterns #-}
{-# OPTIONS_GHC -fwarn-unused-imports #-}

-- | This module contains the code that does the tree normalization and binning.
--   It's the heart of the prgoram.

module Bio.Phylogeny.PhyBin
       ( driver, binthem, normalize, annotateWLabLists, unitTests, acquireTreeFiles,
         deAnnotate )
       where

import qualified Data.Foldable as F
import           Data.Function       (on)
import           Data.List           (delete, minimumBy, sortBy)
import           Data.Maybe          (fromMaybe)
import           Data.Either         (partitionEithers)
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
import           Bio.Phylogeny.PhyBin.Parser (parseNewick, parseNewicks)
import           Bio.Phylogeny.PhyBin.PreProcessor (collapseBranches)
import           Bio.Phylogeny.PhyBin.Visualize (dotToPDF, dotNewickTree, viewNewickTree)
import           Bio.Phylogeny.PhyBin.RFDistance
import           Bio.Phylogeny.PhyBin.Binning
import           Bio.Phylogeny.PhyBin.Util

-- Turn on for extra invariant checking:
debug :: Bool
debug = True

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

    files <- acquireTreeFiles inputs
    bl <- doesDirectoryExist output_dir
    unless bl $ do
      c <- system$ "mkdir -p "++output_dir
      case c of
        ExitSuccess   -> return ()
        ExitFailure c -> error$"Could not create output directory. 'mkdir' command failed with: "++show c
    
    putStrLn$ "Parsing "++show (length files)++" Newick tree files."
    --putStrLn$ "\nFirst ten \n"++ concat (map (++"\n") $ map show $ take 10 files)

    --------------------------------------------------------------------------------
    -- Next, parse the files and do error checking and annotation.
    --------------------------------------------------------------------------------

    (goodFiles,warnings) <- fmap partitionEithers $
      forM files $ \ file -> do
           reg <- is_regular_file file
  	   if reg
             then return$ Left file
             else return$ Right file
    let num_files = length goodFiles
    bstrs <- mapM B.readFile goodFiles
    let (labelTab, fulltrees) = parseNewicks name_hack (zip files bstrs)

    --------------------------------------------------------------------------------

    case branch_collapse_thresh of 
      Just thr -> putStrLn$" !+ Collapsing branches of length less than "++show thr
      Nothing  -> return ()

    let do_one :: FullTree DefDecor -> IO (Int, [NewickTree DefDecor], [(Int, String)])
        do_one (FullTree treename lblAcc parsed) = do 
           let 
               collapser _ _  = (Nothing,0)
               pruned = case branch_collapse_thresh of 
                         Nothing  -> parsed
                         Just thr -> collapseBranches ((< thr) . snd) collapser parsed
--	       annot  = annotateWLabLists pruned
--	       normal = normalize annot
--	       weight = get_weight annot
               numL   = numLeaves pruned
               
           -- TEMPTOGGLE
	   -- when False $ do putStrLn$ "DRAWING TREE"
           --                 viewNewickTree "Annotated"  (FullTree file lblAcc' annot)
           --                 viewNewickTree "Normalized" (FullTree file lblAcc' normal)
	   --      	   putStrLn$ "WEIGHTS OF NORMALIZED' CHILDREN: "++
           --                       show (map get_weight$ get_children normal)

           if numL /= num_taxa
	    then do --putStrLn$ "\n WARNING: file contained an empty or single-node tree: "++ show file
 		    when verbose$ putStrLn$ "\n WARNING: tree contained unexpected number of leaves ("
					    ++ show numL ++"): "++ treename
		    return (0, [], [(numL, treename)])
	    else do 
	     when verbose$ putStr "."
	     return$ (numL, [pruned], [])
             
    -- results contains: label-maps, num-nodes, parsed, warning-files             
    results <- mapM do_one fulltrees
    let (counts::[Int], validtrees, pairs::[[(Int, String)]]) = unzip3 results

    putStrLn$ "\nNumber of input tree files: " ++ show num_files
    when (length warnings > 0) $
      putStrLn$ "Number of bad/unreadable input tree files: " ++ show (length warnings)
    putStrLn$ "Number of VALID trees (correct # of leaves/taxa): " ++ show (length$ concat validtrees)
    putStrLn$ "Total tree nodes contained in valid trees: "++ show (sum counts)

    --------------------------------------------------------------------------------
    -- Do the actual binning:
    --------------------------------------------------------------------------------

    putStrLn$ "Creating equivalence classes (bins)..."

    let classes = --binthem_normed$ zip files $ concat$ map snd3 results
	          binthem$  zip files $ concat validtrees
	binlist = reverse $ sortBy (compare `on` fst3) $
		  map (\ (tr,ls) -> (length (members ls), tr, ls)) $ M.toList classes
	numbins = length binlist
        taxa :: S.Set Int
	taxa = S.unions$ map (S.fromList . all_labels . snd3) binlist
	warnings = concat pairs
	base i size = combine output_dir ("bin" ++ show i ++"_"++ show size)

        binsizes = map fst3 binlist

    putStrLn$ " [outcome] "++show numbins++" bins found, "++show (length$ takeWhile (>1) binsizes)
             ++" non-singleton, top bin sizes: "++show(take 10 binsizes)
    putStrLn$"  ALL bin sizes, excluding singletons:"
    forM_ (zip [1..] binlist) $ \ (ind, (len, tr, BE{bintrees})) -> do
       when (len > 1) $ -- Omit that long tail of single element classes...
          putStrLn$show$ 
           hcat [text ("  * bin#"++show ind++", members "++ show len ++", "), 
                 vcat [text ("avg bootstraps "++show (get_bootstraps$ avg_trees bintrees)++", "),
                       text "all: " <> pPrint (filter (not . null) $ map get_bootstraps bintrees)]]


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

    putStrLn$ "\nTotal unique taxa ("++ show (M.size labelTab) ++"):\n"++ 
	      show (nest 2 $ sep $ map text $ M.elems labelTab)

    putStrLn$ "Final number of tree bins: "++ show (M.size classes)
    let avgs = map (avg_trees . bintrees . thd3) binlist
    forM_ (zip3 [1::Int ..] binlist avgs) $ \ (i, (size, _tr, bentry), avgTree) -> do
       let fullTr = FullTree "fixthis" fixme_HERE avgTree
         
       --putStrLn$ ("  WRITING " ++ combine output_dir ("bin" ++ show i ++"_"++ show size ++".txt"))
       writeFile (base i size ++".txt") (concat$ map (++"\n") (members bentry))
       -- writeFile (base i size ++".tr")  (show (pPrint tr) ++ ";\n")
       -- Printing the average tree instead of the stripped cannonical one:
       when debug$ do
         writeFile (base i size ++".dbg") (show (pPrint avgTree) ++ "\n")
       writeFile   (base i size ++".tr")  (show (displayDefaultTree$ deAnnotate fullTr) ++ ";\n") -- FIXME

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
      forM_ (zip3 [1::Int ..] binlist avgs) $ \ (i, (size, _tr, bentry), avgTree) -> do
         let fullTr = FullTree "fixthis" fixme_HERE avgTree
	 when (size > 1 || numbins < 100) $ do 
           let dot = dotNewickTree ("bin #"++ show i) (1.0 / avg_branchlen (bintrees bentry))
                                   --(annotateWLabLists$ fmap (const 0) tr)
                                   -- TEMP FIXME -- using just ONE representative tree:
                                   ( --trace ("WEIGHTED: "++ show (head$ trees bentry)) $ 
                                     --(head$ trees bentry) 
 				    fullTr)
	   _ <- dotToPDF dot (base i size ++ ".pdf")
	   return ()
      putStrLn$ "[finished] Wrote visual representations of trees to bin<N>_<binsize>.pdf"

    --putStrLn$ "Wrote representative tree to bin<N>_<binsize>.tr"
    putStrLn$ "Finished."
    --------------------------------------------------------------------------------
    -- End driver
    --------------------------------------------------------------------------------

fixme_HERE = M.empty

-- Monadic mapAccum
mapAccumM :: Monad m => (acc -> x -> m (acc,y)) -> acc -> [x] -> m (acc,[y])
mapAccumM fn acc ls = F.foldrM fn' (acc,[]) ls
  where
    fn' x (acc,ls) = do (acc',y) <- fn acc x
                        return (acc',y:ls)

fst3 :: (t, t1, t2) -> t
snd3 :: (t, t1, t2) -> t1
thd3 :: (t, t1, t2) -> t2
fst3 (a,_,_) = a
snd3 (_,b,_) = b
thd3 (_,_,c) = c


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
