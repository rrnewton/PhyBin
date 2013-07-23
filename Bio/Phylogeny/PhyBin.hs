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
import qualified Data.Vector                 as V
import qualified Data.Vector.Unboxed         as U
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
import qualified Data.Clustering.Hierarchical as C

-- For vizualization:
import           Text.PrettyPrint.HughesPJClass hiding (char, Style)
import           Bio.Phylogeny.PhyBin.CoreTypes
import           Bio.Phylogeny.PhyBin.Parser (parseNewick, parseNewicks)
import           Bio.Phylogeny.PhyBin.PreProcessor (collapseBranches)
import           Bio.Phylogeny.PhyBin.Visualize (dotToPDF, dotNewickTree, viewNewickTree)
import           Bio.Phylogeny.PhyBin.RFDistance
import           Bio.Phylogeny.PhyBin.Binning
import           Bio.Phylogeny.PhyBin.Util

import Debug.Trace

-- Turn on for extra invariant checking:
debug :: Bool
debug = True

----------------------------------------------------------------------------------------------------

-- | Driver to put all the pieces together (parse, normalize, bin)
driver :: PhyBinConfig -> IO ()
driver PBC{ verbose, num_taxa, name_hack, output_dir, inputs,
            do_graph, branch_collapse_thresh,
            dist_thresh, clust_mode, print_rfmatrix } =
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
        ExitSuccess     -> return ()
        ExitFailure cde -> error$"Could not create output directory. 'mkdir' command failed with: "++show cde
    
    putStrLn$ "Parsing "++show (length files)++" Newick tree files."
    --putStrLn$ "\nFirst ten \n"++ concat (map (++"\n") $ map show $ take 10 files)

    --------------------------------------------------------------------------------
    -- Next, parse the files and do error checking and annotation.
    --------------------------------------------------------------------------------

    (goodFiles,warnings1) <- fmap partitionEithers $
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

    let do_one :: FullTree DefDecor -> IO (Int, [FullTree DefDecor], [(Int, String)])
        do_one (FullTree treename lblAcc parsed) = do 
           let 
               collapser _ _  = (Nothing,0)
               pruned = case branch_collapse_thresh of 
                         Nothing  -> parsed
                         Just thr -> collapseBranches ((< thr) . snd) collapser parsed
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
	     return$ (numL, [FullTree treename lblAcc pruned], [])

    results <- mapM do_one fulltrees
    let (counts::[Int], validtreess, pairs::[[(Int, String)]]) = unzip3 results
    let validtrees = concat validtreess
        warnings2 = concat pairs
        
    putStrLn$ "\nNumber of input tree files: " ++ show num_files
    when (length warnings2 > 0) $
      putStrLn$ "Number of bad/unreadable input tree files: " ++ show (length warnings2)
    putStrLn$ "Number of VALID trees (correct # of leaves/taxa): " ++ show (length validtrees)
    putStrLn$ "Total tree nodes contained in valid trees: "++ show (sum counts)

    --------------------------------------------------------------------------------
    -- Next, dispatch on the mode and do the actual clustering or binning.
    --------------------------------------------------------------------------------

    classes <- case clust_mode of
      BinThem         -> doBins validtrees 
      ClusterThem lnk -> do
        (mat, dendro) <- doCluster lnk validtrees
        case print_rfmatrix of
          False -> return ()
          True -> do -- treeFiles <- acquireTreeFiles files
                     -- let fn f = do raw <- B.readFile f
                     --               let ls = map (`B.append` (B.pack ";")) $ 
                     --                        B.splitWith (== ';') raw
                     --               return (map (f,) ls)
                     -- trees0 <- concat <$> mapM fn treeFiles
                     -- FIXME: no name_hack here:
                     -- let (lbls, trees) = parseNewicks id trees0 
                     -- putStrLn$ "Read trees! "++show (length trees)
                     -- putStrLn$ "Taxa: "++show (pPrint lbls)
                     -- putStrLn$ "First tree: "++show (displayDefaultTree (head trees))
                     printDistMat mat
        
        case dist_thresh of
          Nothing -> error "Fully hierarchical cluster output is not finished!  Use --editdist."
          Just dstThresh -> do
            let dstThresh' = fromIntegral dstThresh 
            let loop br@(C.Branch dist left right)
                  | dist >= dstThresh' = loop left ++ loop right
                  | otherwise          = [flattenDendro br]
                loop br@(C.Leaf _)     = [flattenDendro br]
            -- Flatten out the dendogram:
            return (clustsToMap $ loop dendro)
    binlist <- reportClusts clust_mode classes
        
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

    unless (null warnings1 && null warnings2) $
	writeFile (combine output_dir (filePrefix ++ "_WARNINGS.txt"))
		  ("This file was generated to record all of the files which WERE NOT incorporated successfully into the results.\n" ++
		   "Each of these files had some kind of problem, likely one of the following:\n"++
		   "  (1) a mismatched number of taxa (leaves) in the tree relative to the rest of the dataset\n"++
		   "  (2) a file that could not be read.\n"++
		   "  (3) a file that could not be parsed.\n\n"++
		   concat (map (\ file -> "Not a regular/readable file: "++ file++"\n")
		           warnings1) ++ 
		   concat (map (\ (n,file) ->
                                 "Wrong number of taxa ("++ show n ++"): "++ file++"\n")
		           warnings2))

    case clust_mode of
      BinThem         -> outputBins binlist output_dir do_graph
      ClusterThem lnk -> outputClusters binlist output_dir do_graph

    --putStrLn$ "Wrote representative tree to bin<N>_<binsize>.tr"
    putStrLn$ "Finished."
    --------------------------------------------------------------------------------
    -- End driver
    --------------------------------------------------------------------------------

-- filePrefix = "bin"
filePrefix = "cluster"    

--------------------------------------------------------------------------------
-- Driver helpers:
--------------------------------------------------------------------------------

-- doBins :: [FullTree DefDecor] -> t1 -> t2 -> IO BinResults
doBins :: [FullTree DefDecor] -> IO (BinResults StandardDecor)
doBins validtrees = do 
    putStrLn$ "Creating equivalence classes (bins)..."
    let classes = --binthem_normed$ zip files $ concat$ map snd3 results
	          binthem validtrees
    return (classes)

doCluster :: C.Linkage -> [FullTree a] -> IO (DistanceMatrix, C.Dendrogram (FullTree a))
doCluster linkage validtrees = do
  putStrLn$ "Clustering using method "++show linkage
  let nwtrees  = map nwtree validtrees
      numtrees = length validtrees 
      mat      = distanceMatrix nwtrees
      ixtrees  = zip [0..] validtrees
      dist (i,t1) (j,t2) | j == i     = 0
--                         | i == numtrees-1 = 0 
                         | j < i      = fromIntegral ((mat V.! i) U.! j)
                         | otherwise  = fromIntegral ((mat V.! j) U.! i)
      dist1 a b = trace ("Taking distance between "++show (fst a, fst b)) $ dist a b
      dendro = fmap snd $ C.dendrogram linkage ixtrees dist
  -- putStrLn$ "Got numtrees ...  "++show numtrees
  -- putStrLn$ "Got the distance matrix ...  "++show (V.length mat)
  return (mat,dendro)
  
reportClusts :: ClustMode -> BinResults StandardDecor -> IO [(Int, StrippedTree, OneCluster StandardDecor)]
reportClusts mode classes = do 
    let binlist = reverse $ sortBy (compare `on` fst3) $
		  map (\ (tr,OneCluster ls) -> (length ls, tr, OneCluster ls)) $ M.toList classes
        taxa :: S.Set Int
	taxa = S.unions$ map (S.fromList . all_labels . snd3) binlist
        binsizes = map fst3 binlist

    putStrLn$ " [outcome] "++show (length binlist)++" bins found, "++show (length$ takeWhile (>1) binsizes)
             ++" non-singleton, top bin sizes: "++show(take 10 binsizes)
    putStrLn$"  ALL bin sizes, excluding singletons:"
    forM_ (zip [1..] binlist) $ \ (ind, (len, tr, OneCluster ftrees)) -> do
       when (len > 1) $ -- Omit that long tail of single element classes...
          putStrLn$show$
           hcat [text ("  * bin#"++show ind++", members "++ show len ++", "), 
                 case mode of
                   BinThem -> vcat [text ("avg bootstraps "++show (get_bootstraps$ avg_trees$ map nwtree ftrees)++", "),
                                    text "all: " <> pPrint (filter (not . null) $ map (get_bootstraps . nwtree) ftrees)]
                   ClusterThem _ -> hcat [] ]
    return binlist

-- | Convert a flat list of clusters into a map from individual trees to clusters.
clustsToMap :: [OneCluster StandardDecor] -> BinResults StandardDecor
clustsToMap clusts = F.foldl' fn M.empty clusts
  where
    fn acc theclust@(OneCluster ftrs) =
      F.foldl' (fn2 theclust) acc ftrs
    fn2 theclust acc (FullTree{nwtree}) =
      M.insert (anonymize_annotated nwtree) theclust acc

flattenDendro :: (C.Dendrogram (FullTree DefDecor)) -> OneCluster StandardDecor
flattenDendro dendro =
  case dendro of
    C.Leaf (FullTree{treename,labelTable,nwtree}) ->
      OneCluster [FullTree treename labelTable (annotateWLabLists nwtree)]
    C.Branch _ left right ->
      -- TODO: fix quadratic append
      flattenDendro left `appendClusts` flattenDendro right
 where
   appendClusts (OneCluster ls1) (OneCluster ls2) = OneCluster (ls1++ls2)
   

--------------------------------------------------------------------------------

outputClusters binlist output_dir do_graph = do
    
    return ()


outputBins binlist output_dir  do_graph = do
    let numbins = length binlist
    let base i size = combine output_dir (filePrefix ++ show i ++"_"++ show size) 
    let avgs = map (avg_trees . map nwtree . clustMembers . thd3) binlist
    forM_ (zip3 [1::Int ..] binlist avgs) $ \ (i, (size, _tr, OneCluster ftrees), avgTree) -> do
       let FullTree fstName labs _ = head ftrees
           fullAvgTr = FullTree fstName labs avgTree
         
       --putStrLn$ ("  WRITING " ++ combine output_dir (filePrefix ++ show i ++"_"++ show size ++".txt"))
       writeFile (base i size ++".txt") (concat$ map ((++"\n") . treename) ftrees)
       -- writeFile (base i size ++".tr")  (show (pPrint tr) ++ ";\n")
       -- Printing the average tree instead of the stripped cannonical one:
       when debug$ do
         writeFile (base i size ++".dbg") (show (pPrint avgTree) ++ "\n")
       writeFile   (base i size ++".tr")  (show (displayDefaultTree$ deAnnotate fullAvgTr) ++ ";\n") -- FIXME

    putStrLn$ "[finished] Wrote contents of each bin to bin<N>_<binsize>.txt"
    putStrLn$ "           Wrote representative trees to bin<N>_<binsize>.tr"  
    when (do_graph) $ do
      putStrLn$ "Next do the time consuming operation of writing out graphviz visualizations:"
      forM_ (zip3 [1::Int ..] binlist avgs) $ \ (i, (size, _tr, OneCluster membs), avgTree) -> do
         let FullTree fstName labs _ = head membs
             fullAvgTr = FullTree fstName labs avgTree 
    	 when (size > 1 || numbins < 100) $ do 
           let dot = dotNewickTree ("cluster #"++ show i) (1.0 / avg_branchlen (map nwtree membs))
                                   --(annotateWLabLists$ fmap (const 0) tr)
                                   -- TEMP FIXME -- using just ONE representative tree:
                                   ( --trace ("WEIGHTED: "++ show (head$ trees bentry)) $ 
                                     --(head$ trees bentry) 
 				    fullAvgTr)
	   _ <- dotToPDF dot (base i size ++ ".pdf")
	   return ()
      putStrLn$ "[finished] Wrote visual representations of trees to bin<N>_<binsize>.pdf"


--------------------------------------------------------------------------------

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
