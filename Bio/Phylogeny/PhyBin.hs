{-# LANGUAGE ScopedTypeVariables, RecordWildCards, TypeSynonymInstances, CPP #-}
{-# LANGUAGE NamedFieldPuns #-}
{-# LANGUAGE OverloadedStrings #-}
{-# OPTIONS_GHC -fwarn-incomplete-patterns #-}
{-# OPTIONS_GHC -fwarn-unused-imports #-}

-- | This module contains the code that does the tree normalization and binning.
--   It's the heart of the program.

module Bio.Phylogeny.PhyBin
       ( driver, binthem, normalize, annotateWLabLists, unitTests, acquireTreeFiles,
         deAnnotate, retrieveHighlights, matchAnyHighlight )
       where

import qualified Data.Foldable as F
import           Data.Function       (on)
import           Data.List           (delete, minimumBy, sortBy, foldl1', foldl', intersperse, isPrefixOf)
import           Data.Maybe          (fromMaybe, catMaybes)
import           Data.Either         (partitionEithers)
import           Data.Time.Clock
import qualified Data.ByteString.Char8 as B
import qualified Data.Map                   as M
import qualified Data.Set                   as S
import qualified Data.Vector                 as V
import qualified Data.Vector.Unboxed         as U
import           Control.Monad       (forM, forM_, filterM, when, unless)
import qualified Control.Concurrent.Async as Async
import           Control.Exception   (evaluate)
import           Control.Applicative ((<$>),(<*>))
import           Control.Concurrent  (Chan)
import           System.FilePath     (combine)
import           System.Directory    (doesFileExist, doesDirectoryExist, createDirectoryIfMissing,
                                      getDirectoryContents, getCurrentDirectory)
import           System.IO           (openFile, hClose, IOMode(..), stdout)
import           System.Info         (os)
import           System.Process      (system)
import           System.Exit         (ExitCode(..))
import           Test.HUnit          ((~:),(~=?),Test,test)
import qualified Data.Clustering.Hierarchical as C

-- For vizualization:
import           Text.PrettyPrint.HughesPJClass hiding (char, Style)
import           Bio.Phylogeny.PhyBin.CoreTypes
import           Bio.Phylogeny.PhyBin.Parser (parseNewick, parseNewicks)
import           Bio.Phylogeny.PhyBin.PreProcessor (collapseBranchLenThresh, collapseBranchBootStrapThresh)
import           Bio.Phylogeny.PhyBin.Visualize (dotToPDF, dotNewickTree, viewNewickTree, dotDendrogram)
import           Bio.Phylogeny.PhyBin.RFDistance
import           Bio.Phylogeny.PhyBin.Binning
import           Bio.Phylogeny.PhyBin.Util

import Debug.Trace
----------------------------------------------------------------------------------------------------

-- Turn on for extra invariant checking:
debug :: Bool
debug = True

#ifdef SEQUENTIALIZE
#warning "SEQUENTIALIING execution.  Disabling all parallelism."
type Async t = t
async a = a
wait  x = return x 
#else
type Async t = Async.Async t
async = Async.async
wait  = Async.wait
#endif

-- | A dendrogram PLUS consensus trees at the intermediate points.
data DendroPlus a = DPLeaf (FullTree a)
                  | DPBranch !Double (FullTree ()) (DendroPlus a) (DendroPlus a)

----------------------------------------------------------------------------------------------------

-- | Driver to put all the pieces together (parse, normalize, bin)
driver :: PhyBinConfig -> IO ()
driver cfg@PBC{ verbose, num_taxa, name_hack, output_dir, inputs=files,
                do_graph, branch_collapse_thresh, bootstrap_collapse_thresh, highlights, 
                dist_thresh, clust_mode, use_hashrf, print_rfmatrix } =
   -- Unused: do_draw
 do 
    --------------------------------------------------------------------------------
    -- First, find out where we are and open the files:
    --------------------------------------------------------------------------------
    cd <- getCurrentDirectory 
    --putStrLn$ "PHYBIN RUNNING IN DIRECTORY: "++ cd

    bl <- doesDirectoryExist output_dir
    unless bl $ createDirectoryIfMissing True output_dir

    if isPrefixOf "mingw" os then
      -- TODO: Need a portable version of this.  'filemanip' would do:
      putStrLn$ "Not cleaning away previous phybin outputs (TODO: port this to Windows)."
     else do 
      putStrLn$ "Cleaning away previous phybin outputs..."
      system$ "rm -f "++output_dir++"/dendrogram.*"
      system$ "rm -f "++output_dir++"/cluster*"
      system$ "rm -f "++output_dir++"/distance_matrix.txt"
      system$ "rm -f "++output_dir++"/WARNINGS.txt"
      return ()

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
          
    highlightTrs <- retrieveHighlights name_hack labelTab highlights

    let allTaxaLs = M.elems labelTab
        totalTaxa = length allTaxaLs
    putStrLn$ "\nTotal unique taxa ("++ show (M.size labelTab) ++"):\n  "++ (unwords allTaxaLs)

    expected_num_taxa <-
      case num_taxa of
        Expected n | n == totalTaxa -> return n
                   | otherwise -> do putStrLn$ " !  Warning, was told to expect "++show n++
                                               " taxa, but there are "++show totalTaxa++" present in the dataset!"
                                     return n 
        _ -> do putStrLn$ "Note: --numtaxa not supplied, defaulting to expecting all "++show totalTaxa++" to be present..."
                return totalTaxa
    --------------------------------------------------------------------------------

    case bootstrap_collapse_thresh of
      Just thr -> putStrLn$" !+ Collapsing branches of bootstrap value less than "++show thr
      Nothing  -> return ()      

    case branch_collapse_thresh of 
      Just thr -> putStrLn$" !+ Collapsing branches of length less than "++show thr
      Nothing  -> return ()

    let do_one :: FullTree DefDecor -> IO (Int, [FullTree DefDecor], Maybe (Int, String))
        do_one (FullTree treename lblAcc parsed) = do 
           let
               -- Important: MUST collapse bootstrap first if we're doing both:
               pruned0 = case bootstrap_collapse_thresh of
                          Nothing  -> parsed
                          Just thr -> collapseBranchBootStrapThresh thr parsed
               pruned1 = case branch_collapse_thresh of 
                          Nothing  -> pruned0
                          Just thr -> collapseBranchLenThresh thr pruned0
               numL   = numLeaves pruned1
               
           -- TEMPTOGGLE
	   -- when False $ do putStrLn$ "DRAWING TREE"
           --                 viewNewickTree "Annotated"  (FullTree file lblAcc' annot)
           --                 viewNewickTree "Normalized" (FullTree file lblAcc' normal)
	   --      	   putStrLn$ "WEIGHTS OF NORMALIZED' CHILDREN: "++
           --                       show (map get_weight$ get_children normal)

           if numL /= expected_num_taxa
	    then do --putStrLn$ "\n WARNING: file contained an empty or single-node tree: "++ show file
 		    when verbose$ putStrLn$ "\n WARNING: tree contained unexpected number of leaves ("
					    ++ show numL ++"): "++ treename
		    return (0, [], Just (numL, treename))
	    else do when verbose$ putStr "."
                    return$ (numL, [FullTree treename lblAcc pruned1], Nothing)

    results <- mapM do_one fulltrees
    let (counts::[Int], validtreess, pairs:: [Maybe (Int, String)]) = unzip3 results
    let validtrees = concat validtreess
        warnings2  = catMaybes pairs
        
    putStrLn$ "\nNumber of input tree files: " ++ show num_files
    when (length warnings2 > 0) $
      putStrLn$ "Number of bad/unreadable input tree files: " ++ show (length warnings2)
    putStrLn$ "Number of VALID trees (correct # of leaves/taxa): " ++ show (length validtrees)
    putStrLn$ "Total tree nodes contained in valid trees: "++ show (sum counts)
    let all_branches = concatMap (F.foldr' (\x ls -> getBranchLen x:ls) []) validtrees    
    unless (null all_branches) $ do
      putStrLn$ "Average branch len over valid trees: "++ show (avg all_branches)
      putStrLn$ "Max/Min branch lengths: "++ show (foldl1' max all_branches,
                                                   foldl1' min all_branches)

    --------------------------------------------------------------------------------
    -- Next, dispatch on the mode and do the actual clustering or binning.
    --------------------------------------------------------------------------------

    (classes,binlist,asyncs) <- case clust_mode of
      BinThem         -> do x <- doBins validtrees
                            -- A list of bins, sorted by size:
                            let binlist = reverse $ sortBy (compare `on` fst) $
                                          map (\ (tr,OneCluster ls) -> (length ls, OneCluster ls)) $
                                          M.toList x
                            return (x,binlist,[])
      ClusterThem{linkage} -> do
        (mat, dendro) <- doCluster use_hashrf expected_num_taxa linkage validtrees        
        when print_rfmatrix $ printDistMat stdout mat
        hnd <- openFile  (combine output_dir ("distance_matrix.txt")) WriteMode
        printDistMat hnd mat
        hClose hnd
        --------------------
        writeFile (combine output_dir ("dendrogram.txt"))
                  (show$ fmap treename dendro)
        putStrLn " [finished] Wrote full dendrogram to file dendrogram.txt"
        let plotIt mnameMap = 
              if True -- do_graph
              then async (do
                putStrLn$ " [async] creating task to plot dendrogram.pdf"
                t0 <- getCurrentTime
                let dot = dotDendrogram cfg "dendrogram" 1.0 dendro mnameMap highlightTrs
                _ <- dotToPDF dot (combine output_dir "dendrogram.pdf") 
                t1 <- getCurrentTime          
                putStrLn$ " [finished] Writing dendrogram diagram ("
                          ++show(diffUTCTime t1 t0)++")")
              else async (return ())
        case dist_thresh of
          Nothing -> do a <- plotIt Nothing
                        return (M.empty,[],[a])
                     -- error "Fully hierarchical cluster output is not finished!  Use --editdist."
          Just dstThresh -> do            
            putStrLn$ "Combining all clusters at distance less than or equal to "++show dstThresh
            let clusts = sliceDendro (fromIntegral dstThresh) dendro
                classes = clustsToMap clusts
                -- Cache the lengths:
                wlens  = [ (length ls, OneCluster ls) | OneCluster ls <- clusts ]
                sorted0 = reverse$ sortBy (compare `on` fst) wlens -- TODO: Parallel sorting?
                nameMap = clustsToNameMap (map snd sorted0)
            gvizAsync <- plotIt (Just nameMap)
            putStrLn$ "After flattening, cluster sizes are: "++show (map fst sorted0)
            -- Flatten out the dendogram:
            return (classes, sorted0, [gvizAsync])

    unless (null binlist)$ reportClusts clust_mode binlist
        
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

--    putStrLn$ "Final number of tree bins: "++ show (M.size classes)

    unless (null warnings1 && null warnings2) $
	writeFile (combine output_dir "WARNINGS.txt")
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

    async2 <- case clust_mode of
                BinThem       -> outputBins binlist output_dir do_graph
                ClusterThem{} -> outputClusters expected_num_taxa binlist output_dir do_graph

    -- Wait on parallel tasks:
    putStrLn$ "Waiting for "++show (length$ async2:asyncs)++" asynchronous tasks to finish..."
    mapM_ wait (async2:asyncs)
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

doCluster :: Bool -> Int -> C.Linkage -> [FullTree a] -> IO (DistanceMatrix, C.Dendrogram (FullTree a))
doCluster use_hashrf expected_num_taxa linkage validtrees = do
  t0 <- getCurrentTime
  when use_hashrf$ putStrLn " Using HashRF-style algorithm..."
  let nwtrees  = map nwtree validtrees
      numtrees = length validtrees
      mat = if use_hashrf 
            then hashRF expected_num_taxa nwtrees 
            else fst (naiveDistMatrix nwtrees)
      ixtrees  = zip [0..] validtrees
      dist (i,t1) (j,t2) | j == i     = 0
--                         | i == numtrees-1 = 0 
                         | j < i      = fromIntegral ((mat V.! i) U.! j)
                         | otherwise  = fromIntegral ((mat V.! j) U.! i)
      dist1 a b = trace ("Taking distance between "++show (fst a, fst b)) $ dist a b
      dendro = fmap snd $ C.dendrogram linkage ixtrees dist
  -- Force the distance matrix:
  V.mapM_ evaluate mat
  t1 <- getCurrentTime
  putStrLn$ "Time to compute distance matrix: "++show(diffUTCTime t1 t0)
  putStrLn$ "Clustering using method "++show linkage
  return (mat,dendro)


reportClusts :: ClustMode -> [(Int, OneCluster StandardDecor)] -> IO ()
reportClusts mode binlist = do 
    let binsizes = map fst binlist

    putStrLn$ " Outcome: "++show (length binlist)++" clusters found, "++show (length$ takeWhile (>1) binsizes)
             ++" non-singleton, top bin sizes: "++show(take 10 binsizes)
    putStrLn$"  Up to first 30 bin sizes, excluding singletons:"
    forM_ (zip [1..30] binlist) $ \ (ind, (len, OneCluster ftrees)) -> do
       when (len > 1) $ -- Omit that long tail of single element classes...
          putStrLn$show$
           hcat [text ("  * cluster#"++show ind++", members "++ show len ++", "), 
                 case mode of
                   BinThem -> vcat [text ("avg bootstraps "++show (get_bootstraps$ avg_trees$ map nwtree ftrees)++", "),
                                    text "all: " <> pPrint (filter (not . null) $ map (get_bootstraps . nwtree) ftrees)]
                   ClusterThem{} -> hcat [] ]

-- | Convert a flat list of clusters into a map from individual trees to clusters.
clustsToMap :: [OneCluster StandardDecor] -> BinResults StandardDecor
clustsToMap clusts = F.foldl' fn M.empty clusts
  where
    fn acc theclust@(OneCluster ftrs) =
      F.foldl' (fn2 theclust) acc ftrs
    fn2 theclust acc (FullTree{nwtree}) =
      M.insert (anonymize_annotated nwtree) theclust acc

-- | Map each tree NAME onto the one-based index (in sorted order) of the cluster it
-- comes from.
clustsToNameMap :: [OneCluster StandardDecor] -> M.Map TreeName Int
clustsToNameMap clusts = foldl' fn M.empty (zip [1..] clusts)
  where
    fn acc (ix,(OneCluster ftrs)) =
      foldl' (\acc' t -> M.insert (treename t) ix acc') acc ftrs

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
   
-- | Turn a hierarchical clustering into a flat clustering.
sliceDendro :: Double -> (C.Dendrogram (FullTree DefDecor)) -> [OneCluster StandardDecor]
sliceDendro dstThresh den = loop den
  where 
   loop br@(C.Branch dist left right)
     -- Too far apart to combine:
     | dist > dstThresh  = loop left ++ loop right
     | otherwise         = [flattenDendro br]
   loop br@(C.Leaf _)    = [flattenDendro br]


--------------------------------------------------------------------------------

-- outputClusters :: (Num a1, Ord a1, Show a1) => Int -> [(a1, t, OneCluster a)] -> String -> Bool -> IO ()
outputClusters :: Int -> [(Int, OneCluster a)] -> String -> Bool -> IO (Async ())
outputClusters num_taxa binlist output_dir do_graph = do
    let numbins = length binlist
    let base i size = combine output_dir (filePrefix ++ show i ++"_"++ show size) 
    let consTrs = [ FullTree "consensus" (labelTable $ head ftrees) nwtr 
                  | (_, OneCluster ftrees) <- binlist
                  , let nwtr :: NewickTree StandardDecor
                        nwtr = annotateWLabLists $ fmap (const (Nothing,0)) $
                               consensusTree num_taxa $ map nwtree ftrees ]
    
    forM_ (zip3 [1::Int ..] binlist consTrs) $ \ (i, (size, OneCluster ftrees), fullConsTr) -> do
       writeFile (base i size ++".txt") (concat$ map ((++"\n") . treename) ftrees)
       writeFile (base i size ++"_consensus.tr") (show (displayStrippedTree fullConsTr) ++ "\n")
       writeFile (base i size ++"_alltrees.tr")
           (unlines [ show (displayStrippedTree ft) | ft <- ftrees ])

    putStrLn$ " [finished] Wrote contents of each cluster to cluster<N>_<size>.txt"
    putStrLn$ " [finished] Wrote representative (consensus) trees to cluster<N>_<size>_consensus.tr"
    if do_graph then do
      putStrLn$ "Next start the time consuming operation of writing out graphviz visualizations:"
      asyncs <- forM (zip3 [1::Int ..] binlist consTrs) $
                \ (i, (size, OneCluster membs), fullConsTr) -> async$ do
    	 when (size > 1 || numbins < 100) $ do
           let dot = dotNewickTree ("cluster #"++ show i) 1.0 fullConsTr
	   _ <- dotToPDF dot (base i size ++ ".pdf")
	   return ()
      async $ do
        mapM_ wait asyncs
        putStrLn$ " [finished] Wrote visual representations of consensus trees to "++filePrefix++"<N>_<size>.pdf"
     else do putStrLn$ "NOT creating processes to build per-cluster visualizations. (Not asked to.)"
             async (return ())

outputBins :: [(Int, OneCluster StandardDecor)] -> String -> Bool -> IO (Async ())
outputBins binlist output_dir  do_graph = do
    let numbins = length binlist
    let base i size = combine output_dir (filePrefix ++ show i ++"_"++ show size) 
    let avgs = map (avg_trees . map nwtree . clustMembers . snd) binlist
    forM_ (zip3 [1::Int ..] binlist avgs) $ \ (i, (size, OneCluster ftrees), avgTree) -> do
       let FullTree fstName labs _ = head ftrees
           fullAvgTr = FullTree fstName labs avgTree
       writeFile (base i size ++".txt") (concat$ map ((++"\n") . treename) ftrees)
       when debug$ do
         writeFile (base i size ++".dbg") (show (pPrint avgTree) ++ "\n")
       -- Printing the average tree instead of the stripped cannonical one:         
       writeFile   (base i size ++"_avg.tr")  (show (displayDefaultTree$ deAnnotate fullAvgTr) ++ "\n") -- FIXME

    putStrLn$ " [finished] Wrote contents of each bin to "++filePrefix++"<N>_<binsize>.txt"
    putStrLn$ "            Wrote representative trees to "++filePrefix++"<N>_<binsize>_avg.tr" 
    if do_graph then do
       putStrLn$ "Next start the time consuming operation of writing out graphviz visualizations:"
       asyncs <- forM (zip3 [1::Int ..] binlist avgs) $ \ (i, (size, OneCluster membs), avgTree) -> async $ do
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
       async $ do
         mapM_ wait asyncs
         putStrLn$ " [finished] Wrote visual representations of trees to "++filePrefix++"<N>_<size>.pdf"
     else async (return ())
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

avg ls = sum ls / fromIntegral (length ls)

-- | Parse extra trees in addition to the main inputs (for --highlight).
retrieveHighlights :: (String->String) -> LabelTable -> [FilePath] -> IO [[NewickTree ()]]
retrieveHighlights name_hack labelTab ls =
  mapM parseHighlight ls
  where
    parseHighlight file = do 
      bs <- B.readFile file
      let (lt2,htr) = parseNewick labelTab name_hack file bs
      unless (lt2 == labelTab) $
        error$"Tree given as --highlight includes taxa not present in main tree set: "++
              show(M.keys$ M.difference lt2 labelTab)            
      return (map (fmap (const())) htr)


-- | Create a predicate that tests trees for consistency with the set of --highlight
-- (consensus) trees.
--
-- Note, tree consistency is not the same as an exact match.  It's
-- like (<=) rather than (==).  All trees are consistent with the
-- "star topology".
matchAnyHighlight :: [[NewickTree ()]] -> NewickTree () -> Bool      
-- matchAnyHighlight :: [[NewickTree ()]] -> NewickTree () -> Maybe Int
-- If there is a match, return the index of the highlight that matched.
matchAnyHighlight highlightTrs =
   let matchers = map mkMatcher highlightTrs
       mkMatcher ls = let fns = map compatibleWith ls -- Multiple predicate functions
                      in \ tr -> -- Does it match any predicate?
                          or$ map (\f -> f tr) fns   
   in \ newtr -> 
       any (\f -> f newtr) matchers
