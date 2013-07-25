{-# LANGUAGE RecordWildCards, TupleSections, NamedFieldPuns #-}
{-# OPTIONS_GHC -fwarn-unused-imports -fwarn-incomplete-patterns #-}

module Main where
import           Data.List (sort, intersperse, foldl')
import qualified Data.ByteString.Lazy.Char8 as B
import qualified Data.Map as M
import qualified Data.Set as S
import qualified Data.IntSet as IS
import           Control.Monad
import           Control.Concurrent    (Chan, readChan, ThreadId, forkIO)
import           System.Environment    (getArgs, withArgs)
import           System.Console.GetOpt (OptDescr(Option), ArgDescr(..), ArgOrder(..), usageInfo, getOpt)
import           System.Exit           (exitSuccess)
import           System.IO             (stdout) 
import           Test.HUnit            (runTestTT, Test, test, (~:))

import Control.Applicative ((<$>))
import qualified Data.Vector                 as V
import           Test.HUnit                  as HU

import Data.GraphViz (runGraphvizCanvas,GraphvizCommand(Dot),GraphvizCanvas(Xlib))
import Text.PrettyPrint.HughesPJClass hiding (char, Style)

import Bio.Phylogeny.PhyBin.CoreTypes          
import Bio.Phylogeny.PhyBin           (driver, binthem, normalize, annotateWLabLists,
                                       unitTests, acquireTreeFiles, deAnnotate)
import Bio.Phylogeny.PhyBin.Parser    (parseNewick, parseNewicks, parseNewickFiles, unitTests)
import Bio.Phylogeny.PhyBin.Visualize (viewNewickTree, dotNewickTree_debug)
import Bio.Phylogeny.PhyBin.RFDistance 
import Bio.Phylogeny.PhyBin.PreProcessor

import qualified Data.Clustering.Hierarchical as C


import Version

----------------------------------------------------------------------------------------------------
-- MAIN script: Read command line options and call the program.
----------------------------------------------------------------------------------------------------

-- Note: ORDER is important here, we process options in this order:
data Flag 
    = Verbose  
    | Version 
    | Output String
    | NumTaxa Int
    | BranchThresh Double      
    | NullOpt
    | Graph | Draw
    | Force 
    | View
    | TabDelimited Int Int

    | HashRF Bool
    | SelfTest
    | RFMatrix | LineSetDiffMode | PrintNorms | PrintReg
    | Cluster C.Linkage
    | BinningMode
    | EditDistThresh Int
    | DendogramOnly

    | NameCutoff String
    | NamePrefix Int
    | NameTable String  -- Must come after Cutoff/Prefix

  deriving (Show, Eq, Ord) -- ORD is really used.


parseTabDelim :: String -> Flag
parseTabDelim _str = 
  TabDelimited 8 9
    
options :: [OptDescr Flag]
options =
     [ Option ['v']     ["verbose"] (NoArg Verbose)    "print WARNINGS and other information (recommended at first)"
     , Option ['V']     ["version"] (NoArg Version)    "show version number"

     , Option ['o']     ["output"]  (ReqArg Output "DIR")  "set directory to contain all output files (default \"./phybin_out/\")"
     , Option []     ["selftest"]   (NoArg SelfTest)   "run internal unit tests"
       
{- -- TODO: FIXME: IMPLEMENT THIS:
     , Option []        []          (NoArg NullOpt)  ""
     , Option ['t']     ["tabbed"]  (ReqArg parseTabDelim "NUM1:NUM2")$  "assume the input is a tab-delimited file with gene names \n"++
		                                                        "in column NUM1 and Newick trees in NUM2"
-}       
     , Option []        []          (NoArg NullOpt)  ""
     , Option []        []  (NoArg$ error "internal problem")  "----------------------------- Clustering Options ------------------------------"

     , Option []    ["bin"]      (NoArg BinningMode)$  "Use simple binning, the cheapest form of 'clustering'"
     , Option []    ["single"]   (NoArg$ Cluster C.SingleLinkage)  $  "Use single-linkage clustering (nearest neighbor)"
     , Option []    ["complete"] (NoArg$ Cluster C.CompleteLinkage)$  "Use complete-linkage clustering (furthest neighbor)"
     , Option []    ["UPGMA"]    (NoArg$ Cluster C.UPGMA)          $  "Use Unweighted Pair Group Method (average linkage)"

     , Option []    ["editdist"]  (ReqArg (EditDistThresh . read) "DIST")$
                                  "Combine all clusters separated by DIST or less.  Report a flat list of clusters."                                  
     , Option []    ["dendogram"] (NoArg DendogramOnly)$ "Report a hierarchical clustering (default)"

     , Option []        []     (NoArg$ error "internal problem")  "  Select Robinson-Foulds (symmetric difference) distance algorithm:"
     , Option []    ["simple"] (NoArg$ HashRF False) "(default) use direct all-to-all comparison"
     , Option []    ["hashrf"] (NoArg$ HashRF True)    "use a variant of the HashRF algorithm for the distance matrix"
       
     , Option []        []          (NoArg NullOpt)  ""
     , Option []        []  (NoArg$ error "internal problem")  "----------------------------- Visualization --------------------------------"
     , Option ['g']     ["graphbins"] (NoArg Graph)  "use graphviz to produce .dot and .pdf output files"
-- TODO: Produce the consensus tree as well as the individual trees.
     , Option ['d']     ["drawbins"]  (NoArg Draw)   "like -g, but open GUI windows to show each bin's tree"

     , Option ['w']     ["view"]    (NoArg View)$  "for convenience, \"view mode\" simply displays input Newick files without binning" 

     , Option []        []          (NoArg NullOpt)  ""
     , Option []        []  (NoArg$ error "internal problem")  "---------------------------- Tree pre-processing -----------------------------"

     , Option ['n']     ["numtaxa"] (ReqArg (NumTaxa . read) "NUM") "expect NUM taxa for this dataset"

     , Option ['b']     ["branchcut"] (ReqArg (BranchThresh . read) "LEN") "collapse branches less than LEN"
              
     , Option []        []          (NoArg NullOpt)  ""
     , Option []        []  (NoArg$ error "internal problem")  "--------------------------- Extracting taxa names ----------------------------"
--     , Option ['n']     ["numtaxa"] (ReqArg (NumTaxa . read) "NUM") "expect NUM taxa for this dataset (otherwise it will guess)"
       --  TODO, FIXME: The "guessing" part doesn't actually work yet -- implement it!!
       --  What's a good algorithm?  Insist they all have the same number?  Take the mode?
       
{- -- TODO: FIXME: IMPLEMENT THIS:
     , Option ['f']     ["force"]   (NoArg Force)    "force phybin to consume and bin trees with different numbers of taxa"
-}
     , Option []        []          (NoArg NullOpt)  ""
     , Option ['p']     ["nameprefix"]  (ReqArg (NamePrefix . read) "NUM") $ 
		  "Leaf names in the input Newick trees can be gene names, not taxa.\n"++
    	    	  "Then it is typical to extract taxa names from genes.  This option extracts\n"++
                  "a prefix of NUM characters to serve as the taxa name."

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
                  
     , Option []        []          (NoArg NullOpt)  ""
     , Option []        []  (NoArg$ error "internal problem")  "--------------------------- Utility Modes ----------------------------"
     , Option [] ["rfdist"]  (NoArg RFMatrix)        "print a Robinson Foulds distance matrix for the input trees"
     , Option [] ["setdiff"] (NoArg LineSetDiffMode) "for convenience, print the set difference between cluster*.txt files"
     , Option [] ["print"]      (NoArg PrintReg)     "simply print out a concise form of each input tree"       
     , Option [] ["printnorms"] (NoArg PrintNorms)   "simply print out a concise and NORMALIZED form of each input tree"
     ]

usage :: String
usage = "\nUsage: phybin [OPTION...] files or directories...\n\n"++

--        "MODE must be one of 'bin' or 'cluster'.\n\n" ++

        "PhyBin takes Newick tree files as input.  Paths of Newick files can\n"++
        "be passed directly on the command line.  Or, if directories are provided,\n"++
        "all files in those directories will be read.  Taxa are named based on the files\n"++
        "containing them.  If a file contains multiple trees, all are read by phybin, and\n"++  
        "the taxa name then includes a suffix indicating the position in the file:\n"++
        " e.g. FILENAME_0, FILENAME_1, etc.\n"++
        "\n"++          

        -- "In binning mode, Phybin output contains, at minimum, files of the form binXX_YY.tr,\n"++
        -- "each containing a list of input file paths that fall into that bin.  XX signifies\n"++
        -- "the rank of the bin (biggest first), and YY is how many trees it contains.\n"++
        -- "\n"++        

        "When clustering trees, Phybin computes a complete all-to-all Robinson-Foulds distance matrix.\n"++
        "If a threshold distance (tree edit distance) is given, then a flat set of clusters\n"++
        "will be produced in files clusterXX_YY.tr.  Otherwise it produces a full dendogram (UNFINISHED).\n"++
        "\n"++  

        "Binning mode provides an especially quick-and-dirty form of clustering.\n"++
        "When running with the --bin option, only exactly equal trees are put in the same cluster.\n"++
        "Tree pre-processing still applies, however: for example collapsing short branches.\n"++
        "\n"++        

	"USAGE NOTES: \n"++
        " * Currently phybin ignores input trees with the wrong number of taxa.\n"++
	" * If given a directory as input phybin will assume all contained files are Newick trees.\n\n"++

	"\nOptions include:\n"

defaultErr :: [String] -> t
defaultErr errs = error $ "ERROR!\n" ++ (concat errs ++ usageInfo usage options)

--------------------------------------------------------------------------------

main :: IO ()
main = 
  do argv <- getArgs 

     (opts,files) <- 
       case getOpt Permute options argv of
	 (o,n,[]  ) -> return (o,n)
         (_,_,errs) -> defaultErr errs

     all_inputs <- acquireTreeFiles files
     let process_opt cfg opt = case opt of 
	   NullOpt -> return cfg
	   Verbose -> return cfg { verbose= True } 
	   Version -> do putStrLn$ "phybin version "++phybin_version; exitSuccess

	   SelfTest -> do _ <- runTestTT allUnitTests; exitSuccess

     	   RFMatrix -> return cfg { print_rfmatrix= True }

           LineSetDiffMode -> do
             bss <- mapM B.readFile files
             case map (S.fromList . B.lines) bss of
               [set1,set2] -> do let [f1,f2] = files
                                 let diff = S.difference set1 set2
                                 putStrLn$" Found "++show(S.size diff)++" lines occuring in "++f1++" but not "++f2
                                 mapM_ B.putStrLn $ S.toList diff
               oth -> error $"Line set difference mode expects two files as input, got "++show(length oth)
             exitSuccess

           PrintNorms -> return cfg
           PrintReg   -> return cfg
           
           Cluster lnk -> return cfg { clust_mode = ClusterThem lnk }
           HashRF  bl  -> return cfg { use_hashrf = bl }
           BinningMode -> return cfg { clust_mode = BinThem }
           EditDistThresh n -> return cfg { dist_thresh = Just n }
           DendogramOnly    -> return cfg { dist_thresh = Nothing }
     
	   Output s -> return cfg { output_dir= s }

	   NumTaxa n -> return cfg { num_taxa= n }
     	   BranchThresh n -> return cfg { branch_collapse_thresh= Just n }
	   Graph     -> return cfg { do_graph= True } 
	   Draw	     -> return cfg { do_draw = True } 
	   View      -> return cfg -- Handled below

	   TabDelimited _ _ -> error "tabbed option not handled yet"
	   Force            -> error "force option not handled yet"


	   NameCutoff str -> let chopper = takeWhile (not . flip S.member (S.fromList str)) 
			     in return cfg { name_hack = chopper . name_hack cfg }
	   NamePrefix n   -> return cfg { name_hack = (take n) . name_hack cfg }

           -- This should always be after cutoff/prefix:
	   NameTable file -> do reader <- name_table_reader file
				return cfg { name_hack = reader . name_hack cfg }

     config <- foldM process_opt default_phybin_config{ inputs= all_inputs } 
	             (sort opts) -- NOTE: options processed in sorted order.

     when (null files) $ do
	defaultErr ["No file arguments!\n"]

     ------------------------------------------------------------
     -- This mode kicks in AFTER config options are processed.
     when (elem PrintReg opts) $ do 
       (_,fts) <- parseNewickFiles (name_hack config) all_inputs
       forM_ fts $ \ ft@(FullTree name _ _) -> do
         putStrLn $ "Tree "++show name
         putStrLn$ show$ displayDefaultTree ft
       exitSuccess
     ------------------------------------------------------------
     when (elem PrintNorms opts) $ do 
       (_,fts) <- parseNewickFiles (name_hack config) all_inputs
       forM_ fts $ \ ft@(FullTree name _ _) -> do
         putStrLn $ "Tree , NORMALIZED:"++show name
         putStrLn$ show$ displayDefaultTree$ deAnnotate $
           liftFT (normalize . annotateWLabLists) ft
       exitSuccess
     ------------------------------------------------------------
     when (View `elem` opts) $ do 
       view_graphs config
       exitSuccess
     ------------------------------------------------------------
     -- Otherwise do the main, normal thing:
     driver config

view_graphs :: PhyBinConfig -> IO ()
view_graphs PBC{..} = 
           do chans <- forM inputs $ \ file -> do 
                putStrLn$ "Drawing "++ file ++"...\n"
		str <- B.readFile file
		putStrLn$ "Parsed: " ++ (B.unpack str)
                let (tbl,tr) = parseNewick M.empty name_hack file str
 	        (chan, _tr) <- viewNewickTree file (FullTree "" tbl (annotateWLabLists (head tr)))
	        return chan
	      forM_ chans readChan 
	      return ()


--------------------------------------------------------------------------------
-- Every dataset it seems needs a new hack on the names!

name_table_reader :: String -> IO (String -> String)
name_table_reader file = 
  do contents <- readFile file
     let mp = M.fromList $ 
	      map (\ls -> case ls of 
		           [a,b] -> (a,b)
		           _ -> error$ "Each line of "++file++" must contain two whitespace free strings: "++ unwords ls) $ 
	      filter (not . null) $
	      map tokenize $ 
	      lines contents
     return$ 
       \ name_to_hack -> 

       case M.lookup name_to_hack mp of -- Could use a trie
	     Just x -> x
	     Nothing -> name_to_hack
  where
    tokenize :: String -> [String]
    tokenize line =
      case words line of
        []        -> []
        [one]     -> error$"Name table contained bad line:  "++ show one
        [one,two] -> [one,two]
        (one:rest) -> [one, unwords rest]

temp :: IO ()
temp = driver default_phybin_config{ num_taxa=7, inputs=["../datasets/test.tr"] }


--------------------------------------------------------------------------------
-- Aggregated Unit Tests
--------------------------------------------------------------------------------

allUnitTests :: Test
-- allUnitTests = unitTests ++
allUnitTests = test 
  [ Bio.Phylogeny.PhyBin.unitTests
  , Bio.Phylogeny.PhyBin.Parser.unitTests
  , "norm/Bip1" ~: (testNorm prob1)
  , "bipTreeConversion" ~: testBipConversion
  , "t3" ~: t3_consensusTest, "t4" ~: t4_consensusTest, "t5" ~: t5_consensusTest
  ]

-- [2013.07.23]      
-- This was INCORRECTLY normalizing to:
--     ((1_, 2_), (7_, (18, 6_)), ((14, 3_), (19, (13, 5_))))
prob1 :: String
prob1 = "(5_, (19, ((3_, 14), ((2_, 1_), (7_, (6_, 18))))), 13);"

-- | Make sure that the normalized version of a tree yields the same bipartitions as
-- the unnormalized one.
testNorm :: String -> IO ()
testNorm str = do
  let (labs,parseds) = parseNewick M.empty id "test" (B.pack str)
      parsed = head parseds
      normed = normalize $ annotateWLabLists parsed
      bips1  = allBips parsed
      bips2  = allBips normed
      added   = S.difference bips2 bips1
      removed = S.difference bips1 bips2
      -- dispBips bip = show$
      --   map (map (labs M.!)) $ 
      --   map IS.toList$ S.toList bip
  unless (bips1 == bips2) $ do
    putStrLn$ "Normalized this: "++show (displayDefaultTree $ FullTree "" labs parsed)
    putStrLn$ "To this        : "++show (displayDefaultTree $ deAnnotate $ FullTree "" labs normed)
    error$ "Normalization added and removed these bipartitions, respectively:\n  "
           ++concat (intersperse " " (map (dispBip labs) (S.toList added))) ++"\n  "
           ++concat (intersperse " " (map (dispBip labs) (S.toList removed)))


-- 112 and 13
rftest :: IO ()
rftest = do 
  (mp,[t1,t2]) <- parseNewickFiles (take 2) ["tests/13.tr", "tests/112.tr"]
  putStrLn$ "Tree 13           : " ++ show (displayDefaultTree t1)
  putStrLn$ "Tree 112          : "++ show (displayDefaultTree t2)

  putStrLn$ "Tree 13 normed    : "++ show (disp t1)
  putStrLn$ "Tree 112 normed   : "++ show (disp t2)

  putStrLn$ "13  collapsed 0.02: " ++show (disp$ liftFT (collapseBranchLenThresh 0.02) t1)
  putStrLn$ "112 collapsed 0.02: " ++show (disp$ liftFT (collapseBranchLenThresh 0.02) t2)

  putStrLn$ "13  collapsed 0.03: " ++show (disp$ liftFT (collapseBranchLenThresh 0.03) t1)
  putStrLn$ "112 collapsed 0.03: " ++show (disp$ liftFT (collapseBranchLenThresh 0.03) t2)  

  let (mat,_) = distanceMatrix [nwtree t1, nwtree t2]
  printDistMat stdout mat
  return ()
 where
  disp (FullTree nm labs tr) =
    let collapsed :: AnnotatedTree 
        collapsed = normalize$ annotateWLabLists tr
    in displayDefaultTree$ deAnnotate $  FullTree nm labs collapsed

-- | This test was done with --editdist 4 --complete
t3_consensusTest :: IO ()
t3_consensusTest = consensusTest "./tests/t3_consensus/cluster1_284_alltrees.tr"
                                 "./tests/t3_consensus/cluster1_284_consensus.tr"

-- | This test was done with --editdist 0 --complete
t4_consensusTest :: IO ()
t4_consensusTest = consensusTest "./tests/t4_consensus/cluster1_16_alltrees.tr"
                                 "./tests/t4_consensus/cluster1_16_consensus.tr"

-- | This test was done with --editdist 1 --complete
t5_consensusTest :: IO ()
t5_consensusTest = consensusTest "./tests/t5_consensus/cluster1_35_alltrees.tr"
                                 "./tests/t5_consensus/cluster1_35_consensus.tr"

consensusTest :: String -> String -> IO ()
consensusTest alltrees consensus = do  
  (_,ctree:ftrees)  <- parseNewickFiles id [consensus,alltrees]
  let num_taxa      = numLeaves (nwtree ctree)
      plainTrs      = map nwtree  ftrees 
      eachbips      = map allBips plainTrs
      totalBips     = foldl' S.union        S.empty eachbips
      intersectBips = foldl' S.intersection S.empty eachbips
      FullTree _ labs _ = ctree
      linesPrnt x = unlines (map (("  "++) . dispBip labs) $ S.toList x)
  putStrLn$ "Bips in each: "++     show (map S.size eachbips)
  putStrLn$ "Total bips in all: "++show (S.size totalBips)  
  putStrLn$ "Bips in common: "++   show (S.size intersectBips)
  putStrLn$ "Bips of first member:\n" ++ linesPrnt (head eachbips)
  putStrLn$ "Some bips in the union that are NOT in the intersection:\n" ++
     linesPrnt (S.fromList$ take 20$ S.toList$ S.difference totalBips intersectBips)
  let cbips = allBips $ nwtree ctree
  putStrLn$ "ConsensusBips ("++show (S.size cbips)++"):\n"++linesPrnt cbips
  putStrLn$"Things in the consensus that should NOT be:\n"++linesPrnt (S.difference cbips intersectBips)
  putStrLn$"Things not in the consensus that SHOULD be:\n"++linesPrnt (S.difference intersectBips cbips)

  putStrLn$ "Now recomputing consensus tree for "++show num_taxa++" taxa"
  let ctree2 = consensusTree num_taxa plainTrs
      cbips2 = allBips ctree2
  putStrLn$ "Freshly recomputed consensusBips ("++show (S.size cbips2)++"):\n"++linesPrnt cbips2
  HU.assertEqual "Consensus tree on disk should match computed one:"
         cbips cbips2 -- (allBips$ fmap (const ()) $ nwtree ctree)         
  
  putStrLn " Partial distance matrix WITHIN this cluster:"
  let (mat,_) = distanceMatrix (map nwtree ftrees)
  printDistMat stdout (V.take 30 mat)
  HU.assertBool "Consensus should only include bicbips2ps in the members" (S.isSubsetOf cbips totalBips)
  HU.assertEqual "Consensus tree matches intersected bips" cbips intersectBips 
  return ()

testBipConversion :: IO ()
testBipConversion = 
  do (_,trs) <- allTestTrees
     mapM_ testOne trs
     putStrLn "All trees passed bip conversion."
  where
    testOne (FullTree{nwtree}) = do
      let sz    = numLeaves nwtree   
          bips  = allBips nwtree
          tr2   = bipsToTree sz bips
          bips2 = allBips tr2
      assertEqual "round trip bips->tree->bips" bips bips2

-- | Read in all test trees which we happen to have put in the repository for testing
-- purposes.
allTestTrees :: IO (LabelTable, [FullTree DefDecor])
allTestTrees =
  parseNewickFiles id $
  [ "./tests/t3_consensus/cluster1_284_alltrees.tr"
  , "./tests/t3_consensus/cluster1_284_consensus.tr"
  , "./tests/t4_consensus/cluster1_16_alltrees.tr"
  , "./tests/t4_consensus/cluster1_16_consensus.tr"
  , "./tests/t5_consensus/cluster1_35_alltrees.tr"
  , "./tests/t5_consensus/cluster1_35_consensus.tr"
  ]
  

----------------------------------------------------------------------------------------------------
-- TODO: expose a command line argument for testing.
-- The below test exposed my normalization bug relating to "deriving Ord".
-- I need to transform it into one or more proper unit tests.

main_test :: IO ()
main_test = 
 withArgs ["-w","~/newton_and_newton_local/datasets/yersinia/yersinia_trees/111.dnd","-m","../datasets/yersinia/name_table_hack_yersinia.txt"]
	  main 

-- [2013.07.22] Disabling for the new Label representation:
{-
pa :: NewickTree DefDecor
pa = set_dec (Nothing,1) $ 
    NTInterior () [NTInterior () [NTLeaf () "RE",NTInterior () [NTLeaf () "SD",NTLeaf () "SM"]],NTInterior () [NTLeaf () "BB",NTLeaf () "BJ"],NTInterior () [NTLeaf () "MB",NTLeaf () "ML"]]

pb :: NewickTree DefDecor
pb = set_dec (Nothing,1) $ 
    NTInterior () [NTInterior () [NTLeaf () "BB",NTLeaf () "BJ"],NTInterior () [NTLeaf () "MB",NTLeaf () "ML"],NTInterior () [NTLeaf () "RE",NTInterior () [NTLeaf () "SD",NTLeaf () "SM"]]]

ls1 :: [(String, NewickTree DefDecor)]
ls1 = [("a",pa),("b",pb)]

-- This is one:
num_binned :: Int
num_binned = M.size $ binthem ls1

a_ :: (String, NewickTree DefDecor)
a_ = ("980.dnd",
      fmap (\x -> (Nothing,x)) $ 
      NTInterior 0.0 [NTInterior 5.697e-2 [NTLeaf 3.95e-2 "SM",NTLeaf 5.977e-2 "SD"],NTLeaf 0.13143 "RE",NTInterior 0.13899 [NTInterior 9.019e-2 [NTLeaf 0.11856 "BB",NTLeaf 0.13592 "BJ"],NTInterior 0.13194 [NTLeaf 0.19456 "MB",NTLeaf 0.16603 "ML"]]])

b_ :: (String, NewickTree DefDecor)
b_ = ("999.dnd",
      fmap (\x -> (Nothing,x)) $ 
      NTInterior 0.0 [NTInterior 6.527e-2 [NTInterior 0.13734 [NTLeaf 2.975e-2 "SM",NTLeaf 3.002e-2 "SD"],NTLeaf 0.18443 "RE"],NTInterior 6.621e-2 [NTLeaf 0.16184 "MB",NTLeaf 0.15233 "ML"],NTInterior 0.23143 [NTLeaf 9.192e-2 "BB",NTLeaf 0.10125 "BJ"]])

-- But THIS is two:  ack!
num2 :: Int
num2 = M.size $ binthem [a_,b_]

-- Here's the test that breaks things:
a_norm :: NewickTree StandardDecor
a_norm = normalize (annotateWLabLists$ snd a_)

b_norm_ :: NewickTree (Double, Int, [Label])
b_norm_ = NTInterior (0.0,7,["BB","BJ","MB","ML","RE","SD","SM"])
         [NTInterior (0.23143,2,["BB","BJ"]) [NTLeaf (9.192e-2,1,["BB"]) "BB",NTLeaf (0.10125,1,["BJ"]) "BJ"],NTInterior (6.621e-2,2,["MB","ML"]) [NTLeaf (0.16184,1,["MB"]) "MB",NTLeaf (0.15233,1,["ML"]) "ML"],NTInterior (6.527e-2,3,["RE","SD","SM"]) [NTLeaf (0.18443,1,["RE"]) "RE",NTInterior (0.13734,2,["SD","SM"]) [NTLeaf (3.002e-2,1,["SD"]) "SD",NTLeaf (2.975e-2,1,["SM"]) "SM"]]]

b_norm :: NewickTree StandardDecor
b_norm = fmap (\ (bl,w,ls) -> StandardDecor bl Nothing w ls) b_norm_

d1 :: IO (Chan (), NewickTree StandardDecor)
d1 = viewNewickTree "" a_norm

d2 :: IO (Chan (), NewickTree StandardDecor)
d2 = viewNewickTree "" b_norm

d1_ :: IO ThreadId
d1_ = forkIO $ do runGraphvizCanvas Dot (dotNewickTree_debug "" a_norm) Xlib; return ()
                  
d2_ :: IO ThreadId
d2_ = forkIO $ do runGraphvizCanvas Dot (dotNewickTree_debug "" b_norm) Xlib; return ()


-- | A of a tree with _____ weights attached to it:
withBootstrap :: String
withBootstrap = "((((A8F330_:0.01131438136322714984,(G0GWK2_:0.00568050636963043226,(Q92FV4_:0.00284163304504484121,((B0BVQ5_:0.00319487112504297311,A8GU65_:0.00000122123005994819)74:0.00279881991324161267,(C3PM27_:0.00560787769333294297,C4K2Z0_:0.00559642713265556899)15:0.00000122123005994819)4:0.00000122123005994819)56:0.00276851661606284868)60:0.00283144414216590342)76:0.00886304965525876697,(A8GQC0_:0.05449879836105625541,(A8F0B2_:0.04736199885985507840,Q4UJN9_:0.02648399728559588939)64:0.00905997055810744446)28:0.00323255855543533657)29:0.02237505187863457132,(Q1RGK5_:0.00000122123005994819,A8GYD7_:0.00000122123005994819)100:0.28299884298270094884)100:0.05776841634437222123,(Q9ZC84_:0.00000122123005994819,D5AYH5_:0.00000122123005994819)99:0.00951976341375833368,Q68VM9_:0.04408933524904214141);"

withBootstrap2 :: String
withBootstrap2 = "((((A8F330_:0.01131438136322714984,(G0GWK2_:0.00568050636963043226,(Q92FV4_:0.00284163304504484121,((B0BVQ5_:0.00319487112504297311,A8GU65_:0.00000122123005994819):0.00279881991324161267[74],(C3PM27_:0.00560787769333294297,C4K2Z0_:0.00559642713265556899):0.00000122123005994819[15]):0.00000122123005994819[4]):0.00276851661606284868[56]):0.00283144414216590342[60]):0.00886304965525876697[76],(A8GQC0_:0.05449879836105625541,(A8F0B2_:0.04736199885985507840,Q4UJN9_:0.02648399728559588939):0.00905997055810744446[64]):0.00323255855543533657[28]):0.02237505187863457132[29],(Q1RGK5_:0.00000122123005994819,A8GYD7_:0.00000122123005994819):0.28299884298270094884[100]):0.05776841634437222123[100],(Q9ZC84_:0.00000122123005994819,D5AYH5_:0.00000122123005994819):0.00951976341375833368[99],Q68VM9_:0.04408933524904214141);"
-}
