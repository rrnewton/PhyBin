{-# LANGUAGE RecordWildCards #-}
{-# OPTIONS_GHC -fwarn-unused-imports #-}

module Main where
import           Data.List (sort)
import qualified Data.ByteString.Lazy.Char8 as B
import qualified Data.Map as M
import qualified Data.Set as S
import           Control.Monad
import           Control.Concurrent    (Chan, readChan, ThreadId, forkIO)
import           System.Environment    (getArgs, withArgs)
import           System.Console.GetOpt (OptDescr(Option), ArgDescr(..), ArgOrder(..), usageInfo, getOpt)
import           System.Exit           (exitSuccess)
import           Test.HUnit            (runTestTT, Test, test)

import Data.GraphViz (runGraphvizCanvas,GraphvizCommand(Dot),GraphvizCanvas(Xlib))
import Bio.Phylogeny.PhyBin.CoreTypes 
        ( NewickTree(..), PhyBinConfig(..), default_phybin_config, DefDecor, StandardDecor(..),
          toLabel, fromLabel, Label, set_dec, map_labels ) 
import Bio.Phylogeny.PhyBin           (driver, binthem, normalize, annotateWLabLists, unitTests)
import Bio.Phylogeny.PhyBin.Parser    (parseNewick, unitTests)
import Bio.Phylogeny.PhyBin.Visualize (viewNewickTree, dotNewickTree_debug)

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

    | SelfTest

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
       --  TODO, FIXME: The "guessing" part doesn't actually work yet -- implement it!!
       --  What's a good algorithm?  Insist they all have the same number?  Take the mode?
       
     , Option ['n']     ["numtaxa"] (ReqArg (NumTaxa . read) "NUM") "expect NUM taxa for this dataset"

     , Option ['b']     ["branchcut"] (ReqArg (BranchThresh . read) "LEN") "collapse branches less than LEN"
       

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

usage :: String
usage = "\nUsage: phybin [OPTION...] files or directories...\n\n"++

        "PhyBin takes Newick tree files as input.  Paths of Newick files can\n"++
        "be passed directly on the command line.  Or, if directories are provided,\n"++
        "all files in those directories will be read.\n\n"++

        "Phybin output contains, at minimum, files of the form binXX_YY.tr,\n"++
        "each containing a list of input file paths that fall into that bin.\n\n"++

	"USAGE NOTES: Currently phybin ignores input trees with the wrong number of taxa.\n"++
	"If given a directory as input phybin will assume all contained files are Newick trees.\n\n"++

	"\nOptions include:\n"

defaultErr :: [String] -> t
defaultErr errs = error $ "ERROR!\n" ++ (concat errs ++ usageInfo usage options)


allUnitTests :: Test
-- allUnitTests = unitTests ++
allUnitTests = test 
  [ Bio.Phylogeny.PhyBin.unitTests
  , Bio.Phylogeny.PhyBin.Parser.unitTests
  ]
--  Bio.Phylogeny.PhyBin.Parser.unitTests

main :: IO ()
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

	   SelfTest -> do _ <- runTestTT allUnitTests; exitSuccess
                          
	   Output s -> return cfg { output_dir= s }

	   NumTaxa n -> return cfg { num_taxa= n }
     	   BranchThresh n -> return cfg { branch_collapse_thresh= Just n }
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
 	        (chan, _tr) <- viewNewickTree file $ 
			       annotateWLabLists$ 
			       map_labels name_hack $ 
			       parseNewick file str
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

----------------------------------------------------------------------------------------------------
-- TODO: expose a command line argument for testing.
-- The below test exposed my normalization bug relating to "deriving Ord".
-- I need to transform it into one or more proper unit tests.

main_test :: IO ()
main_test = 
 withArgs ["-w","~/newton_and_newton_local/datasets/yersinia/yersinia_trees/111.dnd","-m","../datasets/yersinia/name_table_hack_yersinia.txt"]
	  main 

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

--a_norm = NTInterior (0.13899,7,["BB","BJ","MB","ML","RE","SD","SM"]) [NTInterior (0.0,3,["RE","SD","SM"]) [NTLeaf (0.13143,1,["RE"]) "RE",NTInterior (5.697e-2,2,["SD","SM"]) [NTLeaf (5.977e-2,1,["SD"]) "SD",NTLeaf (3.95e-2,1,["SM"]) "SM"]],NTInterior (9.019e-2,2,["BB","BJ"]) [NTLeaf (0.11856,1,["BB"]) "BB",NTLeaf (0.13592,1,["BJ"]) "BJ"],NTInterior (0.13194,2,["MB","ML"]) [NTLeaf (0.19456,1,["MB"]) "MB",NTLeaf (0.16603,1,["ML"]) "ML"]]

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
