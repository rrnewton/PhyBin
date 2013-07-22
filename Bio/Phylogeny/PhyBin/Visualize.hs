{-# LANGUAGE NamedFieldPuns #-}

module Bio.Phylogeny.PhyBin.Visualize
       (dotNewickTree, dotToPDF, viewNewickTree,
        dotNewickTree_debug)
       where
import           Text.Printf        (printf)
import           Data.List          (elemIndex)
import           Data.Maybe         (fromJust)
import           Data.Map           ((!))
import           Data.Text.Lazy     (pack)
import           Control.Concurrent  (Chan, newChan, writeChan, forkIO)
import qualified Data.Graph.Inductive as G  hiding (run)
import qualified Data.GraphViz        as Gv hiding (parse, toLabel)
import qualified Data.GraphViz.Attributes.Complete as GA
import qualified Data.GraphViz.Attributes.Colors   as GC
-- import           Test.HUnit          ((~:),(~=?),Test,test)

import           Bio.Phylogeny.PhyBin.CoreTypes

----------------------------------------------------------------------------------------------------
-- Visualization with GraphViz and FGL:
----------------------------------------------------------------------------------------------------

-- First we need to be able to convert our trees to FGL graphs:
toGraph :: FullTree StandardDecor -> G.Gr String Double
toGraph (FullTree _ tbl tree) = G.run_ G.empty $ loop tree
  where
 fromLabel ix = tbl ! ix
 loop (NTLeaf _ name) = 
    do let str = fromLabel name
       _ <- G.insMapNodeM str
       return str
 loop (NTInterior (StandardDecor{sortedLabels}) ls) =
    do let bigname = concatMap fromLabel sortedLabels
       names <- mapM loop ls
       _ <- G.insMapNodeM bigname
       mapM_ (\x -> G.insMapEdgeM (bigname, x, 0.0)) names
       return bigname

-- This version uses the tree nodes themselves as graph labels.
toGraph2 :: FullTree StandardDecor -> G.Gr (NewickTree StandardDecor) Double       
toGraph2 (FullTree _ tbl tree) = G.run_ G.empty $ loop tree
  where
 loop node@(NTLeaf _  _) =  
    do _ <- G.insMapNodeM node 
       return ()
 loop node@(NTInterior _ ls) =
    do mapM_ loop ls
       _ <- G.insMapNodeM node
       -- Edge weights as just branchLen (not bootstrap):
       mapM_ (\x -> G.insMapEdgeM (node, x, branchLen$ get_dec x)) ls
       return ()

-- | Open a GUI window to displaya tree.
--
--   Fork a thread that then runs graphviz.
--   The channel retuned will carry a single message to signal
--   completion of the subprocess.
viewNewickTree :: String -> FullTree StandardDecor -> IO (Chan (), FullTree StandardDecor)
viewNewickTree title tree@(FullTree{nwtree}) =
  do chan <- newChan
     let dot = dotNewickTree title (1.0 / avg_branchlen [nwtree])
	                     tree
	 runit = do Gv.runGraphvizCanvas default_cmd dot Gv.Xlib
		    writeChan chan ()
     --str <- prettyPrint d
     --putStrLn$ "Generating the following graphviz tree:\n " ++ str
     _ <- forkIO runit
       --do runGraphvizCanvas Dot dot Xlib; return ()
       
     return (chan, tree)


--default_cmd = TwoPi -- Totally ignores edge lengths.
default_cmd :: Gv.GraphvizCommand
default_cmd = Gv.Neato

-- Show a float without scientific notation:
myShowFloat :: Double -> String
-- showFloat weight = showEFloat (Just 2) weight ""
myShowFloat fl = printf "%.4f" fl


dotToPDF :: Gv.DotGraph G.Node -> FilePath -> IO FilePath
dotToPDF dot file =
  Gv.runGraphvizCommand default_cmd dot Gv.Pdf file

-- | Convert a NewickTree to a graphviz Dot graph representation.
dotNewickTree :: String -> Double -> FullTree StandardDecor -> Gv.DotGraph G.Node
dotNewickTree title edge_scale atree@(FullTree _ tbl tree) = 
    --trace ("EDGE SCALE: " ++ show edge_scale) $
    Gv.graphToDot myparams graph
 where 
  graph = toGraph2 atree
  fromLabel ix = tbl ! ix  
  myparams :: Gv.GraphvizParams G.Node (NewickTree StandardDecor) Double () (NewickTree StandardDecor)
  myparams = Gv.defaultParams { Gv.globalAttributes= [Gv.GraphAttrs [GA.Label$ GA.StrLabel$ pack title]],
                                Gv.fmtNode= nodeAttrs, Gv.fmtEdge= edgeAttrs }
  nodeAttrs :: (Int,NewickTree StandardDecor) -> [GA.Attribute]
  nodeAttrs (_num, node) =
    let children = get_children node in 
    [ GA.Label$ GA.StrLabel$ pack$ 
      concatMap fromLabel $
      sortedLabels $ get_dec node
    , GA.Shape (if null children then {-PlainText-} GA.Ellipse else GA.PointShape)
    , GA.Style [GA.SItem GA.Filled []]
    ]

  -- TOGGLE:
  --  edgeAttrs (_,_,weight) = [ArrowHead noArrow, Len (weight * edge_scale + bump), GA.Label (StrLabel$ show (weight))]
  edgeAttrs (_,_,weight) = 
                           let draw_weight = compute_draw_weight weight edge_scale in
                           --trace ("EDGE WEIGHT "++ show weight ++ " drawn at "++ show draw_weight) $
			   [GA.ArrowHead Gv.noArrow,
                            GA.Label$ GA.StrLabel$ pack$ myShowFloat weight] ++ -- TEMPTOGGLE
			   --[ArrowHead noArrow, GA.Label (StrLabel$ show draw_weight)] ++ -- TEMPTOGGLE
			    if weight == 0.0
			    then [GA.Color [weighted$ GA.X11Color Gv.Red], GA.Len minlen]
			    else [GA.Len draw_weight]

  weighted c = GC.WC {GC.wColor=c, GC.weighting=Nothing}
  minlen = 0.7
  maxlen = 3.0
  compute_draw_weight w scale = 
    let scaled = (abs w) * scale + minlen in 
    -- Don't draw them too big or it gets annoying:
    (min scaled maxlen)


-- | This version shows the ordered/rooted structure of the normalized tree.
dotNewickTree_debug :: String -> FullTree StandardDecor -> Gv.DotGraph G.Node
dotNewickTree_debug title atree@(FullTree _ tbl tree) = Gv.graphToDot myparams graph
 where 
  graph = toGraph2 atree
  fromLabel ix = tbl ! ix    
  myparams :: Gv.GraphvizParams G.Node (NewickTree StandardDecor) Double () (NewickTree StandardDecor)
  myparams = Gv.defaultParams { Gv.globalAttributes= [Gv.GraphAttrs [GA.Label$ GA.StrLabel$ pack title]],
			        Gv.fmtNode= nodeAttrs, Gv.fmtEdge= edgeAttrs }
  nodeAttrs :: (Int,(NewickTree StandardDecor)) -> [GA.Attribute]
  nodeAttrs (num,node) =
    let children = get_children node in 
    [ GA.Label (if null children 
  	        then GA.StrLabel$ pack$ concatMap fromLabel $ sortedLabels $ get_dec node
	        else GA.RecordLabel$ take (length children) $ 
                                  -- This will leave interior nodes unlabeled:
	                          map (GA.PortName . GA.PN . pack) $ map show [1..]
		                  -- This version gives some kind of name to interior nodes:
--	                          map (\ (i,ls) -> LabelledTarget (PN$ show i) (fromLabel$ head ls)) $ 
--                                       zip [1..] (map (thd3 . get_dec) children)
               )
    , GA.Shape GA.Record
    , GA.Style [GA.SItem GA.Filled []]
    ]

  edgeAttrs (num1, num2, _weight) = 
    let node1 = fromJust$ G.lab graph num1 
	node2 = fromJust$ G.lab graph num2 	
	ind = fromJust$ elemIndex node2 (get_children node1)
    in [GA.TailPort$ GA.LabelledPort (GA.PN$ pack$ show$ 1+ind) (Just GA.South)]




--------------------------------------------------------------------------------
-- Unit Testing
--------------------------------------------------------------------------------


-- tt0 :: IO (Chan (), FullTree)
-- tt0 = drawNewickTree "tt0" $ annotateWLabLists $ parseNewick "" "(A,(C,D,E),B);"

-- tt3 :: IO (Chan (), FullTree)
-- tt3 = drawNewickTree "tt3" tt

-- tt4 :: IO (Chan (), FullTree)
-- tt4 = drawNewickTree "tt4"$ trace ("FINAL: "++ show (pPrint norm4)) $ norm4

-- tt5 :: IO (Chan (), FullTree)
-- tt5 = drawNewickTree "tt5"$ norm5


-- tt5' :: String
-- tt5' = prettyPrint' $ dotNewickTree "norm5" 1.0 norm5

-- ttall :: IO (Chan (), FullTree)
-- ttall = do tt3; tt4; tt5

-- tt2 :: G.Gr String Double
-- tt2 = toGraph tt


-- unitTests :: Test
-- unitTests = test
--    [ 
--    ]


-- TEMP / HACK:
prettyPrint' :: Show a => a -> String
prettyPrint' = show
