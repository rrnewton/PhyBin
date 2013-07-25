{-# LANGUAGE NamedFieldPuns #-}

module Bio.Phylogeny.PhyBin.Visualize
       (dotNewickTree, dotToPDF, viewNewickTree,
        dotNewickTree_debug,

        -- * Dendrogram visualization
        dendrogramToGraph, dotDendrogram
       )
       where
import           Text.Printf        (printf)
import           Data.List          (elemIndex, isPrefixOf)
import           Data.List.Split    (chunksOf)
import           Data.Maybe         (fromJust, isJust)
import qualified Data.Map           as M
import qualified Data.Set           as S
import qualified Data.Vector        as V
import           Data.Text.Lazy     (pack)
import           Control.Monad      (void)
import           Control.Concurrent  (Chan, newChan, writeChan, forkIO)
import qualified Data.Graph.Inductive as G  hiding (run)
import qualified Data.GraphViz        as Gv hiding (parse, toLabel)
import qualified Data.GraphViz.Attributes.Complete as GA
import qualified Data.GraphViz.Attributes.Colors   as GC
import           Data.GraphViz.Attributes.Colors   (Color(RGB))
-- import           Test.HUnit          ((~:),(~=?),Test,test)

import qualified Data.Clustering.Hierarchical as C

import           Bio.Phylogeny.PhyBin.CoreTypes

----------------------------------------------------------------------------------------------------
-- Visualization with GraphViz and FGL:
----------------------------------------------------------------------------------------------------

-- First we need to be able to convert our trees to FGL graphs:
toGraph :: FullTree StandardDecor -> G.Gr String Double
toGraph (FullTree _ tbl tree) = G.run_ G.empty $ loop tree
  where
 fromLabel ix = tbl M.! ix
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

----------------------------------------------------------------------------------------------------
-- Dendrograms
 ----------------------------------------------------------------------------------------------------

-- | Some duplicated code with dotNewickTree.
dotDendrogram :: PhyBinConfig -> String -> Double -> C.Dendrogram (FullTree a) -> Maybe (M.Map TreeName Int) -> Gv.DotGraph G.Node
dotDendrogram PBC{highlights, show_trees_in_dendro} title edge_scale origDendro mNameMap  =
  Gv.graphToDot myparams (G.nmap uid graph)
 where
  (charsDropped, dendro) = truncateNames origDendro
  -- This is ugly, but we modify the name map to match:
  nameMap' = fmap (M.mapKeys (drop charsDropped)) mNameMap  
  graph :: DendroGraph
  graph = dendrogramToGraph dendro

  uidsToNames = M.fromList $
                map (\NdLabel{uid,tre=Just t,clumpSize} -> (uid,(t,clumpSize))) $ 
                filter (isJust . tre) $ 
                map (fromJust . G.lab graph) $  
                G.nodes graph

  myparams :: Gv.GraphvizParams G.Node String Double () String -- (Either String String)
  myparams = Gv.defaultParams { Gv.globalAttributes= [Gv.GraphAttrs [GA.Label$ GA.StrLabel$ pack title]],
                                Gv.fmtNode= nodeAttrs,
                                Gv.fmtEdge= edgeAttrs
                              }
  nodeAttrs :: (Int, UniqueNodeName) -> [GA.Attribute]
  nodeAttrs (_num, uid) =
    let uid' = if show_trees_in_dendro
               then printed_tree++"\n"++uid
               else uid
        printed_tree =
          case M.lookup uid uidsToNames of
            Nothing -> ""
            Just (ft,sz) -> (if sz>1 then "size "++show sz++"\n" else "")
                            ++ show (displayStrippedTree ft)
        (tag,shp) = -- case eith of
          if isPrefixOf "DUMMY_" uid
          then ("", GA.PointShape)          
          else (uid', GA.Ellipse)
    in 
    [ GA.Label$ GA.StrLabel$ pack tag
    , GA.Shape shp
    , GA.Style [GA.SItem GA.Filled []]
    ] ++
    case (nameMap', M.lookup uid uidsToNames) of
      (Just nm, Just (FullTree{treename=trN},_)) ->
        case M.lookup trN nm of
          Nothing -> []
          -- Here we color the top TEN clusters:
          Just ind | ind <= 10 -> [ GA.Color [weighted$ defaultPalette V.! (ind-1) ] ]
                   | otherwise -> []
      _ -> []
  edgeAttrs = getEdgeAttrs edge_scale


type UniqueNodeName = String
-- type DendroGraph = G.Gr (Maybe TreeName,UniqueNodeName) Double
type DendroGraph = G.Gr NdLabel Double

-- | When we first convert to a graph representation, there is a bunch of information
-- hanging off of each node.
data NdLabel =
  NdLabel
  { uid       :: UniqueNodeName
  , tre       :: Maybe (FullTree ())
  , clumpSize :: !Int
  }
  deriving (Show, Ord, Eq)

-- | Create a graph using TreeNames for node labels and edit-distance for edge labels.
dendrogramToGraph :: C.Dendrogram (FullTree a) -> DendroGraph
dendrogramToGraph orig =
  G.run_ G.empty $ void$
    loop (fmap (fmap (const())) orig)
  where
 loop node@(C.Leaf ft) = G.insMapNodeM (NdLabel (treename ft) (Just$ fmap (const()) ft) 1)
 loop node@(C.Branch 0 left right) = do
   -- As a preprocessing step we collapse clusters that are separated by zero edit distance.
   ----------------------------------------
   let lvs = collapseZeroes left ++ collapseZeroes right
       nms = map (treename) lvs
       lens = map length nms
       total = sum lens
       avg   = total `quot` length nms
       -- The goal here is to make an approximately square arrangement:
       -- goal: avg * perline == total / (avg * perline)
       perline = ceiling$ sqrt (fromIntegral total / ((fromIntegral avg)^2))
       chunked = chunksOf perline nms
       fatname = unlines (map unwords chunked)
   G.insMapNodeM (NdLabel fatname (Just (head lvs)) (length lvs))
   ----------------------------------------   
 loop node@(C.Branch dist left right) =
    do (_,ll@(NdLabel lid _ s1)) <- loop left
       (_,rr@(NdLabel rid _ s2)) <- loop right
       -- Interior nodes do NOT have their names drawn:
       let ndname = "DUMMY_"++(lid++"_"++rid)  -- HACK!
       (midN,mid) <- G.insMapNodeM (NdLabel ndname Nothing (s1+s2))
       G.insMapEdgeM (ll, mid, dist)
       G.insMapEdgeM (rr, mid, dist)
       return (midN,mid)

 collapseZeroes (C.Leaf tr)      = [tr]
 collapseZeroes (C.Branch 0 l r) = collapseZeroes l ++ collapseZeroes r
 collapseZeroes oth = error "dendrogramToGraph: internal error.  Not expecting non-zero branch length here."


-- | The plot looks nicer when the names aren't bloated with repeated stuff. This
-- replaces all tree names with potentially shorter names by removing the common prefix.
-- Returns how many characters were dropped.
truncateNames :: C.Dendrogram (FullTree a) -> (Int, C.Dendrogram (FullTree a))
truncateNames dendro = (prefChars, fmap chopName dendro)
 where 
  chopName ft = ft{ treename= drop prefChars (treename ft) }
  prefChars = length$ commonPrefix$ S.toList$ allNames dendro
  allNames (C.Leaf tr)      = S.singleton (treename tr)
  allNames (C.Branch _ l r) = S.union (allNames l) (allNames r)

----------------------------------------------------------------------------------------------------


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
myShowFloat fl =
  let rnd = round fl in
  if fl == fromIntegral rnd
  then show rnd
  else printf "%.4f" fl

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
  fromLabel ix = tbl M.! ix  
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
  edgeAttrs = getEdgeAttrs edge_scale


-- | Some arbitrarily chosen colors from the X11 set:
defaultPalette :: V.Vector GA.Color
defaultPalette = V.fromList$ concat$ replicate 4 $ map GA.X11Color 
  [ Gv.Aquamarine
  , Gv.PaleVioletRed
  , Gv.MediumPurple
  , Gv.PaleGreen
  , Gv.PapayaWhip
  , Gv.SkyBlue
  , Gv.Yellow
  , Gv.Crimson
  , Gv.Gray
  , Gv.PaleGoldenrod
  ]

-- Grabbed the first couple palettes off a website:
altPalette :: V.Vector GA.Color
altPalette = V.fromList $ concat $ replicate 3 $ 
  -- http://www.colourlovers.com/palette/2962435/Autumn_Roses
  [ RGB 159 74 81, RGB 217 183 173, RGB 149 91 116, RGB 185 138 148
  --  , RGB 101 69 82 -- too dark
  ] ++
  -- http://www.colourlovers.com/palette/2962434/Earthy_warm
  [ RGB 108 74 39, RGB 207 179 83, RGB 180 149 60, RGB 244 242 185
  -- , RGB 61 63 39
  ]
  
getEdgeAttrs :: Double -> (t, t1, Double) -> [GA.Attribute]
getEdgeAttrs edge_scale = edgeAttrs
 where 
  -- TOGGLE:
  --  edgeAttrs (_,_,weight) = [ArrowHead noArrow, Len (weight * edge_scale + bump), GA.Label (StrLabel$ show (weight))]
  edgeAttrs (_,_,weight) = 
                           let draw_weight = compute_draw_weight weight edge_scale in
                           --trace ("EDGE WEIGHT "++ show weight ++ " drawn at "++ show draw_weight) $
			   [GA.ArrowHead Gv.noArrow,
                            GA.Label$ GA.StrLabel$ pack$ myShowFloat weight] ++ -- TEMPTOGGLE
			   --[ArrowHead noArrow, GA.Label (StrLabel$ show draw_weight)] ++ -- TEMPTOGGLE
			    if weight == 0.0
			    then [GA.Color [weighted$ GA.X11Color Gv.Red],
                                  GA.LWidth 3.0, GA.Len minlen]
			    else [GA.Len draw_weight]
  minlen = 0.7
  maxlen = 3.0
  compute_draw_weight w scale = 
    let scaled = (abs w) * scale + minlen in 
    -- Don't draw them too big or it gets annoying:
    (min scaled maxlen)

weighted c = GC.WC {GC.wColor=c, GC.weighting=Nothing}

-- | This version shows the ordered/rooted structure of the normalized tree.
dotNewickTree_debug :: String -> FullTree StandardDecor -> Gv.DotGraph G.Node
dotNewickTree_debug title atree@(FullTree _ tbl tree) = Gv.graphToDot myparams graph
 where 
  graph = toGraph2 atree
  fromLabel ix = tbl M.! ix    
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

-- | Common prefix of a list of lists.
commonPrefix :: Eq a => [[a]] -> [a]
commonPrefix [] = []
commonPrefix ls@(hd:tl)
  | any null ls = []
  | otherwise =
    if   all ((== (head hd)) . head) tl
    then head hd : commonPrefix (map tail ls)
    else           commonPrefix (map tail ls)