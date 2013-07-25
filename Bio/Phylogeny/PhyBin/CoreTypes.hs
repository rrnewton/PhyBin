{-# LANGUAGE NamedFieldPuns, BangPatterns #-}
{-# LANGUAGE OverloadedStrings, TypeSynonymInstances, FlexibleInstances #-}
{-# LANGUAGE CPP #-}
{-# OPTIONS_GHC -fwarn-incomplete-patterns #-}
{-# OPTIONS_GHC -fwarn-unused-imports #-}

module Bio.Phylogeny.PhyBin.CoreTypes
       (
         -- * Tree and tree decoration types
         NewickTree(..), 
         DefDecor, StandardDecor(..), AnnotatedTree, FullTree(..),
         ClustMode(..), TreeName,
         
         -- * Tree operations
         displayDefaultTree, displayStrippedTree, 
         treeSize, numLeaves, liftFT,
         get_dec, set_dec, get_children, 
         map_labels, all_labels, foldIsomorphicTrees,

         -- * Utilities specific to StandardDecor:
         avg_branchlen, get_bootstraps,

         -- * Command line config options
         PhyBinConfig(..), default_phybin_config, 

         -- * General helpers
         Label, LabelTable,
         
         -- * Experimenting with abstracting decoration operations
         HasBranchLen(..)
       )
       where

import qualified Data.Map as M
import Data.Foldable (Foldable(..))
import Data.Maybe (maybeToList)
import Data.Monoid (mappend, mconcat)
import Text.PrettyPrint.HughesPJClass hiding (char, Style)

import qualified Data.Clustering.Hierarchical as C

#define NO_ATOMS
#ifndef NO_ATOMS
import StringTable.Atom
#endif

----------------------------------------------------------------------------------------------------
-- Type definitions
----------------------------------------------------------------------------------------------------

type BranchLen = Double

-- | Even though the Newick format allows it, here we ignore interior node
--   labels. (They are not commonly used.)
--
--   Note that these trees are rooted.  The `normalize` function ensures that a
--   single, canonical rooted representation is chosen.
data NewickTree a = 
   NTLeaf     a {-# UNPACK #-} !Label
 | NTInterior a  [NewickTree a]
 deriving (Show, Eq, Ord)

-- TODO: Ordering maybe shouldn't need to touch the metadata.  At least on the fast
-- path.

{-
-- [2010.09.22] Disabling:
instance NFData Atom where
  rnf a = rnf (fromAtom a :: Int)

instance NFData a => NFData (NewickTree a) where
  rnf (NTLeaf l n)      = rnf (l,n)
  rnf (NTInterior l ls) = rnf (l,ls)
-}

instance Functor NewickTree where 
   fmap fn (NTLeaf dec x)      = NTLeaf (fn dec) x 
   fmap fn (NTInterior dec ls) = NTInterior (fn dec) (map (fmap fn) ls)

instance Foldable NewickTree where
  foldMap f (NTLeaf dec x) = f dec
  foldMap f (NTInterior dec ls) = mappend (f dec) $
                                  mconcat (map (foldMap f) ls)

instance Foldable FullTree where
  foldMap f (FullTree _ _ tr) = foldMap f tr

instance Pretty dec => Pretty (NewickTree dec) where 
 -- I'm using displayDefaultTree for the "prettiest" printing and
 -- replacing pPrint with a whitespace-improved version of show:
 pPrint (NTLeaf dec name)   = "NTLeaf"     <+> pPrint dec <+> text (show name)
 pPrint (NTInterior dec ls) = "NTInterior" <+> pPrint dec <+> pPrint ls

instance Pretty a => Pretty (FullTree a) where
  pPrint (FullTree name mp tr) = 
    "FullTree " <+> text name <+> loop tr
   where
    loop (NTLeaf dec ind)    = "NTLeaf"     <+> pPrint dec <+> text (mp M.! ind)
    loop (NTInterior dec ls) = "NTInterior" <+> pPrint dec <+> pPrint ls

instance (Pretty k, Pretty v) => Pretty (M.Map k v) where
  pPrint mp = pPrint (M.toList mp)

-- | Display a tree WITH the bootstrap and branch lengths.
displayDefaultTree :: FullTree DefDecor -> Doc
displayDefaultTree orig = loop tr <> ";"
  where
    (FullTree _ mp tr) = orig -- normalize orig
    loop (NTLeaf (Nothing,_) name)     = text (mp M.! name)
    loop (NTLeaf _ _)                  = error "WEIRD -- why did a leaf node have a bootstrap value?"
    loop (NTInterior (bootstrap,_) ls) = 
       case bootstrap of
         Nothing -> base
         Just val -> base <> text ":[" <> text (show val) <> text "]"
      where base = parens$ sep$ map_but_last (<>text",") $ map loop ls

displayStrippedTree :: FullTree a -> Doc
displayStrippedTree orig = loop tr <> ";"
  where
    (FullTree _ mp tr) = orig -- normalize orig
    loop (NTLeaf _ name) = text (mp M.! name)
    loop (NTInterior _ ls) = parens$ sep$ map_but_last (<>text",") $ map loop ls

----------------------------------------------------------------------------------------------------
-- Labels
----------------------------------------------------------------------------------------------------


-- Experimental: toggle this to change the representation of labels:
-- Alas I always have problems with the interned string libs (e.g. segfaults)... [2012.11.20]
----------------------------------------
type Label = Int

-- | Map labels back onto meaningful names.
type LabelTable = M.Map Label String

----------------------------------------------------------------------------------------------------
-- Tree metadata (decorators)
----------------------------------------------------------------------------------------------------


-- | The barebones default decorator for NewickTrees contains BOOTSTRAP and
-- BRANCHLENGTH.  The bootstrap values, if present, will range in [0..100]
type DefDecor = (Maybe Int, BranchLen)

-- | Additionally includes some scratch data that is used by the binning algorithm.
type AnnotatedTree = NewickTree StandardDecor

-- | The standard decoration includes everything in `DefDecor` plus
--   some extra cached data:
-- 
--  (1) branch length from parent to "this" node
--  (2) bootstrap values for the node
-- 
--  (3) subtree weights for future use
--      (defined as number of LEAVES, not counting intermediate nodes)
--  (4) sorted lists of labels for symmetry breaking
data StandardDecor = StandardDecor {
  branchLen     :: BranchLen,
  bootStrap     :: Maybe Int,

  -- The rest of these are used by the computations below.  These are
  -- cached (memoized) values that could be recomputed:
  ----------------------------------------
  subtreeWeight :: Int,
  sortedLabels  :: [Label]
 }
 deriving (Show,Read,Eq,Ord)

class HasBranchLen a where
  getBranchLen :: a -> BranchLen

instance HasBranchLen StandardDecor where
  getBranchLen = branchLen

-- This one is kind of sloppy:
instance HasBranchLen DefDecor where
  getBranchLen = snd

-- | A common type of tree contains the standard decorator and also a table for
-- restoring the human-readable node names.
data FullTree a =
  FullTree { treename   :: TreeName
           , labelTable :: LabelTable
           , nwtree     :: NewickTree a 
           }
 deriving (Show)

liftFT :: (NewickTree t -> NewickTree a) -> FullTree t -> FullTree a
liftFT fn (FullTree nm labs x) = FullTree nm labs (fn x)

type TreeName = String

instance Pretty StandardDecor where 
 pPrint (StandardDecor bl bs wt ls) = parens$
    "StandardDecor" <+> hsep [pPrint bl, pPrint bs
--                             , pPrint wt, pPrint ls
                             ]


----------------------------------------------------------------------------------------------------
-- * Configuring and running the command line tool.
----------------------------------------------------------------------------------------------------

-- | Due to the number of configuration options for the driver, we pack them into a record.
data PhyBinConfig = 
  PBC { verbose :: Bool
      , num_taxa :: Int
      , name_hack :: String -> String
      , output_dir :: String
      , inputs :: [String]
      , do_graph :: Bool
      , do_draw :: Bool
      , clust_mode  :: ClustMode
      , highlights  :: [FilePath] -- [FullTree ()]
      , show_trees_in_dendro :: Bool
      , use_hashrf  :: Bool  
      , print_rfmatrix :: Bool
      , dist_thresh :: Maybe Int
      , branch_collapse_thresh :: Maybe Double -- ^ Branches less than this length are collapsed.
      }

-- | The default phybin configuration.
default_phybin_config :: PhyBinConfig
default_phybin_config = 
 PBC { verbose = False
      , num_taxa = error "must be able to determine the number of taxa expected in the dataset.  (Supply it manually with -n.)"
      , name_hack = id -- Default, no transformation of leaf-labels
      , output_dir = "./phybin_out/"
      , inputs = []
      , do_graph = False
      , do_draw = False
      , clust_mode = BinThem
#ifdef USE_HASHRF
      , use_hashrf = True
#else
      , use_hashrf = False
#endif
      , highlights     = []
      , show_trees_in_dendro = False
      , print_rfmatrix = False
      , dist_thresh = Nothing
      , branch_collapse_thresh = Nothing
     }

data ClustMode = BinThem | ClusterThem { linkage :: C.Linkage }

----------------------------------------------------------------------------------------------------
-- * Simple utility functions for the core types:
----------------------------------------------------------------------------------------------------

-- | How many nodes (leaves and interior) are contained in a NewickTree?
treeSize :: NewickTree a -> Int
treeSize (NTLeaf _ _) = 1
treeSize (NTInterior _ ls) = 1 + sum (map treeSize ls)

-- | This counts only leaf nodes, which should include all taxa.
numLeaves :: NewickTree a -> Int
numLeaves (NTLeaf _ _) = 1
numLeaves (NTInterior _ ls) = sum (map numLeaves ls)


map_but_last :: (a -> a) -> [a] -> [a]
map_but_last _ [] = []
map_but_last _ [h] = [h]
map_but_last fn (h:t) = fn h : map_but_last fn t



get_dec :: NewickTree t -> t
get_dec (NTLeaf     dec _) = dec
get_dec (NTInterior dec _) = dec

-- Set all the decorations to a constant:
set_dec :: b -> NewickTree a -> NewickTree b
set_dec d = fmap (const d)
--set_dec d (NTLeaf _ x) = NTLeaf d x
--set_dec d (NTInterior _ ls) = NTInterior d $ map (set_dec d) ls

get_children :: NewickTree t -> [NewickTree t]
get_children (NTLeaf _ _) = []
get_children (NTInterior _ ls) = ls


-- | Average branch length across all branches in all all trees.
avg_branchlen :: HasBranchLen a => [NewickTree a] -> Double
avg_branchlen origls = fst total / snd total
  where
   total = sum_ls $ map sum_tree origls
   sum_ls ls = (sum$ map fst ls, sum$ map snd ls)
   sum_tree (NTLeaf dec _) | getBranchLen dec == 0  = (0,0)
                           | otherwise              = (abs (getBranchLen dec),1)
   sum_tree (NTInterior dec ls) = 
       let branchLen = getBranchLen dec
           (x,y)     = sum_ls$ map sum_tree ls in
       if branchLen == 0 then (x, y) else ((abs branchLen) + x, 1+y)

-- | Retrieve all the bootstraps values actually present in a tree.
get_bootstraps :: NewickTree StandardDecor -> [Int]
get_bootstraps (NTLeaf (StandardDecor{bootStrap}) _) = maybeToList bootStrap
get_bootstraps (NTInterior (StandardDecor{bootStrap}) ls) =
  maybeToList bootStrap ++ concatMap get_bootstraps ls

-- | Apply a function to all the *labels* (leaf names) in a tree.
map_labels :: (Label -> Label) -> NewickTree a -> NewickTree a
map_labels fn (NTLeaf     dec lbl) = NTLeaf dec $ fn lbl
map_labels fn (NTInterior dec ls)  = NTInterior dec$ map (map_labels fn) ls
 
-- | Return all the labels contained in the tree.
all_labels :: NewickTree t -> [Label]
all_labels (NTLeaf     _ lbl) = [lbl]
all_labels (NTInterior _ ls)  = concat$ map all_labels ls



-- | This function allows one to collapse multiple trees while looking
-- only at the "horizontal slice" of all the annotations *at a given
-- position* in the tree.
--
-- "Isomorphic" must apply both to the shape and the name labels or it
-- is an error to apply this function.
foldIsomorphicTrees :: ([a] -> b) -> [NewickTree a] -> NewickTree b
foldIsomorphicTrees _ [] = error "foldIsomorphicTrees: empty list of input trees"
foldIsomorphicTrees fn ls@(hd:_) = fmap fn horiztrees
 where
   -- Preserve the input order:
   horiztrees = Prelude.foldr consTrees (fmap (const []) hd) ls
   -- We use the tree datatype itself as the intermediate data
   -- structure.  This is VERY allocation-expensive, it would be
   -- possible to trade compute for allocation here:
   consTrees a b = case (a,b) of
    (NTLeaf dec nm1, NTLeaf decls nm2) | nm1 /= nm2 -> error$"foldIsomorphicTrees: mismatched names: "++show (nm1,nm2)
                                       | otherwise ->
     NTLeaf (dec : decls) nm1
    (NTInterior dec ls1, NTInterior decls ls2) ->
     NTInterior (dec:decls) $ zipWith consTrees ls1 ls2
    _ -> error "foldIsomorphicTrees: difference in tree shapes"

