{-# LANGUAGE NamedFieldPuns #-}
{-# OPTIONS_GHC -fwarn-incomplete-patterns #-}
{-# OPTIONS_GHC -fwarn-unused-imports #-}

module Bio.Phylogeny.PhyBin.CoreTypes
       (
         NewickTree(..), displayDefaultTree,
         DefDecor, StandardDecor(..),
         get_dec, set_dec, get_children, avg_branchlen,
         
         PhyBinConfig(..), default_phybin_config,

         toLabel, fromLabel, Label,
       )       
       where

import Text.PrettyPrint.HughesPJClass hiding (char, Style)

----------------------------------------------------------------------------------------------------
-- Type definitions
----------------------------------------------------------------------------------------------------

type BranchLen = Double

-- | Even though the Newick format allows it, here we ignore interior node
--   labels. (They are not commonly used.)
data NewickTree a = 
   NTLeaf     a Label
 | NTInterior a [NewickTree a]
 deriving (Show, Eq, Ord)

{-
-- [2010.09.22] Disabling:
instance NFData Atom where
  rnf a = rnf (fromAtom a :: Int)

instance NFData a => NFData (NewickTree a) where
  rnf (NTLeaf l n)      = rnf (l,n)
  rnf (NTInterior l ls) = rnf (l,ls)
-}

instance Pretty (NewickTree dec) where 
 pPrint (NTLeaf _ name)   = text (fromLabel name)
 pPrint (NTInterior _ ls) = 
     --parens$ commasep ls
     (parens$ sep$ map_but_last (<>text",") $ map pPrint ls)

instance Functor NewickTree where 
   fmap fn (NTLeaf dec x)      = NTLeaf (fn dec) x 
   fmap fn (NTInterior dec ls) = NTInterior (fn dec) (map (fmap fn) ls)

----------------------------------------------------------------------------------------------------
-- Labels
----------------------------------------------------------------------------------------------------


-- Experimental: toggle this to change the representation of labels:
----------------------------------------
--type Label = Atom; (toLabel, fromLabel) = (toAtom, fromAtom)
----------------------------------------
type Label = String; (toLabel, fromLabel) = (id, id)
----------------------------------------
fromLabel :: Label -> String
toLabel   :: String -> Label

----------------------------------------------------------------------------------------------------
-- Tree metadata (decorators)
----------------------------------------------------------------------------------------------------


-- | The default decorator for NewickTrees contains BOOTSTRAP and BRANCHLENGTH.
--   The bootstrap values, if present, will range in [0..100]
type DefDecor = (Maybe Int, BranchLen)


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


----------------------------------------------------------------------------------------------------
-- * Configuring and running the command line tool.
----------------------------------------------------------------------------------------------------

-- | Due to the number of configuration options for the driver, we pack them into a record.
data PhyBinConfig = 
  PBC { verbose :: Bool
      , num_taxa :: Int
      , name_hack :: Label -> Label
      , output_dir :: String
      , inputs :: [String]
      , do_graph :: Bool
      , do_draw :: Bool
      }

-- | The default phybin configuration.
default_phybin_config :: PhyBinConfig
default_phybin_config = 
 PBC { verbose = False
      , num_taxa = error "must be able to determine the number of taxa expected in the dataset.  (Supply it manually.)"
      , name_hack = id -- Default, no transformation of leaf-labels
      , output_dir = "./"
      , inputs = []
      , do_graph = False
      , do_draw = False
     }

----------------------------------------------------------------------------------------------------
-- * Simple utility functions for the core types:
----------------------------------------------------------------------------------------------------

-- | Display a tree WITH the bootstrap and branch lengths.
displayDefaultTree :: NewickTree DefDecor -> Doc
displayDefaultTree (NTLeaf (Nothing,_) name)   = text (fromLabel name)
displayDefaultTree (NTLeaf _ _ ) = error "WEIRD -- why did a leaf node have a bootstrap value?"
displayDefaultTree (NTInterior (bootstrap,_) ls) = 
   case bootstrap of
     Nothing -> base
     Just val -> base <> text ":[" <> text (show val) <> text "]"
 where
   base = parens$ sep$ map_but_last (<>text",") $ map pPrint ls

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
avg_branchlen :: [NewickTree StandardDecor] -> Double
avg_branchlen origls = fst total / snd total
  where
   total = sum_ls $ map sum_tree origls
   sum_ls ls = (sum$ map fst ls, sum$ map snd ls)
   sum_tree (NTLeaf (StandardDecor{branchLen=0}) _)    = (0,0)
   sum_tree (NTLeaf (StandardDecor{branchLen}) _)      = (abs branchLen,1)
   sum_tree (NTInterior (StandardDecor{branchLen}) ls) = 
       let (x,y) = sum_ls$ map sum_tree ls in
       if branchLen == 0 then (x, y) else ((abs branchLen) + x, 1+y)
