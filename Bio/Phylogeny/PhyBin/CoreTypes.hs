{-# OPTIONS_GHC -fwarn-incomplete-patterns #-}
{-# OPTIONS_GHC -fwarn-unused-imports #-}

module Bio.Phylogeny.PhyBin.CoreTypes
       (
         NewickTree(..), displayDefaultTree,
         DefDecor, StandardDecor(..),
         
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


