{-# LANGUAGE ScopedTypeVariables, CPP #-}

module Bio.Phylogeny.PhyBin.RFDistance
       where

import           Control.Monad
import           Data.Word
import qualified Data.Vector.Unboxed.Mutable as M
import qualified Data.Vector.Unboxed         as U
import qualified Data.Vector.Unboxed.Bit     as UB
import qualified Data.Bit                    as B
import           Text.PrettyPrint.HughesPJClass hiding (char, Style)

import           Control.LVish
import qualified Data.LVar.Set   as IS
import qualified Data.LVar.SLSet as SL

import           Data.LVar.Map   as IM
import           Data.LVar.NatArray as NA

import           Bio.Phylogeny.PhyBin.CoreTypes
import           Data.BitList
import qualified Data.Set as S
import qualified Data.Foldable as F
import           Data.Monoid

--------------------------------------------------------------------------------

-- | Tree's are identified simply by their order within the list of input trees.
type TreeID = Int

-- | A collection of all observed bipartitons (bips) with a mapping of which trees
-- contain which Bips.
type BipTable s = IMap DenseLabelSet s (SparseLabelSet s)
-- type BipTable = IMap BitList (U.Vector Bool)
-- type BipTable s = IMap BitList s (NA.NatArray s Word8)

-- | Sets of taxa (BiPs) that are expected to be sparse.
type SparseLabelSet s = IS.ISet s Label
-- NA.NatArray s Word8

--------------------------------------------------------------------------------

-- | Dense sets of taxa, aka Bipartitions or BiPs
--   We assume that taxa labels have been mapped onto a dense, contiguous range of integers [0,N). 
type DenseLabelSet = S.Set Label
-- type DenseLabelSet s = BitList
-- type DenseLabelSet = UB.Vector B.Bit

-- M.write vec lab (B.fromBool True)
-- mkEmptyDense size = U.replicate size (B.fromBool False)    

-- markLabel lab set = IS.putInSet lab set 
-- mkEmptyDense _size = IS.newEmptySet

markLabel    :: Label -> DenseLabelSet -> DenseLabelSet
mkEmptyDense :: Int -> DenseLabelSet
denseUnions  :: Int -> [DenseLabelSet] -> DenseLabelSet
bipSize      :: DenseLabelSet -> Int

markLabel lab set  = S.insert lab set 
mkEmptyDense _size = S.empty
denseUnions _size  = S.unions 
bipSize            = S.size

--------------------------------------------------------------------------------

-- The distance matrix is an atomically-bumped matrix of numbers.
-- type DistanceMat s = NA.NatArray s Word32
-- Except... bump isn't supported by our idempotent impl.

-- | Returns a distance matrix as an adjacency matrix encoded as a vector.
distanceMatrix :: [AnnotatedTree] -> IO (U.Vector Word)
distanceMatrix lst = runParIO$ do   

  table :: BipTable s <- IM.newEmptyMap 
  forM_ lst $ \ tree -> do
    insertBips table tree
  
  undefined

  -- runParThenFreeze -- get bip table
  -- followed by ... fill matrix from bip table


insertBips :: BipTable s -> AnnotatedTree -> Par d s ()
insertBips table tree = loop
  where
  loop = undefined

-- | The number of bipartitions implied by a tree is one per EDGE in the tree.  Thus
-- each interior node carries a list of BiPs the same length as its list of children.
labelBips :: NewickTree a -> NewickTree (a, [DenseLabelSet])
labelBips tr = loop tr
  where
    size = treeSize tr    
    zero = mkEmptyDense size
    loop (NTLeaf dec lab) = NTLeaf (dec, [markLabel lab zero]) lab      
    loop (NTInterior dec chlds) =
      let chlds' = map loop chlds
          sets   = map (denseUnions size . snd . get_dec) chlds' in
      NTInterior (dec, sets) chlds'

foldBips :: Monoid m => (DenseLabelSet -> m) -> NewickTree a -> m
foldBips f tr = F.foldMap f' (labelBips tr)
 where
   f' (_,bips) = F.foldMap f bips
  
-- | Get all non-singleton BiPs implied by a tree.
allBips tr = S.filter ((> 1) . bipSize) $ foldBips S.insert tr S.empty


instance Pretty a => Pretty (S.Set a) where
 pPrint s = pPrint (S.toList s)
 
