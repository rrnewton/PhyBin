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

--------------------------------------------------------------------------------

-- type BiPartTable = IMap BitList (U.Vector Bool)
-- type BiPartTable s = IMap BitList s (NA.NatArray s Word8)
type BiPartTable s = IMap DenseLabelSet s (SparseLabelSet s)

type SparseLabelSet s = IS.ISet s Label
-- NA.NatArray s Word8

--------------------------------------------------------------------------------

-- We assume that taxa labels have been mapped onto a dense, contiguous range of integers [0,N).
-- type DenseLabelSet s = BitList
-- type DenseLabelSet = UB.Vector B.Bit
type DenseLabelSet = S.Set Label

-- M.write vec lab (B.fromBool True)
-- mkEmptyDense size = U.replicate size (B.fromBool False)    

-- markLabel lab set = IS.putInSet lab set 
-- mkEmptyDense _size = IS.newEmptySet

markLabel :: Label -> DenseLabelSet -> DenseLabelSet
markLabel lab set  = S.insert lab set 

mkEmptyDense :: Int -> DenseLabelSet
mkEmptyDense _size = S.empty

denseUnions :: Int -> [DenseLabelSet] -> DenseLabelSet
denseUnions _size  = S.unions 

--------------------------------------------------------------------------------

-- The distance matrix is an atomically-bumped matrix of numbers.
-- type DistanceMat s = NA.NatArray s Word32
-- Except... bump isn't supported by our idempotent impl.

-- | Returns a distance matrix as an adjacency matrix encoded as a vector.
distanceMatrix :: [AnnotatedTree] -> IO (U.Vector Word)
distanceMatrix lst = runParIO$ do   

  table::BiPartTable s <- IM.newEmptyMap 
  forM_ lst $ \ tree -> do
    insertBips table tree
  
  undefined

  -- runParThenFreeze -- get bip table
  -- followed by ... fill matrix from bip table


insertBips :: BiPartTable s -> AnnotatedTree -> Par d s ()
insertBips table tree = loop
  where
  loop = undefined


-- labelBips :: AnnotatedTree -> NewickTree (StandardDecor, BitList)
-- labelBips :: AnnotatedTree -> NewickTree (StandardDecor, U.Vector Bool)
labelBips :: AnnotatedTree -> NewickTree (StandardDecor, DenseLabelSet)
labelBips tr = loop tr
  where
    size = treeSize tr    
    zero = mkEmptyDense size
    loop (NTLeaf dec lab) =
      let bitvec = markLabel lab zero in
      NTLeaf (dec,bitvec) lab
      
    loop (NTInterior dec chlds) =
      let chlds' = map loop chlds
          sets   = map (snd . get_dec) chlds' in
      NTInterior (dec, denseUnions size sets) chlds'      

    -- loop2 acc [] = []
    -- loop2 acc (hd:tl) = UB.difference (get_dec hd) acc
    

#if 0
    zero :: U.Vector Bool
    zero = U.replicate size False
    loop (NTLeaf dec lab) =
      let bitvec = U.modify (\v -> M.write v (read lab) True) zero in
      NTLeaf (dec,bitvec) lab
#endif


instance Pretty a => Pretty (S.Set a) where
 pPrint s = pPrint (S.toList s)
 
