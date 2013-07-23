{-# LANGUAGE ScopedTypeVariables, CPP #-}

module Bio.Phylogeny.PhyBin.RFDistance
       (DenseLabelSet, DistanceMatrix, 
        allBips, foldBips, dispBip, 
        distanceMatrix, printDistMat)
       where

import           Control.Monad
import           Data.Word
import qualified Data.Vector                 as V
import qualified Data.Vector.Unboxed.Mutable as MV
import qualified Data.Vector.Unboxed         as U
import qualified Data.Vector.Unboxed.Bit     as UB
import qualified Data.Bit                    as B
import           Text.PrettyPrint.HughesPJClass hiding (char, Style)
import           System.IO      (hPutStrLn, hPutStr, Handle)

-- import           Control.LVish
-- import qualified Data.LVar.Set   as IS
-- import qualified Data.LVar.SLSet as SL

-- import           Data.LVar.Map   as IM
-- import           Data.LVar.NatArray as NA

import           Bio.Phylogeny.PhyBin.CoreTypes
-- import           Data.BitList
import qualified Data.Set as S
import qualified Data.IntSet as SI
import qualified Data.Map.Strict as M
import qualified Data.Foldable as F
import           Data.Monoid
import           Prelude as P
import           Debug.Trace

--------------------------------------------------------------------------------
-- A data structure choice
--------------------------------------------------------------------------------

-- | Dense sets of taxa, aka Bipartitions or BiPs
--   We assume that taxa labels have been mapped onto a dense, contiguous range of integers [0,N).
-- 
--   NORMALIZATION Rule: Bipartitions are really two disjoint sets.  But as long as
--   the parent set (the union of the partitions, aka "all taxa") then a bipartition
--   can be represented just by *one* subset.  Yet we must choose WHICH subset for
--   consistency.  We use the rule that we always choose the SMALLER.  Thus the
--   DenseLabelSet should always be half the size or less, compared to the total
--   number of taxa.
-- 
--   A set that is more than a majority of the taxa can be normalized by "flipping",
--   i.e. taking the taxa that are NOT in that set.
#if 1
-- type DenseLabelSet s = BitList
type DenseLabelSet = UB.Vector B.Bit
markLabel lab = UB.modify (\vec -> MV.write vec lab (B.fromBool True)) 
mkEmptyDense  size = UB.replicate size (B.fromBool False)
mkSingleDense size ind = markLabel ind (mkEmptyDense size)
denseUnions        = UB.unions
bipSize            = U.length
denseDiff          = UB.difference

dispBip labs bip = show$ map (\(ix,_) -> (labs M.! ix)) $
                        filter (\(_,bit) -> B.toBool bit) $
                        zip [0..] (UB.toList bip)
#else
type DenseLabelSet = SI.IntSet
markLabel lab set   = SI.insert lab set 
mkEmptyDense _size  = SI.empty
mkSingleDense _size = SI.singleton
denseUnions _size   = SI.unions 
bipSize             = SI.size
denseDiff           = SI.difference

dispBip labs bip = show$ map (map (labs M.!)) $ 
                         map SI.toList bip
#endif

markLabel    :: Label -> DenseLabelSet -> DenseLabelSet
mkEmptyDense :: Int -> DenseLabelSet
denseUnions  :: Int -> [DenseLabelSet] -> DenseLabelSet
bipSize      :: DenseLabelSet -> Int

-- | Print a BiPartition in a pretty form
dispBip      :: LabelTable -> DenseLabelSet -> String

--------------------------------------------------------------------------------
-- Dirt-simple reference implementation
--------------------------------------------------------------------------------

type DistanceMatrix = V.Vector (U.Vector Int)

-- | Returns a triangular distance matrix encoded as a vector.
distanceMatrix :: [NewickTree a] -> DistanceMatrix
distanceMatrix lst = 
   let sz = P.length lst
       eachbips = V.fromList $ map allBips lst
--   in V.generate (sz-1) $ \ i ->
   in V.generate sz $ \ i ->        
      U.generate i  $ \ j ->
      S.size (S.difference (eachbips V.! i) (eachbips V.! j))
  
-- | The number of bipartitions implied by a tree is one per EDGE in the tree.  Thus
-- each interior node carries a list of BiPs the same length as its list of children.
labelBips :: NewickTree a -> NewickTree (a, [DenseLabelSet])
labelBips tr =
    trace ("labelbips "++show allLeaves++" "++show size) $
    loop tr
  where    
    size = numLeaves tr
    zero = mkEmptyDense size
    loop (NTLeaf dec lab) = NTLeaf (dec, [markLabel lab zero]) lab      
    loop (NTInterior dec chlds) =
      let chlds' = map loop chlds
          sets   = map (normBip . denseUnions size . snd . get_dec) chlds' in
      NTInterior (dec, sets) chlds'

    halfSize = size `quot` 2
    normBip bip =
      let flipped = denseDiff allLeaves bip in
      case compare (bipSize bip) halfSize of
        LT -> bip 
        GT -> flipped -- Flip it
        EQ -> min bip flipped -- This is a painful case, we need a tie-breaker
           
    allLeaves = leafSet tr
    leafSet (NTLeaf _ lab)    = mkSingleDense size lab
    leafSet (NTInterior _ ls) = denseUnions size $ map leafSet ls


foldBips :: Monoid m => (DenseLabelSet -> m) -> NewickTree a -> m
foldBips f tr = F.foldMap f' (labelBips tr)
 where f' (_,bips) = F.foldMap f bips
  
-- | Get all non-singleton BiPs implied by a tree.
allBips :: NewickTree a -> S.Set DenseLabelSet
allBips tr = S.filter ((> 1) . bipSize) $ foldBips S.insert tr S.empty

--------------------------------------------------------------------------------
-- Optimized, LVish version
--------------------------------------------------------------------------------
-- First, necessary types:

#if 0
-- | A collection of all observed bipartitons (bips) with a mapping of which trees
-- contain which Bips.
type BipTable s = IMap DenseLabelSet s (SparseTreeSet s)
-- type BipTable = IMap BitList (U.Vector Bool)
-- type BipTable s = IMap BitList s (NA.NatArray s Word8)

-- | Sets of taxa (BiPs) that are expected to be sparse.
type SparseTreeSet s = IS.ISet s TreeID
-- TODO: make this a set of numeric tree IDs...
-- NA.NatArray s Word8

type TreeID = AnnotatedTree
-- | Tree's are identified simply by their order within the list of input trees.
-- type TreeID = Int
#endif
--------------------------------------------------------------------------------

-- The distance matrix is an atomically-bumped matrix of numbers.
-- type DistanceMat s = NA.NatArray s Word32
-- Except... bump isn't supported by our idempotent impl.

#if 0
-- | Returns a (square) distance matrix encoded as a vector.
distanceMatrix :: [AnnotatedTree] -> IO (U.Vector Word)
distanceMatrix lst = do 
--   IM.IMapSnap (table :: M.Map DenseLabelSet (S.Set TreeID)) <- runParThenFreeze par
--   IM.IMapSnap (table :: M.Map DenseLabelSet (Snapshot IS.ISet TreeID)) <- runParThenFreeze par
   IM.IMapSnap table <- runParThenFreeze par
   let sz = P.length lst
   v <- MV.replicate (sz*sz) (0::Word)
   let fn set () =
         
   F.foldrM 
   undefined
  
  -- runParThenFreeze -- get bip table
  -- followed by ... fill matrix from bip table  
  where
    par = do   
     table <- IM.newEmptyMap 
     forM_ lst (insertBips table)
     return table

insertBips :: BipTable s -> AnnotatedTree -> Par d s ()
insertBips table tree = do
    let bips = allBips tree
        fn bip () = do
          IM.modify table bip (IS.putInSet tree)
          return ()
    F.foldrM fn () bips 
#endif

--------------------------------------------------------------------------------

instance Pretty a => Pretty (S.Set a) where
 pPrint s = pPrint (S.toList s)
 

printDistMat :: Handle -> V.Vector (U.Vector Int) -> IO () 
printDistMat h mat = do
  hPutStrLn h "Robinson-Foulds distance (matrix format):"
  hPutStrLn h "-----------------------------------------"
  V.forM_ mat $ \row -> do 
    U.forM_ row $ \elem -> do
      hPutStr h (show elem)
      hPutStr h " "
    hPutStr h "0\n"          
  hPutStrLn h "-----------------------------------------"