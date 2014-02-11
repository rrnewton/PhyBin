{-# LANGUAGE ScopedTypeVariables, CPP, BangPatterns #-}
{-# LANGUAGE TypeFamilies #-}

module Bio.Phylogeny.PhyBin.RFDistance
       (
         -- * Types
         DenseLabelSet, DistanceMatrix,

         -- * Bipartition (Bip) utilities
         allBips, foldBips, dispBip,
         consensusTree, bipsToTree, filterCompatible, compatibleWith,

         -- * ADT for dense sets
         mkSingleDense, mkEmptyDense, bipSize,
         denseUnions, denseDiff, invertDense, markLabel,
         
        -- * Methods for computing distance matrices
        distanceMatrix, hashRF, 

        -- * Output
        printDistMat)
       where

import           Control.Monad
import           Control.Monad.ST
import           Control.Monad.IO.Class
import           Data.Function       (on)
import           Data.Word
import           Data.Time.Clock
import qualified Data.Vector                 as V
import qualified Data.Vector.Mutable         as MV

-- import qualified Data.Vector.Unboxed.Mutable as MU
-- import qualified Data.Vector.Unboxed         as U

-- TEMP: swapping these out:
import qualified Data.Vector.Storable.Mutable as MU
import qualified Data.Vector.Storable         as U

import qualified Data.Vector.Unboxed.Bit     as UB
import qualified Data.Vector.Storable.Mutable as MS
import qualified Data.Vector.Storable         as SV
import qualified Data.Bit                    as B
import qualified Foreign.ForeignPtr as FP
import           Foreign.Ptr        as Ptr
import           Foreign.Storable (sizeOf)
import qualified Data.Bits.Atomic   as BA
import           Text.PrettyPrint.HughesPJClass hiding (char, Style)
import           System.IO      (hPutStrLn, hPutStr, Handle)
import           System.IO.Unsafe

-- import qualified Control.Monad.Par.Combinator as PC

import           Control.LVish hiding (for_)
import           Control.LVish.DeepFrz (runParThenFreezeIO, Frzn)
#if 0
import qualified Data.LVar.Set   as IS
import           Data.LVar.Map   as IM
#else
#warning "Using skip-list based concurrent data structures for parallel PhyBin."
import qualified Data.LVar.SLSet as IS
import           Data.LVar.SLMap as IM
#endif

import           Data.LVar.NatArray as NA
import           Data.LVar.Counter (newSum, incrSum, freezeSum)

import           Bio.Phylogeny.PhyBin.CoreTypes
-- import           Data.BitList
import qualified Data.Set as S
import qualified Data.List as L
import qualified Data.IntSet as SI
import qualified Data.Map.Strict as M
import qualified Data.Foldable as F
import qualified Data.Traversable as T
import           Data.Monoid
import           Prelude as P
import           Debug.Trace

-- I don't understand WHY, but I seem to get the same answers WITHOUT this.
-- Normalization and symmetric difference do make things somewhat slower (e.g. 1.8
-- seconds vs. 2.2 seconds for 150 taxa / 100 trees)
#define NORMALIZATION
-- define BITVEC_BIPS

--------------------------------------------------------------------------------
-- A data structure choice
--------------------------------------------------------------------------------

-- type DenseLabelSet s = BitList


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
#ifdef BITVEC_BIPS

#  if 1
type DenseLabelSet = UB.Vector B.Bit
markLabel lab = UB.modify (\vec -> MU.write vec lab (B.fromBool True)) 
mkEmptyDense  size = UB.replicate size (B.fromBool False)
mkSingleDense size ind = markLabel ind (mkEmptyDense size)
denseUnions        = UB.unions
bipSize            = UB.countBits
denseDiff          = UB.difference
invertDense size bip = UB.invert bip
dispBip labs bip = show$ map (\(ix,_) -> (labs M.! ix)) $
                        filter (\(_,bit) -> B.toBool bit) $
                        zip [0..] (UB.toList bip)
denseIsSubset a b = UB.or (UB.difference b a)
traverseDense_ fn bip =
  U.ifoldr' (\ ix bit acc ->
              (if B.toBool bit
               then fn ix
               else return ()) >> acc)
        (return ()) bip

#  else
-- TODO: try tracking the size:
data DenseLabelSet = DLS {-# UNPACK #-} !Int (UB.Vector B.Bit)
markLabel lab (DLS _ vec)= DLS (UB.modify (\vec -> return (MU.write vec lab (B.fromBool True))) ) vec
-- ....
#  endif

#else
type DenseLabelSet = SI.IntSet
markLabel lab set   = SI.insert lab set 
mkEmptyDense _size  = SI.empty
mkSingleDense _size = SI.singleton
denseUnions _size   = SI.unions 
bipSize             = SI.size
denseDiff           = SI.difference
denseIsSubset       = SI.isSubsetOf

dispBip labs bip = "[" ++ unwords strs ++ "]"
  where strs = map (labs M.!) $ SI.toList bip
invertDense size bip = loop SI.empty (size-1)
  where -- There's nothing for it but to iterate and test for membership:
    loop !acc ix | ix < 0           = acc
                 | SI.member ix bip = loop acc (ix-1)
                 | otherwise        = loop (SI.insert ix acc) (ix-1)
traverseDense_ fn bip =
  -- FIXME: need guaranteed non-allocating way to do this.
  SI.foldr' (\ix acc ->  fn ix >> acc) (return ()) bip
#endif

markLabel    :: Label -> DenseLabelSet -> DenseLabelSet
mkEmptyDense :: Int -> DenseLabelSet
mkSingleDense :: Int -> Label -> DenseLabelSet
denseUnions  :: Int -> [DenseLabelSet] -> DenseLabelSet
bipSize      :: DenseLabelSet -> Int

-- | Print a BiPartition in a pretty form
dispBip      :: LabelTable -> DenseLabelSet -> String

-- | Assume that total taxa are 0..N-1, invert membership:
invertDense  :: Int -> DenseLabelSet -> DenseLabelSet

traverseDense_ :: Monad m => (Int -> m ()) -> DenseLabelSet -> m ()


--------------------------------------------------------------------------------
-- Dirt-simple reference implementation
--------------------------------------------------------------------------------

type DistanceMatrix = V.Vector (U.Vector Int)

-- | Returns a triangular distance matrix encoded as a vector.
--   Also return the set-of-BIPs representation for each tree.
distanceMatrix :: [NewickTree a] -> (DistanceMatrix, V.Vector (S.Set DenseLabelSet))
distanceMatrix lst = 
   let sz = P.length lst
       eachbips = V.fromList $ map allBips lst
       mat = V.generate sz $ \ i ->        
             U.generate i  $ \ j ->
             let diff1 = S.size (S.difference (eachbips V.! i) (eachbips V.! j))
                 diff2 = S.size (S.difference (eachbips V.! j) (eachbips V.! i))
             in diff1 + diff2
   in (mat, eachbips)

-- | The number of bipartitions implied by a tree is one per EDGE in the tree.  Thus
-- each interior node carries a list of BiPs the same length as its list of children.
labelBips :: NewickTree a -> NewickTree (a, [DenseLabelSet])
labelBips tr =
--    trace ("labelbips "++show allLeaves++" "++show size) $
#ifdef NORMALIZATION  
    fmap (\(a,ls) -> (a,map (normBip size) ls)) $
#endif
    loop tr
  where    
    size = numLeaves tr
    zero = mkEmptyDense size
    loop (NTLeaf dec lab) = NTLeaf (dec, [markLabel lab zero]) lab      
    loop (NTInterior dec chlds) =
      let chlds' = map loop chlds
          sets   = map (denseUnions size . snd . get_dec) chlds' in
      NTInterior (dec, sets) chlds'

    allLeaves = leafSet tr
    leafSet (NTLeaf _ lab)    = mkSingleDense size lab
    leafSet (NTInterior _ ls) = denseUnions size $ map leafSet ls

-- normBip :: DenseLabelSet -> DenseLabelSet -> DenseLabelSet
--    normBip allLeaves bip =
normBip :: Int -> DenseLabelSet -> DenseLabelSet    
normBip totsize bip =
  let -- size     = bipSize allLeaves
      halfSize = totsize `quot` 2
--      flipped  = denseDiff allLeaves bip
      flipped  = invertDense totsize bip 
  in 
  case compare (bipSize bip) halfSize of
    LT -> bip 
    GT -> flipped -- Flip it
    EQ -> min bip flipped -- This is a painful case, we need a tie-breaker
    

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

-- UNFINISHED:
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
-- Alternate way of slicing the problem: HashRF
--------------------------------------------------------------------------------

-- The distance matrix is an atomically-bumped matrix of numbers.
-- type DistanceMat s = NA.NatArray s Word32
-- Except... bump isn't supported by our idempotent impl.

type DistanceMatrix1 = V.Vector (SV.Vector Int)
type DistanceMatrix2 s = V.Vector (SV.Vector (Sum s))

-- | This version slices the problem a different way.  A single pass over the trees
-- populates the table of bipartitions.  Then the table can be processed (locally) to
-- produce (non-localized) increments to a distance matrix.
hashRF :: forall dec . Int -> [NewickTree dec] -> IO DistanceMatrix1
hashRF num_taxa trees = do
    t0  <- getCurrentTime
    bigtable <- getBigtable
    t1  <- getCurrentTime
    res <- runParNonDet $ ingest bigtable
    t2  <- getCurrentTime
    putStrLn$ "hashRF: time spent in first/second runPar and total: "
              ++show (diffUTCTime t1 t0, diffUTCTime t2 t1, diffUTCTime t2 t0)
    return res

  where
    getBigtable :: IO (IMap DenseLabelSet Frzn (IS.ISet Frzn Int))
    getBigtable = runParThenFreezeIO $ isQD par

    par :: forall e s . (HasPut e) => 
           Par e s (IMap DenseLabelSet s (IS.ISet s Int))
    par =  do themp <- newEmptyMap
--              build themp (zip [0..] trees)
              let domain = V.fromList (zip [0..] trees)
              parForTiled Nothing 16 (0,V.length domain) $ \ ix ->
                build themp (domain V.! ix)
              return themp

    num_trees = length trees
    -- First build the table:
    build :: (HasPut e) => 
             IMap DenseLabelSet s (IS.ISet s Int) -> (Int,NewickTree dec) -> Par e s ()
    build acc (ix,hd) = do
      let bips = allBips hd 
          fn bip = IM.gmodify acc bip (IS.insert ix)
      F.traverse_ fn bips

    -- Second, ingest the table to construct the distance matrix:
    ingest :: (HasBump e, HasIO e) =>
              IM.IMap DenseLabelSet Frzn (IS.ISet Frzn Int) -> Par e s (DistanceMatrix2 s)
    ingest bipTable = theST
      where
       theST = do 
        -- Triangular matrix, starting narrow and widening:
        matr <- liftIO$ MV.new num_trees
        -- Too bad MV.replicateM is insufficient.  It should pass index.  
        -- Instead we write this C-style:
        liftIO$ for_ (0,num_trees) $ \ ix -> do 
          row <- MS.replicate ix (0::Int)
          MV.write matr ix row
          return ()

        let bumpMatr i j | j < i     = incr i j
                         | otherwise = incr j i
            incr i j = do -- Not concurrency safe yet:
                          row <- MV.read matr i
                          incrStorable row j
                          return ()
            fn :: S.Set Int -> IO ()
            fn bipMembs =
              -- Here we quadratically consider all pairs of trees and ask whether
              -- their edit distance is increased based on this particular BiP.
              -- Actually, as an optimization, it is sufficient to consider only the
              -- cartesian product of those that have and those that don't.
              let haveIt   = bipMembs
                  -- Depending on how invertDense is written, it could be useful to
                  -- fuse this in and deforest "dontHave".
                  dontHave = invertDense2 num_trees bipMembs
                  fn1 trId = traverseDense_2 (fn2 trId) dontHave
                  fn2 trId1 trId2 = bumpMatr trId1 trId2
              in
--                 trace ("Computed donthave "++ show dontHave) $ 
                 traverseDense_2 fn1 haveIt
--        parForTiled 16 (0,M.size bipTable) $ \ix -> do

#if 0
        -- TODO: Need a proper parallell fold over the SLMap...
        PC.parFor (PC.InclusiveRange 0 (IM.size bipTable - 1)) $ \ix -> do
          -- liftIO$ F.traverse_ fn bipTable
          liftIO$ fn $ IS.fromISet (snd$ M.elemAt ix bipTable)
#else
-- TODO: Restore parallelism ^^:
        liftIO$ F.traverse_ (fn . IS.fromISet) bipTable
#endif

        liftIO$ do
          v1 <- V.unsafeFreeze matr
          T.traverse (SV.unsafeFreeze) v1




-- TEMPORARY:
--------------------------------------------------------------------------------
type DenseLabelSet2 = S.Set Int
markLabel2 lab set   = S.insert lab set 
mkEmptyDense2 _size  = S.empty
mkSingleDense2 _size = S.singleton
denseUnions2 _size   = S.unions 
bipSize2             = S.size
denseDiff2    a b    = S.difference a b
denseIsSubset2 a b   = S.isSubsetOf a b

dispBip2 labs bip = "[" ++ unwords strs ++ "]"
  where strs = map (labs M.!) $ S.toList bip
invertDense2 size bip = loop S.empty (size-1)
  where -- There's nothing for it but to iterate and test for membership:
    loop !acc ix | ix < 0           = acc
                 | S.member ix bip = loop acc (ix-1)
                 | otherwise        = loop (S.insert ix acc) (ix-1)
traverseDense_2 fn bip =
  -- FIXME: need guaranteed non-allocating way to do this.
  S.foldr' (\ix acc ->  fn ix >> acc) (return ()) bip

--------------------------------------------------------------------------------
-- Miscellaneous Helpers
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

-- My own forM for numeric ranges (not requiring deforestation optimizations).
-- Inclusive start, exclusive end.
{-# INLINE for_ #-}
for_ :: Monad m => (Int, Int) -> (Int -> m ()) -> m ()
for_ (start, end) _fn | start > end = error "for_: start is greater than end"
for_ (start, end) fn = loop start
  where
   loop !i | i == end  = return ()
           | otherwise = do fn i; loop (i+1)

-- | Which of a set of trees are compatible with a consensus?
filterCompatible :: NewickTree a -> [NewickTree b] -> [NewickTree b]
filterCompatible consensus trees =
    let cbips = allBips consensus in
    [ tr | tr <- trees
         , cbips `S.isSubsetOf` allBips tr ]

-- | Is a tree compatible with a consensus?
--   This is more efficient if partially applied then used repeatedly.
compatibleWith :: NewickTree a -> NewickTree b -> Bool
compatibleWith consensus newTr =
  S.isSubsetOf (allBips consensus) (allBips newTr)

-- Consensus between two trees, which may even have different label maps.
consensusTreeFull (FullTree n1 l1 t1) (FullTree n2 l2 t2) =
  error "FINISHME - consensusTreeFull"

-- | Take only the bipartitions that are agreed on by all trees.
consensusTree :: Int -> [NewickTree a] -> NewickTree ()
consensusTree _ [] = error "Cannot take the consensusTree of the empty list"
consensusTree num_taxa (hd:tl) = bipsToTree num_taxa intersection
  where
    intersection = L.foldl' S.intersection (allBips hd) (map allBips tl)
--     intersection = loop (allBips hd) tl
--     loop :: S.Set DenseLabelSet -> [NewickTree a] -> S.Set DenseLabelSet
--     loop !remain []      = remain
--     -- Was attempting to use foldBips here as an optimization:
-- --     loop !remain (hd:tl) = loop (foldBips S.delete hd remain) tl
--     loop !remain (hd:tl) = loop (S.difference remain (allBips hd)) tl    
      
-- | Convert from bipartitions BACK to a single tree.
bipsToTree :: Int -> S.Set DenseLabelSet -> NewickTree ()
bipsToTree num_taxa origbip =
--  trace ("Doing bips in order: "++show sorted++"\n") $ 
  loop lvl0 sorted
  where
    -- We consider each subset in increasing size order.
    -- FIXME: If we tweak the order on BIPs, then we can just use S.toAscList here:
    sorted = L.sortBy (compare `on` bipSize) (S.toList origbip)

    lvl0 = [ (mkSingleDense num_taxa ix, NTLeaf () ix)
           | ix <- [0..num_taxa-1] ]

    -- VERY expensive!  However, due to normalization issues this is necessary for now:
    -- TODO: in the future make it possible to definitively denormalize.
    -- isMatch bip x = denseIsSubset x bip || denseIsSubset x (invertDense num_taxa bip)
    isMatch bip x = denseIsSubset x bip 

    -- We recursively glom together subtrees until we have a complete tree.
    -- We only process larger subtrees after we have processed all the smaller ones.
    loop !subtrees [] =
      case subtrees of
        []    -> error "bipsToTree: internal error"
        [(_,one)] -> one
        lst   -> NTInterior () (map snd lst)
    loop !subtrees (bip:tl) =
--      trace (" -> looping, subtrees "++show subtrees) $ 
      let (in_,out) = L.partition (isMatch bip. fst) subtrees in
      case in_ of
        [] -> error $"bipsToTree: Internal error!  No match for bip: "++show bip
              ++" out is\n "++show out++"\n and remaining bips "++show (length tl)
              ++"\n when processing orig bip set:\n  "++show origbip
          -- loop out tl
        _ -> 
         -- Here all subtrees that match the current bip get merged:
         loop ((denseUnions num_taxa (map fst in_),
                NTInterior ()        (map snd in_)) : out) tl

{-# INLINE incrStorable #-}
incrStorable :: MS.MVector RealWorld Int -> Int -> IO ()
incrStorable (MS.MVector _ rowFP) offset =
  FP.withForeignPtr rowFP $ \ (ptr :: Ptr Int) -> do 
--    let ptr' = (P.ptrToIntPtr ptr) + offset
    let ptr' = ptr `plusPtr` (offset * sizeOf (0::Int))
    BA.fetchAndAdd ptr' (1::Int)
    return ()
  
