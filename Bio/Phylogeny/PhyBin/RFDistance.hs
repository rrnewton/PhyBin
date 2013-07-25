{-# LANGUAGE ScopedTypeVariables, CPP, BangPatterns #-}

module Bio.Phylogeny.PhyBin.RFDistance
       (
         -- * Types
         DenseLabelSet, DistanceMatrix,

         -- * Bipartition (Bip) utilities
         allBips, foldBips, dispBip,
         consensusTree, bipsToTree,

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
import           Data.Function       (on)
import           Data.Word
import qualified Data.Vector                 as V
import qualified Data.Vector.Mutable         as MV
import qualified Data.Vector.Unboxed.Mutable as MU
import qualified Data.Vector.Unboxed         as U
import qualified Data.Vector.Unboxed.Bit     as UB
import qualified Data.Bit                    as B
import           Text.PrettyPrint.HughesPJClass hiding (char, Style)
import           System.IO      (hPutStrLn, hPutStr, Handle)
import           System.IO.Unsafe

-- import           Control.LVish
-- import qualified Data.LVar.Set   as IS
-- import qualified Data.LVar.SLSet as SL

-- import           Data.LVar.Map   as IM
-- import           Data.LVar.NatArray as NA

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
                 (q,rem) = (diff1 + diff2) `quotRem` 2
#            ifdef NORMALIZATION          
             in q + rem
               -- if rem==0
               -- then q
               -- else -- trace "Warning, when dividing symmetric difference by two, there was a remainder!"
               --      q
#            else
             in diff1
#            endif
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

-- | This version slices the problem a different way.  A single pass over the trees
-- populates the table of bipartitions.  Then the table can be processed (locally) to
-- produce (non-localized) increments to a distance matrix.
hashRF :: Int -> [NewickTree a] -> DistanceMatrix
hashRF num_taxa trees = build M.empty (zip [0..] trees)
  where
    num_trees = length trees
    -- First build the table:
    build acc [] = ingest acc
    build acc ((ix,hd):tl) =
      let bips = allBips hd
          acc' = S.foldl' fn acc bips
          fn acc bip = M.alter fn2 bip acc
          fn2 (Just membs) = Just (markLabel ix membs)
          fn2 Nothing      = Just (mkSingleDense num_taxa ix)
      in      
      build acc' tl

    -- Second, ingest the table to construct the distance matrix:
    ingest :: M.Map DenseLabelSet DenseLabelSet -> DistanceMatrix
    ingest bipTable = runST theST
      where
       theST :: forall s0 . ST s0 DistanceMatrix
       theST = do 
        -- Triangular matrix, starting narrow and widening:
        matr <- MV.new num_trees
        -- Too bad MV.replicateM is insufficient.  It should pass index.  
        -- Instead we write this C-style:
        for_ (0,num_trees) $ \ ix -> do 
          row <- MU.replicate ix (0::Int)
          MV.write matr ix row
          return ()

        unsafeIOToST$ putStrLn$" Built matrix for dim "++show num_trees

        let bumpMatr i j | j < i     = incr i j
                         | otherwise = incr j i
            incr :: Int -> Int -> ST s0 ()
            incr i j = do -- Not concurrency safe yet:
--                          unsafeIOToST$ putStrLn$" Reading at position "++show(i,j)
                          row <- MV.read matr i
                          elm <- MU.read row j
                          MU.write row j (elm+1)
                          return ()
            fn bipMembs =
              -- Here we quadratically consider all pairs of trees and ask whether
              -- their edit distance is increased based on this particular BiP.
              -- Actually, as an optimization, it is sufficient to consider only the
              -- cartesian product of those that have and those that don't.
              let haveIt   = bipMembs
                  -- Depending on how invertDense is written, it could be useful to
                  -- fuse this in and deforest "dontHave".
                  dontHave = invertDense num_trees bipMembs
                  fn1 trId = traverseDense_ (fn2 trId) dontHave
                  fn2 trId1 trId2 = bumpMatr trId1 trId2
              in
--                 trace ("Computed donthave "++ show dontHave) $ 
                 traverseDense_ fn1 haveIt
        F.traverse_ fn bipTable
        v1 <- V.unsafeFreeze matr
        T.traverse (fmap (U.map (`quot` 2)) . U.unsafeFreeze) v1



#if 0
-- | Returns a (square) distance matrix encoded as a vector.
distanceMatrix :: [AnnotatedTree] -> IO (U.Vector Word)
distanceMatrix lst = do 
--   IM.IMapSnap (table :: M.Map DenseLabelSet (S.Set TreeID)) <- runParThenFreeze par
--   IM.IMapSnap (table :: M.Map DenseLabelSet (Snapshot IS.ISet TreeID)) <- runParThenFreeze par
   IM.IMapSnap table <- runParThenFreeze par
   let sz = P.length lst
   v <- MU.replicate (sz*sz) (0::Word)
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
  

-- | Take only the bipartitions that are agreed on by all trees.
consensusTree :: Int -> [NewickTree a] -> NewickTree ()
consensusTree _ [] = error "Cannot take the consensusTree of the empty list"
consensusTree num_taxa (hd:tl) = bipsToTree num_taxa intersection
  where
    intersection = L.foldl1' S.intersection (map allBips (hd:tl))
--     intersection = loop (allBips hd) tl
--     loop :: S.Set DenseLabelSet -> [NewickTree a] -> S.Set DenseLabelSet
--     loop !remain []      = remain
--     -- Was attempting to use foldBips here as an optimization:
-- --     loop !remain (hd:tl) = loop (foldBips S.delete hd remain) tl
--     loop !remain (hd:tl) = loop (S.difference remain (allBips hd)) tl    
      
-- | Convert from bipartitions BACK to a single tree.
bipsToTree :: Int -> S.Set DenseLabelSet -> NewickTree ()
bipsToTree num_taxa bip =
  loop lvl0 sorted
  where
    -- We consider each subset in increasing size order.
    -- FIXME: If we tweak the order on BIPs, then we can just use S.toAscList here:
    sorted = L.sortBy (compare `on` bipSize) (S.toList bip)

    lvl0 = [ (mkSingleDense num_taxa ix, NTLeaf () ix)
           | ix <- [0..num_taxa-1] ]

    -- We recursively glom together subtrees until we have a complete tree.
    -- We only process larger subtrees after we have processed all the smaller ones.
    loop !subtrees [] =
      case subtrees of
        []    -> error "bipsToTree: internal error"
        [(_,one)] -> one
        lst   -> NTInterior () (map snd lst)
    loop !subtrees (bip:tl) =
      let (in_,out) = L.partition ((denseIsSubset bip) . fst) subtrees in
      -- Here all subtrees that match the current bip get merged:
      loop ((denseUnions num_taxa (map fst in_),
             NTInterior ()        (map snd in_)) : out) tl

    -- sizeChop [] = []
    -- sizeChop (hd:tl) =
    --   let sz        = bipSize hd 
    --       (ths,rst) = span ((sz ==) . bipSize) tl
    --   in (hd:ths) : sizeChop rst

