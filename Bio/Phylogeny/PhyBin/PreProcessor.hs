{-# LANGUAGE ScopedTypeVariables #-}

-- | This is a preprocessor which can (optionally) be invoked to
-- collapse short branch lengths.

module Bio.Phylogeny.PhyBin.PreProcessor
       ( collapseBranches,
         collapseBranchLenThresh, collapseBranchBootStrapThresh,
         pruneTreeLeaves
       )
       where

import qualified Data.Set as S
import Data.Maybe (catMaybes)
import  Bio.Phylogeny.PhyBin.CoreTypes 


-- | Prune the leaves of the tree to only those leaves in the provided set.
-- 
--   If ALL leaves are pruned from the set, this function returns nothing.
pruneTreeLeaves :: S.Set Label -> NewickTree a -> Maybe (NewickTree a)
pruneTreeLeaves set tr = loop tr
 where
   loop orig@(NTLeaf _ lab)
     | S.member lab set = Just orig
     | otherwise        = Nothing
   loop (NTInterior dec ls) =
     case catMaybes $ map loop ls of
       []    -> Nothing
       [one] -> Just one
       ls'   -> Just (NTInterior dec ls')


-- | Removes branches that do not meet a predicate, leaving a shallower, "bushier"
--   tree.  This does NOT change the set of leaves (taxa), it only removes interior
--   nodes.
--
--   `collapseBranches pred collapser tree` uses `pred` to test the meta-data to see
--   if collapsing the intermediate node below the branch is necessary, and if it is,
--   it uses `collapser` to reduce all the metadata for the collapsed branches into a
--   single piece of metadata.
collapseBranches :: forall a . (a -> Bool) -> (a -> a -> a) -> NewickTree a -> NewickTree a
collapseBranches isCollapsable collapse origtr = final
  where    
    (_,_, final) = loop origtr
    
    -- This loop returns:
    --  (1) a list of leaf "floaters" that can still move upwards,
    --  (2) immovable subtrees that can't
    --  (3) a final node IF the result is the root.
    loop :: NewickTree a -> ([(a,Label)], [NewickTree a], NewickTree a)
    loop lf@(NTLeaf dec lb) | isCollapsable dec = ([(dec,lb)], [],  lf)
                            | otherwise         = ([],        [lf], lf)
    loop (NTInterior dec children) =
      let (floats, anchors, _) = unzip3 $ map loop children 
          thenode = NTInterior dec $ concat anchors ++ 
                                     map (uncurry NTLeaf) (concat floats)
      in      
      if isCollapsable dec then
        -- If we are collapsable we keep floating upwards.
        -- We combine our metadata (on the left) into the floatees:
        (map (\ (a,l) -> (collapse dec a, l))
             (concat floats),
         concat anchors,
         thenode)
      else
        -- Otherwise this is the end of the road for these floaters:
        ([], [thenode], thenode)


-- | A common configuration.  Collapse branches based on a length threshold.
collapseBranchLenThresh :: Double -> NewickTree DefDecor -> NewickTree DefDecor
-- collapseBranchLenThresh :: HasBranchLen a => Double -> NewickTree a -> NewickTree a    
collapseBranchLenThresh thr tr =
  collapseBranches ((< thr) . getBranchLen) collapser tr
  where
    -- We REMOVE BootStraps as part of this process, they are not well-defined after this point.
    collapser _intermediate@(_bt1, len1) _floatee@(_bt2, len2) =
      (Nothing, len1 + len2)

-- | A common configuration.  Collapse branches based on bootstrap values.
collapseBranchBootStrapThresh :: Int -> NewickTree DefDecor -> NewickTree DefDecor
-- collapseBranchLenThresh :: HasBranchLen a => Double -> NewickTree a -> NewickTree a    
collapseBranchBootStrapThresh thr tr =
  collapseBranches ((< thr) . getBoot) collapser tr
  where
    getBoot (Nothing,_)    = error$"collapseBranchBootStrapThresh: bootstrap value missing on tree node!\n"++
                                   "They must be present if --minbootstrap is used."
    getBoot (Just boot,_)  = boot
    -- This had better happen BEFORE branch-length based collapsing is done:
    collapser (_,len1) (_,len2)  = (Nothing, len1+len2)



