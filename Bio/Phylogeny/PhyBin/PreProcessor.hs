{-# LANGUAGE ScopedTypeVariables #-}

-- | This is a preprocessor which can (optionally) be invoked to
-- collapse short branch lengths.

module Bio.Phylogeny.PhyBin.PreProcessor
       where

import  Bio.Phylogeny.PhyBin.CoreTypes (NewickTree(..), DefDecor, Label)

-- | Removes branches that do not meet a predicate, leaving a
--   shallower, "bushier" tree.
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
