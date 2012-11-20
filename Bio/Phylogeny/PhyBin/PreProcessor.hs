

-- | This is a preprocessor which can (optionally) be invoked to
-- collapse short branch lengths.

module Bio.Phylogeny.PhyBin.PreProcessor
       where

import  Bio.Phylogeny.PhyBin.CoreTypes (NewickTree(..), DefDecor)

-- | Removes branches with length less than the threshold.
collapseBranches :: Double -> NewickTree DefDecor -> NewickTree DefDecor
collapseBranches thresh tr =
  tr

