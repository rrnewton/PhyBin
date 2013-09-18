{-# LANGUAGE TemplateHaskell #-}
{-# LANGUAGE ScopedTypeVariables, BangPatterns, ParallelListComp, OverloadedStrings  #-}
module Main  where

import Data.Vector         as V
import Data.Vector.Unboxed as U
import Data.ByteString.Char8 as B
import Bio.Phylogeny.PhyBin
import Bio.Phylogeny.PhyBin.CoreTypes
import Bio.Phylogeny.PhyBin.RFDistance
import Bio.Phylogeny.PhyBin.Parser
import System.IO
import Prelude as P

import qualified Test.HUnit as HU
-- import Test.Framework
import Test.Framework.Providers.HUnit
import Test.Framework.TH (defaultMainGenerator)

-- Here we compare a tree against the same tree with 18 removed:
-- Or the same tree with one intermediate node removed.
(_,t30) = parseNewicks id $
          [ ("t30_0", "(2_, (((14, 3_), (19, (5_, 13))), (18, (6_, 7_))), 1_);")
          , ("t30_1", "(2_, (((14, 3_), (19, (5_, 13))), (6_, 7_)), 1_);")
          , ("t30_2", "(2_, (((14, 3_), (19, (5_, 13))), (18, 6_, 7_)), 1_);")              
          ]

t30' = naiveDistMatrix$ P.map nwtree t30
 -- "[[],[1]]"

case_t30 = HU.assertEqual "3-tree distance matrix"
              "[[],[0],[1,0]]" (showMat t30')

-- m = printDistMat stdout (fst d)


-- Simple show for a distance matrix:
showMat (m,_) = show$ V.toList$ V.map U.toList m

main = $(defaultMainGenerator)

