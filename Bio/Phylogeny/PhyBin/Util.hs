{-# LANGUAGE ScopedTypeVariables #-}
-- RecordWildCards, TypeSynonymInstances, CPP
-- {-# LANGUAGE NamedFieldPuns #-}
-- {-# LANGUAGE OverloadedStrings #-}
-- {-# OPTIONS_GHC -fwarn-incomplete-patterns #-}
-- {-# OPTIONS_GHC -fwarn-unused-imports #-}

-- | This module contains misc bits used by (multiple) other modules.

module Bio.Phylogeny.PhyBin.Util
       ( 
         is_regular_file, acquireTreeFiles,
         safePrintDendro, sanityCheck
       )
       where

import qualified Data.Foldable as F
import           Data.Function       (on)
import           Data.List           (delete, minimumBy, sortBy, insertBy, intersperse, sort)
import           Data.Maybe          (fromMaybe, catMaybes)
import qualified Data.Map                   as M
import qualified Data.Set                   as S
import           Control.Monad       (forM, forM_, filterM, when, unless)
import           Control.Exception   (evaluate)
import           Control.Applicative ((<$>),(<*>))
import           Control.Concurrent  (Chan)
import           System.FilePath     (combine)
import           System.Directory    (doesFileExist, doesDirectoryExist,
                                      getDirectoryContents, getCurrentDirectory)
import           System.IO           (openFile, hClose, IOMode(ReadMode), stderr,
                                      hPutStr, hPutStrLn)
import           System.Exit         (ExitCode(..))
import           System.Timeout      (timeout)
import           Test.HUnit          ((~:),(~=?),Test,test)

-- For vizualization:
import           Text.PrettyPrint.HughesPJClass hiding (char, Style)
import           Bio.Phylogeny.PhyBin.CoreTypes
import           Bio.Phylogeny.PhyBin.Parser (parseNewick)
import           Bio.Phylogeny.PhyBin.Visualize (dotToPDF, dotNewickTree, viewNewickTree)
import           Bio.Phylogeny.PhyBin.RFDistance

import qualified Data.Clustering.Hierarchical as C
import qualified Data.Graph.Inductive as G
import qualified Data.GraphViz        as Gv
import  Data.GraphViz.Types.Canonical (nodeStmts, graphStatements)

----------------------------------------------------------------------------------------------------
-- OS specific bits:
----------------------------------------------------------------------------------------------------
-- #ifdef WIN32
-- is_regular_file = undefined
-- is_directory path = 
--   getFileAttributes
-- --getFileInformationByHandle
-- --    bhfiFileAttributes
-- file_exists = undefined
-- #else
-- is_regular_file :: FilePath -> IO Bool
-- is_regular_file file = 
--   do stat <- getFileStatus file; 
--      -- Hmm, this is probably bad practice... hard to know its exhaustive:
--      return$ isRegularFile stat || isNamedPipe stat || isSymbolicLink stat
-- is_directory :: FilePath -> IO Bool
-- is_directory path = 
--   do stat <- getFileStatus path
--      return (isDirectory stat)
-- file_exists = fileExist
-- #endif

-- Here we ASSUME it exists, then these functions are good enough:
is_regular_file :: FilePath -> IO Bool
is_regular_file = doesFileExist

is_directory :: FilePath -> IO Bool
is_directory = doesDirectoryExist 

file_exists :: FilePath -> IO Bool
file_exists path = 
  do f <- doesFileExist path
     d <- doesDirectoryExist path
     return (f || d)

--------------------------------------------------------------------------------

-- | Expand out directories to find all the tree files.
acquireTreeFiles :: [String] -> IO [String]
acquireTreeFiles inputs = do 
    all :: [[String]] <- forM inputs $ \ path -> do
      exists <- file_exists path 

      --stat   <- if exists then getFileStatus path else return (error "internal invariant")
      -- [2010.09.23] This is no longer really necessary:
      if not exists then do
        error$ "No file or directory found at this path!: "++path
	 -- hPutStr stderr$ "Input not a file/directory, assuming wildcard, using 'find' for expansion"
	 -- entries <- HSH.run$ "find " ++ path	 
	 -- hPutStrLn stderr$ "("++show (length entries)++" files found):  "++ show path
	 -- return entries
       else do
	 isdir <- is_directory path
	 reg  <- is_regular_file path
	 if isdir then do 
	    hPutStr stderr$ "Input is a directory, reading all regular files contained "
	    children <- getDirectoryContents path
	    filtered <- filterM is_regular_file $ map (combine path) children
	    hPutStrLn stderr$ "("++show (length filtered)++" regular files found):  "++ show path
	    return$ filtered
          else if reg then do 
	    return [path]
	  else error$ "phybin: Unhandled input path: " ++ path

    return (concat all)

--------------------------------------------------------------------------------

-- Detect cycles:
-- safePrintDendro :: C.Dendrogram (FullTree a) -> IO (Maybe String)
safePrintDendro :: Gv.DotGraph G.Node -> IO (Maybe String)
safePrintDendro dotg= do 
--  putStrLn$ "Dendrogram graph size: "++ show (F.foldl' (\a _ -> a+1) 0 dotg)
  mx <- timeout (2 * 1000 * 1000) $ do
--        putStrLn$ "Dendrogram graph, is directed?: "++ show (Gv.directedGraph dotg)
        putStrLn$ "Dendrogram graph size: "++ show (length $ nodeStmts $ graphStatements dotg)
        let str = show dotg
        evaluate (length str)
        return str
  case mx of
    Nothing -> do putStrLn "WARNING: DotGraph appears to be a cyclic structure.  This is probably a bug."
                  return Nothing
    _ -> return mx

sanityCheck :: C.Dendrogram (FullTree DefDecor) -> IO ()
sanityCheck dendro = do 
  let fn seen elm | S.member (treename elm) seen =
                       error$"Dendrogram failed sanity check!  Tree name occurs multiple times: "++(treename elm)
                  | otherwise = S.insert (treename elm) seen
      sz = S.size $ F.foldl' fn S.empty dendro
  putStrLn$ "Sanity checked dendrogram of size: "++show sz

