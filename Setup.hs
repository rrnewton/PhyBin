#!/usr/bin/env runhaskell

import Distribution.Simple
import Distribution.Simple.PreProcess
import Distribution.PackageDescription	
import Distribution.Simple.LocalBuildInfo	
import System.Cmd(system) 
import System.Exit

main :: IO () 
main = do putStrLn$ "Running Setup.hs ..."
	  defaultMainWithHooks simpleUserHooks 
