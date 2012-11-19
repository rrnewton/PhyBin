

-- | This is a quick hack to split up a table in a tab-delimited file of
-- the format Irene produces.  The result is a simple two-column,
-- white-space delimited table of the kind phybin expects.


import System.Environment (getArgs)

import Control.Monad
import Data.List
import Data.List.Split

import System.IO


isDataLine :: String -> Bool
isDataLine l =
  case words l of
    ("Roundup":"Orthology":_) -> False
    (_:"results":"found":_)   -> False
    ("Gene":"Cluster":_) -> False
    ("Id":"Genome":_)    -> False
    []                   -> False
    _                    -> True

main = do
  args <- getArgs
  let file = case args of
              [f] -> f
              _   -> error "Expects one argument!! [filename]"

  raw <- readFile file
  let lns  = lines raw
      filt = filter isDataLine lns
      toks = map (splitOn ['\t']) filt
      put  = hPutStrLn stderr
  
  put$"Read "++show (length lns)++" lines from file "++show file
  put$"     "++show (length filt)++" contain data"
  put$"  Distinct #toks found on data lines: " ++ show(nub (map length toks))
  put$"  A sample of ten parsed lines :"
  mapM_ (put .  ("    "++) . show) $ take 10 $ toks

  put$"  Echoing columns one and two in space-separated form:"

  forM_ toks $ \ (one:two:_) ->
    putStrLn$ one ++"  "++ two

  return ()
