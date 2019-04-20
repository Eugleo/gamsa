module Main where

import           Input
import           MultipleSeqAlignment

main :: IO ()
main = runIO inp >>= print
