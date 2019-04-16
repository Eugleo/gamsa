module Main where

import           Input
import           Lib

main :: IO ()
main = runIO inp >>= print
