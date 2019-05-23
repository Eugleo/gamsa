module Main where

import           Control.Monad        (mapM_)
import           Data.Random          (StdRandom (..), runRVar)
import           Input
import           Model                (Alignment (..))
import           MultipleSeqAlignment (run)
import           Utils                (fill)

main :: IO ()
main = runIO inp >>= print

runIO :: Alignment -> IO Alignment
runIO a = do
  result <- runRVar (run a) StdRandom
  mapM_ putStrLn . fill . aProteins $ result
  return result
