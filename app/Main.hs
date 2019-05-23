module Main where

import           Control.Monad        (mapM_)
import           Data.Random          (StdRandom (..), runRVar)
import           Model                (Alignment (..))
import           MultipleSeqAlignment (run)
import           System.Environment
import           Utils                (fill, mkAlignment)

main :: IO ()
main = do
  args <- getArgs
  al <- mkAlignment . lines <$> readFile (head args)
  result <- runRVar (run al) StdRandom
  mapM_ putStrLn . fill . aProteins $ result
