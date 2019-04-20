module Scoring
  ( score
  ) where

import           Blosum          (blosum62)
import           Data.List       (sortOn)
import qualified Data.Map.Strict as Map (lookup)
import           Data.Maybe      (fromMaybe)
import           Data.Vector     (Vector, slice)
import qualified Data.Vector     as V (zip)
import           Model           (Gap, Protein (..))

fstGap :: Int
fstGap = -4

anotherGap :: Int
anotherGap = -1

-- TODO: Solve indices in scoreAlignmentBetween
score :: Protein -> Protein -> Int
score Protein {pSeq = s1, pGaps = g1} Protein {pSeq = s2, pGaps = g2} =
  go 0 0 . sortOn fst $ g1 ++ g2
  where
    minLength = min (length s1) (length s2)
    go :: Int -> Int -> [Gap] -> Int
    go i acc [] = acc + scoreAlignmentBetween i (minLength - 1) s1 s2
    go i acc (g@(gi, gl):gs) =
      let gapCost = scoreGap g
          scoreBeforeGap = scoreAlignmentBetween i (min minLength gi) s1 s2
       in go (gi + gl + 1) (acc + scoreBeforeGap + gapCost) gs

scoreAlignmentBetween :: Int -> Int -> Vector Char -> Vector Char -> Int
scoreAlignmentBetween i1 i2 s1 s2
  | i1 < i2 = sum (needlemanWunsch <$> V.zip slice1 slice2)
  | otherwise = 0
  where
    slice1 = slice i1 (i2 - i1) s1
    slice2 = slice i1 (i2 - i1) s2

scoreGap :: Gap -> Int
scoreGap (_, l) = fstGap + (l - 1) * anotherGap

needlemanWunsch :: (Char, Char) -> Int
needlemanWunsch k =
  fromMaybe (error $ "Key not in DB: " ++ show k) (Map.lookup k blosum62)
