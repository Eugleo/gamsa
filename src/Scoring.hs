{-# LANGUAGE TupleSections #-}

module Scoring
  ( scorePair
  ) where

import           Blosum          (anotherGap, blosum62, fstGap)
import           Data.List       (sortOn)
import qualified Data.Map.Strict as Map (lookup)
import           Data.Maybe      (fromMaybe)
import           Data.Vector     (Vector, slice)
import qualified Data.Vector     as V (zip)
import           Debug.Trace     (traceShow)
import           Model           (Gap, Protein (..))

data ProteinID
  = A
  | B
  deriving (Eq)

scorePair :: Protein -> Protein -> Int
scorePair Protein {pSeq = s1, pGaps = g1} Protein {pSeq = s2, pGaps = g2} =
  go 0 0 0 0 . sortOn (fst . fst) $
  label A (helper 0 $ sortOn fst g1) ++ label B (helper 0 $ sortOn fst g2)
  where
    helper _ []            = []
    helper n ((gi, gl):gs) = (gi + n, gl) : helper (n + gl) gs
    label l = map (, l)
    go :: Int -> Int -> Int -> Int -> [(Gap, ProteinID)] -> Int
    go i adjA adjB acc [] =
      let sliceLength = minLength - i
          minLength = min (length s1 - adjA) (length s2 - adjB)
       in acc + scoreProteinSlice (i + adjA) (i + adjB) sliceLength s1 s2
    go i adjA adjB acc ((g@(gi, gl), name):gs) =
      let iA = i + adjA
          iB = i + adjB
          (newA, newB) =
            if name == A
              then (adjA - gl, adjB)
              else (adjA, adjB - gl)
          gapCost = scoreGap g
          minLength = min (length s1 - adjA) (length s2 - adjB)
          newIndex = max i (gi + gl)
          scoreBeforeGap = scoreProteinSlice iA iB (min gi minLength - i) s1 s2
       in if i > minLength
            then acc + gapCost
            else go newIndex newA newB (acc + gapCost + scoreBeforeGap) gs

scoreProteinSlice :: Int -> Int -> Int -> Vector Char -> Vector Char -> Int
scoreProteinSlice i1 i2 len s1 s2
  | len > 0 = sum $ needlemanWunsch <$> V.zip slice1 slice2
  | otherwise = 0
  where
    slice1 = slice i1 len s1
    slice2 = slice i2 len s2

scoreGap :: Gap -> Int
scoreGap (_, l) = fstGap + (l - 1) * anotherGap

needlemanWunsch :: (Char, Char) -> Int
needlemanWunsch k =
  fromMaybe (error $ "Key not in DB: " ++ show k) (Map.lookup k blosum62)
