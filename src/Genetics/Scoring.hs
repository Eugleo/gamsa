{-# LANGUAGE TupleSections #-}

module Genetics.Scoring
  ( scoreProteins
  ) where

import           Data.List             (sortOn, tails)
import qualified Data.Map.Strict       as Map (lookup)
import           Data.Maybe            (fromMaybe)
import           Data.Vector           (Vector, slice)
import qualified Data.Vector           as V (zip)

import           Genetics.ScoringTable (anotherGap, blosum62, fstGap)
import           Model                 (Gap, Protein (..))

data ProteinID
  = A
  | B
  deriving (Eq)

scoreProteins :: [Protein] -> Int
scoreProteins al = sum [scorePair a b | (a:rest) <- tails al, b <- rest]

scorePair :: Protein -> Protein -> Int
scorePair Protein {pSeq = s1, pGaps = g1} Protein {pSeq = s2, pGaps = g2} =
  go 0 0 0 0 . sortOn (fst . fst) $
  label A (helper 0 $ sortOn fst g1) ++ label B (helper 0 $ sortOn fst g2)
  where
    helper _ []            = []
    helper n ((gi, gl):gs) = (gi + n, gl) : helper (n + gl) gs
    label l = map (, l)
    minLength = min (length s1 + gapLength g1) (length s2 + gapLength g2)
    gapLength = sum . map snd
    go :: Int -> Int -> Int -> Int -> [(Gap, ProteinID)] -> Int
    -- No gaps left, score the rest of the sequences
    go i a b acc [] =
      acc + scoreProteinSlice (i + a) (i + b) (minLength - i) s1 s2
    -- There's a gap g remaining
    go i a b acc ((g@(gi, gl), pID):gs)
      -- We're beyond the end of the shorter sequence
      -- That means we can now score the remaining gaps, which are all in the longer seq
      | i >= minLength =
        foldr (\(gp, _) part -> part + scoreGap gp) (acc + scoreGap g) gs
      | otherwise =
        let (newA, newB) =
              if pID == A
                then (a - gl, b)
                else (a, b - gl)
            newIndex = max i (gi + gl)
            scoreBeforeGap =
              scoreProteinSlice (i + a) (i + b) (min gi minLength - i) s1 s2
         in go newIndex newA newB (acc + scoreGap g + scoreBeforeGap) gs

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
