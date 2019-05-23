module Crossover
  ( recombineV
  , recombineH
  ) where

import           Data.List   (sortOn)
import           Data.Random (RVar, stdUniform, uniform)
import           Model       (Alignment (..), Gap, Protein (..))
import           Scoring     (scoreProteins)

recombineH :: Alignment -> Alignment -> RVar Alignment
recombineH a@Alignment {aProteins = protA} Alignment {aProteins = protB}
  | protA == protB = return a
  | otherwise = do
    newProt <- mapM recombineProt (zip protA protB)
    return $ Alignment newProt (scoreProteins newProt)
  where
    recombineProt (a, Protein {pGaps = gapsB}) = do
      coin <- stdUniform
      if coin
        then return a
        else return $ Protein (pSeq a) gapsB (pMeanGapCount a)

recombineV :: Alignment -> Alignment -> RVar Alignment
recombineV a@Alignment {aProteins = protA} Alignment {aProteins = protB} = do
  i <- breakpoint
  newProt <- mapM (recombineProt i) (zip protA protB)
  return $ Alignment newProt (scoreProteins newProt)
  where
    seqLength s = length (pSeq s) -- + foldr (\g acc -> snd g + acc) 0 (pGaps s)
    minLength = minimum $ map seqLength protA
    breakpoint = uniform 1 (minLength - 1)
    recombineProt i (p1, p2) = do
      let (g1A, g2A) = splitGaps i (sortOn fst $ pGaps p1)
      let (g1B, g2B) = splitGaps i (sortOn fst $ pGaps p2)
      coin <- stdUniform
      if coin
        then return $ Protein (pSeq p1) (g1A ++ g2B) (pMeanGapCount p1)
        else return $ Protein (pSeq p1) (g1B ++ g2A) (pMeanGapCount p1)

splitGaps :: Int -> [Gap] -> ([Gap], [Gap])
splitGaps i = helper []
  where
    helper acc [] = (acc, [])
    helper acc (g@(ga, _):gs)
      | ga < i = helper (g : acc) gs
      | otherwise = (acc, gs)
