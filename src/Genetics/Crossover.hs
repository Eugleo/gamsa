module Genetics.Crossover
  ( recombineV
  , recombineH
  ) where

import Data.List        (partition)
import Data.Random      (RVar, stdUniform, uniform)

import Genetics.Scoring (scoreProteins)
import Model            (Alignment (..), Protein (..))

-- | Udělá horizontální crossover dvou alignmentů (A1 a A2), které se (mezerami) liší,
-- ale jejichž proteiny jsou stejné. Vytvoří nový alignment,
-- v němž budou v každém jednotlivém proteinu mezery náhodně buďto z A1, nebo A2.
recombineH :: Alignment -> Alignment -> RVar Alignment
recombineH al@Alignment {aProteins = proteinsA} Alignment {aProteins = proteinsB}
  | proteinsA == proteinsB = return al -- nemusíme nic kombinovat, alignmenty jsou stejné
  | otherwise = do
    newProteins <- mapM recombineH' (zip proteinsA proteinsB)
    return $ Alignment newProteins (scoreProteins newProteins)

-- | Pomocná funkce k recombineH, provádí samotnou rekombinaci dvou proteinů.
recombineH' :: (Protein, Protein) -> RVar Protein
recombineH' (p1, p2) = do
  coin <- stdUniform
  if coin
    -- p1 i p2 mají stejné hodnoty pSeq i pMeanGapCount, liší se pouze v mezerách
    then return p1
    else return p2

-- | Udělá vertikální crossover dvou alignmentů (A1 a A2), které se (mezerami) liší,
-- ale jejichž proteiny jsou stejné. Vytvoří nový alignment,
-- v němž budou u všech proteinů mezery před k-tou aminokyselinou z A1 (A2) a za k-tou
-- aminokyselinou budou všechny mezery z A2 (A1).
recombineV :: Alignment -> Alignment -> RVar Alignment
recombineV al@Alignment {aProteins = proteinsA} Alignment {aProteins = proteinsB}
  | proteinsA == proteinsB = return al -- nemusíme nic kombinovat, alignmenty jsou stejné
  | otherwise = do
    breakpoint <- uniform 1 (minLength - 1) -- generujeme index AK, podle které proteiny rozdělíme
    newProteins <- mapM (recombineV' breakpoint) (zip proteinsA proteinsB)
    return $ Alignment newProteins (scoreProteins newProteins)
  where
    minLength = minimum $ map (length . pSeq) proteinsA -- délka nejkratšího proteinu (bez mezer)

-- | Pomocná funkce k recombineV, provádí samotnou rekombinaci dvou proteinů.
recombineV' :: Int -> (Protein, Protein) -> RVar Protein
recombineV' i (Protein {pSeq = seq, pGaps = gaps1, pMeanGapCount = mgc}, Protein {pGaps = gaps2}) = do
  let (g1, g2) = partition ((> i) . fst) gaps1
  let (h1, h2) = partition ((> i) . fst) gaps2
  coin <- stdUniform
  -- oba proteiny mají stejné hodnoty pSeq i pMeanGapCount, liší se pouze v mezerách
  return $
    Protein
      seq
      (if coin
         then g1 ++ h2
         else h1 ++ g2)
      mgc
