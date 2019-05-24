{-# LANGUAGE TupleSections #-}

module Genetics.Scoring
  ( scoreProteins
  ) where

import           Data.List              (sortOn, tails)
import qualified Data.Map.Strict        as Map (lookup)
import           Data.Maybe             (fromMaybe)
import           Data.Vector            (Vector, slice)
import qualified Data.Vector            as V (zip)

import           Debug.Trace            (traceShow)

import           Genetics.ScoringMatrix (anotherGap, blosum62, fstGap)
import           Model                  (Gap, Protein (..))

data ProteinID
  = A
  | B
  deriving (Eq)

-- | Oskóruje seznam proteinů pomocí "sum of pairs" metody. Spočítá skóre každého páru a výsledky sečte.
scoreProteins :: [Protein] -> Int
scoreProteins prots =
  sum
    [ scorePair maxLength a b
    | (a:rest) <- tails prots -- chceme oskórovat každý pár jen jednou
    , b <- rest -- b /= a, b je vždy v seznamu proteinů až po a
    ]
  where
    maxLength =
      maximum . map (\p -> length (pSeq p) + sum (map snd $ pGaps p)) $ prots

-- | Oskóruje dvojici proteinů podle tabulky BLOSUM62.
scorePair :: Int -> Protein -> Protein -> Int
scorePair maxLen Protein {pSeq = s1, pGaps = g1} Protein {pSeq = s2, pGaps = g2} =
  scorePairGo (s1, s2, minLength, maxLen) 0 0 0 0 . sortOn (fst . fst) $
  -- sloučíme mezery z obou proteinů, ale předtím si označíme, v jakém proteinu byly
  -- kromě toho je také seřadíme a převedeme na alternativní systém souřadnic
  label A (helper 0 $ sortOn fst g1) ++ label B (helper 0 $ sortOn fst g2)
  where
    label l = map (, l)
    -- převedeme souřadnice mezer z (číslo AK, délka) na (místo v sekvenci, délka)
    -- neboli souřadnice mezer dále v sekvenci bude ovlivněna tím, kolik mezer je před nimi
    -- mezery musí být seřazené podle souřadnic
    helper _ []            = []
    helper n ((gi, gl):gs) = (gi + n, gl) : helper (n + gl) gs
    -- délka nejkratší proteinu (včetně délky mezer)
    minLength = min (length s1 + gapLength g1) (length s2 + gapLength g2)
    gapLength = sum . map snd -- sečteme délky všech mezer

-- | Pomocná funkce scorePair, dělá vlastně to samé, ale provádí už samotnou rekurzi. Postupně
-- prochází postupně všemi mezerami v alignmentu, skóruje úseky mezi nimi, a k tomu vždy přičítá
-- skóre úseku s mezerou. Dává pozor na to, abychom omylem neskórovali úseky, kde mají mezeru oba proteiny.
scorePairGo ::
     ( Vector Char -- ^ sekvence prvního proteinu
     , Vector Char -- ^ sekvence druhého proteinu
     , Int -- ^ délka nejkratší sekvence (včetně mezer)
     , Int -- ^ délka celého alignmentu
      )
  -> Int -- ^ index současné aminokyseliny vzhledem k alignmentu
  -> Int -- ^ korekce indexu současné aminokyseliny vzhledem k prvnímu proteinu
  -> Int -- ^ korekce indexu současné aminokyseliny vzhledem k druhému proteinu
  -> Int -- ^ akumulátor skóre
  -> [(Gap, ProteinID)] -- ^ seznam zbývajících mezer
  -> Int
scorePairGo (s1, s2, minLength, maxLength) i a b acc [] -- nezbývají mezery, stačí nám oskórovat zbytek sekvencí
 =
  acc + scoreProteinSlice (i + a) (i + b) (minLength - i) s1 s2 +
  scoreGap (0, maxLength - minLength)
scorePairGo t@(s1, s2, minLength, maxLength) i a b acc ((g@(gi, gl), pID):gs)
  | i >= minLength -- prošli jsme celý kratší protein, stačí nám oskórovat zbylé mezery v delším proteinu
   =
    scoreGap (0, maxLength - i) +
    foldr (\(gp, _) part -> part + scoreGap gp) (acc + scoreGap g) gs
  | otherwise =
    let (newA, newB) =
          if pID == A -- updatujeme korekce
            then (a - gl, b) -- index současné AK v A nyní bude o gl posunutý, nebo
            else (a, b - gl) -- index současné AK v B nyní bude o gl posunutý
        newIndex = max i (gi + gl) -- updatujeme součsný index
        -- oskórujeme proteiny až do místa začátku současné mezery
        scoreBeforeGap =
          scoreProteinSlice (i + a) (i + b) (min gi minLength - i) s1 s2
     in scorePairGo t newIndex newA newB (acc + scoreGap g + scoreBeforeGap) gs

-- | Oskóruje úsek dvou proteinů. Předpokládá, že na skórovaném úseku nejsou
-- ani v jednom proteinu žádné mezery.
scoreProteinSlice ::
     Int -- ^ index počáteční aminokyseliny v prvním proteinu
  -> Int -- ^ index počáteční aminokyseliny v druhém proteinu
  -> Int -- ^ délka úseku
  -> Vector Char -- ^ sekvence prvního proteinu
  -> Vector Char -- ^ sekvence druhého proteinu
  -> Int -- ^ skóre
scoreProteinSlice i1 i2 len s1 s2
  | len > 0 = sum $ scoreAA <$> V.zip slice1 slice2
  | otherwise = 0
  where
    slice1 = slice i1 len s1
    slice2 = slice i2 len s2

-- | Oskóruje mezeru, podle hodnot udaných ve ScoringMatrix.
scoreGap :: Gap -> Int
scoreGap (_, l) = fstGap + (l - 1) * anotherGap

-- | Oskóruje pár aminokyselin, podle tabulky BLOSUM62.
scoreAA :: (Char, Char) -> Int
scoreAA k = fromMaybe (error "Key not in DB") (Map.lookup k blosum62)
