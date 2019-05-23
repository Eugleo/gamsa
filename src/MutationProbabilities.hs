{-# LANGUAGE NamedFieldPuns #-}

module MutationProbabilities
  ( nextGenProbabilities
  , Probabilities(..)
  , Stats(..)
  , MutationState
  , pick
  , listToS
  ) where

import Control.Monad.Trans.Class      (lift)
import Control.Monad.Trans.State.Lazy (StateT (..), get)
import Data.Random                    (RVar, stdUniform)

import Utils                          (between)

-- | Popisuje současný stav mutačních operací. Během běhu generací jsou
-- sbírány statistiky o průběhu různých druhů mutací, které jsou dále využívány k dynamickému
-- přizpůsobování toho, s jakou pravděpodobností bude vybrán konkrétní druh mutace (ciz níže).
-- Tyto pravděpodobnosti se mění pouze napříč populacemi, v rámci jedné generace jsou konstantní.
type MutationState = (Probabilities, Stats)

-- | Pravděpodobnosti mutačních operacích, jejich součet by měl vždy být 1.0
data Probabilities = P
  { pis :: Double -- ^ insert
  , pic :: Double -- ^ increase
  , pdc :: Double -- ^ decrease
  , pdl :: Double -- ^ delete
  , phf :: Double -- ^ shift
  }

-- | Informace o mutačních operacích; pro každou z nich je zde dvojice (A, B), kde
-- A udává, kolikrát byla operace v současné generaci použita ("total")
-- B je součet rozdílů skóre alignmentu před provedením operace a po ní ("dif")
data Stats = S
  { sis :: [Double] -- ^ insert
  , sic :: [Double] -- ^ increase
  , sdc :: [Double] -- ^ decrease
  , sdl :: [Double] -- ^ delete
  , shf :: [Double] -- ^ shift
  }

listToS :: [[Double]] -> Stats
listToS [a, b, c, d, e] = S a b c d e
listToS _               = error "Wrong number of list items"

-- | Tato funkce je volána vždy po dokončení jedné generace; zpracuje nasbíraná data
-- a na základě nich spočítá nové hodnoty pravděpodobností jendotlivých operací.
-- Jsou používány komlikované vzorce z originálního paperu (DOI: 10.1155/2009/963150),
-- víceméně se jedná o to, že operace s pozitivním dopadem budou mít vyšší pravděpodobnost.
nextGenProbabilities :: MutationState -> Probabilities
nextGenProbabilities (P {pis, pic, pdc, pdl, phf}, S {sis, sic, sdc, sdl, shf}) =
  let fin = map (/ newProbSum) newProbList -- chceme, aby finální pravděpodobnosti měly součet 1.0
   in P (fin !! 0) (fin !! 1) (fin !! 2) (fin !! 3) (fin !! 4)
  where
    probList = [pis, pic, pdc, pdl, phf] -- se seznamy se lépe pracuje
    statList = [sis, sic, sdc, sdl, shf]
    tDSOs = map sum statList
    normalizedTDSOs = map (/ maximum tDSOs) tDSOs
    newProbList = zipWith (\p t -> p + p * 0.1 * t) probList normalizedTDSOs -- 0.1 je zvolená konstanta
    newProbSum = sum newProbList

-- | Total attributed difference by suboperator, udává, jak pozitivní/negativní dopad měl daný operátor
-- vzhledem k tomu, jak často byl používán v poměru s ostatními.
tDSO :: Int -> (Int, Int) -> Double
tDSO mutationsCount (total, dif) =
  fromIntegral dif * fromIntegral total / fromIntegral mutationsCount

-- | Dostane číslo od 0 do 1 a na základě něj vybere s danou pravděpodobností jednu z mutačních
-- operací (respektive, vrátí její index, protože k samotné funkci nemá přístup)
pick :: StateT MutationState RVar Int
pick = do
  (P {pis, pic, pdc, pdl, phf}, _) <- get
  num <- lift stdUniform
  return $ go num 0 0 pis [pis, pic, pdc, pdl, phf] -- num
  where
    go num n _ _ [] = n - 1
    go num n l r (p:ps)
      | num `between` (l, r) = n
      | otherwise = go num (n + 1) r (r + p) ps
