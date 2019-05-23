{-# LANGUAGE NamedFieldPuns #-}

module MutationProbabilities
  ( nextGenProbabilities
  , Probabilities(..)
  , Stats(..)
  , Mutations
  , pick
  ) where

import Utils (between)

-- | Popisuje současný stav mutačních operací. Během běhu generací jsou
-- sbírány statistiky o průběhu různých druhů mutací, které jsou dále využívány k dynamickému
-- přizpůsobování toho, s jakou pravděpodobností bude vybrán konkrétní druh mutace (ciz níže).
-- Tyto pravděpodobnosti se mění pouze napříč populacemi, v rámci jedné generace jsou konstantní.
type Mutations = (Probabilities, Stats)

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
  { sis :: (Int, Int) -- ^ insert
  , sic :: (Int, Int) -- ^ increase
  , sdc :: (Int, Int) -- ^ decrease
  , sdl :: (Int, Int) -- ^ delete
  , shf :: (Int, Int) -- ^ shift
  }

-- | Tato funkce je volána vždy po dokončení jedné generace; zpracuje nasbíraná data
-- a na základě nich spočítá nové hodnoty pravděpodobností jendotlivých operací.
-- Jsou používány komlikované vzorce z originálního paperu (DOI: 10.1155/2009/963150),
-- víceméně se jedná o to, že operace s pozitivním dopadem budou mít vyšší pravděpodobnost.
nextGenProbabilities :: Mutations -> Probabilities
nextGenProbabilities (P {pis, pic, pdc, pdl, phf}, S {sis, sic, sdc, sdl, shf}) =
  let fin = map (/ newProbSum) newProbList -- chceme, aby finální pravděpodobnosti měly součet 1.0
   in P (fin !! 0) (fin !! 1) (fin !! 2) (fin !! 3) (fin !! 4)
  where
    probList = [pis, pic, pdc, pdl, phf] -- se seznamy se lépe pracuje
    statList = [sis, sic, sdc, sdl, shf]
    mutationsCount = sum . map fst $ statList -- suma přes počet použití všech mutačních operací
    tDSOs = map (tDSO mutationsCount) statList
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
pick :: Probabilities -> Double -> Int
pick P {pis, pic, pdc, pdl, phf} num = go 0 0 pis [pis, pic, pdc, pdl, phf]
  where
    go n _ _ [] = n - 1
    go n l r (p:ps)
      | num `between` (l, r) = n
      | otherwise = go (n + 1) r (r + p) ps
