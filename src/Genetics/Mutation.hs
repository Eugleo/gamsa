{-# LANGUAGE NamedFieldPuns #-}

module Genetics.Mutation
  ( mutate
  ) where

import Control.Monad                        (replicateM)
import Control.Monad.Trans.Class            (lift)
import Control.Monad.Trans.State.Lazy       (StateT (..), get, modify, put,
                                             runStateT)
import Data.List                            (partition)
import Data.Random                          (RVar, stdUniform)
import Data.Random.Distribution.Exponential (exponential)

import Genetics.Scoring                     (scoreProteins)
import Model
import MutationProbabilities                (MutationState, Stats (..), listToS,
                                             pickOperationIndex)
import Utils                                (choose, chooseI, updateAt)

-- | Pomocný typ popsiující funkci, která na základě nějakého MutationState
-- budě umět zmutovat Alignment
type Mutator = StateT MutationState RVar

-- | Vratí Mutator, který bude n-krát náhodně mutovat alignment,
-- přičemž průměr přes všechna n (pro všechny alignmenty) bude 1.8 (zvolená konstanta)
mutate :: Alignment -> Mutator Alignment
mutate al@Alignment {aScore} = do
  count <- lift (exponential 1.8 :: RVar Float) -- count může být i 0
  (opIndices, newAl) <- runStateT (replicateM (round count) mutateOnce) al
  let dif = scoreProteins (aProteins newAl) - aScore -- rozdíl mezi skóre před mutacemi a po nich
  modify $ updateMutationUseCounts opIndices dif -- musíme výsledky provedených mutací uložit
  return newAl

-- | Poměry jednotlivých mutačních operací se adaptují tomu, jak "úspěšné" jednotlivé oprerace jsou.
-- Tato funkce pro provedeném souboru mutací spočítá jejich úspěšnost (dSO, viz níže)
updateMutationUseCounts ::
     [Int] -- ^ seznamů indexů použitých operací (indexy se mohou opakovat)
  -> Int -- ^ rozdíl skóre před provedením mutací a po nich
  -> MutationState -- ^ starý MutationState
  -> MutationState -- ^ nový MutationState
updateMutationUseCounts indices dif (p, S {sis, sic, sdc, sdl, shf}) =
  (p, listToS statsList) -- neměníme pravděpodobonsti, pouze měříme statistiky
  where
    go :: [Int] -> (Int, Int, Int, Int, Int) -> (Int, Int, Int, Int, Int)
    go [] tp = tp
    go (i:is) (a, b, c, d, e) -- spočítáme, kolikrát byla jaká operace v současném souboru použita
     =
      case i of
        0 -> go is (a + 1, b, c, d, e)
        1 -> go is (a, b + 1, c, d, e)
        2 -> go is (a, b, c + 1, d, e)
        3 -> go is (a, b, c, d + 1, e)
        4 -> go is (a, b, c, d, e + 1)
        _ -> error "Wrong number was picked"
    -- metrika udávající, jak pozitivní dopad operace má vzhledem k tomu, kolikrát byla použita
    dSO x =
      if usesTotal == 0
        then 0
        else fromIntegral dif * fromIntegral x / fromIntegral usesTotal
    usesList =
      let (a, b, c, d, e) = go indices (0, 0, 0, 0, 0)
       in [a, b, c, d, e]
    usesTotal = sum usesList
    -- k už nasbíraným dSO přidáme dSO se současného souboru mutací
    statsList =
      zipWith3
        (\uses dso dsos ->
           if uses == 0
             then dsos
             else dso : dsos)
        usesList
        (map dSO usesList)
        [sis, sic, sdc, sdl, shf]

-- | Mutuje alignment pouze jednou, mutační operaci volí na základě pravděpodobností v MutationState.
-- Vrací mimo jiné i Int, který udává, jaká oprace byla použita, a je později
-- využit k počítání pravděpodobnosti užití dané operace v další generaci.
mutateOnce :: StateT Alignment Mutator Int
mutateOnce = do
  Alignment {aProteins, aStartingGapSize} <- get
  operationIndex <- lift pickOperationIndex -- vybereme jednu z operací
  let operation = operations aStartingGapSize !! operationIndex
  (proteinIndex, protein) <- lift . lift . choose $ aProteins -- vybereme jeden z proteinů
  newProtein <- lift . lift . operation $ protein -- vybranou operaci aplikujeme na vybraný protein
  let newProteins = updateAt proteinIndex [newProtein] aProteins
  put $ Alignment newProteins (scoreProteins newProteins) aStartingGapSize -- updatujeme alignment
  return operationIndex
  where
    operations x = [insert x, increase, decrease, delete, shift]

-- | Vloží do proteinu na náhodné místo novou mezeru s náhodnou délkou.
-- Délka je určena exponenciální distribucí pravděpodobnosti, s průměrem rovným průměrné délce
-- mezer v počátečním alignmentu.
insert :: Double -> Protein -> RVar Protein
insert lambda Protein {pSeq, pGaps} = do
  gl <- exponential lambda
  gs <- chooseI pSeq
  -- vymažeme případnou mezeru, která už na gs začíná
  let newGaps = (gs, ceiling gl) : filter ((/= gs) . fst) pGaps
  return $ Protein pSeq newGaps

-- | Prodlouží v proteinu náhodnou mezeru o náhodnou délku.
-- Délka je určena exponenciální distribucí pravděpodobnosti, s průměrem 1.5.
increase :: Protein -> RVar Protein
increase prot@Protein {pGaps = []} = return prot
increase Protein {pSeq, pGaps} = do
  (index, (gs, gl)) <- choose pGaps
  lengthAdd <- round <$> (exponential 1.5 :: RVar Float)
  let newGaps = updateAt index [(gs, gl + lengthAdd)] pGaps -- přepíšeme stávající mezeru
  return $ Protein pSeq newGaps

-- | Zkrátí v proteinu náhodnou mezeru o 1. Pokud je délka mezery 1, vymaže ji.
decrease :: Protein -> RVar Protein
decrease prot@Protein {pGaps = []} = return prot
decrease Protein {pSeq, pGaps} = do
  (index, (gs, gl)) <- choose pGaps
  return $
    Protein pSeq $
    if gl - 1 == 0
      then updateAt index [] pGaps -- vymažeme stávající mezeru
      else updateAt index [(gs, gl - 1)] pGaps -- přepíšeme stávající mezeru

-- | Vybere v proteinu náhodnou pozici P a náhodnou mezeru M. Pokud na P již mezera je,
-- prohodí délka P a M. Pokud na P mezera není, přesune M na P.
shift :: Protein -> RVar Protein
shift prot@Protein {pGaps = []} = return prot
shift Protein {pSeq, pGaps} = do
  (_, (gs, gl)) <- choose pGaps
  targetIndex <- chooseI pSeq
  -- rozdělíme mezery na ty, které na gs začínají a na zbytek
  -- těch prvních by správně mělo být 0 nebo 1
  let (atTargetIndex, rest) =
        partition ((targetIndex ==) . fst) $ filter ((/=) gs . fst) pGaps
  return $
    Protein pSeq $
    case atTargetIndex of
      []         -> (targetIndex, gl) : filter ((/=) gs . fst) rest -- přesouváme mezeru
      [(hs, hl)] -> (hs, gl) : (gs, hl) : filter ((/=) gs . fst) rest -- prohazujeme délky
      _          -> error "Found gaps with duplicate starting points"

-- | Vymaže v proteinu náhodnou mezeru. Šance na to, že bude vybraná mezera opravdu smazána,
-- je rovna [2 * 1 / (délka mezery)], je tedy složitější vymazat dlouhé mezery.
delete :: Protein -> RVar Protein
delete s@Protein {pGaps = []} = return s
delete Protein {pSeq, pGaps} = do
  (index, (_, gl)) <- choose pGaps
  newGaps <-
    do coin <- stdUniform :: RVar Float
       return $
         if coin <= (2 / fromIntegral gl)
           then updateAt index [] pGaps
           else pGaps
  return $ Protein pSeq newGaps
