{-# LANGUAGE NamedFieldPuns #-}

module Genetics.Mutation
  ( mutate
  ) where

import Control.Monad.Trans.Class            (lift)
import Control.Monad.Trans.State.Lazy       (StateT (..), modify)
import Data.List                            (partition)
import Data.Random                          (RVar, stdUniform)
import Data.Random.Distribution.Exponential (exponential)

import Genetics.Scoring                     (scoreProteins)
import Model
import MutationProbabilities                (MutationState, Stats (..), listToS,
                                             pickOperationIndex)
import Utils                                (choose, chooseI, updateAt)

type Mutator = StateT MutationState RVar Alignment

mutate :: Alignment -> Mutator
mutate al@Alignment {aScore} = do
  count <- lift (exponential (1 / 1.8) :: RVar Float)
  (opIndices, newAl) <- applyNTimes (round count) mutate' al
  let dif = scoreProteins (aProteins newAl) - aScore
  modify $ updateMutationUseCounts opIndices dif
  return newAl

updateMutationUseCounts :: [Int] -> Int -> MutationState -> MutationState
updateMutationUseCounts indices dif (p, S {sis, sic, sdc, sdl, shf}) =
  (p, listToS statsList)
  where
    go :: [Int] -> (Int, Int, Int, Int, Int) -> (Int, Int, Int, Int, Int)
    go [] tp = tp
    go (i:is) (a, b, c, d, e) =
      case i of
        0 -> go is (a + 1, b, c, d, e)
        1 -> go is (a, b + 1, c, d, e)
        2 -> go is (a, b, c + 1, d, e)
        3 -> go is (a, b, c, d + 1, e)
        4 -> go is (a, b, c, d, e + 1)
        _ -> error "Wrong number was picked"
    dSO x =
      if usesTotal == 0
        then 0
        else fromIntegral dif * fromIntegral x / fromIntegral usesTotal
    usesList =
      let (a, b, c, d, e) = go indices (0, 0, 0, 0, 0)
       in [a, b, c, d, e]
    usesTotal = sum usesList
    statsList =
      zipWith3
        (\uses dso dsos ->
           if uses == 0
             then dsos
             else dso : dsos)
        usesList
        (map dSO usesList)
        [sis, sic, sdc, sdl, shf]

applyNTimes ::
     Int
  -> (Alignment -> StateT MutationState RVar (Int, Alignment))
  -> Alignment
  -> StateT MutationState RVar ([Int], Alignment)
applyNTimes 0 _ a = return ([], a)
applyNTimes k f a = do
  (index, al) <- f a
  (indices, finalAl) <- applyNTimes (k - 1) f al
  return (index : indices, finalAl)

mutate' :: Alignment -> StateT MutationState RVar (Int, Alignment)
mutate' Alignment {aProteins, aStartingGapSize} = do
  operationIndex <- pickOperationIndex
  let operation = operations !! operationIndex
  (proteinIndex, protein) <- lift $ choose aProteins
  newProtein <- lift $ operation protein
  let newProteins = updateAt proteinIndex [newProtein] aProteins
  let newAl = Alignment newProteins (scoreProteins newProteins) aStartingGapSize
  return (operationIndex, newAl)
  where
    operations = [insert aStartingGapSize, increase, decrease, delete, shift]

insert :: Double -> Protein -> RVar Protein
insert lambda Protein {pSeq, pGaps} = do
  l <- exponential lambda
  i <- chooseI pSeq
  let newGaps = (i, ceiling l) : filter ((/= i) . fst) pGaps
  return $ Protein pSeq newGaps

increase :: Protein -> RVar Protein
increase prot@Protein {pGaps = []} = return prot
increase Protein {pSeq, pGaps} = do
  (index, (gs, gl)) <- choose pGaps
  lengthAdd <- round <$> (exponential (1 / 1.5) :: RVar Float)
  let newGaps = updateAt index [(gs, gl + lengthAdd)] pGaps
  return $ Protein pSeq newGaps

decrease :: Protein -> RVar Protein
decrease prot@Protein {pGaps = []} = return prot
decrease Protein {pSeq, pGaps} = do
  (index, (gs, gl)) <- choose pGaps
  return $
    Protein pSeq $
    if gl - 1 == 0
      then updateAt index [] pGaps
      else updateAt index [(gs, gl - 1)] pGaps

shift :: Protein -> RVar Protein
shift prot@Protein {pGaps = []} = return prot
shift Protein {pSeq, pGaps} = do
  (_, (gs, gl)) <- choose pGaps
  targetIndex <- chooseI pSeq
  let (atTargetIndex, rest) =
        partition ((targetIndex ==) . fst) $ filter ((/=) gs . fst) pGaps
  return $
    Protein pSeq $
    case atTargetIndex of
      []         -> (targetIndex, gl) : filter ((/=) gs . fst) rest
      [(hs, hl)] -> (hs, gl) : (gs, hl) : filter ((/=) gs . fst) rest
      _          -> error "Found gaps with duplicate starting points"

delete :: Protein -> RVar Protein
delete s@Protein {pGaps = []} = return s
delete Protein {pSeq, pGaps} = do
  (index, (_, gl)) <- choose pGaps
  newGaps <-
    do coin <- stdUniform :: RVar Float
       return $
         if coin <= (1 / fromIntegral gl)
           then updateAt index [] pGaps
           else pGaps
  return $ Protein pSeq newGaps
