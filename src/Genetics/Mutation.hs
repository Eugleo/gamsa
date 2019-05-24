{-# LANGUAGE NamedFieldPuns #-}

module Genetics.Mutation
  ( mutate
  ) where

import Control.Monad.Trans.Class            (lift)
import Control.Monad.Trans.State.Lazy       (StateT (..), modify)
import Data.List                            (group, partition, sort)
import Data.Random                          (RVar, uniform)
import Data.Random.Distribution.Exponential (exponential)

import Genetics.Scoring                     (scoreProteins)
import Model
import MutationProbabilities                (MutationState, Stats (..), listToS,
                                             pickOperationIndex)
import Utils                                (choose, updateAt)

type Mutator = StateT MutationState RVar Alignment

mutate :: Alignment -> Mutator
mutate al@Alignment {aScore} = do
  count <- lift (exponential (1 / 1.8) :: RVar Float)
  (opIndices, newAl) <- applyNTimes 1 mutate' al -- TODO: N times
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
    usesTotal =
      let (a, b, c, d, e) = usesTup
       in a + b + c + d + e
    dSO x =
      if usesTotal == 0
        then 0
        else fromIntegral dif * fromIntegral x / fromIntegral usesTotal
    usesTup = go indices (0, 0, 0, 0, 0)
    dSOList =
      let (a, b, c, d, e) = usesTup
       in [dSO a, dSO b, dSO c, dSO d, dSO e]
    statsList =
      let (a, b, c, d, e) = usesTup
       in zipWith3
            (\x xs uses ->
               if uses == 0
                 then xs
                 else x : xs)
            dSOList
            [sis, sic, sdc, sdl, shf]
            [a, b, c, d, e]

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
mutate' Alignment {aProteins} = do
  operationIndex <- pickOperationIndex
  let operation = operations !! operationIndex
  (proteinIndex, protein) <- lift $ choose aProteins
  newProtein <- lift $ operation protein
  let newProteins = updateAt proteinIndex [newProtein] aProteins
  let newAl = Alignment newProteins (scoreProteins newProteins)
  return (operationIndex, newAl)
  where
    operations = [insert, increase, decrease, delete, shift]

-- TODO: Add pattern matching and simpler updates
insert :: Protein -> RVar Protein
insert sq = do
  let s = pGaps sq
  l <- exponential (pMeanGapCount sq)
  i <- uniform 0 (length (pSeq sq) - 1)
  let newS = (i, ceiling l) : filter ((/= i) . fst) s
  return $ Protein (pSeq sq) newS (pMeanGapCount sq)

increase :: Protein -> RVar Protein
increase s@Protein {pGaps = []} = return s
increase Protein {pGaps = s, pSeq = aa, pMeanGapCount = seedGapMean} = do
  i <- uniform 0 (length s - 1)
  let (start, l) = s !! i
  let gs = updateAt i [(start, l + 1)] s
  return $ Protein aa gs seedGapMean

decrease :: Protein -> RVar Protein
decrease s@Protein {pGaps = []} = return s
decrease sq = do
  let s = pGaps sq
  i <- uniform 0 (length s - 1)
  let (start, l) = s !! i
  if l - 1 == 0
    then let newS = updateAt i [] s
          in return $ Protein (pSeq sq) newS (pMeanGapCount sq)
    else let newS = updateAt i [(start, l - 1)] s
          in return $ Protein (pSeq sq) newS (pMeanGapCount sq)

shift :: Protein -> RVar Protein
shift s@Protein {pGaps = []} = return s
shift sq = do
  let s = pGaps sq
  gapI <- uniform 0 (length s - 1)
  let (gs, gl) = s !! gapI
  i <- uniform 0 (length (pSeq sq) - 1)
  let (ms, nms) = partition ((i ==) . fst) $ filter ((/=) gs . fst) s
  if null ms
    then return $
         Protein
           (pSeq sq)
           ((i, gl) : filter ((/=) gs . fst) nms)
           (pMeanGapCount sq)
    else let [(hs, hl)] = ms
          in return $
             Protein
               (pSeq sq)
               ((hs, gl) : (gs, hl) : filter ((/=) gs . fst) nms)
               (pMeanGapCount sq)

delete :: Protein -> RVar Protein
delete s@Protein {pGaps = []} = return s
delete sq = do
  let s = pGaps sq
  i <- uniform 0 (length s - 1)
  let newS = updateAt i [] s
  return $ Protein (pSeq sq) newS (pMeanGapCount sq)
