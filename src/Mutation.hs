module Mutation
  ( mutate
  ) where

import           Control.Monad                        (foldM, replicateM)
import           Data.List                            (partition)
import           Data.Random                          (RVar, stdUniform,
                                                       uniform)
import           Data.Random.Distribution.Exponential (exponential)
import           Model
import           MutationProbabilities
import           Scoring                              (scoreProteins)

mutate :: State -> Alignment -> RVar (Alignment, State)
mutate st@(S a b c d e) Alignment {aProteins = seqs, aScore = oldScore} = do
  p <- replicateM 1 stdUniform
  let opIndices = map (pick st) p
  let op = (ops !!) <$> opIndices
  let sts = (ps !!) <$> opIndices
  seqIndex <- uniform 0 (length seqs - 1)
  let s = seqs !! seqIndex
  newS <- foldM (\acc o -> o acc) s op
  let newSequences = updateAt seqIndex [newS] seqs
  let newScore = scoreProteins newSequences
  let newSt = go st sts opIndices (newScore - oldScore)
  return (Alignment newSequences newScore, newSt)
  where
    ops = [insert, increase, decrease, delete, shift]
    ps = [a, b, c, d, e]
    go stt _ [] _ = stt
    go stt [] _ _ = stt
    go stt ((prob, tot, d1):muts) (i:is) diff =
      let s = updateMutation stt (prob, tot, d1 + diff) i
       in go s muts is diff

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

updateAt :: Int -> [a] -> [a] -> [a]
updateAt i new xs = take i xs ++ new ++ drop (i + 1) xs

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

updateMutation :: State -> (Double, Int, Int) -> Int -> State
updateMutation st (prob, tot, diff) i =
  let mut = divideState st
   in createState $ updateAt i [(prob, tot + 1, diff)] mut

delete :: Protein -> RVar Protein
delete s@Protein {pGaps = []} = return s
delete sq = do
  let s = pGaps sq
  i <- uniform 0 (length s - 1)
  let newS = updateAt i [] s
  return $ Protein (pSeq sq) newS (pMeanGapCount sq)
