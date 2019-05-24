{-# LANGUAGE MultiWayIf     #-}
{-# LANGUAGE NamedFieldPuns #-}

module MultipleSeqAlignment
  ( run
  , defaultConfig
  ) where

import Control.Monad                  (replicateM)
import Control.Monad.Trans.Class      (lift)
import Control.Monad.Trans.State.Lazy (StateT (..), evalStateT, get, modify,
                                       put)
import Data.Function                  (on)
import Data.List                      (sortBy)
import Data.Ord                       (Down (..))
import Data.Random                    (RVar, stdUniform)
import Debug.Trace                    (traceShow)

import Genetics.Crossover             (recombineH, recombineV)
import Genetics.Mutation              (mutate)
import Model                          (Alignment (..), Generation)
import MutationProbabilities          (MutationState, Probabilities (..),
                                       blankStats, nextGenProbabilities,
                                       sToList)
import Utils                          (between, choose)

data Config = Config
  { generationCount       :: Int
  , startingProbabilities :: Probabilities
  , tournamentSize        :: Int
  , tournamentCount       :: Int
  , eliteCount            :: Int
  }

defaultConfig :: Config
defaultConfig = Config 5000 (P 0.2 0.2 0.2 0.2 0.2) 3 33 1

run :: Config -> Alignment -> RVar Alignment
run config@Config { generationCount
                  , startingProbabilities
                  , tournamentSize
                  , tournamentCount
                  , eliteCount
                  } al = do
  let state = (startingProbabilities, blankStats)
  let populationSize = tournamentSize * tournamentCount + eliteCount
  generation <- mkPopulation populationSize al state
  fin <- evalStateT (repeatNGenerations config generation generationCount) state
  return $ maximum fin

mkPopulation :: Int -> Alignment -> MutationState -> RVar Generation
mkPopulation size al = evalStateT (go size [])
  where
    go 0 acc = return acc
    go n acc = do
      oldState <- get
      newAl <- mutate al
      put oldState
      go (n - 1) (newAl : acc)

repeatNGenerations ::
     Config -> Generation -> Int -> StateT MutationState RVar Generation
repeatNGenerations _ g 0 = return g
repeatNGenerations config g n = do
  (prob, _) <- get
  let message =
        "Run " ++
        show n ++
        ", top 3: " ++
        show (map aScore $ take 3 tops) ++ ", p: " ++ show (stateMsg prob)
  if n `mod` 100 == 0
    then traceShow message result
    else result
  where
    stateMsg P {pis, pic, pdc, pdl, phf} =
      map (round . (* 100)) [pis, pic, pdc, pdl, phf]
    result =
      nextGeneration config g >>= \newGeneration ->
        repeatNGenerations config newGeneration (n - 1)
    tops = sortBy (compare `on` Down) g

nextGeneration :: Config -> Generation -> StateT MutationState RVar Generation
nextGeneration Config {tournamentSize, tournamentCount, eliteCount} g = do
  let topN = top eliteCount g
  alignments <- concat <$> replicateM tournamentCount go
  modify (\st@(_, s) -> (nextGenProbabilities st, blankStats))
  return $ topN ++ alignments
  where
    go = do
      elite <- lift $ tournament tournamentSize g
      mapM mutate elite

top :: Int -> Generation -> [Alignment]
top n = take n . sortBy (compare `on` Down)

tournament :: Int -> Generation -> RVar [Alignment]
tournament size g = do
  contenders <- fmap (map snd) . replicateM size $ choose g
  let first = maximum contenders
  coin <- stdUniform :: RVar Float
  if | coin `between` (0, 0.3) -> mapM (recombineH first) contenders
     | coin `between` (0.3, 0.8) -> mapM (recombineV first) contenders
     | otherwise -> return contenders
