{-# LANGUAGE MultiWayIf #-}

module MultipleSeqAlignment
  ( run
  ) where

import           Control.Monad                        (replicateM)
import           Crossover                            (recombineH, recombineV)
import           Data.Function                        (on)
import           Data.List                            (sortBy)
import           Data.Ord                             (Down (..))
import           Data.Random                          (RVar, stdUniform,
                                                       uniform)
import           Data.Random.Distribution.Exponential (exponential)
import           Debug.Trace                          (traceShow)
import           Model                                (Alignment (..),
                                                       Generation)
import           Mutation
import           MutationProbabilities
import           Utils                                (between)

run :: Alignment -> RVar Alignment
run a = do
  let state = S (0.3, 0, 0) (0.2, 0, 0) (0.2, 0, 0) (0.1, 0, 0) (0.2, 0, 0)
  startingG <- populate state a
  fin <- doRuns state startingG 20000
  return $ maximum fin

populate :: State -> Alignment -> RVar Generation
populate st a = go 100 []
  where
    go 0 acc = return acc
    go n acc = mutate st a >>= \(na, _) -> go (n - 1) (na : acc)

doRuns :: State -> Generation -> Int -> RVar Generation
doRuns _ g 0 = return g
doRuns st g n =
  if n `mod` 100 == 0
    then traceShow
           ( "Remaining runs: " ++ show n
           , (map aScore (take 3 tops), map aScore (drop 97 tops)))
           result
    else result
  where
    result = nextGen st g >>= \(ng, ns) -> doRuns st ng (n - 1)
    tops = sortBy (compare `on` Down) g

-- TODO: Find out if we should mutate the top5 as well
nextGen :: State -> Generation -> RVar (Generation, State)
nextGen st g = do
  let top = top5 g
  (als, ns) <- go st 19 []
  return (top ++ als, newState ns) -- newState ns
  where
    go :: State -> Int -> [Alignment] -> RVar (Generation, State)
    go s 0 acc = return (acc, s)
    go s n acc = do
      best <- tournament' g
      k <- exponential (1 / 1.8) :: RVar Float
      more5 <- mapM (applyNTimes (round k) s mutate) best
      go s (n - 1) (map fst more5 ++ acc) -- get new state

applyNTimes ::
     Monad m
  => Int
  -> State
  -> (State -> a -> m (a, State))
  -> a
  -> m (a, State)
applyNTimes 0 st _ a = return (a, st)
applyNTimes k st f a = do
  (result, newSt) <- f st a
  applyNTimes (k - 1) st f result -- zde použít newst

top5 :: Generation -> [Alignment]
top5 gen = take 5 tops
  where
    tops = sortBy (compare `on` Down) gen

tournament' :: Generation -> RVar [Alignment]
tournament' g = do
  indices <- replicateM 5 $ uniform 0 99
  let contenders = (g !!) <$> indices
  let first = maximum contenders
  coin <- stdUniform :: RVar Float
  if | coin `between` (0, 0.3) -> mapM (recombineH first) contenders
     | coin `between` (0.3, 0.8) -> mapM (recombineV first) contenders
     | otherwise -> return contenders
