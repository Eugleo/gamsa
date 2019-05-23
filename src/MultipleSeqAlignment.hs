module MultipleSeqAlignment where

import           Control.Monad (foldM, mapM_, replicateM)
import           Data.Function (on)
import           Data.List     (maximumBy, partition, sortBy, sortOn)
import           Data.Ord      (Down (..))
import           Data.Random   (RVar, StdRandom (..), runRVar, stdUniform,
                                uniform)
import           Data.Vector   (toList)
import           Debug.Trace   (traceShow)
import           Model         (Alignment (..), Gap, Generation, Protein (..),
                                State (..))
import           Scoring       (scorePair)

exponential :: Double -> RVar Double
exponential lambda = do
  x <- stdUniform
  if x == 0
    then exponential lambda
    else return $ negate lambda * log x

getRandomR :: (Int, Int) -> RVar Int
getRandomR (a, b) = uniform a b

runIO :: Alignment -> IO Alignment
runIO a = do
  result <- runRVar (run a) StdRandom
  mapM_ putStrLn . fill . aProteins $ result
  return result

testrecomb :: Alignment -> Alignment -> IO Alignment
testrecomb a b = do
  result <- runRVar (recombineV a b) StdRandom
  mapM_ putStrLn . fill . aProteins $ result
  return result

fill :: [Protein] -> [String]
fill al = map (fill' maxSequenceLength) al
  where
    maxSequenceLength = seqLength . maximumOn seqLength $ al
    seqLength s = length (pSeq s) + foldr (\g acc -> snd g + acc) 0 (pGaps s)

fill' :: Int -> Protein -> String
fill' n Protein {pSeq = aa, pGaps = gaps} =
  let result = go (sortOn fst gaps) (zip [0 ..] $ toList aa) ""
   in result ++ replicate (n - length result) '-'
  where
    go :: [(Int, Int)] -> [(Int, Char)] -> String -> String
    go _ [] acc = reverse acc
    go [] str acc = reverse acc ++ map snd str
    go gps@((start, leng):gs) ((i, c):xs) acc
      | start <= i = go gs xs (c : replicate leng '-' ++ acc)
      | otherwise = go gps xs (c : acc)

run :: Alignment -> RVar Alignment
run a = do
  let state = S (0.2, 0, 0) (0.2, 0, 0) (0.2, 0, 0) (0.2, 0, 0) (0.2, 0, 0)
  startingG <- populate state a
  fin <- doRuns state startingG 100000
  return $ maximum fin

maximumOn :: Ord b => (a -> b) -> [a] -> a
maximumOn f = maximumBy (\a b -> compare (f a) (f b))

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
           (nextGen st g >>= \(ng, ns) -> doRuns st ng (n - 1))
    else nextGen st g >>= \(ng, ns) -> doRuns st ng (n - 1)
  where
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
      more5 <- mapM (mutate s) best
      go s (n - 1) (map fst more5 ++ acc) -- get new state

top5 :: Generation -> [Alignment]
top5 gen = take 5 tops
  where
    tops = sortBy (compare `on` Down) gen

tournament :: Generation -> RVar Alignment
tournament g = do
  indices <- replicateM 5 $ getRandomR (0, 99)
  return $ maximum $ (g !!) <$> indices

tournament' :: Generation -> RVar [Alignment]
tournament' g = do
  indices <- replicateM 5 $ getRandomR (0, 99)
  let contenders = (g !!) <$> indices
  let first = maximum contenders
  mapM (recombineH first) contenders

recombineH :: Alignment -> Alignment -> RVar Alignment
recombineH a@Alignment {aProteins = protA} Alignment {aProteins = protB}
  | protA == protB = return a
  | otherwise = do
    newProt <- mapM recombineProt (zip protA protB)
    return $ Alignment newProt (aScore a)
  where
    recombineProt (a, Protein {pGaps = gapsB}) = do
      coin <- stdUniform
      if coin
        then return a
        else return $ Protein (pSeq a) gapsB (pMeanGapCount a)

-- TODO: Implement better safety measures than simple <
scoreProteins :: [Protein] -> Int
scoreProteins al = sum [scorePair a b | a <- al, b <- al, a > b]

-- TODO: Implement 1.2 mutation probability
mutate :: State -> Alignment -> RVar (Alignment, State)
mutate st@(S a b c d e) al@Alignment {aProteins = seqs, aScore = oldScore} = do
  check <- getRandomR (1, 100)
  if check >= 2
    then do
      p <- replicateM 1 stdUniform
      let opIndices = map (pick st) p
      let op = (ops !!) <$> opIndices
      let sts = (ps !!) <$> opIndices
      seqIndex <- getRandomR (0, length seqs - 1)
      let s = seqs !! seqIndex
      newS <- foldM (\acc o -> o acc) s op
      let newSequences = updateAt seqIndex [newS] seqs
      let newScore = scoreProteins newSequences
      let newSt = go st sts opIndices (newScore - oldScore)
      return (Alignment newSequences newScore, newSt)
    else return (al, st)
  where
    ops = [insert, increase, decrease, delete, shift]
    ps = [a, b, c, d, e]
    go stt _ [] _ = stt
    go stt [] _ _ = stt
    go stt ((prob, tot, d1):muts) (i:is) diff =
      let s = updateMutation stt (prob, tot, d1 + diff) i
       in go s muts is diff

divideState :: State -> [(Double, Int, Int)]
divideState (S a b c d e) = [a, b, c, d, e]

updateMutation :: State -> (Double, Int, Int) -> Int -> State
updateMutation st (prob, tot, diff) i =
  let mut = divideState st
   in createState $ updateAt i [(prob, tot + 1, diff)] mut

createState :: [(Double, Int, Int)] -> State
createState [a, b, c, d, e] = S a b c d e
createState _               = undefined

newState :: State -> State
newState (S a b c d e) =
  S (suPr a, 0, 0) (suPr b, 0, 0) (suPr c, 0, 0) (suPr d, 0, 0) (suPr e, 0, 0)
  where
    every = [a, b, c, d, e]
    total = sum . map (\(_, x, _) -> x) $ every
    maxTso = maximum . map (abs . tso) $ every
    tso (_, tot, diff) =
      fromIntegral diff * fromIntegral tot / fromIntegral total
    pr x@(p, _, _) = p + p * 0.1 * (tso x / maxTso)
    totalPr = sum . map pr $ every
    suPr x = pr x / totalPr

pick :: State -> Double -> Int
pick (S (ins, _, _) (inc, _, _) (dec, _, _) (del, _, _) (shf, _, _)) num =
  go 0 0 ins [ins, inc, dec, del, shf]
  where
    go n _ _ [] = n - 1
    go n l r (p:ps)
      | (l <= num) && (num <= r) = n
      | otherwise = go (n + 1) r (r + p) ps

-- TODO: Add pattern matching and simpler updates
insert :: Protein -> RVar Protein
insert sq = do
  let s = pGaps sq
  l <- exponential (pMeanGapCount sq)
  i <- getRandomR (0, length (pSeq sq) - 1)
  let newS = (i, ceiling l) : filter ((/= i) . fst) s
  return $ Protein (pSeq sq) newS (pMeanGapCount sq)

increase :: Protein -> RVar Protein
increase s@Protein {pGaps = []} = return s
increase Protein {pGaps = s, pSeq = aa, pMeanGapCount = seedGapMean} = do
  i <- getRandomR (0, length s - 1)
  let (start, l) = s !! i
  let gs = updateAt i [(start, l + 1)] s
  return $ Protein aa gs seedGapMean

decrease :: Protein -> RVar Protein
decrease s@Protein {pGaps = []} = return s
decrease sq = do
  let s = pGaps sq
  i <- getRandomR (0, length s - 1)
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
  gapI <- getRandomR (0, length s - 1)
  let (gs, gl) = s !! gapI
  i <- getRandomR (0, length (pSeq sq) - 1)
  let (ms, nms) = partition ((i ==) . fst) $ filter ((/=) gs . fst) s
  if null ms
    then return $
         Protein
           (pSeq sq)
           ((i, gl) : filter ((/=) gs . fst) nms)
           (pMeanGapCount sq)
    else let [(hs, hl)] = traceShow ms ms
          in return $
             Protein
               (pSeq sq)
               ((hs, gl) : (gs, hl) : filter ((/=) gs . fst) nms)
               (pMeanGapCount sq)

delete :: Protein -> RVar Protein
delete s@Protein {pGaps = []} = return s
delete sq = do
  let s = pGaps sq
  i <- getRandomR (0, length s - 1)
  let newS = updateAt i [] s
  return $ Protein (pSeq sq) newS (pMeanGapCount sq)
