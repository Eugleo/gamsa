module Lib where

import           Control.Monad                        (foldM, mapM_, replicateM)
import           Data.List                            (find, group, groupBy,
                                                       maximumBy, partition,
                                                       sort, sortOn)
import           Data.Ord                             (Down (..))
import           Data.Random
import           Data.Random.Distribution.Exponential (exponential)
import           Debug.Trace                          (traceShow)

type Generation = [Alignment]

type Alignment = [Seq]

data State =
  S (Double, Int, Double)
    (Double, Int, Double)
    (Double, Int, Double)
    (Double, Int, Double)
    (Double, Int, Double)
  deriving (Show, Eq, Ord)

data Seq = Seq
  { len         :: Int
  , aa          :: String
  , gaps        :: [Gap]
  , seedGapMean :: Double -- TODO: Move this to alignment or someplace else
  } deriving (Show, Eq)

type Gap = (Int, Int)

getRandomR :: (Int, Int) -> RVar Int
getRandomR (a, b) = uniform a b

runIO :: Alignment -> IO Alignment
runIO a = do
  result <- runRVar (run a) StdRandom
  mapM_ putStrLn . fill $ result
  return result

run :: Alignment -> RVar Alignment
run a = do
  let state = S (0.2, 0, 0) (0.2, 0, 0) (0.2, 0, 0) (0.2, 0, 0) (0.2, 0, 0)
  startingG <- populate state a
  fin <- doRuns state startingG 2000
  return $ maximumOn score fin

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
  traceShow
    ("Remaining runs: " ++ show n)
    (nextGen st g >>= \(ng, ns) -> doRuns ns ng (n - 1))

-- TODO: Find out if we should mutate the top5 as well
nextGen :: State -> Generation -> RVar (Generation, State)
nextGen st g = do
  let top = top5 g
  (als, ns) <- go st 95 []
  return (als ++ top, newState ns)
  where
    go :: State -> Int -> [Alignment] -> RVar ([Alignment], State)
    go s 0 acc = return (acc, s)
    go s n acc = do
      best <- tournament g
      (ng, ns) <- mutate s best
      go ns (n - 1) (ng : acc)

top5 :: Generation -> [Alignment]
top5 = take 5 . sortOn (Down . score)

tournament :: Generation -> RVar Alignment
tournament g = do
  indices <- replicateM 5 $ getRandomR (0, length g - 1)
  return $ head $ sortOn (Down . score) $ map (g !!) indices

score :: Alignment -> Double
score al =
  1000 * meanColumnHomogenicity al + 2000 * gapBlocks al -
  2 * columnsIncrement al

meanColumnHomogenicity :: Alignment -> Double
meanColumnHomogenicity al =
  mean $ map (columnHomogenicity (fill al)) [0 .. len (head al) - 1]
  where
    mean :: [Double] -> Double
    mean ns = sum ns / fromIntegral (length ns)

fill :: Alignment -> [String]
fill al = fill' (maximum $ map (\s -> len s + sum (map snd (gaps s))) al) al
  where
    fill' _ [] = []
    fill' n (x:xs) =
      let res = go x
          l = length res
       in (if l < n
             then res ++ replicate (n - l) '-'
             else res) :
          fill' n xs
    go Seq {aa = aa, len = _, gaps = gs} = foldr (helper gs) "" (zip aa [1 ..])
    helper :: [Gap] -> (Char, Int) -> String -> String
    helper gs (c, i) acc =
      case find ((== i) . fst) gs of
        Nothing     -> c : acc
        Just (_, l) -> c : replicate l '-' ++ acc

columnHomogenicity :: [String] -> Int -> Double
columnHomogenicity al n = fromIntegral numerator / fromIntegral denominator
  where
    numerator :: Int
    numerator =
      sum . map (sqr . countWithoutGaps) . group . sort . map (!! n) $ al
    denominator :: Int
    denominator =
      sqr . sum . map countWithoutGaps . group . sort . map (!! n) $ al
    countWithGaps :: [a] -> Int
    countWithGaps = length
    countWithoutGaps ('-':_) = 0
    countWithoutGaps xs      = length xs
    sqr x = x * x

gapBlocks :: Alignment -> Double
gapBlocks seqs =
  if gb == 0
    then 0
    else 1 / gb
  where
    gb = sum (map (fromIntegral . length . gaps) seqs)

columnsIncrement :: Alignment -> Double
columnsIncrement seqs =
  if mSeq <= 1
    then 0
    else fromIntegral newMax / fromIntegral mSeq - 1
  where
    mSeq = maximum $ map len seqs
    newMax = maximum $ map (\s -> len s + sum (map snd (gaps s))) seqs

mutate :: State -> Alignment -> RVar (Alignment, State)
mutate st@(S a b c d e) al = do
  check <- getRandomR (1, 10)
  if check >= 2
    then do
      p <- replicateM 1 stdUniform
      let oldScore = score al
      let opIndices = map (pick st) p
      let op = (ops !!) <$> opIndices
      let sts = (ps !!) <$> opIndices
      seqIndex <- getRandomR (0, length al - 1)
      let s = al !! seqIndex
      newS <- foldM (\acc o -> o acc) s op
      let newA = take seqIndex al ++ [newS] ++ drop (seqIndex + 1) al
      let newScore = score newA
      let newSt = go st sts opIndices (newScore - oldScore)
      return (newA, newSt)
    else return (al, st)
  where
    ops = [insert, increase, decrease, delete, shift]
    ps = [a, b, c, d, e]
    go stt _ [] _ = stt
    go stt [] _ _ = stt
    go stt ((prob, tot, d1):muts) (i:is) diff =
      let s = updateMutation stt (prob, tot, d1 + diff) i
       in go s muts is diff

divideState :: State -> [(Double, Int, Double)]
divideState (S a b c d e) = [a, b, c, d, e]

updateMutation :: State -> (Double, Int, Double) -> Int -> State
updateMutation st (prob, tot, diff) i =
  let mut = divideState st
   in createState $ take i mut ++ [(prob, tot + 1, diff)] ++ drop (i + 1) mut

createState :: [(Double, Int, Double)] -> State
createState [a, b, c, d, e] = S a b c d e
createState xs              = undefined

newState :: State -> State
newState (S a b c d e) =
  S (suPr a, 0, 0) (suPr b, 0, 0) (suPr c, 0, 0) (suPr d, 0, 0) (suPr e, 0, 0)
  where
    every = [a, b, c, d, e]
    total = sum . map (\(_, x, _) -> x) $ every
    maxTso = maximum . map (abs . tso) $ every
    tso (_, tot, diff) = diff * fromIntegral tot / fromIntegral total
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
insert :: Seq -> RVar Seq
insert sq = do
  let s = gaps sq
  l <- round <$> exponential (seedGapMean sq)
  i <- getRandomR (1, len sq - 1)
  let newS = (i, l) : filter ((/= i) . fst) s
  return $ Seq (len sq) (aa sq) newS (seedGapMean sq)

increase :: Seq -> RVar Seq
increase s@Seq {gaps = []} = return s
increase sq = do
  let s = gaps sq
  i <- getRandomR (0, length s - 1)
  let (start, l) = s !! i
  let gs = take i s ++ [(start, l + 1)] ++ drop (i + 1) s
  return $ Seq (len sq) (aa sq) gs (seedGapMean sq)

decrease :: Seq -> RVar Seq
decrease s@Seq {gaps = []} = return s
decrease sq = do
  let s = gaps sq
  i <- getRandomR (0, length s - 1)
  let (start, l) = s !! i
  if l - 1 == 0
    then let newS = take i s ++ drop (i + 1) s
          in return $ Seq (len sq) (aa sq) newS (seedGapMean sq)
    else let newS = take i s ++ [(start, l - 1)] ++ drop (i + 1) s
          in return $ Seq (len sq) (aa sq) newS (seedGapMean sq)

shift :: Seq -> RVar Seq
shift s@Seq {gaps = []} = return s
shift sq = do
  let s = gaps sq
  gapI <- getRandomR (0, length s - 1)
  let (gs, gl) = s !! gapI
  i <- getRandomR (0, len sq - 1)
  let (ms, nms) = partition ((i ==) . fst) $ filter ((/=) gs . fst) s
  if null ms
    then return $
         Seq
           (len sq)
           (aa sq)
           ((i, gl) : filter ((/=) gs . fst) nms)
           (seedGapMean sq)
    else let [(hs, hl)] = ms
          in return $
             Seq
               (len sq)
               (aa sq)
               ((hs, gl) : (gs, hl) : filter ((/=) gs . fst) nms)
               (seedGapMean sq)

delete :: Seq -> RVar Seq
delete s@Seq {gaps = []} = return s
delete sq = do
  let s = gaps sq
  i <- getRandomR (0, length s - 1)
  let newS = take i s ++ drop (i + 1) s
  return $ Seq (len sq) (aa sq) newS (seedGapMean sq)
