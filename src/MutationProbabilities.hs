module MutationProbabilities where

-- (prob, tot, diff)
data State =
  S (Double, Int, Int)
    (Double, Int, Int)
    (Double, Int, Int)
    (Double, Int, Int)
    (Double, Int, Int)
  deriving (Show, Eq)

divideState :: State -> [(Double, Int, Int)]
divideState (S a b c d e) = [a, b, c, d, e]

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
