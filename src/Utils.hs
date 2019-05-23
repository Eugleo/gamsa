module Utils
  ( maximumOn
  , between
  , fill
  , mkAlignment
  ) where

import           Data.List        (groupBy, maximumBy, sortOn)
import           Data.Vector      (fromList, toList)
import           Genetics.Scoring (scoreProteins)
import           Model            (Alignment (..), Protein (..))

maximumOn :: Ord b => (a -> b) -> [a] -> a
maximumOn f = maximumBy (\a b -> compare (f a) (f b))

between :: Ord a => a -> (a, a) -> Bool
between x (l, r) = (l <= x) && (x <= r)

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

mkAlignment :: [String] -> Alignment
mkAlignment input = Alignment seqs (scoreProteins seqs)
  where
    seqs = map mkSeq input
    mkSeq str = Protein (aa str) (gps str) (1 / seed)
    seed =
      if denominator == 0
        then 1
        else numerator / denominator
    numerator = fromIntegral (sum (map snd (concatMap gps input)))
    denominator = fromIntegral (length (concatMap gps input))
    aa = fromList . filter (not . isGap)
    gps = go 0 [] . groupBy (\a b -> isGap a && isGap b)
    isGap = (==) '-'
    go i acc (x:xs) =
      go
        (if isGap (head x)
           then i
           else i + 1)
        (if isGap (head x) && not (null xs)
           then (i, length x) : acc
           else acc)
        xs
    go _ acc [] = acc
