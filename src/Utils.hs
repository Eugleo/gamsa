module Utils
  ( maximumOn
  , between
  , fill
  ) where

import           Data.List   (maximumBy, sortOn)
import           Data.Vector (toList)
import           Model       (Protein (..))

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
