module Align where

import qualified Data.Map.Strict as Map

dna :: Map.Map (Char, Char) Int
dna =
  Map.fromList $
  [(('A', 'A'), 5), (('C', 'C'), 5), (('G', 'G'), 5), (('T', 'T'), 5)] ++
  [((a, b), -4) | a <- "ACGT", b <- "ACGT", a /= b]

nw :: String -> String -> Int
nw = undefined
