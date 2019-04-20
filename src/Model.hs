module Model
  ( Generation
  , Alignment(..)
  , State(..)
  , Protein(..)
  , Gap
  ) where

import           Data.Vector (Vector)

type Generation = [Alignment]

data Alignment = Alignment
  { aProteins :: [Protein]
  , aScore    :: Int
  } deriving (Show)

instance Eq Alignment where
  Alignment {aScore = s1} == Alignment {aScore = s2} = s1 == s2

instance Ord Alignment where
  Alignment {aScore = s1} `compare` Alignment {aScore = s2} = s1 `compare` s2

data Protein = Protein
  { pSeq          :: Vector Char
  , pGaps         :: [Gap]
  , pMeanGapCount :: Double -- TODO: Move this to alignment or someplace else
  } deriving (Show, Eq, Ord)

-- (prob, tot, diff)
data State =
  S (Double, Int, Int)
    (Double, Int, Int)
    (Double, Int, Int)
    (Double, Int, Int)
    (Double, Int, Int)
  deriving (Show, Eq)

type Gap = (Int, Int)
