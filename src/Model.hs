module Model
  ( Generation
  , Alignment(..)
  , Protein(..)
  , Gap
  ) where

import Data.Vector (Vector)

type Generation = [Alignment]

-- | Alignment je nějaké seřazení sekvencí pod sebe, tj. vlastně nějaké seskupení mezer.
data Alignment = Alignment
  { aProteins        :: [Protein] -- ^ seznam proteinů
  , aScore           :: Int -- ^ skóre celého alignmentu
  , aStartingGapSize :: Double -- počáteční průměrná délka mezer (je konstantní napříč všemi alignmenty)
  } deriving (Show)

instance Eq Alignment where
  Alignment {aScore = s1} == Alignment {aScore = s2} = s1 == s2

instance Ord Alignment where
  Alignment {aScore = s1} `compare` Alignment {aScore = s2} = s1 `compare` s2

-- | Protein je kombinace sekvence a mezer; sekvence samotná je původní textová sekvence aminokyselin (AK).
-- Sekvence zlstávají během alignmentu konstantí, zatímco mezery se mění.
data Protein = Protein
  { pSeq  :: Vector Char
  , pGaps :: [Gap]
  } deriving (Show, Eq, Ord)

-- | Gap popisuje mezeru v proteinu pomocí dvojice souřadnic: (pozice v sekvenci, délka mezery).
-- | To znamená, že souřadnice mezery se nezmění, když před ní v rámci proteinu vložíme další mezeru,
-- protože *pozice v sekvenci* je brána vůči konstantímu seznamu AK, a ne vzhledem k proteinu i s mezerami.
type Gap = (Int, Int)
