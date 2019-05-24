import Data.Random          (StdRandom (..), runRVar)
import Test.Hspec

import Model                (Alignment (..))
import MultipleSeqAlignment (defaultConfig, run)
import Utils                (mkAlignment)

testAl = mkAlignment ["KMMEEAABBGHGHI", "-MMABGHI", "EAABB-H-HI"]

testAl2 = mkAlignment ["AA--", "--AA"]

main :: IO ()
main =
  hspec $
  describe "Hlavní část programu" $ do
    it "zvýší původní skóre" $ do
      let originalScore = aScore testAl
      result <- runRVar (run defaultConfig testAl) StdRandom
      let newScore = aScore result
      originalScore `shouldSatisfy` (>=) newScore
    it "počítá správně mezery i na konci" $
      aScore testAl2 `shouldBe` (-6 :: Int)
