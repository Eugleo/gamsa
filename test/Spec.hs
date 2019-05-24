import Data.Random          (StdRandom (..), runRVar)
import Test.Hspec

import Model                (Alignment (..))
import MultipleSeqAlignment (defaultConfig, run)
import Utils                (mkAlignment)

testAl = mkAlignment ["KMMEEAABBGHGHI", "-MMABGHI", "EAABB-H-HI"]

main :: IO ()
main =
  hspec $
  describe "Hlavní část programu" $
  it "zvýší původní skóre" $ do
    let originalScore = aScore testAl
    result <- runRVar (run defaultConfig testAl) StdRandom
    let newScore = aScore result
    originalScore `shouldSatisfy` (>=) newScore
