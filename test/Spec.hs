import           Data.Random          (StdRandom (..), runRVar)
import           Model
import           MultipleSeqAlignment (run)
import           Test.Hspec
import           Utils                (mkAlignment)

testAl = mkAlignment ["KMMEEAABBGHGHI", "-MMABGHI", "EAABB-H-HI"]

main :: IO ()
main =
  hspec $
  describe "Hlavní část programu" $
  it "zvýší původní skóre" $ do
    let originalScore = aScore testAl
    result <- runRVar (run testAl) StdRandom
    let newScore = aScore result
    originalScore `shouldSatisfy` (>=) newScore
