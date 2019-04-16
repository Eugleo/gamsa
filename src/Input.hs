module Input where

import           Data.List   (groupBy)
import           Data.Vector (fromList)
import           Lib

mkAl :: [String] -> Alignment
mkAl input = Alignment (fromList $ map mkSeq input) 0
  where
    mkSeq str = Seq (length (aa str)) (fromList $ aa str) (gps str) seed
    seed :: Double
    seed =
      fromIntegral (sum (map snd (concatMap gps input))) /
      fromIntegral (length (concatMap gps input))
    aa = filter (not . isGap)
    gps = go 0 [] . groupBy (\a b -> isGap a && isGap b)
    isGap = (==) '-'
    go i acc (x:xs) =
      go
        (i + 1)
        (if isGap (head x)
           then (i, length x) : acc
           else acc)
        xs
    go _ acc [] = acc

inp = mkAl [aa1, aa2, aa3, aa4, aa5, aa6, aa7]

aa1 =
  "---------------------------------------------------------------------------------------------------------------------------------------MIQNFRVYYRDSRD--PVWKGPAKLL-----------------------WKGEGAV----VIQDNSDIKVVPRRKAK---IIRD---------------"

aa2 =
  "-----------------------------------------------------------------------------------------------------------ADRKLCADQECSHPIS------------MAVALQDYMAPDCRFLTIHRGQVVYV------------FSKLKGRGRLFWGGSVQGDY--YGDLAARLGYFPSSIVREDQTLKPGKVDVKTDKWDFYCQ"

aa3 =
  "-----------------------------------------------------------------------------------------------------------------GSSGSSGEIA--------------QVTSAYVASGSEQLSLAPGQLILI---------------LKKNTSGWWQGELQA-----RGKKRQKGWFPASHVK---LLGPSSERASGPSSG----"

aa4 =
  "------------------------------------------------------------------------------------------------------------------------NLF--------------VALYDFVASGDNTLSITKGEKLRV---------------LGYNHNGEWCEAQTK---------NGQGWVPSNYIT---PVNS---------------"

aa5 =
  "---------------------------------------------------------------------------------------------------------------------AEGYQY--------------RALYDYKKEREEDIDLHLGDILTVNKGSLVALGFSDGQEARPEEIGWLNGYNET--------TGERGDFPGTYVE---YIGRKKISPP---------"

aa6 =
  "PLALLLDSSLEGEFDLVQRIIYEVDDPSLPNDEGITALHNAVCAGHTEIVKFLVQFGVNVNAADSDGWTPLHCAASCNNVQVCKFLVESGAAVFAMTYSDMQTAADKCEEMEEGYTQCSQFLYGVQEKMGIMNKGVIYALWDYEPQNDDELPMKEGDCMTI------------IHREDEDEIEWWWARLND----------KEGYVPRNLLG----LYP---------------"

aa7 =
  "--------------------------------------------------------------------------------------------------------------TTGRLDLPPGFMF------------KVQAQHDYTATDTDELQLKAGDVVLV----------IPFQNPEEQDEGWLMGVKESDWNQHKELEKCRGVFPENFTE---RVQ----------------"
