# Kono-et-al.-2023
ℹ We do not guarantee the integrity of these programs in their existing versions. We also assume no responsibility for any disadvantages that may result from the reader's use of these programs.

## Image analysis pipelines for CellProfiler ver. 3.1.9
These contain several modifications to [Kono-et-al.-2022](http://doi.org/10.1083/jcb.202201024)[^1] by Arata Komatsubara, which measure the mean fluorescence intensity in the nucleoplasm and the nuclear lamina.

**NP_NL_intensity_fixed.cppipe** recognizes nuclei from DNA images. 3 pixels from the outside of the nucleus are regarded as the nuclear lamina. More than 10 pixels inside from the rim of the nucleus are regarded as the nucleoplasm.  
When nuclei cannot be detected properly using Hoechst images, use **NP_NL_intensity_fixed_EX.cppipe** which switches from `two classes` to `three classes` thresholding distinguishing brighter foci and nuclei from cytoplasm and background.

**NP_NL_intensity_live.cppipe** recognizes nuclei from NLS-sfCherry images. 3 pixels from the outside of the nucleus are regarded as the nuclear lamina. More than 10 pixels inside from the rim of the nucleus are regarded as the nucleoplasm.  
When nuclei cannot be detected properly using NLS images, use **NP_NL_intensity_live_EX.cppipe** which switches from `two classes` to `three classes` thresholding distinguishing brighter foci and nuclei from cytoplasm and background.

## Statistic codes for R ver. 4.2.2
These contain Games-Howell post-hoc multiple comparison test. Used in R Commander ver. 2.8-0 with EZR plugin ver. 1.61.

**Aoki_all.R** is forked from http://aoki2.si.gunma-u.ac.jp/R/src/all.R (last updated at Feb 01, 2019). Code is programmed by Shigenobu AOKI (Professor Emeritus, Gunma University). Encoding has been changed from the original `EUC-JP` to `UTF-8` for R 4.2.0 and later versions.  
This program is provided as a draft for users to freely rewrite and use. I encourage readers to refer to his book[^2] or [website](http://aoki2.si.gunma-u.ac.jp/R/) (in Japanese, but easy to translate) for further information and cite his work appropriately in accordance with international academic publishing guidelines to use this.

**Hatano_socialStatisticsBasic.R** is forked from http://kyoto-edu.sakura.ne.jp/weblesson/statistics/socialStatisticsBasic.R (Ver.1.3; Nov 24, 2021). Code is programmed by Shinsuke HATANO (Ryukoku University).  
I have obtained permission to redistribute these codes on GitHub. I encourage readers to refer to [his website](http://kyoto-edu.sakura.ne.jp/ryukoku/?&course=statistics&content=R01_source) (in Japanese, but easy to translate) for further information. In accordance with international academic publishing guidelines, appropriate citation should be provided to use this.

ℹ To perform the Games-Howell test using R software, the **rosetta**, **PMCMRplus** and ***rstatix*** packages are now available. (The **userfriendlyscience** package was provided for older versions than R 4.0.5).  
* However, the **rosetta** package does not show the digit after p=0.05 even if p=0.055, so it is not possible to determine whether p<0.05 or p>0.05.  
* **PMCMRplus** shows p=0.055 as it displays up to three decimal places, but p=0.0499 is shown as p=0.050, so it is unclear whether p<0.05 or not.  
* ***rstatix*** shows asterisks to indicate significant differences, but up to two decimal places (e.g. p=0.05*, even if p=0.0499).  
* On the other hand, **Aoki's** and **Hatano's codes** can be calculated to eight decimal places (e.g. p=0.04994141), which allow determination of p<0.05 and indicate accurate p-values.

## Modified EZR ver. 1.61a for Windows
This contains several modifications by Yohei Kono to the [EZR ver. 1.61](http://www.jichi.ac.jp/saitama-sct/SaitamaHP.files/statmedEN.html) programmed by Yoshinobu Kanda, which is a customized plugin of the [R Commander ver. 2.8-0](http://socialsciences.mcmaster.ca/jfox/Misc/Rcmdr/) developed by John Fox[^3][^4][^5].
These patches enable user-friendly manipulation of Games-Howell tests on **Aoki_all.R**, **Hatano_socialStatisticsBasic.R** and ***rstatix*** package, and other useful statistic analyses.  
If you use the EZR in a scientific paper, you should cite his report[^6] as a reference.

ℹ If there are 10 more results from the Games-Howell test on ***rstatix*** package like below,  
> "… with 18 more rows"  
> "Use `print(n = ...)` to see more rows"  

, then write `print (res, n = 28)` in R script window and run it.

### References
[^1]:Kono Y, et al. (2022). Nucleoplasmic lamin C rapidly accumulates at sites of nuclear envelope rupture with BAF and cGAS. [Journal of Cell Biology, 221(12), e202201024.](http://doi.org/10.1083/jcb.202201024)  
[^2]:Aoki S (2009). [Statistical analysis by R] R ni yoru Toukei Kaiseki (in Japanese), [Ohmsha, Ltd. Tokyo Japan](http://shop.ohmsha.co.jp/shopdetail/000000001819/)  
[^3]:Fox J (2005). “The R Commander: A Basic Statistics Graphical User Interface to R.” [Journal of Statistical Software, 14(9), 1–42.](http://doi.org/10.18637/jss.v014.i09)  
[^4]:Fox J (2017). Using the R Commander: A Point-and-Click Interface for R. [Chapman and Hall/CRC Press, Boca Raton FL.](http://socialsciences.mcmaster.ca/jfox/Books/RCommander/)  
[^5]:Fox J, Bouchet-Valat M (2022). Rcmdr: R Commander. R package version 2.8-0, http://socialsciences.mcmaster.ca/jfox/Misc/Rcmdr/  
[^6]:Kanda Y (2013). Investigation of the freely available easy-to-use software ‘EZR’ for medical statistics. [Bone Marrow Transplantation, 48,452-458.](http://doi.org/10.1038/bmt.2012.244)
