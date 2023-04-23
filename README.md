# Kono-et-al.-2023
# Image analysis pipelines for CellProfiler ver. 3.1.9
These contain several modifications to Kono-et-al.-2022 by Arata Komatsubara, which measure the mean fluorescence intensity in the nucleoplasm and the nuclear lamina.

NP_NL_intensity_fixed.cppipe recognizes nuclei from DNA images. 3 pixels from the outside of the nucleus are regarded as the nuclear lamina. More than 10 pixels inside from the rim of the nucleus are regarded as the nucleoplasm. 
When nuclei cannot be detected properly using Hoechst images, use NP_NL_intensity_fixed_EX.cppipe which switches from 'two classes' to 'three classes' thresholding distinguishing brighter foci and nuclei from cytoplasm and background.

NP_NL_intensity_live.cppipe recognizes nuclei from NLS-sfCherry images. 3 pixels from the outside of the nucleus are regarded as the nuclear lamina. More than 10 pixels inside from the rim of the nucleus are regarded as the nucleoplasm. 
When nuclei cannot be detected properly using NLS images, use NP_NL_intensity_live_EX.cppipe which switches from 'two classes' to 'three classes' thresholding distinguishing brighter foci and nuclei from cytoplasm and background.

# Statistic codes for R ver. 4.2.2
These contain Games-Howell post-hoc multiple comparison test. Used in R Commander ver. 2.8-0 with EZR plugin ver. 1.61.

Aoki_all.R is forked from http://aoki2.si.gunma-u.ac.jp/R/src/all.R (last updated at Feb 01, 2019). Code is programmed by Shigenobu AOKI (Professor Emeritus, Gunma University). Encoding has been changed from the original "EUC-JP" to "UTF-8" for R 4.2.0 and later versions.

Hatano_socialStatisticsBasic.R is forked from http://kyoto-edu.sakura.ne.jp/weblesson/statistics/socialStatisticsBasic.R (Ver.1.3; Nov 24, 2021). Code is programmed by Shinsuke HATANO (Ryukoku University).

#To perform the Games-Howell test using R software, the "rosetta", "PMCMRplus" and "rstatix" packages are now available. (The "userfriendlyscience" package was provided in older versions than R 4.0.5).

#However, the rosetta package does not show the digit after p=0.05 even if p=0.055, so it is not possible to determine whether p<0.05 or p>0.05.

#PMCMRplus shows p=0.055 as it displays up to three decimal places, but p=0.0499 is shown as p=0.050, so it is unclear whether p<0.05 or not.

#rstatix shows asterisks to indicate significant differences, but up to two decimal places (e.g. p=0.05*, even if p=0.0499).

#On the other hand, the Aoki and Hatano's codes can be shown to eight decimal places (e.g. p=0.04994141), which allow decisions of p<0.05 and indicate accurate p-values.

# Modified EZR ver. 1.61 for Windows
This contains several modifications by Yohei Kono to the EZR ver. 1.61 programmed by Yoshinobu Kanda, which is a customized plugin of the R Commander ver. 2.8-0 developed by John Fox.
These patches enable user-friendly operation of the Games-Howell test of Aoki_all.R and rstatix package, and other useful statistic analyses.

ℹ Games-Howell test on rstatix package

#If there are 10 more results like below,
"… with 18 more rows"
"Use `print(n = ...)` to see more rows"
, then write `print (res, n = 28)` in R script window and run it.
