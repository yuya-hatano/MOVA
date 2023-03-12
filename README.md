# MOVA
a method for evaluating the pathogenicity of missense variants using AlphaFold2
***
MOVA algorithms are described in detail in:  
Yuya Hatano, Tomohiko Ishihara,  Osamu Onodera. Accuracy of a machine learning method based on structural and locational information from AlphaFold2 for predicting the pathogenicity of TARDBP and FUS gene variants in ALS. bioRxiv 2022.07.07.499092; doi: https://doi.org/10.1101/2022.07.07.499092
## installation
- MOVA has been confirmed to work with R 4.2.2.
- data.table, ggplot2, LowMACA, randomForest, ROCR, bio3d, seqinr, dplyr, tidyverse, cvAUC, splitstackshape, kernlab, xgboost, Matrix, gridExtra installation is required.

***
```
install.packages("data.table")
install.packages("ggplot2")
install.packages("LowMACA")
install.packages("randomForest")
install.packages("ROCR")
install.packages("bio3d")
install.packages("seqinr")
install.packages("dplyr")
install.packages("tidyverse")
install.packages("cvAUC")
install.packages("splitstackshape")
install.packages("kernlab")
install.packages("xgboost")
install.packages("Matrix")
install.packages("gridExtra")
```
***

