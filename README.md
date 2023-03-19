# MOVA
a method for evaluating the pathogenicity of missense variants using AlphaFold2
***
MOVA algorithms are described in detail in:  
Yuya Hatano, Tomohiko Ishihara,  Osamu Onodera. Accuracy of a machine learning method based on structural and locational information from AlphaFold2 for predicting the pathogenicity of TARDBP and FUS gene variants in ALS. bioRxiv 2022.07.07.499092; doi: https://doi.org/10.1101/2022.07.07.499092
## installation
- MOVA has been confirmed to work with R 4.2.2.
- devtools installation is required.

***
```
install.packages("devtools", dependencies = TRUE)
```
***

- To install MOVA, execute the following two lines

***
```
library(devtools)
install_github("yuya-hatano/MOVA")
```
***
