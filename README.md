# MOVA
a method for evaluating the pathogenicity of missense variants using AlphaFold2
***
MOVA algorithms are described in detail in:  
Yuya Hatano, Tomohiko Ishihara,  Osamu Onodera. Accuracy of a machine learning method based on structural and locational information from AlphaFold2 for predicting the pathogenicity of TARDBP and FUS gene variants in ALS. bioRxiv 2022.07.07.499092; doi: https://doi.org/10.1101/2022.07.07.499092
## installation
- MOVA has been confirmed to work with R 4.2.2.
- devtools and LowMACA installation is required.

```
install.packages("devtools", dependencies = TRUE)
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("LowMACA")
```

- To install MOVA, execute the following two lines

```
library(devtools)
install_github("yuya-hatano/MOVA")
```
***
## Tutorial

- To load the package, enter the following in your R session:
```
library(LowMACA)
library(MOVA)
```
- Create a csv file from the HGMD data by manually creating four columns as shown in the figure below.
- Specifically, based on HGMD data, list amino acid changes (column: HGVSprotein), VCF (column: VCF), DM if DM in HGMD, DMp if DM? (column: possible), Pathogenic if pathogenic but not a target phenotype, and Target if target phenotype(column: Type) for each variant as shown in the figure below and save it. 
![readme_pic1](https://user-images.githubusercontent.com/108056381/226163655-e107e999-751a-4c11-8eb8-a7b7e8afc4a6.png)
- Perform the following in R.
```
Hgmd_divide("created csv file")
```
- Search for target genes in gnomAD v3.1.2.
- Select Ensembl canonical transcript.
![readme_pic2](https://user-images.githubusercontent.com/108056381/226165662-1cf3d48a-5c20-44a6-a800-8da242bccdb1.png)
- Please select only Missense / inframe indel.
- Deselect Indels and select SNVs.
![readme_pic3](https://user-images.githubusercontent.com/108056381/226166217-7bf41609-458e-4cf6-9717-cfcf4b1bc3d2.png)
- Select 'Export variants to CSV' and download.
![readme_pic4](https://user-images.githubusercontent.com/108056381/226166553-c767cbeb-64f7-405c-a172-cc2ad332a5df.png)
- Manually remove all but the missense variants from the downloaded file.
