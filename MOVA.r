############################################################
# MOVA(a method of variant pathogenicity using alphafold2) #
############################################################

# Requires the data.table, ggplot2, LowMACA, randomForest, ROCR, bio3d, seqinr, dplyr
#  packages  (needs to be installed first)
library(data.table)
library(ggplot2)
library(LowMACA)
library(randomForest)
library(ROCR)
library(bio3d)
library(seqinr)
library(dplyr)

MOVA <- function(fasta_file_name, protein_name,  pdb_file_name, variant_file){
# Import BLOSUM62
 data("BLOSUM62")
 amino_code1 <- c(
    'A', 'R', 'N', 'D',
    'C', 'Q', 'E', 'G',
    'I', 'L', 'K', 'M',
    'F', 'P', 'S', 'T',
    'W', 'Y', 'V', 'H')
# Import fasta file
 Protein <- read.fasta(file=fasta_file_name, seqtype = "AA")
 Proteinr <- Protein$Protein[1:length(Protein$Protein)]
# Import pdb file
 Protein3D <- read.pdb(pdb_file_name)
# Get information related to atom records in pdb file
 Protein3Dr <- as.data.frame(Protein3D$atom) 
# Summarize the atom's information by amino acid residue.
 Protein3Dp <- Protein3Dr %>%
  group_by(resno) %>%
  summarise(x = mean(x),y = mean(y), z = mean(z), b = mean(b))
# Save the data of atoms in the pdb file summarized by amino acid residue.
 fwrite(Protein3Dp,  paste(protein_name, "_3D.csv", sep = ""))
# Load the list of variants
 D <- fread(variant_file)
 D <- D[(D$possible!="DMp")|(is.na(D$possible)), ]
#
 D <- D[(D$Type=="ALSFTD")|(D$Type=="Ctrl"), ]
#
# Display the number of variants.
 paste("positive_variant(n)", nrow(D[D$Type=="ALSFTD",]), sep = ":")
 paste("negative_variant(n)", nrow(D[D$Type=="Ctrl",]), sep = ":")
#Extract the variant positions and substituted residues.
 D[, possible:=NULL]
 D$Pos <- as.integer(substr(D$change, 2, nchar(D$change)-1))
 D$ss <- substr(D$change, nchar(D$change), nchar(D$change))
# Initialize data.frame.
 all_aa <- data.frame(matrix(rep(NA, 10), nrow=1))[numeric(0), ]
 colnames(all_aa) <- c("change", "Type", "Pos", "ss", "con", "x", "y", "z", "b", "Type2")
 fpredict <- data.frame(matrix(rep(NA, 3), nrow=1))[numeric(0), ]
 colnames(fpredict) <- c("change", "predict", "n")
 imp.rfa <- data.frame(matrix(rep(NA, 3), nrow=1))[numeric(0), ]
 colnames(imp.rfa) <- c("name", "IncNodePurity", "n")
# Obtain information on each variant of the teacher data
 for (i in 1:nrow(D)){
#(BLOSUM62 score of amino acid residues in WT) - (BLOSUM62 score of amino acid residues in MT) (con)
  D$con[i] <- BLOSUM62[Proteinr[D$Pos[i]],Proteinr[D$Pos[i]]]-BLOSUM62[Proteinr[D$Pos[i]],D$ss[i]]
#x, y, z coordinates (x, y, z)
  D$x[i] <- Protein3Dp$x[D$Pos[i]]
  D$y[i] <- Protein3Dp$y[D$Pos[i]]
  D$z[i] <- Protein3Dp$z[D$Pos[i]]
#plddt score (b)
  D$b[i] <- Protein3Dp$b[D$Pos[i]]
 }
# Obtain information on each variant for final predict
 for (i in 1:length(Proteinr)){
  for(i2 in 1:length(amino_code1)){
   aa <- data.frame( change =  paste(Proteinr[i], i, amino_code1[i2], sep = ""), Type = "unknown", Pos = i, ss = amino_code1[i2], predict = 0, con = BLOSUM62[Proteinr[i],Proteinr[i]]-BLOSUM62[Proteinr[i],amino_code1[i2]], x = Protein3Dp$x[i], y = Protein3Dp$y[i], z = Protein3Dp$z[i], b = Protein3Dp$b[i], Type2 = 3)
   all_aa <- rbind(all_aa, aa)
  }
 }

 df <- data.frame(matrix(rep(NA, 4), nrow=1))[numeric(0), ]
 colnames(df) <- c("ID", "result", "predict",  "change")
 D$label <- 1:nrow(D)
#randomForest
 for (i in 1:nrow(D)){
#leave-one-out. Separate training data (dat) from test data (val)
  dat <- D[D$label != i,]
  val <- D[D$label == i,]
  val_type <- D[D$label == i,]$Type2
#Model construction with randomForest.
#Type2 is the objective variable and con+x+y+z+b is the explanatory variable.
  rf<-randomForest(Type2~con+x+y+z+b,dat, importance=TRUE)
#Measure variable importance.
  imp.rf <- importance(rf)
#Store the predicted values of the test data in data.frame.
  f <-data.frame(ID = val$ID, result = val$Type2, predict = predict(rf,val), change = val$change)
  df <- rbind(df, f)
#Store the variable importance in data.frame.
  f <- data.frame(name = rownames(imp.rf), IncNodePurity = imp.rf, n = i)
  imp.rfa <- rbind(imp.rfa,f)
 }
#Store the predicted values of all variants for final predict in data.frame.
 for (i in 1:100){
  rf <-randomForest(Type2~con+x+y+z+b,D)
  f <-data.frame(change = all_aa$change, predict = predict(rf,all_aa), n = i) 
  fpredict <- rbind(fpredict, f)
 }
#Take the average of the predictions to make the final prediction.
 fpredictb <- fpredict %>%
  group_by(change) %>%
  summarise(predict = mean(predict))
#Take the average of variable importance
 imp.rfb <- imp.rfa %>%
  group_by(name) %>%
  summarise(IncMSE = mean(IncNodePurity..IncMSE), IncNodePurity = mean(IncNodePurity.IncNodePurity))
#If the amino acid residues do not change, the predicted value is 0.
 for (i in 1:nrow(fpredictb)){
  if(substr(fpredictb$change[i], nchar(fpredictb$change[i]), nchar(fpredictb$change[i]))==substr(fpredictb$change[i], 1, 1)){
   fpredictb$predict[i] = 0
  }
 }

 fpredictb$ref <- substr(fpredictb$change, 1, 1)
 fpredictb$Pos <- as.integer(substr(fpredictb$change, 2, nchar(fpredictb$change)-1))
 fpredictb$alt <- substr(fpredictb$change, nchar(fpredictb$change), nchar(fpredictb$change))
#Save the final prediction.
 fwrite(fpredictb, paste(protein_name, "finalpredict.csv", sep = ""))
#Save the average of variable importance.
 fwrite(imp.rfb, paste(protein_name, "importance.csv", sep = ""))
#Plot ROC curve
 pred <- prediction(df$predict, df$result)
 perf <- performance(pred, "tpr", "fpr")
 plot(perf, col = "red")
 auc.tmp <- performance(pred,"auc")
 auc <- as.numeric(auc.tmp@y.values)
 tab <- data.frame(Cutoff=unlist(pred@cutoffs),
  TP=unlist(pred@tp), FP=unlist(pred@fp),
  FN=unlist(pred@fn), TN=unlist(pred@tn),
  Sensitivity=unlist(pred@tp)/(unlist(pred@tp)+unlist(pred@fn)),
  Specificity=unlist(pred@tn)/(unlist(pred@fp)+unlist(pred@tn)),
  Accuracy=((unlist(pred@tp)+unlist(pred@tn))/nrow(df)),
  Precision=(unlist(pred@tp)/(unlist(pred@tp)+unlist(pred@fp)))
  )
 tab$Youden <- tab$Sensitivity+tab$Specificity-1
 YI <- tab[order(tab$Youden, decreasing=T), ]
 YI <- cbind(YI[1,], data.frame(AUC = c(auc)), data.frame(positive_variant_num = c(nrow(D[D$Type=="ALSFTD",]))), data.frame(negative_variant_num = c(nrow(D[D$Type=="Ctrl",]))))
 fwrite(YI, paste(protein_name, "_result.csv", sep = ""))
 fwrite(merge(D,df), paste(protein_name, "_predict.csv", sep = ""))
 return(YI)
}

PolyPhen_MOVA <- function(file_name){
 df <- fread(file_name)
 #
 df <- df[(df$Type=="ALSFTD")|(df$Type=="Ctrl"), ]
 #
 df$predict <- df$pph2_prob
 df$result <- df$Type2
 #
 pred <- prediction(df$predict, df$result)
 perf <- performance(pred, "tpr", "fpr")
 plot(perf, col = "orange", add=TRUE)
 auc.tmp <- performance(pred,"auc")
 auc <- as.numeric(auc.tmp@y.values)
 return(auc)
}

AlphScore_MOVA <- function(file_name){
 df <- fread(file_name)
 #
 df <- df[(df$Type=="ALSFTD")|(df$Type=="Ctrl"), ]
 df$predict <- df$AlphScore
 df$result <- df$Type2
 #
 pred <- prediction(df$predict, df$result)
 perf <- performance(pred, "tpr", "fpr")
 plot(perf, col = "green", add=TRUE)
 auc.tmp <- performance(pred,"auc")
 auc <- as.numeric(auc.tmp@y.values)
 return(auc)
}

CADD_MOVA <- function(file_name){
 df <- fread(file_name)
 #A
 df <- df[(df$Type=="ALSFTD")|(df$Type=="Ctrl"), ]
 #
 df$predict <- df$CADD_raw
 df$result <- df$Type2
 #
 pred <- prediction(df$predict, df$result)
 perf <- performance(pred, "tpr", "fpr")
 plot(perf, col = "blue", add=TRUE)
 auc.tmp <- performance(pred,"auc")
 auc <- as.numeric(auc.tmp@y.values)
 return(auc)
}

Com_CADD_MOVA <- function(predict_file_name){
 df <- fread(predict_file_name)
#
 df <- df[(df$Type=="ALSFTD")|(df$Type=="Ctrl"), ]
#
 df$loo_MOVA <- df$predict
 df$result <- df$Type2
#MOVA
 pred <- prediction(df$loo_MOVA, df$result)
 perf <- performance(pred, "tpr", "fpr")
 plot(perf, col = "blue")
 auc.tmp <- performance(pred,"auc")
 auc <- as.numeric(auc.tmp@y.values)
 auc_bind <- data.frame(MOVA = c(auc))
#CADD
 df$predict <- df$CADD_raw
 df$result <- df$Type2
 pred <- prediction(df$predict, df$result)
 perf <- performance(pred, "tpr", "fpr")
 plot(perf, col = "black", add=TRUE)
 auc.tmp <- performance(pred,"auc")
 auc <- as.numeric(auc.tmp@y.values)
 auc_bind <- cbind(auc_bind, data.frame(CADD = c(auc)))
#Alph_Score+CADD
 df$predict <- df$glm_AlphCadd
 df$result <- df$Type2
#
 pred <- prediction(df$predict, df$result)
 perf <- performance(pred, "tpr", "fpr")
 plot(perf, col = "gray", add=TRUE)
 auc.tmp <- performance(pred,"auc")
 auc <- as.numeric(auc.tmp@y.values)
 auc_bind <- cbind(auc_bind, data.frame(com_Alph_CADD = c(auc)))
#MOVA+CADD
 df$label <- 1:nrow(df)
 df2 <- data.frame(matrix(rep(NA, 4), nrow=1))[numeric(0), ]
 colnames(df2) <- c("ID","result", "CADDMOVA", "change")
 for (i in 1:nrow(df)){
  dat <- df[df$label != i,]
  val <- df[df$label == i,]
  glmals <-glm(Type2~loo_MOVA+CADD_raw,data=dat, family=binomial)
  f <-data.frame(ID= val$ID, result = val$Type2, CADDMOVA = predict(glmals,newdata = val), change = val$change)
  df2 <- rbind(df2, f)
 }
 df <- merge(df,df2)
 df$predict <- df$CADDMOVA
 df$result <- df$Type2
#
 pred <- prediction(df$predict, df$result)
 perf <- performance(pred, "tpr", "fpr")
 plot(perf, col = "red", add=TRUE)
 auc.tmp <- performance(pred,"auc")
 auc <- as.numeric(auc.tmp@y.values)
 auc_bind <- cbind(auc_bind, data.frame(com_MOVA_CADD = c(auc)))
 return(auc_bind)
}


REVEL_MOVA <- function(file_name){
 df <- fread(file_name)
#
 df <- df[(df$Type=="ALSFTD")|(df$Type=="Ctrl"), ]
#
 df$predict <- df$REVEL_score
 df$result <- df$Type2
#
 pred <- prediction(df$predict, df$result)
 perf <- performance(pred, "tpr", "fpr")
 plot(perf, col = "black", add=TRUE)
 auc.tmp <- performance(pred,"auc")
 auc <- as.numeric(auc.tmp@y.values)
 return(auc)
}

Com_REVEL_MOVA <- function(predict_file_name){
 df <- fread(predict_file_name)
#
 df <- df[(df$Type=="ALSFTD")|(df$Type=="Ctrl"), ]
#MOVA
 df$loo_MOVA <- df$predict
 df$result <- df$Type2
#
 pred <- prediction(df$loo_MOVA, df$result)
 perf <- performance(pred, "tpr", "fpr")
 plot(perf, col = "blue")
 auc.tmp <- performance(pred,"auc")
 auc <- as.numeric(auc.tmp@y.values)
 auc_bind <- data.frame(MOVA = c(auc))
#REVEL
 df$predict <- df$REVEL_score
 df$result <- df$Type2
#
 pred <- prediction(df$predict, df$result)
 perf <- performance(pred, "tpr", "fpr")
 plot(perf, col = "black", add=TRUE)
 auc.tmp <- performance(pred,"auc")
 auc <- as.numeric(auc.tmp@y.values)
 auc_bind <- cbind(auc_bind, data.frame(REVEL = c(auc)))
#AlphScore+REVEL
 df$predict <- df$glm_AlphRevel
 df$result <- df$Type2
#
 pred <- prediction(df$predict, df$result)
 perf <- performance(pred, "tpr", "fpr")
 plot(perf, col = "gray", add=TRUE)
 auc.tmp <- performance(pred,"auc")
 auc <- as.numeric(auc.tmp@y.values)
 auc_bind <- cbind(auc_bind, data.frame(com_Alph_REVEL = c(auc)))
#MOVA+REVEL
 df$label <- 1:nrow(df)
 df2 <- data.frame(matrix(rep(NA, 4), nrow=1))[numeric(0), ]
 colnames(df2) <- c("ID","result", "REVELMOVA", "change")
 for (i in 1:nrow(df)){
  dat <- df[df$label != i,]
  val <- df[df$label == i,]
  glmals <-glm(Type2~loo_MOVA+REVEL_score,data=dat, family=binomial)
  f <-data.frame(ID= val$ID, result = val$Type2, REVELMOVA = predict(glmals,newdata = val), change = val$change)
  df2 <- rbind(df2, f)
 }
 df <- merge(df,df2)
 df$predict <- df$REVELMOVA
 df$result <- df$Type2
#
 pred <- prediction(df$predict, df$result)
 perf <- performance(pred, "tpr", "fpr")
 plot(perf, col = "red", add=TRUE)
 auc.tmp <- performance(pred,"auc")
 auc <- as.numeric(auc.tmp@y.values)
 auc_bind <- cbind(auc_bind, data.frame(com_MOVA_REVEL = c(auc)))
 return(auc_bind)
}

EVE_MOVA <- function(file_name, EVE_file_name){
 df <- fread(file_name)
 y <- fread(EVE_file_name)
#
 df <- df[(df$Type=="ALSFTD")|(df$Type=="Ctrl"), ]
#
 df[, EVE:=NA]
 df[, EVE_class:=NA]
 df$Pos <- as.integer(substr(df$change, 2, nchar(df$change)-1))
 for (i in 1:nrow(df)){
  df$EVE[i] <- y[(y$position==df$Pos[i])&(y$mt_aa==df$aaalt[i])]$EVE_scores_ASM
  df$EVE_class[i] <- y[(y$position==df$Pos[i])&(y$mt_aa==df$aaalt[i])]$EVE_classes_75_pct_retained_ASM
 }
 df$predict <- df$EVE
 df$result <- df$Type2
 df <- subset(df, !(is.na(df$predict)))
#
 pred <- prediction(df$predict, df$result)
 perf <- performance(pred, "tpr", "fpr")
 plot(perf, col = "gray", add=TRUE)
 auc.tmp <- performance(pred,"auc")
 auc <- as.numeric(auc.tmp@y.values)
 return(auc)
}

#Be sure to check the contents of the file yourself and delete all but the missense mutations.
Edit_gnomAD_file <- function(input_file_name, export_file_name){
 amino_code1 <- c(
    'A', 'R', 'N', 'D',
    'C', 'Q', 'E', 'G',
    'I', 'L', 'K', 'M',
    'F', 'P', 'S', 'T',
    'W', 'Y', 'V', 'H')
 amino_code3 = c(
    'Ala', 'Arg', 'Asn', 'Asp',
    'Cys', 'Gln', 'Glu', 'Gly',
    'Ile', 'Leu', 'Lys', 'Met',
    'Phe', 'Pro', 'Ser', 'Thr',
    'Trp', 'Tyr', 'Val', 'His')
 names(amino_code1)<-amino_code3
 x <- fread(input_file_name)
 colnames(x)[10] <- "Protein_Consequence"
 x$ref <- substr(x$Protein_Consequence, 3, 5)
 x$Pos <- substr(x$Protein_Consequence, 6, nchar(x$Protein_Consequence)-3)
 x$alt <- substr(x$Protein_Consequence, nchar(x$Protein_Consequence)-2, nchar(x$Protein_Consequence))
 for (i in 1:nrow(x)){
  x$change[i] <- paste(amino_code1[x$ref[i]], x$Pos[i], amino_code1[x$alt[i]], sep = "")
 }
 y <- data.frame(Chromosome = x$Chromosome, Position = x$Position, Reference = x$Reference, Alternate = x$Alternate, change = x$change)
 fwrite(y, export_file_name)
}

#Hgmd_files must be created from the HGMD data yourself.
Edit_variant_data <- function(AlphScore_final_file, hgmd_file_name, gnomad_file_name, uniprot_name, gene_name, export_file_name){
#https://zenodo.org/record/6288139
 x <- fread(AlphScore_final_file)
 fwrite(x[x$Uniprot_acc_split==uniprot_name,], paste(gene_name, "_alph.csv", sep = ""))

 x <- fread(hgmd_file_name)
 y <- fread(gnomad_file_name)
 D <- merge(x,y,all=T)
 D <- D[(D$possible!="DMp")|(is.na(D$possible)), ]
 D$Type <- replace(D$Type, which(is.na(D$Type)), "Ctrl")
 for (i in 1:nrow(D)){
  D$change[i] <- replace(D$change[i], which(is.na(D$change[i])), D$HGVSprotein[i])
 }
 D[, HGVSprotein:=NULL]
 fwrite(D, paste(gene_name, "_variant.csv", sep = ""))
 x <- fread(paste(gene_name, "_alph.csv", sep = ""))
 y <- fread(paste(gene_name, "_variant.csv", sep = ""))
 colnames(x)[1] <- "Chromosome"
 colnames(x)[2] <- "Position"
 colnames(x)[3] <- "Reference"
 colnames(x)[4] <- "Alternate"
 D <- merge(x,y)
 D$Type2 <- 0
 D$Type3 <- 0
 D[D$Type=="ALSFTD",]$Type2 <- 1
 D[D$Type=="ALSFTD",]$Type3 <- 1
 D[D$Type=="Pathogenic",]$Type2 <- 0
 D[D$Type=="Pathogenic",]$Type3 <- 1
 fwrite(D, export_file_name)
}

