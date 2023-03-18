##########################################################################################
# MOVA((a method for evaluating the pathogenicity of missense variants using AlphaFold2) #
##########################################################################################

# Requires the data.table, ggplot2, LowMACA, randomForest, ROCR, bio3d, seqinr, dplyr, tidyverse, 
#  cvAUC, splitstackshape, kernlab, xgboost, Matrix, gridExtra
#  packages  (needs to be installed first)
library(data.table)
library(ggplot2)
library(LowMACA)
library(randomForest)
library(ROCR)
library(bio3d)
library(seqinr)
library(dplyr)
library(tidyverse)
library(cvAUC)
library(splitstackshape)
library(kernlab)
library(xgboost)
library(Matrix)
library(gridExtra)

MOVA <- function(fasta_file_name, protein_name,  pdb_file_name, final_predict_input_file, variant_file, phenotype = "Target"){
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
# Load the list of variants
 D <- fread(variant_file)
 D$ID <- paste(D$Chromosome, D$Position, D$Reference, D$Alternate, sep=":")
 D <- D[(D$possible!="DMp")|(is.na(D$possible)), ]
#
 if(phenotype == "Target"){
  D <- D[(D$Type=="Target")|(D$Type=="Ctrl"), ]
 }
#
#Extract the variant positions and substituted residues.
 D[, possible:=NULL]
 D$Pos <- as.integer(substr(D$change, 2, nchar(D$change)-1))
 D$ss <- substr(D$change, nchar(D$change), nchar(D$change))
# Initialize data.frame.
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
 final_data <- fread(final_predict_input_file)
 final_data$ID <- paste(final_data$Chromosome, final_data$Position, final_data$Reference, final_data$Alternate, final_data$aaref, final_data$Pos, final_data$aaalt,  sep=":")
 for (i in 1:nrow(final_data)){
  final_data$con[i] <- BLOSUM62[Proteinr[final_data$Pos[i]],Proteinr[final_data$Pos[i]]]-BLOSUM62[Proteinr[final_data$Pos[i]],final_data$aaalt[i]]
#x, y, z coordinates (x, y, z)
  final_data$x[i] <- Protein3Dp$x[final_data$Pos[i]]
  final_data$y[i] <- Protein3Dp$y[final_data$Pos[i]]
  final_data$z[i] <- Protein3Dp$z[final_data$Pos[i]]
#plddt score (b)
  final_data$b[i] <- Protein3Dp$b[final_data$Pos[i]]
 }

 df <- data.frame(matrix(rep(NA, 6), nrow=1))[numeric(0), ]
 colnames(df) <- c("ID", "result", "predict",  "change", "iter1", "iter2")

 if(phenotype == "Target"){
  D$result <- D$Type2
 } else {
  D$result <- D$Type3
 }
#randomForest
 oldw <- getOption("warn")
 options(warn=-1)
 D <- rowid_to_column(D, var = "id2")
 for(i2 in 1:5){## split data into train and test using stratified sampling
  D$id3 <- D$id2
  j <- 5
  dat1 <- D %>% stratified(., group = "result", size = 1/j)
  dat1$iter <- 1
  datzan <- D[-dat1$id2, ]
  for(i in 2:j){
   datzan[ , c('id2')] <- list(NULL)
   datzan <- rowid_to_column(datzan, var = "id2")
   if(i != j){
    dat2 <- datzan %>% stratified(., group = "result", size = 1/(j-i+1))
   } else {
    dat2 <- datzan
   }
   dat2$iter <- i
   datzan <- datzan[-dat2$id2, ]
   dat1 <- rbind(dat1, dat2)
  }
  dat1[ , c('id2')] <- list(NULL)
  dat1 <- rowid_to_column(dat1, var = "id2")
  for (i in 1:j){
   val <- dat1[dat1$iter == i,]
   train.data <- dat1[-val$id2,]
#Model construction with randomForest.
#Type2 is the objective variable and con+x+y+z+b is the explanatory variable.
   rf<-randomForest(result~con+x+y+z+b,train.data, importance=TRUE)
#Measure variable importance.
   imp.rf <- importance(rf)
#Store the predicted values of the test data in data.frame.
    f <-data.frame(ID = val$ID, result = val$result, predict = predict(rf,val), change = val$change, iter1 = i, iter2 =i2)
   df <- rbind(df, f)
#Store the variable importance in data.frame.
   f <- data.frame(name = rownames(imp.rf), IncNodePurity = imp.rf, n = i)
   imp.rfa <- rbind(imp.rfa,f)
  }
 }
#Store the predicted values of all variants for final predict in data.frame.
 for (i in 1:30){
   rf <-randomForest(result~con+x+y+z+b,D)
  f <-data.frame(ID = final_data$ID, predict = predict(rf,final_data), n = i) 
  fpredict <- rbind(fpredict, f)
 }
 options(warn = oldw)
#Take the average of the predictions to make the final prediction.
 fpredictb <- fpredict %>%
  group_by(ID) %>%
  summarise(MOVA_predict = mean(predict))
#Take the average of variable importance
 imp.rfb <- imp.rfa %>%
  group_by(name) %>%
  summarise(IncMSE = mean(IncNodePurity..IncMSE), IncNodePurity = mean(IncNodePurity.IncNodePurity))
#If the amino acid residues do not change, the predicted value is 0.
 ID <- str_split(fpredictb$ID, pattern = ":", simplify = TRUE)
 for (i in 1:nrow(fpredictb)){
  if(ID[i, 5]==ID[i, 7]){
   fpredictb$predict[i] = 0
  }
 }

 final_data <- merge(fpredictb, final_data)
#Save the final prediction.
 fwrite(final_data, paste(protein_name, phenotype, "finalpredict.csv", sep = "_"))
#Save the average of variable importance.
 fwrite(imp.rfb, paste(protein_name, phenotype, "importance.csv", sep = "_"))
#Plot ROC curve
 df$iter3 <- df$iter1 + (df$iter2 * 5)
 out <- cvAUC(df$predict, df$result, label.ordering = NULL, folds = df$iter3)
 plot(out$perf, col = "red", avg = "vertical")
 cvauc <- out$cvAUC
 YI2 <- data.frame(matrix(rep(NA, 5), nrow=1))[numeric(0), ]
 colnames(YI2) <- c("Cutoff", "positive_variant_num", "negative_variant_num", "AUC", "cvauc")
 for(i in 6:30){
  pred <- prediction(df[df$iter3==i,]$predict, df[df$iter3==i,]$result)
  auc.tmp <- performance(pred,"auc")
  auc <- as.numeric(auc.tmp@y.values)
  tab <- data.frame(Cutoff=unlist(pred@cutoffs),
   TP=unlist(pred@tp), FP=unlist(pred@fp),
   FN=unlist(pred@fn), TN=unlist(pred@tn),
   Sensitivity=unlist(pred@tp)/(unlist(pred@tp)+unlist(pred@fn)),
   Specificity=unlist(pred@tn)/(unlist(pred@fp)+unlist(pred@tn)),
   Accuracy=((unlist(pred@tp)+unlist(pred@tn))/nrow(df)),
   Precision=(unlist(pred@tp)/(unlist(pred@tp)+unlist(pred@fp))))
  tab$Youden <- tab$Sensitivity+tab$Specificity-1
  YI <- tab[order(tab$Youden, decreasing=T), ]
  YI <- data.frame(Cutoff = YI$Cutoff[1], positive_variant_num = c(nrow(df[(df$result==1)&(df$iter2==1),])), negative_variant_num = c(nrow(df[(df$result==0)&(df$iter2==1),])),AUC = c(auc),cvauc = c(cvauc))
  YI2 <- rbind(YI2, YI)
 }
 fwrite(YI2, paste(protein_name, phenotype, "result.csv", sep = "_"))
 df2 <- df %>%
  group_by(ID) %>%
  summarise(MOVA_predict = mean(predict))
 D[ , c('result')] <- list(NULL)
 fwrite(merge(D,df2), paste(protein_name, phenotype, "predict.csv", sep = "_"))
 fwrite(df, paste(protein_name, phenotype, "predict_orig.csv", sep = "_"))
 return(YI2)
}

MOVA_LOPOV<- function(phenotype = "Target"){
 read_name <- c(
   "TARDBP",
   "FUS",
   "SETX",
   "TBK1",
   "OPTN",
   "SOD1",
   "VCP",
   "SQSTM1",
   "ANG",
   "UBQLN2",
   "DCTN1",
   "CCNF")
 D <- fread(paste(read_name[1], phenotype, "predict.csv", sep="_"))
 final_data <- fread(paste(read_name[1], phenotype, "finalpredict.csv", sep="_"))
 D$gene <- read_name[1]
 D$gene_number <- 1
 D$glm_CaddDeogenRevel <- 0
 D[, glm_CaddDeogenRevel:=NULL]
 final_data$gene <- read_name[1]
 final_data$gene_number <- 1
 final_data$glm_CaddDeogenRevel <- 0
 final_data[, glm_CaddDeogenRevel:=NULL]
 D$t_distance <- 0
 if(phenotype == "Target"){
  D$result <- D$Type2
 } else {
  D$result <- D$Type3
 }
 x <- D[D$result==1,]$x
 y <- D[D$result==1,]$y
 z <- D[D$result==1,]$z
 for(i3 in 1:nrow(D)){
    t_distance <- sqrt((x-D[i3]$x)*(x-D[i3]$x)+(y-D[i3]$y)*(y-D[i3]$y)+(z-D[i3]$z)*(z-D[i3]$z))
    t_distance <- sort(t_distance,decreasing=F)
    if(D[i3]$result == 1){
       D[i3]$t_distance <- t_distance[2]
    } else {
       D[i3]$t_distance <- t_distance[1]
    }

 }
 for(i in 2:12){
  df2 <- fread(paste(read_name[i], phenotype, "predict.csv", sep="_"))
  final_data2 <-  fread(paste(read_name[i], phenotype, "finalpredict.csv", sep="_"))
  df2$gene <- read_name[i]
  df2$gene_number <- i
  df2$glm_CaddDeogenRevel <- 0
  df2[, glm_CaddDeogenRevel:=NULL]
  if(phenotype == "Target"){
   df2$result <- df2$Type2
  } else {
   df2$result <- df2$Type3
  }

  df2$t_distance <- 0
  x <-df2[df2$result==1,]$x
  y <- df2[df2$result==1,]$y
  z <- df2[df2$result==1,]$z
  for(i3 in 1:nrow(df2)){
    t_distance <- sqrt((x-df2[i3]$x)*(x-df2[i3]$x)+(y-df2[i3]$y)*(y-df2[i3]$y)+(z-df2[i3]$z)*(z-df2[i3]$z))
    t_distance <- sort(t_distance,decreasing=F)
    if(df2[i3]$result == 1){
       df2[i3]$t_distance <- t_distance[2]
    } else {
       df2[i3]$t_distance <- t_distance[1]
    }

  }

  D <- rbind(D, df2)
  final_data2$gene <- read_name[i]
  final_data2$gene_number <- i
  final_data2$glm_CaddDeogenRevel <- 0
  final_data2[, glm_CaddDeogenRevel:=NULL]
  final_data <- rbind(final_data, final_data2)
 }
 df <- data.frame(matrix(rep(NA, 6), nrow=1))[numeric(0), ]
 colnames(df) <- c("ID", "result", "predict",  "change", "gene_number", "gene")
 fpredict <- data.frame(matrix(rep(NA, 3), nrow=1))[numeric(0), ]
 colnames(fpredict) <- c("change", "predict", "n")
 imp.rfa <- data.frame(matrix(rep(NA, 3), nrow=1))[numeric(0), ]
 colnames(imp.rfa) <- c("name", "IncNodePurity", "n")

#randomForest
 oldw <- getOption("warn")
 options(warn=-1)
  for (i in 1:12){

   val <- D[D$gene_number == i,]
   train.data <- D[D$gene_number != i,]
#Model construction with randomForest.
#Type2 is the objective variable and con+x+y+z+b is the explanatory variable.
    rf<-randomForest(result~con+t_distance+b,train.data, importance=TRUE)
#Measure variable importance.
   imp.rf <- importance(rf)
#Store the predicted values of the test data in data.frame.
    f <-data.frame(ID = val$ID, result = val$result, predict = predict(rf,val), change = val$change, gene_number = i, gene =val$gene)
   df <- rbind(df, f)
   f <-data.frame(ID = final_data[final_data$gene_number==i,]$ID, MOVA_LOPOV_predict  = predict(rf,final_data[final_data$gene_number==i,]), gene_number = i) 
   fpredict <- rbind(fpredict, f)
#Store the variable importance in data.frame.
   f <- data.frame(name = rownames(imp.rf), IncNodePurity = imp.rf, n = i)
   imp.rfa <- rbind(imp.rfa,f)
  }

 options(warn = oldw)
#Take the average of the predictions to make the final prediction.
 fpredictb <- fpredict
#Take the average of variable importance
 imp.rfb <- imp.rfa %>%
  group_by(name) %>%
  summarise(IncMSE = mean(IncNodePurity..IncMSE), IncNodePurity = mean(IncNodePurity.IncNodePurity))

 final_data <- merge(fpredictb, final_data)
#Save the final prediction.
 fwrite(final_data, paste("LOPOV", phenotype, "finalpredict.csv", sep = "_"))
#Save the average of variable importance.
 fwrite(imp.rfb, paste("LOPOV", phenotype, "importance.csv", sep = "_"))
#Plot ROC curve
 out <- cvAUC(df$predict, df$result, label.ordering = NULL, folds = df$gene_number)
 plot(out$perf, col = "red", avg = "vertical")
 cvauc <- out$cvAUC
 YI2 <- data.frame(matrix(rep(NA, 5), nrow=1))[numeric(0), ]
 colnames(YI2) <- c("Cutoff", "positive_variant_num", "negative_variant_num", "AUC", "cvauc")
 for(i in 1:12){
  pred <- prediction(df[df$gene_number==i,]$predict, df[df$gene_number==i,]$result)
  auc.tmp <- performance(pred,"auc")
  auc <- as.numeric(auc.tmp@y.values)
  tab <- data.frame(Cutoff=unlist(pred@cutoffs),
   TP=unlist(pred@tp), FP=unlist(pred@fp),
   FN=unlist(pred@fn), TN=unlist(pred@tn),
   Sensitivity=unlist(pred@tp)/(unlist(pred@tp)+unlist(pred@fn)),
   Specificity=unlist(pred@tn)/(unlist(pred@fp)+unlist(pred@tn)),
   Accuracy=((unlist(pred@tp)+unlist(pred@tn))/nrow(df)),
   Precision=(unlist(pred@tp)/(unlist(pred@tp)+unlist(pred@fp))))
  tab$Youden <- tab$Sensitivity+tab$Specificity-1
  YI <- tab[order(tab$Youden, decreasing=T), ]
  YI <- data.frame(Cutoff = YI$Cutoff[1], positive_variant_num = c(nrow(df[(df$result==1)&(df$gene_number==1),])), negative_variant_num = c(nrow(df[(df$result==0)&(df$gene_number==1),])),AUC = c(auc),cvauc = c(cvauc), gene = c(read_name[i]))
  YI2 <- rbind(YI2, YI)
 }
 fwrite(YI2, paste("LOPOV", phenotype, "result.csv", sep = "_"))
 df2 <- df %>%
  group_by(ID) %>%
  summarise(MOVA_LOPOV_predict = mean(predict))
 fwrite(merge(D,df2), paste("LOPOV", phenotype, "predict.csv", sep = "_"))
 fwrite(df, paste("LOPOV", phenotype, "predict_orig.csv", sep = "_"))
 return(YI2)
 
}

MOVA_redraw <- function(MOVA_predict_orig_file, col, add=FALSE){
 df <- fread(MOVA_predict_orig_file)
 out <- cvAUC(df$predict, df$result, label.ordering = NULL, folds = df$iter3)
 plot(out$perf, col = col, avg = "vertical", add=add)
}

MOVA_without_3d_coordinates <- function(protein_name, MOVA_final_predict_file, MOVA_predict_file, phenotype = "Target"){
 D <- fread(MOVA_predict_file)
 D$without_3d_coordinates_predict <- 0
 D[ , c('without_3d_coordinates_predict')] <- list(NULL)
#
 if(phenotype == "Target"){
  D <- D[(D$Type=="Target")|(D$Type=="Ctrl"), ]
 }
#
# Initialize data.frame.
 fpredict <- data.frame(matrix(rep(NA, 3), nrow=1))[numeric(0), ]
 colnames(fpredict) <- c("change", "predict", "n")
# Obtain information on each variant for final predict
 final_data <- fread(MOVA_final_predict_file)
 final_data$without_3d_coordinates_predict<- 0
 final_data[ , c('without_3d_coordinates_predict')] <- list(NULL)
 df <- data.frame(matrix(rep(NA, 6), nrow=1))[numeric(0), ]
 colnames(df) <- c("ID", "result", "predict",  "change", "iter1", "iter2")
 D$id2 <- 0
 D[ , c('id2')] <- list(NULL)
 if(phenotype == "Target"){
  D$result <- D$Type2
 } else {
  D$result <- D$Type3
 }
 oldw <- getOption("warn")
 options(warn=-1)
 D <- rowid_to_column(D, var = "id2")
 for(i2 in 1:5){## split data into train and test using stratified sampling
  D$id3 <- D$id2
  j <- 5
  dat1 <- D %>% stratified(., group = "result", size = 1/j)
  dat1$iter <- 1
  datzan <- D[-dat1$id2, ]
  for(i in 2:j){
   datzan[ , c('id2')] <- list(NULL)
   datzan <- rowid_to_column(datzan, var = "id2")
   if(i != j){
    dat2 <- datzan %>% stratified(., group = "result", size = 1/(j-i+1))
   } else {
    dat2 <- datzan
   }
   dat2$iter <- i
   datzan <- datzan[-dat2$id2, ]
   dat1 <- rbind(dat1, dat2)
  }
  dat1[ , c('id2')] <- list(NULL)
  dat1 <- rowid_to_column(dat1, var = "id2")
  for (i in 1:j){
   val <- dat1[dat1$iter == i,]
   train.data <- dat1[-val$id2,]#Model construction with randomForest.
#Type2 is the objective variable and con+x+y+z+b is the explanatory variable.
    rf<-randomForest(result~con+b,train.data, importance=TRUE)
#Store the predicted values of the test data in data.frame
    f <-data.frame(ID = val$ID, result = val$result, predict = predict(rf,val), change = val$change, iter1 = i, iter2 =i2)
   df <- rbind(df, f)
  }
 }
#Store the predicted values of all variants for final predict in data.frame.
 for (i in 1:30){
   rf <-randomForest(result~con+b,D)
  f <-data.frame(ID = final_data$ID, predict = predict(rf,final_data), n = i) 
  fpredict <- rbind(fpredict, f)
 }
 options(warn = oldw)
#Take the average of the predictions to make the final prediction.
 fpredictb <- fpredict %>%
  group_by(ID) %>%
  summarise(without_3d_coordinates_predict = mean(predict))
#If the amino acid residues do not change, the predicted value is 0.
 ID <- str_split(fpredictb$ID, pattern = ":", simplify = TRUE)
 for (i in 1:nrow(fpredictb)){
  if(ID[i, 5]==ID[i, 7]){
   fpredictb$predict[i] = 0
  }
 }

 final_data <- merge(fpredictb, final_data)
#Save the final prediction.
 fwrite(final_data, MOVA_final_predict_file)
#Plot ROC curve
 df$iter3 <- df$iter1 + (df$iter2 * 5)
 out <- cvAUC(df$predict, df$result, label.ordering = NULL, folds = df$iter3)
 plot(out$perf, col = "black", avg = "vertical", add=TRUE)
 cvauc <- out$cvAUC
 YI2 <- data.frame(matrix(rep(NA, 5), nrow=1))[numeric(0), ]
 colnames(YI2) <- c("Cutoff", "positive_variant_num", "negative_variant_num", "AUC", "cvauc")
 for(i in 6:30){
  pred <- prediction(df[df$iter3==i,]$predict, df[df$iter3==i,]$result)
  auc.tmp <- performance(pred,"auc")
  auc <- as.numeric(auc.tmp@y.values)
  tab <- data.frame(Cutoff=unlist(pred@cutoffs),
   TP=unlist(pred@tp), FP=unlist(pred@fp),
   FN=unlist(pred@fn), TN=unlist(pred@tn),
   Sensitivity=unlist(pred@tp)/(unlist(pred@tp)+unlist(pred@fn)),
   Specificity=unlist(pred@tn)/(unlist(pred@fp)+unlist(pred@tn)),
   Accuracy=((unlist(pred@tp)+unlist(pred@tn))/nrow(df)),
   Precision=(unlist(pred@tp)/(unlist(pred@tp)+unlist(pred@fp))))
  tab$Youden <- tab$Sensitivity+tab$Specificity-1
  YI <- tab[order(tab$Youden, decreasing=T), ]
  YI <- data.frame(Cutoff = YI$Cutoff[1], positive_variant_num = c(nrow(df[(df$result==1)&(df$iter2==1),])), negative_variant_num = c(nrow(df[(df$result==0)&(df$iter2==1),])),AUC = c(auc),cvauc = c(cvauc))
  YI2 <- rbind(YI2, YI)
 }
 fwrite(YI2, paste(protein_name, phenotype, "result_without3d.csv", sep = "_"))
 df2 <- df %>%
  group_by(ID) %>%
  summarise(without_3d_coordinates_predict = mean(predict))
 D[ , c('result')] <- list(NULL)
 fwrite(merge(D,df2), MOVA_predict_file)
 fwrite(df, paste(protein_name, phenotype, "predict_orig_without3d.csv", sep = "_"))
 return(YI2)
}


MOVA_only_pLDDT <- function(protein_name, MOVA_final_predict_file, MOVA_predict_file, phenotype = "Target"){
 D <- fread(MOVA_predict_file)
 D$only_pLDDT_predict<- 0
 D[ , c('only_pLDDT_predict')] <- list(NULL)
#
 if(phenotype == "Target"){
  D <- D[(D$Type=="Target")|(D$Type=="Ctrl"), ]
 }
#
# Initialize data.frame.
 fpredict <- data.frame(matrix(rep(NA, 3), nrow=1))[numeric(0), ]
 colnames(fpredict) <- c("change", "predict", "n")
# Obtain information on each variant for final predict
 final_data <- fread(MOVA_final_predict_file)
 final_data$only_pLDDT_predict<- 0
 final_data[ , c('only_pLDDT_predict')] <- list(NULL)
 df <- data.frame(matrix(rep(NA, 6), nrow=1))[numeric(0), ]
 colnames(df) <- c("ID", "result", "predict",  "change", "iter1", "iter2")
 D$id2 <- 0
 D[ , c('id2')] <- list(NULL)
 if(phenotype == "Target"){
  D$result <- D$Type2
 } else {
  D$result <- D$Type3
 }
 oldw <- getOption("warn")
 options(warn=-1)
 D <- rowid_to_column(D, var = "id2")
 for(i2 in 1:5){## split data into train and test using stratified sampling
  D$id3 <- D$id2
  j <- 5
  dat1 <- D %>% stratified(., group = "result", size = 1/j)
  dat1$iter <- 1
  datzan <- D[-dat1$id2, ]
  for(i in 2:j){
   datzan[ , c('id2')] <- list(NULL)
   datzan <- rowid_to_column(datzan, var = "id2")
   if(i != j){
    dat2 <- datzan %>% stratified(., group = "result", size = 1/(j-i+1))
   } else {
    dat2 <- datzan
   }
   dat2$iter <- i
   datzan <- datzan[-dat2$id2, ]
   dat1 <- rbind(dat1, dat2)
  }
  dat1[ , c('id2')] <- list(NULL)
  dat1 <- rowid_to_column(dat1, var = "id2")
  for (i in 1:j){
   val <- dat1[dat1$iter == i,]
   train.data <- dat1[-val$id2,]#Model construction with randomForest.
#Type2 is the objective variable and con+x+y+z+b is the explanatory variable.
   rf<-randomForest(result~b,train.data, importance=TRUE)
#Store the predicted values of the test data in data.frame.
   f <-data.frame(ID = val$ID, result = val$result, predict = predict(rf,val), change = val$change, iter1 = i, iter2 =i2)
   df <- rbind(df, f)
  }
 }
#Store the predicted values of all variants for final predict in data.frame.
 for (i in 1:30){
   rf <-randomForest(result~b,D)
  f <-data.frame(ID = final_data$ID, predict = predict(rf,final_data), n = i) 
  fpredict <- rbind(fpredict, f)
 }
 options(warn = oldw)
#Take the average of the predictions to make the final prediction.
 fpredictb <- fpredict %>%
  group_by(ID) %>%
  summarise(only_pLDDT_predict = mean(predict))
#If the amino acid residues do not change, the predicted value is 0.
 ID <- str_split(fpredictb$ID, pattern = ":", simplify = TRUE)
 for (i in 1:nrow(fpredictb)){
  if(ID[i, 5]==ID[i, 7]){
   fpredictb$predict[i] = 0
  }
 }

 final_data <- merge(fpredictb, final_data)
#Save the final prediction.
 fwrite(final_data, MOVA_final_predict_file)
#Plot ROC curve
 df$iter3 <- df$iter1 + (df$iter2 * 5)
 out <- cvAUC(df$predict, df$result, label.ordering = NULL, folds = df$iter3)
 plot(out$perf, col = "blue", avg = "vertical", add=TRUE)
 cvauc <- out$cvAUC
 YI2 <- data.frame(matrix(rep(NA, 5), nrow=1))[numeric(0), ]
 colnames(YI2) <- c("Cutoff", "positive_variant_num", "negative_variant_num", "AUC", "cvauc")
 for(i in 6:30){
  pred <- prediction(df[df$iter3==i,]$predict, df[df$iter3==i,]$result)
  auc.tmp <- performance(pred,"auc")
  auc <- as.numeric(auc.tmp@y.values)
  tab <- data.frame(Cutoff=unlist(pred@cutoffs),
   TP=unlist(pred@tp), FP=unlist(pred@fp),
   FN=unlist(pred@fn), TN=unlist(pred@tn),
   Sensitivity=unlist(pred@tp)/(unlist(pred@tp)+unlist(pred@fn)),
   Specificity=unlist(pred@tn)/(unlist(pred@fp)+unlist(pred@tn)),
   Accuracy=((unlist(pred@tp)+unlist(pred@tn))/nrow(df)),
   Precision=(unlist(pred@tp)/(unlist(pred@tp)+unlist(pred@fp))))
  tab$Youden <- tab$Sensitivity+tab$Specificity-1
  YI <- tab[order(tab$Youden, decreasing=T), ]
  YI <- data.frame(Cutoff = YI$Cutoff[1], positive_variant_num = c(nrow(df[(df$result==1)&(df$iter2==1),])), negative_variant_num = c(nrow(df[(df$result==0)&(df$iter2==1),])),AUC = c(auc),cvauc = c(cvauc))
  YI2 <- rbind(YI2, YI)
 }
 fwrite(YI2, paste(protein_name, phenotype, "result_only_pLDDT.csv", sep = "_"))
 df2 <- df %>%
  group_by(ID) %>%
  summarise(only_pLDDT_predict = mean(predict))
 D[ , c('result')] <- list(NULL)
 fwrite(merge(D,df2), MOVA_predict_file)
 fwrite(df, paste(protein_name, phenotype, "predict_orig_only_pLDDT.csv", sep = "_"))
 return(YI2)
}

MOVA_only_BLOSUM62 <- function(protein_name, MOVA_final_predict_file, MOVA_predict_file, phenotype = "Target"){
 D <- fread(MOVA_predict_file)
 D$only_BLOSUM62_predict<- 0
 D[ , c('only_BLOSUM62_predict')] <- list(NULL)
#
 if(phenotype == "Target"){
  D <- D[(D$Type=="Target")|(D$Type=="Ctrl"), ]
 }
#
# Initialize data.frame.
 fpredict <- data.frame(matrix(rep(NA, 3), nrow=1))[numeric(0), ]
 colnames(fpredict) <- c("change", "predict", "n")
# Obtain information on each variant for final predict
 final_data <- fread(MOVA_final_predict_file)
 final_data$only_BLOSUM62_predict<- 0
 final_data[ , c('only_BLOSUM62_predict')] <- list(NULL)
 df <- data.frame(matrix(rep(NA, 6), nrow=1))[numeric(0), ]
 colnames(df) <- c("ID", "result", "predict",  "change", "iter1", "iter2")
 D$id2 <- 0
 D[ , c('id2')] <- list(NULL)
 if(phenotype == "Target"){
  D$result <- D$Type2
 } else {
  D$result <- D$Type3
 }
 oldw <- getOption("warn")
 options(warn=-1)
 D <- rowid_to_column(D, var = "id2")
 for(i2 in 1:5){## split data into train and test using stratified sampling
  D$id3 <- D$id2
  j <- 5
  dat1 <- D %>% stratified(., group = "result", size = 1/j)
  dat1$iter <- 1
  datzan <- D[-dat1$id2, ]
  for(i in 2:j){
   datzan[ , c('id2')] <- list(NULL)
   datzan <- rowid_to_column(datzan, var = "id2")
   if(i != j){
    dat2 <- datzan %>% stratified(., group = "result", size = 1/(j-i+1))
   } else {
    dat2 <- datzan
   }
   dat2$iter <- i
   datzan <- datzan[-dat2$id2, ]
   dat1 <- rbind(dat1, dat2)
  }
  dat1[ , c('id2')] <- list(NULL)
  dat1 <- rowid_to_column(dat1, var = "id2")
  for (i in 1:j){
   val <- dat1[dat1$iter == i,]
   train.data <- dat1[-val$id2,]#Model construction with randomForest.
#Type2 is the objective variable and con+x+y+z+b is the explanatory variable.
   rf<-randomForest(result~con,train.data, importance=TRUE)
#Store the predicted values of the test data in data.frame.
   f <-data.frame(ID = val$ID, result = val$result, predict = predict(rf,val), change = val$change, iter1 = i, iter2 =i2)
   df <- rbind(df, f)
  }
 }
#Store the predicted values of all variants for final predict in data.frame.
 for (i in 1:30){
  rf <-randomForest(result~con,D)
  f <-data.frame(ID = final_data$ID, predict = predict(rf,final_data), n = i) 
  fpredict <- rbind(fpredict, f)
 }
 options(warn = oldw)
#Take the average of the predictions to make the final prediction.
 fpredictb <- fpredict %>%
  group_by(ID) %>%
  summarise(only_BLOSUM62_predict = mean(predict))
#If the amino acid residues do not change, the predicted value is 0.
 ID <- str_split(fpredictb$ID, pattern = ":", simplify = TRUE)
 for (i in 1:nrow(fpredictb)){
  if(ID[i, 5]==ID[i, 7]){
   fpredictb$predict[i] = 0
  }
 }

 final_data <- merge(fpredictb, final_data)
#Save the final prediction.
 fwrite(final_data, MOVA_final_predict_file)
#Plot ROC curve
 df$iter3 <- df$iter1 + (df$iter2 * 5)
 out <- cvAUC(df$predict, df$result, label.ordering = NULL, folds = df$iter3)
 plot(out$perf, col = "yellow", avg = "vertical", add=TRUE)
 cvauc <- out$cvAUC
 YI2 <- data.frame(matrix(rep(NA, 5), nrow=1))[numeric(0), ]
 colnames(YI2) <- c("Cutoff", "positive_variant_num", "negative_variant_num", "AUC", "cvauc")
 for(i in 6:30){
  pred <- prediction(df[df$iter3==i,]$predict, df[df$iter3==i,]$result)
  auc.tmp <- performance(pred,"auc")
  auc <- as.numeric(auc.tmp@y.values)
  tab <- data.frame(Cutoff=unlist(pred@cutoffs),
   TP=unlist(pred@tp), FP=unlist(pred@fp),
   FN=unlist(pred@fn), TN=unlist(pred@tn),
   Sensitivity=unlist(pred@tp)/(unlist(pred@tp)+unlist(pred@fn)),
   Specificity=unlist(pred@tn)/(unlist(pred@fp)+unlist(pred@tn)),
   Accuracy=((unlist(pred@tp)+unlist(pred@tn))/nrow(df)),
   Precision=(unlist(pred@tp)/(unlist(pred@tp)+unlist(pred@fp))))
  tab$Youden <- tab$Sensitivity+tab$Specificity-1
  YI <- tab[order(tab$Youden, decreasing=T), ]
  YI <- data.frame(Cutoff = YI$Cutoff[1], positive_variant_num = c(nrow(df[(df$result==1)&(df$iter2==1),])), negative_variant_num = c(nrow(df[(df$result==0)&(df$iter2==1),])),AUC = c(auc),cvauc = c(cvauc))
  YI2 <- rbind(YI2, YI)
 }
 fwrite(YI2, paste(protein_name, phenotype, "result_only_BLOSUM62.csv", sep = "_"))
 df2 <- df %>%
  group_by(ID) %>%
  summarise(only_BLOSUM62_predict = mean(predict))
 D[ , c('result')] <- list(NULL)
 fwrite(merge(D,df2), MOVA_predict_file)
 fwrite(df, paste(protein_name, phenotype, "predict_orig_only_BLOSUM62.csv", sep = "_"))
 return(YI2)
}

MOVA_only_3d_coordinates <- function(protein_name, MOVA_final_predict_file, MOVA_predict_file, phenotype = "Target"){
 D <- fread(MOVA_predict_file)
 D$only_3d_coordinates_predict<- 0
 D[ , c('only_3d_coordinates_predict')] <- list(NULL)
#
 if(phenotype == "Target"){
  D <- D[(D$Type=="Target")|(D$Type=="Ctrl"), ]
 }
#
# Initialize data.frame.
 fpredict <- data.frame(matrix(rep(NA, 3), nrow=1))[numeric(0), ]
 colnames(fpredict) <- c("change", "predict", "n")
# Obtain information on each variant for final predict
 final_data <- fread(MOVA_final_predict_file)
 final_data$only_3d_coordinates_predict<- 0
 final_data[ , c('only_3d_coordinates_predict')] <- list(NULL)
 df <- data.frame(matrix(rep(NA, 6), nrow=1))[numeric(0), ]
 colnames(df) <- c("ID", "result", "predict",  "change", "iter1", "iter2")
 D$id2 <- 0
 D[ , c('id2')] <- list(NULL)
 if(phenotype == "Target"){
  D$result <- D$Type2
 } else {
  D$result <- D$Type3
 }
 oldw <- getOption("warn")
 options(warn=-1)
 D <- rowid_to_column(D, var = "id2")
 for(i2 in 1:5){## split data into train and test using stratified sampling
  D$id3 <- D$id2
  j <- 5
  dat1 <- D %>% stratified(., group = "result", size = 1/j)
  dat1$iter <- 1
  datzan <- D[-dat1$id2, ]
  for(i in 2:j){
   datzan[ , c('id2')] <- list(NULL)
   datzan <- rowid_to_column(datzan, var = "id2")
   if(i != j){
    dat2 <- datzan %>% stratified(., group = "result", size = 1/(j-i+1))
   } else {
    dat2 <- datzan
   }
   dat2$iter <- i
   datzan <- datzan[-dat2$id2, ]
   dat1 <- rbind(dat1, dat2)
  }
  dat1[ , c('id2')] <- list(NULL)
  dat1 <- rowid_to_column(dat1, var = "id2")
  for (i in 1:j){
   val <- dat1[dat1$iter == i,]
   train.data <- dat1[-val$id2,]#Model construction with randomForest.
#Type2 is the objective variable and con+x+y+z+b is the explanatory variable.
    rf<-randomForest(result~x+y+z,train.data, importance=TRUE)
#Store the predicted values of the test data in data.frame.
    f <-data.frame(ID = val$ID, result = val$result, predict = predict(rf,val), change = val$change, iter1 = i, iter2 =i2)
   df <- rbind(df, f)
  }
 }
#Store the predicted values of all variants for final predict in data.frame.
 for (i in 1:30){
  rf <-randomForest(result~x+y+z,D)
  f <-data.frame(ID = final_data$ID, predict = predict(rf,final_data), n = i) 
  fpredict <- rbind(fpredict, f)
 }
 options(warn = oldw)
#Take the average of the predictions to make the final prediction.
 fpredictb <- fpredict %>%
  group_by(ID) %>%
  summarise(only_3d_coordinates_predict = mean(predict))
#If the amino acid residues do not change, the predicted value is 0.
 ID <- str_split(fpredictb$ID, pattern = ":", simplify = TRUE)
 for (i in 1:nrow(fpredictb)){
  if(ID[i, 5]==ID[i, 7]){
   fpredictb$predict[i] = 0
  }
 }

 final_data <- merge(fpredictb, final_data)
#Save the final prediction.
 fwrite(final_data, MOVA_final_predict_file)
#Plot ROC curve
 df$iter3 <- df$iter1 + (df$iter2 * 5)
 out <- cvAUC(df$predict, df$result, label.ordering = NULL, folds = df$iter3)
 plot(out$perf, col = "green", avg = "vertical", add=TRUE)
 cvauc <- out$cvAUC
 YI2 <- data.frame(matrix(rep(NA, 5), nrow=1))[numeric(0), ]
 colnames(YI2) <- c("Cutoff", "positive_variant_num", "negative_variant_num", "AUC", "cvauc")
 for(i in 6:30){
  pred <- prediction(df[df$iter3==i,]$predict, df[df$iter3==i,]$result)
  auc.tmp <- performance(pred,"auc")
  auc <- as.numeric(auc.tmp@y.values)
  tab <- data.frame(Cutoff=unlist(pred@cutoffs),
   TP=unlist(pred@tp), FP=unlist(pred@fp),
   FN=unlist(pred@fn), TN=unlist(pred@tn),
   Sensitivity=unlist(pred@tp)/(unlist(pred@tp)+unlist(pred@fn)),
   Specificity=unlist(pred@tn)/(unlist(pred@fp)+unlist(pred@tn)),
   Accuracy=((unlist(pred@tp)+unlist(pred@tn))/nrow(df)),
   Precision=(unlist(pred@tp)/(unlist(pred@tp)+unlist(pred@fp))))
  tab$Youden <- tab$Sensitivity+tab$Specificity-1
  YI <- tab[order(tab$Youden, decreasing=T), ]
  YI <- data.frame(Cutoff = YI$Cutoff[1], positive_variant_num = c(nrow(df[(df$result==1)&(df$iter2==1),])), negative_variant_num = c(nrow(df[(df$result==0)&(df$iter2==1),])),AUC = c(auc),cvauc = c(cvauc))
  YI2 <- rbind(YI2, YI)
 }
 fwrite(YI2, paste(protein_name, phenotype, "result_only_3d_coordinates.csv", sep = "_"))
 df2 <- df %>%
  group_by(ID) %>%
  summarise(only_3d_coordinates_predict = mean(predict))
 D[ , c('result')] <- list(NULL)
 fwrite(merge(D,df2), MOVA_predict_file)
 fwrite(df, paste(protein_name, phenotype, "predict_orig_only_3d_coordinates.csv", sep = "_"))
 return(YI2)
}


MOVA_only_3d_distance <- function(protein_name, MOVA_final_predict_file, MOVA_predict_file, phenotype = "Target"){
 D <- fread(MOVA_predict_file)
 D$t_distance<- 0
 D[ , c('t_distance')] <- list(NULL)
#
 if(phenotype == "Target"){
  D <- D[(D$Type=="Target")|(D$Type=="Ctrl"), ]
 }
#
# Initialize data.frame.
 fpredict <- data.frame(matrix(rep(NA, 2), nrow=1))[numeric(0), ]
 colnames(fpredict) <- c("ID", "distance")
# Obtain information on each variant for final predict
 final_data <- fread(MOVA_final_predict_file)
 df <- data.frame(matrix(rep(NA, 6), nrow=1))[numeric(0), ]
 colnames(df) <- c("ID", "re_result", "distance",  "change", "iter1", "iter2")
 D$id2 <- 0
 D[ , c('id2')] <- list(NULL)
 D$re_result <- 0
 D[ , c('re_result')] <- list(NULL)
 D$t_distance <- 0
 D[ , c('t_distance')] <- list(NULL)
 D$iter1 <- 0
 D[ , c('iter1')] <- list(NULL)
 D$iter2 <- 0
 D[ , c('iter2')] <- list(NULL)
 if(phenotype == "Target"){
  D$result <- D$Type2
 } else {
  D$result <- D$Type3
 }
 oldw <- getOption("warn")
 options(warn=-1)
 D <- rowid_to_column(D, var = "id2")
 for(i2 in 1:5){## split data into train and test using stratified sampling
  D$id3 <- D$id2
  j <- 5
  dat1 <- D %>% stratified(., group = "result", size = 1/j)
  dat1$iter <- 1
  datzan <- D[-dat1$id2, ]
  for(i in 2:j){
   datzan[ , c('id2')] <- list(NULL)
   datzan <- rowid_to_column(datzan, var = "id2")
   if(i != j){
    dat2 <- datzan %>% stratified(., group = "result", size = 1/(j-i+1))
   } else {
    dat2 <- datzan
   }
   dat2$iter <- i
   datzan <- datzan[-dat2$id2, ]
   dat1 <- rbind(dat1, dat2)
  }
  dat1[ , c('id2')] <- list(NULL)
  dat1 <- rowid_to_column(dat1, var = "id2")
  for (i in 1:j){
   val <- dat1[dat1$iter == i,]
   train.data <- dat1[-val$id2,]#Model construction with randomForest.
#Type2 is the objective variable and con+x+y+z+b is the explanatory variable.
    x <- train.data[train.data$result==1,]$x
    y <- train.data[train.data$result==1,]$y
    z <- train.data[train.data$result==1,]$z
    for(i3 in 1:nrow(val)){
      f <-data.frame(ID = val[i3]$ID, re_result = abs(val[i3]$result-1), t_distance = min(sqrt((x-val[i3]$x)*(x-val[i3]$x)+(y-val[i3]$y)*(y-val[i3]$y)+(z-val[i3]$z)*(z-val[i3]$z))), change = val[i3]$change, iter1 = i, iter2 =i2)
     df <- rbind(df, f)
    }
  }
 }

 options(warn = oldw)
#Take the average of the predictions to make the final prediction.
#If the amino acid residues do not change, the predicted value is 0.

#Save the final prediction.
#Plot ROC curve
 df$iter3 <- df$iter1 + (df$iter2 * 5)
 out <- cvAUC(df$t_distance, df$re_result, label.ordering = NULL, folds = df$iter3)
 plot(out$perf, col = "black", avg = "vertical", add=TRUE)
 cvauc <- out$cvAUC
 YI2 <- data.frame(matrix(rep(NA, 5), nrow=1))[numeric(0), ]
 colnames(YI2) <- c("Cutoff", "positive_variant_num", "negative_variant_num", "AUC", "cvauc")
 for(i in 6:30){
  pred <- prediction(df[df$iter3==i,]$t_distance, df[df$iter3==i,]$re_result)
  auc.tmp <- performance(pred,"auc")
  auc <- as.numeric(auc.tmp@y.values)
  tab <- data.frame(Cutoff=unlist(pred@cutoffs),
   TP=unlist(pred@tp), FP=unlist(pred@fp),
   FN=unlist(pred@fn), TN=unlist(pred@tn),
   Sensitivity=unlist(pred@tp)/(unlist(pred@tp)+unlist(pred@fn)),
   Specificity=unlist(pred@tn)/(unlist(pred@fp)+unlist(pred@tn)),
   Accuracy=((unlist(pred@tp)+unlist(pred@tn))/nrow(df)),
   Precision=(unlist(pred@tp)/(unlist(pred@tp)+unlist(pred@fp))))
  tab$Youden <- tab$Sensitivity+tab$Specificity-1
  YI <- tab[order(tab$Youden, decreasing=T), ]
  YI <- data.frame(Cutoff = YI$Cutoff[1], positive_variant_num = c(nrow(df[(df$re_result==0)&(df$iter2==1),])), negative_variant_num = c(nrow(df[(df$re_result==1)&(df$iter2==1),])),AUC = c(auc),cvauc = c(cvauc))
  YI2 <- rbind(YI2, YI)
 }
 fwrite(YI2, paste(protein_name, phenotype, "result_only_3d_distance.csv", sep = "_"))
 df2 <- df %>%
  group_by(ID) %>%
  summarise(t_distance = mean(t_distance))
 D[ , c('result')] <- list(NULL)
 fwrite(merge(D,df2), MOVA_predict_file)
 fwrite(df, paste(protein_name, phenotype, "predict_orig_only_3d_distance.csv", sep = "_"))
 return(YI2)
}


MOVA_predistance <- function(protein_name, MOVA_final_predict_file, MOVA_predict_file, phenotype = "Target"){
 D <- fread(MOVA_predict_file)
#
 if(phenotype == "Target"){
  D <- D[(D$Type=="Target")|(D$Type=="Ctrl"), ]
 }
 final_data <- fread(MOVA_final_predict_file)
  if(phenotype == "Target"){
   x <- D[Type2==1,]$x
   y <- D[Type2==1,]$y
   z <- D[Type2==1,]$z
   D$t_distance <- 0
    final_data$t_distance <- 0
   for(i3 in 1:nrow(final_data)){
     t_distance <- sqrt((x-final_data[i3]$x)*(x-final_data[i3]$x)+(y-final_data[i3]$y)*(y-final_data[i3]$y)+(z-final_data[i3]$z)*(z-final_data[i3]$z))
      t_distance <- sort(t_distance,decreasing=F)
       final_data[i3]$t_distance <- t_distance[1]

    }
  } else {
   x <- D[Type3==1,]$x
   y <- D[Type3==1,]$y
   z <- D[Type3==1,]$z
    D$t_distance <- 0
    final_data$t_distance <- 0
   for(i3 in 1:nrow(final_data)){
     t_distance <- sqrt((x-final_data[i3]$x)*(x-final_data[i3]$x)+(y-final_data[i3]$y)*(y-final_data[i3]$y)+(z-final_data[i3]$z)*(z-final_data[i3]$z))
      t_distance <- sort(t_distance,decreasing=F)
       final_data[i3]$t_distance <- t_distance[1]
    }
  }
 fwrite(final_data, MOVA_final_predict_file)
}


MOVA_3d_distance <- function(protein_name, MOVA_final_predict_file, MOVA_predict_file, phenotype = "Target"){
 D <- fread(MOVA_predict_file)
 D$MOVA_3d_distance_predict<- 0
 D[ , c('MOVA_3d_distance_predict')] <- list(NULL)
#
 if(phenotype == "Target"){
  D <- D[(D$Type=="Target")|(D$Type=="Ctrl"), ]
 }
#
# Initialize data.frame.
 fpredict <- data.frame(matrix(rep(NA, 3), nrow=1))[numeric(0), ]
 colnames(fpredict) <- c("change", "predict", "n")
 imp.rfa <- data.frame(matrix(rep(NA, 3), nrow=1))[numeric(0), ]
 colnames(imp.rfa) <- c("name", "IncNodePurity", "n")
# Obtain information on each variant for final predict
 final_data <- fread(MOVA_final_predict_file)
 final_data$MOVA_3d_distance_predict<- 0
 final_data[ , c('MOVA_3d_distance_predict')] <- list(NULL)
 df <- data.frame(matrix(rep(NA, 6), nrow=1))[numeric(0), ]
 colnames(df) <- c("ID", "re_result", "distance",  "change", "iter1", "iter2")
 D$id2 <- 0
 D[ , c('id2')] <- list(NULL)
 if(phenotype == "Target"){
  D$result <- D$Type2
 } else {
  D$result <- D$Type3
 }
 oldw <- getOption("warn")
 options(warn=-1)
 D <- rowid_to_column(D, var = "id2")
 for(i2 in 1:5){## split data into train and test using stratified sampling
  D$id3 <- D$id2
  j <- 5
  dat1 <- D %>% stratified(., group = "result", size = 1/j)
  dat1$iter <- 1
  datzan <- D[-dat1$id2, ]
  for(i in 2:j){
   datzan[ , c('id2')] <- list(NULL)
   datzan <- rowid_to_column(datzan, var = "id2")
   if(i != j){
    dat2 <- datzan %>% stratified(., group = "result", size = 1/(j-i+1))
   } else {
    dat2 <- datzan
   }
   dat2$iter <- i
   datzan <- datzan[-dat2$id2, ]
   dat1 <- rbind(dat1, dat2)
  }
  dat1[ , c('id2')] <- list(NULL)
  dat1 <- rowid_to_column(dat1, var = "id2")
  for (i in 1:j){
   val <- dat1[dat1$iter == i,]
   train.data <- dat1[-val$id2,]
#Model construction with randomForest.
#Type2 is the objective variable and con+x+y+z+b is the explanatory variable.
    x <- train.data[train.data$result==1,]$x
    y <- train.data[train.data$result==1,]$y
    z <- train.data[train.data$result==1,]$z
    train.data$t_distance <- 0
    for(i3 in 1:nrow(train.data)){
      t_distance <- sqrt((x-train.data[i3]$x)*(x-train.data[i3]$x)+(y-train.data[i3]$y)*(y-train.data[i3]$y)+(z-train.data[i3]$z)*(z-train.data[i3]$z))
      t_distance <- sort(t_distance,decreasing=F)
      if(train.data[i3]$result == 1){
        train.data[i3]$t_distance <- t_distance[2]
      } else {
        train.data[i3]$t_distance <- t_distance[1]
      }
    }
    rf<-randomForest(result~con+t_distance+b,train.data, importance=TRUE)
#Measure variable importance.
   imp.rf <- importance(rf)
#Store the predicted values of the test data in data.frame.
    x <- train.data[train.data$result==1,]$x
    y <- train.data[train.data$result==1,]$y
    z <- train.data[train.data$result==1,]$z
    val$t_distance <- 0
    for(i3 in 1:nrow(val)){
      val[i3]$t_distance <- min(sqrt((x-val[i3]$x)*(x-val[i3]$x)+(y-val[i3]$y)*(y-val[i3]$y)+(z-val[i3]$z)*(z-val[i3]$z)))
    }
    f <-data.frame(ID = val$ID, result = val$result, predict = predict(rf,val), change = val$change, iter1 = i, iter2 =i2)
   df <- rbind(df, f)
#Store the variable importance in data.frame.
  }
 }

#Store the variable importance in data.frame.
   f <- data.frame(name = rownames(imp.rf), IncNodePurity = imp.rf, n = i)
   imp.rfa <- rbind(imp.rfa,f)
#Store the predicted values of all variants for final predict in data.frame.
   x <- D[D$result==1,]$x
   y <- D[D$result==1,]$y
   z <- D[D$result==1,]$z
   D$t_distance <- 0
   for(i3 in 1:nrow(D)){
     t_distance <- sqrt((x-D[i3]$x)*(x-D[i3]$x)+(y-D[i3]$y)*(y-D[i3]$y)+(z-D[i3]$z)*(z-D[i3]$z))
      t_distance <- sort(t_distance,decreasing=F)
     if(D[i3]$result == 1){
       D[i3]$t_distance <- t_distance[2]
     } else {
       D[i3]$t_distance <- t_distance[1]
     }

   }
 for (i in 1:30){
   rf <-randomForest(result~con+t_distance+b,D)
  f <-data.frame(ID = final_data$ID, predict = predict(rf,final_data), n = i) 
  fpredict <- rbind(fpredict, f)
 }

 options(warn = oldw)
 fpredictb <- fpredict %>%
  group_by(ID) %>%
  summarise(MOVA_3d_distance_predict = mean(predict))
#Take the average of variable importance
 imp.rfb <- imp.rfa %>%
  group_by(name) %>%
  summarise(IncMSE = mean(IncNodePurity..IncMSE), IncNodePurity = mean(IncNodePurity.IncNodePurity))


 final_data <- merge(fpredictb, final_data)
#Save the final prediction.
 fwrite(final_data, MOVA_final_predict_file)
#Save the average of variable importance.
 fwrite(imp.rfb, paste(protein_name, phenotype, "importance_3d_distance.csv", sep = "_"))
#Plot ROC curve
 df$iter3 <- df$iter1 + (df$iter2 * 5)
 out <- cvAUC(df$predict, df$result, label.ordering = NULL, folds = df$iter3)
 plot(out$perf, col = "blue", avg = "vertical", add=TRUE)
 cvauc <- out$cvAUC
 YI2 <- data.frame(matrix(rep(NA, 5), nrow=1))[numeric(0), ]
 colnames(YI2) <- c("Cutoff", "positive_variant_num", "negative_variant_num", "AUC", "cvauc")
 for(i in 6:30){
  pred <- prediction(df[df$iter3==i,]$predict, df[df$iter3==i,]$result)
  auc.tmp <- performance(pred,"auc")
  auc <- as.numeric(auc.tmp@y.values)
  tab <- data.frame(Cutoff=unlist(pred@cutoffs),
   TP=unlist(pred@tp), FP=unlist(pred@fp),
   FN=unlist(pred@fn), TN=unlist(pred@tn),
   Sensitivity=unlist(pred@tp)/(unlist(pred@tp)+unlist(pred@fn)),
   Specificity=unlist(pred@tn)/(unlist(pred@fp)+unlist(pred@tn)),
   Accuracy=((unlist(pred@tp)+unlist(pred@tn))/nrow(df)),
   Precision=(unlist(pred@tp)/(unlist(pred@tp)+unlist(pred@fp))))
  tab$Youden <- tab$Sensitivity+tab$Specificity-1
  YI <- tab[order(tab$Youden, decreasing=T), ]
  YI <- data.frame(Cutoff = YI$Cutoff[1], positive_variant_num = c(nrow(df[(df$result==1)&(df$iter2==1),])), negative_variant_num = c(nrow(df[(df$result==0)&(df$iter2==1),])),AUC = c(auc),cvauc = c(cvauc))
  YI2 <- rbind(YI2, YI)
 }
 fwrite(YI2, paste(protein_name, phenotype, "result_3d_distance.csv", sep = "_"))
 df2 <- df %>%
  group_by(ID) %>%
  summarise(MOVA_3d_distance_predict = mean(predict))
 D[ , c('result')] <- list(NULL)
 fwrite(merge(D,df2),  MOVA_predict_file)
 fwrite(df, paste(protein_name, phenotype, "predict_orig_3d_distance.csv", sep = "_"))
 return(YI2)
}


MOVA_3d_distance_log <- function(protein_name, MOVA_final_predict_file, MOVA_predict_file, phenotype = "Target"){
 D <- fread(MOVA_predict_file)
 D$MOVA_3d_distance_log_predict<- 0
 D[ , c('MOVA_3d_distance_log_predict')] <- list(NULL)
#
 if(phenotype == "Target"){
  D <- D[(D$Type=="Target")|(D$Type=="Ctrl"), ]
 }
#
# Initialize data.frame.
 fpredict <- data.frame(matrix(rep(NA, 3), nrow=1))[numeric(0), ]
 colnames(fpredict) <- c("change", "predict", "n")
 imp.rfa <- data.frame(matrix(rep(NA, 3), nrow=1))[numeric(0), ]
 colnames(imp.rfa) <- c("name", "IncNodePurity", "n")
# Obtain information on each variant for final predict
 final_data <- fread(MOVA_final_predict_file)
 final_data$MOVA_3d_distance_log_predict<- 0
 final_data[ , c('MOVA_3d_distance_log_predict')] <- list(NULL)
 df <- data.frame(matrix(rep(NA, 6), nrow=1))[numeric(0), ]
 colnames(df) <- c("ID", "re_result", "distance",  "change", "iter1", "iter2")
 D$id2 <- 0
 D[ , c('id2')] <- list(NULL)
 if(phenotype == "Target"){
  D$result <- D$Type2
 } else {
  D$result <- D$Type3
 }
 oldw <- getOption("warn")
 options(warn=-1)
 D <- rowid_to_column(D, var = "id2")
 for(i2 in 1:5){## split data into train and test using stratified sampling
  D$id3 <- D$id2
  j <- 5
  dat1 <- D %>% stratified(., group = "result", size = 1/j)
  dat1$iter <- 1
  datzan <- D[-dat1$id2, ]
  for(i in 2:j){
   datzan[ , c('id2')] <- list(NULL)
   datzan <- rowid_to_column(datzan, var = "id2")
   if(i != j){
    dat2 <- datzan %>% stratified(., group = "result", size = 1/(j-i+1))
   } else {
    dat2 <- datzan
   }
   dat2$iter <- i
   datzan <- datzan[-dat2$id2, ]
   dat1 <- rbind(dat1, dat2)
  }
  dat1[ , c('id2')] <- list(NULL)
  dat1 <- rowid_to_column(dat1, var = "id2")
  for (i in 1:j){
   val <- dat1[dat1$iter == i,]
   train.data <- dat1[-val$id2,]
#Model construction with randomForest.
#Type2 is the objective variable and con+x+y+z+b is the explanatory variable.
    x <- train.data[train.data$result==1,]$x
    y <- train.data[train.data$result==1,]$y
    z <- train.data[train.data$result==1,]$z
    train.data$t_distance <- 0
    for(i3 in 1:nrow(train.data)){
      t_distance <- sqrt((x-train.data[i3]$x)*(x-train.data[i3]$x)+(y-train.data[i3]$y)*(y-train.data[i3]$y)+(z-train.data[i3]$z)*(z-train.data[i3]$z))
      t_distance <- sort(t_distance,decreasing=F)
      if(train.data[i3]$result == 1){
        train.data[i3]$t_distance <- t_distance[2]
      } else {
        train.data[i3]$t_distance <- t_distance[1]
      }
    }
    glm<-glm(result~con+t_distance+b,train.data, family = "binomial")
#Measure variable importance.
#Store the predicted values of the test data in data.frame.
    x <- train.data[train.data$result==1,]$x
    y <- train.data[train.data$result==1,]$y
    z <- train.data[train.data$result==1,]$z
    val$t_distance <- 0
    for(i3 in 1:nrow(val)){
      val[i3]$t_distance <- min(sqrt((x-val[i3]$x)*(x-val[i3]$x)+(y-val[i3]$y)*(y-val[i3]$y)+(z-val[i3]$z)*(z-val[i3]$z)))
    }
    f <-data.frame(ID = val$ID, result = val$result, predict = predict(glm,newdata=val, type="response"), change = val$change, iter1 = i, iter2 =i2)
   df <- rbind(df, f)
  }
 }


#Store the predicted values of all variants for final predict in data.frame.
   x <- D[D$result==1,]$x
   y <- D[D$result==1,]$y
   z <- D[D$result==1,]$z
   D$t_distance <- 0
   for(i3 in 1:nrow(D)){
     t_distance <- sqrt((x-D[i3]$x)*(x-D[i3]$x)+(y-D[i3]$y)*(y-D[i3]$y)+(z-D[i3]$z)*(z-D[i3]$z))
      t_distance <- sort(t_distance,decreasing=F)
     if(D[i3]$result == 1){
       D[i3]$t_distance <- t_distance[2]
     } else {
       D[i3]$t_distance <- t_distance[1]
     }

   }
 for (i in 1:30){
   glm<-glm(result~con+t_distance+b,D, family = "binomial")
  f <-data.frame(ID = final_data$ID, predict = predict(glm,newdata=final_data, type="response"), n = i) 
  fpredict <- rbind(fpredict, f)
 }

 options(warn = oldw)
 fpredictb <- fpredict %>%
  group_by(ID) %>%
  summarise(MOVA_3d_distance_log_predict = mean(predict))

 final_data <- merge(fpredictb, final_data)
#Save the final prediction.
 fwrite(final_data, MOVA_final_predict_file)
#Plot ROC curve
 df$iter3 <- df$iter1 + (df$iter2 * 5)
 out <- cvAUC(df$predict, df$result, label.ordering = NULL, folds = df$iter3)
 plot(out$perf, col = "blue", avg = "vertical", add=TRUE)
 cvauc <- out$cvAUC
 YI2 <- data.frame(matrix(rep(NA, 5), nrow=1))[numeric(0), ]
 colnames(YI2) <- c("Cutoff", "positive_variant_num", "negative_variant_num", "AUC", "cvauc")
 for(i in 6:30){
  pred <- prediction(df[df$iter3==i,]$predict, df[df$iter3==i,]$result)
  auc.tmp <- performance(pred,"auc")
  auc <- as.numeric(auc.tmp@y.values)
  tab <- data.frame(Cutoff=unlist(pred@cutoffs),
   TP=unlist(pred@tp), FP=unlist(pred@fp),
   FN=unlist(pred@fn), TN=unlist(pred@tn),
   Sensitivity=unlist(pred@tp)/(unlist(pred@tp)+unlist(pred@fn)),
   Specificity=unlist(pred@tn)/(unlist(pred@fp)+unlist(pred@tn)),
   Accuracy=((unlist(pred@tp)+unlist(pred@tn))/nrow(df)),
   Precision=(unlist(pred@tp)/(unlist(pred@tp)+unlist(pred@fp))))
  tab$Youden <- tab$Sensitivity+tab$Specificity-1
  YI <- tab[order(tab$Youden, decreasing=T), ]
  YI <- data.frame(Cutoff = YI$Cutoff[1], positive_variant_num = c(nrow(df[(df$result==1)&(df$iter2==1),])), negative_variant_num = c(nrow(df[(df$result==0)&(df$iter2==1),])),AUC = c(auc),cvauc = c(cvauc))
  YI2 <- rbind(YI2, YI)
 }
 fwrite(YI2, paste(protein_name, phenotype, "result_3d_distance_log.csv", sep = "_"))
 df2 <- df %>%
  group_by(ID) %>%
  summarise(MOVA_3d_distance_log_predict = mean(predict))
 D[ , c('result')] <- list(NULL)
 fwrite(merge(D,df2),  MOVA_predict_file)
 fwrite(df, paste(protein_name, phenotype, "predict_orig_3d_distance_log.csv", sep = "_"))
 return(YI2)
}



MOVA_plus_3d_distance <- function(protein_name, MOVA_final_predict_file, MOVA_predict_file, phenotype = "Target"){
 D <- fread(MOVA_predict_file)
 D$MOVA_plus_3d_distance_predict<- 0
 D[ , c('MOVA_plus_3d_distance_predict')] <- list(NULL)
#
 if(phenotype == "Target"){
  D <- D[(D$Type=="Target")|(D$Type=="Ctrl"), ]
 }
#
# Initialize data.frame.
 fpredict <- data.frame(matrix(rep(NA, 3), nrow=1))[numeric(0), ]
 colnames(fpredict) <- c("change", "predict", "n")
 imp.rfa <- data.frame(matrix(rep(NA, 3), nrow=1))[numeric(0), ]
 colnames(imp.rfa) <- c("name", "IncNodePurity", "n")
# Obtain information on each variant for final predict
 final_data <- fread(MOVA_final_predict_file)
 final_data$MOVA_plus_3d_distance_predict<- 0
 final_data[ , c('MOVA_plus_3d_distance_predict')] <- list(NULL)
 df <- data.frame(matrix(rep(NA, 6), nrow=1))[numeric(0), ]
 colnames(df) <- c("ID", "re_result", "distance",  "change", "iter1", "iter2")
 D$id2 <- 0
 D[ , c('id2')] <- list(NULL)
 if(phenotype == "Target"){
  D$result <- D$Type2
 } else {
  D$result <- D$Type3
 }
 oldw <- getOption("warn")
 options(warn=-1)
 D <- rowid_to_column(D, var = "id2")
 for(i2 in 1:5){## split data into train and test using stratified sampling
  D$id3 <- D$id2
  j <- 5
  dat1 <- D %>% stratified(., group = "result", size = 1/j)
  dat1$iter <- 1
  datzan <- D[-dat1$id2, ]
  for(i in 2:j){
   datzan[ , c('id2')] <- list(NULL)
   datzan <- rowid_to_column(datzan, var = "id2")
   if(i != j){
    dat2 <- datzan %>% stratified(., group = "result", size = 1/(j-i+1))
   } else {
    dat2 <- datzan
   }
   dat2$iter <- i
   datzan <- datzan[-dat2$id2, ]
   dat1 <- rbind(dat1, dat2)
  }
  dat1[ , c('id2')] <- list(NULL)
  dat1 <- rowid_to_column(dat1, var = "id2")
  for (i in 1:j){
   val <- dat1[dat1$iter == i,]
   train.data <- dat1[-val$id2,]
#Model construction with randomForest.
#Type2 is the objective variable and con+x+y+z+b is the explanatory variable.
    x <- train.data[train.data$result==1,]$x
    y <- train.data[train.data$result==1,]$y
    z <- train.data[train.data$result==1,]$z
    train.data$t_distance <- 0
    for(i3 in 1:nrow(train.data)){
      t_distance <- sqrt((x-train.data[i3]$x)*(x-train.data[i3]$x)+(y-train.data[i3]$y)*(y-train.data[i3]$y)+(z-train.data[i3]$z)*(z-train.data[i3]$z))
      t_distance <- sort(t_distance,decreasing=F)
      if(train.data[i3]$result == 1){
        train.data[i3]$t_distance <- t_distance[2]
      } else {
        train.data[i3]$t_distance <- t_distance[1]
      }
    }
    rf<-randomForest(result~con+x+y+z+t_distance+b,train.data, importance=TRUE)
#Measure variable importance.
   imp.rf <- importance(rf)
#Store the predicted values of the test data in data.frame.
    x <- train.data[train.data$result==1,]$x
    y <- train.data[train.data$result==1,]$y
    z <- train.data[train.data$result==1,]$z
    val$t_distance <- 0
    for(i3 in 1:nrow(val)){
      val[i3]$t_distance <- min(sqrt((x-val[i3]$x)*(x-val[i3]$x)+(y-val[i3]$y)*(y-val[i3]$y)+(z-val[i3]$z)*(z-val[i3]$z)))
    }
    f <-data.frame(ID = val$ID, result = val$result, predict = predict(rf,val), change = val$change, iter1 = i, iter2 =i2)
   df <- rbind(df, f)
#Store the variable importance in data.frame.
  }
 }

#Store the variable importance in data.frame.
   f <- data.frame(name = rownames(imp.rf), IncNodePurity = imp.rf, n = i)
   imp.rfa <- rbind(imp.rfa,f)
#Store the predicted values of all variants for final predict in data.frame.
   x <- D[D$result==1,]$x
   y <- D[D$result==1,]$y
   z <- D[D$result==1,]$z
   D$t_distance <- 0
   for(i3 in 1:nrow(D)){
     t_distance <- sqrt((x-D[i3]$x)*(x-D[i3]$x)+(y-D[i3]$y)*(y-D[i3]$y)+(z-D[i3]$z)*(z-D[i3]$z))
      t_distance <- sort(t_distance,decreasing=F)
     if(D[i3]$result == 1){
       D[i3]$t_distance <- t_distance[2]
     } else {
       D[i3]$t_distance <- t_distance[1]
     }

   }
 for (i in 1:30){
   rf <-randomForest(result~con+x+y+z+t_distance+b,D)
  f <-data.frame(ID = final_data$ID, predict = predict(rf,final_data), n = i) 
  fpredict <- rbind(fpredict, f)
 }

 options(warn = oldw)
 fpredictb <- fpredict %>%
  group_by(ID) %>%
  summarise(MOVA_plus_3d_distance_predict = mean(predict))
#Take the average of variable importance
 imp.rfb <- imp.rfa %>%
  group_by(name) %>%
  summarise(IncMSE = mean(IncNodePurity..IncMSE), IncNodePurity = mean(IncNodePurity.IncNodePurity))


 final_data <- merge(fpredictb, final_data)
#Save the final prediction.
 fwrite(final_data, MOVA_final_predict_file)
#Save the average of variable importance.
 fwrite(imp.rfb, paste(protein_name, phenotype, "importance_plus_3d_distance.csv", sep = "_"))
#Plot ROC curve
 df$iter3 <- df$iter1 + (df$iter2 * 5)
 out <- cvAUC(df$predict, df$result, label.ordering = NULL, folds = df$iter3)
 plot(out$perf, col = "blue", avg = "vertical", add=TRUE)
 cvauc <- out$cvAUC
 YI2 <- data.frame(matrix(rep(NA, 5), nrow=1))[numeric(0), ]
 colnames(YI2) <- c("Cutoff", "positive_variant_num", "negative_variant_num", "AUC", "cvauc")
 for(i in 6:30){
  pred <- prediction(df[df$iter3==i,]$predict, df[df$iter3==i,]$result)
  auc.tmp <- performance(pred,"auc")
  auc <- as.numeric(auc.tmp@y.values)
  tab <- data.frame(Cutoff=unlist(pred@cutoffs),
   TP=unlist(pred@tp), FP=unlist(pred@fp),
   FN=unlist(pred@fn), TN=unlist(pred@tn),
   Sensitivity=unlist(pred@tp)/(unlist(pred@tp)+unlist(pred@fn)),
   Specificity=unlist(pred@tn)/(unlist(pred@fp)+unlist(pred@tn)),
   Accuracy=((unlist(pred@tp)+unlist(pred@tn))/nrow(df)),
   Precision=(unlist(pred@tp)/(unlist(pred@tp)+unlist(pred@fp))))
  tab$Youden <- tab$Sensitivity+tab$Specificity-1
  YI <- tab[order(tab$Youden, decreasing=T), ]
  YI <- data.frame(Cutoff = YI$Cutoff[1], positive_variant_num = c(nrow(df[(df$result==1)&(df$iter2==1),])), negative_variant_num = c(nrow(df[(df$result==0)&(df$iter2==1),])),AUC = c(auc),cvauc = c(cvauc))
  YI2 <- rbind(YI2, YI)
 }
 fwrite(YI2, paste(protein_name, phenotype, "result_plus_3d_distance.csv", sep = "_"))
 df2 <- df %>%
  group_by(ID) %>%
  summarise(MOVA_plus_3d_distance_predict = mean(predict))
 D[ , c('result')] <- list(NULL)
 fwrite(merge(D,df2),  MOVA_predict_file)
 fwrite(df, paste(protein_name, phenotype, "predict_orig_plus_3d_distance.csv", sep = "_"))
 return(YI2)
}


MOVA_plus_3d_distance_SVM <- function(protein_name, MOVA_final_predict_file, MOVA_predict_file, phenotype = "Target"){
 D <- fread(MOVA_predict_file)
 D$MOVA_plus_3d_distance_SVM_predict<- 0
 D[ , c('MOVA_plus_3d_distance_SVM_predict')] <- list(NULL)
#
 if(phenotype == "Target"){
  D <- D[(D$Type=="Target")|(D$Type=="Ctrl"), ]
 }
#
# Initialize data.frame.
 fpredict <- data.frame(matrix(rep(NA, 3), nrow=1))[numeric(0), ]
 colnames(fpredict) <- c("change", "predict", "n")
 imp.rfa <- data.frame(matrix(rep(NA, 3), nrow=1))[numeric(0), ]
 colnames(imp.rfa) <- c("name", "IncNodePurity", "n")
# Obtain information on each variant for final predict
 final_data <- fread(MOVA_final_predict_file)
 final_data$MOVA_plus_3d_distance_SVM_predict<- 0
 final_data[ , c('MOVA_plus_3d_distance_SVM_predict')] <- list(NULL)
 df <- data.frame(matrix(rep(NA, 6), nrow=1))[numeric(0), ]
 colnames(df) <- c("ID", "re_result", "distance",  "change", "iter1", "iter2")
 D$id2 <- 0
 if(phenotype == "Target"){
  D$result <- D$Type2
 } else {
  D$result <- D$Type3
 }
 D[ , c('id2')] <- list(NULL)
 oldw <- getOption("warn")
 options(warn=-1)
 D <- rowid_to_column(D, var = "id2")
 for(i2 in 1:5){## split data into train and test using stratified sampling
  D$id3 <- D$id2
  j <- 5
  dat1 <- D %>% stratified(., group = "result", size = 1/j)
  dat1$iter <- 1
  datzan <- D[-dat1$id2, ]
  for(i in 2:j){
   datzan[ , c('id2')] <- list(NULL)
   datzan <- rowid_to_column(datzan, var = "id2")
   if(i != j){
    dat2 <- datzan %>% stratified(., group = "result", size = 1/(j-i+1))
   } else {
    dat2 <- datzan
   }
   dat2$iter <- i
   datzan <- datzan[-dat2$id2, ]
   dat1 <- rbind(dat1, dat2)
  }
  dat1[ , c('id2')] <- list(NULL)
  dat1 <- rowid_to_column(dat1, var = "id2")
  for (i in 1:j){
   val <- dat1[dat1$iter == i,]
   train.data <- dat1[-val$id2,]
#Model construction with randomForest.
#Type2 is the objective variable and con+x+y+z+b is the explanatory variable.
    x <- train.data[train.data$result==1,]$x
    y <- train.data[train.data$result==1,]$y
    z <- train.data[train.data$result==1,]$z
    train.data$t_distance <- 0
    for(i3 in 1:nrow(train.data)){
      t_distance <- sqrt((x-train.data[i3]$x)*(x-train.data[i3]$x)+(y-train.data[i3]$y)*(y-train.data[i3]$y)+(z-train.data[i3]$z)*(z-train.data[i3]$z))
      t_distance <- sort(t_distance,decreasing=F)
      if(train.data[i3]$result == 1){
        train.data[i3]$t_distance <- t_distance[2]
      } else {
        train.data[i3]$t_distance <- t_distance[1]
      }
    }
    svm<-ksvm(result~con+x+y+z+t_distance+b,train.data, importance=TRUE)
#Measure variable importance.
#Store the predicted values of the test data in data.frame.
    x <- train.data[train.data$result==1,]$x
    y <- train.data[train.data$result==1,]$y
    z <- train.data[train.data$result==1,]$z
    val$t_distance <- 0
    for(i3 in 1:nrow(val)){
      val[i3]$t_distance <- min(sqrt((x-val[i3]$x)*(x-val[i3]$x)+(y-val[i3]$y)*(y-val[i3]$y)+(z-val[i3]$z)*(z-val[i3]$z)))
    }
    f <-data.frame(ID = val$ID, result = val$result, predict = predict(svm,val), change = val$change, iter1 = i, iter2 =i2)
   df <- rbind(df, f)
  }
 }

#Store the predicted values of all variants for final predict in data.frame.
   x <- D[D$result==1,]$x
   y <- D[D$result==1,]$y
   z <- D[D$result==1,]$z
   D$t_distance <- 0
   for(i3 in 1:nrow(D)){
     t_distance <- sqrt((x-D[i3]$x)*(x-D[i3]$x)+(y-D[i3]$y)*(y-D[i3]$y)+(z-D[i3]$z)*(z-D[i3]$z))
      t_distance <- sort(t_distance,decreasing=F)
     if(D[i3]$result == 1){
       D[i3]$t_distance <- t_distance[2]
     } else {
       D[i3]$t_distance <- t_distance[1]
     }

   }
 for (i in 1:30){
   svm<-ksvm(result~con+x+y+z+t_distance+b,D)
  f <-data.frame(ID = final_data$ID, predict = predict(svm,final_data), n = i) 
  fpredict <- rbind(fpredict, f)
 }

 options(warn = oldw)
 fpredictb <- fpredict %>%
  group_by(ID) %>%
  summarise(MOVA_plus_3d_distance_SVM_predict = mean(predict))


 final_data <- merge(fpredictb, final_data)
#Save the final prediction.
 fwrite(final_data, MOVA_final_predict_file)
#Plot ROC curve
 df$iter3 <- df$iter1 + (df$iter2 * 5)
 out <- cvAUC(df$predict, df$result, label.ordering = NULL, folds = df$iter3)
 plot(out$perf, col = "blue", avg = "vertical", add=TRUE)
 cvauc <- out$cvAUC
 YI2 <- data.frame(matrix(rep(NA, 5), nrow=1))[numeric(0), ]
 colnames(YI2) <- c("Cutoff", "positive_variant_num", "negative_variant_num", "AUC", "cvauc")
 for(i in 6:30){
  pred <- prediction(df[df$iter3==i,]$predict, df[df$iter3==i,]$result)
  auc.tmp <- performance(pred,"auc")
  auc <- as.numeric(auc.tmp@y.values)
  tab <- data.frame(Cutoff=unlist(pred@cutoffs),
   TP=unlist(pred@tp), FP=unlist(pred@fp),
   FN=unlist(pred@fn), TN=unlist(pred@tn),
   Sensitivity=unlist(pred@tp)/(unlist(pred@tp)+unlist(pred@fn)),
   Specificity=unlist(pred@tn)/(unlist(pred@fp)+unlist(pred@tn)),
   Accuracy=((unlist(pred@tp)+unlist(pred@tn))/nrow(df)),
   Precision=(unlist(pred@tp)/(unlist(pred@tp)+unlist(pred@fp))))
  tab$Youden <- tab$Sensitivity+tab$Specificity-1
  YI <- tab[order(tab$Youden, decreasing=T), ]
  YI <- data.frame(Cutoff = YI$Cutoff[1], positive_variant_num = c(nrow(df[(df$result==1)&(df$iter2==1),])), negative_variant_num = c(nrow(df[(df$result==0)&(df$iter2==1),])),AUC = c(auc),cvauc = c(cvauc))
  YI2 <- rbind(YI2, YI)
 }
 fwrite(YI2, paste(protein_name, phenotype, "result_plus_3d_distance_SVM.csv", sep = "_"))
 df2 <- df %>%
  group_by(ID) %>%
  summarise(MOVA_plus_3d_distance_SVM_predict = mean(predict))
 D[ , c('result')] <- list(NULL)
 fwrite(merge(D,df2),  MOVA_predict_file)
 fwrite(df, paste(protein_name, phenotype, "predict_orig_plus_3d_distance_SVM.csv", sep = "_"))
 return(YI2)
}


MOVA_plus_3d_distance_xgboost <- function(protein_name, MOVA_final_predict_file, MOVA_predict_file, phenotype = "Target"){
 D <- fread(MOVA_predict_file)
 D$MOVA_plus_3d_distance_xgboost_predict<- 0
 D[ , c('MOVA_plus_3d_distance_xgboost_predict')] <- list(NULL)
#
 if(phenotype == "Target"){
  D <- D[(D$Type=="Target")|(D$Type=="Ctrl"), ]
 }
#
# Initialize data.frame.
 fpredict <- data.frame(matrix(rep(NA, 3), nrow=1))[numeric(0), ]
 colnames(fpredict) <- c("change", "predict", "n")
 imp.rfa <- data.frame(matrix(rep(NA, 3), nrow=1))[numeric(0), ]
 colnames(imp.rfa) <- c("name", "IncNodePurity", "n")
# Obtain information on each variant for final predict
 final_data <- fread(MOVA_final_predict_file)
 final_data$MOVA_plus_3d_distance_xgboost_predict<- 0
 final_data[ , c('MOVA_plus_3d_distance_xgboost_predict')] <- list(NULL)
 df <- data.frame(matrix(rep(NA, 6), nrow=1))[numeric(0), ]
 colnames(df) <- c("ID", "re_result", "distance",  "change", "iter1", "iter2")
 D$id2 <- 0
 if(phenotype == "Target"){
  D$result <- D$Type2
 } else {
  D$result <- D$Type3
 }
 D[ , c('id2')] <- list(NULL)
 oldw <- getOption("warn")
 options(warn=-1)
 D <- rowid_to_column(D, var = "id2")
 for(i2 in 1:5){## split data into train and test using stratified sampling
  D$id3 <- D$id2
  j <- 5
  dat1 <- D %>% stratified(., group = "result", size = 1/j)
  dat1$iter <- 1
  datzan <- D[-dat1$id2, ]
  for(i in 2:j){
   datzan[ , c('id2')] <- list(NULL)
   datzan <- rowid_to_column(datzan, var = "id2")
   if(i != j){
    dat2 <- datzan %>% stratified(., group = "result", size = 1/(j-i+1))
   } else {
    dat2 <- datzan
   }
   dat2$iter <- i
   datzan <- datzan[-dat2$id2, ]
   dat1 <- rbind(dat1, dat2)
  }
  dat1[ , c('id2')] <- list(NULL)
  dat1 <- rowid_to_column(dat1, var = "id2")
  for (i in 1:j){
   val <- dat1[dat1$iter == i,]
   train.data <- dat1[-val$id2,]
#Model construction with randomForest.
#Type2 is the objective variable and con+x+y+z+b is the explanatory variable.
    x <- train.data[train.data$result==1,]$x
    y <- train.data[train.data$result==1,]$y
    z <- train.data[train.data$result==1,]$z
    train.data$t_distance <- 0
    for(i3 in 1:nrow(train.data)){
      t_distance <- sqrt((x-train.data[i3]$x)*(x-train.data[i3]$x)+(y-train.data[i3]$y)*(y-train.data[i3]$y)+(z-train.data[i3]$z)*(z-train.data[i3]$z))
      t_distance <- sort(t_distance,decreasing=F)
      if(train.data[i3]$result == 1){
        train.data[i3]$t_distance <- t_distance[2]
      } else {
        train.data[i3]$t_distance <- t_distance[1]
      }
    }
    train.data.mx <- sparse.model.matrix(result~.,data.frame(result = train.data$result, con = train.data$con, x = train.data$x, y =train.data$y, z = train.data$z, b =train.data$b, t_distance =train.data$t_distance))
    train.data.dm <- xgb.DMatrix(train.data.mx,label=train.data$result)
    xgb.result <- xgb.train( data=train.data.dm, label=train.data$result, objective="binary:logistic", booster="gbtree", nrounds=100, verbose=1)
#Measure variable importance.
#Store the predicted values of the test data in data.frame.
    x <- train.data[train.data$result==1,]$x
    y <- train.data[train.data$result==1,]$y
    z <- train.data[train.data$result==1,]$z
    val$t_distance <- 0
    for(i3 in 1:nrow(val)){
      val[i3]$t_distance <- min(sqrt((x-val[i3]$x)*(x-val[i3]$x)+(y-val[i3]$y)*(y-val[i3]$y)+(z-val[i3]$z)*(z-val[i3]$z)))
    }
    test.data.mx <- sparse.model.matrix(result~.,data.frame(result = val$result, con = val$con, x = val$x, y =val$y, z = val$z, b =val$b, t_distance = val$t_distance))
    test.data.dm <- xgb.DMatrix(test.data.mx,label=val$result) 
    f <-data.frame(ID = val$ID, result = val$result, predict = predict(object=xgb.result,newdata=test.data.dm), change = val$change, iter1 = i, iter2 =i2)
   df <- rbind(df, f)
  }
 }

#Store the predicted values of all variants for final predict in data.frame.
   x <- D[D$result==1,]$x
   y <- D[D$result==1,]$y
   z <- D[D$result==1,]$z
   D$t_distance <- 0
   for(i3 in 1:nrow(D)){
     t_distance <- sqrt((x-D[i3]$x)*(x-D[i3]$x)+(y-D[i3]$y)*(y-D[i3]$y)+(z-D[i3]$z)*(z-D[i3]$z))
     t_distance <- sort(t_distance,decreasing=F)
     if(D[i3]$result == 1){
       D[i3]$t_distance <- t_distance[2]
     } else {
       D[i3]$t_distance <- t_distance[1]
     }

   }
  D.mx <- sparse.model.matrix(result~.,data.frame(result = D$result, con = D$con, x = D$x, y =D$y, z = D$z, b =D$b, t_distance = D$t_distance))
  D.dm <- xgb.DMatrix(D.mx,label=D$result)
 final_data.mx <- sparse.model.matrix(dammy~.,data.frame(dammy = final_data$dammy, con =  final_data$con, x =  final_data$x, y = final_data$y, z =  final_data$z, b = final_data$b, t_distance = final_data$t_distance))
 final_data.dm <- xgb.DMatrix(final_data.mx,label=final_data$dammy)
 for (i in 1:30){
   xgb.result <- xgb.train( data=D.dm, label=D$result, objective="binary:logistic", booster="gbtree", nrounds=100, verbose=1)
  f <-data.frame(ID = final_data$ID, predict = predict(object=xgb.result,newdata=final_data.dm), n = i) 
  fpredict <- rbind(fpredict, f)
 }

 options(warn = oldw)
 fpredictb <- fpredict %>%
  group_by(ID) %>%
 summarise(MOVA_plus_3d_distance_xgboost_predict = mean(predict))


 final_data <- merge(fpredictb, final_data)
#Save the final prediction.
 fwrite(final_data, MOVA_final_predict_file)
#Plot ROC curve
 df$iter3 <- df$iter1 + (df$iter2 * 5)
 out <- cvAUC(df$predict, df$result, label.ordering = NULL, folds = df$iter3)
 plot(out$perf, col = "black", avg = "vertical", add=TRUE)
 cvauc <- out$cvAUC
 YI2 <- data.frame(matrix(rep(NA, 5), nrow=1))[numeric(0), ]
 colnames(YI2) <- c("Cutoff", "positive_variant_num", "negative_variant_num", "AUC", "cvauc")
 for(i in 6:30){
  pred <- prediction(df[df$iter3==i,]$predict, df[df$iter3==i,]$result)
  auc.tmp <- performance(pred,"auc")
  auc <- as.numeric(auc.tmp@y.values)
  tab <- data.frame(Cutoff=unlist(pred@cutoffs),
   TP=unlist(pred@tp), FP=unlist(pred@fp),
   FN=unlist(pred@fn), TN=unlist(pred@tn),
   Sensitivity=unlist(pred@tp)/(unlist(pred@tp)+unlist(pred@fn)),
   Specificity=unlist(pred@tn)/(unlist(pred@fp)+unlist(pred@tn)),
   Accuracy=((unlist(pred@tp)+unlist(pred@tn))/nrow(df)),
   Precision=(unlist(pred@tp)/(unlist(pred@tp)+unlist(pred@fp))))
  tab$Youden <- tab$Sensitivity+tab$Specificity-1
  YI <- tab[order(tab$Youden, decreasing=T), ]
  YI <- data.frame(Cutoff = YI$Cutoff[1], positive_variant_num = c(nrow(df[(df$result==1)&(df$iter2==1),])), negative_variant_num = c(nrow(df[(df$result==0)&(df$iter2==1),])),AUC = c(auc),cvauc = c(cvauc))
  YI2 <- rbind(YI2, YI)
 }
 fwrite(YI2, paste(protein_name, phenotype, "result_plus_3d_distance_xgboost.csv", sep = "_"))
 df2 <- df %>%
  group_by(ID) %>%
  summarise(MOVA_plus_3d_distance_xgboost_predict = mean(predict))
 D[ , c('result')] <- list(NULL)
 fwrite(merge(D,df2),  MOVA_predict_file)
 fwrite(df, paste(protein_name, phenotype, "predict_orig_plus_3d_distance_xgboost.csv", sep = "_"))
 return(YI2)
}


MOVA_SVM <- function(protein_name, MOVA_final_predict_file, MOVA_predict_file, phenotype = "Target"){
 D <- fread(MOVA_predict_file)
 D$SVM_predict<- 0
 D[ , c('SVM_predict')] <- list(NULL)
#
 if(phenotype == "Target"){
  D <- D[(D$Type=="Target")|(D$Type=="Ctrl"), ]
 }
#
# Initialize data.frame.
 fpredict <- data.frame(matrix(rep(NA, 3), nrow=1))[numeric(0), ]
 colnames(fpredict) <- c("change", "predict", "n")
# Obtain information on each variant for final predict
 final_data <- fread(MOVA_final_predict_file)
 final_data$SVM_predict<- 0
 final_data[ , c('SVM_predict')] <- list(NULL)
 df <- data.frame(matrix(rep(NA, 6), nrow=1))[numeric(0), ]
 colnames(df) <- c("ID", "result", "predict",  "change", "iter1", "iter2")
 D$id2 <- 0
 D[ , c('id2')] <- list(NULL)
 if(phenotype == "Target"){
  D$result <- D$Type2
 } else {
  D$result <- D$Type3
 }
#SVM
 oldw <- getOption("warn")
 options(warn=-1)
 D <- rowid_to_column(D, var = "id2")
 for(i2 in 1:5){## split data into train and test using stratified sampling
  D$id3 <- D$id2
  j <- 5
  dat1 <- D %>% stratified(., group = "result", size = 1/j)
  dat1$iter <- 1
  datzan <- D[-dat1$id2, ]
  for(i in 2:j){
   datzan[ , c('id2')] <- list(NULL)
   datzan <- rowid_to_column(datzan, var = "id2")
   if(i != j){
    dat2 <- datzan %>% stratified(., group = "result", size = 1/(j-i+1))
   } else {
    dat2 <- datzan
   }
   dat2$iter <- i
   datzan <- datzan[-dat2$id2, ]
   dat1 <- rbind(dat1, dat2)
  }
  dat1[ , c('id2')] <- list(NULL)
  dat1 <- rowid_to_column(dat1, var = "id2")
  for (i in 1:j){
   val <- dat1[dat1$iter == i,]
   train.data <- dat1[-val$id2,]
#Model construction with randomForest.
#Type2 is the objective variable and con+x+y+z+b is the explanatory variable.
    svm<-ksvm(result~con+x+y+z+b,data=train.data)
#Store the predicted values of the test data in data.frame.
    f <-data.frame(ID = val$ID, result = val$result, predict = predict(svm,val), change = val$change, iter1 = i, iter2 =i2)
   df <- rbind(df, f)
  }
 }
#Store the predicted values of all variants for final predict in data.frame.
 for (i in 1:30){
   svm<-ksvm(result~con+x+y+z+b,data=D)
  f <-data.frame(ID = final_data$ID, predict = predict(svm,final_data), n = i) 
  fpredict <- rbind(fpredict, f)
 }
 options(warn = oldw)
#Take the average of the predictions to make the final prediction.
 fpredictb <- fpredict %>%
  group_by(ID) %>%
  summarise(SVM_predict = mean(predict))
#If the amino acid residues do not change, the predicted value is 0.
 ID <- str_split(fpredictb$ID, pattern = ":", simplify = TRUE)
 for (i in 1:nrow(fpredictb)){
  if(ID[i, 5]==ID[i, 7]){
   fpredictb$predict[i] = 0
  }
 }

 final_data <- merge(fpredictb, final_data)
#Save the final prediction.
 fwrite(final_data, MOVA_final_predict_file)
#Plot ROC curve
 df$iter3 <- df$iter1 + (df$iter2 * 5)
 out <- cvAUC(df$predict, df$result, label.ordering = NULL, folds = df$iter3)
 plot(out$perf, col = "black", avg = "vertical", add=TRUE)
 cvauc <- out$cvAUC
 YI2 <- data.frame(matrix(rep(NA, 5), nrow=1))[numeric(0), ]
 colnames(YI2) <- c("Cutoff", "positive_variant_num", "negative_variant_num", "AUC", "cvauc")
 for(i in 6:30){
  pred <- prediction(df[df$iter3==i,]$predict, df[df$iter3==i,]$result)
  auc.tmp <- performance(pred,"auc")
  auc <- as.numeric(auc.tmp@y.values)
  tab <- data.frame(Cutoff=unlist(pred@cutoffs),
   TP=unlist(pred@tp), FP=unlist(pred@fp),
   FN=unlist(pred@fn), TN=unlist(pred@tn),
   Sensitivity=unlist(pred@tp)/(unlist(pred@tp)+unlist(pred@fn)),
   Specificity=unlist(pred@tn)/(unlist(pred@fp)+unlist(pred@tn)),
   Accuracy=((unlist(pred@tp)+unlist(pred@tn))/nrow(df)),
   Precision=(unlist(pred@tp)/(unlist(pred@tp)+unlist(pred@fp))))
  tab$Youden <- tab$Sensitivity+tab$Specificity-1
  YI <- tab[order(tab$Youden, decreasing=T), ]
  YI <- data.frame(Cutoff = YI$Cutoff[1], positive_variant_num = c(nrow(df[(df$result==1)&(df$iter2==1),])), negative_variant_num = c(nrow(df[(df$result==0)&(df$iter2==1),])),AUC = c(auc),cvauc = c(cvauc))
  YI2 <- rbind(YI2, YI)
 }
 fwrite(YI2, paste(protein_name, phenotype, "result_svm.csv", sep = "_"))
 df2 <- df %>%
  group_by(ID) %>%
  summarise(SVM_predict = mean(predict))
 D[ , c('result')] <- list(NULL)
 fwrite(merge(D,df2), MOVA_predict_file)
 fwrite(df, paste(protein_name, phenotype, "predict_orig_SVM.csv", sep = "_"))
 return(YI2)
}


MOVA_xgboost <- function(protein_name, MOVA_final_predict_file, MOVA_predict_file, phenotype = "Target"){
 D <- fread(MOVA_predict_file)
 D$xgboost_predict<- 0
 D[ , c('xgboost_predict')] <- list(NULL)
#
 if(phenotype == "Target"){
  D <- D[(D$Type=="Target")|(D$Type=="Ctrl"), ]
 }
#
# Initialize data.frame.
 fpredict <- data.frame(matrix(rep(NA, 3), nrow=1))[numeric(0), ]
 colnames(fpredict) <- c("change", "predict", "n")
# Obtain information on each variant for final predict
 final_data <- fread(MOVA_final_predict_file)
 final_data$xgboost_predict<- 0
 final_data[ , c('xgboost_predict')] <- list(NULL)
 final_data$dammy <- 0
 df <- data.frame(matrix(rep(NA, 6), nrow=1))[numeric(0), ]
 colnames(df) <- c("ID", "result", "predict",  "change", "iter1", "iter2")
 D$id2 <- 0
 D[ , c('id2')] <- list(NULL)
 if(phenotype == "Target"){
  D$result <- D$Type2
 } else {
  D$result <- D$Type3
 }
#SVM
 oldw <- getOption("warn")
 options(warn=-1)
 D <- rowid_to_column(D, var = "id2")
 for(i2 in 1:5){## split data into train and test using stratified sampling
  D$id3 <- D$id2
  j <- 5
  dat1 <- D %>% stratified(., group = "result", size = 1/j)
  dat1$iter <- 1
  datzan <- D[-dat1$id2, ]
  for(i in 2:j){
   datzan[ , c('id2')] <- list(NULL)
   datzan <- rowid_to_column(datzan, var = "id2")
   if(i != j){
    dat2 <- datzan %>% stratified(., group = "result", size = 1/(j-i+1))
   } else {
    dat2 <- datzan
   }
   dat2$iter <- i
   datzan <- datzan[-dat2$id2, ]
   dat1 <- rbind(dat1, dat2)
  }
  dat1[ , c('id2')] <- list(NULL)
  dat1 <- rowid_to_column(dat1, var = "id2")
  for (i in 1:j){
   val <- dat1[dat1$iter == i,]
   train.data <- dat1[-val$id2,]
#Model construction with randomForest.
#Type2 is the objective variable and con+x+y+z+b is the explanatory variable.
    train.data.mx <- sparse.model.matrix(result~.,data.frame(result = train.data$result, con = train.data$con, x = train.data$x, y =train.data$y, z = train.data$z, b =train.data$b))
    train.data.dm <- xgb.DMatrix(train.data.mx,label=train.data$result)
    xgb.result <- xgb.train( data=train.data.dm, label=train.data$result, objective="binary:logistic", booster="gbtree", nrounds=100, verbose=1)
#Store the predicted values of the test data in data.frame.
    test.data.mx <- sparse.model.matrix(result~.,data.frame(result = val$result, con = val$con, x = val$x, y =val$y, z = val$z, b =val$b))
    test.data.dm <- xgb.DMatrix(test.data.mx,label=val$result) 
    f <-data.frame(ID = val$ID, result = val$result, predict = predict(object=xgb.result,newdata=test.data.dm), change = val$change, iter1 = i, iter2 =i2)
   df <- rbind(df, f)
#Store the variable importance in data.frame.
  }
 }
#Store the predicted values of all variants for final predict in data.frame.
 final_data.mx <- sparse.model.matrix(dammy~.,data.frame(dammy = final_data$dammy, con =  final_data$con, x =  final_data$x, y = final_data$y, z =  final_data$z, b = final_data$b))
 final_data.dm <- xgb.DMatrix(final_data.mx,label=final_data$dammy)
  D.mx <- sparse.model.matrix(result~.,data.frame(result = D$result, con = D$con, x = D$x, y =D$y, z = D$z, b =D$b))
  D.dm <- xgb.DMatrix(D.mx,label=D$result)
 for (i in 1:30){
   xgb.result <- xgb.train( data=D.dm, label=D$result, objective="binary:logistic", booster="gbtree", nrounds=100, verbose=1)
  f <-data.frame(ID = final_data$ID, predict = predict(object=xgb.result,newdata=final_data.dm), n = i) 
  fpredict <- rbind(fpredict, f)
 }
 options(warn = oldw)
#Take the average of the predictions to make the final prediction.
 fpredictb <- fpredict %>%
  group_by(ID) %>%
  summarise(xgboost_predict = mean(predict))


 final_data <- merge(fpredictb, final_data)
#Save the final prediction.
 fwrite(final_data, MOVA_final_predict_file)
#Plot ROC curve
 df$iter3 <- df$iter1 + (df$iter2 * 5)
 out <- cvAUC(df$predict, df$result, label.ordering = NULL, folds = df$iter3)
 plot(out$perf, col = "blue", avg = "vertical", add=TRUE)
 cvauc <- out$cvAUC
 YI2 <- data.frame(matrix(rep(NA, 5), nrow=1))[numeric(0), ]
 colnames(YI2) <- c("Cutoff", "positive_variant_num", "negative_variant_num", "AUC", "cvauc")
 for(i in 6:30){
  pred <- prediction(df[df$iter3==i,]$predict, df[df$iter3==i,]$result)
  auc.tmp <- performance(pred,"auc")
  auc <- as.numeric(auc.tmp@y.values)
  tab <- data.frame(Cutoff=unlist(pred@cutoffs),
   TP=unlist(pred@tp), FP=unlist(pred@fp),
   FN=unlist(pred@fn), TN=unlist(pred@tn),
   Sensitivity=unlist(pred@tp)/(unlist(pred@tp)+unlist(pred@fn)),
   Specificity=unlist(pred@tn)/(unlist(pred@fp)+unlist(pred@tn)),
   Accuracy=((unlist(pred@tp)+unlist(pred@tn))/nrow(df)),
   Precision=(unlist(pred@tp)/(unlist(pred@tp)+unlist(pred@fp))))
  tab$Youden <- tab$Sensitivity+tab$Specificity-1
  YI <- tab[order(tab$Youden, decreasing=T), ]
  YI <- data.frame(Cutoff = YI$Cutoff[1], positive_variant_num = c(nrow(df[(df$result==1)&(df$iter2==1),])), negative_variant_num = c(nrow(df[(df$result==0)&(df$iter2==1),])),AUC = c(auc),cvauc = c(cvauc))
  YI2 <- rbind(YI2, YI)
 }
 fwrite(YI2, paste(protein_name, phenotype, "result_xgboost.csv", sep = "_"))
 df2 <- df %>%
  group_by(ID) %>%
  summarise(xgboost_predict = mean(predict))
 D[ , c('result')] <- list(NULL)
 fwrite(merge(D,df2), MOVA_predict_file)
 fwrite(df, paste(protein_name, phenotype, "predict_orig_xgboost.csv", sep = "_"))
 return(YI2)
}



PolyPhen_MOVA <- function(file_name, phenotype = "Target", col = "orange"){
 df <- fread(file_name)
 #
 if(phenotype == "Target"){
  df <- df[(df$Type=="Target")|(df$Type=="Ctrl"), ]
 }
 #
 df$predict <- df$pph_prob
 if(phenotype == "Target"){
  df$result <- df$Type2
 } else {
  df$result <- df$Type3
 }
 #
 pred <- prediction(df$predict, df$result)
 perf <- performance(pred, "tpr", "fpr")
 plot(perf, col = col, add=TRUE)
 auc.tmp <- performance(pred,"auc")
 auc <- as.numeric(auc.tmp@y.values)
 return(auc)
}

AlphScore_MOVA <- function(file_name, phenotype = "Target", col = "green"){
 df <- fread(file_name)
 #
 if(phenotype == "Target"){
 df <- df[(df$Type=="Target")|(df$Type=="Ctrl"), ]
 }
 df$predict <- df$AlphScore
 if(phenotype == "Target"){
  df$result <- df$Type2
 } else {
  df$result <- df$Type3
 }
 #
 pred <- prediction(df$predict, df$result)
 perf <- performance(pred, "tpr", "fpr")
 plot(perf, col = col, add=TRUE)
 auc.tmp <- performance(pred,"auc")
 auc <- as.numeric(auc.tmp@y.values)
 return(auc)
}

CADD_MOVA <- function(file_name, phenotype = "Target", col = "blue"){
 df <- fread(file_name)
 #
 if(phenotype == "Target"){
  df <- df[(df$Type=="Target")|(df$Type=="Ctrl"), ]
 }
 #
 df$predict <- df$CADD_raw
 if(phenotype == "Target"){
  df$result <- df$Type2
 } else {
  df$result <- df$Type3
 }
 #
 pred <- prediction(df$predict, df$result)
 perf <- performance(pred, "tpr", "fpr")
 plot(perf, col = col, add=TRUE)
 auc.tmp <- performance(pred,"auc")
 auc <- as.numeric(auc.tmp@y.values)
 return(auc)
}

CADD_AlphScore <- function(file_name, phenotype = "Target", col = "gray"){
 df <- fread(file_name)
 #
 if(phenotype == "Target"){
  df <- df[(df$Type=="Target")|(df$Type=="Ctrl"), ]
 }
 #
 df$predict <- df$glm_AlphCadd
 if(phenotype == "Target"){
  df$result <- df$Type2
 } else {
  df$result <- df$Type3
 }
 #
 pred <- prediction(df$predict, df$result)
 perf <- performance(pred, "tpr", "fpr")
 plot(perf, col = col, add=TRUE)
 auc.tmp <- performance(pred,"auc")
 auc <- as.numeric(auc.tmp@y.values)
 return(auc)
}

Com_CADD_MOVA <- function(protein_name, MOVA_final_predict_file, MOVA_predict_orig_file, MOVA_predict_file, alphscore = 0, phenotype = "Target"){
 final_data <- fread(MOVA_final_predict_file)
#MOVA
 df <- fread(MOVA_predict_orig_file)
 out <- cvAUC(df$predict, df$result, label.ordering = NULL, folds = df$iter3)
 plot(out$perf, col = "blue", avg = "vertical")
 auc_bind <- data.frame(MOVA = c(out$cvAUC))

 df <- fread(MOVA_predict_file)
 df$CADDMOVA_predict <- 0
 df[ , c('CADDMOVA_predict')] <- list(NULL)
 if(phenotype == "Target"){
  df <- df[(df$Type=="Target")|(df$Type=="Ctrl"), ]
 }
 if(phenotype == "Target"){
  df$result <- df$Type2
 } else {
  df$result <- df$Type3
 }
#CADD
 df$predict <- df$CADD_raw
 pred <- prediction(df$predict, df$result)
 perf <- performance(pred, "tpr", "fpr")
 plot(perf, col = "black", add=TRUE)
 auc.tmp <- performance(pred,"auc")
 auc <- as.numeric(auc.tmp@y.values)
 auc_bind <- cbind(auc_bind, data.frame(CADD = c(auc)))
 if(alphscore == 1){
#Alph_Score+CADD
  df$predict <- df$glm_AlphCadd
#
  pred <- prediction(df$predict, df$result)
  perf <- performance(pred, "tpr", "fpr")
  plot(perf, col = "gray", add=TRUE)
  auc.tmp <- performance(pred,"auc")
  auc <- as.numeric(auc.tmp@y.values)
  auc_bind <- cbind(auc_bind, data.frame(com_Alph_CADD = c(auc)))
 }
#MOVA+CADD
 D <- df
 df <- data.frame(matrix(rep(NA, 6), nrow=1))[numeric(0), ]
 colnames(df) <- c("ID", "result", "CADDMOVA",  "change", "iter1", "iter2")
 D$id2 <- 0
 D[ , c('id2')] <- list(NULL)
 oldw <- getOption("warn")
 options(warn=-1)
 D <- rowid_to_column(D, var = "id2")
 for(i2 in 1:5){## split data into train and test using stratified sampling
  D$id3 <- D$id2
  j <- 5
  dat1 <- D %>% stratified(., group = "result", size = 1/j)
  dat1$iter <- 1
  datzan <- D[-dat1$id2, ]
  for(i in 2:j){
   datzan[ , c('id2')] <- list(NULL)
   datzan <- rowid_to_column(datzan, var = "id2")
   if(i != j){
    dat2 <- datzan %>% stratified(., group = "result", size = 1/(j-i+1))
   } else {
    dat2 <- datzan
   }
   dat2$iter <- i
   datzan <- datzan[-dat2$id2, ]
   dat1 <- rbind(dat1, dat2)
  }
  dat1[ , c('id2')] <- list(NULL)
  dat1 <- rowid_to_column(dat1, var = "id2")
  for (i in 1:j){
   val <- dat1[dat1$iter == i,]
   train.data <- dat1[-val$id2,]
   glmals <-glm(result~MOVA_predict+CADD_raw,data=train.data, family=binomial)
#Store the predicted values of the test data in data.frame.
    f <-data.frame(ID = val$ID, result = val$result, predict = predict(glmals,newdata = val), change = val$change, iter1 = i, iter2 =i2)
   df <- rbind(df, f)
  }
 }
   glmals <-glm(result~MOVA_predict+CADD_raw,data=D, family=binomial)
  f <-data.frame(ID = final_data$ID, CADDMOVA = predict(glmals,newdata = final_data))

 options(warn = oldw)

#Take the average of variable importance


 final_data <- merge(f, final_data)
#Save the final prediction.
 fwrite(final_data, MOVA_final_predict_file)
#Plot ROC curve
 df$iter3 <- df$iter1 + (df$iter2 * 5)
 out <- cvAUC(df$predict, df$result, label.ordering = NULL, folds = df$iter3)
 plot(out$perf, col = "red", avg = "vertical", add=TRUE)
 cvauc <- out$cvAUC
 auc_bind <- cbind(auc_bind, data.frame(com_MOVA_CADD = c(cvauc)))
 YI2 <- data.frame(matrix(rep(NA, 5), nrow=1))[numeric(0), ]
 colnames(YI2) <- c("Cutoff", "positive_variant_num", "negative_variant_num", "AUC", "cvauc")
 for(i in 6:30){
  pred <- prediction(df[df$iter3==i,]$predict, df[df$iter3==i,]$result)
  auc.tmp <- performance(pred,"auc")
  auc <- as.numeric(auc.tmp@y.values)
  tab <- data.frame(Cutoff=unlist(pred@cutoffs),
   TP=unlist(pred@tp), FP=unlist(pred@fp),
   FN=unlist(pred@fn), TN=unlist(pred@tn),
   Sensitivity=unlist(pred@tp)/(unlist(pred@tp)+unlist(pred@fn)),
   Specificity=unlist(pred@tn)/(unlist(pred@fp)+unlist(pred@tn)),
   Accuracy=((unlist(pred@tp)+unlist(pred@tn))/nrow(df)),
   Precision=(unlist(pred@tp)/(unlist(pred@tp)+unlist(pred@fp))))
  tab$Youden <- tab$Sensitivity+tab$Specificity-1
  YI <- tab[order(tab$Youden, decreasing=T), ]
  YI <- data.frame(Cutoff = YI$Cutoff[1], positive_variant_num = c(nrow(df[(df$result==1)&(df$iter2==1),])), negative_variant_num = c(nrow(df[(df$result==0)&(df$iter2==1),])),AUC = c(auc))
  YI <- cbind(YI, auc_bind)
  YI2 <- rbind(YI2, YI)
 }
 fwrite(YI2, paste(protein_name, phenotype, "result_CADDMOVA.csv", sep = "_"))
 df2 <- df %>%
  group_by(ID) %>%
  summarise(CADDMOVA_predict = mean(CADDMOVA))
 D[ , c('result')] <- list(NULL)
 fwrite(merge(D,df2),  MOVA_predict_file)
 fwrite(df, paste(protein_name, phenotype, "predict_orig_CADDMOVA.csv", sep = "_"))
 return(YI2)
}


REVEL_MOVA <- function(file_name, phenotype = "Target", col = "black"){
 df <- fread(file_name)
#
 if(phenotype == "Target"){
  df <- df[(df$Type=="Target")|(df$Type=="Ctrl"), ]
 }
#
 df$predict <- df$REVEL_score
 df <- df[df$predict!=".",]
 df$predict <- as.numeric(df$predict)
 if(phenotype == "Target"){
  df$result <- df$Type2
 } else {
  df$result <- df$Type3
 }
#
 pred <- prediction(df$predict, df$result)
 perf <- performance(pred, "tpr", "fpr")
 plot(perf, col = col, add=TRUE)
 auc.tmp <- performance(pred,"auc")
 auc <- as.numeric(auc.tmp@y.values)
 return(auc)
}

REVEL_AlphScore <- function(file_name, phenotype = "Target", col = "gray"){
 df <- fread(file_name)
#
 if(phenotype == "Target"){
  df <- df[(df$Type=="Target")|(df$Type=="Ctrl"), ]
 }
#
 df$predict <- df$glm_AlphRevel
 df <- df[df$predict!=".",]
 df$predict <- as.numeric(df$predict)
 if(phenotype == "Target"){
  df$result <- df$Type2
 } else {
  df$result <- df$Type3
 }
#
 pred <- prediction(df$predict, df$result)
 perf <- performance(pred, "tpr", "fpr")
 plot(perf, col = col, add=TRUE)
 auc.tmp <- performance(pred,"auc")
 auc <- as.numeric(auc.tmp@y.values)
 return(auc)
}


Com_REVEL_MOVA <- function(protein_name, MOVA_final_predict_file, MOVA_predict_orig_file, MOVA_predict_file, alphscore = 0, phenotype = "Target"){
 final_data <- fread(MOVA_final_predict_file)
#MOVA
 df <- fread(MOVA_predict_orig_file)
 out <- cvAUC(df$predict, df$result, label.ordering = NULL, folds = df$iter3)
 plot(out$perf, col = "blue", avg = "vertical")
 auc_bind <- data.frame(MOVA = c(out$cvAUC))

 df <- fread(MOVA_predict_file)
 df$REVELMOVA_predict <- 0
 df[ , c('REVELMOVA_predict')] <- list(NULL)

 if(phenotype == "Target"){
  df <- df[(df$Type=="Target")|(df$Type=="Ctrl"), ]
 }
 if(phenotype == "Target"){
  df$result <- df$Type2
 } else {
  df$result <- df$Type3
 }
#REVEL
 df$predict <- df$REVEL_score
 pred <- prediction(df$predict, df$result)
 perf <- performance(pred, "tpr", "fpr")
 plot(perf, col = "black", add=TRUE)
 auc.tmp <- performance(pred,"auc")
 auc <- as.numeric(auc.tmp@y.values)
 auc_bind <- cbind(auc_bind, data.frame(REVEL = c(auc)))
 if(alphscore == 1){
#Alph_Score+REVEL
  df$predict <- df$glm_AlphRevel
#
  pred <- prediction(df$predict, df$result)
  perf <- performance(pred, "tpr", "fpr")
  plot(perf, col = "gray", add=TRUE)
  auc.tmp <- performance(pred,"auc")
  auc <- as.numeric(auc.tmp@y.values)
  auc_bind <- cbind(auc_bind, data.frame(com_Alph_REVEL = c(auc)))
 }
#MOVA+REVEL
 D <- df
 df <- data.frame(matrix(rep(NA, 6), nrow=1))[numeric(0), ]
 colnames(df) <- c("ID", "result", "REVELMOVA",  "change", "iter1", "iter2")
 D$id2 <- 0
 D[ , c('id2')] <- list(NULL)
 oldw <- getOption("warn")
 options(warn=-1)
 D <- rowid_to_column(D, var = "id2")
 for(i2 in 1:5){## split data into train and test using stratified sampling
  D$id3 <- D$id2
  j <- 5
  dat1 <- D %>% stratified(., group = "result", size = 1/j)
  dat1$iter <- 1
  datzan <- D[-dat1$id2, ]
  for(i in 2:j){
   datzan[ , c('id2')] <- list(NULL)
   datzan <- rowid_to_column(datzan, var = "id2")
   if(i != j){
    dat2 <- datzan %>% stratified(., group = "result", size = 1/(j-i+1))
   } else {
    dat2 <- datzan
   }
   dat2$iter <- i
   datzan <- datzan[-dat2$id2, ]
   dat1 <- rbind(dat1, dat2)
  }
  dat1[ , c('id2')] <- list(NULL)
  dat1 <- rowid_to_column(dat1, var = "id2")
  for (i in 1:j){
   val <- dat1[dat1$iter == i,]
   train.data <- dat1[-val$id2,]
   glmals <-glm(result~MOVA_predict+REVEL_score,data=train.data, family=binomial)
#Store the predicted values of the test data in data.frame.
    f <-data.frame(ID = val$ID, result = val$result, predict = predict(glmals,newdata = val), change = val$change, iter1 = i, iter2 =i2)
   df <- rbind(df, f)
  }
 }
   glmals <-glm(result~MOVA_predict+REVEL_score,data=D, family=binomial)
  f <-data.frame(ID = final_data$ID, REVELMOVA = predict(glmals,newdata = final_data))

 options(warn = oldw)

#Take the average of variable importance


 final_data <- merge(f, final_data)
#Save the final prediction.
 fwrite(final_data, MOVA_final_predict_file)
#Plot ROC curve
 df$iter3 <- df$iter1 + (df$iter2 * 5)
 out <- cvAUC(df$predict, df$result, label.ordering = NULL, folds = df$iter3)
 plot(out$perf, col = "red", avg = "vertical", add=TRUE)
 cvauc <- out$cvAUC
 auc_bind <- cbind(auc_bind, data.frame(com_MOVA_REVEL = c(cvauc)))
 YI2 <- data.frame(matrix(rep(NA, 5), nrow=1))[numeric(0), ]
 colnames(YI2) <- c("Cutoff", "positive_variant_num", "negative_variant_num", "AUC", "cvauc")
 for(i in 6:30){
  pred <- prediction(df[df$iter3==i,]$predict, df[df$iter3==i,]$result)
  auc.tmp <- performance(pred,"auc")
  auc <- as.numeric(auc.tmp@y.values)
  tab <- data.frame(Cutoff=unlist(pred@cutoffs),
   TP=unlist(pred@tp), FP=unlist(pred@fp),
   FN=unlist(pred@fn), TN=unlist(pred@tn),
   Sensitivity=unlist(pred@tp)/(unlist(pred@tp)+unlist(pred@fn)),
   Specificity=unlist(pred@tn)/(unlist(pred@fp)+unlist(pred@tn)),
   Accuracy=((unlist(pred@tp)+unlist(pred@tn))/nrow(df)),
   Precision=(unlist(pred@tp)/(unlist(pred@tp)+unlist(pred@fp))))
  tab$Youden <- tab$Sensitivity+tab$Specificity-1
  YI <- tab[order(tab$Youden, decreasing=T), ]
  YI <- data.frame(Cutoff = YI$Cutoff[1], positive_variant_num = c(nrow(df[(df$result==1)&(df$iter2==1),])), negative_variant_num = c(nrow(df[(df$result==0)&(df$iter2==1),])),AUC = c(auc))
  YI <- cbind(YI, auc_bind)
  YI2 <- rbind(YI2, YI)
 }
 fwrite(YI2, paste(protein_name, phenotype, "result_REVELMOVA.csv", sep = "_"))
 df2 <- df %>%
  group_by(ID) %>%
  summarise(REVELMOVA_predict = mean(REVELMOVA))
 D[ , c('result')] <- list(NULL)
 fwrite(merge(D,df2),  MOVA_predict_file)
 fwrite(df, paste(protein_name, phenotype, "predict_orig_REVELMOVA.csv", sep = "_"))
 return(YI2)
}


EVE_Final <- function(file_name, EVE_file_name){
 df <- fread(file_name)
 y <- fread(EVE_file_name)
 df[, EVE:=NA]
 df[, EVE_class:=NA]
 df$Pos <- as.integer(substr(df$change, 2, nchar(df$change)-1))
 for (i in 1:nrow(df)){
  df$EVE[i] <- y[(y$position==df$Pos[i])&(y$mt_aa==df$aaalt[i])]$EVE_scores_ASM
  df$EVE_class[i] <- y[(y$position==df$Pos[i])&(y$mt_aa==df$aaalt[i])]$EVE_classes_75_pct_retained_ASM
 }
 fwrite(df, file_name)
}

EVE_MOVA <- function(file_name, EVE_file_name, phenotype = "Target"){
 df <- fread(file_name)
 y <- fread(EVE_file_name)
#
 if(phenotype == "Target"){
  df <- df[(df$Type=="Target")|(df$Type=="Ctrl"), ]
 }
#
 df[, EVE:=NA]
 df[, EVE_class:=NA]
 df$Pos <- as.integer(substr(df$change, 2, nchar(df$change)-1))
 for (i in 1:nrow(df)){
  df$EVE[i] <- y[(y$position==df$Pos[i])&(y$mt_aa==df$aaalt[i])]$EVE_scores_ASM
  df$EVE_class[i] <- y[(y$position==df$Pos[i])&(y$mt_aa==df$aaalt[i])]$EVE_classes_75_pct_retained_ASM
 }
 fwrite(df, file_name)
 df$predict <- df$EVE
 if(phenotype == "Target"){
  df$result <- df$Type2
 } else {
  df$result <- df$Type3
 }
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

Edit_final_variant_file <- function(input_file_name, export_file_name){
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
 x$ref <- substr(x$HGVSp_VEP_split, 3, 5)
 x$Pos <- substr(x$HGVSp_VEP_split, 6, nchar(x$HGVSp_VEP_split)-3)
 x$alt <- substr(x$HGVSp_VEP_split, nchar(x$HGVSp_VEP_split)-2, nchar(x$HGVSp_VEP_split))
 for (i in 1:nrow(x)){
  x$change[i] <- paste(amino_code1[x$ref[i]], x$Pos[i], amino_code1[x$alt[i]], sep = "")
 }
 fwrite(x, export_file_name)
}

Edit_variant_data <- function(AlphScore_final_file, hgmd_file_name, gnomad_file_name, uniprot_name, gene_name, export_file_name){
#https://zenodo.org/record/6288139
 x <- fread(AlphScore_final_file)
 x <- x[x$Uniprot_acc_split==uniprot_name,]
 colnames(x)[1] <- "Chromosome"
 colnames(x)[2] <- "Position"
 colnames(x)[3] <- "Reference"
 colnames(x)[4] <- "Alternate"
 x[, DEOGEN2_score:=NULL]
 x[, b_factor:=NULL]
 x[, SOLVENT_ACCESSIBILITY_core:=NULL]
 x[, in_gnomad_train:=NULL]
 x[, in_clinvar_ds:=NULL]
 x[, glm_AlphCaddDeogen:=NULL]
 x[, glm_AlphDeogenRevel:=NULL]
 x[, glm_DeogenRevel:=NULL]
 x[, glm_CaddDeogen:=NULL]
 x[, glm_AlphDeogen:=NULL]
 x[, glm_AlphRevelCadd:=NULL]
 x[, glm_CaddDeogenRevel:=NULL]
 fwrite(x, paste(gene_name, "_alph.csv", sep = ""))

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
 D <- merge(x,y)
 D$Type2 <- 0
 D$Type3 <- 1
 D[D$Type=="Target",]$Type2 <- 1
 D[D$Type=="Target",]$Type3 <- 1
 D[D$Type=="Ctrl",]$Type2 <- 0
 D[D$Type=="Ctrl",]$Type3 <- 0
 fwrite(D, export_file_name)
}

Edit_polyphen_data <- function(dbNSFP4.3a_final_file, input_file_name1,  input_file_name2,  uniprot_name, position1, position2, export_file_name1, export_file_name2){
 D <- fread(dbNSFP4.3a_final_file)
 head(D)
 colnames(D)[1] <- "Chromosome"
 colnames(D)[2] <- "Position"
 colnames(D)[3] <- "Reference"
 colnames(D)[4] <- "Alternate"
 dren <- D[(D$Position<=position2)&(D$Position>=position1),]
 dt <- data.table(Chromosome=dren$Chromosome, Position=dren$Position, Reference=dren$Reference, Alternate=dren$Alternate, Uniprot_acc=dren$Uniprot_acc, Polyphen_prob=dren$Polyphen2_HDIV_score)

 Uniplot_v <- str_split(dt$Uniprot_acc, pattern = ";", simplify = TRUE)
 Polyphen2_HDIV_v <- str_split(dt$Polyphen_prob, pattern = ";", simplify = TRUE)
 dt$pph_prob <- ""
 i <- 1
 while (i <nrow(dt)){
  CHm <- match(uniprot_name, Uniplot_v[i,])
  if(is.na(CHm)==FALSE){
   dt$pph_prob[i] <- Polyphen2_HDIV_v[i, CHm]
   i<-i+1
  } else {
   dt <- dt[-i, ] 
   Uniplot_v <- Uniplot_v[-i,]
   Polyphen2_HDIV_v <- Polyphen2_HDIV_v[-i, ]
  }
 }
 dt[, Uniprot_acc:=NULL]
 dt[, Polyphen_prob:=NULL]

 dat <- fread(input_file_name1)
 pph_data <- merge(dat, dt)
 fwrite(pph_data, export_file_name1)
 dat <- fread(input_file_name2)
 pph_data <- merge(dat, dt)
 fwrite(pph_data, export_file_name2)
}

Hgmd_divide <- function(hgmd_file_name){
 x <- fread(hgmd_file_name)
 x <- x %>% separate(VCF, c("buf1", "buf2"), sep=":")
 x <- x %>% separate(buf2, c("buf3", "Alternate"), sep="/")
 x$Chromosome <- str_sub(x$buf1, 4, -1)
 x$Position <- str_sub(x$buf3, 1,-2)
 x$Reference <- str_sub(x$buf3, -1, -1)
 x[, buf1:=NULL]
 x[, buf3:=NULL]
 fwrite(x, hgmd_file_name)
}

variant_predict_plot <- function(file_name, legend_type=FALSE){
 df <- fread( file_name)
 df$Typevariant <- "negative"
 df[df$Type=="Target", Typevariant := "positive"]
 f <- data.table(variant=df$Typevariant, predict=df$MOVA_predict, position=df$Pos, type="MOVA")
 f2 <- data.table(variant=df$Typevariant, predict=df$pph_prob, position=df$Pos, type="Polyphen-2")
 df2 <- rbind(f,f2)
 df2 <- transform(df2, variant= factor(variant, levels = c("positive","negative")))
 g <- ggplot(data = df2) 
 g<-g+theme(axis.text=element_text(size=12))
 g<-g+theme(legend.text=element_text(size=12))
 g<-g+geom_point(aes(x=position, y = predict, colour = variant))
 g<-g+facet_grid(rows = vars(type), scales="free_y")
 if(legend_type==TRUE){ 
  g<-g+theme(legend.position="bottom")
 } else {
  g<-g+theme(legend.position="none")
 }
 g
}