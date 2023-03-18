\name{MOVA_SVM}
\alias{MOVA_SVM}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
The training was performed in SVM.
%%  ~~function to do ... ~~
}
\description{
This function is trained in a Support vector machine with ΔBLOSUM62, 3D location information, and pLDDT as explanatory variables and the pathogenicity of the variant as objective variable.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
MOVA_SVM(protein_name, MOVA_final_predict_file, MOVA_predict_file, phenotype = "Target")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{protein_name}{
       Name of target protein/gene. Used to name the file to be exported.
%%     ~~Describe \code{protein_name} here~~
}
  \item{MOVA_final_predict_file}{
       The "gene name_Target or Pathogenic_finalpredict.csv" file output by the MOVA function. The final MOVA_3d_distance predicted values will be in this file.
%%     ~~Describe \code{MOVA_final_predict_file} here~~
}
  \item{MOVA_predict_file}{
       The "gene name_Target or Pathogenic_predict.csv" file output by the MOVA function.
%%     ~~Describe \code{MOVA_predict_file} here~~
}
  \item{phenotype}{
       Specify "Target" if you want the positive variant to be the Target variant only, or "Pathogenic" if you want to include the Pathogenic variant as well. The default is "Target".
%%     ~~Describe \code{phenotype} here~~
}
}
\details{
        Draw ROC curve for MOVA_SVM (black). MOVA_SVM is evaluated using the Stratified 5-fold cross validation method, and the models are repeated five more times. The average of the predicted values is added to the "SVM_predict" column of the MOVA_predict_file. The final MOVA_SVM predicted value are added to the "SVM_predict" column in MOVA_final_predict_file. "protein/gene name_Target or Pathogenic_result_svm.csv" contains the Cutoff value (Youden index) for each fold (Column: Cutoff), the number of positive variants for each gene used in the analysis (Column: positive_variant_num), number of negative variants (Column: negative_variant_num), AUC for each fold of MOVA_SVM (Column: AUC), cvAUC for MOVA_SVM (Column: cvauc). The file required for redrawing with the MOVA_redraw function of MOVA_SVM is output in "protein/gene name_Target or Pathogenic_predict_orig_SVM.csv".
%%  ~~ If necessary, more details than the description above ~~
}
\value{
 A data.frame containing the Cutoff value (Youden index) for each fold (Column: Cutoff), the number of positive variants for each gene used in the analysis (Column: positive_variant_num), number of negative variants (Column: negative_variant_num), AUC for each fold of MOVA_SVM (Column: AUC), cvAUC for MOVA_SVM (Column: cvauc).
The same data is output to "protein/gene name_Target or Pathogenic_result_svm.csv".
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Yuya Hatano
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
Hgmd_divide("./source/TARDBP.csv")
Edit_gnomAD_file("./source/gnomAD_v3.1.2_ENST00000240185_2023_02_14_13_23_48.csv", "./source/TARDBP_gnomAD.csv")
Edit_variant_data("../CADD_REVEL/AlphScore_final.tsv", "./source/TARDBP.csv", "./source/TARDBP_gnomAD.csv", "Q13148", "TARDBP", "./source/TARDBPvariantdata.csv")
Edit_polyphen_data("../dbNSFP/dbNSFP4.3a/dbNSFP4.3a_variant.chr1","./source/TARDBPvariantdata.csv", "TARDBP_alph.csv",  "Q13148",  11012654, 11025492 ,"./source/TARDBPvariantdatapol.csv", "./source/TARDBP_alphpol.csv")
Edit_final_variant_file("./source/TARDBP_alphpol.csv","./source/TARDBP_alphpol2.csv")
MOVA("./source/Q13148.fa", "TARDBP", "./source/AF-Q13148-F1-model_v2.pdb","./source/TARDBP_alphpol2.csv", "./source/TARDBPvariantdatapol.csv")
MOVA_SVM("TARDBP", "TARDBP_Target_finalpredict.csv","TARDBP_Target_predict.csv")
MOVA_xgboost("TARDBP", "TARDBP_Target_finalpredict.csv","TARDBP_Target_predict.csv")

## The function is currently defined as
function (protein_name, MOVA_final_predict_file, MOVA_predict_file, 
    phenotype = "Target") 
{
    D <- fread(MOVA_predict_file)
    D$SVM_predict <- 0
    D[, c("SVM_predict")] <- list(NULL)
    if (phenotype == "Target") {
        D <- D[(D$Type == "Target") | (D$Type == "Ctrl"), ]
    }
    fpredict <- data.frame(matrix(rep(NA, 3), nrow = 1))[numeric(0), 
        ]
    colnames(fpredict) <- c("change", "predict", "n")
    final_data <- fread(MOVA_final_predict_file)
    final_data$SVM_predict <- 0
    final_data[, c("SVM_predict")] <- list(NULL)
    df <- data.frame(matrix(rep(NA, 6), nrow = 1))[numeric(0), 
        ]
    colnames(df) <- c("ID", "result", "predict", "change", "iter1", 
        "iter2")
    D$id2 <- 0
    D[, c("id2")] <- list(NULL)
    if (phenotype == "Target") {
        D$result <- D$Type2
    }
    else {
        D$result <- D$Type3
    }
    oldw <- getOption("warn")
    options(warn = -1)
    D <- rowid_to_column(D, var = "id2")
    for (i2 in 1:5) {
        D$id3 <- D$id2
        j <- 5
        dat1 <- D \%>\% stratified(., group = "result", size = 1/j)
        dat1$iter <- 1
        datzan <- D[-dat1$id2, ]
        for (i in 2:j) {
            datzan[, c("id2")] <- list(NULL)
            datzan <- rowid_to_column(datzan, var = "id2")
            if (i != j) {
                dat2 <- datzan \%>\% stratified(., group = "result", 
                  size = 1/(j - i + 1))
            }
            else {
                dat2 <- datzan
            }
            dat2$iter <- i
            datzan <- datzan[-dat2$id2, ]
            dat1 <- rbind(dat1, dat2)
        }
        dat1[, c("id2")] <- list(NULL)
        dat1 <- rowid_to_column(dat1, var = "id2")
        for (i in 1:j) {
            val <- dat1[dat1$iter == i, ]
            train.data <- dat1[-val$id2, ]
            svm <- ksvm(result ~ con + x + y + z + b, data = train.data)
            f <- data.frame(ID = val$ID, result = val$result, 
                predict = predict(svm, val), change = val$change, 
                iter1 = i, iter2 = i2)
            df <- rbind(df, f)
        }
    }
    for (i in 1:30) {
        svm <- ksvm(result ~ con + x + y + z + b, data = D)
        f <- data.frame(ID = final_data$ID, predict = predict(svm, 
            final_data), n = i)
        fpredict <- rbind(fpredict, f)
    }
    options(warn = oldw)
    fpredictb <- fpredict \%>\% group_by(ID) \%>\% summarise(SVM_predict = mean(predict))
    ID <- str_split(fpredictb$ID, pattern = ":", simplify = TRUE)
    for (i in 1:nrow(fpredictb)) {
        if (ID[i, 5] == ID[i, 7]) {
            fpredictb$predict[i] = 0
        }
    }
    final_data <- merge(fpredictb, final_data)
    fwrite(final_data, MOVA_final_predict_file)
    df$iter3 <- df$iter1 + (df$iter2 * 5)
    out <- cvAUC(df$predict, df$result, label.ordering = NULL, 
        folds = df$iter3)
    plot(out$perf, col = "black", avg = "vertical", add = TRUE)
    cvauc <- out$cvAUC
    YI2 <- data.frame(matrix(rep(NA, 5), nrow = 1))[numeric(0), 
        ]
    colnames(YI2) <- c("Cutoff", "positive_variant_num", "negative_variant_num", 
        "AUC", "cvauc")
    for (i in 6:30) {
        pred <- prediction(df[df$iter3 == i, ]$predict, df[df$iter3 == 
            i, ]$result)
        auc.tmp <- performance(pred, "auc")
        auc <- as.numeric(auc.tmp@y.values)
        tab <- data.frame(Cutoff = unlist(pred@cutoffs), TP = unlist(pred@tp), 
            FP = unlist(pred@fp), FN = unlist(pred@fn), TN = unlist(pred@tn), 
            Sensitivity = unlist(pred@tp)/(unlist(pred@tp) + 
                unlist(pred@fn)), Specificity = unlist(pred@tn)/(unlist(pred@fp) + 
                unlist(pred@tn)), Accuracy = ((unlist(pred@tp) + 
                unlist(pred@tn))/nrow(df)), Precision = (unlist(pred@tp)/(unlist(pred@tp) + 
                unlist(pred@fp))))
        tab$Youden <- tab$Sensitivity + tab$Specificity - 1
        YI <- tab[order(tab$Youden, decreasing = T), ]
        YI <- data.frame(Cutoff = YI$Cutoff[1], positive_variant_num = c(nrow(df[(df$result == 
            1) & (df$iter2 == 1), ])), negative_variant_num = c(nrow(df[(df$result == 
            0) & (df$iter2 == 1), ])), AUC = c(auc), cvauc = c(cvauc))
        YI2 <- rbind(YI2, YI)
    }
    fwrite(YI2, paste(protein_name, phenotype, "result_svm.csv", 
        sep = "_"))
    df2 <- df \%>\% group_by(ID) \%>\% summarise(SVM_predict = mean(predict))
    D[, c("result")] <- list(NULL)
    fwrite(merge(D, df2), MOVA_predict_file)
    fwrite(df, paste(protein_name, phenotype, "predict_orig_SVM.csv", 
        sep = "_"))
    return(YI2)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
