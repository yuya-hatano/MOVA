\name{Com_CADD_MOVA}
\alias{Com_CADD_MOVA}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Combining CADD and MOVA
%%  ~~function to do ... ~~
}
\description{
MOVA was combined with CADD using logistic regression as implemented in the R-function glm with the option family=binomial  and evaluated using the stratified 5-fold cross validation method. MOVA + CADD was compared to MOVA, CADD, and CADD + AlphScore. This function also outputs the final MOVA + CADD predicted value.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
Com_CADD_MOVA(protein_name, MOVA_final_predict_file, MOVA_predict_orig_file, MOVA_predict_file, alphscore = 0, phenotype = "Target")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{protein_name}{
       Name of the target protein/gene. Used for the name of the file generation.
%%     ~~Describe \code{protein_name} here~~
}
  \item{MOVA_final_predict_file}{
       The "gene name_Target or Pathogenic_finalpredict.csv" file output by the MOVA function. The final MOVA + CADD predicted values will be in this file.
%%     ~~Describe \code{MOVA_final_predict_file} here~~
}
  \item{MOVA_predict_orig_file}{
       The "gene name_Target or Pathogenic_predict_orig.csv" file output by the MOVA function. Used to redraw the ROC curve for MOVA.
%%     ~~Describe \code{MOVA_predict_orig_file} here~~
}
  \item{MOVA_predict_file}{
       The "gene name_Target or Pathogenic_predict.csv" file output by the MOVA function.
%%     ~~Describe \code{MOVA_predict_file} here~~
}
  \item{alphscore}{
       When alphscore is 1, CADD + alphscore is drawn; when alphscore is 0, it is not drawn. The default is 0.
%%     ~~Describe \code{alphscore} here~~
}
  \item{phenotype}{
             Specify "Target" if you want the positive variant to be the Target variant only, or "Pathogenic" if you want to include the Pathogenic variant as well. The default is "Target".
%%     ~~Describe \code{phenotype} here~~
}
}
\details{
       Draw ROC curve for MOVA + CADD (red), MOVA (blue), CADD (black), (and optional: CADD + AlphScore (gray)). MOVA, CADD, CADD + AlphScore use predetermined values. So, it is necessary to perform the MOVA function in advance. The predicted value of MOVA + CADD is determined by this function. The final MOVA + CADD predicted value are added to the "CADDMOVA" column in the MOVA_final_predict_file. "protein/gene name_Target or Pathogenic_result_CADDMOVA.csv" contains the Cutoff value (Youden index) for each fold (Column: Cutoff), the number of positive variants for each gene used in the analysis (Column: positive_variant_num), number of negative variants (Column: negative_variant_num), AUC for each fold of MOVA + CADD (Column: AUC), cvAUC for MOVA (Column: MOVA), AUC for CADD (Column: CADD), AUC for CADD + AlphScore (Column: com_Alph_CADD), and cvAUC for MOVA + CADD (Column: com_MOVA_CADD). The file required for redrawing with the redraw function of MOVA + CADD is output in "protein/gene name_Target or Pathogenic_predict_orig_CADDMOVA.csv".


%%  ~~ If necessary, more details than the description above ~~
}
\value{
 A data.frame containing the Cutoff value (Youden index) for each fold (Column: Cutoff), the number of positive variants for each gene used in the analysis (Column: positive_variant_num), number of negative variants (Column: negative_variant_num), AUC for each fold of MOVA + CADD (Column: AUC), cvAUC for MOVA (Column: MOVA), AUC for CADD (Column: CADD), AUC for CADD + AlphScore (Column: com_Alph_CADD), and cvAUC for MOVA + CADD (Column: com_MOVA_CADD).
The same data is output to "protein/gene name_Target or Pathogenic_result_CADDMOVA.csv".
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
Hgmd_divide("./source/OPTN_hgmd.csv")
Edit_gnomAD_file("./source/gnomAD_v3.1.2_ENST00000378748_2023_02_15_13_19_36.csv", "./source/OPTN_gnomAD.csv")
Edit_variant_data("../CADD_REVEL/AlphScore_final.tsv", "./source/OPTN_hgmd.csv", "./source/OPTN_gnomAD.csv", "Q96CV9", "OPTN", "./source/OPTNvariantdata.csv")
Edit_polyphen_data("../dbNSFP/dbNSFP4.3a/dbNSFP4.3a_variant.chr10","./source/OPTNvariantdata.csv", "OPTN_alph.csv", "Q96CV9",13100173,13138291,  "./source/OPTNvariantdatapol.csv", "./source/OPTN_alphpol.csv")
Edit_final_variant_file("./source/OPTN_alphpol.csv","./source/OPTN_alphpol2.csv")
MOVA("./source/Q96CV9.fa", "OPTN", "./source/AF-Q96CV9-F1-model_v2.pdb","./source/OPTN_alphpol2.csv", "./source/OPTNvariantdatapol.csv")
Com_CADD_MOVA("OPTN", "OPTN_Target_finalpredict.csv", "OPTN_Target_predict_orig.csv", "OPTN_Target_predict.csv", alphscore = 1)

## The function is currently defined as
function (protein_name, MOVA_final_predict_file, MOVA_predict_orig_file, 
    MOVA_predict_file, alphscore = 0, phenotype = "Target") 
{
    final_data <- fread(MOVA_final_predict_file)
    df <- fread(MOVA_predict_orig_file)
    out <- cvAUC(df$predict, df$result, label.ordering = NULL, 
        folds = df$iter3)
    plot(out$perf, col = "blue", avg = "vertical")
    auc_bind <- data.frame(MOVA = c(out$cvAUC))
    df <- fread(MOVA_predict_file)
    df$CADDMOVA_predict <- 0
    df[, c("CADDMOVA_predict")] <- list(NULL)
    if (phenotype == "Target") {
        df <- df[(df$Type == "Target") | (df$Type == "Ctrl"), 
            ]
    }
    if (phenotype == "Target") {
        df$result <- df$Type2
    }
    else {
        df$result <- df$Type3
    }
    df$predict <- df$CADD_raw
    pred <- prediction(df$predict, df$result)
    perf <- performance(pred, "tpr", "fpr")
    plot(perf, col = "black", add = TRUE)
    auc.tmp <- performance(pred, "auc")
    auc <- as.numeric(auc.tmp@y.values)
    auc_bind <- cbind(auc_bind, data.frame(CADD = c(auc)))
    if (alphscore == 1) {
        df$predict <- df$glm_AlphCadd
        pred <- prediction(df$predict, df$result)
        perf <- performance(pred, "tpr", "fpr")
        plot(perf, col = "gray", add = TRUE)
        auc.tmp <- performance(pred, "auc")
        auc <- as.numeric(auc.tmp@y.values)
        auc_bind <- cbind(auc_bind, data.frame(com_Alph_CADD = c(auc)))
    }
    D <- df
    df <- data.frame(matrix(rep(NA, 6), nrow = 1))[numeric(0), 
        ]
    colnames(df) <- c("ID", "result", "CADDMOVA", "change", "iter1", 
        "iter2")
    D$id2 <- 0
    D[, c("id2")] <- list(NULL)
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
            glmals <- glm(result ~ MOVA_predict + CADD_raw, data = train.data, 
                family = binomial)
            f <- data.frame(ID = val$ID, result = val$result, 
                predict = predict(glmals, newdata = val), change = val$change, 
                iter1 = i, iter2 = i2)
            df <- rbind(df, f)
        }
    }
    glmals <- glm(result ~ MOVA_predict + CADD_raw, data = D, 
        family = binomial)
    f <- data.frame(ID = final_data$ID, CADDMOVA = predict(glmals, 
        newdata = final_data))
    options(warn = oldw)
    final_data <- merge(f, final_data)
    fwrite(final_data, MOVA_final_predict_file)
    df$iter3 <- df$iter1 + (df$iter2 * 5)
    out <- cvAUC(df$predict, df$result, label.ordering = NULL, 
        folds = df$iter3)
    plot(out$perf, col = "red", avg = "vertical", add = TRUE)
    cvauc <- out$cvAUC
    auc_bind <- cbind(auc_bind, data.frame(com_MOVA_CADD = c(cvauc)))
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
            0) & (df$iter2 == 1), ])), AUC = c(auc))
        YI <- cbind(YI, auc_bind)
        YI2 <- rbind(YI2, YI)
    }
    fwrite(YI2, paste(protein_name, phenotype, "result_CADDMOVA.csv", 
        sep = "_"))
    df2 <- df \%>\% group_by(ID) \%>\% summarise(CADDMOVA_predict = mean(CADDMOVA))
    D[, c("result")] <- list(NULL)
    fwrite(merge(D, df2), MOVA_predict_file)
    fwrite(df, paste(protein_name, phenotype, "predict_orig_CADDMOVA.csv", 
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
