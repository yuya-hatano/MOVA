\name{MOVA_3d_distance}
\alias{MOVA_3d_distance}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  MOVA with 'distance to the known pathogenic variant' learned instead of position coordinates
%%  ~~function to do ... ~~
}
\description{
MOVA with 'distance to the known pathogenic variant' learned instead of position coordinates.
MOVA_predistance function must be performed before this function.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
MOVA_3d_distance(protein_name, MOVA_final_predict_file, MOVA_predict_file, phenotype = "Target")
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
        Draw ROC curve for MOVA_3d_distance (blue). MOVA_3d_distance is the MOVA features minus the location information and plus 'distance to the known pathogenic variant'. MOVA_3d_distance is evaluated using the Stratified 5-fold cross validation method, and the model is repeated five more times. The average of the predicted values is added to the "MOVA_3d_distance_predict" column of the MOVA_predict_file. The final MOVA_3d_distance predicted value are added to the "MOVA_3d_distance_predict" column in MOVA_final_predict_file. "protein/gene name_Target or Pathogenic_result_3d_distance.csv" contains the Cutoff value (Youden index) for each fold (Column: Cutoff), the number of positive variants for each gene used in the analysis (Column: positive_variant_num), number of negative variants (Column: negative_variant_num), AUC for each fold of MOVA_3d_distance (Column: AUC), cvAUC for MOVA_3d_distance (Column: cvauc). The file required for redrawing with the MOVA_redraw function of MOVA_3d_distance is output in "protein/gene name_Target or Pathogenic_predict_orig_3d_distance.csv".
%%  ~~ If necessary, more details than the description above ~~
}
\value{
 A data.frame containing the Cutoff value (Youden index) for each fold (Column: Cutoff), the number of positive variants for each gene used in the analysis (Column: positive_variant_num), number of negative variants (Column: negative_variant_num), AUC for each fold of MOVA_3d_distance (Column: AUC), cvAUC for MOVA_3d_distance (Column: cvauc).
The same data is output to "protein/gene name_Target or Pathogenic_result_3d_distance.csv".
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
MOVA_predistance("TARDBP", "TARDBP_Target_finalpredict.csv","TARDBP_Target_predict.csv")
MOVA_3d_distance("TARDBP", "TARDBP_Target_finalpredict.csv","TARDBP_Target_predict.csv")

## The function is currently defined as
function (protein_name, MOVA_final_predict_file, MOVA_predict_file, 
    phenotype = "Target") 
{
    D <- fread(MOVA_predict_file)
    D$MOVA_3d_distance_predict <- 0
    D[, c("MOVA_3d_distance_predict")] <- list(NULL)
    if (phenotype == "Target") {
        D <- D[(D$Type == "Target") | (D$Type == "Ctrl"), ]
    }
    fpredict <- data.frame(matrix(rep(NA, 3), nrow = 1))[numeric(0), 
        ]
    colnames(fpredict) <- c("change", "predict", "n")
    imp.rfa <- data.frame(matrix(rep(NA, 3), nrow = 1))[numeric(0), 
        ]
    colnames(imp.rfa) <- c("name", "IncNodePurity", "n")
    final_data <- fread(MOVA_final_predict_file)
    final_data$MOVA_3d_distance_predict <- 0
    final_data[, c("MOVA_3d_distance_predict")] <- list(NULL)
    df <- data.frame(matrix(rep(NA, 6), nrow = 1))[numeric(0), 
        ]
    colnames(df) <- c("ID", "re_result", "distance", "change", 
        "iter1", "iter2")
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
            x <- train.data[train.data$result == 1, ]$x
            y <- train.data[train.data$result == 1, ]$y
            z <- train.data[train.data$result == 1, ]$z
            train.data$t_distance <- 0
            for (i3 in 1:nrow(train.data)) {
                t_distance <- sqrt((x - train.data[i3]$x) * (x - 
                  train.data[i3]$x) + (y - train.data[i3]$y) * 
                  (y - train.data[i3]$y) + (z - train.data[i3]$z) * 
                  (z - train.data[i3]$z))
                t_distance <- sort(t_distance, decreasing = F)
                if (train.data[i3]$result == 1) {
                  train.data[i3]$t_distance <- t_distance[2]
                }
                else {
                  train.data[i3]$t_distance <- t_distance[1]
                }
            }
            rf <- randomForest(result ~ con + t_distance + b, 
                train.data, importance = TRUE)
            imp.rf <- importance(rf)
            x <- train.data[train.data$result == 1, ]$x
            y <- train.data[train.data$result == 1, ]$y
            z <- train.data[train.data$result == 1, ]$z
            val$t_distance <- 0
            for (i3 in 1:nrow(val)) {
                val[i3]$t_distance <- min(sqrt((x - val[i3]$x) * 
                  (x - val[i3]$x) + (y - val[i3]$y) * (y - val[i3]$y) + 
                  (z - val[i3]$z) * (z - val[i3]$z)))
            }
            f <- data.frame(ID = val$ID, result = val$result, 
                predict = predict(rf, val), change = val$change, 
                iter1 = i, iter2 = i2)
            df <- rbind(df, f)
        }
    }
    f <- data.frame(name = rownames(imp.rf), IncNodePurity = imp.rf, 
        n = i)
    imp.rfa <- rbind(imp.rfa, f)
    x <- D[D$result == 1, ]$x
    y <- D[D$result == 1, ]$y
    z <- D[D$result == 1, ]$z
    D$t_distance <- 0
    for (i3 in 1:nrow(D)) {
        t_distance <- sqrt((x - D[i3]$x) * (x - D[i3]$x) + (y - 
            D[i3]$y) * (y - D[i3]$y) + (z - D[i3]$z) * (z - D[i3]$z))
        t_distance <- sort(t_distance, decreasing = F)
        if (D[i3]$result == 1) {
            D[i3]$t_distance <- t_distance[2]
        }
        else {
            D[i3]$t_distance <- t_distance[1]
        }
    }
    for (i in 1:30) {
        rf <- randomForest(result ~ con + t_distance + b, D)
        f <- data.frame(ID = final_data$ID, predict = predict(rf, 
            final_data), n = i)
        fpredict <- rbind(fpredict, f)
    }
    options(warn = oldw)
    fpredictb <- fpredict \%>\% group_by(ID) \%>\% summarise(MOVA_3d_distance_predict = mean(predict))
    imp.rfb <- imp.rfa \%>\% group_by(name) \%>\% summarise(IncMSE = mean(IncNodePurity..IncMSE), 
        IncNodePurity = mean(IncNodePurity.IncNodePurity))
    final_data <- merge(fpredictb, final_data)
    fwrite(final_data, MOVA_final_predict_file)
    fwrite(imp.rfb, paste(protein_name, phenotype, "importance_3d_distance.csv", 
        sep = "_"))
    df$iter3 <- df$iter1 + (df$iter2 * 5)
    out <- cvAUC(df$predict, df$result, label.ordering = NULL, 
        folds = df$iter3)
    plot(out$perf, col = "blue", avg = "vertical", add = TRUE)
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
    fwrite(YI2, paste(protein_name, phenotype, "result_3d_distance.csv", 
        sep = "_"))
    df2 <- df \%>\% group_by(ID) \%>\% summarise(MOVA_3d_distance_predict = mean(predict))
    D[, c("result")] <- list(NULL)
    fwrite(merge(D, df2), MOVA_predict_file)
    fwrite(df, paste(protein_name, phenotype, "predict_orig_3d_distance.csv", 
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
