\name{MOVA_only_3d_distance}
\alias{MOVA_only_3d_distance}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
A simple baseline approach based on distance to the known pathogenic variants
%%  ~~function to do ... ~~
}
\description{
This function evaluates whether pathogenicity can be predicted by the distance between the variant of interest and a known pathogenic variant using a Stratified 5-fold cross validation method.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
MOVA_only_3d_distance(protein_name, MOVA_final_predict_file, MOVA_predict_file, phenotype = "Target")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{protein_name}{
       Name of target protein/gene. Used to name the file to be exported.
%%     ~~Describe \code{protein_name} here~~
}
  \item{MOVA_final_predict_file}{
       The "gene name_Target or Pathogenic_finalpredict.csv" file output by the MOVA function. The final MOVA_only_3d_distance predicted values will be in this file.
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
        Draw ROC curve for MOVA_only_3d_distance (black).  MOVA_only_3d_distance is evaluated using the Stratified 5-fold cross validation method, and the model is repeated five more times. The average of the distance to the known pathogenic variants is added to the "t_distance" column of the MOVA_predict_file. "protein/gene name_Target or Pathogenic_result_only_3d_distance.csv" contains the Cutoff value (Youden index) for each fold (Column: Cutoff), the number of positive variants for each gene used in the analysis (Column: positive_variant_num), number of negative variants (Column: negative_variant_num), AUC for each fold of MOVA_only_3d_distance (Column: AUC), cvAUC for MOVA_only_3d_distance (Column: cvauc). The file required for redrawing with the MOVA_redraw function of MOVA_only_3d_distance is output in "protein/gene name_Target or Pathogenic_predict_orig_only_3d_distance.csv".
%%  ~~ If necessary, more details than the description above ~~
}
\value{
 A data.frame containing the Cutoff value (Youden index) for each fold (Column: Cutoff), the number of positive variants for each gene used in the analysis (Column: positive_variant_num), number of negative variants (Column: negative_variant_num), AUC for each fold of MOVA_only_3d_distance (Column: AUC), cvAUC for MOVA_only_3d_distance (Column: cvauc).
The same data is output to "protein/gene name_Target or Pathogenic_result_only_3d_distance.csv".
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
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (protein_name, MOVA_final_predict_file, MOVA_predict_file, 
    phenotype = "Target") 
{
    D <- fread(MOVA_predict_file)
    D$t_distance <- 0
    D[, c("t_distance")] <- list(NULL)
    if (phenotype == "Target") {
        D <- D[(D$Type == "Target") | (D$Type == "Ctrl"), ]
    }
    fpredict <- data.frame(matrix(rep(NA, 2), nrow = 1))[numeric(0), 
        ]
    colnames(fpredict) <- c("ID", "distance")
    final_data <- fread(MOVA_final_predict_file)
    df <- data.frame(matrix(rep(NA, 6), nrow = 1))[numeric(0), 
        ]
    colnames(df) <- c("ID", "re_result", "distance", "change", 
        "iter1", "iter2")
    D$id2 <- 0
    D[, c("id2")] <- list(NULL)
    D$re_result <- 0
    D[, c("re_result")] <- list(NULL)
    D$t_distance <- 0
    D[, c("t_distance")] <- list(NULL)
    D$iter1 <- 0
    D[, c("iter1")] <- list(NULL)
    D$iter2 <- 0
    D[, c("iter2")] <- list(NULL)
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
            for (i3 in 1:nrow(val)) {
                f <- data.frame(ID = val[i3]$ID, re_result = abs(val[i3]$result - 
                  1), t_distance = min(sqrt((x - val[i3]$x) * 
                  (x - val[i3]$x) + (y - val[i3]$y) * (y - val[i3]$y) + 
                  (z - val[i3]$z) * (z - val[i3]$z))), change = val[i3]$change, 
                  iter1 = i, iter2 = i2)
                df <- rbind(df, f)
            }
        }
    }
    options(warn = oldw)
    df$iter3 <- df$iter1 + (df$iter2 * 5)
    out <- cvAUC(df$t_distance, df$re_result, label.ordering = NULL, 
        folds = df$iter3)
    plot(out$perf, col = "black", avg = "vertical", add = TRUE)
    cvauc <- out$cvAUC
    YI2 <- data.frame(matrix(rep(NA, 5), nrow = 1))[numeric(0), 
        ]
    colnames(YI2) <- c("Cutoff", "positive_variant_num", "negative_variant_num", 
        "AUC", "cvauc")
    for (i in 6:30) {
        pred <- prediction(df[df$iter3 == i, ]$t_distance, df[df$iter3 == 
            i, ]$re_result)
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
        YI <- data.frame(Cutoff = YI$Cutoff[1], positive_variant_num = c(nrow(df[(df$re_result == 
            0) & (df$iter2 == 1), ])), negative_variant_num = c(nrow(df[(df$re_result == 
            1) & (df$iter2 == 1), ])), AUC = c(auc), cvauc = c(cvauc))
        YI2 <- rbind(YI2, YI)
    }
    fwrite(YI2, paste(protein_name, phenotype, "result_only_3d_distance.csv", 
        sep = "_"))
    df2 <- df \%>\% group_by(ID) \%>\% summarise(t_distance = mean(t_distance))
    D[, c("result")] <- list(NULL)
    fwrite(merge(D, df2), MOVA_predict_file)
    fwrite(df, paste(protein_name, phenotype, "predict_orig_only_3d_distance.csv", 
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
