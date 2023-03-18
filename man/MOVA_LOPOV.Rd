\name{MOVA_LOPOV}
\alias{MOVA_LOPOV}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  do a leave one protein out validation (using 11 proteins' variations for training, 1 protein's for testing)
%%  ~~function to do ... ~~
}
\description{
　The Random forest analysis with 'distance to the known pathogenic variant', ΔBLOSUM62,pLDDT as explanatory variables and the pathogenicity of the variant as objective variable.
  The study was trained with 11 genes, one out of 12 ALS-related genes (TARDBP, FUS, SETX, TBK1, OPTN, SOD1, VCP, SQSTM1, ANG, UBQLN2, DCTN1, CCNF), and validated with the remaining one gene. The process was repeated with different validation genes to ensure that all genes were validated.
MOVA_predistance function must be performed before this function.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
MOVA_LOPOV(phenotype = "Target")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{phenotype}{
      Specify "Target" if you want the positive variant to be the Target variant only, or "Pathogenic" if you want to include the Pathogenic variant as well. The default is "Target".
%%     ~~Describe \code{phenotype} here~~
}
}
\details{
       To perform this function, MOVA fuction must be performed on all genes in TARDBP, FUS, SETX, TBK1, OPTN, SOD1, VCP, SQSTM1, ANG, UBQLN2, DCTN1, CCNF. MOVA_LOPOV is a random forest analysis with 'distance to the known pathogenic variant', ΔBLOSUM62, pLDDT as explanatory variables and the pathogenicity of the variant as objective variable. MOVA_LOPOV is evaluated using the leave one protein out validation. The MOVA_LOPOV predicted value are added to the "MOVA_LOPOV_predict" column in "LOPOV_Target or Pathogenic_finalpredict.csv". "LOPOV_Target or Pathogenic_result.csv" contains the Cutoff value (Youden index) for each gene (Column: Cutoff), the number of positive variants for each gene used in the analysis (Column: positive_variant_num), number of negative variants (Column: negative_variant_num), AUC for each gene of MOVA_LOPOV (Column: AUC), cvAUC for MOVA_LOPOV (Column: cvauc). The file required for redrawing with the redraw function of MOVA_LOPOV is output in "LOPOV_Target or Pathogenic_predict_orig.csv".
%%  ~~ If necessary, more details than the description above ~~
}
\value{
 A data.frame containing the Cutoff value (Youden index) for each gene (Column: Cutoff), the number of positive variants for each gene used in the analysis (Column: positive_variant_num), number of negative variants (Column: negative_variant_num), AUC for each gene of MOVA_LOPOV (Column: AUC), cvAUC for MOVA_LOPOV (Column: cvauc).
The same data is output to "LOPOV_Target or Pathogenic_result.csv".
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
MOVA_LOPOV()

## The function is currently defined as
function (phenotype = "Target") 
{
    read_name <- c("TARDBP", "FUS", "SETX", "TBK1", "OPTN", "SOD1", 
        "VCP", "SQSTM1", "ANG", "UBQLN2", "DCTN1", "CCNF")
    D <- fread(paste(read_name[1], phenotype, "predict.csv", 
        sep = "_"))
    final_data <- fread(paste(read_name[1], phenotype, "finalpredict.csv", 
        sep = "_"))
    D$gene <- read_name[1]
    D$gene_number <- 1
    D$glm_CaddDeogenRevel <- 0
    D[, `:=`(glm_CaddDeogenRevel, NULL)]
    final_data$gene <- read_name[1]
    final_data$gene_number <- 1
    final_data$glm_CaddDeogenRevel <- 0
    final_data[, `:=`(glm_CaddDeogenRevel, NULL)]
    D$t_distance <- 0
    if (phenotype == "Target") {
        D$result <- D$Type2
    }
    else {
        D$result <- D$Type3
    }
    x <- D[D$result == 1, ]$x
    y <- D[D$result == 1, ]$y
    z <- D[D$result == 1, ]$z
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
    for (i in 2:12) {
        df2 <- fread(paste(read_name[i], phenotype, "predict.csv", 
            sep = "_"))
        final_data2 <- fread(paste(read_name[i], phenotype, "finalpredict.csv", 
            sep = "_"))
        df2$gene <- read_name[i]
        df2$gene_number <- i
        df2$glm_CaddDeogenRevel <- 0
        df2[, `:=`(glm_CaddDeogenRevel, NULL)]
        if (phenotype == "Target") {
            df2$result <- df2$Type2
        }
        else {
            df2$result <- df2$Type3
        }
        df2$t_distance <- 0
        x <- df2[df2$result == 1, ]$x
        y <- df2[df2$result == 1, ]$y
        z <- df2[df2$result == 1, ]$z
        for (i3 in 1:nrow(df2)) {
            t_distance <- sqrt((x - df2[i3]$x) * (x - df2[i3]$x) + 
                (y - df2[i3]$y) * (y - df2[i3]$y) + (z - df2[i3]$z) * 
                (z - df2[i3]$z))
            t_distance <- sort(t_distance, decreasing = F)
            if (df2[i3]$result == 1) {
                df2[i3]$t_distance <- t_distance[2]
            }
            else {
                df2[i3]$t_distance <- t_distance[1]
            }
        }
        D <- rbind(D, df2)
        final_data2$gene <- read_name[i]
        final_data2$gene_number <- i
        final_data2$glm_CaddDeogenRevel <- 0
        final_data2[, `:=`(glm_CaddDeogenRevel, NULL)]
        final_data <- rbind(final_data, final_data2)
    }
    df <- data.frame(matrix(rep(NA, 6), nrow = 1))[numeric(0), 
        ]
    colnames(df) <- c("ID", "result", "predict", "change", "gene_number", 
        "gene")
    fpredict <- data.frame(matrix(rep(NA, 3), nrow = 1))[numeric(0), 
        ]
    colnames(fpredict) <- c("change", "predict", "n")
    imp.rfa <- data.frame(matrix(rep(NA, 3), nrow = 1))[numeric(0), 
        ]
    colnames(imp.rfa) <- c("name", "IncNodePurity", "n")
    oldw <- getOption("warn")
    options(warn = -1)
    for (i in 1:12) {
        val <- D[D$gene_number == i, ]
        train.data <- D[D$gene_number != i, ]
        rf <- randomForest(result ~ con + t_distance + b, train.data, 
            importance = TRUE)
        imp.rf <- importance(rf)
        f <- data.frame(ID = val$ID, result = val$result, predict = predict(rf, 
            val), change = val$change, gene_number = i, gene = val$gene)
        df <- rbind(df, f)
        f <- data.frame(ID = final_data[final_data$gene_number == 
            i, ]$ID, MOVA_LOPOV_predict = predict(rf, final_data[final_data$gene_number == 
            i, ]), gene_number = i)
        fpredict <- rbind(fpredict, f)
        f <- data.frame(name = rownames(imp.rf), IncNodePurity = imp.rf, 
            n = i)
        imp.rfa <- rbind(imp.rfa, f)
    }
    options(warn = oldw)
    fpredictb <- fpredict
    imp.rfb <- imp.rfa \%>\% group_by(name) \%>\% summarise(IncMSE = mean(IncNodePurity..IncMSE), 
        IncNodePurity = mean(IncNodePurity.IncNodePurity))
    final_data <- merge(fpredictb, final_data)
    fwrite(final_data, paste("LOPOV", phenotype, "finalpredict.csv", 
        sep = "_"))
    fwrite(imp.rfb, paste("LOPOV", phenotype, "importance.csv", 
        sep = "_"))
    out <- cvAUC(df$predict, df$result, label.ordering = NULL, 
        folds = df$gene_number)
    plot(out$perf, col = "red", avg = "vertical")
    cvauc <- out$cvAUC
    YI2 <- data.frame(matrix(rep(NA, 5), nrow = 1))[numeric(0), 
        ]
    colnames(YI2) <- c("Cutoff", "positive_variant_num", "negative_variant_num", 
        "AUC", "cvauc")
    for (i in 1:12) {
        pred <- prediction(df[df$gene_number == i, ]$predict, 
            df[df$gene_number == i, ]$result)
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
            1) & (df$gene_number == 1), ])), negative_variant_num = c(nrow(df[(df$result == 
            0) & (df$gene_number == 1), ])), AUC = c(auc), cvauc = c(cvauc), 
            gene = c(read_name[i]))
        YI2 <- rbind(YI2, YI)
    }
    fwrite(YI2, paste("LOPOV", phenotype, "result.csv", sep = "_"))
    df2 <- df \%>\% group_by(ID) \%>\% summarise(MOVA_LOPOV_predict = mean(predict))
    fwrite(merge(D, df2), paste("LOPOV", phenotype, "predict.csv", 
        sep = "_"))
    fwrite(df, paste("LOPOV", phenotype, "predict_orig.csv", 
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
