\name{MOVA}
\alias{MOVA}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Derivation of predictions and ROC curves by MOVA
%%  ~~function to do ... ~~
}
\description{
Conduct MOVA. The training was performed in a random forest, and the evaluation was performed using the Stratified 5-fold cross validation method, which was repeated five more times to reduce the variability from one execution to another. The ROC curve was used to evaluate the results. In addition, final predictions were calculated for possible nonsynonymous variants.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
MOVA(fasta_file_name, protein_name, pdb_file_name, final_predict_input_file, variant_file, phenotype = "Target")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{fasta_file_name}{
       fasta file that records the amino acid sequence information of the target protein. The first line should describe >Protein, and the amino acid sequence should begin on the second line.
%%     ~~Describe \code{fasta_file_name} here~~
}
  \item{protein_name}{
       Name of target protein/gene. Used to name the file to be exported.
%%     ~~Describe \code{protein_name} here~~
}
  \item{pdb_file_name}{
      AlphaFold2 prediction results pdb file for the target protein. We use data from the AlphaFold2 Protein Structure Database. (https://alphafold.ebi.ac.uk/)
%%     ~~Describe \code{pdb_file_name} here~~
}
  \item{final_predict_input_file}{
       File output by the Edit_final_variant_file function
%%     ~~Describe \code{final_predict_input_file} here~~
}
  \item{variant_file}{
       export_file_name1 by the Edit_polyphen_data function
%%     ~~Describe \code{variant_file} here~~
}
  \item{phenotype}{
       Specify "Target" if you want the positive variant to be the Target variant only, or "Pathogenic" if you want to include the Pathogenic variant as well. The default is "Target".
%%     ~~Describe \code{phenotype} here~~
}
}
\details{
       Draw ROC curve for MOVA (red). MOVA is evaluated using the Stratified 5-fold cross validation method, and the model is repeated five more times. The average of the predicted values (5F-MV value: 5-fold MOVA value) is added to the "MOVA_predict" column of the "gene name_Target or Pathogenic_predict.csv" file. The final MOVA predicted value (MOVA value) are added to the "MOVA_predict" column in "protein/gene name_Target or Pathogenic_finalpredict.csv". "protein/gene name_Target or Pathogenic_result.csv" contains the Cutoff value (Youden index) for each fold (Column: Cutoff), the number of positive variants for each gene used in the analysis (Column: positive_variant_num), number of negative variants (Column: negative_variant_num), AUC for each fold of MOVA (Column: AUC), cvAUC for MOVA (Column: cvauc). The file required for redrawing with the MOVA_redraw function of MOVA is output in "protein/gene name_Target or Pathogenic_predict_orig.csv".
%%  ~~ If necessary, more details than the description above ~~
}
\value{
 A data.frame containing the Cutoff value (Youden index) for each fold (Column: Cutoff), the number of positive variants for each gene used in the analysis (Column: positive_variant_num), number of negative variants (Column: negative_variant_num), AUC for each fold of MOVA (Column: AUC), cvAUC for MOVA (Column: cvauc).
The same data is output to "protein/gene name_Target or Pathogenic_result.csv".
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

## The function is currently defined as
function (fasta_file_name, protein_name, pdb_file_name, final_predict_input_file, 
    variant_file, phenotype = "Target") 
{
    data("BLOSUM62")
    amino_code1 <- c("A", "R", "N", "D", "C", "Q", "E", "G", 
        "I", "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V", 
        "H")
    Protein <- read.fasta(file = fasta_file_name, seqtype = "AA")
    Proteinr <- Protein$Protein[1:length(Protein$Protein)]
    Protein3D <- read.pdb(pdb_file_name)
    Protein3Dr <- as.data.frame(Protein3D$atom)
    Protein3Dp <- Protein3Dr \%>\% group_by(resno) \%>\% summarise(x = mean(x), 
        y = mean(y), z = mean(z), b = mean(b))
    D <- fread(variant_file)
    D$ID <- paste(D$Chromosome, D$Position, D$Reference, D$Alternate, 
        sep = ":")
    D <- D[(D$possible != "DMp") | (is.na(D$possible)), ]
    if (phenotype == "Target") {
        D <- D[(D$Type == "Target") | (D$Type == "Ctrl"), ]
    }
    D[, `:=`(possible, NULL)]
    D$Pos <- as.integer(substr(D$change, 2, nchar(D$change) - 
        1))
    D$ss <- substr(D$change, nchar(D$change), nchar(D$change))
    fpredict <- data.frame(matrix(rep(NA, 3), nrow = 1))[numeric(0), 
        ]
    colnames(fpredict) <- c("change", "predict", "n")
    imp.rfa <- data.frame(matrix(rep(NA, 3), nrow = 1))[numeric(0), 
        ]
    colnames(imp.rfa) <- c("name", "IncNodePurity", "n")
    for (i in 1:nrow(D)) {
        D$con[i] <- BLOSUM62[Proteinr[D$Pos[i]], Proteinr[D$Pos[i]]] - 
            BLOSUM62[Proteinr[D$Pos[i]], D$ss[i]]
        D$x[i] <- Protein3Dp$x[D$Pos[i]]
        D$y[i] <- Protein3Dp$y[D$Pos[i]]
        D$z[i] <- Protein3Dp$z[D$Pos[i]]
        D$b[i] <- Protein3Dp$b[D$Pos[i]]
    }
    final_data <- fread(final_predict_input_file)
    final_data$ID <- paste(final_data$Chromosome, final_data$Position, 
        final_data$Reference, final_data$Alternate, final_data$aaref, 
        final_data$Pos, final_data$aaalt, sep = ":")
    for (i in 1:nrow(final_data)) {
        final_data$con[i] <- BLOSUM62[Proteinr[final_data$Pos[i]], 
            Proteinr[final_data$Pos[i]]] - BLOSUM62[Proteinr[final_data$Pos[i]], 
            final_data$aaalt[i]]
        final_data$x[i] <- Protein3Dp$x[final_data$Pos[i]]
        final_data$y[i] <- Protein3Dp$y[final_data$Pos[i]]
        final_data$z[i] <- Protein3Dp$z[final_data$Pos[i]]
        final_data$b[i] <- Protein3Dp$b[final_data$Pos[i]]
    }
    df <- data.frame(matrix(rep(NA, 6), nrow = 1))[numeric(0), 
        ]
    colnames(df) <- c("ID", "result", "predict", "change", "iter1", 
        "iter2")
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
            rf <- randomForest(result ~ con + x + y + z + b, 
                train.data, importance = TRUE)
            imp.rf <- importance(rf)
            f <- data.frame(ID = val$ID, result = val$result, 
                predict = predict(rf, val), change = val$change, 
                iter1 = i, iter2 = i2)
            df <- rbind(df, f)
            f <- data.frame(name = rownames(imp.rf), IncNodePurity = imp.rf, 
                n = i)
            imp.rfa <- rbind(imp.rfa, f)
        }
    }
    for (i in 1:30) {
        rf <- randomForest(result ~ con + x + y + z + b, D)
        f <- data.frame(ID = final_data$ID, predict = predict(rf, 
            final_data), n = i)
        fpredict <- rbind(fpredict, f)
    }
    options(warn = oldw)
    fpredictb <- fpredict \%>\% group_by(ID) \%>\% summarise(MOVA_predict = mean(predict))
    imp.rfb <- imp.rfa \%>\% group_by(name) \%>\% summarise(IncMSE = mean(IncNodePurity..IncMSE), 
        IncNodePurity = mean(IncNodePurity.IncNodePurity))
    ID <- str_split(fpredictb$ID, pattern = ":", simplify = TRUE)
    for (i in 1:nrow(fpredictb)) {
        if (ID[i, 5] == ID[i, 7]) {
            fpredictb$predict[i] = 0
        }
    }
    final_data <- merge(fpredictb, final_data)
    fwrite(final_data, paste(protein_name, phenotype, "finalpredict.csv", 
        sep = "_"))
    fwrite(imp.rfb, paste(protein_name, phenotype, "importance.csv", 
        sep = "_"))
    df$iter3 <- df$iter1 + (df$iter2 * 5)
    out <- cvAUC(df$predict, df$result, label.ordering = NULL, 
        folds = df$iter3)
    plot(out$perf, col = "red", avg = "vertical")
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
    fwrite(YI2, paste(protein_name, phenotype, "result.csv", 
        sep = "_"))
    df2 <- df \%>\% group_by(ID) \%>\% summarise(MOVA_predict = mean(predict))
    D[, c("result")] <- list(NULL)
    fwrite(merge(D, df2), paste(protein_name, phenotype, "predict.csv", 
        sep = "_"))
    fwrite(df, paste(protein_name, phenotype, "predict_orig.csv", 
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
