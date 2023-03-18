\name{CADD_AlphScore}
\alias{CADD_AlphScore}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function to compare MOVA and CADD + AlphScore.
%%  ~~function to do ... ~~
}
\description{
This function plots the ROC curve in CADD + AlphScore for each gene and returns the AUC.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
CADD_AlphScore(file_name, phenotype = "Target", col = "gray")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{file_name}{
       The "gene name_Target or Pathogenic_predict.csv" file output by the MOVA function
%%     ~~Describe \code{file_name} here~~
}
  \item{phenotype}{
       Specify "Target" if you want the positive variant to be the Target variant only, or "Pathogenic" if you want to include the Pathogenic variant as well. The default is "Target".
%%     ~~Describe \code{phenotype} here~~
}
  \item{col}{
       Specifies the color of the ROC curve. Default is gray.
%%     ~~Describe \code{col} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
Return the AUC for each gene in CADD + AlphScore.
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
MOVA_redraw("OPTN_Target_predict_orig.csv",col="blue")
CADD_MOVA("OPTN_Target_predict.csv", col="black")
CADD_AlphScore("OPTN_Target_predict.csv", col="gray")
MOVA_redraw("OPTN_Target_predict_orig_CADDMOVA.csv",col="red", add=TRUE)


## The function is currently defined as
function (file_name, phenotype = "Target", col = "gray") 
{
    df <- fread(file_name)
    if (phenotype == "Target") {
        df <- df[(df$Type == "Target") | (df$Type == "Ctrl"), 
            ]
    }
    df$predict <- df$glm_AlphCadd
    if (phenotype == "Target") {
        df$result <- df$Type2
    }
    else {
        df$result <- df$Type3
    }
    pred <- prediction(df$predict, df$result)
    perf <- performance(pred, "tpr", "fpr")
    plot(perf, col = col, add = TRUE)
    auc.tmp <- performance(pred, "auc")
    auc <- as.numeric(auc.tmp@y.values)
    return(auc)
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
