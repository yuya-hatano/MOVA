\name{PolyPhen_MOVA}
\alias{PolyPhen_MOVA}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function to compare MOVA and PolyPhen-2.
%%  ~~function to do ... ~~
}
\description{
This function plots the ROC curve in PolyPhen-2 for each gene and returns the AUC.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
PolyPhen_MOVA(file_name, phenotype = "Target", col = "orange")
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
       Specifies the color of the ROC curve. Default is orange.
%%     ~~Describe \code{col} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
Return the AUC for each gene in PolyPhen-2.
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
PolyPhen_MOVA("TARDBP_Target_predict.csv")
AlphScore_MOVA("TARDBP_Target_predict.csv")
CADD_MOVA("TARDBP_Target_predict.csv")
REVEL_MOVA("TARDBP_Target_predict.csv")
EVE_MOVA("TARDBP_Target_predict.csv", "../dementia/EVE_all_data/variant_files/TADBP_HUMAN.csv")
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (file_name, phenotype = "Target", col = "orange") 
{
    df <- fread(file_name)
    if (phenotype == "Target") {
        df <- df[(df$Type == "Target") | (df$Type == "Ctrl"), 
            ]
    }
    df$predict <- df$pph_prob
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
