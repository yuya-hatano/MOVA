\name{MOVA_predistance}
\alias{MOVA_predistance}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
 Calculate the 'distance to the known pathogenic variant' for each variant in MOVA_final_predict_file.
%%  ~~function to do ... ~~
}
\description{
  Calculate the 'distance to the known pathogenic variant' for each variant in MOVA_final_predict_file. This function should be done in advance of the machine learning function with 'distance to the known pathogenic variant' as a feature.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
MOVA_predistance(protein_name, MOVA_final_predict_file, MOVA_predict_file, phenotype = "Target")
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
       This function must be performed before the following functions are performed.
       *MOVA_3d_distance
       *MOVA_3d_distance_log
       *MOVA_LOPOV
       *MOVA_plus_3d_distance
       *MOVA_plus_3d_distance_SVM
       *MOVA_plus_3d_distance_xgboost
%%  ~~ If necessary, more details than the description above ~~
}
\value{
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
MOVA_plus_3d_distance("TARDBP", "TARDBP_Target_finalpredict.csv","TARDBP_Target_predict.csv")

## The function is currently defined as
function (protein_name, MOVA_final_predict_file, MOVA_predict_file, 
    phenotype = "Target") 
{
    D <- fread(MOVA_predict_file)
    if (phenotype == "Target") {
        D <- D[(D$Type == "Target") | (D$Type == "Ctrl"), ]
    }
    final_data <- fread(MOVA_final_predict_file)
    if (phenotype == "Target") {
        x <- D[Type2 == 1, ]$x
        y <- D[Type2 == 1, ]$y
        z <- D[Type2 == 1, ]$z
        D$t_distance <- 0
        final_data$t_distance <- 0
        for (i3 in 1:nrow(final_data)) {
            t_distance <- sqrt((x - final_data[i3]$x) * (x - 
                final_data[i3]$x) + (y - final_data[i3]$y) * 
                (y - final_data[i3]$y) + (z - final_data[i3]$z) * 
                (z - final_data[i3]$z))
            t_distance <- sort(t_distance, decreasing = F)
            final_data[i3]$t_distance <- t_distance[1]
        }
    }
    else {
        x <- D[Type3 == 1, ]$x
        y <- D[Type3 == 1, ]$y
        z <- D[Type3 == 1, ]$z
        D$t_distance <- 0
        final_data$t_distance <- 0
        for (i3 in 1:nrow(final_data)) {
            t_distance <- sqrt((x - final_data[i3]$x) * (x - 
                final_data[i3]$x) + (y - final_data[i3]$y) * 
                (y - final_data[i3]$y) + (z - final_data[i3]$z) * 
                (z - final_data[i3]$z))
            t_distance <- sort(t_distance, decreasing = F)
            final_data[i3]$t_distance <- t_distance[1]
        }
    }
    fwrite(final_data, MOVA_final_predict_file)
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
