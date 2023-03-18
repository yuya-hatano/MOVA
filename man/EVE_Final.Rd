\name{EVE_Final}
\alias{EVE_Final}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
 Add the predicted value of EVE for each variant.
%%  ~~function to do ... ~~
}
\description{
   Add the predicted value of EVE for each variant in "gene name_Target or Pathogenic_finalpredict.csv".
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
EVE_Final(file_name, EVE_file_name)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{file_name}{
       Path to "gene name_Target or Pathogenic_finalpredict.csv" file output by the MOVA function
%%     ~~Describe \code{file_name} here~~
}
  \item{EVE_file_name}{
       EVE variant files for targeted genes. These files are available at https://evemodel.org/.
%%     ~~Describe \code{EVE_file_name} here~~
}
}
\details{
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
EVE_Final("TARDBP_Target_finalpredict.csv", "../dementia/EVE_all_data/variant_files/TADBP_HUMAN.csv")

## The function is currently defined as
function (file_name, EVE_file_name) 
{
    df <- fread(file_name)
    y <- fread(EVE_file_name)
    df[, `:=`(EVE, NA)]
    df[, `:=`(EVE_class, NA)]
    df$Pos <- as.integer(substr(df$change, 2, nchar(df$change) - 
        1))
    for (i in 1:nrow(df)) {
        df$EVE[i] <- y[(y$position == df$Pos[i]) & (y$mt_aa == 
            df$aaalt[i])]$EVE_scores_ASM
        df$EVE_class[i] <- y[(y$position == df$Pos[i]) & (y$mt_aa == 
            df$aaalt[i])]$EVE_classes_75_pct_retained_ASM
    }
    fwrite(df, file_name)
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
