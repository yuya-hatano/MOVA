\name{MOVA_redraw}
\alias{MOVA_redraw}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Redraw the resulting ROC curve once analyzed.
%%  ~~function to do ... ~~
}
\description{
Redraw the resulting ROC curve once analyzed.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
MOVA_redraw(MOVA_predict_orig_file, col, add = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{MOVA_predict_orig_file}{
              "protein/gene name_Target or Pathogenic_predict_orig.csv" output by the MOVA function, or csv file for MOVA_redraw output by each function
%%     ~~Describe \code{MOVA_predict_orig_file} here~~
}
  \item{col}{
       Color of ROC curve
%%     ~~Describe \code{col} here~~
}
  \item{add}{
             When TRUE, it overlays an existing graph; when FALSE, it creates a new graph.
%%     ~~Describe \code{add} here~~
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
MOVA_SVM("TARDBP", "TARDBP_Target_finalpredict.csv","TARDBP_Target_predict.csv")
MOVA_xgboost("TARDBP", "TARDBP_Target_finalpredict.csv","TARDBP_Target_predict.csv")
MOVA_redraw("TARDBP_Target_predict_orig.csv",col="red")
PolyPhen_MOVA("TARDBP_Target_predict.csv")
AlphScore_MOVA("TARDBP_Target_predict.csv")
CADD_MOVA("TARDBP_Target_predict.csv")
REVEL_MOVA("TARDBP_Target_predict.csv")
EVE_MOVA("TARDBP_Target_predict.csv", "../dementia/EVE_all_data/variant_files/TADBP_HUMAN.csv")
MOVA_redraw("TARDBP_Target_predict_orig.csv",col="red")
MOVA_redraw("TARDBP_Target_predict_orig_SVM.csv",col="black", add=TRUE)
MOVA_redraw("TARDBP_Target_predict_orig_xgboost.csv",col="blue", add=TRUE)

##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (MOVA_predict_orig_file, col, add = FALSE) 
{
    df <- fread(MOVA_predict_orig_file)
    out <- cvAUC(df$predict, df$result, label.ordering = NULL, 
        folds = df$iter3)
    plot(out$perf, col = col, avg = "vertical", add = add)
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
