\name{variant_predict_plot}
\alias{variant_predict_plot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
The predicted probability of pathogenicity was compared between MOVA and PolyPhen-2.
%%  ~~function to do ... ~~
}
\description{
Create a graph corresponding to Figure 3 in the paper.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
variant_predict_plot(file_name, legend_type = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{file_name}{
       The "gene name_Target or Pathogenic_predict.csv" file output by the MOVA function
%%     ~~Describe \code{file_name} here~~
}
  \item{legend_type}{
             Specify "TRUE" if you want the legend, or "FALSE" if you don't want the legend. The default is "FALSE".
%%     ~~Describe \code{legend_type} here~~
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
variant_predict_plot("TARDBP_Target_predict.csv")
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (file_name, legend_type = FALSE) 
{
    df <- fread(file_name)
    df$Typevariant <- "negative"
    df[df$Type == "Target", `:=`(Typevariant, "positive")]
    f <- data.table(variant = df$Typevariant, predict = df$MOVA_predict, 
        position = df$Pos, type = "MOVA")
    f2 <- data.table(variant = df$Typevariant, predict = df$pph_prob, 
        position = df$Pos, type = "Polyphen-2")
    df2 <- rbind(f, f2)
    df2 <- transform(df2, variant = factor(variant, levels = c("positive", 
        "negative")))
    g <- ggplot(data = df2)
    g <- g + theme(axis.text = element_text(size = 12))
    g <- g + theme(legend.text = element_text(size = 12))
    g <- g + geom_point(aes(x = position, y = predict, colour = variant))
    g <- g + facet_grid(rows = vars(type), scales = "free_y")
    if (legend_type == TRUE) {
        g <- g + theme(legend.position = "bottom")
    }
    else {
        g <- g + theme(legend.position = "none")
    }
    g
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
