\name{REVEL_AlphScore}
\alias{REVEL_AlphScore}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function to compare MOVA and REVEL + AlphScore.
%%  ~~function to do ... ~~
}
\description{
This function plots the ROC curve in REVEL + AlphScore for each gene and returns the AUC.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
REVEL_AlphScore(file_name, phenotype = "Target", col = "gray")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{file_name}{
       The "gene name_Target or Pathogenic_predict.csv" file output by the MOVA function.
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
Return the AUC for each gene in REVEL + AlphScore.
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
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (file_name, phenotype = "Target", col = "gray") 
{
    df <- fread(file_name)
    if (phenotype == "Target") {
        df <- df[(df$Type == "Target") | (df$Type == "Ctrl"), 
            ]
    }
    df$predict <- df$glm_AlphRevel
    df <- df[df$predict != ".", ]
    df$predict <- as.numeric(df$predict)
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
