\name{Hgmd_divide}
\alias{Hgmd_divide}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Creation of Hgmd files
%%  ~~function to do ... ~~
}
\description{
Create a positive variant file from the Hgmd data.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
Hgmd_divide(hgmd_file_name)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{hgmd_file_name}{
      Path to hgmd file name
%%     ~~Describe \code{hgmd_file_name} here~~
}
}
\details{
    Based on HGMD data, list amino acid changes (column: HGVSprotein), VCF (column: VCF), DM if DM in HGMD, DMp  if DM? (column: possible), Pathogenic if pathogenic but not a target phenotype, and Target if target phenotype(column: Type) for each variant and save it. This is processed by Hgmd_divide and can be processed by Edit_variant_data.
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

## The function is currently defined as
function (hgmd_file_name) 
{
    x <- fread(hgmd_file_name)
    x <- x \%>\% separate(VCF, c("buf1", "buf2"), sep = ":")
    x <- x \%>\% separate(buf2, c("buf3", "Alternate"), sep = "/")
    x$Chromosome <- str_sub(x$buf1, 4, -1)
    x$Position <- str_sub(x$buf3, 1, -2)
    x$Reference <- str_sub(x$buf3, -1, -1)
    x[, `:=`(buf1, NULL)]
    x[, `:=`(buf3, NULL)]
    fwrite(x, hgmd_file_name)
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
