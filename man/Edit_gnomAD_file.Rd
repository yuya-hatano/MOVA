\name{Edit_gnomAD_file}
\alias{Edit_gnomAD_file}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function to retrieve the required information from the downloaded gnomAD file
%%  ~~function to do ... ~~
}
\description{
Function to retrieve the required information from the downloaded gnomAD file.
Before applying this function, be sure to check the contents of the file yourself and delete all but the missense mutations.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
Edit_gnomAD_file(input_file_name, export_file_name)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{input_file_name}{
       The name of downloaded gnomAD file
%%     ~~Describe \code{input_file_name} here~~
}
  \item{export_file_name}{
       Output file name
%%     ~~Describe \code{export_file_name} here~~
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

## The function is currently defined as
function (input_file_name, export_file_name) 
{
    amino_code1 <- c("A", "R", "N", "D", "C", "Q", "E", "G", 
        "I", "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V", 
        "H")
    amino_code3 = c("Ala", "Arg", "Asn", "Asp", "Cys", "Gln", 
        "Glu", "Gly", "Ile", "Leu", "Lys", "Met", "Phe", "Pro", 
        "Ser", "Thr", "Trp", "Tyr", "Val", "His")
    names(amino_code1) <- amino_code3
    x <- fread(input_file_name)
    colnames(x)[10] <- "Protein_Consequence"
    x$ref <- substr(x$Protein_Consequence, 3, 5)
    x$Pos <- substr(x$Protein_Consequence, 6, nchar(x$Protein_Consequence) - 
        3)
    x$alt <- substr(x$Protein_Consequence, nchar(x$Protein_Consequence) - 
        2, nchar(x$Protein_Consequence))
    for (i in 1:nrow(x)) {
        x$change[i] <- paste(amino_code1[x$ref[i]], x$Pos[i], 
            amino_code1[x$alt[i]], sep = "")
    }
    y <- data.frame(Chromosome = x$Chromosome, Position = x$Position, 
        Reference = x$Reference, Alternate = x$Alternate, change = x$change)
    fwrite(y, export_file_name)
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
