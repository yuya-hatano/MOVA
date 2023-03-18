\name{Edit_polyphen_data}
\alias{Edit_polyphen_data}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
This function adds PolyPhen-2 predictions to the MOVA input file.
%%  ~~function to do ... ~~
}
\description{
This function adds PolyPhen-2 predictions to the MOVA input file. This function is required to create the input files needed for MOVA functions.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
Edit_polyphen_data(dbNSFP4.3a_final_file, input_file_name1, input_file_name2, uniprot_name, position1, position2, export_file_name1, export_file_name2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dbNSFP4.3a_final_file}{
       Path of the dbNSFP4.3a_variant file for the Chromosome in which the gene resides.
%%     ~~Describe \code{dbNSFP4.3a_final_file} here~~
}
  \item{input_file_name1}{
       Path of export_file_name output by the Edit_variant_data function.
%%     ~~Describe \code{input_file_name1} here~~
}
  \item{input_file_name2}{
       Path of "gene name_alph.csv" output by the Edit_variant_data function.
%%     ~~Describe \code{input_file_name2} here~~
}
  \item{uniprot_name}{
       uniprot name of target protein
%%     ~~Describe \code{uniprot_name} here~~
}
  \item{position1}{
       Lower limit position of the targeted gene (hg38)
%%
}
  \item{position2}{
       Upper limit position of the targeted gene (hg38)
%%     ~~Describe \code{position2} here~~
}
  \item{export_file_name1}{
       Name of output file 1
%%     ~~Describe \code{export_file_name1} here~~
}
  \item{export_file_name2}{
       Name of output file 2
%%     ~~Describe \code{export_file_name2} here~~
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
function (dbNSFP4.3a_final_file, input_file_name1, input_file_name2, 
    uniprot_name, position1, position2, export_file_name1, export_file_name2) 
{
    D <- fread(dbNSFP4.3a_final_file)
    head(D)
    colnames(D)[1] <- "Chromosome"
    colnames(D)[2] <- "Position"
    colnames(D)[3] <- "Reference"
    colnames(D)[4] <- "Alternate"
    dren <- D[(D$Position <= position2) & (D$Position >= position1), 
        ]
    dt <- data.table(Chromosome = dren$Chromosome, Position = dren$Position, 
        Reference = dren$Reference, Alternate = dren$Alternate, 
        Uniprot_acc = dren$Uniprot_acc, Polyphen_prob = dren$Polyphen2_HDIV_score)
    Uniplot_v <- str_split(dt$Uniprot_acc, pattern = ";", simplify = TRUE)
    Polyphen2_HDIV_v <- str_split(dt$Polyphen_prob, pattern = ";", 
        simplify = TRUE)
    dt$pph_prob <- ""
    i <- 1
    while (i < nrow(dt)) {
        CHm <- match(uniprot_name, Uniplot_v[i, ])
        if (is.na(CHm) == FALSE) {
            dt$pph_prob[i] <- Polyphen2_HDIV_v[i, CHm]
            i <- i + 1
        }
        else {
            dt <- dt[-i, ]
            Uniplot_v <- Uniplot_v[-i, ]
            Polyphen2_HDIV_v <- Polyphen2_HDIV_v[-i, ]
        }
    }
    dt[, `:=`(Uniprot_acc, NULL)]
    dt[, `:=`(Polyphen_prob, NULL)]
    dat <- fread(input_file_name1)
    pph_data <- merge(dat, dt)
    fwrite(pph_data, export_file_name1)
    dat <- fread(input_file_name2)
    pph_data <- merge(dat, dt)
    fwrite(pph_data, export_file_name2)
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
