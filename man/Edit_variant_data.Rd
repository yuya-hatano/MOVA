\name{Edit_variant_data}
\alias{Edit_variant_data}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Integrate HGMD file and gnomAD file.  
%%  ~~function to do ... ~~
}
\description{
Integrate the HGMD file processed by the Hgmd_divide function and the gnomAD file processed by the Edit_gnomAD_file function. In addition, this function also adds information on predictions obtained from existing in silico algorithms other than PolyPhen-2, a necessary step for MOVA execution.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
Edit_variant_data(AlphScore_final_file, hgmd_file_name, gnomad_file_name, uniprot_name, gene_name, export_file_name)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{AlphScore_final_file}{
       Path to AlphScore_final.tsv. AlphScore_final.tsv is available at https://zenodo.org/record/6288139#.ZA1Zo3bP3rc.
%%     ~~Describe \code{AlphScore_final_file} here~~
}
  \item{hgmd_file_name}{
       The HGMD file processed by the Hgmd_divide function
%%     ~~Describe \code{hgmd_file_name} here~~
}
  \item{gnomad_file_name}{
       The gnomAD file processed by the Edit_gnomAD_file function
%%     ~~Describe \code{gnomad_file_name} here~~
}
  \item{uniprot_name}{
       Target protein uniprot name
%%     ~~Describe \code{uniprot_name} here~~
}
  \item{gene_name}{
       Target gene name
%%     ~~Describe \code{gene_name} here~~
}
  \item{export_file_name}{
       Name of the file to be exported.
%%     ~~Describe \code{export_file_name} here~~
}
}
\details{
  In addition to the main output file (export_file_name), "gene name_alph.csv" is also output. edit_polyphen_data function must be applied to the output file of the Edit_variant_data function. The main output file should be specified in argument: input_file_name1, and "gene name_alph.csv" in argument: input_file_name2.
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
function (AlphScore_final_file, hgmd_file_name, gnomad_file_name, 
    uniprot_name, gene_name, export_file_name) 
{
    x <- fread(AlphScore_final_file)
    x <- x[x$Uniprot_acc_split == uniprot_name, ]
    colnames(x)[1] <- "Chromosome"
    colnames(x)[2] <- "Position"
    colnames(x)[3] <- "Reference"
    colnames(x)[4] <- "Alternate"
    x[, `:=`(DEOGEN2_score, NULL)]
    x[, `:=`(b_factor, NULL)]
    x[, `:=`(SOLVENT_ACCESSIBILITY_core, NULL)]
    x[, `:=`(in_gnomad_train, NULL)]
    x[, `:=`(in_clinvar_ds, NULL)]
    x[, `:=`(glm_AlphCaddDeogen, NULL)]
    x[, `:=`(glm_AlphDeogenRevel, NULL)]
    x[, `:=`(glm_DeogenRevel, NULL)]
    x[, `:=`(glm_CaddDeogen, NULL)]
    x[, `:=`(glm_AlphDeogen, NULL)]
    x[, `:=`(glm_AlphRevelCadd, NULL)]
    x[, `:=`(glm_CaddDeogenRevel, NULL)]
    fwrite(x, paste(gene_name, "_alph.csv", sep = ""))
    x <- fread(hgmd_file_name)
    y <- fread(gnomad_file_name)
    D <- merge(x, y, all = T)
    D <- D[(D$possible != "DMp") | (is.na(D$possible)), ]
    D$Type <- replace(D$Type, which(is.na(D$Type)), "Ctrl")
    for (i in 1:nrow(D)) {
        D$change[i] <- replace(D$change[i], which(is.na(D$change[i])), 
            D$HGVSprotein[i])
    }
    D[, `:=`(HGVSprotein, NULL)]
    fwrite(D, paste(gene_name, "_variant.csv", sep = ""))
    x <- fread(paste(gene_name, "_alph.csv", sep = ""))
    y <- fread(paste(gene_name, "_variant.csv", sep = ""))
    D <- merge(x, y)
    D$Type2 <- 0
    D$Type3 <- 1
    D[D$Type == "Target", ]$Type2 <- 1
    D[D$Type == "Target", ]$Type3 <- 1
    D[D$Type == "Ctrl", ]$Type2 <- 0
    D[D$Type == "Ctrl", ]$Type3 <- 0
    fwrite(D, export_file_name)
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
