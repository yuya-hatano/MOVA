source("MOVA.r")

#TARDBP
#Be sure to check the contents of the gnomAD file yourself and delete all but the missense mutations.
Edit_gnomAD_file("gnomAD_v3.1.2_ENST00000240185.csv", "TARDBP_gnomAD.csv")
#Please download AlphScore_final.tsv from the following URL.
#https://zenodo.org/record/6288139
#TARDBP_hgmd.csv must be created from the HGMD data yourself.
Edit_variant_data("../CADD_REVEL/AlphScore_final.tsv", "TARDBP_hgmd.csv", "TARDBP_gnomAD.csv", "Q13148", "TARDBP", "TARDBPvariantdata.csv")

#Please download AlphaFold2 pdb files from the following URL.
#https://alphafold.ebi.ac.uk/
MOVA("Q13148.fa",  "TARDBP", "AF-Q13148-F1-model_v2.pdb", "TARDBPvariantdata.csv")
#Please create TARDBPpolyphen2.csv by yourself referring to the following URL.
#http://genetics.bwh.harvard.edu/pph2/bgi.shtml
PolyPhen_MOVA("TARDBPpolyphen2.csv")
AlphScore_MOVA("TARDBPvariantdata.csv")
CADD_MOVA("TARDBPvariantdata.csv")
REVEL_MOVA("TARDBPvariantdata.csv")
EVE_MOVA("TARDBPvariantdata.csv", "TADBP_HUMAN.csv")
Com_CADD_MOVA("TARDBP_Target_predict.csv")
Com_REVEL_MOVA("TARDBP_Target_predict.csv")


#FUS
Edit_gnomAD_file("gnomAD_v3.1.2_ENST00000254108.csv", "FUS_gnomAD.csv")
Edit_variant_data("../CADD_REVEL/AlphScore_final.tsv", "FUS_hgmd.csv", "FUS_gnomAD.csv", "P35637", "FUS", "FUSvariantdata.csv")


MOVA("P35637.fa",  "FUS", "AF-P35637-F1-model_v2.pdb", "FUSvariantdata.csv")
PolyPhen_MOVA("FUSpolyphen2.csv")
AlphScore_MOVA("FUSvariantdata.csv")
CADD_MOVA("FUSvariantdata.csv")
REVEL_MOVA("FUSvariantdata.csv")
Com_CADD_MOVA("FUS_Target_predict.csv")
Com_REVEL_MOVA("FUS_Target_predict.csv")


#SETX
Edit_gnomAD_file("gnomAD_v3.1.2_ENST00000224140.csv", "SETX_gnomAD.csv")
Edit_variant_data("../CADD_REVEL/AlphScore_final.tsv", "SETX_hgmd.csv", "SETX_gnomAD.csv", "Q7Z333", "SETX", "SETXvariantdata.csv")


MOVA("Q7Z333.fa",  "SETX", "AF-Q7Z333-F1-model_v2.pdb", "SETXvariantdata.csv")
PolyPhen_MOVA("SETXpolyphen2.csv")
AlphScore_MOVA("SETXvariantdata.csv")
CADD_MOVA("SETXvariantdata.csv")
REVEL_MOVA("SETXvariantdata.csv")
EVE_MOVA("SETXvariantdata.csv", "SETX_HUMAN.csv")
Com_CADD_MOVA("SETX_Target_predict.csv")
Com_REVEL_MOVA("SETX_Target_predict.csv")

#SETX_ALL_Phenotype
MOVA("Q7Z333.fa",  "SETX", "AF-Q7Z333-F1-model_v2.pdb", "SETXvariantdata.csv", phenotype = "All")
PolyPhen_MOVA("SETXpolyphen2.csv", phenotype = "All")
AlphScore_MOVA("SETXvariantdata.csv", phenotype = "All")
CADD_MOVA("SETXvariantdata.csv", phenotype = "All")
REVEL_MOVA("SETXvariantdata.csv", phenotype = "All")
EVE_MOVA("SETXvariantdata.csv", "SETX_HUMAN.csv", phenotype = "All")



#OPTN
Edit_gnomAD_file("gnomAD_v3.1.2_ENST00000378748.csv", "OPTN_gnomAD.csv")
Edit_variant_data("../CADD_REVEL/AlphScore_final.tsv", "OPTN_hgmd.csv", "OPTN_gnomAD.csv", "Q96CV9", "OPTN", "OPTNvariantdata.csv")


MOVA("Q96CV9.fa",  "OPTN", "AF-Q96CV9-F1-model_v2.pdb", "OPTNvariantdata.csv")
PolyPhen_MOVA("OPTNpolyphen2.csv")
AlphScore_MOVA("OPTNvariantdata.csv")
CADD_MOVA("OPTNvariantdata.csv")
REVEL_MOVA("OPTNvariantdata.csv")
EVE_MOVA("OPTNvariantdata.csv", "OPTN_HUMAN.csv")
Com_CADD_MOVA("OPTN_Target_predict.csv")
Com_REVEL_MOVA("OPTN_Target_predict.csv")



#TBK1
Edit_gnomAD_file("gnomAD_v3.1.2_ENST00000331710.csv", "TBK1_gnomAD.csv")
Edit_variant_data("../CADD_REVEL/AlphScore_final.tsv", "TBK1_hgmd.csv", "TBK1_gnomAD.csv", "Q9UHD2", "TBK1", "TBK1variantdata.csv")


MOVA("Q9UHD2.fa",  "TBK1", "AF-Q9UHD2-F1-model_v2.pdb", "TBK1variantdata.csv")
PolyPhen_MOVA("TBK1polyphen2.csv")
AlphScore_MOVA("TBK1variantdata.csv")
CADD_MOVA("TBK1variantdata.csv")
REVEL_MOVA("TBK1variantdata.csv")
Com_CADD_MOVA("TBK1_Target_predict.csv")
Com_REVEL_MOVA("TBK1_Target_predict.csv")



#SOD1
Edit_gnomAD_file("gnomAD_v3.1.2_ENST00000270142.csv", "SOD1_gnomAD.csv")
Edit_variant_data("../CADD_REVEL/AlphScore_final.tsv", "SOD1_hgmd.csv", "SOD1_gnomAD.csv", "P00441", "SOD1", "SOD1variantdata.csv")


MOVA("P00441.fa",  "SOD1", "AF-P00441-F1-model_v2.pdb", "SOD1variantdata.csv")
PolyPhen_MOVA("SOD1polyphen2.csv")
AlphScore_MOVA("SOD1variantdata.csv")
CADD_MOVA("SOD1variantdata.csv")
REVEL_MOVA("SOD1variantdata.csv")
EVE_MOVA("SOD1variantdata.csv", "SODC_HUMAN.csv")
Com_CADD_MOVA("SOD1_Target_predict.csv")
Com_REVEL_MOVA("SOD1_Target_predict.csv")

