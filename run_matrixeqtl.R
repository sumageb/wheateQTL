library(MatrixEQTL)
useModel = modelLINEAR
snps_file = paste("genotypes.txt", sep= "\t")
snps_location_file = paste("final_snps_position.txt", sep="\t")
expression_file = paste("tmm_ematrix_filtered.txt", sep="\t")
gene_location_file = paste("final_annofile_filtered.txt", sep= "\t")
covariates = paste("covariates.txt", sep="\t")

output_file_name_cis = NULL
output_file_name_tra = NULL

#Only associations significant at this level will be saved
pvOutputThreshold_cis = 0.00005
pvOutputThreshold_tra = 0.00005

# Distance for local gene-SNP pairs
cisDist = 1e6           #maximum distance at which gene-SNP pair is considered local.
#loading genotype data
snps = SlicedData$new();
snps$fileDelimiter = "\t";      # the TAB character
snps$fileOmitCharacters = "NA"; # denote missing values;
snps$fileSkipRows = 1;          # one row of column labels
snps$fileSkipColumns = 1;       # two column of row labels
snps$fileSliceSize = 2000;      # read file in pieces of 2,000 rows
snps$LoadFile(snps_file)

#loading gene expression data
gene = SlicedData$new();
gene$fileDelimiter = "\t";      # the TAB character
gene$fileOmitCharacters = "NA"; # denote missing values;
gene$fileSkipRows = 1;          # one row of column labels
gene$fileSkipColumns = 1;       # two column of row labels
gene$fileSliceSize = 2000;      # read file in slices of 2,000 rows
gene$LoadFile(expression_file)

## Load covariates
cvrt = SlicedData$new();
cvrt$fileDelimiter = "\t";      # the TAB character
cvrt$fileOmitCharacters = "NA"; # denote missing values;
cvrt$fileSkipRows = 1;          # one row of column labels
cvrt$fileSkipColumns = 1;       # one column of row labels
if(length(covariates)>0) {
cvrt$LoadFile(covariates);
}

## Run the analysis
snpspos = read.table(snps_location_file, header = TRUE, stringsAsFactors = FALSE);
genepos = read.table(gene_location_file, header = TRUE, stringsAsFactors = FALSE)

matrix_eqtl = Matrix_eQTL_main(
snps = snps,
gene = gene,
cvrt = cvrt,
output_file_name = output_file_name_tra,
pvOutputThreshold = pvOutputThreshold_tra,
useModel = useModel,
verbose = TRUE,
output_file_name.cis = output_file_name_cis,
pvOutputThreshold.cis = pvOutputThreshold_cis,
snpspos = snpspos,
genepos = genepos,
cisDist = cisDist,
min.pv.by.genesnp = TRUE,
noFDRsaveMemory = FALSE)
save(matrix_eqtl, file="meqtl.Rdata")
