##############################################################################################
# This script simulates GWAS stats from hapmap data using both linear and logistic regression
# It has no connection to genomic regions or LD
# REQUIRED PREINSTALLED LIBRARIES:
# - data.table
# - S4Vectors
# - SummarizedExperiment
# - getopt
##############################################################################################
# INPUT
# - afsource: hapmap data, for example "data/hapmap3_r3_b36_chr22.raw.shrunk"
# - outlin: file name of simulated linear regression data (default: linear_testStats.txt)
# - outlog: file name of simulated logistic regression data (default: logistic_testStats.txt)

#Collect script info
initial.options <- commandArgs(trailingOnly = F)
file.arg.name <- "--file="
script.name <- sub(file.arg.name, "", initial.options[grep(file.arg.name, initial.options)])
script.basename <- dirname(script.name)

suppressMessages(library('getopt'))

#input params
#0=no argument required
#1=argument required
#2=argument optional

spec = matrix(c(
'afsource', 'a', 1, "character",
'outlin' , 'b', 1, "character",
'outlog' , 'c', 1, "character",
'help' , 'h', 0, "logical"
), byrow=TRUE, ncol=4)

opt = getopt(spec)
# if help was asked for print a friendly message
# and exit with a non-zero error code
if ( !is.null(opt$help) ) {
  cat(getopt(spec, usage=TRUE))
  q(status=1)
}
#print(opt$afsource)
#print(opt$outlin)
#print(opt$outlog)
if(is.null(opt$afsource)){stop("source has to be specified")}else{ afsource <- opt$afsource }
if(is.null(opt$outlin)){outlin <- "linear_testStats.txt"}else{ outlin <- opt$outlin }
if(is.null(opt$outlog)){outlog <- "logistic_testStats.txt"}else{ outlog <- opt$outlog }
#print("hej")
#print(afsource)
#print(outlin)
#print(outlog)

# Use data.table to read and and write fast using fread and fwrite
suppressMessages(library(data.table))
# Use S4Vectors for a more advanced and efficient infrastructure
suppressMessages(library(S4Vectors))
# Use SummarizedExperiment to create a container structure using the S4Vector framework
suppressMessages(library(SummarizedExperiment))

# Select all SNPs in the file
dataraw <- fread(afsource)
# prepare genotype part of SE object, snps as rows and individuals as columns
# Note, in the df object snps are as columns
df <- DataFrame(dataraw)

# pick out info on rsid and effect allele
header <- colnames(df)[7:ncol(df)]
strspl <- tstrsplit(header, split="_")
rsid <- strspl[[1]]
effectAllele <- strspl[[2]]

#pick out the actual genotype information
geno <- t(as.matrix(df[,7:ncol(df)]))
rownames(geno) <- NULL
colnames(geno) <- df[,"FID"]
pheno <- df[,1:6]
rownames(pheno) <- df[,"FID"]

# Use chr, pos info here (but make up genomic ranges as they are missing in this dataset)
gr <- GRanges(seqnames=1:nrow(geno), ranges=IRanges(start=1:nrow(geno), end=1:nrow(geno)))
mcols(gr)["RSID"] <- rsid
mcols(gr)["effectAllele"] <- effectAllele

# Create a SummarizedExperiments object, which uses copy by replacement (= gentle to memory, and faster)
se <- SummarizedExperiment(assays=SimpleList(geno=geno),
                          rowData=NULL, rowRanges=gr,
                          colData=pheno,
                          metadata=list())

# remove initial data formatting to save space, keeping only the se
rm(dataraw, df, geno, pheno, gr, header, strspl, rsid, effectAllele)

# replace all NAs with zeros
assays(se)[["geno"]][is.na(assays(se)[["geno"]])] <- 0


# Simulate some betas one for each variant (= one each row)
rowData(se)[["B.true"]] <- 0
# Simulate betas different from 0. Use three random values from the normal distribution.
rowData(se)[["B.true"]][sample(length(se), 3)] <- rnorm(3)

# Multiply the matrix (Remember, nrows of left matrix has to be equal to ncols of right matrix).
colData(se)[["g.pheno"]] <- as.vector(t(assays(se)[["geno"]]) %*% rowData(se)[["B.true"]])
# Scale and center the data around the mean, but why? And what happens if we don't?
colData(se)[["g.pheno"]] <- (colData(se)[["g.pheno"]] - mean(colData(se)[["g.pheno"]]))/sd(colData(se)[["g.pheno"]])
# Sqrt(0.1) reduces the explained genotypic variance for the phenotype to 10%.
colData(se)[["g.pheno"]] <- sqrt(0.1)*colData(se)[["g.pheno"]]
# Make a corresponding environmental effect that is normal distributed random explaining 90% of the phenotype.
colData(se)[["e.pheno"]] <- rnorm(ncol(se))
colData(se)[["e.pheno"]] <- (colData(se)[["e.pheno"]] - mean(colData(se)[["e.pheno"]]))/sd(colData(se)[["e.pheno"]])
colData(se)[["e.pheno"]] <- sqrt(0.9)*colData(se)[["e.pheno"]]
# Combine the genotype and environment variance components
colData(se)[["q.pheno"]] <- colData(se)[["g.pheno"]] + colData(se)[["e.pheno"]]
# Generate binary phenotype using the liability model with a threshold for disease at 0.1 at the upper end of the normal distribution
colData(se)[["b.pheno"]] <- 1*( colData(se)[["q.pheno"]] >= qnorm( 0.1, lower.tail=F ) )

####
# Linear regression
####

# This is a time critical step
se2 <- se
rowData(se2)[["B"]] <- NA
rowData(se2)[["SE"]] <- NA
rowData(se2)[["Z"]] <- NA
rowData(se2)[["P"]] <- NA

# There are no packages that have vectorized the lm function, and therefore we have to use a for loop
for (i in 1:nrow(se2)){
  temp <- summary(lm(colData(se2)[,"q.pheno"] ~ as.vector(assays(se2[i])[["geno"]])))$coefficients
  rowData(se2[i])[["B"]] <- temp[2,1]
  rowData(se2[i])[["SE"]] <- temp[2,2]
  rowData(se2[i])[["Z"]] <- temp[2,3]
  rowData(se2[i])[["P"]] <- temp[2,4]
}

# Add number of individuals and allele frequency
rowData(se2)[["N"]] <- rep(ncol(se2), nrow(se2))
rowData(se2)[["AFREQ"]] <- rowSums(assays(se2)[["geno"]])/(2*ncol(se2))

#write simulated sumstats to file (skip true beta)
rowData(se2)[["B.true"]] <- NULL
write.table(rowData(se2), outlin, sep="\t", quote=FALSE, row.names=FALSE)

####
# Logisitc regression
####

#re-initialize se2 stats
se2 <- se

rowData(se2)[["B"]] <- NA
rowData(se2)[["SE"]] <- NA
rowData(se2)[["OR"]] <- NA 
rowData(se2)[["OR_u95"]] <- NA 
rowData(se2)[["OR_l95"]] <- NA  
rowData(se2)[["Z"]] <- NA
rowData(se2)[["P"]] <- NA

# Even slower than the linear regression above
for (i in 1:nrow(se2)){
  #Here there is a warnings that we suppress for the sake of the markdown document, but we need to fix this later
  ### Warning: glm.fit: fitted probabilities numerically 0 or 1 occurred
  ### Warning in regularize.values(x, y, ties, missing(ties)): collapsing to unique 'x' values
  suppressWarnings(temp <- glm(colData(se2)[,"b.pheno"] ~ as.vector(assays(se2[i])[["geno"]]), family="binomial"))
  temp2 <- coef(summary(temp))
  rowData(se2[i])[["B"]] <- temp2[2,1]
  rowData(se2[i])[["SE"]] <- temp2[2,2]
  rowData(se2[i])[["Z"]] <- temp2[2,3]
  rowData(se2[i])[["P"]] <- temp2[2,4]
  suppressWarnings(rowData(se2[i])[["OR"]] <- exp(coef(temp))[2])
  suppressWarnings(suppressMessages(rowData(se2[i])[["OR_u95"]] <- exp(confint(temp) )[2,2]))
  suppressWarnings(suppressMessages(rowData(se2[i])[["OR_l95"]] <- exp(confint(temp) )[2,1]))
}

# Add number of individuals and allele frequency
rowData(se2)[["N"]] <- rep(ncol(se2), nrow(se2))
rowData(se2)[["Ncase"]] <- rep(table(colData(se2)[,"b.pheno"])[2], nrow(se2))
rowData(se2)[["Ncont"]] <- rep(table(colData(se2)[,"b.pheno"])[1], nrow(se2))
rowData(se2)[["AFREQ"]] <- rowSums(assays(se2)[["geno"]])/(2*ncol(se2))

#write sumstats to file (skip true beta)
rowData(se2)[["B.true"]] <- NULL
write.table(rowData(se2), outlog, sep="\t", quote=FALSE, row.names=FALSE)
