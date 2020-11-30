---
output:
  html_document:
    keep_md: true
---

# A tutorial on how to simulate GWAS summary statistics

The only prerequisite is to download a file with allele frequencies e.g., hapmap3.
Here, a subset of chr22 is used to demonstrate.

```
# Use awk to shrink the hapmap dataset (used to generate the data included with this repo)
awk '{out="";for(i=1;i<=125;i++){out=out" "$i}; print out}' hapmap3_r3_b36_chr22.raw > hapmap3_r3_b36_chr22.raw.shrunk

```

Load important R packages, set paths and make output directory


To simulate gwas summary statistics we use existing allele frequencies from a real dataset, and in this example the real data come from the public hapmap3 data. In this tutorial, we have tried to use speed and memory efficient code to allow the simulation to scale better in case we decide toe.g., use all 80 million SNPs from the 1000 genomes project. To quickly import the data we use ```fread()``` and then to allow for better memory management, and indirectly makes processing faster as well, we use the ```SummarizedExperiments``` object class to store data and results.


As we want to make three associations in the dataset, we need to make three betas that are not 0. We takes these betas from the normal distribution using the ```rnorm()``` function. (Andrew: what happens if we don't use rnorm, but just select three random numbers?).



```r
# Simulate some betas one for each variant (= one each row)
rowData(se)[["B.true"]] <- 0
# Simulate betas different from 0. Use three random values from the normal distribution.
rowData(se)[["B.true"]][sample(length(se), 3)] <- rnorm(3)
```

Now we will create our quantitative traits based on the genotypes(0, 1 or 2) and add the effect of the simulated "true" betas.


```r
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
```

## Generate summary statistics - Linear GWAS




```r
# Apply linear regression model and collect statistics
# This is a time critical step, so in this example we only use 100 rows
se2 <- se[1:100]
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
```



```r
# Add number of individuals and allele frequency
rowData(se2)[["N"]] <- rep(ncol(se2), nrow(se2))
rowData(se2)[["AFREQ"]] <- rowSums(assays(se2)[["geno"]])/(2*ncol(se2))

#write sumstats to file (skip true beta)
write.table(rowData(se2)[,-1], "out/linear_testStats.txt", sep="\t", quote=FALSE)
```

## Generate summary statistics - Logistic GWAS


```r
#re-initialize se2
se2 <- se[1:100]

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
write.table(rowData(se2)[,-1], "out/logistic_testStats.txt", sep="\t", quote=FALSE)
```

This html document can be regenerated by using rmarkdown


```r
## From within R
# Use data.table to read and and write fast using fread and fwrite
library(data.table)
# Use S4Vectors for a more advanced and efficient infrastructure
library(S4Vectors)
# Use SummarizedExperiment to create a container structure using the S4Vector framework
library(SummarizedExperiment)
# Use the rmarkdown package to generate/render the markdown/html document
library(rmarkdown)
render("gwas_simulate_sumstats_in_R.rmd")

#If it doesn't work, it might be because of an outdated pandoc
#If so try download and install a new one (here an example for a debian machine)
wget https://github.com/jgm/pandoc/releases/download/2.11.1.1/pandoc-2.11.1.1-1-amd64.deb
sudo dpkg -i pandoc-2.11.1.1-1-amd64.deb
pandoc --version
```

