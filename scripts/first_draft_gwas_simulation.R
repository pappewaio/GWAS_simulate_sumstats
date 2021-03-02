######
## Stat conversion formulas to be implemented in sumstats cleaner
######

## Output Columns
##	0 CHR POS RSID EffectAllele OtherAllele B SE OR OR_u95 OR_l95 Z P N	N_inf N_case N_Control INFO AFREQ	

library( data.table )
setwd( "/Users/AndrewSchork/Desktop/Projects/SumStats_Pipeline/ConversionFormulas" )

## Generate Test Data

## Load test data

# Select 1000 SNPs, fill in missing data
data <- fread( 'hapmap3_r3_b36_chr22.raw' )
data <- data[ ,7:1006 ]
data[ is.na( data ) ] <- 0

# Simulate some betas
beta.true <- rep( 0, 1000 )
beta.true[ sample( 1000, 3 ) ] <- rnorm( 3 )

# Convert to quantative and binary traits
g.pheno <- as.matrix( data ) %*% beta.true
g.pheno <- sqrt( 0.1 )*( ( g.pheno - mean( g.pheno) ) / sd( g.pheno ) )
e.pheno <- rnorm( 1397 )
e.pheno <- sqrt( 0.9 )*( ( e.pheno - mean( e.pheno) ) / sd( e.pheno ) )
q.pheno <- g.pheno + e.pheno
b.pheno <- 1*( q.pheno >= qnorm( 0.1, lower.tail=F ) )

## Generate summary statistics - Linear GWAS

sumstats <- data.table( names( data) )
names( sumstats ) <- "RSID"

sumstats$B <- NA 
sumstats$SE <- NA 
sumstats$Z <- NA 
sumstats$P <- NA 
sumstats$N <- NA 
sumstats$AFREQ <- NA 

for ( i in 1:1000 ) {
	
	temp <- lm( q.pheno ~ data[[ i ]] )
	sumstats$B[[ i ]] <- summary( temp )$coefficients[ 2,1 ]
	sumstats$SE[[ i ]] <- summary( temp )$coefficients[ 2,2 ]
	sumstats$Z[[ i ]] <- summary( temp )$coefficients[ 2,3 ]
	sumstats$P[[ i ]] <- summary( temp )$coefficients[ 2,4 ]

}

sumstats$N <- rep( 1397, 1000 )
sumstats$AFREQ <- colSums( data ) / ( 2*1397 )

fwrite( sumstats, "linear_testStats.txt")

xaxis <- -log10( seq( 1,1000 ) / 1001 )
yaxis <- -log10( sumstats$P[ order( sumstats$P ) ] )
plot( xaxis, yaxis )
abline(0,1)

## Generate summary statistics - Logistic GWAS

sumstats <- data.table( names( data) )
names( sumstats ) <- "RSID"

sumstats$B <- NA 
sumstats$SE <- NA
sumstats$OR <- NA 
sumstats$OR_u95 <- NA 
sumstats$OR_l95 <- NA  
sumstats$Z <- NA 
sumstats$P <- NA 
sumstats$N <- NA
sumstats$Ncase <- NA 
sumstats$Ncont <- NA  
sumstats$AFREQ <- NA 

for ( i in 1:1000 ) {
	
	temp <- glm( b.pheno ~ data[[ i ]], family="binomial" )
	sumstats$B[[ i ]] <- summary( temp )$coefficients[ 2,1 ]
	sumstats$SE[[ i ]] <- summary( temp )$coefficients[ 2,2 ]
	sumstats$OR[[ i ]] <- exp( coef( temp ) )[ 2 ]
	sumstats$OR_u95[[ i ]] <- exp( confint( temp ) )[ 2,2 ]
	sumstats$OR_l95[[ i ]] <- exp( confint( temp ) )[ 2,1 ]
	sumstats$Z[[ i ]] <- summary( temp )$coefficients[ 2,3 ]
	sumstats$P[[ i ]] <- summary( temp )$coefficients[ 2,4 ]

}

sumstats$N <- rep( 1397, 1000 )
sumstats$Ncase <- rep( table( b.pheno )[ 2 ], 1000 )
sumstats$Ncont <- rep( table( b.pheno )[ 1 ], 1000 )
sumstats$AFREQ <- colSums( data ) / ( 2*1397 )

fwrite( sumstats, "logistic_testStats.txt")

xaxis <- -log10( seq( 1,1000 ) / 1001 )
yaxis <- -log10( sumstats$P[ order( sumstats$P ) ] )
plot( xaxis, yaxis )
abline(0,1)

save( data, g.pheno, e.pheno, q.pheno, b.pheno, beta.true, file="simData.RData")




