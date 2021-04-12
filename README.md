# GWAS_simulate_sumstats

Access all required software from this conda environment
```
conda env create -f my_env.yml
conda activate gwas_stat_simulation_2
```

There is a script to produce gwas stats from hapmap data, see below for example of how to run
```
Rscript scripts/gwas_stat_simulation_script.R --afsource data/hapmap3_r3_b36_chr22.raw.shrunk --outlin out/linear_testStats.txt --outlog out/logistic_testStats.txt
```

To create a complete GWAS simulated dataset using stats from hapmap and annotation from dbsnp, please use the workflow executed by this command (takes about 14 minutes to run for these ~4000 variants)
```
# Run the exampled data that mimics the source (splits = nr. of cpus)
nextflow run main.nf \
       	--dbsnp151vcf data/dbsnp151/All_20180418_example_data.vcf.gz \
       	--hapmap3pedmap data/hapmap3_r3_b36/original_input_reduced/test_merge \
       	--splits 4 \
       	--out out/nf-data-gen
```

All intermediate files are saved in the output folder, and the final result is located in 'out/nf-data-gen/concatenate_chunks'. For the original source see:

```
## Download from web the whole set of variants from dbSNP
wget ftp://ftp.ncbi.nlm.nih.gov/snp/organisms/human_9606_b151_GRCh38p7/VCF/All_20180418.vcf

## Download from web the whole set of variants from HapMap3
wget -r --no-parent https://ftp.ncbi.nlm.nih.gov/hapmap/genotypes/hapmap3_r3/plink_format/*
```


_There is a markdown document generatd from rmarkdown which describes more in depth how the simulated data is created, see the tutorial folder._

