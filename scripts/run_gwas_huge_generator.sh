## Download from web the whole set from
#wget ftp://ftp.ncbi.nlm.nih.gov/snp/organisms/human_9606_b151_GRCh38p7/VCF/All_20180418.vcf

# Run the exampled data that mimics the source
nextflow run scripts/workflow_to_create_huge_simulated_gwas_dataset.nf \
       	--dbsnp151vcf data/dbsnp151/All_20180418_example_data.vcf.gz \
       	--hapmap3pedmap data/hapmap3_r3_b36/original_input_reduced/test_merge \
       	--out out/nf-data-gen

