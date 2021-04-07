## Download from web
#wget ftp://ftp.ncbi.nlm.nih.gov/snp/organisms/human_9606_b151_GRCh38p7/VCF/All_20180418.vcf

nextflow run scripts/workflow_to_create_huge_simulated_gwas_dataset.nf --dbsnp151vcf
