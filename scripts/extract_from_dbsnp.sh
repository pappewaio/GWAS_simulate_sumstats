## Download from web
#wget ftp://ftp.ncbi.nlm.nih.gov/snp/organisms/human_9606_b151_GRCh38p7/VCF/All_20180418.vcf

# Select required fields 1-5 and split into files based on chr 1-22, MT, X and Y
mkdir -p simulated_gwas_coordinates
zcat All_20180418.vcf.gz | grep -v "#" | awk '{chr=$1; print $1,$2,$3,$4,$5 > "simulated_gwas_coordinates/sim_"chr }'

cd simulated_gwas_coordinates

seedval=1337
openssl enc -aes-256-ctr -pass pass:"$seedval" -nosalt </dev/zero 2>/dev/null | head -1000 > random_seed_file_source

# As example file in this tutorial in this repo we only need a subset of these chromosomes
for chr in {1..22} MT X Y; do \
 ( \
 echo "$chr starting ..."; \
 sort -R --random-source=random_seed_file_source sim_${chr} | head -n 100 > dbsnp_chr${chr}_chunk
 echo "$chr done ..."; \
 ) & \
done; wait

# Move these chunks into the data folder to be used as example for the merge script

