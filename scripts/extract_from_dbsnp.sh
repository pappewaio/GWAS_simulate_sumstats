## Download from web
#wget ftp://ftp.ncbi.nlm.nih.gov/snp/organisms/human_9606_b151_GRCh38p7/VCF/All_20180418.vcf

# Select required fields 1-5 and split into files based on chr 1-22, MT, X and Y
zcat All_20180418.vcf.gz | grep -v "#" | awk -vOFS="\t" '{chr=$1; print $1,$2,$3,$4,$5 > "sim_"chr }'

seedval=1337
openssl enc -aes-256-ctr -pass pass:"$seedval" -nosalt </dev/zero 2>/dev/null | head -10000 > random_seed_file_source

# As example file in this tutorial in this repo we only need a subset of these chromosomes
for chr in {15..22} MT X Y; do \
for chr in {1..14}; do \
 ( \
 echo "$chr starting ..."; \
 echo -e "CHR\tPOS\tRSID\tEA\tOA" > dbsnp_chr${chr}_chunk
 sort -R --random-source=random_seed_file_source sim_${chr} | head -n 100 >> dbsnp_chr${chr}_chunk
 echo "$chr done ..."; \
 ) & \
done; wait

for chr in {1..22} MT X Y; do \
 echo "$chr starting ..."; \
 echo -e "CHR\tPOS\tRSID\tEA\tOA" > dbsnp_chr${chr}_chunk
 sort -R --random-source=random_seed_file_source sim_${chr} | head -n 100 >> dbsnp_chr${chr}_chunk
 echo "$chr done ..."; \
done

# Move these chunks into the data folder to be used as example for the merge script


# Generate vcf example dataset
zcat All_20180418.vcf.gz | head -n10000 | grep "#" > All_20180418_example_data.vcf
for chr in {3..22} MT X Y; do \
 ( \
 echo "$chr starting ..."; \
  zcat All_20180418.vcf.gz | grep -v "#" | awk -vchr=${chr} -vOFS="\t" '{if(chr==$1){print $0}}' | head -n200 >> "All_20180418_example_data.vcf" 
 echo "$chr done ..."; \
 ) & \
done


