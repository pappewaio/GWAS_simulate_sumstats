# Prepare example data for this repo
#download plink ped and map files from
# wget -r --no-parent https://ftp.ncbi.nlm.nih.gov/hapmap/genotypes/hapmap3_r3/plink_format/*
 
# convert to plink bfile format
plink --ped hapmap3_r3_b36_fwd.consensus.qc.poly.ped --map hapmap3_r3_b36_fwd.consensus.qc.poly.map --make-bed --out hapmap3_r3_b36

# convert to a format suitable to read from R
for chr in {1..22} X; do
 echo "chromosome $chr"
 plink2 -bfile hapmap3_r3_b36 --chr ${chr} --export A --out hapmap3_r3_b36_chr${chr} --threads 10
done


# Use awk to shrink the hapmap dataset (used to generate the data included with this repo)
size=125
for chr in {1..22} X; do
 echo "chromosome $chr"
 awk -vsize=${size} '{out="";for(i=1;i<=size;i++){out=out" "$i}; print out}' hapmap3_r3_b36_chr${chr}.raw > hapmap3_r3_b36_chr${chr}.raw.shrunk
done


