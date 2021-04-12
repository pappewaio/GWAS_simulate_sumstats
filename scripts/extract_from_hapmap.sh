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


####
# Make plink ped and map file to mimic the source, so that we can use it in our example data
####

# make a list of snps to subset
for chr in {1..22} X; do
  echo "chr ${chr} done"
  if [ "${chr}" == "X" ];  then    chr2=23; else    chr2=${chr}; fi
  cat hapmap3_r3_b36_fwd.consensus.qc.poly.map | awk -vchr=${chr2} -vOFS="\t" '{if(chr==$1){print $2}}' | head -n200 > source_example_data/to_extract_${chr}
done

for chr in {1..22} X; do \
 ( \
 echo "$chr starting ..."; \
  if [ "${chr}" == "X" ];  then    chr2=23; else    chr2=${chr}; fi
  plink --chr ${chr2} --ped hapmap3_r3_b36_fwd.consensus.qc.poly.ped --map hapmap3_r3_b36_fwd.consensus.qc.poly.map --extract source_example_data/to_extract_${chr} --make-bed --out source_example_data/test_${chr}
 echo "$chr done ..."; \
 ) & \
done

# Merge files 
rm source_example_data/test_merge_list
#for chr in {2..22} X; do
#skip X chromosome due to problems with the haploied heterygous markers
for chr in {2..22}; do
  echo "source_example_data/test_${chr}" >> source_example_data/test_merge_list
done
plink --bfile source_example_data/test_1 --merge-list source_example_data/test_merge_list --out source_example_data/test_merge

# Convert merged file back to ped map
plink --bfile source_example_data/test_merge --recode --out source_example_data/test_merge


