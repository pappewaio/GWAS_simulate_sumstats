
for chr in {1..22} X; do \
 ( \
 echo "$chr starting ..."; \
 Rscript scripts/gwas_stat_simulation_script.R --afsource data/hapmap3_r3_b36/hapmap3_r3_b36_chr${chr}.raw.shrunk --outlin chr${chr}.lin --outlog chr${chr}.log
 echo "$chr done ..."; \
 ) & \
done; wait

