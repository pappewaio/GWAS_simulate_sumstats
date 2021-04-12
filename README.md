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

There is a markdown document generatd from rmarkdown which describes more in depth how the simulated data is created

