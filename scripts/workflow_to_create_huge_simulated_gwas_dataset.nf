nextflow.enable.dsl=2

process generate_random_source_file {
    publishDir "${params.out}/generate_random_source_file", mode: 'rellink', overwrite: true

    output:
      path 'random_seed_file_source'
    script:
      seedval=1337
      """
      openssl enc -aes-256-ctr -pass pass:"$seedval" -nosalt </dev/zero 2>/dev/null | head -10000 > random_seed_file_source
      """
}


process split_dbsnp_into_smaller_chunks {
    publishDir "${params.out}/split_of_chunks", mode: 'rellink', overwrite: true

    input:
      path dbsnp151vcf
      val dbsnpsplits
    output:
      path 'chunk_*'
    script:
      """
      zcat $dbsnp151vcf | grep -v "#" | awk -vOFS="\t" '{print \$1,\$2,\$3,\$4,\$5}' > five-first-columns-no-header

      #split into dbsnpsplit number of unix split files
      split -dn ${dbsnpsplits} five-first-columns-no-header chunk_

      """
}

process select_random_snps {
    publishDir "${params.out}/select_random_snps", mode: 'rellink', overwrite: true
    input:
      tuple val(chunkname), path(chunk), path(random_seed_file_source)
      val(rsize)
    output:
      tuple val(chunkname), path(chunk_random)
    script:
      """
      sort -R --random-source=${random_seed_file_source} ${chunk} | head -n ${rsize} > chunk_random
      """
}

process convert_to_bfile_format {
    publishDir "${params.out}/convert_to_bfile_format", mode: 'rellink', overwrite: true
    input:
      path(pedfile)
      path(mapfile)
    output:
      tuple path('hm3_b36.bed'), path('hm3_b36.bim'), path('hm3_b36.fam'), emit: bfiles
    script:
      """
      plink --ped ${pedfile} --map ${mapfile} --make-bed --out hm3_b36
      """
}

process create_split_selection {
    publishDir "${params.out}/create_split_selection", mode: 'rellink', overwrite: true
   publishDir "${params.out}/convert_to_bfile_format", mode: 'rellink', overwrite: true
    input:
      tuple path('hm3_b36.bed'), path('hm3_b36.bim'), path('hm3_b36.fam')
      val hapmapsplits
    output:
      path 'chunk_*', emit: rsid_chunk
      env(size), emit: size_chunk
    script:
      """
      awk '{print \$2}' hm3_b36.bim > ids_to_split
      split -dn ${hapmapsplits} ids_to_split chunk_
      size=\$(wc -l chunk_00 | awk '{print \$1}')
      """
}

process extract_and_transform_hapmap_genotypes {
    publishDir "${params.out}/extract_and_transform_hapmap_genotypes", mode: 'rellink', overwrite: true
    input:
      tuple path('hm3_b36.bed'), path('hm3_b36.bim'), path('hm3_b36.fam'), val(id), path('to_extract')
    output:
      tuple val(id), path('chunk.raw')
    script:
      """
      plink2 --bfile hm3_b36 --extract to_extract --export A --out chunk
      """
}

process simulate_gwas_stats_from_hapmap_genotypes {
    publishDir "${params.out}/simulate_gwas_stats_from_hapmap_genotypes", mode: 'rellink', overwrite: true
    input:
      tuple val(id), path('hapmap_chunk.raw')
      path('Rsimscript.R')
    output:
      tuple val(id), path('logistic_testStats.txt'), path('linear_testStats.txt')
    script:
      """
      Rscript Rsimscript.R --afsource hapmap_chunk.raw --outlin linear_testStats.txt --outlog logistic_testStats.txt

      """
}

process merge_dbsnp_and_hapmap_chunks {
    publishDir "${params.out}/merge_dbsnp_and_hapmap_chunks", mode: 'rellink', overwrite: true
    input:
      tuple val(id), path('logistic_testStats.txt'), path('linear_testStats.txt'), path('dbsnpchunk')
    output:
      path('logistic_merged.txt'), emit: logistic
      path('linear_merged.txt'), emit: linear
    script:
      """
      echo -en "CHR\tBP\tRSID\tEA\tOA\t" > logistic_merged.txt
      head -n1 logistic_testStats.txt | awk '{\$1=""; \$2=""; \$0=\$0;\$1=\$1; print \$0}' >> logistic_merged.txt
      awk -vOFS="\t" -vlargerFile=dbsnpchunk 'NR>1{\$1=""; \$2=""; \$0=\$0;\$1=\$1; getline x<largerFile; print x, \$0}' logistic_testStats.txt >> logistic_merged.txt

      echo -en "CHR\tBP\tRSID\tEA\tOA\t" > linear_merged.txt
      head -n1 linear_testStats.txt | awk '{\$1=""; \$2=""; \$0=\$0;\$1=\$1; print \$0}' >> linear_merged.txt
      awk -vOFS="\t" -vlargerFile=dbsnpchunk 'NR>1{\$1=""; \$2=""; \$0=\$0;\$1=\$1; getline x<largerFile; print x, \$0}' linear_testStats.txt >> linear_merged.txt
      """
}

process concatenate_chunks  {
    publishDir "${params.out}/concatenate_chunks", mode: 'rellink', overwrite: true

      input:
        tuple val(type), path('test_chunks')
      output:
        path("${type}_tests_GWAS_GRCh38.txt")
      script:
        """
        # Concatenate tests
        for chunk in ${test_chunks}
        do
          cat \${chunk} >> ${type}_tests_GWAS_GRCh38.txt
        done
        """
}

workflow {

  Rsimscript=file("scripts/gwas_stat_simulation_script.R")

  println("reading file")
  println("${params.dbsnp151vcf}")
  println("${params.hapmap3pedmap}")

  // The workflow is designed to be able to use the same number of splits for dbsnp and hapmap
  // The number of splits should be increased +1 for each additional cpu used (default: 4)
  if (params.splits) { splits = params.splits }

  // Prepare hapmap part
  if(params.hapmap3pedmap) mapfile = file("${params.hapmap3pedmap}.map")
  if(params.hapmap3pedmap) pedfile = file("${params.hapmap3pedmap}.ped")

  // Prepare dbsnp part
  if(params.dbsnp151vcf) dbsnp151_channel = channel.fromPath("${params.dbsnp151vcf}")

  // Prepare general
  generate_random_source_file()

  // Run hapmap part
  convert_to_bfile_format(pedfile, mapfile)
  create_split_selection(convert_to_bfile_format.out, splits)
  create_split_selection.out.rsid_chunk
    .map { file -> tuple(file.baseName, file) }
    .transpose()
    .set { create_split_selection2 }
  convert_to_bfile_format.out.bfiles
    .combine(create_split_selection2)
    .set { ready_to_split_hapmap }
  extract_and_transform_hapmap_genotypes(ready_to_split_hapmap)

  // Simulate the GWAS from hapmap data
  simulate_gwas_stats_from_hapmap_genotypes(extract_and_transform_hapmap_genotypes.out, Rsimscript)
  
  // Run dbsnp part (creating an equally sized set of chunks)
  split_dbsnp_into_smaller_chunks(dbsnp151_channel, splits)
  split_dbsnp_into_smaller_chunks.out
    .map { file -> tuple(file.baseName, file) }
    .transpose()
    .set { split_dbsnp_into_smaller_chunks2 }
  split_dbsnp_into_smaller_chunks2.combine(generate_random_source_file.out).set { ch_random_selection }
  select_random_snps(ch_random_selection, create_split_selection.out.size_chunk)

  // Join dbsnp and sim gwas stats
  simulate_gwas_stats_from_hapmap_genotypes.out
    .join(select_random_snps.out)
    .set { to_merge }
  merge_dbsnp_and_hapmap_chunks(to_merge)

  // Concat all joined chunks
  merge_dbsnp_and_hapmap_chunks.out.linear.collect().map { chunks -> tuple('linear', chunks) }
    .mix(merge_dbsnp_and_hapmap_chunks.out.logistic.collect().map { chunks -> tuple('logistic', chunks) })
    .set { to_concatenate }
  concatenate_chunks(to_concatenate)

}

