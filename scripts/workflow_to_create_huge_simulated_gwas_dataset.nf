nextflow.enable.dsl=2


process split_dbsnp_into_smaller_chunks {
    input:
      path dbsnp151vcf
    output:
      path 'chunk_*'
    script:
      dbsnpsplits=3
      """
      zcat $dbsnp151vcf | grep -v "#" | awk -vOFS="\t" '{print \$1,\$2,\$3,\$4,\$5}' > five-first-columns-no-header

      #split into dbsnpsplit number of unix split files
      split -dn ${dbsnpsplits} five-first-columns-no-header chunk_

      """
}

//process bar {
//    input:
//      path x
//    output:
//      path 'bar.txt'
//    script:
//      """
//      cat $x > bar.txt
//      """
//}

workflow {

  println("reading file")
  println("${params.dbsnp151vcf}")
  dbsnp151_channel = channel.fromPath("${params.dbsnp151vcf}")

  split_dbsnp_into_smaller_chunks(dbsnp151_channel)
  //bar(split_dbsnp_into_smaller_chunks.out)
}

