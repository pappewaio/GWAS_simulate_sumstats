# Run using
./scripts/merge_simgwas_and_dbsnp_chunks_randomly.sh data/dbsnp151/dbsnp_chr22_chunk data/gwas_sim_output/chr22.lin | head

test_res=$1
dbsnp_ann=$2

#Add test_res on the left of dbsnp_ann. The stats will be totally disconnected from the annotation, that is important to remember. This data is only built to create test data for the purpose of software tests.

#check size of both data
s1=$(wc -l ${test_res})
s2=$(wc -l ${dbsnp_ann})

if [ ${s1} == ${s2} ]; do

  paste ${test_res} ${dbsnp_ann}

else

  echo "Not same file sizem, which is required. file1: ${1} versus file2: ${2}"
  exit 1

done


