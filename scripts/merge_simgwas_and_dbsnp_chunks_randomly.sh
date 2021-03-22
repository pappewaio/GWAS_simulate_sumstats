# Run using
#./scripts/merge_simgwas_and_dbsnp_chunks_randomly.sh data/dbsnp151/dbsnp_chr22_chunk data/gwas_sim_output/chr22.lin test.out

dbsnp_ann=$1
test_res=$2
outfile=$3
sep="\t"

#Add test_res on the left of dbsnp_ann. The stats will be totally disconnected from the annotation, that is important to remember. This data is only built to create test data for the purpose of software tests.

#check size of both data
echo "checking size of ${test_res}"
s1="$(awk 'END{print NR}' ${test_res})"
echo "checking size of ${dbsnp_ann}"
s2="$(awk 'END{print NR}' ${dbsnp_ann})"
if [ "${s1}" -gt "${s2}" ]; then
  # awk join the two files
  echo "merging the two files into a file of size: ${s2}"
  awk -vOFS="${sep}" -vlargerFile=${test_res} '{getline x<largerFile; print $0, x}' ${dbsnp_ann} > ${outfile}
elif [ "${s1}" -lt "${s2}" ]; then
  # awk join the two files
  echo "merging the two files into a file of size: ${s1}"
  awk -vOFS="${sep}" -vlargerFile=${dbsnp_ann} '{getline x<largerFile; print $0, x}' ${test_res} > ${outfile}
elif [ "${s1}" -eq "${s2}" ]; then
  # awk join the two files
  echo "merging the two files into a file of size: ${s1}"
  awk -vOFS="${sep}" -vlargerFile=${dbsnp_ann} '{getline x<largerFile; print $0, x}' ${test_res} > ${outfile}
else
  echo "something is wrong with the file size check. file1: ${1}, file2: ${2}"
  exit 1

fi


