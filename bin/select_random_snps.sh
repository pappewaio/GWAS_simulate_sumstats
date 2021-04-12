chunk=$1
random_seed_file_source=$2
rsize=$3
sort -R --random-source=${random_seed_file_source} ${chunk} | head -n ${rsize}

