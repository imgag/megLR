#!/bin/bash
# Change samplename for all files in a raw data folder (excluding .fast5)

old_sample=21073LRa141_02787  # Give full old samplename!
new_sample=21102LRa041_02782  # Full new samplename

old_project=21073
new_project=21102

find . -type f -not -name '*.fast5' -not -name '*.fastq.gz' -exec sed -i "s/$old_sample/$new_sample/g" {} +
find . -type f -not -name '*.fast5' -not -name '*.fastq.gz' -exec sed -i "s/$old_project/$new_project/g" {} +

change_fastq () {
    tmp=$(mktemp)
    echo "Modifying file $1"
    cp $1 $tmp
    gunzip $tmp 
    sed -i "s/$old_sample/$new_sample/g" $tmp
    sed -i "s/$old_project/$new_project/g" $tmp
    gzip $tmp > $1
    rm $tmp
}

export old
export new
export -f change_fastq

find . -type f -name '*.fastq.gz' -exec bash -c 'change_fastq "{}"' \;