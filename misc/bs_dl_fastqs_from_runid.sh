#!/bin/bash
run_name=$1
config=$2

## get run id
run_id=$(bs --filter-field=ExperimentName --filter-term=$run_name -c $config -f csv | cut -f1 | tail -1 | cut -d, -f2)

## get dataaset IDs
./bs list datasets --input-run=$run_id -F Name -F Id -F AppSession.Application.Name -F AppSession.DateCreated -c $config -f csv | tail -n +2 | sort -t, -k4r | grep -v "Undetermined" > files

## get most recent dataset ids in case of multiple fastq generations
most_recent=$(cut -d, -f4 files | uniq | head -n 1)
grep "$most_recent" files | cut -d, -f 2 > dataset_ids

## download fastqs
while read line; do bs download dataset -i $line --extension=.fastq.gz; done < dataset_ids

## clean up
rm files
rm dataset_ids