#!/bin/bash

## Set paths to code binaries
code_dir="/home/ec2-user/Workspace/testingOpenMP"
code_00="${code_dir}/00_splitChromosome/00_splitChromosome"
code_01="${code_dir}/01_modScores/01_modScores"

## Set inputs such as input correlation file, max length of LD block, and max number of variants in a chunk for parallel processing
input_file="$1"
model="$3"
max_length="5000000"
max_vars="100000"

## Set output directory and create it if it doesn't exist already
output_dir="$2"

## Do not change the name of this file -- partitions is the name of the file the C++ script will output
partitions="${output_dir}partitions.txt"
mkdir -p ${output_dir}

echo "Running chunking script"

## Run script to preprocess input file and identify chunk ranges for parallel or sequential processing
${code_00} -i ${input_file} -m ${model} -l ${max_length} -v ${max_vars} -o ${output_dir} 

echo "Chunking script done."
echo "Beginning modularity score script..."

## Run script to calculate modularity scores over each chunk
counter=1
last_line=$(wc -l < "$partitions")

## This loop reads lines in from the partition file to get chunk ranges
while read line; do
	wordarray=($line)
	range_start=${wordarray[0]}
	range_stop=${wordarray[1]}
	buffer_stop=${wordarray[2]}
	
	## If it is the final chunk then backtracking should be performed and the -f flag is needed
	if [ ${counter} -eq ${last_line} ]
	then
		${code_01} -i ${input_file} -o ${output_dir} -m ${model} -l ${max_length} -s ${range_start} -e ${range_stop} -b ${buffer_stop} -f
	
	## Otherwise the -f flag is not needed yet
	else
		${code_01} -i ${input_file} -o ${output_dir} -m ${model} -l ${max_length} -s ${range_start} -e ${range_stop} -b ${buffer_stop}
		echo "Range ${range_start} to ${range_stop} complete" 
		echo " "
		counter=$((counter+1))
	fi
 
done < $partitions

