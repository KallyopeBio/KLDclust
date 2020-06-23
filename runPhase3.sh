#!/bin/bash

## Set input folder
input_folder="/home/ec2-user/input/"

## Download file of individual IDs from FTP 
mkdir -p "${input_folder}raw/"
cd "${input_folder}raw/"
wget -q "ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/integrated_call_samples_v3.20130502.ALL.panel" 

## Create lists of individuals by population (we remove the first line which is a header)
rm -rf "${input_folder}pop/"
mkdir -p "${input_folder}pop/"
panel_file="${input_folder}raw/integrated_call_samples_v3.20130502.ALL.panel"
tail --lines=+2 ${panel_file} > "${input_folder}pop/integrated_call_samples_v3.20130502.ALL.panel.clean"
while read line; do
	wordarray=($line)
	sample=${wordarray[0]}
	super_pop=${wordarray[2]}
        
	echo "${sample}" >> "${input_folder}pop/${super_pop}.txt"
	
done < "${input_folder}pop/integrated_call_samples_v3.20130502.ALL.panel.clean"

## Set paths to code binaries
code_dir="/home/ec2-user/Workspace"
code_00="${code_dir}/testingOpenMP/00_calculateCorrelation/00_calculateCorrelation"

## Set inputs such as input correlation file, max length of LD block, and max number of variants in a chunk for parallel processing
max_length="5000000"
minimum_r2="0.00001"
maf="0.05"
model="constant"

## For each chromosome we will run each population (including ALL)
for ((i=21; i<=21; i++)); do

	# Get input file 
	input_file="${input_folder}raw/ALL.chr${i}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf"
	
	if [ ! -e ${input_file} ]; then
		cd "${input_folder}raw/"
		wget -q ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr${i}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes*
	fi

	echo "Unzipping ${input_file}"
	gunzip "${input_file}.gz"
	echo "Done upzipping file"

	for pop_file in ${input_folder}pop/*.txt; do

		## Start timer
		start=`date +%s`

		## Get super population and create output folder for this chromosome and population 
		pop=$(basename "${pop_file}" .txt)

                ## Sample data randomly to same size
		sample_size=100
		pop_sample=${input_folder}pop/${pop}_${sample_size}.txt
                shuf -n ${sample_size} ${pop_file} > ${pop_sample}

		output_folder="/home/ec2-user/output/modularity_sample/chr${i}/${pop}/"
		correlation_folder="${output_folder}Correlations/"
		clustering_folder="${output_folder}Clustering/"
		mkdir -p ${output_folder}
		mkdir -p ${correlation_folder}
		mkdir -p ${clustering_folder}
		
		if [ -e ${correlation_folder}r2.txt ]
		then
			echo "Correlation script for Chromosome ${i} and ${pop} already calculated"

		else
			echo "Running correlation script for Chromosome ${i} and ${pop}"

			## Run script to preprocess input file and identify chunk ranges for parallel or sequential processing
			${code_00} -i ${input_file} -g ${pop_file} -l ${max_length} -f ${maf} -r2 ${minimum_r2} -o ${correlation_folder} -p

			echo "Correlation script for Chromosome ${i} and ${pop} done"
			mv ${correlation_folder}penalty.txt ${clustering_folder}penalty.txt
			cd ${correlation_folder}

			## Combine edges hashed to the same bucket by different threads into one file
			for ((j=0; j<96; j++)); do
				cat *_${j}.txt > ${j}.txt
				rm *_${j}.txt
			done
		fi

		## Now run clustering algorithm on output
		cd ${code_dir}
		./runKLDclust.sh ${correlation_folder} ${clustering_folder} ${model}

		echo "Chromosome ${i} ${pop} Duration: $((($(date +%s)-${start})/60)) minutes"

		## Zip the correlation folder
		cd ${output_folder}
		tar -I pigz -cf correlation_archive.tar.gz ${correlation_folder}

		## Upload max correlation file to S3 
		aws s3 cp "${output_folder}correlation_archive.tar.gz" "s3://ktmp/kldclust/${model}/chr${i}/${pop}/R2_MAF05.tar.gz" --quiet

		## Remove correlation file to save space
		rm ${output_folder}correlation_archive.tar.gz
		rm -r ${correlation_folder}

	done 

	## Remove input file
	rm "${input_file}"
done
