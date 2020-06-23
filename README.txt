1) To recompile the source code do the following:

	cd ./src/00_calculateCorrelation/
	make clean
	make

	cd ./src/00_splitChromosome/
	make clean
	make

	cd ./src/01_modScores/
	make clean 
	make

2) To run the C++ code on 1000 Genomes Phase 3 data, use runPhase3.sh. 
This code will automatically download the individual ID numbers and VCF data from the FTP site. 
It loops over the provided chromosome number range and superpopulation code. 

For each chromosome and superpopulation...
	a) It first calls 00_calculateCorrelation to compute the correlation input data and chi square penalty. 
	b) It then calls runKLDclust.sh which computes the clustering on the correlation files. 
	c) It finally can gzip and copy the correlation files to be backed up on AWS.

Some values that you can set...
        max_length: Max length of an LD block in BP (suggested value 5000000)
        max_vars: Max number of variants in a chunk for parallelism and should be based on RAM
                  Does not influence output unless too small to be possible (it will fail with relevant message)
	model: This refers to what penalty term should be used. The options are modularity, background, chi, or constant. 
		Modularity will use modularity score, background will calculate correlations longer than 2,000,000 BP, chi
		uses the chi square penalty at 5%, and constant will just allow you to manually provide a penalty. 


Some things you may need to configure or install to be able to run this code:

	To use AWS, you need to configure your AWS profile prior to running the code.

	You may want to install pigz which does gzip in parallel. 
	https://zlib.net/pigz/
	sudo yum install pigz

	You may need to increase your ulimit in linux because depending on how many threads you use, you may have many files open at once.
	I found the instructions on the following page to work well for me on AWS
	https://gist.github.com/somandubey/52bff8c7cc8639292629



3) The script that actually runs the clustering code is runKLDclust.sh. Right now, it takes 3 arguments as input from the command line:

	input_file or input_folder: An input correlation file with three columns (may be space or tab deliminated)
		Column 1: First variant
		Column 2: Second variant
		Column 3: Correlation between first and second variants
	*OR* a folder containing such files. 
	Note that providing a folder containing correlation files has the benefit that they can be be read in parallel.
	In order for this to be thread safe, the left endpoint (i.e. the smaller variant location in BP) must have all of its edges in the same file. 
	The provided code already does this when it hashes out the edges in the correlation calculation script, but something to be aware of if you try to provide your own folder of files generated elsewhere.    

	output_dir: A path to a directory that can store intermediate steps as well as final output. If the Chi Square or constant penalty is to be used, there should be a file in this folder called penalty.txt containing that penalty.

	model: See above for discussion of the model options.  

	There are additional hardcoded parameters in the script (e.g. chunk size) that you may want to modify but do not change with each run of the Phase 3 script so they do not get input from the command line. 
