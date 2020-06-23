/*
 * 00_splitChromosome.cpp
 *
 *  Created on: Oct 17, 2018
 *      Author: sarah
 */

#include "FileReader.h"
#include "DictionaryK.h"
#include "Chromosome.h"
#include "./Penalty_Class/Penalty.h"
#include "./Penalty_Class/ThresholdPenalty.h"
#include "./Penalty_Class/DictionaryPenalty.h"
#include <assert.h>
#include <ctime>
#include <cmath>
#include <iostream>
#include <math.h>
#include <string>
#include <vector>
#include <unordered_set>
#include <algorithm>

static void show_usage(std::string name)
{
    std::cerr << "Usage: " << name << std::endl;
    std::cerr << "\t-i/--input INPUT FILE OR FOLDER (required) with three tab or space separated columns and no header (Location_1 Location_2 Correlation)." << std::endl;
    std::cerr << "\t-o/--output OUTPUT FOLDER (required) where output and intermediate results can be stored." << std::endl;
    std::cerr << "\t-m/--model NULL HYOPTHESIS MODEL (required) used to normalize correlation values." << std::endl;
    std::cerr << "\t-l/--length MAXIMUM LENGTH (required) of LD block given in base pair units (bp)." << std::endl;
    std::cerr << "\t-v/--variants (APPROXIMATE) MAXIMUM VARIANTS (required) allowed in a chunk for processing efficiency/RAM management." << std::endl;
    return;
}

bool file_exists(std::string fileName)
{
    return FileReader::is_file(fileName);
}

int main(int argc, char** argv) {

	// Print usage message if the proper number of arguments have not been provided
	if (argc < 10) {
		show_usage(argv[0]);
		return 1;
	}

	// Parse input parameters
	std::string filePath;
	std::string outputFolder;
	std::string model;
	unsigned int maxLength(0);
	unsigned int maxVars(0);

	for (int i = 1; i < argc; ++i) {
		std::string arg = argv[i];
	    if ((arg == "-h") || (arg == "--help")) {
	    	show_usage(argv[0]);
	            return 0;
	        }
	    else if ((arg == "-i") || (arg == "--input")) {
	    	// Make sure argument was supplied for this flag
	        if (i + 1 < argc) {filePath = argv[++i];}
	        else { // Otherwise no argument was provided with this option
	                  std::cerr << "Check that arguments were provided for all flag options." << std::endl;
	                return 1;
	        }
	    }
	    else if ((arg == "-o") || (arg == "--output")) {
	    	if (i + 1 < argc) {outputFolder = argv[++i];}
	    		else {std::cerr << "Check that arguments were provided for all flag options." << std::endl;
	    	    	return 1;
	    	    }
	    }
	    else if ((arg == "-m") || (arg == "--model")) {
	    	    	if (i + 1 < argc) {model = argv[++i];}
	    	    		else {std::cerr << "Check that arguments were provided for all flag options." << std::endl;
	    	    	    	return 1;
	    	    	    }
	    }
	    else if ((arg == "-l") || (arg == "--length")) {
	    	if (i + 1 < argc) {
	    		maxLength = std::stoul(argv[++i]);}
	    		else {std::cerr << "Check that arguments were provided for all flag options." << std::endl;
	    	    	return 1;
	    	    }
	    }
	    else if ((arg == "-v") || (arg == "--variants") || (arg == "--variant")) {
	    	if (i + 1 < argc) {maxVars = std::stoul(argv[++i]);}
	    		else {std::cerr << "Check that arguments were provided for all flag options." << std::endl;
	    	    	return 1;
	    	    }
	    }
	    else {
	    	std::cerr << "Flag provided is not recognized. Please type -h or --help for details on input format." << std::endl;
	    	return 1;
	    }
	}

	// Print parsed input values
	std::cout << "\nWelcome to LD Clustering: Step 1 Chromosome Chunking and Preprocessing!" << std::endl;
	std::cout << "\nReading from the following file path: " << filePath << std::endl;
	std::cout << "Output will be written into the following folder: " << outputFolder << std::endl;
	std::cout << "Null hypothesis model selected: " << model << std::endl;
	std::cout << "Maximum length of LD block: " << maxLength << std::endl;
	std::cout << "(Approximate) Maximum number of variants in a chunk for processing: " << maxVars << std::endl;

	// Initialize timers
	time_t timer1;
	time_t timer2;

	// Report on long range correlations for reference
	std::cout << "\nBeginning calculation of long-range correlation summary statistics" << std::endl;
	time(&timer1);
	std::vector<std::string> inputFiles(FileReader::get_files(filePath));
	long double backgroundPenalty = ThresholdPenalty::background(inputFiles, 2000000);
	time(&timer2);
	std::cout << "\nCalculation of long-range correlation summary statistics complete and took " << difftime(timer2,timer1) << " seconds" << std::endl;
	std::cout << "\nTwo standard deviations above mean will be used as threshold penalty: " << backgroundPenalty << std::endl;

	// MODULARITY SCORE PREPROCESSING
	if (model=="modularity"){
		// Check if preprocessing files have not yet been calculated
		if (!file_exists(outputFolder+"variants.txt") || !file_exists(outputFolder+"m.txt") || !file_exists(outputFolder+"normDictK.txt")){

			// Create Dictionary K Normalized by sqrt 2m
			// Each key is a variant location (in BP) and each value is the degree of the corresponding variant node
			// To see how this is used, see definition of modularity score (i.e. terms Ki Kj)
			std::cout << "\nBeginning calculation of the dictionary K" << std::endl;
			time(&timer1);
			DictionaryK dictK(inputFiles);
			time(&timer2);
			std::cout << "Calculation of the normalized dictionary K complete and took " << difftime(timer2,timer1) << " seconds" << std::endl;
			std::cout << "\t Total weight (m) in implicit input graph is " << dictK.getM() << std::endl;

			// Create Chromosome object
			Chromosome chrom(dictK.getK());
			std::vector<unsigned int> chromVect = chrom.getVector();

			// Print items to output folder for later programs to use
			dictK.printK(outputFolder+"normDictK.txt");
			dictK.printM(outputFolder+"m.txt");
			chrom.printChrom(outputFolder+"variants.txt");
		}
	}

	// CONSTANT LONG-RANGE BACKGROUND THRESHOLD PREPROCESSING
	else if (model=="background"){

		if (!file_exists(outputFolder+"penalty.txt")){

			// Create vector of input files (whether from folder or file)
			std::cout << "\nBeginning calculation of long-range correlation summary statistics" << std::endl;
			time(&timer1);
			std::vector<std::string> inputFiles(FileReader::get_files(filePath));
			long double backgroundPenalty = ThresholdPenalty::background(inputFiles, 2000000);
			time(&timer2);

			std::cout << "\nCalculation of long-range correlation summary statistics complete and took " << difftime(timer2,timer1) << " seconds" << std::endl;
			std::cout << "\nTwo standard deviations above mean will be used as threshold penalty: " << backgroundPenalty << std::endl;
			std::ofstream ofm(outputFolder+"penalty.txt");
			ofm << backgroundPenalty << std::endl;
			ofm.close();
		}
	}

	// FOR ALL MODELS A VARIANT LIST IS REQUIRED (e.g. FISHER MODEL NO ADDITIONAL PREPROCESSING NEEDED HERE)
	if (!file_exists(outputFolder+"variants.txt")){
		// Create vector of input files (whether from folder or file)
		std::vector<std::string> inputFiles(FileReader::get_files(filePath));
		std::cout << "The number of files detected for variant list calculation is: " << inputFiles.size() << std::endl;

		// Create Chromosome object
		std::cout << "\nBeginning calculation of variant list" << std::endl;
		time(&timer1);
		Chromosome chrom = Chromosome::correlationFiles(inputFiles);
		chrom.printChrom(outputFolder+"variants.txt");
		time(&timer2);
		std::cout << "\nCalculation of variant set complete and took " << difftime(timer2,timer1) << " seconds" << std::endl;

	}

	// Reconstruct chromosome from file
	Chromosome chrom = Chromosome::parsedFile(outputFolder+"variants.txt");
	std::vector<unsigned int> chromVect(chrom.getVector());
	unsigned int numVar(chrom.getLengthCount());
	std::cout << "Number of variants in input file: " << numVar << std::endl;
	std::cout << "Chromosome Start: " << chrom.getStart() << " Chromosome End: " << chrom.getEnd() << std::endl;

	// Get intervals over this chromosome (startBP, endBP, bufferEndBP)
	std::cout << "\nBeginning calculation of chunks"<< std::endl;
	std::vector<std::tuple<unsigned int, unsigned int, unsigned int>> intervalVector(chrom.getIntervals(maxLength, maxVars));

	// Print intervals to file
	std::ofstream ofs(outputFolder+"partitions.txt");
	for (int i = (intervalVector.size()-1); i >= 0; --i){

		// Write this interval to output file with following form: startBP endBP bufferEndBP
		ofs << std::get<0>(intervalVector.at(i)) << " " << std::get<1>(intervalVector[i]) << " " << std::get<2>(intervalVector[i]) << std::endl;

	}
	ofs.close();

	return 0;
}


