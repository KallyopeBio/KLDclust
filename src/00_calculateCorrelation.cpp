/*
 * 00_calculateCorrelation.cpp
 *
 *  Created on: Nov 13, 2018
 *      Author: sarah
 */


#include "GeneticMap.h"
#include "VCF.h"
#include "Chromosome.h"
#include "./Penalty_Class/Penalty.h"
#include "./Penalty_Class/ThresholdPenalty.h"
#include "./Penalty_Class/DictionaryPenalty.h"
#include <assert.h>
#include <ctime>
#include <cstring>

static void show_usage(std::string name)
{
    std::cerr << "Usage: " << name << std::endl;
    std::cerr << "\t-i/--input INPUT VCF FILE (required)." << std::endl;
    std::cerr << "\t-o/--output OUTPUT FOLDER (required) where output and intermediate results can be stored." << std::endl;
    std::cerr << "\t-g/--group SAMPLE IDS OF INTEREST (required) file with list of sample IDs to be included (one per line)." << std::endl;
    std::cerr << "\t-l/--length MAXIMUM LENGTH (required) of LD block given in base pair units (bp)." << std::endl;
    std::cerr << "\t-m/--map GENETIC MAP (required for shrinkage correction to correlation) file containing interpolated genetic map positions (SNP, Location BP, Genetic Position)." << std::endl;
    std::cerr << "\t-s/--size EFFECTIVE POPULATION SIZE (required for shrinkage correction to correlation) for shrinkage correction." << std::endl;
    std::cerr << "\t-f/--maf MINIMUM ALLELE FREQUENCY (optional) threshold to print variant correlations to output file." << std::endl;
    std::cerr << "\t-r2/--r2 MINIMUM CORRELATION R2 VALUE (optional) threshold to print to output file." << std::endl;
    std::cerr << "\t-p/--phased PHASED VCF DATA FLAG (optional) add this flag if VCF data should be treated as phased." << std::endl;
    return;
}

int main(int argc, char** argv) {

	// Print usage message if the proper number of arguments have not been provided
	if (argc < 6) {
		show_usage(argv[0]);
		return 1;
	}

	// Parse input parameters
	std::string filePath;
	std::string outputFolder;
	std::string mapFile("missing");
	std::string indFile;
	unsigned int maxLength(0);
	unsigned int popSize(0);
	double minR2(0);
	double maf(0);
	bool phased(false);

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
		else if ((arg == "-m") || (arg == "--map")) {
			// Make sure argument was supplied for this flag
			if (i + 1 < argc) {mapFile = argv[++i];}
		    else { // Otherwise no argument was provided with this option
		    	std::cerr << "Check that arguments were provided for all flag options." << std::endl;
		    	return 1;
		    }
		}
		else if ((arg == "-g") || (arg == "--group")) {
			// Make sure argument was supplied for this flag
			if (i + 1 < argc) {indFile = argv[++i];}
			else { // Otherwise no argument was provided with this option
				std::cerr << "Check that arguments were provided for all flag options." << std::endl;
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
		else if ((arg == "-s") || (arg == "--size")) {
			if (i + 1 < argc) {popSize = std::stoul(argv[++i]);}
		    else {std::cerr << "Check that arguments were provided for all flag options." << std::endl;
		    	return 1;
		    }
		}
		else if ((arg == "-r2") || (arg == "--r2")) {
			if (i + 1 < argc) {minR2 = std::stod(argv[++i]);}
		    else {std::cerr << "Check that arguments were provided for all flag options." << std::endl;
		    	return 1;
		    }
		}
		else if ((arg == "-f") || (arg == "--maf")) {
			if (i + 1 < argc) {maf = std::stod(argv[++i]);}
		    else {std::cerr << "Check that arguments were provided for all flag options." << std::endl;
		    	return 1;
		    }
		}
		else if ((arg == "-p") || (arg == "--phased")) {phased=true;}
		else {
			std::cerr << "Flag provided is not recognized. Please type -h or --help for details on input format." << std::endl;
		    return 1;
		}
	}

	// Print parsed input values
	std::cout << "Welcome to LD Clustering: Preprocessing Step to Calculating Correlation Values!" << std::endl;
	std::cout << "Reading from the following VCF file path: " << filePath << std::endl;
	if ((mapFile != "missing") && (popSize != 0)) {
		std::cout << "Reading from the following Genetic Map file path: " << mapFile << std::endl;
	}
	std::cout << "Subsetting to individuals with IDs in the following file path: " << indFile << std::endl;
	std::cout << "Output will be written into the following folder: " << outputFolder << std::endl;
	std::cout << "\t The maximum allowed window (in BP) is " << maxLength << std::endl;
	if ((mapFile != "missing") && (popSize != 0)) {
		std::cout << "\t The effective population size is " << popSize << std::endl;
	}
	std::cout << "\t The minimum R2 threshold is " << minR2 << std::endl;
	std::cout << "\t The minimum minor allele frequency threshold is " << maf << std::endl;
	if (phased) {std::cout << "\t The input data will be treated as PHASED" << std::endl;}
	if (!phased) {std::cout << "\t The input data will be treated as UNPHASED" << std::endl;}

	// Initialize timers
	time_t timer1;
	time_t timer2;

	// Get the number of lines in the VCF file
	std::fstream in(filePath);
	int lines = 0;
	char endline_char = '\n';
	while (in.ignore(std::numeric_limits<std::streamsize>::max(), in.widen(endline_char))){++lines;}
	in.close();
	std::cout << "\nThere are " << lines << " lines in the VCF file" << std::endl;

	// Parse VCF file to initialize member variables
	// (VCF file, estimated number of variants, max window length, phased flag)
	std::cout << "\nParsing VCF file" << std::endl;
	time(&timer1);
	VCF x(filePath, lines, maxLength, phased, indFile);
	time(&timer2);
	std::cout << "Parsing of VCF file complete and took " << difftime(timer2,timer1) << " seconds" << std::endl;

	// Get critical value for this population size
	unsigned int obsCount(x.getObsCount(phased));
	long double chiThreshold(ThresholdPenalty::chiSquare(obsCount, 0.05));
	std::ofstream ofm(outputFolder+"penalty.txt");
	ofm << chiThreshold << std::endl;
	ofm.close();
	std::cout << "The Chi Squared Test Threshold at 5% Significance is " << chiThreshold << std::endl;

	// Calculate correlations and write result to files (1 per thread)
	time(&timer1);
	if ((mapFile != "missing") && (popSize != 0)){

		std::cout << "\nCalculating shrinkage correlations" << std::endl;

		// Create Genetic Map dictionary
		GeneticMap gMap(mapFile);

		// Run shrinkage correction to correlation
		x.calculateShrinkageCorrelation(outputFolder, gMap, popSize, minR2);
	}
	else{
		std::cout << "\nCalculating correlations with MAF threshold of " << maf << std::endl;

		// Just filter out variants below maf threshold
		x.calculateMAFCorrelation(outputFolder, minR2, maf);
	}
	time(&timer2);
	std::cout << "Correlation calculations complete and took " << difftime(timer2,timer1) << " seconds \n" << std::endl;

	return 0;

}

