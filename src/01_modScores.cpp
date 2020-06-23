//============================================================================
// Name        : LD_Clustering.cpp
// Author      : Sarah Christensen
// Version     :
// Copyright   : Your copyright notice
// Description : Perform LD Clustering on data
//============================================================================

#include "FileReader.h"
#include "DictionaryK.h"
#include "Chromosome.h"
#include "ModScores.h"
#include "Partition.h"
#include "./Penalty_Class/Penalty.h"
#include "./Penalty_Class/ThresholdPenalty.h"
#include "./Penalty_Class/DictionaryPenalty.h"
#include <assert.h>
#include <ctime>
#include <cstring>
#include <climits>

static void show_usage(std::string name)
{
    std::cerr << "Usage: " << name << std::endl;
    std::cerr << "\t-i/--input INPUT FILE OF FOLDER (required) files with three tab or space separated columns and no header (Location_1 Location_2 Correlation)." << std::endl;
    std::cerr << "\t-o/--output OUTPUT FOLDER (required) where output and intermediate results can be stored." << std::endl;
    std::cerr << "\t-m/--model NULL HYOPTHESIS MODEL (required) used to normalize correlation values." << std::endl;
    std::cerr << "\t-l/--length MAXIMUM LENGTH (required) of LD block given in base pair units (bp)." << std::endl;
    std::cerr << "\t-s/--start STARTING LOCATION (optional) in base pair units(bp) for parallel job." << std::endl;
    std::cerr << "\t-e/--end ENDING LOCATION (optional) in base pair units(bp) for parallel job." << std::endl;
    std::cerr << "\t-b/--buffer BUFFER ENDING LOCATION (optional) in base pair units(bp) for parallel job." << std::endl;
    std::cerr << "\t-f/--final FINAL CHUNK (optional) add this flag to perform backtracking if it is the final (or only) chunk." << std::endl;
    return;
}

int main(int argc, char** argv) {

	// Print usage message if the proper number of arguments have not been provided
	if (argc < 8) {
		show_usage(argv[0]);
		return 1;
	}

	// Parse input parameters
	std::string filePath;
	std::string outputFolder;
	std::string model;
	unsigned int maxLength(0);
	unsigned int startBP(0);
	unsigned int endBP(UINT_MAX);
	unsigned int bufferEndBP(UINT_MAX);
	bool final(false);

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
		else if ((arg == "-s") || (arg == "--start")) {
			if (i + 1 < argc) {startBP = std::stoul(argv[++i]);}
		    else {std::cerr << "Check that arguments were provided for all flag options." << std::endl;
		    	return 1;
		    }
		}
		else if ((arg == "-e") || (arg == "--end")) {
			if (i + 1 < argc) {endBP = std::stoul(argv[++i]);}
		    else {std::cerr << "Check that arguments were provided for all flag options." << std::endl;
		    	return 1;
		    }
		}
		else if ((arg == "-b") || (arg == "--buffer")) {
			if (i + 1 < argc) {bufferEndBP = std::stoul(argv[++i]);}
		    else {std::cerr << "Check that arguments were provided for all flag options." << std::endl;
		    	return 1;
		    }
		}
		else if ((arg == "-f") || (arg == "--final")) {final=true;}
		else {
			std::cerr << "Flag provided is not recognized. Please type -h or --help for details on input format." << std::endl;
		    return 1;
		}
	}

	// Print parsed input values
	std::cout << "Welcome to LD Clustering: Step 2 Modularity Score Calculation!" << std::endl;
	std::cout << "Reading input from the following path: " << filePath << std::endl;
	std::cout << "Output will be written into the following folder: " << outputFolder << std::endl;
	std::cout << "\t The maximum allowed window (in BP) is " << maxLength << std::endl;
	std::cout << "\t The chromosome starts at position " << startBP << std::endl;
	std::cout << "\t The chromosome ends at position " << endBP << std::endl;
	std::cout << "\t The buffer ends at position " << bufferEndBP << std::endl;

	// Initialize timers
	time_t timer1;
	time_t timer2;


	// DYNAMIC PROGRAMMING I: Get community modularity scores

	// Items we will continue to need in DP II
	// This vector is intentionally made with new so that it will be allocated from the heap
	// Since the dictionaries are used throughout the entire program, we do not delete the vector, which would otherwise be a bottleneck
	// This is much faster than slow, sequential deletion from the stack
	std::vector<std::unordered_map<unsigned int, double>>* dictC;
	dictC= new std::vector<std::unordered_map<unsigned int, double>>;

	// Items in following scope only needed in DP I
	{
	// Get chromosome on relevant subset
	Chromosome chrom = Chromosome::parsedFile(outputFolder+"variants.txt");
	chrom.subsetChromosome(startBP, bufferEndBP);
	std::vector<unsigned int> chromVect(chrom.getVector());

	// Set values in case was not set in input and relying on default values of 0 and MAX
	startBP=chrom.getStart();
	bufferEndBP=chrom.getEnd();
	endBP=std::min(endBP, bufferEndBP);

	if (model=="modularity"){

		// Get Dictionary K on relevant subset
		DictionaryK K(outputFolder+"normDictK.txt", outputFolder+"m.txt", startBP, bufferEndBP);
		std::unordered_map<unsigned int, double> dictNormK(K.getK());

		// Now do community modularity score calculation
		std::cout << "Beginning calculation of modularity scores" << std::endl;
		time(&timer1);
		Chromosome chrom = Chromosome::parsedFile(outputFolder+"variants.txt");
		double coeff(chrom.getLengthBP());
		coeff /= maxLength;
		DictionaryPenalty penalty(dictNormK, coeff);
		ModScores mScores(filePath, chromVect, penalty, maxLength, startBP, bufferEndBP);
		time(&timer2);
		double seconds = difftime(timer2,timer1);
		std::cout << "Calculation of dictionary S complete and took " << seconds << " seconds" << std::endl;
		time(&timer1);
		mScores.getScores(*dictC);
		time(&timer2);
		seconds = difftime(timer2,timer1);
		std::cout << "Calculation of all modularity scores complete and took " << seconds << " seconds" << std::endl;

	}

	else if (model=="chi" || model=="background" || model=="constant"){

		// Get penalty term
		std::cout << "Reading in threshold penalty" << std::endl;
		double threshold(0);
		std::ifstream m_in(outputFolder+"penalty.txt");
		m_in >> threshold;
		m_in.close();
		std::cout << "Threshold penalty: " << threshold << std::endl;

		// Now do community modularity score calculation
		std::cout << "Beginning calculation of threshold scores" << std::endl;
		time(&timer1);
		ThresholdPenalty penalty(threshold);
		ModScores mScores(filePath, chromVect, penalty, maxLength, startBP, bufferEndBP);
		time(&timer2);
		double seconds = difftime(timer2,timer1);
		std::cout << "Calculation of dictionary S complete and took " << seconds << " seconds" << std::endl;
		time(&timer1);
		mScores.getScores(*dictC);
		time(&timer2);
		seconds = difftime(timer2,timer1);
		std::cout << "Calculation of all threshold scores complete and took " << seconds << " seconds" << std::endl;
	}

	else {
		throw std::runtime_error("Model type is not recognized");
	}
	// Remove inaccurate items (buffer indices) from DictC
	for (int i=((*dictC).size()-1); i>=0; i--){
		if (chromVect.at(i)>endBP){
			(*dictC).pop_back();
		}
		else{
			break;
		}
	}

	}

	// DYNAMIC PROGRAMMING II: Find best modularity score
	std::cout << "Beginning calculation of optimal modularity score" << std::endl;
	Chromosome chrom = Chromosome::parsedFile(outputFolder+"variants.txt");
	time(&timer1);
	Partition breaks(chrom, maxLength, startBP, endBP, *dictC, outputFolder);
	time(&timer2);
	double seconds = difftime(timer2,timer1);
	std::cout << "Calculation of optimal modularity score complete and took " << seconds << " seconds" << std::endl;

	// Save to file
	time(&timer1);
	breaks.printScores(outputFolder+"allOptScores.txt");
	breaks.printBreaks(outputFolder+"allOptBreaks.txt");
	time(&timer2);
	seconds = difftime(timer2,timer1);
	std::cout << "Intermediate results saved to file and took " << seconds << " seconds" << std::endl;

	// BACKTRACKING: Return a partition that achieves optimal modularity score
	if (final){

		std::cout << "\nBeginning backtracking" << std::endl;

		double optScore(breaks.getOpt());
		std::ofstream ofm(outputFolder+"FINAL_optScore.txt");
		ofm << optScore << std::endl;
		ofm.close();
		std::cout << "Optimal score: " << optScore << std::endl;

		std::vector<std::pair<unsigned int, unsigned int>> optPartition;
		if (!breaks.getOptPartition(optPartition)) {return 1;}

		// Save optimal partition
		int previousEnd(chrom.getStart());
		int totalUncovered(0);
		unsigned int maxGap(0);
		std::ofstream ofs(outputFolder+"FINAL_optBreaks.txt");
		for(auto it = optPartition.begin(); it != optPartition.end(); ++it){
			// Write this entry to file
			ofs << std::get<0>(*it) << " " << std::get<1>(*it) << std::endl;

			// Add length of Chromosome that is not covered (in BP)
			totalUncovered = totalUncovered + (std::get<0>(*it)-previousEnd-1);

			// Find max length of gap
			if (std::get<0>(*it)-previousEnd > maxGap){
				maxGap = std::get<0>(*it)-previousEnd;
			}

			// Update endpoint
			previousEnd = std::get<1>(*it);
		}
		// Close file
		ofs.close();

		std::cout << "Number of intervals in the optimal partition: " << optPartition.size() << std::endl;
		std::cout << "The chromosome range in this data starts at " << chrom.getStart() << " and ends at " << chrom.getEnd() << std::endl;
		std::cout << "This translates to a length of " << chrom.getEnd()-chrom.getStart() << " BP" << std::endl;
		std::cout << "There are " << totalUncovered << " BP in this range not assigned to an LD block, which is " << ((totalUncovered*100.0)/(chrom.getLengthBP()*1.0)) << " percent of the chromosome" << std::endl;
		std::cout << "The maximum gap is " << maxGap-1 << " BP" << std::endl;
		std::cout << "\nLD Clustering complete!" << std::endl;
		std::cout << "=^_^=" << std::endl;

	}

	return 0;

}

