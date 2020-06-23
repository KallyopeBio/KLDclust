/*
 * 00_calculateCorrelation.cpp
 *
 *  Created on: Nov 13, 2018
 *      Author: sarah
 */

#include "FileReader.h"
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <unordered_map>
#include <assert.h>
#include <ctime>
#include <limits>
#include <sys/stat.h>
#include "FileReader.h"
#include <omp.h>
#include <boost/functional/hash.hpp>


struct pairhash {
public:
  std::size_t operator()(const std::pair<unsigned int, unsigned int> &x) const
  {
    return std::hash<unsigned int>{}(x.first) + std::hash<unsigned int>{}(x.second);
  }
};

static void show_usage(std::string name)
{
    std::cerr << "Usage: " << name << std::endl;
    std::cerr << "\t-i/--input INPUT FILES OR FOLDERS (required) in a comma separated list each with three columns of the form: LOC_A LOC_B CORRELATION." << std::endl;
    std::cerr << "\t-o/--output OUTPUT FOLDER (required) where output and intermediate results can be stored." << std::endl;
    return;
}

int main(int argc, char** argv) {

	// Print usage message if the proper number of arguments have not been provided
	if (argc < 4) {
		show_usage(argv[0]);
		return 1;
	}

	// Parse input parameters
	std::vector<std::string> input;
	std::string outputFolder;

	for (int i = 1; i < argc; ++i) {
		std::string arg = argv[i];
		if ((arg == "-h") || (arg == "--help")) {
		    show_usage(argv[0]);
		    return 0;
		}
		else if ((arg == "-i") || (arg == "--input")) {
			// Make sure argument was supplied for this flag
			if (i + 1 < argc) {
				std::stringstream tmp(argv[++i]);

				while( tmp.good() )
				{
				    std::string substr;
				    std::getline( tmp, substr, ',' );
				    input.push_back( substr );
				}
			}
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
		else {
			std::cerr << "Flag provided is not recognized. Please type -h or --help for details on input format." << std::endl;
		    return 1;
		}
	}

	// Print parsed input values
	std::cout << "Welcome to LD Clustering: Preprocessing Step to Combine Correlation Values!" << std::endl;
	std::cout << "Reading from the following input path(s): " << std::endl;
	for (unsigned int i=0; i<input.size(); i++){
		std::cout << "\t" << input[i] << std::endl;
	}
	std::cout << "Output will be written into the following folder: " << outputFolder << std::endl;

	// Initialize timers
	time_t timer1;
	time_t timer2;

	// Check if input is files or folders
	std::vector<std::vector<std::string>> inputFiles;
	inputFiles.resize(input.size());

	for (unsigned int i=0; i<input.size(); i++){
		inputFiles[i]=FileReader::get_files(input[i]);
	}

	// Check that the number of files in each folder was equal
	unsigned int numFiles = inputFiles[0].size();
	for (unsigned int i=0; i<inputFiles.size(); i++){
		if (inputFiles[i].size() != numFiles){
			throw std::runtime_error("There are not an equal number of files in each folder");
		}
	}

	// Create an array of output streams
	unsigned int outThreads(std::ceil(numFiles));
	std::vector<std::ofstream> stream;
	stream.resize(outThreads);
	for (unsigned int i = 0; i < outThreads; ++i){
		std::string outFile(outputFolder+std::to_string(i)+".txt");
		stream[i].open(outFile);
		assert(stream[i].good());
	}

	// Loop over each chunk
	std::cout << "\nProcessing each chunk in parallel" << std::endl;
	time(&timer1);
	#pragma omp parallel for schedule(static) num_threads(outThreads)
	for (unsigned int j=0; j<numFiles; j++){

		// Get thread number
		unsigned int t = omp_get_thread_num();

		// Create dictionary for just this chunk (will avoid fitting all correlations in memory at once)
		std::unordered_map<std::pair<unsigned int,unsigned int>, double, pairhash> corrDict;

		// Process this chunk from each correlation folder
		for (unsigned int i=0; i<inputFiles.size(); i++){

			FileReader fileStream(inputFiles[i][j]);
			unsigned int loc_a(0);
			unsigned int loc_b(0);
			double r_score(0);

			// Loop until the next read is the end of the file
			while(fileStream.getNextRead(loc_a, loc_b, r_score)){

				// Initialize this crazy type which is what will be returned by dictionary.insert()
				std::pair< std::unordered_map<std::pair<unsigned int,unsigned int>, double, pairhash>::iterator, bool > f;

				// Try to add correlation to the dictionary
				f = corrDict.insert(std::make_pair(std::make_pair(loc_a, loc_b), r_score));

				// Check if the insertion failed because the key was already present
				if (!f.second){
					// If already present, keep max correlation
					double tmp(std::max(f.first->second, r_score));
					f.first->second = tmp;
				}
			}
		}

		// Write dictionary to file
		for(auto& it : corrDict){
			// Write this entry to file
			stream[t] << std::get<0>(it.first) << " " << std::get<1>(it.first) << " " << it.second << std::endl;
		}
	}

	time(&timer2);
	std::cout << "Processing chunks complete and took " << difftime(timer2,timer1) << " seconds" << std::endl;

	// Close streams
	for (unsigned int i = 0; i < numFiles; ++i){
		assert(stream[i].good());
		stream[i].close();}

	return 0;

}
