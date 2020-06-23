/*
 * Chromosome.cpp
 *
 *  Created on: Sep 14, 2018
 *      Author: Sarah Christensen
 */

#include "Chromosome.h"

// Instantiate a new chromosome instance from file path (two possible types of files, correlation file and variant list file)
Chromosome::Chromosome(const std::string file, const bool calculate)
{
	// If the chromosome needs to be calculated from the input file
	if (calculate) {
		// Import file as stream
		FileReader fileStream(file);

		// Initialize dummy variables to store lines
		unsigned int loc_a(0);
		unsigned int loc_b(0);
		double r_score(0);

		// Initialize set to only keep unique variants in sorted order
		std::set<unsigned int> variantsSet;

		// Keep track of the number of lines being read
		unsigned int i=0;

		// Loop until the next read is the end of the file
		while(fileStream.getNextRead(loc_a, loc_b, r_score)){

			// Add first variant location to dictionary
			variantsSet.emplace(loc_a);

			// Add second variant location to dictionary
			variantsSet.emplace(loc_b);

			// Increment i counter
			++i;
		}

		if (i==0)
			std::cerr << "ERROR: No lines read from file when attempting to create chromosome" << std::endl;

		// Convert set to vector for constant time lookup
		m_variants.assign(variantsSet.begin(), variantsSet.end());
	}

	// Otherwise, the file contains an already processed chromosome so we can read it in directly
	else {
		unsigned int loc;
		std::ifstream in(file);
		while(true){
		  in >> loc;
		  if( in.eof() ) break;
		  m_variants.push_back(loc);
		}
		in.close();
	}

	// Set start and end location of this chromosome vector
	m_start = m_variants.front();
	m_end = m_variants.back();
}

// Instantiate a new chromosome instance from a file just containing a list of variants
Chromosome Chromosome::parsedFile(std::string filepath){
	return Chromosome(filepath, false);
}

// Instantiate a new chromosome instance from a correlation file with 3 columns
Chromosome Chromosome::correlationFile(std::string filepath){
	return Chromosome(filepath, true);
}

// Instantiate a new chromosome instance from vector of file paths
Chromosome Chromosome::correlationFiles(const std::vector<std::string> inputFiles){

	// Initialize values that will be reduced across threads
	std::vector<std::unordered_set<unsigned int>> variantSets;
	variantSets.resize(omp_get_max_threads());

	#pragma omp parallel for
	for (unsigned int i=0; i<inputFiles.size(); i++){

		// Get thread number
		unsigned int t = omp_get_thread_num();

		// Initialize values
		unsigned int loc_a(0);
		unsigned int loc_b(0);
		double r_score(0);

		FileReader fileStream(inputFiles[i]);
		while(fileStream.getNextRead(loc_a, loc_b, r_score)){
			// Add variant locations to set for this thread
			variantSets[t].insert(loc_a);
			variantSets[t].insert(loc_b);
		}
	}

	// Reduce variants sets into one dictionary
	std::unordered_map<unsigned int, double> variantDict;
	for (unsigned int i=0; i<variantSets.size(); i++){
		for (const auto& elem: variantSets[i]) {
			variantDict[elem] = 0;
		}
	}

	return Chromosome(variantDict);
}

// Instantiate a new chromosome instance from vector
Chromosome::Chromosome(const std::vector<unsigned int>& chromVect): m_variants(chromVect)
{
	// Set start and end location of this chromosome vector
	m_start = m_variants.front();
	m_end = m_variants.back();
}

// Instantiate a new chromosome instance from double dictionary
Chromosome::Chromosome(const std::unordered_map<unsigned int, double>& dictK)
{
	// Keep track of the number of lines being read
	unsigned int i=0;

	// Loop over dictionary
	for (auto it = dictK.begin(); it != dictK.end(); it++)
	{
		// Add key to variant list
		unsigned int key = it->first;
		m_variants.push_back(key);
		++i;
	}

	if (i==0)
	{
		std::cerr << "ERROR: No lines read when attempting to create chromosome" << std::endl;
	}

	// Sort vector
	std::sort (m_variants.begin(), m_variants.end());

	// Set start and end location of this chromosome vector
	m_start = m_variants.front();
	m_end = m_variants.back();
}

// Instantiate a new chromosome instance from unsigned int dictionary
Chromosome::Chromosome(const std::unordered_map<unsigned int, unsigned int>& dict)
{
	// Keep track of the number of lines being read
	unsigned int i=0;

	// Loop over dictionary
	for (auto it = dict.begin(); it != dict.end(); it++)
	{
		// Add key to variant list
		unsigned int key = it->first;
		m_variants.push_back(key);
		++i;
	}

	if (i==0)
	{
		std::cerr << "ERROR: No lines read from Dictionary K when attempting to create chromosome" << std::endl;
	}

	// Sort vector
	std::sort (m_variants.begin(), m_variants.end());

	// Set start and end location of this chromosome vector
	m_start = m_variants.front();
	m_end = m_variants.back();
}

//= Get chromosome variants returned as ordered set
std::vector<unsigned int> Chromosome::getVector(){return m_variants;}

// Fill dictionary with key (BP location) and value (index location)
void Chromosome::fillBP2IndexMap(std::unordered_map<unsigned int, unsigned int>  & BP2IndexMap){
	// Loop over vector of variants
	for(std::vector<unsigned int>::size_type i = 0; i != m_variants.size(); i++) {
		BP2IndexMap.insert(std::make_pair(m_variants[i], i));
	}
}

// Subset this chromosome to specified range
void Chromosome::subsetChromosome(const unsigned int subsetStart, const unsigned int subsetEnd){
	unsigned int numVar(m_variants.size());
	std::vector<unsigned int> tmp;

	for (unsigned int i = 0; i < numVar; i++) {
		if ((m_variants[i] >= subsetStart) && (m_variants[i] <= subsetEnd)){
			tmp.push_back(m_variants[i]);
		}
	}

	// Update chromosome variants to be this new vector
	m_variants = tmp;

	// Update chromosome start and end location
	m_start = m_variants.front();
	m_end = m_variants.back();
}

//= Get beginning of the chromosome
unsigned int Chromosome::getStart(){return m_start;}

//= Get end of the chromosome
unsigned int Chromosome::getEnd(){return m_end;}

//= Get length of the chromosome in BP
unsigned int Chromosome::getLengthBP(){
	return m_end-m_start;
}

//= Get length of the chromosome in number of variants
unsigned int Chromosome::getLengthCount(){
	return m_variants.size();
}

//= Get farthest variant within x BP away to the right from variant at some index location
unsigned int Chromosome::getRightNeighborIndex(const unsigned int locBP, const unsigned int distanceBP)
{
	// Largest possible chromosome location that could be included in a window with locBP
	unsigned int candidateUpper=locBP+distanceBP;

	// Returns an iterator pointing to first element in set greater than candidateUB (i.e. first element outside of neighbor window)
	auto itup = upper_bound(m_variants.cbegin(), m_variants.cend(), candidateUpper);

	// Need to actually find first element NOT greater than candidate UB
	unsigned int upperIndex;

	// .end() is an rvalue expression so decrementing it is not guaranteed to work. We handle this separately.
	if (itup==m_variants.cend()) {
		upperIndex = m_variants.size()-1;  // Because indexing starts at 0 we subtract 1
	}
	else
	{
		// Move left one in iterator to find first max eligible neighbor in set of variants
		itup--;
		upperIndex = itup-m_variants.cbegin();
	}

	// Return this rightmost neighbor as index location
	return upperIndex;
}

//= Get farthest variant within x BP away to the left from variant at some index location
unsigned int Chromosome::getLeftNeighborIndex(const unsigned int locBP, const unsigned int distanceBP)
{
	// Smallest possible chromosome location that could be included in a window with locBP
	int candidateLower=locBP-distanceBP;
	unsigned int absCandidateLower= abs(candidateLower);

	// Check if candidate is out of range
	if ((absCandidateLower<= m_start) || (candidateLower<0)){
		return 0;
	}

	// Otherwise
	else{
		// Lower_bound does a binary search for the first element <= the query.
		auto itlow = lower_bound(m_variants.cbegin(), m_variants.cend(), candidateLower);

		int index = itlow - m_variants.cbegin();

		// Return this leftmost neighbor as index location
		return index;
	}
}

//= Get farthest variant within x BP away to the right from variant at some BP location
unsigned int Chromosome::getRightNeighborBP(const unsigned int locBP, const unsigned int distanceBP)
{
	// Get this rightmost neighbor as index location
	int index = Chromosome::getRightNeighborIndex(locBP, distanceBP);

	// Return this leftmost neighbor as BP location
	return m_variants[index];
}

//= Get farthest variant within x BP away to the left from variant at some BP location
unsigned int Chromosome::getLeftNeighborBP(const unsigned int locBP, const unsigned int distanceBP)
{
	// Get this leftmost neighbor as index location
	int index = Chromosome::getLeftNeighborIndex(locBP, distanceBP);

	// Return this leftmost neighbor as BP location
	return m_variants[index];
}

//= Get a vector of tuples containing intervals with (approximately) at most maxVars variants and will be correct for LD blocks of at most maxLength
std::vector<std::tuple<unsigned int, unsigned int, unsigned int>> Chromosome::getIntervals(
		const unsigned int maxLength,
		const unsigned int maxVars)
{

	unsigned int numVar(m_variants.size());
	std::vector<std::tuple<unsigned int, unsigned int, unsigned int>> vectIntervalsBP;
	unsigned int i = 0;

	while (i < numVar) {
		// i is an index (e.g., the ith variant present in the data, as opposed to the variant in location i)
		// Initialize start and end in terms of variant count
		unsigned int startVar(i);
		unsigned int bufferEndVar(
				std::min(startVar + (maxVars - 1), numVar - 1));

		// Convert these variants to BP locations
		unsigned int startBP(m_variants.at(startVar));
		unsigned int bufferEndBP(m_variants.at(bufferEndVar));

		// Find end of accurate window in this range
		unsigned int endVar(0);
		// endVar will exclude the "buffer region".
		if (bufferEndVar == numVar - 1){ //If the buffer reaches the end then everything up to the buffer will be accurate
			endVar = bufferEndVar;
		}
		else{
			endVar = getLeftNeighborIndex(bufferEndBP, maxLength); // Otherwise, we make overlap
		}
		unsigned int endBP(m_variants.at(endVar));

		// Since get left and get right var may not be symmetric, we need to update the right buffer endpoint
		bufferEndVar = getRightNeighborIndex(endBP, maxLength);
		bufferEndBP = m_variants.at(bufferEndVar);

		std::cout << "\nStart Variant: " << startVar << " End Variant: "
				<< endVar << " Buffer End Variant: " << bufferEndVar
				<< std::endl;
		std::cout << "Start BP: " << startBP << " End BP: " << endBP
				<< " Buffer End Variant: " << bufferEndBP << std::endl;

		// Print the actual number of variants that will be in this range
		std::cout << "Actual chunk size (including overlap buffer): "
				<< bufferEndVar - startVar + 1 << std::endl;
		std::cout << "Accurate chunk size (excluding overlap buffer): "
				<< endVar - startVar + 1 << std::endl;

		// Check that range is valid (if not, then max var range must be increased to be feasible)
		if (startVar > endVar){
			throw std::runtime_error("The maximum variants per chunk must be increased to guarantee finding the correct optimal partition.");
		}

		// Write this entry to output file
		vectIntervalsBP.push_back(std::make_tuple(startBP, endBP, bufferEndBP));

		// Now i gets updated to be the next variant after the end of the accurate window
		i = endVar + 1;

	}

	// Return overlapping intervals covering the chromosome
	return vectIntervalsBP;
}


// Subset correlation file based on chromosome intervals
void Chromosome::subsetCorrelationFile(const std::string file,
		std::vector<std::tuple<unsigned int, unsigned int, unsigned int>> intervals,
		const std::string outputFolder, const unsigned int numOutFiles){

	// Stream over correlation input file
	#pragma omp parallel for schedule(static)
	for (unsigned int i=0; i<intervals.size(); i++){

		// Get start and End buffer for each interval
		unsigned int start = std::get<0>(intervals.at(i));
		unsigned int end = std::get<2>(intervals.at(i));

		// Make a file for each interval in output folder
		std::string tmp(outputFolder+std::to_string(start)+"_"+std::to_string(end));
		const char* intervalFolder(tmp.c_str());
		if (mkdir(intervalFolder, 0777) == -1){
			std::cerr << "ERROR: When creating directory  " << strerror(errno) << std::endl;
		}

		// Each thread will have a vector of streams
		std::vector<std::ofstream> outStreams;
		outStreams.resize(numOutFiles);

		// Open multiple streams per interval to speed up read in later
		for (unsigned int j=0; j<numOutFiles; j++){
			outStreams[j].open(tmp+"/"+std::to_string(j)+".txt");
		}

		// Initialize dummy variables to store lines
		unsigned int loc_a(0);
		unsigned int loc_b(0);
		double r_score(0);

		FileReader fileStream(file);
		while(fileStream.getNextRead(loc_a, loc_b, r_score)){


			// Place this input edge in the corresponding file
			if ((loc_a >= start) && (loc_b<= end)){

				// Take mod of starting index so that it is always in same file
				unsigned int j= loc_a % numOutFiles;
				outStreams[j] << loc_a << " " << loc_b << " " << r_score << std::endl;
			}
		}

		// Close all file streams
		for (unsigned int j=0; j<numOutFiles; j++){
			outStreams[j].close();
		}
	}
}

void Chromosome::printChrom(const std::string file){
	// Print chromosome locations to file
	std::ofstream ofp(file);
	for (std::vector<unsigned int>::iterator it = m_variants.begin(); it != m_variants.end(); ++it){
		// Write variant location to output file
		ofp << *it << std::endl;
	}
	// Close file
	ofp.close();
}


