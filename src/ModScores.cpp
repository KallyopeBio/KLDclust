/*
 * ModTerms.cpp
 *
 *  Created on: Sep 17, 2018
 *      Author: Sarah Christensen
 */
#include "ModScores.h"
#include "sys/types.h"

// Instantiate new ModScores instance with dictionary of penalties
ModScores::ModScores(std::string fof, const std::vector<unsigned int>& chromVect,
		const Penalty& penalty, const unsigned int maxLength,
		const unsigned int bufferStartBP, const unsigned int bufferEndBP):
m_chromVect(chromVect), m_maxLength(maxLength), m_penalty(penalty), m_bufferStartBP(bufferStartBP), m_bufferEndBP(bufferEndBP)
{
	// Get chromosome as subset to desired range
	Chromosome chrom(chromVect);
	m_varCount=m_chromVect.size();
	m_dictS.resize(m_varCount);

	// Get map to convert BP locations to indices
	chrom.fillBP2IndexMap(m_BP2IndexMap);

	// Calculate S Dictionary

	// Get vector containing files (will just have single file if a file was provided and not a folder)
	std::vector<std::string> fileNames(FileReader::get_files(fof));

	std::cout << "There are " << fileNames.size() << " files that will be read as input to calculate S" << std::endl;

	// Fill out dictionary S in parallel
	// Each file has all edges with loc a so no dictionary in vector can be written to at the same time
	#pragma omp parallel for schedule(static)
	for (unsigned int i = 0; i < fileNames.size(); i++){

		// Get input string
		std::string threadFile(fileNames[i]);
		FileReader fileStream(threadFile);
		ModScores::calculateS(fileStream, penalty);
	}
}

//= Get dictionary containing S score for each edge in the input data between specified range
void ModScores::calculateS(FileReader& fileStream, const Penalty& penalty)
{
	// Initialize dummy variables to store lines
	unsigned int loc_a(0);
	unsigned int loc_b(0);
	double r_score(0);

	// Loop until the next read is the end of the file
	while(fileStream.getNextRead(loc_a, loc_b, r_score))
	{
		// Only keep edge if it is in relevant range
		if ((loc_a >= m_bufferStartBP) && (loc_b <= m_bufferEndBP)){
			// Initialize start and end points
			unsigned int startIndex;
			unsigned int endIndex;

			// Calculate value
			double value = r_score - (penalty.getPenalty(loc_a, loc_b));

			// Order matters in m_DictS. It is indexed [startLoc][endLoc] -> value
			// Order is cleaned in file reader so we know loc_a <= loc_b
			startIndex=m_BP2IndexMap.at(loc_a);
			endIndex=m_BP2IndexMap.at(loc_b);

			// Try to add to dictionary
			std::pair< std::unordered_map<unsigned int, double>::iterator, bool > f=
					m_dictS[startIndex].insert(std::pair<unsigned int, double>(endIndex, value));

			// If the insertion failed it is because the key was already present
			if (!f.second) {

				// We only allow one correlation between each location
				throw std::runtime_error("ERROR: Input file contained two different values for correlation between same locations on chromosome.");
			}
		}
	}
}

//= Calculate modularity scores and add them to queue
bool ModScores::getScores(std::vector<std::unordered_map<unsigned int, double>>& dictC)
{
	// State flag
	bool bError=false;

	// Resize input array of dictionaries to proper size
	dictC.resize(m_varCount);

	// Get maximum number of threads available to do work and amount of work each thread should do
	unsigned int maxThreads = omp_get_max_threads();

	std::cout << "Running new parallel job with " << maxThreads << " threads" << std::endl;

	// BASE CASE (window size 0)

	// Length 0: Since self-loops removed, all r scores are 0
	#pragma omp parallel for schedule(static) num_threads(maxThreads)
	for (unsigned int i = 0; i < m_varCount; ++i){
		unsigned int loc = m_chromVect.at(i);
		dictC[i][i] = -(m_penalty.getPenalty(loc, loc));
	}

	// ITERATE OVER WINDOW OF SIZE > 0

	// Initialize window size
	unsigned int windowSize=1;

	// Initialize counter to keep track of how many variants have a window up to size x (in terms of number of variants)
	unsigned int count(1);

	//  Once counter remains 0 for full iteration then all variants have reached the max window size
	while (count>0){

		// Reset count of active variants for this iteration
		count = 0;

		// We do the following loop twice-- once over even and then once over odd indices
		// We do this to avoid needing to use locks to keep the parallel implementation thread safe
		// Even indices only ever look up values from odd indices in dictC and vice versa
		for (unsigned int j = 0; j <=1; j++) {

			// Create iterator over chromosome and use Open MP to parallelize
			#pragma omp parallel for schedule(static) num_threads(maxThreads) reduction(+:count)
			for (unsigned int i = j; i < m_varCount; i=i+2)
			{

				// Check if ending location past the end of the chromosome
				if ((i+windowSize)<m_varCount){

					// Get first chromosome location
					unsigned int communityStart = m_chromVect.at(i);

					// Get potential ending location
					unsigned int communityEnd = m_chromVect.at(i+windowSize);

					// Community is valid only if the length in BP is less than the max allowed
					if (communityEnd-communityStart <= m_maxLength){

						// Formula for prior column
						double a = dictC.at(i).at(i+windowSize-1);

						// Formula for next row
						double b = dictC.at(i+1).at(i+windowSize);

						// Formula for diagonal (if undefined contributes 0)
						double c(0);
						if (m_chromVect.at(i+1) <= m_chromVect.at(i+windowSize-1))
							c = dictC.at(i+1).at(i+windowSize-1);

						// Formula for new pair
						double d;
						// Check if the pair was in input and score already calculated
						if (m_dictS.at(i).count(i+windowSize)>0) {
							d = m_dictS.at(i).at(i+windowSize);
						}
						else{
							// Otherwise the pair was not in input and has an implicit correlation of 0
							d = -m_penalty.getPenalty(communityStart, communityEnd);
						}

						// Final score calculated as follows and added into modularity score dictionary
						dictC[i][i+windowSize]=a+b-c+d;

						// Increment count because this starting location could potentially be part of an even larger window
						++count;
					}
				}
			}
		}

		// Increment window size in terms of number of variants up by 1
		++windowSize;
	}

	return !bError;
}
