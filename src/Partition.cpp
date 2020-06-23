/*
 * Partition.cpp
 *
 *  Created on: Sep 19, 2018
 *      Author: Sarah Christensen
 */

#include "Partition.h"
#include <limits>
#include <algorithm>

// Instantiate new partition instance
Partition::Partition(Chromosome& chrom, const unsigned int maxLength,
		const unsigned int startBP, const unsigned int endBP,
		std::vector<std::unordered_map<unsigned int, double> >& dictC, std::string outputFolder):
m_chrom(chrom), m_maxLength(maxLength)
{
	// Load in opt scores and breaks that have already been computed (if any)
	std::ifstream ins(outputFolder+"allOptScores.txt");
	if (ins.good()){
		unsigned int loc;
		double val;
		while(true){
			ins >> loc;
			ins >> val;
			if( ins.eof() ) break;
			m_optScores[loc]=val;
		}
		ins.close();
	}

	std::ifstream inb(outputFolder+"allOptBreaks.txt");
	if (inb.good()){
		unsigned int loc;
		double val;
		while(true){
			inb >> loc;
			inb >> val;
			if( inb.eof() ) break;
			m_optBreaks[loc]=val;
		}
		inb.close();
	}

	// Get chromosome as a vector
	m_chromVect = chrom.getVector();

	// Convert BP locations to indices in the DictC dictionary (which does not span the whole chromosome)
	m_startCIndex = find(m_chromVect.begin(), m_chromVect.end(), startBP) - m_chromVect.begin();
	int endCIndex = find(m_chromVect.begin(), m_chromVect.end(), endBP) - m_chromVect.begin();
	int varCount = endCIndex-m_startCIndex;

	// For each location in the chromosome (starting with the largest)
	for (int i = endCIndex; i >= m_startCIndex; i--){

		// Get location of index
		unsigned int initLoc = m_chromVect.at(i);

		// Get rightmost location index that can be in window with location of interest
		unsigned int rightBound = m_chrom.getRightNeighborIndex(initLoc, m_maxLength);

		// Add optimal partition and score to member variables for this location
		if (!calculatePartition(i, initLoc, rightBound, dictC)){
			throw std::runtime_error("Error when calculating partition at site " + std::to_string(i));
		}

		// Track progress through this phase (based on number of accurate variants in dictC)
		unsigned int cIndex = i-m_startCIndex;
		unsigned int x(round(cIndex * 10 / varCount));
		unsigned int y(round((cIndex+1) * 10 / varCount));
		if (x == 8 && y == 9) {
					std::cout << "[====                ]20%\r"<< std::endl;
				}
				else if (x == 6 && y == 7) {
					std::cout << "[========            ]40%\r"<< std::endl;
				}
				else if (x == 4 && y == 5) {
					std::cout << "[============        ]60%\r"<< std::endl;
				}
				else if (x == 2 && y == 3) {
					std::cout << "[================    ]80%\r"<< std::endl;
				}
				else if (x == 1 && y == 2) {
					std::cout << "[====================]100%" << std::endl;
				}
	}
}

//= Calculates optimal partition for a given index
bool Partition::calculatePartition(unsigned int initIndex, unsigned int initLoc, unsigned int rightBound, std::vector<std::unordered_map<unsigned int, double> >& dictC)
{
	// State flag
	bool bError=false;

	// Create vector that stores minimum scores and corresponding location for each thread
	unsigned int maxThreads = omp_get_max_threads();
	std::vector<double> threadBestScore;
	std::vector<unsigned int> threadBestBreakEnd;

	// Initialize best score vector to lowest possible score
	threadBestScore.resize(maxThreads);
	threadBestBreakEnd.resize(maxThreads);
	for (unsigned int i = 0; i < maxThreads; i++) {threadBestScore[i]=std::numeric_limits<double>::lowest();}

	// Get number of variants
	unsigned int varCount=m_chromVect.size();

	// Get init index converted to dictC which is on a subset of the chromosome
	unsigned int cInitIndex = initIndex-m_startCIndex;

	// Check all possible windows that could contain this index for the best score
	#pragma omp parallel for schedule(static) num_threads(maxThreads)
	for (unsigned int candidateIndex = initIndex; candidateIndex <= rightBound; ++candidateIndex)
	{
		// Get thread number
		unsigned int t = omp_get_thread_num();

		// Initialize candidate score
		double candidateScore = 0.0;

		// If window reaches the end of the chromosome
		if (candidateIndex == (varCount-1))
			try
			{
				// Score is just the window itself (need to reindex since dictC does not span whole chromosome!!)
				candidateScore = dictC[cInitIndex].at(candidateIndex-m_startCIndex);
			}
			catch (...)
			{
				bError = true;
				std::cerr << "ERROR: Community was not found in community score dictionary." << std::endl;
			}
		// Otherwise,
		else
			try
			{
				// The candidate score is the sum of the "recursive" call plus the new window
				// Note that we need to reindex dictC since it does not span the whole chromosome (initIndex and candidateIndex are based on whole chromosome)
				candidateScore = m_optScores.at(m_chromVect[candidateIndex+1]) + dictC[cInitIndex][candidateIndex-m_startCIndex];
			}
			catch (...)
			{
				bError = true;
				std::cerr << "ERROR: Either community was not found in community score dictionary or recursive partition not found." << std::endl;
			}

		// Compare candidate score to existing best score. Keep if it is bigger.
		if (candidateScore > threadBestScore[t])
		{
			// Keep best score
			threadBestScore[t] = candidateScore;

			// Keep best location (in BP)
			threadBestBreakEnd[t] = m_chromVect[candidateIndex];
		}
	}

	// Join parallel threads
	#pragma omp barrier

	// Reduce results from each parallel thread
	double bestScore = std::numeric_limits<double>::lowest();
	unsigned int bestBreakEnd = initLoc;
	for (unsigned int i = 0; i < maxThreads; i++)
	{
	    if (threadBestScore[i] > bestScore)
	    {
	        bestScore = threadBestScore[i];
	        bestBreakEnd = threadBestBreakEnd[i];
	    }
	}

	// Add best score and partition ending at index to dictionaries
	m_optScores[initLoc] = bestScore;
	m_optBreaks[initLoc] = bestBreakEnd;

	return !bError;
}

//= Return optimal modularity score
double Partition::getOpt()
{
	// Initialize opt score
	double optScore(0);

	// Find start of chromosome
	unsigned int x = m_chrom.getStart();

	// Get score corresponding to the start of the chromosome
	try {optScore = m_optScores.at(x);}
	catch (...)
	{
		std::cerr << "ERROR: Cannot locate chromosome location " << x << " in dictionary of optimal scores." << std::endl;
	}

	return optScore;
}

//= Return partition that achieves optimal modularity score
bool Partition::getOptPartition(std::vector<std::pair<unsigned int, unsigned int> >& optPartition)
{
	// Set state flag
	bool bError=false;

	// Start with community containing the leftmost variant
	unsigned int l = m_chrom.getStart();

	// Initialize l to be 0 (not a valid chromosome location)
	unsigned int r=0;

	// Get end of chromosome
	unsigned int chromEnd = m_chrom.getEnd();

	// Continue until the last location of the chromosome is reached
	while (l <= chromEnd)
	{
		// Get optimal right end point
		try {r = m_optBreaks.at(l);}
		catch (...)
		{
			bError = true;
			std::cerr << "ERROR: Cannot locate chromosome location " << l << " in dictionary of optimal partitions." << std::endl;
		}

		// Add to output
		optPartition.push_back(std::make_pair(l,r));

		// If right end of interval is end of chromosome then you are done (this is a bit redundant with the while loop)
		if (r==chromEnd) {break;}

		// Otherwise, find the left end point (one bigger than r) of the next window implied by break (log time)
		// If this is bottleneck, potentially can be done in constant time by different data structures
		std::vector<unsigned int>::iterator it;
		it = std::upper_bound (m_chromVect.begin(), m_chromVect.end(), r);
		l=*it;

	}
	return !bError;
}

// Print optimal scores to file
void Partition::printScores(const std::string file){
	std::ofstream of(file);
	for (std::unordered_map<unsigned int, double>::iterator it = m_optScores.begin(); it != m_optScores.end(); ++it)
		of << it->first << " " << it->second << std::endl;
	of.close();
}

// Print optimal breaks to file
void Partition::printBreaks(const std::string file){
	std::ofstream of(file);
	for (std::unordered_map<unsigned int, unsigned int>::iterator it = m_optBreaks.begin(); it != m_optBreaks.end(); ++it)
		of << it->first << " " << it->second << std::endl;
	of.close();
}
