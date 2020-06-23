/*
 * ModScores.h
 *
 *  Created on: Sep 14, 2018
 *      Author: Sarah Christensen
 * Description: Class takes as input a PLINK file, chromosome, and K dictionary. Returns a dictionary containing modularity scores for all possible communities.
 */

#ifndef MODSCORES_H_
#define MODSCORES_H_

#include "Chromosome.h"
#include "DictionaryK.h"
#include "FileReader.h"
#include "./Penalty_Class/Penalty.h"
#include "./Penalty_Class/ThresholdPenalty.h"
#include "./Penalty_Class/DictionaryPenalty.h"
#include <queue>
#include <omp.h>
#include <algorithm>
#include <cmath>
#include <sys/stat.h>
#include <stdlib.h>
#include <random>
#include <chrono>

class ModScores
{
public:
	// Instantiate a new ModScore from file of folder path and dictionary of normalized K values (can optionally specify subset start and stop values)
	ModScores(std::string fof, const std::vector<unsigned int>& chromVect, const Penalty& penalty, const unsigned int maxLength, const unsigned int bufferStartBP, const unsigned int bufferEndBP);

	//= Calculate modularity scores and add them to queue
	bool getScores(std::vector<std::unordered_map<unsigned int, double>>& dictC);

protected:
	// Vector of variants in this chromosome in BP
	std::vector<unsigned int> m_chromVect;

	// Map where Key is BP of variant and Value is Index of variant
	std::unordered_map<unsigned int, unsigned int> m_BP2IndexMap;

	// Max length of an allowed LD block in BP
	unsigned int m_maxLength;

	// Penalty function
	const Penalty& m_penalty;

	// Start of the region on which mod scores are to be calculated (in BP)
	unsigned int m_bufferStartBP;

	// End of the region on which mod scores are to be calculated (in BP)
	unsigned int m_bufferEndBP;

	// Total number of variants
	unsigned int m_varCount;

	// Dictionary of S scores
	// Each vector position is the start point of an edge
	// The dictionary located at each vector position contains all endpoints as its keys
	// Values are the S value for the edge from start point (vector index) to endpoint (dictionary key)
	std::vector<std::unordered_map<unsigned int, double>> m_dictS;

	// Calculate dictionary of S scores with modularity penalty
	void calculateS(FileReader& fileStream, const Penalty& penalty);

};


#endif /* MODSCORES_H_ */
