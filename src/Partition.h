/*
 * Partition.h
 *
 *  Created on: Sep 18, 2018
 *      Author: Sarah Christensen
 *      Description: Calculates the optimal modularity score corresponding to some partition. Valid partitions must correspond to cuts along the chromosome.
 */

#ifndef PARTITION_H_
#define PARTITION_H_

#include "ModScores.h"
#include "Chromosome.h"
#include <omp.h>
#include <algorithm>

class Partition
{
public:
	// Instantiate new partition from chromosome object and community modularity score dictionary
	Partition(Chromosome& chrom, const unsigned int maxLength, const unsigned int startBP, const unsigned int endBP, std::vector<std::unordered_map<unsigned int, double> >& dictC, std::string outputFolder);

	//= Return optimal modularity score
	double getOpt();

	//= Return partition that achieves optimal modularity score
	bool getOptPartition(std::vector<std::pair<unsigned int, unsigned int> >& optPartition);

	// Print optimal scores to file
	void printScores(const std::string file);

	// Print optimal breaks to file
	void printBreaks(const std::string file);

protected:
	// Chromosome object
	Chromosome m_chrom;

	// Max length of allowed LD block
	unsigned int m_maxLength;

	// Vector containing variants in sorted order
	std::vector<unsigned int> m_chromVect;

	// For each variant location, the optimal modularity score ending at that point
	std::unordered_map<unsigned int, double> m_optScores;

	// For each variant location, a pointer back to the start position for the break point that achieved the optimal score
	std::unordered_map<unsigned int, unsigned int> m_optBreaks;

	// This is an unfortunate variable to help with reindexing
	// The cDict vector/dictionary is only over a chunk of the chromosome while the other dictionaries are over the full chromosome
	int m_startCIndex;

	// Calculates optimal partition for a given index
	bool calculatePartition(unsigned int initIndex, unsigned int initLoc, unsigned int leftBound, std::vector<std::unordered_map<unsigned int, double> >& dictC);

};



#endif /* PARTITION_H_ */
