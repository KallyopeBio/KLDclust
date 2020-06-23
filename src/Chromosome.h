/*
 * Chromosome.h
 *
 *  Created on: Sep 14, 2018
 *      Author: Sarah Christensen
 */

#ifndef CHROMOSOME_H_
#define CHROMOSOME_H_

#include <math.h>
#include <unordered_map>
#include <set>
#include <limits>
#include <omp.h>
#include <bits/stdc++.h>
#include <iostream>
#include <sys/stat.h>
#include <sys/types.h>
#include "FileReader.h"


class Chromosome
{
public:
	// Instantiate a new chromosome instance from correlation (3 column) file path
	static Chromosome correlationFile(std::string filepath);

	// Instantiate a new chromosome instance from vector of file paths
	static Chromosome correlationFiles(const std::vector<std::string> chromVect);

	// Instantiate a new chromosome instance from a file just containing a list of variants
	static Chromosome parsedFile(std::string filepath);

	// Instantiate a new chromosome instance from vector
	Chromosome(const std::vector<unsigned int>& chromVect);

	// Instantiate a new chromosome instance from double dictionary
	Chromosome(const std::unordered_map<unsigned int, double>& dictK);

	// Instantiate a new chromosome instance from unsigned int dictionary
	Chromosome(const std::unordered_map<unsigned int, unsigned int>& dict);

	//= Get chromosome variants returned as a sorted vector
	std::vector<unsigned int> getVector();

	// Fill dictionary with key (BP location) and value (index location)
	void fillBP2IndexMap(std::unordered_map<unsigned int, unsigned int> & BP2IndexMap);

	// Subset this chromosome to specified range
	void subsetChromosome(const unsigned int subsetStart, const unsigned int subsetEnd);

	//= Get beginning of the chromosome
	unsigned int getStart();

	//= Get end of the chromosome
	unsigned int getEnd();

	//= Get length of the chromosome in BP
	unsigned int getLengthBP();

	//= Get length of the chromosome in number of variants
	unsigned int getLengthCount();

	//= Get farthest variant within x BP away to the right from variant at some index location
	unsigned int getRightNeighborIndex(const unsigned int locBP, const unsigned int distanceBP);

	//= Get farthest variant within x BP away to the left from variant at some index location
	unsigned int getLeftNeighborIndex(const unsigned int locBP, const unsigned int distanceBP);

	//= Get farthest variant within x BP away to the right from variant at some BP location
	// Returns BP of this variant
	unsigned int getRightNeighborBP(const unsigned int locBP, const unsigned int distanceBP);

	//= Get farthest variant within x BP away to the left from variant at some BP location
	// Returns BP of this variant
	unsigned int getLeftNeighborBP(const unsigned int locBP, const unsigned int distanceBP);

	//= Get a vector of tuples containing intervals with (approximately) at most maxVars variants and will be correct for LD blocks of at most maxLength
	std::vector<std::tuple<unsigned int, unsigned int, unsigned int>> getIntervals(const unsigned int maxLength, const unsigned int maxVars);

	// Subset correlation file based on chromosome intervals
	static void subsetCorrelationFile(const std::string file, std::vector<std::tuple<unsigned int, unsigned int, unsigned int>> intervals, const std::string outputFolder, const unsigned int numOutFiles);

	// Print chromosome vector to file
	void printChrom(const std::string file);

protected:
	// Start of chromosome in BP
	unsigned int m_start;

	// End of chromosome in BP
	unsigned int m_end;

	// Vector of varients in chromosome in BP in sorted order
	std::vector<unsigned int> m_variants;

	// Constructor for file paths
	Chromosome(const std::string file, const bool calculate);

};




#endif /* CHROMOSOME_H_ */
