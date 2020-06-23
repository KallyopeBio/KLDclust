/*
 * VCF.h
 *
 *  Created on: Nov 6, 2018
 *      Author: sarah
 */

#ifndef VCF_H_
#define VCF_H_

#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <tuple>
#include <vector>
#include <cmath>
#include <algorithm>
#include <unordered_map>
#include <unordered_set>
#include <omp.h>
#include <limits>
#include "Chromosome.h"
#include "GeneticMap.h"

class VCF
{
public:
	// Instantiate a new VCF file instance from file path
	VCF(const std::string file, const unsigned int size, const unsigned int maxLength, const bool phased, const std::string indFile);

	// Calculate correlations between input variants above some threshold level
	void calculateShrinkageCorrelation(std::string outputDir, GeneticMap gmap, unsigned int NE, double minCorr);

	// Calculate correlations between input variants above some threshold level
	void calculateMAFCorrelation(std::string outputDir, double minCorr, double maf);

	// Get critical value for R2 at 5% level
	unsigned int getObsCount(bool phased);

protected:

	// Max window length in BP
	unsigned int m_maxLength;

	// Set containing the column numbers for individuals of interest
	std::unordered_set<int> m_individualsSubset;

	// Vector containing the list of variants (note there may be duplicates if multiple SNPs at same site)
	std::vector<unsigned int> m_variants;

	// Set containing the variants that have been eliminated for one of the following reasons
	// All missing for individuals of interest
	// No variation
	// Multiple SNPs at same location
	std::unordered_set<int> m_variantsRemoved;

	// A dictionary where the variant location (in BP) are keys and the values are the number of SNPs that have that location
	std::unordered_map<unsigned int, unsigned int> m_variantCounts;

	// Vector where each position corresponds to a variant
	// Each position contains a tuple containing various counts (# score 0, # score 1, # score 2, # missing, # non-missing)
	std::vector<std::tuple<unsigned int, unsigned int, unsigned int, unsigned int, unsigned int>> m_counts;

	// Vector where each position corresponds to a variant
	// Each position contains the average score for this variant
	std::vector<long double> m_avgScores;

	// Vector where each position corresponds to a variant
	// Each position contains a dictionary of individuals (key based on column number) who have non-zero scores for that variant (value)
	// A score of -1 indicates that the individual was missing data for that variant
	std::vector<std::unordered_map<int, int>> m_indScores;

	// Add variants to set m_variantsRemoved that should be removed from analysis
	void fillVariantsRemoved();

	// Parse individuals included in sample
	void parseIndividuals(const std::string indFile, const std::string line);

	// Parse chromosome location in data line
	void parseChromLoc(const unsigned int loc, const unsigned int index);

	// Parse input file as phased data
	void parseLine(const std::string line, const unsigned int index, bool phased);

	// Parse input file as phased data
	void parsePhasedData(const std::string line, const unsigned int index, const unsigned int column);

	// Parse input file as unphased data
	void parseUnphasedData(const std::string line, const unsigned int index, const unsigned int column);

	// Calculate theta parameter based on Li and Stephens (2003); EQ 2.8 in Wen and Stephens 2010
	double calculateTheta();

	// Calculate empirical covariance
	std::tuple<unsigned int, unsigned int> calculateInteraction(unsigned int i, unsigned int j);

};




#endif /* VCF_H_ */
