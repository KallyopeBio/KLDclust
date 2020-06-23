/*
 * DictionaryK.h
 *
 *  Created on: Sep 13, 2018
 *      Author: Sarah Christensen
 * Description: Class takes as input a PLINK file reader and can calculate certain inputs to the modularity score from it
 */

#ifndef DICTIONARYK_H_
#define DICTIONARYK_H_

#include <unordered_map>
#include <cmath>
#include <iostream>
#include <omp.h>
#include "FileReader.h"

class DictionaryK
{

public:
	// Instantiate dictionary where keys are variants (in BP) and values are degree of variant node
	// Created from correlation file (loc_a, loc_b, score)
	DictionaryK(const std::string file);

	// Instantiate dictionary where keys are variants (in BP) and values are degree of variant node
	// Allow denominator of null hypothesis term to be scaled by factor c
	DictionaryK(const std::string file, double scaler);

	// Instantiate dictionary where keys are variants (in BP) and values are degree of variant node
	// Created from a vector of correlation files
	DictionaryK(const std::vector<std::string> inputFiles);

	// Instantiate (subset) of dictionary K from file already containing dictionary K values
	DictionaryK(const std::string dictK_file, const std::string m_file, const unsigned int start, const unsigned int stop);

	//= Get m
	double getM();

	//= Get K
	std::unordered_map<unsigned int, double> const & getK();

	// Print dictionary K to file
	void printK(const std::string file);

	// Print value m to file
	void printM(const std::string file);

protected:
	// Dictionary of K values (normalized in constructor by square root of 2m)
	std::unordered_map<unsigned int, double> m_dictK;

	// m is the total weight in the input graph (see modularity score equation)
	double m_M;

	// Get dictionary containing degree K_i for each variant i in the input data
	void calculateK(FileReader& fileStream, std::unordered_map<unsigned int, double>& dictK);

	// Calculate the total weight in the graph
	void calculateM();

	// Normalize K dictionary
	void normalize(double scaler);

};



#endif /* DICTIONARYK_H_ */
