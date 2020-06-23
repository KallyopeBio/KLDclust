/*
 * GeneticMap.h
 *
 *  Created on: Nov 13, 2018
 *      Author: sarah
 */

#ifndef GENETICMAP_H_
#define GENETICMAP_H_

#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <tuple>
#include <vector>
#include <cmath>
#include <algorithm>
#include <unordered_map>

class GeneticMap
{
public:
	// Instantiate a new VCF file instance from file path
	GeneticMap(const std::string file);

	// Given a variant in BP, return genetic map value
	double getValue(unsigned int loc);

protected:
std::unordered_map<unsigned int, double> m_locMap;

};



#endif /* GENETICMAP_H_ */
