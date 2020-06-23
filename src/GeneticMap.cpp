/*
 * GeneticMap.cpp
 *
 *  Created on: Nov 13, 2018
 *      Author: sarah
 */

#include "GeneticMap.h"

// Instantiate a new file from file path
GeneticMap::GeneticMap(const std::string file)
{
	std::cout << "Reading in genetic map" << std::endl;
	std::string SNP;
	unsigned int loc;
	double val;
	std::ifstream in(file);
	while(true){
		in >> SNP;
		in >> loc;
		in >> val;
		if( in.eof() ) break;
		// Add conversion to dictionary
		m_locMap[loc]=val;
	}
	in.close();
}

// Given a variant in BP, return genetic map value
double GeneticMap::getValue(unsigned int loc){return m_locMap.at(loc);}
