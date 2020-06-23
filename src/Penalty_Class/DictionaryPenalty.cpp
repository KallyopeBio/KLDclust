/*
 * DictionaryPenalty.cpp
 *
 *  Created on: Nov 29, 2018
 *      Author: sarah
 */

#include "DictionaryPenalty.h"

// Create constructor for this penalty
DictionaryPenalty::DictionaryPenalty(const std::unordered_map<unsigned int, double>& dictionary): m_dictionary(dictionary), m_coeff(1.0){


}

// Create constructor for this penalty
DictionaryPenalty::DictionaryPenalty(const std::unordered_map<unsigned int, double>& dictionary, const double coeff): m_dictionary(dictionary), m_coeff(coeff){

}

// Returns penalty for this derived penalty class
double DictionaryPenalty::getPenalty(unsigned int loc_a, unsigned int loc_b) const{
	return (m_coeff*m_dictionary.at(loc_a)*m_dictionary.at(loc_b));
}


