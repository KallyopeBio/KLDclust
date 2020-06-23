/*
 * DictionaryPenalty.h
 *
 *  Created on: Nov 29, 2018
 *      Author: sarah
 */

#ifndef PENALTY_CLASS_DICTIONARYPENALTY_H_
#define PENALTY_CLASS_DICTIONARYPENALTY_H_

#include <unordered_map>
#include "../Penalty_Class/Penalty.h"

class DictionaryPenalty : public Penalty
{
public:
	// Create constructor
	DictionaryPenalty(const std::unordered_map<unsigned int, double>& dictionary);

	// Create constructor
	DictionaryPenalty(const std::unordered_map<unsigned int, double>& dictionary, const double coeff);

	// Return penalty
	double getPenalty(unsigned int loc_a, unsigned int loc_b) const;

private:
	const std::unordered_map<unsigned int, double>& m_dictionary;

	const double m_coeff;

};



#endif /* PENALTY_CLASS_DICTIONARYPENALTY_H_ */
