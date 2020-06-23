/*
 * penalty.h
 *
 *  Created on: Nov 29, 2018
 *      Author: sarah
 */

#ifndef PENALTY_CLASS_PENALTY_H_
#define PENALTY_CLASS_PENALTY_H_

class Penalty
{
public:
	// Return penalty
	virtual double getPenalty(unsigned int loc_a, unsigned int loc_b) const = 0;

	// Create virtual destructor
	virtual ~Penalty(){};

};



#endif /* PENALTY_CLASS_PENALTY_H_ */
