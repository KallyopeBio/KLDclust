/*
 * ThresholdPenalty.h
 *
 *  Created on: Nov 29, 2018
 *      Author: sarah
 */

#ifndef PENALTY_CLASS_THRESHOLDPENALTY_H_
#define PENALTY_CLASS_THRESHOLDPENALTY_H_

#include <vector>
#include <string>
#include <math.h>
#include <omp.h>
#include "../FileReader.h"
#include "../Penalty_Class/PwrChiSq.h"
#include "../Penalty_Class/Penalty.h"

class ThresholdPenalty : public Penalty
{
public:
	// Create constructor
	ThresholdPenalty(const double threshold);

	// Return penalty (locations in this case are ignored since it is constant)
	double getPenalty(unsigned int loc_a, unsigned int loc_b) const;

	// Return penalty
	double getPenalty() const;

	// Calculate background penalty
	static long double background(std::vector<std::string> inputFiles, unsigned int distance);

	// Calculate Fisher penalty
	static long double fisher(unsigned int obsCount);

	// Calculate Chi Square P value penalty
	static long double chiSquare(unsigned int obsCount, double alpha);

	// Calculate Chi Square Power penalty
	static long double chiSquarePower(unsigned int obsCount, double effectSize, double power);

private:
	double m_threshold;

};



#endif /* PENALTY_CLASS_THRESHOLDPENALTY_H_ */
