/*
 * ThresholdPenalty.cpp
 *
 *  Created on: Nov 29, 2018
 *      Author: sarah
 */

#include "ThresholdPenalty.h"

// Create constructor for this penalty
ThresholdPenalty::ThresholdPenalty(const double threshold){
	// Set member variable
	m_threshold=threshold;
}

// Returns penalty for this derived penalty class
double ThresholdPenalty::getPenalty() const {
	return m_threshold;
}

// Returns penalty for this derived penalty class
double ThresholdPenalty::getPenalty(unsigned int loc_a, unsigned int loc_b) const {
	return m_threshold;
}

// Calculate background penalty
long double ThresholdPenalty::background(std::vector<std::string> inputFiles, unsigned int distance){

	// Initialize values that will be reduced across threads
	unsigned int edgesCount(0);
	long double sum(0);
	long double squareSum(0);
	long double n(0);

	#pragma omp parallel for reduction(+:edgesCount) reduction(+:sum) reduction(+:squareSum) reduction(+:n)
	for (unsigned int i=0; i<inputFiles.size(); i++){

		// Get summary statistics on long range correlations
		unsigned int loc_a(0);
		unsigned int loc_b(0);
		double r_score(0);

		FileReader fileStream(inputFiles[i]);
		while(fileStream.getNextRead(loc_a, loc_b, r_score)){
			// Count all edges
			edgesCount += 1;

			// Check if distance is longer that allowed LD block in BP
			if ((loc_b-loc_a)> distance){
				n += 1;
				sum += r_score;
				squareSum += (r_score*r_score);
			}
		}
	}

	if (n==0){
		throw std::runtime_error("No correlation observations were greater than the specified distance so the threshold cannot be calculated.");
	}

	// Calculate standard deviation
	long double variance = (squareSum*n)-(sum*sum);
	variance /= (n*n);
	long double std = std::sqrt(variance);
	long double mean = sum/n;


	// Report values
	std::cout << "\t Number of input correlations: " << edgesCount << std::endl;
	std::cout << "\t Number of long-range correlations: " << n << std::endl;
	std::cout << "\t Mean: " << mean << std::endl;
	std::cout << "\t Standard deviation: " << std << std::endl;

	return mean+(2*std);
}

// Calculate Fisher penalty
long double ThresholdPenalty::fisher(unsigned int obsCount){

	// Compute r threshold at the 5% level of significance which corresponds 1.96 critical value
	long double thresholdR = 1.96/(std::sqrt(obsCount-3.0));

	// The Fisher transformation is just the inverse hyperbolic function
	// We take the inverse of F to solve for r, so we get back tanh
	thresholdR = std::tanh(thresholdR);

	// We convert this R threshold to R2
	long double thresholdRR = thresholdR*thresholdR;

	return thresholdRR;
}

// Calculate Fisher penalty
long double ThresholdPenalty::chiSquare(unsigned int obsCount, double alpha){

	// rr threshold at the 5% level of significance with DF (2-1)(2-1) = 1 corresponds to 3.841 chi squared critical value
	// Note the relationship between chi square and the phi correlation X^2/N = r^2
	long double thresholdRR = PwrChiSq::solve_cv(1, alpha)/obsCount;

	return thresholdRR;
}

// Calculate Chi Square Power penalty
long double ThresholdPenalty::chiSquarePower(unsigned int obsCount, double effectSize, double power){
	double alpha = PwrChiSq::solve_alpha(obsCount, effectSize, 1, power);
	long double thresholdRR = PwrChiSq::solve_cv(1, alpha)/obsCount;

	std::cout << "The Chi Square Power alpha value is: " << alpha << std::endl;

	return thresholdRR;
}
