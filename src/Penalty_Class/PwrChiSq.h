/*
 * PwrChiSq.h
 *
 * An implementation of the function pwr.chisq.test in the R package pwr,
 * for the scenarios in which we wish to calculate w, the effect size, or
 * alpha, the significance level.
 *
 * Created on: Dec 17, 2018
 *     Author: Justin McManus, PhD
 */

#ifndef PWR_CHI_SQ_H_
#define PWR_CHI_SQ_H_

#include <iostream>

class PwrChiSq
{
public:
	//= the minimum effect size, w, that can be detected at a given
	//  power, for a specified sample size, N, and significance level,
	//  alpha
	static double solve_w(std::size_t n, unsigned int df, double alpha, double power);

	//= the significance level, alpha, at which the null hypothesis can be
	//  rejected with the given power, for a specified sample size, N, and
	//  effect size, w
	static double solve_alpha(std::size_t n, double w, unsigned int df, double power);

	//= the critical value, r, such that Pr(X > r) = alpha, where X is
	//  a Chi-squared random variable with the given degrees of freedom
	static double solve_cv(unsigned int df, double alpha);

private:
	// A functor to solve for w
	class FuncW {
	public:
		FuncW(std::size_t n, double df, double alpha, double power);
		double operator()(double w);
	private:
		double m_n;
		double m_df;
		double m_alpha;
		double m_power;
	};

	// A functor to solve for alpha
	class FuncAlpha {
	public:
		FuncAlpha(std::size_t n, double w, double df, double power);
		double operator()(double alpha);
	private:
		double m_n;
		double m_w;
		double m_df;
		double m_power;
	};
};

#endif
