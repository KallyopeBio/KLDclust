/*
 * PwrChiSq.cpp
 *
 * An implementation of the function pwr.chisq.test in the R package pwr,
 * for the scenarios in which we wish to calculate w, the effect size, or
 * alpha, the significance level.
 *
 * Created on: Dec 17, 2018
 *     Author: Justin McManus, PhD
 */

#include "PwrChiSq.h"
#include <boost/math/distributions.hpp>
#include <boost/math/distributions/chi_squared.hpp>
#include <boost/math/distributions/non_central_chi_squared.hpp>
#include <boost/math/tools/toms748_solve.hpp>

// A functor to solve for the effect size, w
PwrChiSq::FuncW::FuncW(std::size_t n, double df, double alpha, double power):
m_n((double)n),
m_df(df),
m_alpha(alpha),
m_power(power)
{}

double PwrChiSq::FuncW::operator()(double w) 
{	
	boost::math::chi_squared_distribution<double> chisq_dist(m_df);
	double k= quantile(complement(chisq_dist, m_alpha));

	boost::math::non_central_chi_squared_distribution<double> chisq_nc_dist(m_df, m_n*w*w);
	return cdf(complement(chisq_nc_dist, k)) - m_power;
}

// A functor to solve for the significance threshold, alpha
PwrChiSq::FuncAlpha::FuncAlpha(std::size_t n, double w, double df, double power):
m_n((double)n),
m_w(w),
m_df(df),
m_power(power)
{}

double PwrChiSq::FuncAlpha::operator()(double alpha)
{
	boost::math::chi_squared_distribution<double> chisq_dist(m_df);
	double k= quantile(complement(chisq_dist, alpha));

	boost::math::non_central_chi_squared_distribution<double> chisq_nc_dist(m_df, m_n*m_w*m_w);
	return cdf(complement(chisq_nc_dist, k)) - m_power;
}

//= the minimum effect size, w, that can be detected at a given
//  power, for a specified sample size, N, and significance level,
//  alpha
double PwrChiSq::solve_w(std::size_t n, unsigned int df, double alpha, double power) 
{
	FuncW f(n, df, alpha, power);
	boost::math::tools::eps_tolerance<double> tol(64);
	boost::uintmax_t miter= 1000000;
	std::pair<double, double> bracket= toms748_solve(f, 1e-10, 1.0 - 1e-10, tol, miter);
	return (bracket.first + bracket.second)/2.0;
}

//= the significance level, alpha, at which the null hypothesis can be
//  rejected with the given power, for a specified sample size, N, and
//  effect size, w
double PwrChiSq::solve_alpha(std::size_t n, double w, unsigned int df, double power)
{
	FuncAlpha f(n, w, df, power);
	boost::math::tools::eps_tolerance<double> tol(64);
	boost::uintmax_t miter= 1000000;
	std::pair<double, double> bracket= toms748_solve(f, 1e-10, 1.0 - 1e-10, tol, miter);
	return (bracket.first + bracket.second)/2.0;
}

//= the critical value, r, such that Pr(X > r) = alpha, where X is
//  a Chi-squared random variable with the given degrees of freedom
double PwrChiSq::solve_cv(unsigned int df, double alpha)
{
	boost::math::chi_squared_distribution<double> chisq_dist((double)df);
	return quantile(complement(chisq_dist, alpha));
}

