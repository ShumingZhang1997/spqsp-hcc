#include "RNG.h"
#include <iostream>
#include <algorithm>
#include <numeric>

RNG::RNG()
	:_rng()
	, _unif_01_gen(_rng, UNI_01_DIST())
	, _norm_std_gen(_rng, BOOST_NORMAL_DIST(0,1))
	, _exponential_1_gen(_rng, BOOST_EXP_DIST(1.0))
{
}


RNG::~RNG()
{
}
//
void RNG::seed(unsigned int s){
	_rng.seed(s);
}

int RNG::getInt(int n) //Return a random number from 0 to n - 1
{
	return _unif_01_gen() * n;
}

double RNG::get_unif_01(){
	return _unif_01_gen();
}

double RNG::get_norm_std(){
	return _norm_std_gen();
}

double RNG::get_normal(double mean, double sd){
	double x = get_norm_std();
	return x*sd + mean;
}

double RNG::get_exponential(double mean){
	return _exponential_1_gen() * mean;
}

/*! \param[in]: cdf: vector of double, size n
	cumulative probability of n events. sorted low to high
*/
int RNG::sample_cdf(const std::vector<double>& cdf){
	double p = get_unif_01();
	return std::lower_bound(cdf.begin(), cdf.end(), p) - cdf.begin();
}

std::vector<double> RNG::pdf_to_cdf(const std::vector<double>& pdf) {
	double sum = std::accumulate(pdf.begin(), pdf.end(), 0.0);
	std::vector<double> result;
	result.reserve(pdf.size());
	double partial_sum = 0;
	for (double const& val : pdf) {
		partial_sum += val;
		result.push_back(partial_sum / sum);
	}
	return result;
}
