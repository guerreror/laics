/*
 * ran_mk
 *  InvCoal
 *  Version:  labp_v17
 *
 * Routines for generating various kinds of random numbers.
 *
 * Started by Shane Pope & Mark Kirkpatrick, August 2008
 *
 * Modified by Rafael Guerrero 2009-2016
 * Included exponentially distributed numbers and standard-library random distributions.
 *
 */

#ifndef RANDOM_
#define RANDOM_


#include <random>

extern std::mt19937_64 gen;	// Declared externally in main.cpp so that it can be initialized there and called in random.cpp

inline std::uniform_int_distribution<>::result_type randint(int lower_bound, int upper_bound)
{
//
// randint(a, b) returns a uniformly distributed random integer on [a, b]
//
    return std::uniform_int_distribution<>(lower_bound, upper_bound)(gen);
}



inline std::uniform_real_distribution<>::result_type randreal(double lower_bound = 0, double upper_bound = 1)
{
//
// Returns a uniformly distributed random deviate.  
//	randreal() is uniform on [0, 1]
//	randreal(a, b) is uniform on [a, b].

    return std::uniform_real_distribution<>(lower_bound, upper_bound)(gen);
}



inline std::exponential_distribution<>::result_type randexp(double lambda)
{
//
// Returns a random deviate from an exponential distribution with parameter lambda
//
    return std::exponential_distribution<>(lambda)(gen);
}

inline int randpois( double n ) {
    using dist_t = std::poisson_distribution<int>;
    static dist_t d{n};
    return d(gen);
}

inline std::binomial_distribution<>::result_type randbinom(unsigned int trials, double p)
{
    //
    // Returns a random deviate from an poisson distribution with parameter lambda
    //
    return std::binomial_distribution<>(trials, p)(gen);
}

#endif //RANDOM_
