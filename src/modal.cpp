// Copyright (c) 2018-2020  Robert J. Hijmans
//
// This file is part of the "spat" library.
//
// spat is free software: you can redistribute it and/or modify it
// under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 2 of the License, or
// (at your option) any later version.
//
// spat is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with spat. If not, see <http://www.gnu.org/licenses/>.

#include <vector>
#include <algorithm>
#include <cmath>
#include <random>

#include "math_utils.h"

/*
unsigned random_value(std::mt19937 &rng, std::uniform_int_distribution<std::mt19937::result_type>& dist) {
	return dist(rng);
}
std::mt19937 rng(seed);
unsigned rndmax = 10000
std::uniform_int_distribution<std::mt19937::result_type> dist(0, rndmax);
double x = random_value(rng, dist);
*/



double modal_value(std::vector<double> values, unsigned ties, bool narm, std::default_random_engine rgen, std::uniform_real_distribution<double> dist) {
	
	if (narm) {
		na_omit(values);
	}
	size_t n = values.size();
	if (n == 0) return (NAN);
	if (n == 1) return (values[0]);		
    std::vector<unsigned> counts(n, 0);

	if (ties < 3) {
		std::sort(values.begin(), values.end());
	}

	
    for (size_t i=0; i<n; ++i) {
        counts[i] = 0;
        size_t j = 0;
        while ((j < i) && (values[i] != values[j])) {
            ++j;
        }
        ++(counts[j]);
    }
	
    size_t maxCount = 0;
	// first (lowest due to sorting)
	if (ties == 0) {
		for (size_t i = 1; i < n; ++i) {
			if (counts[i] > counts[maxCount]) {
				maxCount = i;
			}
		}
	// last	
	} else if (ties == 1) {
		for (size_t i = 1; i < n; ++i) {
			if (counts[i] >= counts[maxCount]) {
				maxCount = i;
			}
		}

	// dont care (first, but not sorted)
	} else if (ties == 2) {
		for (size_t i = 1; i < n; ++i) {
			if (counts[i] > counts[maxCount]) {
				maxCount = i;
			}
		}

	// random
	} else if (ties == 3) {
		size_t tieCount = 1;
		for (size_t i = 1; i < n; ++i) {
			if (counts[i] > counts[maxCount]) {
				maxCount = i;
				tieCount = 1;
			} else if (counts[i] == counts[maxCount]) {
				tieCount++;
				double rand = dist(rgen);
				if (rand < (1 / tieCount)) {
					maxCount = i;
				}			
			}
		}		
	} else {
		size_t tieCount = 1;
		for (size_t i = 1; i < n; ++i) {
			if (counts[i] > counts[maxCount]) {
				maxCount = i;
				tieCount = 1;
			} else if (counts[i] == counts[maxCount]) {
				tieCount++;
			}
		}
		if (tieCount > 1 ) {
			return(NAN);
		}
	}
	
    return values[maxCount];
}

