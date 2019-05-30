// Copyright (c) 2018-2019  Robert J. Hijmans
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

double modal_value(std::vector<double> values, int ties) {
	int n = values.size();
    std::vector<unsigned> counts(n, 0);

	if (ties < 3) {
		std::sort(values.begin(), values.end());
	}
	
    for (int i = 0; i < n; ++i) {
        counts[i] = 0;
        int j = 0;
        while ((j < i) && (values[i] != values[j])) {
            ++j;
        }
        ++(counts[j]);
    }
	
    int maxCount = 0;
	// first (lowest due to sorting)
	if (ties == 0) {
		for (int i = 1; i < n; ++i) {
			if (counts[i] > counts[maxCount]) {
				maxCount = i;
			}
		}
	// last	
	} else if (ties == 1) {
		for (int i = 1; i < n; ++i) {
			if (counts[i] >= counts[maxCount]) {
				maxCount = i;
			}
		}

	// dont care (first, but not sorted)
	} else if (ties == 2) {
		for (int i = 1; i < n; ++i) {
			if (counts[i] > counts[maxCount]) {
				maxCount = i;
			}
		}

	// random
	/* } 
		else if (ties == 3) {
		int tieCount = 1;
		for (int i = 1; i < n; ++i) {
			if (counts[i] > counts[maxCount]) {
				maxCount = i;
				tieCount = 1;
			} else if (counts[i] == counts[maxCount]) {
				tieCount++;
				if (R::runif(0,1) < (1.0 / tieCount)) {
					maxCount = i;
				}			
			}
		}
	*/	
   // NA		
	} else {
		int tieCount = 1;
		for (int i = 1; i < n; ++i) {
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

