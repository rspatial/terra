#include <vector>
#include <algorithm>
#include <cmath>
using namespace std;

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

