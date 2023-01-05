// Copyright (c) 2018-2023  Robert J. Hijmans
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

#include <cmath>
#include <limits>
#include <vector>
#include <algorithm>
#include <random>

/*

std::vector<double> mean2d(const std::vector<std::vector<double>> &x) {
	size_t n = x[0].size();
	size_t nn = x.size();
	std::vector<double> out(n, NAN);
	size_t d;
	double v;
	for (size_t i=0; i<n; i++) {
		v = 0;
		d = 0;
		for (size_t j=0; j<nn; j++) {
			if (!std::isnan(x[i][j])) {
				v += x[i][j];
				d++;
			}
		}
		if (d > 0) {
			out[i] = v / d;
		}
	}
	return out;
}

*/


void na_omit(std::vector<double> &x) {
	x.erase(std::remove_if(std::begin(x), std::end(x),
        [](const double& value) { return std::isnan(value); }),
        std::end(x));
}


void vector_minmax(std::vector<double> v, double &min, int &imin, double &max, int &imax) {
    std::vector<double>::size_type p=0;
    imax = -1; imin=-1;
    min = std::numeric_limits<double>::max();
    max = std::numeric_limits<double>::lowest();
    for (auto &val : v) {
		if (!std::isnan(val)) {
			if (val > max) {
				imax = p;
				max = val;
			}
			if (val < min) {
				imin = p;
				min = val;
			}
		}
        p++;
    }
	if (imax == -1) {
		max = NAN;
		min = NAN;
	}
}


double roundn(double x, int n){
	double d = pow(10.0, n);
	return std::round(x * d) / d;
}

double signif(double x, unsigned n) {
	double b = x;
	unsigned i;
	for (i = 0; b >= 1; ++i) {
		b = b / 10;
	}
	int d = n-i;
	return roundn(x, d);
}

bool is_equal(double a, double b, double tolerance=10.0) {

	double tol = std::max(tolerance, std::abs(std::min(a,b))) * std::numeric_limits<double>::epsilon();
	return ((a==b) || (std::abs(a-b) < tol) );
}

bool about_equal(double a, double b, double tolerance) {
	return ((a==b) || (std::abs(a-b) < tolerance));
}

bool is_equal_relative(double a, double b, double tolerance) {
	tolerance = std::max(fabs(a), fabs(b)) * tolerance;
    return about_equal(a, b, tolerance);
}

bool is_equal_range(double x, double y, double range, double tolerance) {
	return (fabs(x - y) / range) < tolerance ;
}



double median(const std::vector<double>& v) {
	size_t n = v.size();
	std::vector<double> vv;
	vv.reserve(n);
	for (size_t i=0; i<n; i++) {
        if (!std::isnan(v[i])) {
            vv.push_back(v[i]);
        }
	}
	n = vv.size();
	if (n == 0) {
		return(NAN);
	}
	size_t n2 = n / 2;
	std::nth_element(vv.begin(), vv.begin()+n2, vv.end());
	double med = vv[n2];
	return med;
}



std::vector<double> movingMedian(const std::vector<double> &x, size_t n) {
	std::vector<double> out(x.size());
	std::vector<double> d(n, NAN);
	size_t half = (n/2);
	size_t half1 = half+1;
	// fill left side
	for (size_t i=0; i<half; i++) {
		for (size_t j=0; j< (half1+i); j++) {
			d[j] = x[j];
		}
		out[i] = median(d);
	}
	// middle
	size_t maxn = out.size() - half;
	std::vector<double> v;
	for (size_t i=half; i<maxn; i++) {
		v = std::vector<double>(x.begin()+i-half, x.begin()+i+half1);
 		out[i] = median(v);
	}
	// right side
	int j=0;
	for (size_t i=maxn; i<out.size(); i++) {
		v[j++] = NAN;
		out[i] = median(v);
	}
	return(out);
}



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

