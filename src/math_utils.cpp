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

#include <cmath>
#include <limits>
#include <vector>


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
  double d = pow(10, n);
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


bool is_equal(double a, double b, double error_factor=1.0) {
	return a==b || std::abs(a-b)<std::abs(std::min(a,b))*std::numeric_limits<double>::epsilon()*error_factor;
}

bool is_equal_range(double x, double y, double range, double tolerance) {
	return (fabs(x - y) / range) < tolerance ;
}

