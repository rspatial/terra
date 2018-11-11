// Copyright (c) 2018  Robert J. Hijmans
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

bool is_equal(double a, double b, double error_factor=1.0);
bool is_equal_range(double x, double y, double range, double tolerance);
void vector_minmax(std::vector<double> v, double &min, int &imin, double &max, int &imax);
double roundn(double x, int n);

template <typename Iterator>
void minmax(Iterator start, Iterator end, double &min, double &max) {
    min = std::numeric_limits<double>::max();
    max = std::numeric_limits<double>::lowest();
    bool none = true;
	for (Iterator v = start; v !=end; ++v) {
		if (!std::isnan(*v)) {
			if (*v > max) {
				max = *v;
                none = false;
			}
			if (*v < min) {
				min = *v;
			}
		}
    }
    if (none) {
        min = NAN;
        max = NAN;
    }
}

