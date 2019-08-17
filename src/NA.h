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


template <class T> 
class NA {
  public:
    static constexpr T value = std::is_floating_point<T>::value ? NAN : std::numeric_limits<T>::min();
};

// bool has no NA
template <> class NA<bool> {
public:
  static constexpr unsigned value = false;
};

template <typename T>
struct is_string {
  static const bool value = false;
};

template <class T, class Traits, class Alloc>
struct is_string<std::basic_string<T, Traits, Alloc>> {
  static const bool value = true;
};

/*
template <> class NAvalue<unsigned> {
public:
    static constexpr unsigned value = std::numeric_limits<unsigned>::max();
};
*/

template <typename T> 
bool is_NA(const T v) {
    if (std::is_floating_point<T>::value) {
        return std::isnan(v);
    } else {
        bool b = v == (NA<T>::value);
        return b;
  	}
}

template <> class NA<unsigned> {
public:
    static constexpr unsigned value = std::numeric_limits<unsigned>::max();
};


template <typename T>
void set_NA(std::vector<T> &v, double naflag) {
	if (!std::isnan(naflag)) {
		T flag = naflag;
		T navalue = NA<T>::value;
		std::replace(v.begin(), v.end(), flag, navalue);
	}
}

/*
class NA_long {
public:
    static constexpr long value = std::numeric_limits<long>::min();
};

class NA_unsigned {
public:
    static constexpr unsigned value = std::numeric_limits<unsigned>::max();
};

class NA_double {
public:
    static constexpr double value = NAN;
};

class NA_float {
public:
    static constexpr float value = NAN;
};


bool is_NAN(unsigned v) {
	return (v == NA_unsigned::value);
}

bool is_NAN(long v) {
    return (v == NA_long::value);
}

bool is_NAN(double v) {
    return std::isnan(v);
}

bool is_NAN(float v) {
    return std::isnan(v);
}

*/

