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

#ifndef SPATVARIANT_GUARD
#define SPATVARIANT_GUARD

#include <vector>
#include <string>
#include "date.h"

#include <vecmath>

class SpatFactor {
public:
	std::vector<unsigned> v;
	std::vector<unsigned> levels;
	std::vector<std::string> labels;
  
	void compute_levels() {
		levels = vunique(v);
		labels = vtostring(levels);
	}

	void set_values(std::vector<unsigned> _values) {
		v = _values;
		levels.resize(0);
		labels.resize(0);		
	}

	void set_levels(std::vector<unsigned> _levels, std::vector<std::string> _labels) {
		levels = _levels;
		labels = _labels;
		labels.resize(levels.size());
	}
	
	void set_labels(std::vector<std::string> _labels) {
		labels = _labels;
		labels.resize(levels.size());
	}
}


class SpatVariant {
public:
  unsigned type=0;
  std::vector<double> d;
  std::vector<int> i;
  std::vector<std::string> s;
  std::vector<unsigned> u;
  std::vector<bool> b;
  std::vector<SpatDate> date;
  std::vector<SpatFactor> fact;
  
  const unsigned ntypes = 7;
  
  SpatVariant();
  SpatVariant(const unsigned _type) {
    set_type(_type);
  }  
  SpatVariant(unsigned _type, const size_t _size) {
    set_type(_type);
 	resize(_size);
  };

  void set_type(const unsigned _type){
    if (_type < ntypes) type = _type;
  }
  
  size_t size() {
    if (type == 0) return d.size();
    if (type == 1) return i.size();
    if (type == 2) return s.size();
    if (type == 3) return u.size();
    if (type == 4) return b.size();
    if (type == 5) return date.size();
    if (type == 6) return fact.size();
	return 0;
  }

  void resize(size_t _size) {
    if (type == 0) d.resize(_size);
    else if (type == 1) i.resize(_size);
    else if (type == 2) s.resize(_size);
    else if (type == 3) u.resize(_size);
    else if (type == 4) b.resize(_size);
    else if (type == 5) date.resize(_size);
    if (type == 6) return fact.resize(_size);
    else vdate.resize(_size);
  }
  
};


