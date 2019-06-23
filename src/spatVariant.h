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

class SpatVariant {
public:
  unsigned type=0;
  std::vector<double> vd;
  std::vector<int> vi;
  std::vector<std::string> vs;
  std::vector<unsigned> vu;
  std::vector<bool> vb;
  std::vector<SpatDate> vdate;
  
  const unsigned ntypes = 5;
  
  SpatVariant();
  SpatVariant(const unsigned _type) {
    set_type(_type);
  }  
  SpatVariant(unsigned _type, const unsigned _size) {
    set_type(_type);
 	resize(_size);
  };

  void set_type(const unsigned _type){
    if (_type < ntypes) type = _type;
  }
  
  unsigned size() {
    if (type == 0) return vd.size();
    if (type == 1) return vi.size();
    if (type == 2) return vs.size();
    if (type == 3) return vu.size();
    if (type == 4) return vb.size();
    return vdate.size();
  }

  void resize(unsigned _size) {
    if (type == 0) vd.resize(_size);
    else if (type == 1) vi.resize(_size);
    else if (type == 2) vs.resize(_size);
    else if (type == 3) vu.resize(_size);
    else if (type == 4) vb.resize(_size);
    else vdate.resize(_size);
  }
  
};
