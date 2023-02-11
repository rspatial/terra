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

#ifndef SPATFACTOR_GUARD
#define SPATFACTOR_GUARD


class SpatFactor {
public:
	virtual ~SpatFactor(){}
	SpatFactor(){} ;
	SpatFactor(size_t _size, unsigned _value) {
		v.resize(_size, _value);
	};

	SpatFactor(std::vector<unsigned> _values, std::vector<std::string> _labels);
	SpatFactor(std::vector<unsigned> _values);
	SpatFactor(std::vector<std::string> _values);

	
	std::vector<unsigned> v;
	//std::vector<unsigned> levels;
	std::vector<std::string> labels;
  
	size_t size() { return v.size(); }
	bool empty() { return v.empty(); }
	
	//void compute_levels();
	void push_back(unsigned x) { v.push_back(x); }

	
	bool set_labels(std::vector<std::string> _labels);
	
	void reserve(size_t n) { v.reserve(n); }
	void resize(size_t n) { v.resize(n); }
	void resize(size_t n, unsigned x) {v.resize(n, x);}	
	
//	template <typename T>
//	  SpatFactor(std::vector<T> _v) {
//	   set_values(_v);
//	}
	
	SpatFactor subset(std::vector<unsigned> i);
	std::string getLabel(size_t i); 
	std::vector<std::string> getLabels();
	
};

#endif

