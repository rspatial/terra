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

#include "spatRaster.h"
#include <limits>
#include <set>


std::map<double, unsigned> table(std::vector<double> &v) {
	std::map<double, unsigned> count;
	for_each( v.begin(), v.end(), [&count]( double val ){ 
			if(!std::isnan(val)) count[val]++; 
		} 
	);	
	return count;
}


std::map<double, unsigned> ctable(std::map<double, unsigned> &x, std::map<double, unsigned> &y) {
	for(auto p : y) {
		x[p.first] += p.second;
	}	
	return(x);
}


std::vector<double> vtable(std::map<double, unsigned> &x) {
	std::vector<std::vector<double>> out(2);
	for( auto p : x ) {
		out[0].push_back(p.first);
		out[1].push_back(p.second);
	}  
	out[0].insert(out[0].end(), out[1].begin(), out[1].end());
	return out[0];
}



std::vector<std::vector<double>> SpatRaster::freq(bool bylayer) {
	std::vector<std::vector<double>> out;
	if (!hasValues()) return out;
	BlockSize bs = getBlockSize(4);
	unsigned nc = ncol();
	unsigned nl = nlyr();
	readStart();
	if (bylayer) {
		out.resize(nl);
		std::vector<std::map<double, unsigned>> tabs(nl);
		for (size_t i = 0; i < bs.n; i++) {
			unsigned n = bs.nrows[i] * nc;
			std::vector<double> v = readValues(bs.row[i], bs.nrows[i], 0, nc);
			for (size_t lyr=0; lyr<nl; lyr++) {
				unsigned off = lyr*n;
				std::vector<double> vv(v.begin()+off, v.begin()+off+n);
				std::map<double, unsigned> tab = table(vv);
				tabs[lyr] = ctable(tabs[lyr], tab);
			}
		}
		for (size_t lyr=0; lyr<nl; lyr++) {
			out[lyr] = vtable(tabs[lyr]);
		}
	} else {
		out.resize(1);
		std::map<double, unsigned> tabs;
		for (size_t i = 0; i < bs.n; i++) {
			std::vector<double> v = readValues(bs.row[i], bs.nrows[i], 0, nc);
			std::map<double, unsigned> tab = table(v);
			tabs = ctable(tabs, tab);
		}
		out[0] = vtable(tabs);
	}
	readStop();
	return(out);
}

