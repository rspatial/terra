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
#include "vecmath.h"



void do_stats(std::vector<double> &v, std::string fun, bool narm, double &stat, unsigned &n) {
	double s;
	if (fun == "sum") {
		s = vsum(v, narm);
		stat = stat + s;
	} else if (fun == "mean") {
		stat = vsum(v, narm);
		std::vector<bool> nna = visnotna(v);
		for (size_t i=0; i<nna.size(); i++) {
			n += nna[i];
		}
	} else if (fun == "min") {
		s = vmin(v, narm);
		stat = std::min(stat, s);
	} else if (fun == "max") {
		s = vmax(v, narm);		
		stat = std::max(stat, s);
	}
}


SpatDataFrame SpatRaster::global(std::string fun, bool narm) {

	SpatDataFrame out;
	std::vector<std::string> f {"sum", "mean", "min", "max"};
	if (std::find(f.begin(), f.end(), fun) == f.end()) {
		out.setError("not a valid function");
		return(out);
	}

	if (!hasValues()) {
		out.setError("SpatRaster has no values");
		return(out);
	}

	std::vector<double> stats(nlyr());
	std::vector<unsigned> n(nlyr());
	readStart();
	BlockSize bs = getBlockSize(2);
	for (size_t i=0; i<bs.n; i++) {
		std::vector<double> v =   readValues(bs.row[i], bs.nrows[i], 0, ncol());
		unsigned off = bs.nrows[i] * ncol() ;
		for (size_t lyr=0; lyr<nlyr(); lyr++) {
			unsigned offset = lyr * off;
			std::vector<double> vv = {  v.begin()+offset,  v.begin()+offset+off };
			do_stats(vv, fun, narm, stats[lyr], n[lyr]);
		}
	}
	readStop();


	if (fun=="mean") {
		for (size_t lyr=0; lyr<nlyr(); lyr++) {
			if (n[lyr] > 0) {
				stats[lyr] = stats[lyr] / n[lyr];
			} else {
				stats[lyr] = NAN;
			}
		}
	}


	out.add_column(stats, fun);
	return(out);
}

