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
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURP0OSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with spat. If not, see <http://www.gnu.org/licenses/>.

#include "spatRaster.h"
#include "recycle.h"
#include "vecmath.h"


SpatRaster SpatRaster::apply(std::vector<unsigned> ind, std::string fun, bool narm, std::vector<std::string> nms, SpatOptions &opt) {


	recycle(ind, nlyr());
	std::vector<unsigned> ui = vunique(ind);
	unsigned nl = ui.size();		
	SpatRaster out = geometry(nl);
	recycle(nms, nl);
	out.setNames(nms);

	std::vector<std::string> f {"sum", "mean", "min", "max", "prod", "any", "all"};
	if (std::find(f.begin(), f.end(), fun) == f.end()) {
		out.setError("unknown apply function");
		return out;
	}
	
	if (!hasValues()) return(out);
 	if (!out.writeStart(opt)) { return out; }
	BlockSize bs = getBlockSize(8);
	readStart();
	std::vector<std::vector<double>> v(nl);
	std::vector<unsigned> ird(ind.size()); 
	std::vector<unsigned> jrd(ind.size()); 
	for (size_t i=0; i<nl; i++) {
		for (size_t j=0; j<ind.size(); j++) {
			if (ui[i] == ind[j]) {
				v[i].push_back(0);
				ird[j] = i;
				jrd[j] = v[i].size()-1;				
			}
		}
	}
	
	for (size_t i=0; i<bs.n; i++) {
        std::vector<double> a = readBlock(bs, i);
		unsigned nc = out.bs.nrows[i] * out.ncol();
		std::vector<double> b(nc * nl);
		for (size_t j=0; j<nc; j++) {
			for (size_t k=0; k<ird.size(); k++) {
				v[ird[k]][jrd[k]] = a[j+k*nc];
			}
			for (size_t k=0; k<ui.size(); k++) {
				size_t off = k * nc;
				if (fun == "sum") {
					b[off+j] = vsum(v[k], narm);
				} else if (fun == "mean") {
					b[off+j] = vmean(v[k], narm);
				} else if (fun == "prod") {
					b[off+j] = vprod(v[k], narm);
				} else if (fun == "min") {
					b[off+j] = vmin(v[k], narm);
				} else if (fun == "max") {
					b[off+j] = vmax(v[k], narm);
				} else if (fun == "any") {
					b[off+j] = vany(v[k], narm);
				} else if (fun == "all") {
					b[off+j] = vall(v[k], narm);
				}
			}
		
		}
		if (!out.writeValues(b, out.bs.row[i], out.bs.nrows[i], 0, ncol())) return out;
	}
	readStop();
	out.writeStop();
	return(out);
}

