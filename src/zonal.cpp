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


void jointstats(const std::vector<double> &u, const std::vector<double> &v, const std::vector<double> &z, std::string fun, bool narm, std::vector<double>& out, std::vector<double> &cnt) {

	std::vector<double> cmp;
	std::vector<bool> done(z.size(), false);
	//recycle(v, z);

	for (size_t j=0; j<u.size(); j++) {
		cmp.resize(0);
		for (size_t k=0; k<v.size(); k++) {
			if (j==0) {
				if (std::isnan(z[k])) {
					done[k] = true;
				}
			}
			if (!done[k]) {
				if (z[k] == u[j]) {
					if (!(narm & std::isnan(v[k]))) {
						cmp.push_back(v[k]);
						done[k]=true;
					}
				}
			}
		}

		if (fun=="sum") {
			double s = vsum(cmp, narm);
			out[j] = s + out[j];
		} else if (fun=="mean") {
			double s = vsum(cmp, narm);
			if (narm) {
				for (size_t k=1; k<cmp.size(); k++) {
					cnt[j] += !std::isnan(cmp[k]);
				}
			} else {
				cnt[j] += cmp.size();
			}
			out[j] = s + out[j];
		} else if (fun == "min") {
			double m = vmin(cmp, narm);
			if (narm) {
				if (!std::isnan(m)) {
					if (cnt[j] == 0) {
						out[j] = m;
						cnt[j] = 1;
					} else {
						out[j] = std::min(m, out[j]);
					}
				}
			} else {
				if (cnt[j] == 0) {
					out[j] = m;
					cnt[j] = 1;
				} else {
					out[j] = std::min(m, out[j]);
				}
			}
		} else if (fun == "max") {
			double m = vmax(cmp, narm);
			if (narm) {
				if (!std::isnan(m)) {
					if (cnt[j] == 0) {
						out[j] = m;
						cnt[j] = 1;
					} else {
						out[j] = std::max(m, out[j]);
					}
				}
			} else {
				if (cnt[j] == 0) {
					out[j] = m;
					cnt[j] = 1;
				} else {
					out[j] = std::max(m, out[j]);
				}
			}
		}


	}
}




SpatDataFrame SpatRaster::zonal(SpatRaster z, std::string fun, bool narm) {

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
	if (!z.hasValues()) {
		out.setError("zonal SpatRaster has no values");
		return(out);
	}
	if (!compare_geom(z, false, true)) {
		out.setError("dimensions and/or extent do not match");
		return(out);
	}

	if (z.nlyr() > 1) {
		SpatOptions opt;
		std::vector<unsigned> lyr = {0};
		z = z.subset(lyr, opt);
	}

	std::vector<std::vector<double>> uq = z.unique(true);
	std::vector<double> u = uq[0];
	std::vector<std::vector<double>> stats(nlyr(), std::vector<double>(u.size()));
	std::vector<std::vector<double>> cnt(nlyr(), std::vector<double>(u.size()));

	readStart();
	z.readStart();
	BlockSize bs = getBlockSize(8);
	for (size_t i=0; i<bs.n; i++) {
		std::vector<double> v =   readValues(bs.row[i], bs.nrows[i], 0, ncol());
		std::vector<double> zv = z.readValues(bs.row[i], bs.nrows[i], 0, ncol());
		unsigned off = bs.nrows[i] * ncol() ;
		for (size_t lyr=0; lyr<nlyr(); lyr++) {
			unsigned offset = lyr * off;
			std::vector<double> vv = {  v.begin()+offset,  v.begin()+offset+off };
			jointstats(u, vv, zv, fun, narm, stats[lyr], cnt[lyr]);
		}
	}
	readStop();
	z.readStop();

	if (fun=="mean") {
		for (size_t lyr=0; lyr<nlyr(); lyr++) {
			for (size_t j=0; j<u.size(); j++) {
				if (cnt[lyr][j] > 0) {
					stats[lyr][j] = stats[lyr][j] / (1+cnt[lyr][j]);
				} else {
					stats[lyr][j] = NAN;
				}
			}
		}
	}

	out.add_column(u, "zone");
	std::vector<std::string> nms = getNames();
	for (size_t i=0; i<nlyr(); i++) {
		out.add_column(stats[i], nms[i]);
	}
	return(out);
}



