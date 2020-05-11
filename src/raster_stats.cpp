// Copyright (c) 2018-2020  Robert J. Hijmans
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
#include <cmath>
#include <algorithm>
#include <map>

#include "vecmath.h"
#include "math_utils.h"
#include "string_utils.h"

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


static inline double interpolate(double x, double y1, double y2, unsigned x1, unsigned x2) {
	double denom = (x2-x1);
	return y1 + (x-x1) * (y2-y1)/denom;
}


static inline std::vector<double> vquantile(std::vector<double> v, const std::vector<double>& probs, bool narm) {
	size_t n = v.size();
    if (n==0) {
        return std::vector<double>(probs.size(), NAN);
    }
    if (n == 1) {
        return std::vector<double>(probs.size(), v[0]);
    }
	na_omit(v);
	if ((!narm) & (v.size() < n)) {
        return std::vector<double>(probs.size(), NAN);
	}
	n = v.size();
    std::sort(v.begin(), v.end());

	size_t pn = probs.size();
	std::vector<double> q(pn);

    for (size_t i = 0; i < pn; ++i) {
		double x = probs[i] * (n-1);
		unsigned x1 = std::floor(x);
		unsigned x2 = std::ceil(x);
		if (x1 == x2) {
			q[i] = v[x1];
		} else {
			q[i] = interpolate(x, v[x1], v[x2], x1, x2);
		}
    }
    return q;
}


SpatRaster SpatRaster::quantile(std::vector<double> probs, bool narm, SpatOptions &opt) {
	size_t n = probs.size();
	if (n == 0) {
		SpatRaster out = geometry(1);
		out.setError("no probs");
		return out;
	}
	double pmin = vmin(probs, false);
	double pmax = vmin(probs, false);
	if ((std::isnan(pmin)) | (std::isnan(pmax)) | (pmin < 0) | (pmax > 1)) {
		SpatRaster out = geometry(1);
		out.setError("intvalid probs");
		return out;
	}
	SpatRaster out = geometry(probs.size());
	out.source[0].names = double_to_string(probs, "q");
  	if (!hasValues()) { return out; }

  	if (!out.writeStart(opt)) { return out; }
	readStart();
	unsigned nl = nlyr();
	std::vector<double> v(nl);

	for (size_t i = 0; i < out.bs.n; i++) {
		std::vector<double> a = readBlock(out.bs, i);
		unsigned nc = out.bs.nrows[i] * out.ncol();
		std::vector<double> b(nc * n);
		for (size_t j=0; j<nc; j++) {
			for (size_t k=0; k<nl; k++) {
				v[k] = a[j+k*nc];
			}
			std::vector<double> p = vquantile(v, probs, narm);
			for (size_t k=0; k<n; k++) {
				b[j+(k*nc)] = p[k];
			}
		}
		if (!out.writeValues(b, out.bs.row[i], out.bs.nrows[i], 0, ncol())) return out;
	}
	out.writeStop();
	readStop();
	return(out);
}





void unique_values_alt(std::vector<double> &d) {
	d.erase(std::remove_if(d.begin(), d.end(),
            [](const double& value) { return std::isnan(value); }), d.end());
	std::sort(d.begin(), d.end());
	d.erase(std::unique(d.begin(), d.end()), d.end());
}


void unique_values(std::vector<double> &d) {
	d.erase(std::remove_if(d.begin(), d.end(),
            [](const double& value) { return std::isnan(value); }), d.end());
	std::set<double> u { d.begin(), d.end()};
	std::copy(u.begin(), u.end(), d.begin());
	d.erase(d.begin()+u.size(), d.end());
}


std::vector<std::vector<double>> SpatRaster::unique(bool bylayer) {

	std::vector<std::vector<double>> out;
	if (!hasValues()) return out;

	constexpr double lowest_double = std::numeric_limits<double>::lowest();

	BlockSize bs = getBlockSize(4);
	unsigned nc = ncol();
	unsigned nl = nlyr();
	readStart();
	out.resize(nl);

	if (nl == 1) bylayer = true;
	if (bylayer) {
		for (size_t i = 0; i < bs.n; i++) {
			unsigned n = bs.nrows[i] * nc;
			std::vector<double> v = readValues(bs.row[i], bs.nrows[i], 0, nc);
			for (size_t lyr=0; lyr<nl; lyr++) {
				unsigned off = lyr*n;
				out[lyr].insert(out[lyr].end(), v.begin()+off, v.begin()+off+n);
				unique_values(out[lyr]);
			}
		}
	} else {
		std::vector<std::vector<double>> temp;
		for (size_t i = 0; i < bs.n; i++) {
			unsigned n = bs.nrows[i] * nc;
			std::vector<std::vector<double>> m(n, std::vector<double>(nl));
			std::vector<double> v = readValues(bs.row[i], bs.nrows[i], 0, nc);
			for (size_t j = 0; j < v.size(); j++) {
				if (std::isnan(v[j])) v[j] = lowest_double;
			}

			for (size_t lyr=0; lyr<nl; lyr++) {
				unsigned off = lyr*n;
				for (size_t j=0; j<n; j++) {
					m[j][lyr] = v[off+j];
				}
			}
			std::sort(m.begin(), m.end());
			m.erase(std::unique(m.begin(), m.end()), m.end());
			for (size_t j=0; j<m.size(); j++) {
				temp.insert(temp.end(), m[j]);
			}
		}
		std::sort(temp.begin(), temp.end());
		temp.erase(std::unique(temp.begin(), temp.end()), temp.end());
		for (size_t i = 0; i < temp.size(); i++) {
			for (size_t j = 0; j < temp[0].size(); j++) {
				out[j].resize(temp.size());
				if (temp[i][j] == lowest_double) {
					out[j][i] = NAN;
				} else {
					out[j][i] = temp[i][j];
				}
			}
		}
	}
	readStop();

	return(out);
}




void jointstats(const std::vector<double> &u, const std::vector<double> &v, const std::vector<double> &z, std::string fun, bool narm, std::vector<double>& out, std::vector<double> &cnt) {

	std::vector<double> cmp;
	//recycle(v, z);

	for (size_t j=0; j<u.size(); j++) {
		cmp.resize(0);
		cmp.reserve(v.size() / u.size());

		for (size_t k=0; k<v.size(); k++) {
			if (z[k] == u[j]) {
				if (!(narm & std::isnan(v[k]))) {
					cmp.push_back(v[k]);
				}
			}
		}
		if (cmp.size() == 0) continue;
		if (fun=="sum") {
			double s = vsum(cmp, narm);
			out[j] = s + out[j];
		} else if (fun=="mean") {
			double s = vsum(cmp, narm);
			if (narm) {
				for (size_t k=0; k<cmp.size(); k++) {
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
	BlockSize bs = getBlockSize(12);
	for (size_t i=0; i<bs.n; i++) {
		std::vector<double> v =    readValues(bs.row[i], bs.nrows[i], 0, ncol());
		std::vector<double> zv = z.readValues(bs.row[i], bs.nrows[i], 0, ncol());
		unsigned off = bs.nrows[i] * ncol() ;
		for (size_t lyr=0; lyr<nlyr(); lyr++) {
			unsigned offset = lyr * off;
			std::vector<double> vx( v.begin()+offset,  v.begin()+offset+off);
			std::vector<double> vz(zv.begin()+offset, zv.begin()+offset+off);
			jointstats(u, vx, vz, fun, narm, stats[lyr], cnt[lyr]);
		}
	}
	readStop();
	z.readStop();

	if (fun=="mean") {
		for (size_t lyr=0; lyr<nlyr(); lyr++) {
			for (size_t j=0; j<u.size(); j++) {
				if (cnt[lyr][j] > 0) {
					stats[lyr][j] = stats[lyr][j] / cnt[lyr][j];
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



