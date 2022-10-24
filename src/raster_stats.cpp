// Copyright (c) 2018-2022  Robert J. Hijmans
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

std::map<double, unsigned long long> table(std::vector<double> &v) {
	std::map<double, unsigned long long> count;
	for_each( v.begin(), v.end(), [&count]( double val ){
			if(!std::isnan(val)) count[val]++;
		}
	);
	return count;
}


std::map<double, unsigned long long int> ctable(std::map<double, unsigned long long int> &x, std::map<double, unsigned long long int> &y) {
	for(auto p : y) {
		x[p.first] += p.second;
	}
	return(x);
}


std::vector<double> vtable(std::map<double, unsigned long long int> &x) {
	std::vector<std::vector<double>> out(2);
	for( auto p : x ) {
		out[0].push_back(p.first);
		out[1].push_back(p.second);
	}
	out[0].insert(out[0].end(), out[1].begin(), out[1].end());
	return out[0];
}



std::vector<std::vector<double>> SpatRaster::freq(bool bylayer, bool round, int digits, SpatOptions &opt) {
	std::vector<std::vector<double>> out;
	if (!hasValues()) return out;
	BlockSize bs = getBlockSize(opt);
	unsigned nc = ncol();
	unsigned nl = nlyr();
	if (!readStart()) {
		return(out);
	}

	if (bylayer) {
		out.resize(nl);
		std::vector<std::map<double, unsigned long long int>> tabs(nl);
		for (size_t i = 0; i < bs.n; i++) {
			unsigned nrc = bs.nrows[i] * nc;
			std::vector<double> v;
			readValues(v, bs.row[i], bs.nrows[i], 0, nc);
			if (round) {
				for(double& d : v) d = roundn(d, digits);
			}
			for (size_t lyr=0; lyr<nl; lyr++) {
				unsigned off = lyr*nrc;
				std::vector<double> vv(v.begin()+off, v.begin() + off + nrc);
				std::map<double, unsigned long long int> tab = table(vv);
				tabs[lyr] = ctable(tabs[lyr], tab);
			}
		}
		for (size_t lyr=0; lyr<nl; lyr++) {
			out[lyr] = vtable(tabs[lyr]);
		}
	} else {
		out.resize(1);
		std::map<double, long long unsigned> tabs;
		for (size_t i = 0; i < bs.n; i++) {
			std::vector<double> v;
			readValues(v, bs.row[i], bs.nrows[i], 0, nc);
			if (round) {
				for (double& d : v) d = roundn(d, digits);
			}
			std::map<double, long long unsigned> tab = table(v);
			tabs = ctable(tabs, tab);
		}
		out[0] = vtable(tabs);
	}
	readStop();
	return(out);
}


std::vector<size_t> SpatRaster::count(double value, bool bylayer, bool round, int digits, SpatOptions &opt) {
	std::vector<size_t> out;
	if (!hasValues()) return out;
	BlockSize bs = getBlockSize(opt);
	unsigned nc = ncol();
	unsigned nl = nlyr();
	if (!readStart()) {
		return(out);
	}

	if (bylayer) {
		out.resize(nl);
		for (size_t i = 0; i < bs.n; i++) {
			unsigned nrc = bs.nrows[i] * nc;
			std::vector<double> v;
			readValues(v, bs.row[i], bs.nrows[i], 0, nc);
			if (round) {
				for(double& d : v) d = roundn(d, digits);
			}
			if (std::isnan(value)) {
				for (size_t lyr=0; lyr<nl; lyr++) {
					unsigned off = lyr*nrc;
					out[lyr] += count_if(v.begin()+off, v.begin()+off+nrc,
							[](double d){return std::isnan(d);});
				}
			} else {
				for (size_t lyr=0; lyr<nl; lyr++) {
					unsigned off = lyr*nrc;
					out[lyr] += std::count(v.begin()+off, v.begin()+off+nrc, value);
				}
			}
		}
	} else {
		out.resize(1);
		for (size_t i = 0; i < bs.n; i++) {
			std::vector<double> v;
			readValues(v, bs.row[i], bs.nrows[i], 0, nc);
			if (round) {
				for (double& d : v) d = roundn(d, digits);
			}
			if (std::isnan(value)) {
				out[0] += count_if(v.begin(), v.end(),
								[](double d){return std::isnan(d);});
			} else {
				out[0] += std::count(v.begin(), v.end(), value);
			}
		}
	}
	readStop();
	return(out);
}



SpatRaster SpatRaster::quantile(std::vector<double> probs, bool narm, SpatOptions &opt) {

	SpatRaster out = geometry(1);
	size_t n = probs.size();

	if (n == 0) {
		out.setError("no probs");
		return out;
	} else if (nlyr() < 2) {
		out.setError("more than one layer needed to compute quantiles");
		return out;
	}

	double pmin = vmin(probs, false);
	double pmax = vmin(probs, false);
	if ((std::isnan(pmin)) || (std::isnan(pmax)) || (pmin < 0) || (pmax > 1)) {
		SpatRaster out = geometry(1);
		out.setError("intvalid probs");
		return out;
	}
	out = geometry(probs.size());
	out.source[0].names = double_to_string(probs, "q");
  	if (!hasValues()) { return out; }

	if (!readStart()) {
		out.setError(getError());
		return(out);
	}
  	if (!out.writeStart(opt, filenames())) {
		readStop();
		return out;
	}
	unsigned nl = nlyr();
	std::vector<double> v(nl);

	for (size_t i = 0; i < out.bs.n; i++) {
		std::vector<double> a;
		readBlock(a, out.bs, i);
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
		if (!out.writeBlock(b, i)) return out;
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


void unique_values(std::vector<double> &d, bool narm) {
	if (narm) {
		d.erase(std::remove_if(d.begin(), d.end(),
				[](const double& value) { return std::isnan(value); }), d.end());
		std::set<double> u { d.begin(), d.end()};
		std::copy(u.begin(), u.end(), d.begin());
		d.erase(d.begin()+u.size(), d.end());
	} else {
		size_t s = d.size();
		d.erase(std::remove_if(d.begin(), d.end(),
				[](const double& value) { return std::isnan(value); }), d.end());
		bool addNAN = s > d.size();
		std::set<double> u { d.begin(), d.end()};
		std::copy(u.begin(), u.end(), d.begin());
		d.erase(d.begin()+u.size(), d.end());
		if (addNAN) d.push_back(NAN);
	}
}


std::vector<std::vector<double>> SpatRaster::unique(bool bylayer, bool narm, SpatOptions &opt) {

	std::vector<std::vector<double>> out;
	if (!hasValues()) return out;

	constexpr double lowest_double = std::numeric_limits<double>::lowest();

	BlockSize bs = getBlockSize(opt);
	unsigned nc = ncol();
	unsigned nl = nlyr();
	if (!readStart()) {
		return(out);
	}
	out.resize(nl);

	if (nl == 1) bylayer = true;
	if (bylayer) {
		for (size_t i = 0; i < bs.n; i++) {
			unsigned n = bs.nrows[i] * nc;
			std::vector<double> v;
			readValues(v, bs.row[i], bs.nrows[i], 0, nc);
			for (size_t lyr=0; lyr<nl; lyr++) {
				unsigned off = lyr*n;
				out[lyr].insert(out[lyr].end(), v.begin()+off, v.begin()+off+n);
				unique_values(out[lyr], narm);
			}
		}
	} else {
		std::vector<std::vector<double>> temp;
		for (size_t i = 0; i < bs.n; i++) {
			unsigned n = bs.nrows[i] * nc;
			std::vector<std::vector<double>> m(n, std::vector<double>(nl));
			std::vector<double> v;
			readValues(v, bs.row[i], bs.nrows[i], 0, nc);
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

	if (fun=="sum") {
		if (narm) {
			for (size_t i=0; i<z.size(); i++) {
				if ((!std::isnan(z[i])) && (!std::isnan(v[i]))) {
					out[z[i]] += v[i];
				}
			}
		} else {
			for (size_t i=0; i<z.size(); i++) {
				if ((!std::isnan(z[i]))) {
					out[z[i]] += v[i];
				}
			}
		}
	} else if (fun=="mean") {
		if (narm) {
			for (size_t i=0; i<z.size(); i++) {
				if ((!std::isnan(z[i])) && (!std::isnan(v[i]))) {
					out[z[i]] += v[i];
					cnt[z[i]]++;
				}
			}
		} else {
			for (size_t i=0; i<z.size(); i++) {
				if ((!std::isnan(z[i]))) {
					out[z[i]] += v[i];
					cnt[z[i]]++;
				}
			}
		}
	} else if (fun=="min") {
		if (narm) {
			for (size_t i=0; i<z.size(); i++) {
				if ((!std::isnan(z[i])) && (!std::isnan(v[i]))) {
					out[z[i]] = std::min(out[z[i]], v[i]);
				}
			}
		} else {
			for (size_t i=0; i<z.size(); i++) {
				if ((!std::isnan(z[i]))) {
					out[z[i]] = std::min(out[z[i]], v[i]);
				}
			}
		}
	} else if (fun=="max") {
		if (narm) {
			for (size_t i=0; i<z.size(); i++) {
				if ((!std::isnan(z[i])) && (!std::isnan(v[i]))) {
					out[z[i]] = std::max(out[z[i]], v[i]);
				}
			}
		} else {
			for (size_t i=0; i<z.size(); i++) {
				if ((!std::isnan(z[i]))) {
					out[z[i]] = std::max(out[z[i]], v[i]);
				}
			}
		}
	}
}



SpatDataFrame SpatRaster::zonal_old(SpatRaster z, std::string fun, bool narm, SpatOptions &opt) {

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
	if (!compare_geom(z, false, true, opt.get_tolerance())) {
		out.setError("dimensions and/or extent do not match");
		return(out);
	}

	if (z.nlyr() > 1) {
		SpatOptions xopt(opt);
		std::vector<unsigned> lyr = {0};
		z = z.subset(lyr, xopt);
		out.addWarning("only the first zonal layer is used");
	}

	size_t nl = nlyr();
	std::vector<std::vector<double>> uq = z.unique(true, true, opt);
	std::vector<double> u = uq[0];
	double initv = 0;
	double posinf = std::numeric_limits<double>::infinity();
	double neginf = -posinf;
	if (fun == "max") initv = neginf;
	if (fun == "min") initv = posinf;
	std::vector<std::vector<double>> stats(nl, std::vector<double>(u.size(), initv));
	std::vector<std::vector<double>> cnt(nl, std::vector<double>(u.size(), 0));
	if (!readStart()) {
		out.setError(getError());
		return(out);
	}
	if (!z.readStart()) {
		out.setError(z.getError());
		return(out);
	}
	opt.ncopies = 6;
	BlockSize bs = getBlockSize(opt);
	for (size_t i=0; i<bs.n; i++) {
		std::vector<double> v, zv;
		readValues(v, bs.row[i], bs.nrows[i], 0, ncol());
		z.readValues(zv, bs.row[i], bs.nrows[i], 0, ncol());
		std::vector<double> zvr(zv.size());
		for (size_t j=0; j<zvr.size(); j++)	 {
			if (std::isnan(zv[j])) {
				zvr[j] = NAN;
			} else {
				for (size_t k=0; k<u.size(); k++) {
					if (zv[j] == u[k]) {
						zvr[j] = k;
						continue;
					}
				}
			}
		}
		zv.resize(0);
		unsigned off = bs.nrows[i] * ncol() ;
		if (nl > 1) {
			for (size_t lyr=0; lyr<nl; lyr++) {
				unsigned offset = lyr * off;
				std::vector<double> vx( v.begin()+offset,  v.begin()+offset+off);
				jointstats(u, vx, zvr, fun, narm, stats[lyr], cnt[lyr]);
			}
		} else {
			jointstats(u, v, zvr, fun, narm, stats[0], cnt[0]);
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

	else if (fun == "min") {
		for (size_t lyr=0; lyr<nlyr(); lyr++) {
			for (size_t j=0; j<u.size(); j++) {
				if (stats[lyr][j] == posinf) {
					stats[lyr][j] = NAN;
				}
			}
		}
	} else if (fun == "max") {
		for (size_t lyr=0; lyr<nlyr(); lyr++) {
			for (size_t j=0; j<u.size(); j++) {
				if (stats[lyr][j] == neginf) {
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




SpatDataFrame SpatRaster::zonal(SpatRaster z, std::string fun, bool narm, SpatOptions &opt) {

	SpatDataFrame out;
	std::vector<std::string> f {"sum", "mean", "min", "max", "isNA", "notNA"};
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
	if (!compare_geom(z, false, true, opt.get_tolerance(), true)) {
		out.setError(getError());
		return(out);
	}
	if (hasWarning()) {
		std::vector<std::string> w = getWarnings();
		for (size_t i=0; i<w.size(); i++) {
			out.addWarning(w[i]);
		}
	}

	if (z.nlyr() > 1) {
		SpatOptions xopt(opt);
		std::vector<unsigned> lyr = {0};
		z = z.subset(lyr, xopt);
		out.addWarning("only the first zonal layer is used");
	}

	double posinf = std::numeric_limits<double>::infinity();
	double neginf = -posinf;
	if (!readStart()) {
		out.setError(getError());
		return(out);
	}
	if (!z.readStart()) {
		out.setError(z.getError());
		return(out);
	}
	opt.ncopies = 6;
	BlockSize bs = getBlockSize(opt);

	size_t nl = nlyr();
	size_t nc = ncol();
	std::vector<std::map<double, double>> m(nl);
	std::vector<std::map<double, size_t>> cnt(nl);

	for (size_t i=0; i<bs.n; i++) {
		unsigned nrc = bs.nrows[i] * nc;
		std::vector<double> vv, zv;
		readValues(vv, bs.row[i], bs.nrows[i], 0, ncol());
		z.readValues(zv, bs.row[i], bs.nrows[i], 0, ncol());
		if (fun == "sum") {
			for (size_t j=0; j<nl; j++) {
				size_t off = j*nrc;
				std::vector<double> v(vv.begin()+off, vv.begin() + off + nrc);
				for (size_t k=0; k<nrc; k++) {
					if (std::isnan(zv[k])) {
						continue;
					}
					if (narm && std::isnan(v[k])) {
						if (m[j].find(zv[k]) == m[j].end()) {
							m[j][zv[k]] = 0;
							cnt[j][zv[k]] = 0;
						}
					} else if (m[j].find(zv[k]) == m[j].end()) {
						m[j][zv[k]] = v[k];
						cnt[j][zv[k]] = 1;
					} else {
						m[j][zv[k]] += v[k];
						cnt[j][zv[k]] = 1; // may be necessary if the first case was NAN
					}
				}
			}
		} else if (fun == "mean") {
			for (size_t j=0; j<nl; j++) {
				size_t off = j*nrc;
				std::vector<double> v(vv.begin()+off, vv.begin() + off + nrc);
				for (size_t k=0; k<nrc; k++) {
					if (std::isnan(zv[k])) {
						continue;
					}
					if (narm && std::isnan(v[k])) {
						if (m[j].find(zv[k]) == m[j].end()) {
							m[j][zv[k]] = 0;
							cnt[j][zv[k]] = 0;
						}
					} else if (m[j].find(zv[k]) == m[j].end()) {
						m[j][zv[k]] = v[k];
						cnt[j][zv[k]] = 1;
					} else {
						m[j][zv[k]] += v[k];
						cnt[j][zv[k]]++;
					}
				}
			}
		} else if (fun == "min") {
			for (size_t j=0; j<nl; j++) {
				size_t off = j*nrc;
				std::vector<double> v(vv.begin()+off, vv.begin() + off + nrc);
				for (size_t k=0; k<nrc; k++) {
					if (std::isnan(zv[k])) {
						continue;
					}
					if (narm && std::isnan(v[k])) {
						if (m[j].find(zv[k]) == m[j].end()) {
							m[j][zv[k]] = posinf;
							cnt[j][zv[k]] = 0;
						}
					} else if (m[j].find(zv[k]) == m[j].end()) {
						m[j][zv[k]] = v[k];
						cnt[j][zv[k]] = 1;
					} else {
						m[j][zv[k]] = std::min(v[k], m[j][zv[k]]);
						cnt[j][zv[k]] = 1;
					}
				}
			}
		} else if (fun == "max") {
			for (size_t j=0; j<nl; j++) {
				size_t off = j*nrc;
				std::vector<double> v(vv.begin()+off, vv.begin() + off + nrc);
				for (size_t k=0; k<nrc; k++) {
					if (std::isnan(zv[k])) {
						continue;
					}
					if (narm && std::isnan(v[k])) {
						if (m[j].find(zv[k]) == m[j].end()) {
							m[j][zv[k]] = neginf;
							cnt[j][zv[k]] = 0;
						}
					} else if (m[j].find(zv[k]) == m[j].end()) {
						m[j][zv[k]] = v[k];
						cnt[j][zv[k]] = 1;
					} else {
						m[j][zv[k]] = std::max(v[k], m[j][zv[k]]);
						cnt[j][zv[k]] = 1;
					}
				}
			}
		} else if (fun == "isNA") {
			for (size_t j=0; j<nl; j++) {
				size_t off = j*nrc;
				std::vector<double> v(vv.begin()+off, vv.begin() + off + nrc);
				for (size_t k=0; k<nrc; k++) {
					if (std::isnan(zv[k])) {
						continue;
					}
					if (!std::isnan(v[k])) {
						if (cnt[j].find(zv[k]) == cnt[j].end()) {
							cnt[j][zv[k]] = 0;
						}
					} else {
						if (cnt[j].find(zv[k]) == cnt[j].end()) {
							cnt[j][zv[k]] = 1;
						} else {
							cnt[j][zv[k]]++;
						}
					}
				}
			}
		} else if (fun == "notNA") {
			for (size_t j=0; j<nl; j++) {
				size_t off = j*nrc;
				std::vector<double> v(vv.begin()+off, vv.begin() + off + nrc);
				for (size_t k=0; k<nrc; k++) {
					if (std::isnan(zv[k])) {
						continue;
					}
					if (std::isnan(v[k])) {
						if (cnt[j].find(zv[k]) == cnt[j].end()) {
							cnt[j][zv[k]] = 0;
						}
					} else {
						if (cnt[j].find(zv[k]) == cnt[j].end()) {
							cnt[j][zv[k]] = 1;
						} else {
							cnt[j][zv[k]]++;
						}
					}
				}
			}
		}
	}
	readStop();
	z.readStop();

	std::vector<double> zone;
	std::vector<std::string> nms = getNames();

	if ((fun == "notNA") || (fun == "isNA")){
		size_t n = cnt[0].size();
		zone.reserve(n);
		for (size_t i=0; i<nl; i++) {
			std::vector<double> value;
			value.reserve(n);
			if (i==0) {
				for (auto& it : cnt[0]) {
					zone.push_back(it.first);
					value.push_back(it.second);
				}
				out.add_column(zone, "zone");
			} else {
				for (auto& it : cnt[i]) {
					value.push_back(it.second);
				}
			}
			out.add_column(value, nms[i]);
		}
	} else {
		size_t n = m[0].size();
		zone.reserve(n);
		for (size_t i=0; i<nl; i++) {
			std::vector<double> value;
			value.reserve(n);
			if (i==0) {
				for (auto& it : m[0]) {
					zone.push_back(it.first);
					value.push_back(it.second);
				}
				out.add_column(zone, "zone");
			} else {
				for (auto& it : m[i]) {
					value.push_back(it.second);
				}
			}
			size_t j = 0;
			if (fun == "mean") {
				for (auto& it : cnt[i]) {
					double d = (double)it.second;
					if (d > 0) {
						value[j] /= d;
					} else {
						value[j] = NAN;
					}
					j++;
				}
			} else {
				for (auto& it : cnt[i]) {
					if (it.second == 0) {
						value[j] = NAN;
					}
					j++;
				}
			}
			out.add_column(value, nms[i]);
		}
	}
//	std::vector<std::string> nms = getNames();
//	for (size_t i=0; i<nlyr(); i++) {
//		out.add_column(stats[i], nms[i]);
//	}
	return(out);
}
