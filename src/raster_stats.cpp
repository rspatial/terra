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

#include "spatRaster.h"
#include <limits>
#include <set>
#include <cmath>
#include <algorithm>
#include <map>

#include "vecmath.h"
#include "vecmathse.h"

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

	double pmin = min_se(probs, 0, probs.size());
	double pmax = max_se(probs, 0, probs.size());
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


std::vector<std::vector<double>> SpatRaster::unique(bool bylayer, double digits, bool narm, SpatOptions &opt) {

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
			if (!std::isnan(digits)) {
				int dig = digits;
				for(double& d : v) d = roundn(d, dig);
			}
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
			if (!std::isnan(digits)) {
				int dig = digits;
				for (size_t j = 0; j < v.size(); j++) {
					if (std::isnan(v[j])) {
						v[j] = lowest_double;
					} else {
						v[j] = roundn(v[j], dig);
					}
				}
			} else {
				for (size_t j = 0; j < v.size(); j++) {
					if (std::isnan(v[j])) v[j] = lowest_double;
				}
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


/*
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


*/

void zonalsum(const std::vector<double> &v, 
			const std::vector<double> &zv, 
			std::vector<std::map<double, double>> &m, 
			std::vector<std::map<double, size_t>> &cnt, 
			size_t nl, unsigned &nrc, bool narm){
				
	for (size_t i=0; i<nl; i++) {
		size_t off = i*nrc;
		for (size_t j=0; j<nrc; j++) {
			double z=zv[j];
			size_t k=j+off;
			if (std::isnan(z)) continue;
			if (narm && std::isnan(v[k])) {
				if (m[i].find(z) == m[i].end()) {
					m[i][z] = 0;
					cnt[i][z] = 0;
				}
			} else if (m[i].find(z) == m[i].end()) {
				m[i][z] = v[k];
				cnt[i][z] = 1;
			} else {
				m[i][z] += v[k];
				cnt[i][z] = 1; // may be necessary if the first case was NAN
			}
		}
	}
}


void zonalsumgroup(const std::vector<double> &v, 
		const std::vector<double> &zv, const std::vector<double> &gv, 
		std::vector<std::map<double, std::map<double, double>>> &m, 
		std::vector<std::map<double, std::map<double, size_t>>> &cnt, 
		size_t nl, unsigned &nrc, bool narm){

	for (size_t i=0; i<nl; i++) {
		size_t off = i*nrc;
		for (size_t j=0; j<nrc; j++) {
			double z = zv[j];
			size_t k = j+off;
			if (std::isnan(z) || std::isnan(gv[j])) continue;
			size_t g = gv[j];						
			if (narm && std::isnan(v[k])) {
				if (m[i].find(g) == m[i].end() || m[i][g].find(z) == m[i][g].end()) {
					m[i][g][z] = 0;
					cnt[i][g][z] = 0;
				}
			} else if (m[i].find(g) == m[i].end() || m[i][g].find(z) == m[i][g].end()) {
				m[i][g][z] = v[k];
				cnt[i][g][z] = 1;
			} else {
				m[i][g][z] += v[k];
				cnt[i][g][z] = 1; // may be necessary if the first case was NAN
			}
		}
	}
}


void zonalmean(const std::vector<double> &v, 
		const std::vector<double> &zv, 
		std::vector<std::map<double, double>> &m, 
		std::vector<std::map<double, size_t>> &cnt, 
		size_t nl, unsigned &nrc, bool narm){

	for (size_t i=0; i<nl; i++) {
		size_t off = i*nrc;
		for (size_t j=0; j<nrc; j++) {
			double z = zv[j];
			size_t k = j+off;
			if (std::isnan(z)) continue;
			if (narm && std::isnan(v[k])) {
				if (m[i].find(z) == m[i].end()) {
					m[i][z] = 0;
					cnt[i][z] = 0;
				}
			} else if (m[i].find(z) == m[i].end()) {
				m[i][z] = v[k];
				cnt[i][z] = 1;
			} else {
				m[i][z] += v[k];
				cnt[i][z]++;
			}
		}
	}
}


void zonalmeangroup(const std::vector<double> &v, 
		const std::vector<double> &zv, const std::vector<double> &gv, 
		std::vector<std::map<double, std::map<double, double>>> &m, 
		std::vector<std::map<double, std::map<double, size_t>>> &cnt, 
		size_t nl, unsigned &nrc, bool narm){

	for (size_t i=0; i<nl; i++) {
		size_t off = i*nrc;
		for (size_t j=0; j<nrc; j++) {
			double z = zv[j];
			size_t k = j+off;
			if (std::isnan(z) || std::isnan(gv[j])) continue;
			size_t g = gv[j];

			if (narm && std::isnan(v[k])) {
				if (m[i].find(g) == m[i].end() || m[i][g].find(z) == m[i][g].end()) {
					m[i][g][z] = 0;
					cnt[i][g][z] = 0;
				}
			} else if (m[i].find(g) == m[i].end() || m[i][g].find(z) == m[i][g].end()) {
				m[i][g][z] = v[k];
				cnt[i][g][z] = 1;
			} else {
				m[i][g][z] += v[k];
				cnt[i][g][z]++;
			}
		}
	}
}


void zonalmin(const std::vector<double> &v, 
		const std::vector<double> &zv, 
		std::vector<std::map<double, double>> &m, 
		std::vector<std::map<double, size_t>> &cnt, 
		size_t nl, unsigned &nrc, bool narm){

	double posinf = std::numeric_limits<double>::infinity();

	for (size_t i=0; i<nl; i++) {
		size_t off = i*nrc;
		for (size_t j=0; j<nrc; j++) {
			double z = zv[j];			
			size_t k = j+off;
			if (std::isnan(z)) continue;
			if (narm && std::isnan(v[k])) {
				if (m[i].find(z) == m[i].end()) {
					m[i][z] = posinf;
					cnt[i][z] = 0;
				}
			} else if (m[i].find(z) == m[i].end()) {
				m[i][z] = v[k];
				cnt[i][z] = 1;
			} else {
				m[i][z] = std::min(v[k], m[i][z]);
				cnt[i][z] = 1;
			}
		}
	}
}

void zonalmingroup(const std::vector<double> &v, 
		const std::vector<double> &zv, const std::vector<double> &gv, 
		std::vector<std::map<double, std::map<double, double>>> &m, 
		std::vector<std::map<double, std::map<double, size_t>>> &cnt, 
		size_t nl, unsigned &nrc, bool narm){

	double posinf = std::numeric_limits<double>::infinity();

	for (size_t i=0; i<nl; i++) {
		size_t off = i*nrc;
		for (size_t j=0; j<nrc; j++) {
			double z = zv[j];
			size_t k = j+off;
			if (std::isnan(z) || std::isnan(gv[j])) continue;
			size_t g = gv[j];
			if (narm && std::isnan(v[k])) {
				if (m[i].find(g) == m[i].end() || m[i][g].find(z) == m[i][g].end()) {
					m[i][g][z] = posinf;
					cnt[i][g][z] = 0;
				}
			} else if (m[i].find(g) == m[i].end() || m[i][g].find(z) == m[i][g].end()) {
				m[i][g][z] = v[k];
				cnt[i][g][z] = 1;
			} else {
				m[i][g][z] = std::min(v[k], m[i][g][z]);
				cnt[i][g][z] = 1;
			}
		}
	}
}


void zonalmax(const std::vector<double> &v, 
		const std::vector<double> &zv, 
		std::vector<std::map<double, double>> &m, 
		std::vector<std::map<double, size_t>> &cnt, 
		size_t nl, unsigned &nrc, bool narm){

	double neginf = -std::numeric_limits<double>::infinity();

	for (size_t i=0; i<nl; i++) {
		size_t off = i*nrc;
		for (size_t j=0; j<nrc; j++) {
			double z = zv[j];
			size_t k = j+off;
			if (std::isnan(z)) continue;
			if (narm && std::isnan(v[k])) {
				if (m[i].find(z) == m[i].end()) {
					m[i][z] = neginf;
					cnt[i][z] = 0;
				}
			} else if (m[i].find(z) == m[i].end()) {
				m[i][z] = v[k];
				cnt[i][z] = 1;
			} else {
				m[i][z] = std::max(v[k], m[i][z]);
				cnt[i][z] = 1;
			}
		}
	}
}


void zonalmaxgroup(const std::vector<double> &v, 
		const std::vector<double> &zv, const std::vector<double> &gv, 
		std::vector<std::map<double, std::map<double, double>>> &m, 
		std::vector<std::map<double, std::map<double, size_t>>> &cnt, 
		size_t nl, unsigned &nrc, bool narm){

	double neginf = -std::numeric_limits<double>::infinity();

	for (size_t i=0; i<nl; i++) {
		size_t off = i*nrc;
		for (size_t j=0; j<nrc; j++) {
			double z = zv[j];
			size_t k = j+off;
			if (std::isnan(z) || std::isnan(gv[j])) continue;
			size_t g = gv[j];
			if (narm && std::isnan(v[k])) {
				if (m[i].find(g) == m[i].end() || m[i][g].find(z) == m[i][g].end()) {
					m[i][g][z] = neginf;
					cnt[i][g][z] = 0;
				}
			} else if (m[i].find(g) == m[i].end() || m[i][g].find(z) == m[i][g].end()) {
				m[i][g][z] = v[k];
				cnt[i][g][z] = 1;
			} else {
				m[i][g][z] = std::max(v[k], m[i][g][z]);
				cnt[i][g][z] = 1;
			}
		}
	}
}


void zonalisna(const std::vector<double> &v, 
		const std::vector<double> &zv, 
		std::vector<std::map<double, double>> &m, 
		std::vector<std::map<double, size_t>> &cnt, 
		size_t nl, unsigned &nrc, bool narm){

	for (size_t i=0; i<nl; i++) {
		size_t off = i*nrc;
		for (size_t j=0; j<nrc; j++) {
			double z = zv[j];			
			size_t k = j+off;
			if (std::isnan(z)) continue;
			if (!std::isnan(v[k])) {
				if (cnt[i].find(z) == cnt[i].end()) {
					cnt[i][z] = 0;
				}
			} else {
				if (cnt[i].find(z) == cnt[i].end()) {
					cnt[i][z] = 1;
				} else {
					cnt[i][z]++;
				}
			}
		}
	}
}


void zonalisnagroup(const std::vector<double> &v, 
		const std::vector<double> &zv, const std::vector<double> &gv, 
		std::vector<std::map<double, std::map<double, double>>> &m, 
		std::vector<std::map<double, std::map<double, size_t>>> &cnt, 
		size_t nl, unsigned &nrc, bool narm){

	for (size_t i=0; i<nl; i++) {
		size_t off = i*nrc;
		for (size_t j=0; j<nrc; j++) {
			double z = zv[j];
			size_t k = j+off;
			if (std::isnan(z) || std::isnan(gv[j])) continue;
			size_t g = gv[j];

			if (!std::isnan(v[k])) {
				if (cnt[i].find(g) == cnt[i].end() || cnt[i][g].find(z) == cnt[i][g].end()) {
					cnt[i][g][z] = 0;
				}
			} else {
				if (cnt[i].find(g) == cnt[i].end() || cnt[i][g].find(z) == cnt[i][g].end()) {
					cnt[i][g][z] = 1;
				} else {
					cnt[i][g][z]++;
				}
			}
		}
	}
}
			

void zonalnotna(const std::vector<double> &v, 
		const std::vector<double> &zv, 
		std::vector<std::map<double, double>> &m, 
		std::vector<std::map<double, size_t>> &cnt, 
		size_t nl, unsigned &nrc, bool narm){


	for (size_t i=0; i<nl; i++) {
		size_t off = i*nrc;
		for (size_t j=0; j<nrc; j++) {
			double z = zv[j];
			size_t k = j+off;
			if (std::isnan(z)) continue;
			if (std::isnan(v[k])) {
				if (cnt[i].find(z) == cnt[i].end()) {
					cnt[i][z] = 0;
				}
			} else {
				if (cnt[i].find(z) == cnt[i].end()) {
					cnt[i][z] = 1;
				} else {
					cnt[i][z]++;
				}
			}
		}
	}
}

void zonalnotnagroup(const std::vector<double> &v, 
		const std::vector<double> &zv, const std::vector<double> &gv, 
		std::vector<std::map<double, std::map<double, double>>> &m, 
		std::vector<std::map<double, std::map<double, size_t>>> &cnt, 
		size_t nl, unsigned &nrc, bool narm) {

	for (size_t i=0; i<nl; i++) {
		size_t off = i*nrc;
		for (size_t j=0; j<nrc; j++) {
			double z = zv[j];
			size_t k = j+off;
			if (std::isnan(z) || std::isnan(gv[j])) continue;
			size_t g = gv[j];
			if (std::isnan(v[k])) {
				if (cnt[i].find(g) == cnt[i].end() || cnt[i][g].find(z) == cnt[i][g].end()) {
					cnt[i][g][z] = 0;
				}
			} else {
				if (cnt[i].find(g) == cnt[i].end() || cnt[i][g].find(z) == cnt[i][g].end()) {
					cnt[i][g][z] = 1;
				} else {
					cnt[i][g][z]++;
				}
			}
		}
	}
}


SpatDataFrame SpatRaster::zonal(SpatRaster z, SpatRaster g, std::string fun, bool narm, SpatOptions &opt) {

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
	bool groups = false;
	std::vector<size_t> groupid;
	if (g.hasValues()) {
		if (g.nlyr() > 1) {
			SpatOptions xopt(opt);
			g = g.subset({0}, xopt);
			out.addWarning("only the first grouping layer is used");
		}
		groups = true;
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
		// this is not very efficient. Should deal with multiple z layers below.
		SpatOptions xopt(opt);
		SpatDataFrame spout;
		std::vector<std::string> nms = z.getNames();
		make_unique_names(nms);
		z.setNames(nms);
		
		for (unsigned i=0; i<z.nlyr(); i++) {
			std::vector<unsigned> lyr = {i};
			SpatRaster zz = z.subset(lyr, xopt);
			SpatDataFrame spd = zonal(zz, g, fun, narm, xopt);
			std::vector<long> id(spd.nrow(), i);
			spd.add_column(id, "zlyr");
			if (i == 0) {
				spout = spd;
			} else {
				spout.rbind(spd);
			}
		}
		return spout;
	}

	if (!readStart()) {
		out.setError(getError());
		return(out);
	}
	if (!z.readStart()) {
		out.setError(z.getError());
		return(out);
	}
	opt.ncopies = 6;
	std::vector<SpatCategories> cats;
	if (groups) {
		if (!g.readStart()) {
			out.setError(g.getError());
			return(out);
		}
		opt.ncopies = 8;
	}
	
	BlockSize bs = getBlockSize(opt);

	size_t nl = nlyr();
	size_t nc = ncol();
	std::vector<std::map<double, double>> m;
	std::vector<std::map<double, size_t>> cnt;
	std::vector<std::map<double, std::map<double, double>>> gm;
	std::vector<std::map<double, std::map<double, size_t>>> gcnt;
	if (groups) {
		gm.resize(nl);
		gcnt.resize(nl);
	} else {
		m.resize(nl);
		cnt.resize(nl);
	}


	for (size_t i=0; i<bs.n; i++) {
		unsigned nrc = bs.nrows[i] * nc;
		std::vector<double> vv, zv, gv;
		readValues(vv, bs.row[i], bs.nrows[i], 0, ncol());
		z.readValues(zv, bs.row[i], bs.nrows[i], 0, ncol());
		if (groups) {
			g.readValues(gv, bs.row[i], bs.nrows[i], 0, ncol());
			if (fun == "sum") {
				zonalsumgroup(vv, zv, gv, gm, gcnt, nl, nrc, narm);				
			} else if (fun == "mean") {
				zonalmeangroup(vv, zv, gv, gm, gcnt, nl, nrc, narm);
			} else if (fun == "min") {
				zonalmingroup(vv, zv, gv, gm, gcnt, nl, nrc, narm);
			} else if (fun == "max") {
				zonalmaxgroup(vv, zv, gv, gm, gcnt, nl, nrc, narm);	
			} else if (fun == "isNA") {
				zonalisnagroup(vv, zv, gv, gm, gcnt, nl, nrc, narm);	
			} else if (fun == "notNA") {
				zonalnotnagroup(vv, zv, gv, gm, gcnt, nl, nrc, narm);	
			}
		} else {
			if (fun == "sum") {
				zonalsum(vv, zv, m, cnt, nl, nrc, narm);
			} else if (fun == "mean") {
				zonalmean(vv, zv, m, cnt, nl, nrc, narm);
			} else if (fun == "min") {
				zonalmin(vv, zv, m, cnt, nl, nrc, narm);			
			} else if (fun == "max") {
				zonalmax(vv, zv, m, cnt, nl, nrc, narm);			
			} else if (fun == "isNA") {
				zonalisna(vv, zv, m, cnt, nl, nrc, narm);			
			} else if (fun == "notNA") {
				zonalnotna(vv, zv, m, cnt, nl, nrc, narm);			
			}
		}
	}
	readStop();
	z.readStop();
	if (groups) g.readStop();
	
	std::vector<double> zone;
	std::vector<std::string> nms = getNames();

	if (groups) {
		size_t n1 = gm.size();
		size_t n2 = gm[0].size();
		size_t n = n1*n2;
		std::vector<double> layer;
		std::vector<double> group;
		std::vector<double> value;
		std::vector<size_t> cnter;
		layer.reserve(n);
		group.reserve(n);
		zone.reserve(n);
		value.reserve(n);			
		cnter.reserve(n);
		if ((fun == "notNA") || (fun == "isNA")){
			for (size_t i=0; i<nl; i++) {
				for (auto& it1:gcnt[i]) {
					std::map<double, size_t> mcnt = it1.second;
					for (auto& it2:mcnt) {
						layer.push_back(i);
						group.push_back(it1.first);
						zone.push_back(it2.first);
						value.push_back(it2.second);
					}
				}
			}
		} else {	
			for (size_t i=0; i<nl; i++) {
				for (auto& it1:gcnt[i]) {
					std::map<double, size_t> mcnt = it1.second;
					for (auto& it2:mcnt) {
						cnter.push_back(it2.second);
					}
				}
				for (auto& it1:gm[i]) {
					std::map<double, double> mg = it1.second;
					for (auto& it2:mg) {
						layer.push_back(i);
						group.push_back(it1.first);
						zone.push_back(it2.first);
						value.push_back(it2.second);
					}
				}
			} 
			if (fun == "mean")  {
				for (size_t i=0; i<cnter	.size(); i++) {
					if (cnter[i] > 0) {
						value[i] /= cnter[i];
					} else {
						value[i] = NAN;								
					}
				}
			} else {
				for (size_t i=0; i<cnter.size(); i++) {
					if (cnter[i] < 1) {
						value[i] = NAN;								
					}
				}
			}
		}
		out.add_column(layer, "layer");
		out.add_column(zone, "zone");
		out.add_column(group, "group");
		out.add_column(value, "value");
	} else {
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
	}
//	std::vector<std::string> nms = getNames();
//	for (size_t i=0; i<nlyr(); i++) {
//		out.add_column(stats[i], nms[i]);
//	}
	return(out);
}




SpatDataFrame SpatRaster::zonal_weighted(SpatRaster z, SpatRaster w, bool narm, SpatOptions &opt) {

	SpatDataFrame out;
	if (!hasValues()) {
		out.setError("SpatRaster has no values");
		return(out);
	}
	if (!z.hasValues()) {
		out.setError("zonal SpatRaster has no values");
		return(out);
	}
	if (!w.hasValues()) {
		out.setError("weights SpatRaster has no values");
		return(out);
	}
	if (!compare_geom(z, false, true, opt.get_tolerance(), true)) {
		out.setError(getError());
		return(out);
	}
	if (!compare_geom(w, false, true, opt.get_tolerance(), true)) {
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
	if (w.nlyr() > 1) {
		SpatOptions xopt(opt);
		std::vector<unsigned> lyr = {0};
		w = w.subset(lyr, xopt);
		out.addWarning("only the first weights layer is used");
	}

	if (!readStart()) {
		out.setError(getError());
		return(out);
	}
	if (!z.readStart()) {
		out.setError(z.getError());
		return(out);
	}
	if (!w.readStart()) {
		out.setError(z.getError());
		return(out);
	}
	opt.ncopies = 8;
	BlockSize bs = getBlockSize(opt);

	size_t nl = nlyr();
	size_t nc = ncol();
	std::vector<std::map<double, double>> m(nl);
	std::vector<std::map<double, size_t>> wsum(nl);

	for (size_t i=0; i<bs.n; i++) {
		unsigned nrc = bs.nrows[i] * nc;
		std::vector<double> vv, zv, wv;
		readValues(vv, bs.row[i], bs.nrows[i], 0, ncol());
		z.readValues(zv, bs.row[i], bs.nrows[i], 0, ncol());
		w.readValues(wv, bs.row[i], bs.nrows[i], 0, ncol());

		for (size_t j=0; j<nl; j++) {
			size_t off = j*nrc;
			std::vector<double> v(vv.begin()+off, vv.begin() + off + nrc);
			for (size_t k=0; k<nrc; k++) {
				if (std::isnan(zv[k]) || std::isnan(wv[k]) || (wv[k] < 0)) {
					continue;
				}
				if (narm && std::isnan(v[k])) {
					if (m[j].find(zv[k]) == m[j].end()) {
						m[j][zv[k]] = 0;
						wsum[j][zv[k]] = 0;
					}
				} else if (m[j].find(zv[k]) == m[j].end()) {
					m[j][zv[k]] = v[k] * wv[k];
					wsum[j][zv[k]] = wv[k];
				} else {
					m[j][zv[k]] += v[k] * wv[k];
					wsum[j][zv[k]] += wv[k];
				}
			}
		}
	}

	readStop();
	z.readStop();
	w.readStop();

	std::vector<double> zone;
	std::vector<std::string> nms = getNames();

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
		for (auto& it : wsum[i]) {
			double d = (double)it.second;
			if (d > 0) {
				value[j] /= d;
			} else {
				value[j] = NAN;
			}
			j++;
		}
		out.add_column(value, nms[i]);
	}
	return(out);
}


SpatDataFrame SpatRaster::zonal_poly(SpatVector x, std::string fun, bool weights, bool exact, bool touches,bool narm, SpatOptions &opt) {

	SpatDataFrame out;
	std::string gtype = x.type();
	if (gtype != "polygons") {
		out.setError("SpatVector must have polygon geometry");
		return out;
	}
	
	if (!hasValues()) {
		out.setError("raster has no values");
		return out;
	}

	if ((weights || exact)) {
		if ((fun != "mean") && (fun!="min") && (fun!="max")) {
			out.setError("fun should be 'min', 'max' or 'mean' when using weights/exact");
			return out;			
		}
	}

	if (!haveseFun(fun)) {
		out.setError("Unknown function");
		return out;
	}
	std::function<double(std::vector<double>&, size_t, size_t)> zfun;
	if (!getseFun(zfun, fun, narm)) {
		out.setError("Unknown function");
		return out;
	}

    unsigned nl = nlyr();
    unsigned ng = x.size();

	std::vector<std::vector<double>> zv(nl, std::vector<double>(ng));
	
    SpatRaster r = geometry(1);
    for (size_t i=0; i<ng; i++) {
		SpatGeom g = x.getGeom(i);
		SpatVector p(g);
		p.srs = x.srs;
		std::vector<double> cell, wgt;
		if (weights) {
			rasterizeCellsWeights(cell, wgt, p, opt);
		} else if (exact) {
			rasterizeCellsExact(cell, wgt, p, opt);
		} else {
			cell = rasterizeCells(p, touches, opt);
        }
		
		std::vector<std::vector<double>> e = extractCell(cell);
 		if ((weights || exact) && fun == "mean") {
			if (narm) {
				for (size_t j=0; j<nl; j++) {
					double wsum = 0;
					double vsum = 0;
					for (size_t k=0; k<e[j].size(); k++) {
						if (!std::isnan(e[j][k])) {
							wsum += wgt[k];
							vsum += (e[j][k] * wgt[k]);  
						}
					}
					zv[j][i] = vsum / wsum;
				}
			} else {
				for (size_t j=0; j<nl; j++) {
					double wsum = 0;
					double vsum = 0;
					for (size_t k=0; k<e[j].size(); k++) {
						wsum += wgt[k];
						vsum += (e[j][k] * wgt[k]);  
					}
					zv[j][i] = vsum / wsum;
				}
			}
		} else {
			for (size_t j=0; j<nl; j++) {
				zv[j][i] = zfun(e[j], 0, e[j].size());
			}
		}
	}
	std::vector<std::string> nms = getNames();	
	for (size_t j=0; j<nl; j++) {
		out.add_column(zv[j], nms[j]);
	}
	
	return out;
}


std::vector<double> tabfun(std::vector<double> x, std::vector<double> w) {
//	if (w.size() == 0) {
		std::map<double, long long unsigned> tab = table(x);
		return vtable(tab);
//	} else {
		
//	}
}


std::vector<std::vector<double>> SpatRaster::zonal_poly_table(SpatVector x, bool weights, bool exact, bool touches,bool narm, SpatOptions &opt) {

	std::vector<std::vector<double>> out;
	std::string gtype = x.type();
	if (gtype != "polygons") {
		setError("SpatVector must have polygon geometry");
		return out;
	}
	
	if (!hasValues()) {
		setError("raster has no values");
		return out;
	}

    unsigned nl = nlyr();
	if (nl > 1) {
		SpatOptions ops(opt);
		SpatRaster r = subset({0}, ops);
		out = r.zonal_poly_table(x, weights, exact, touches, narm, opt);
		addWarning("only the first layer of the raster is used");		
		return out;
	}

    unsigned ng = x.size();
	std::vector<std::vector<double>> zv(nl, std::vector<double>(ng));
	out.resize(ng);
    SpatRaster r = geometry(1);
    for (size_t i=0; i<ng; i++) {
		SpatGeom g = x.getGeom(i);
		SpatVector p(g);
		p.srs = x.srs;
		std::vector<double> cell, wgt;
		if (weights) {
			rasterizeCellsWeights(cell, wgt, p, opt);
		} else if (exact) {
			rasterizeCellsExact(cell, wgt, p, opt);
		} else {
			cell = rasterizeCells(p, touches, opt);
        }	
		std::vector<std::vector<double>> e = extractCell(cell);
		out[i] = tabfun(e[0], wgt);
	}

	return out;
}


SpatDataFrame SpatRaster::zonal_poly_weighted(SpatVector x, SpatRaster w, bool weights, bool exact, bool touches, bool narm, SpatOptions &opt) {

	SpatDataFrame out;
	std::string gtype = x.type();
	if (gtype != "polygons") {
		out.setError("SpatVector must have polygon geometry");
		return out;
	}
	
	if (!compare_geom(w, false, true, opt.get_tolerance(), true)) {
		out.setError(getError());
		return(out);
	}
	if (!hasValues()) {
		out.setError("raster has no values");
		return out;
	}
	if (!w.hasValues()) {
		out.setError("raster has no values");
		return out;
	}

    unsigned nl = nlyr();
    unsigned ng = x.size();

	std::vector<std::vector<double>> zv(nl, std::vector<double>(ng));
	
    SpatRaster r = geometry(1);
    for (size_t i=0; i<ng; i++) {
		SpatGeom g = x.getGeom(i);
		SpatVector p(g);
		p.srs = x.srs;
		std::vector<double> cell, wgt;
		if (weights) {
			rasterizeCellsWeights(cell, wgt, p, opt);
		} else if (exact) {
			rasterizeCellsExact(cell, wgt, p, opt);
		} else {
			cell = rasterizeCells(p, touches, opt);
        }
		
		std::vector<std::vector<double>> e = extractCell(cell);
		std::vector<std::vector<double>> we = w.extractCell(cell);
		
 		if (weights || exact) {
			if (narm) {
				for (size_t j=0; j<nl; j++) {
					double wsum = 0;
					double vsum = 0;
					for (size_t k=0; k<e[j].size(); k++) {
						if (!std::isnan(e[j][k])) {
							wsum += we[0][k] * wgt[k];
							vsum += (e[j][k] * we[0][k] * wgt[k]);  
						}
					}
					zv[j][i] = vsum / wsum;
				}
				for (size_t j=0; j<nl; j++) {
					double wsum = 0;
					double vsum = 0;
					for (size_t k=0; k<e[j].size(); k++) {
						if ((!std::isnan(e[j][k])) && (!std::isnan(we[0][k]))) {
							wsum += we[0][k] * wgt[k];
							vsum += (e[j][k] * we[0][k] * wgt[k]);  
						}
					}
					zv[j][i] = vsum / wsum;
				}				


			} else {
				for (size_t j=0; j<nl; j++) {
					double wsum = 0;
					double vsum = 0;
					for (size_t k=0; k<e[j].size(); k++) {
						wsum += we[0][k];
						vsum += (e[j][k] * we[0][k]);  
					}
					zv[j][i] = vsum / wsum;
				}
			}
		} else {
			if (narm) {
				for (size_t j=0; j<nl; j++) {
					double wsum = 0;
					double vsum = 0;
					for (size_t k=0; k<e[j].size(); k++) {
						if ((!std::isnan(e[j][k])) && (!std::isnan(we[0][k]))) {
							wsum += we[0][k];
							vsum += (e[j][k] * we[0][k]);  
						}
					}
					zv[j][i] = vsum / wsum;
				}				
			} else {
				for (size_t j=0; j<nl; j++) {
					double wsum = 0;
					double vsum = 0;
					for (size_t k=0; k<e[j].size(); k++) {
						wsum += we[0][k];
						vsum += (e[j][k] * we[0][k]);  
					}
					zv[j][i] = vsum / wsum;
				}
			}
		}
	}
	std::vector<std::string> nms = getNames();	
	for (size_t j=0; j<nl; j++) {
		out.add_column(zv[j], nms[j]);
	}
	
	return out;
}

