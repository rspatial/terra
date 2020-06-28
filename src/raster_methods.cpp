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

//#include <vector>
#include "spatRasterMultiple.h"
#include "recycle.h"
#include "vecmathfun.h"
//#include "vecmath.h"
#include <cmath>
#include "math_utils.h"


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
	out.bs = getBlockSize(opt.get_blocksizemp(), opt.get_memfrac(), opt.get_steps());
    #ifdef useRcpp
	out.pbar = new Progress(out.bs.n+2, opt.do_progress(bs.n));
	out.pbar->increment();
	#endif
	
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

	std::function<double(std::vector<double>&, bool)> theFun;
	theFun = getFun(fun);

	for (size_t i=0; i<out.bs.n; i++) {
        std::vector<double> a = readBlock(out.bs, i);
		unsigned nc = out.bs.nrows[i] * ncol();
		std::vector<double> b(nc * nl);
		for (size_t j=0; j<nc; j++) {
			for (size_t k=0; k<ird.size(); k++) {
				v[ird[k]][jrd[k]] = a[j+k*nc];
			}
			for (size_t k=0; k<ui.size(); k++) {
				size_t off = k * nc + j;
				b[off] = theFun(v[k], narm);
			}
		}
		if (!out.writeValues(b, out.bs.row[i], out.bs.nrows[i], 0, ncol())) return out;
	}
	readStop();
	out.writeStop();
	return(out);
}



SpatRaster SpatRaster::mask(SpatRaster x, bool inverse, double maskvalue, double updatevalue, SpatOptions &opt) {

	unsigned nl = std::max(nlyr(), x.nlyr());
	SpatRaster out = geometry(nl);

	if (!out.compare_geom(x, false, true, false, true, true, false)) {
		return(out);
	}

	readStart();
	x.readStart();
  	if (!out.writeStart(opt)) { return out; }
	std::vector<double> v, m;
	for (size_t i = 0; i < out.bs.n; i++) {
		v = readValues(out.bs.row[i], out.bs.nrows[i], 0, ncol());
		m = x.readValues(out.bs.row[i], out.bs.nrows[i], 0, ncol());
		recycle(v, m);
		if (inverse) {
			if (std::isnan(maskvalue)) {
				for (size_t i=0; i < v.size(); i++) {
					if (!std::isnan(m[i])) {
						v[i] = updatevalue;
					}
				}
			} else {
				for (size_t i=0; i < v.size(); i++) {
					if (m[i] != maskvalue) {
						v[i] = updatevalue;
					}
				}
			}
		} else {
			if (std::isnan(maskvalue)) {
				for (size_t i=0; i < v.size(); i++) {
					if (std::isnan(m[i])) {
						v[i] = updatevalue;
					}
				}
			} else {
				for (size_t i=0; i < v.size(); i++) {
					if (m[i] == maskvalue) {
						v[i] = updatevalue;
					}
				}
			}
		}
		if (!out.writeValues(v, out.bs.row[i], out.bs.nrows[i], 0, ncol())) return out;

	}
	out.writeStop();
	readStop();
	x.readStop();
	return(out);
}


SpatRaster SpatRaster::mask(SpatVector x, bool inverse, double updatevalue, SpatOptions &opt) {

//return grasterize(x, "", {updatevalue}, NAN, true, false, !inverse, opt);
// gdal_rasterize with inverse=true does not work well with overlapping polygons 
// also can't use NA as update value, it appears
// looks like GDAL bug
//	eturn grasterize(x, "", {updatevalue}, NAN, true, false, !inverse, opt);
// so do it in two steps
	SpatRaster out;
	if (!hasValues()) {
		out.setError("SpatRaster has no values");
		return out;
	}
	std::string filename = opt.get_filename();
	opt.set_filename("");
	SpatRaster m = rasterize(x, "", {1}, 0, false, false, false, opt);
	opt.set_filename(filename);
	out = mask(m, inverse, 0, updatevalue, opt);
	return(out);

/*
	std::string filename = opt.get_filename();
	opt.set_filename("");
	std::vector<double> feats(x.size(), 1) ;
	SpatRaster m = rasterize(x, feats, 0, false, opt);
	opt.set_filename(filename);
	SpatRaster out = mask(m, inverse, 0, updatevalue, opt);
	return(out);
*/
}



SpatRaster SpatRaster::transpose(SpatOptions &opt) {

	SpatRaster out = geometry();
	SpatExtent eold = getExtent();
	SpatExtent enew = getExtent();
	enew.xmin = eold.ymin;
	enew.xmax = eold.ymax;
	enew.ymin = eold.xmin;
	enew.ymax = eold.xmax;
	out.setExtent(enew, false, "");
	out.source[0].ncol = nrow();
	out.source[0].nrow = ncol();
	if (!hasValues()) return out;
 	if (!out.writeStart(opt)) { return out; }
	readStart();
	for (size_t i=0; i < out.bs.n; i++) {
		unsigned nr = nrow();
		unsigned nc = out.bs.nrows[i];
		std::vector<double> v = readValues(0, nr, out.bs.row[i], nc);
		std::vector<double> vv(v.size());
		for (size_t lyr=0; lyr<nlyr(); lyr++) {
			size_t off = lyr*ncell();
			for (size_t r = 0; r < nr; r++) {
				size_t rnc = off + r * nc;
				for (size_t c = 0; c < nc; c++) {
					vv[c*nr+r+off] = v[rnc+c];
				}
			}
		}
		if (!out.writeValues(vv, out.bs.row[i], out.bs.nrows[i], 0, out.ncol())) return out;
	}
	out.writeStop();
	readStop();
	return(out);
}




SpatRaster SpatRaster::trim(unsigned padding, SpatOptions &opt) {

	long nrl = nrow() * nlyr();
	long ncl = ncol() * nlyr();

	std::vector<double> v;
	unsigned r;
	for (r=0; r<nrow(); r++) {
		v = readValues(r, 1, 0, ncol());
		if (std::count_if( v.begin(), v.end(), [](double d) { return std::isnan(d); } ) < ncl) {
			break;
		}
	}

	if (r == nrow()) { //stop('only NA values found')
	}
	unsigned firstrow = std::min(std::max(r - padding, unsigned(0)), nrow());

	for (r=nrow()-1; r>firstrow; r--) {
		v = readValues(r, 1, 0, ncol());
		if (std::count_if(v.begin(), v.end(), [](double d) { return std::isnan(d); } ) < ncl) {
			break;
		}
	}

	unsigned lastrow = std::max(std::min(r+padding, nrow()), unsigned(0));

	unsigned tmp;
	if (lastrow < firstrow) {
		tmp = firstrow;
		firstrow = lastrow;
		lastrow = tmp;
	}
	unsigned c;
	for (c=0; c<ncol(); c++) {
		v = readValues(0, nrow(), c, 1);
		if (std::count_if( v.begin(), v.end(), [](double d) { return std::isnan(d); } ) < nrl) {
			break;
		}
	}
	unsigned firstcol = std::min(std::max(c-padding, unsigned(0)), ncol());


	for (c=ncol()-1; c>firstcol; c--) {
		v = readValues(0, nrow(), c, 1);
		if (std::count_if( v.begin(), v.end(), [](double d) { return std::isnan(d); } ) < nrl) {
			break;
		}
	}
	unsigned lastcol = std::max(std::min(c+padding, ncol()), unsigned(0));

	if (lastcol < firstcol) {
		tmp = firstcol;
		firstcol = lastcol;
		lastcol = tmp;
	}

	std::vector<double> res = resolution();
	double xr = res[0];
	double yr = res[1];
	SpatExtent e = SpatExtent(xFromCol(firstcol)-0.5*xr, xFromCol(lastcol)+0.5*xr, yFromRow(lastrow)-0.5*yr, yFromRow(firstrow)+0.5*yr);

	return( crop(e, "near", opt) ) ;
}





void clamp_vector(std::vector<double> &v, double low, double high, bool usevalue) {
	size_t n = v.size();
	if (usevalue) {
		for (size_t i=0; i<n; i++) {
			if ( v[i] < low ) {
				v[i] = low;
			} else if ( v[i] > high ) {
				v[i] = high;
			}
		}
	} else {
		for (size_t i=0; i<n; i++) {
			if ( (v[i] < low )| (v[i] > high)) {
				v[i] = NAN;
			}
		}
	}
}



SpatRaster SpatRaster::clamp(double low, double high, bool usevalue, SpatOptions &opt) {

	SpatRaster out = geometry(nlyr());
	if (low > high) {
		out.setError("lower clamp value cannot be larger than the higher clamp value");
		return out;
	}
	if (!hasValues()) {
		out.setError("cannot clamp a raster with no values");
		return out;
	}

  	if (!out.writeStart(opt)) { return out; }
	readStart();
	for (size_t i = 0; i < out.bs.n; i++) {
		std::vector<double> v = readBlock(out.bs, i);
		clamp_vector(v, low, high, usevalue);
		if (!out.writeValues(v, out.bs.row[i], out.bs.nrows[i], 0, ncol())) return out;
	}
	readStop();
	out.writeStop();
	return(out);
}




SpatRaster SpatRaster::selRange(SpatRaster x, int z, SpatOptions &opt) {

	int nl = nlyr();
	z = std::max(1, std::min(z, nl));
	SpatRaster out = geometry(z);
	if (!out.compare_geom(x, false, false)) {
		return(out);
	}
	if (!hasValues()) return(out);

	if (x.nlyr() > 1) {
		out.setError("index raster must have only one layer");
		return out;
	}
	if (!x.hasValues()) {
		out.setError("index raster has no values");
		return out;
	}

 	if (!out.writeStart(opt)) { return out; }
	readStart();
	x.readStart();
	for (size_t i=0; i<out.bs.n; i++) {
		std::vector<double> v = readBlock(out.bs, i);
		std::vector<double> idx = x.readBlock(out.bs, i);
		size_t is = idx.size();
		std::vector<double> vv(is*z, NAN);
		size_t ncell = out.bs.nrows[i] * ncol();

		for (size_t j=0; j<is; j++) {
			int start = idx[j] - 1;
			if ((start >= 0) && (start < nl)) {
				int zz = std::min(nl-start, z);
				for (int i = 0; i<zz; i++){
					size_t offin = (start + i) * ncell + j;
					size_t offout = i * ncell + j;
					vv[offout] = v[offin];   
				}
			}
		}
	
		//for (size_t j=0; j<is; j++) {
		//	int index = idx[j] - 1;
		//	if ((index >= 0) && (index < nl)) {
		//		vv[j] = v[j + index * ncell];
		//	}
		//}
		if (!out.writeValues(vv, out.bs.row[i], out.bs.nrows[i], 0, ncol())) return out;
	}
	readStop();
	x.readStop();
	out.writeStop();
	return(out);
}



SpatRaster SpatRaster::rapply(SpatRaster x, std::string fun, bool narm, SpatOptions &opt) {

	SpatRaster out = geometry(1);
	if (!out.compare_geom(x, false, false)) {
		return(out);
	}
	if (!hasValues()) return(out);
	if (!x.hasValues()) {
		out.setError("index raster has no values");
		return out;
	}
	if (x.nlyr() != 2) {
		out.setError("index raster must have two layers");
		return out;
	}
	
	std::vector<std::string> f {"sum", "mean", "min", "max", "prod", "any", "all"};
	if (std::find(f.begin(), f.end(), fun) == f.end()) {
		out.setError("unknown apply function");
		return out;
	}
	std::function<double(std::vector<double>&, bool)> theFun;
	theFun = getFun(fun);

	int nl = nlyr();
 	if (!out.writeStart(opt)) { return out; }
	readStart();
	x.readStart();
	for (size_t i=0; i<out.bs.n; i++) {
		std::vector<double> v = readBlock(out.bs, i);
		std::vector<double> idx = x.readBlock(out.bs, i);
		size_t ncell = out.bs.nrows[i] * ncol();
		std::vector<double> vv(ncell, NAN);
		for (size_t j=0; j<ncell; j++) {
			int start = idx[j] - 1;
			int end   = idx[j+ncell];
			if ((start >= 0) && (end <= nl) && (end >= start)) {
				std::vector<double> se;
				se.reserve(end-start+1);
				for (int i = start; i<end; i++){
					size_t off = i * ncell + j;
					se.push_back(v[off]);   
				}
				vv[j] = theFun(se, narm);
			}
		}
		if (!out.writeValues(vv, out.bs.row[i], out.bs.nrows[i], 0, ncol())) return out;
	}
	readStop();
	x.readStop();
	out.writeStop();
	return(out);
}


std::vector<std::vector<double>> SpatRaster::rappvals(SpatRaster x, size_t startrow, size_t nrows) {

	std::vector<std::vector<double>> r;
	if (!compare_geom(x, false, false)) {
		return(r);
	}
	if (!hasValues()) {
		return r;
	}
	if (!x.hasValues()) {
		return r;
	}
	if (x.nlyr() != 2) {
		return r;
	}
	
	int nl = nlyr();
	readStart();
	x.readStart();
	std::vector<double> v = readValues(startrow, nrows, 0, ncol());
	std::vector<double> idx = x.readValues(startrow, nrows, 0, ncol());
	size_t ncell = nrows * ncol();
	r.resize(ncell);
	for (size_t j=0; j<ncell; j++) {
		int start = idx[j] - 1;
		int end   = idx[j+ncell];
		if ((start >= 0) && (end <= nl) && (end >= start)) {
			r[j].reserve(end-start+1);
			for (int i = start; i<end; i++){
				size_t off = i * ncell + j;
				r[j].push_back(v[off]);   
			}
		}
	}
	readStop();
	x.readStop();
	return(r);
}


bool disaggregate_dims(std::vector<unsigned> &fact, std::string &message ) {

	unsigned fs = fact.size();
	if ((fs > 3) | (fs == 0)) {
		message = "argument 'fact' should have length 1, 2, or 3";
		return false;
	}
	auto min_value = *std::min_element(fact.begin(),fact.end());
	if (min_value < 1) {
		message = "values in argument 'fact' should be > 0";
		return false;
	}
	auto max_value = *std::max_element(fact.begin(),fact.end());
	if (max_value == 1) {
		message = "all values in argument 'fact' are 1, nothing to do";
		return false;
	}

	fact.resize(3);
	if (fs == 1) {
		fact[1] = fact[0];
	}
	fact[2] = 1;
	return true;
}



SpatRaster SpatRaster::disaggregate(std::vector<unsigned> fact, SpatOptions &opt) {

    SpatRaster out = geometry();


	std::string message = "";
	bool success = disaggregate_dims(fact, message);
	if (!success) {
		out.setError(message);
		return out;
	}

    out.source[0].nrow = out.source[0].nrow * fact[0];
    out.source[0].ncol = out.source[0].ncol * fact[1];
    out.source[0].nlyr = out.source[0].nlyr * fact[2];

    if (!hasValues()) {
        return out;
    }

	unsigned bsmp = opt.get_blocksizemp()*fact[0]*fact[1]*fact[2];
	BlockSize bs = getBlockSize(bsmp, opt.get_memfrac());
	//opt.set_blocksizemp();
	std::vector<double> v, vout;
	unsigned nc = ncol();
	unsigned nl = nlyr();
	std::vector<double> newrow(nc*fact[1]);
  	readStart();

  	if (!out.writeStart(opt)) { return out; }
	for (size_t i = 0; i < bs.n; i++) {
		v = readValues(bs.row[i], bs.nrows[i], 0, nc);
		vout.resize(0);
		vout.reserve(v.size() * fact[0] * fact[1] * fact[2]);

		for (size_t lyr=0; lyr<nl; lyr++) {
			for (size_t row=0; row<bs.nrows[i]; row++) {
				unsigned rowoff = row*nc + lyr*nc*bs.nrows[i];
				// for each new column
				unsigned jfact = 0;
				for (size_t j=0; j<nc; j++) {
					unsigned coloff = rowoff + j;
					for (size_t k=0; k<fact[1]; k++) {
						newrow[jfact+k] = v[coloff];
					}
					jfact += fact[1];
				}
				// for each new row
				for (size_t j=0; j<fact[0]; j++) {
					vout.insert(vout.end(), newrow.begin(), newrow.end());
				}
			}
		}
		if (!out.writeValues(vout, bs.row[i]*fact[0], bs.nrows[i]*fact[0], 0, out.ncol())) return out;
	}
	vout.resize(0);
	out.writeStop();
	readStop();
	return(out);
}




SpatRaster SpatRaster::init(std::string value, bool plusone, SpatOptions &opt) {

	SpatRaster out = geometry();

	std::vector<std::string> f {"row", "col", "cell", "x", "y", "chess"};
	bool test = std::find(f.begin(), f.end(), value) == f.end();
	if (test) {
		out.setError("not a valid init option");
		return out;
	}

	size_t nr = nrow();
	size_t steps = nr; // for the pbar
	if (value == "chess") {
		steps = steps / 2;
	}
	opt.set_steps(steps);
 	if (!out.writeStart(opt)) { return out; }

	if (value == "row") {
		std::vector<double> v(ncol());
		for (size_t i = 0; i < nr; i++) {
			std::fill(v.begin(), v.end(), i+plusone);
			if (!out.writeValues(v, i, 1, 0, ncol())) return out;
		}
	} else if (value == "col") {
		std::vector<double> col(ncol());
		double start = plusone ? 1 : 0;
		std::iota(col.begin(), col.end(), start);
		for (unsigned i = 0; i < nr; i++) {
			if (!out.writeValues(col, i, 1, 0, ncol())) return out;
		}
	} else if (value == "cell") {
		std::vector<long> col(ncol());
		std::iota(col.begin(), col.end(), 0);
		std::vector<long> row(1);
		for (unsigned i = 0; i < nr; i++) {
			row[0] = i;
			std::vector<double> v = cellFromRowCol(row, col);
			if (plusone) for(double& d : v) d=d+1;
			if (!out.writeValues(v, i, 1, 0, ncol())) return out;
		}
	} else if (value == "x") {
		std::vector<long> col(ncol());
		std::iota(col.begin(), col.end(), 0);
		std::vector<double> x = xFromCol(col);
		for (unsigned i = 0; i < nr; i++) {
			if (!out.writeValues(x, i, 1, 0, ncol())) return out;
		}
	} else if (value == "y") {
		std::vector<double> v(ncol());
		for (unsigned i = 0; i < nr; i++) {
			double y = yFromRow(i);
			std::fill(v.begin(), v.end(), y);
			if (!out.writeValues(v, i, 1, 0, ncol())) return out;
		}
	} else if (value == "chess") {
		std::vector<double> a(ncol());
		std::vector<double> b(ncol());
		size_t nr = nrow();
		size_t nc = ncol();
		a[0] = 1;
		b[0] = 0;
		for (size_t i=1; i<nc; i++) {
			bool test = i%2 == 0;
			a[i] = test;
			b[i] = !test;
		}
		out.bs.n = nr/2; // for the pbar
		for (unsigned i=0; i<(nr-1); i=i+2) {
			if (!out.writeValues(a, i, 1, 0, ncol())) return out;
			if (!out.writeValues(b, i+1, 1, 0, ncol())) return out;
		}
		if (nr%2 == 0) {
			if (!out.writeValues(a, nr-2, 1, 0, ncol())) return out;
			if (!out.writeValues(b, nr-1, 1, 0, ncol())) return out;
		} else {
			if (!out.writeValues(a, nr-1, 1, 0, ncol())) return out;
		}
	}

	out.writeStop();
	return(out);
}



SpatRaster SpatRaster::init(double value, SpatOptions &opt) {
	SpatRaster out = geometry();
 	if (!out.writeStart(opt)) { return out; }
	unsigned nc = ncol();
	std::vector<double> v(out.bs.nrows[0]*nc, value);
	for (size_t i = 0; i < out.bs.n; i++) {
		if (i > 0 && i == (out.bs.n-1)) {
			v.resize(bs.nrows[i]*nc);
		}
		if (!out.writeValues(v, i, 1, 0, ncol())) return out;
	}
	out.writeStop();
	return(out);
}


SpatRaster SpatRaster::isnan(SpatOptions &opt) {
	SpatRaster out = geometry();
    if (!hasValues()) return out;

	if (!out.writeStart(opt)) { return out; }
	readStart();
	for (size_t i=0; i<out.bs.n; i++) {
		std::vector<double> v = readBlock(out.bs, i);
		for (double &d : v) d = std::isnan(d);
		if (!out.writeValues(v, out.bs.row[i], out.bs.nrows[i], 0, ncol())) return out;
	}
	readStop();
	out.writeStop();
	return(out);
}


SpatRaster SpatRaster::isnotnan(SpatOptions &opt) {
	SpatRaster out = geometry();
    if (!hasValues()) return out;

	if (!out.writeStart(opt)) { return out; }
	readStart();
	for (size_t i=0; i<out.bs.n; i++) {
		std::vector<double> v = readBlock(out.bs, i);
		for (double &d : v) d = ! std::isnan(d);
		if (!out.writeValues(v, out.bs.row[i], out.bs.nrows[i], 0, ncol())) return out;
	}
	readStop();
	out.writeStop();
	return(out);
}


SpatRaster SpatRaster::isfinite(SpatOptions &opt) {
	SpatRaster out = geometry();
    if (!hasValues()) return out;

	if (!out.writeStart(opt)) { return out; }
	readStart();
	for (size_t i=0; i<out.bs.n; i++) {
		std::vector<double> v = readBlock(out.bs, i);
		for (double &d : v) d = std::isfinite(d);
		if (!out.writeValues(v, out.bs.row[i], out.bs.nrows[i], 0, ncol())) return out;
	}
	readStop();
	out.writeStop();
	return(out);
}


SpatRaster SpatRaster::isinfinite(SpatOptions &opt) {
	SpatRaster out = geometry();
    if (!hasValues()) return out;

	if (!out.writeStart(opt)) { return out; }
	readStart();
	for (size_t i=0; i<out.bs.n; i++) {
		std::vector<double> v = readBlock(out.bs, i);
		for (double &d : v) d = std::isinf(d);
		if (!out.writeValues(v, out.bs.row[i], out.bs.nrows[i], 0, ncol())) return out;
	}
	readStop();
	out.writeStop();
	return(out);
}


SpatRaster SpatRaster::rotate(bool left, SpatOptions &opt) {

	unsigned nc = ncol();
	unsigned nl = nlyr();
	unsigned hnc = (nc / 2);
	double addx = hnc * xres();
	if (left) {
		addx = -addx;
	}
	SpatRaster out = geometry();
	out.extent.xmin = out.extent.xmin + addx;
	out.extent.xmax = out.extent.xmax + addx;

	if (!hasValues()) return out;

 	if (!out.writeStart(opt)) { return out; }
	readStart();
	std::vector<double> b;
	for (size_t i=0; i < out.bs.n; i++) {
		std::vector<double> a = readBlock(out.bs, i);
		for (size_t j=0; j < nl; j++) {
			for (size_t r=0; r < out.bs.nrows[i]; r++) {
				unsigned s1 = j * out.bs.nrows[i] * nc + r * nc;
				unsigned e1 = s1 + hnc;
				unsigned s2 = e1;
				unsigned e2 = s1 + nc;
				b.insert(b.end(), a.begin()+s2, a.begin()+e2);
				b.insert(b.end(), a.begin()+s1, a.begin()+e1);
			}
		}
		if (!out.writeValues(b, out.bs.row[i], nrow(), 0, ncol())) return out;
		b.resize(0);
	}
	out.writeStop();
	readStop();
	return(out);
}



SpatRaster SpatRaster::extend(SpatExtent e, SpatOptions &opt) {

	SpatRaster out = geometry(nlyr());
	e = out.align(e, "near");
	e.unite(extent);
	if (extent.compare(e, "==", 1000)) {
		out = deepCopy();
		return out;
	}

	out.setExtent(e, true);
	if (!hasValues()) return(out);

 	if (!out.writeStart(opt)) { return out; }
	out.fill(NAN);
	BlockSize bs = getBlockSize(4, opt.get_memfrac());
	readStart();
	for (size_t i=0; i<bs.n; i++) {
        std::vector<double> v = readValues(bs.row[i], bs.nrows[i], 0, ncol());
        unsigned row1 = out.rowFromY(yFromRow(bs.row[i]));
        unsigned row2 = out.rowFromY(yFromRow(bs.row[i]+bs.nrows[i]-1));
        unsigned col1 = out.colFromX(xFromCol(0));
        unsigned col2 = out.colFromX(xFromCol(ncol()-1));
        if (!out.writeValues(v, row1, row2-row1+1, col1, col2-col1+1)) return out;
	}
	readStop();
	out.writeStop();
	return(out);
}


/*
SpatRaster SpatRaster::filler(SpatRaster x, SpatOptions &opt) {

	unsigned nl = std::max(nlyr(), x.nlyr());
	SpatRaster out = geometry(nl);
	SpatExtent e = getExtent();
	SpatExtent xe = x.getExtent();
	
	e.intersect(xe);
	if (!e.valid()) {
		out.setError("SpatRasters do not overlap");
		return out;
	}

	if (!shared_basegeom(x, 0.1, true)) {
		out.setError("raster dimensions do not match");
		return(out);
	}
	// x is larger
	if (e.compare(xe, "<=")) {
		return x.crop(e, "near", opt);
	}
	// x is not within
	if (!e.compare(xe, ">=")) {
		SpatOptions xopt(opt);
		x = x.crop(e, "near", xopt);
		xe = x.getExtent();
	}

	std::vector<unsigned> rc = rowcolFromExtent(e);


	if (!x.hasValues()) {
		return *this;
	}
	if (!hasValues()) {
		if (rmatch) {
			return x.deepCopy();
		} else {
			SpatExtent e = getExtent();
			return x.extend(e, opt);
		}
	}
	
	readStart();
	x.readStart();
  	if (!out.writeStart(opt)) { return out; }
	std::vector<double> v, m;
	for (size_t i = 0; i < out.bs.n; i++) {
		v = readValues(out.bs.row[i], out.bs.nrows[i], 0, ncol());
		m = x.readValues(out.bs.row[i], out.bs.nrows[i], 0, ncol());
		recycle(v, m);
		if (std::isnan(value)) {
			for (size_t i=0; i < v.size(); i++) {
				if (std::isnan(v[i])) {
					v[i] = m[i];
				}
			}
		} else {
			for (size_t i=0; i < v.size(); i++) {
				if (v[i] == value) {
					v[i] = m[i];
				}
			}
		}
		if (!out.writeValues(v, out.bs.row[i], out.bs.nrows[i], 0, ncol())) return out;
	}
	out.writeStop();
	readStop();
	x.readStop();
	return(out);
}
*/

bool SpatRaster::shared_basegeom(SpatRaster &x, double tol, bool test_overlap) {
	if (!compare_origin(x.origin(), tol)) return false;	
	if (!about_equal(xres(), x.xres(), xres() * tol)) return false; 
	if (!about_equal(yres(), x.yres(), yres() * tol)) return false; 
	if (test_overlap) {
		SpatExtent e = x.getExtent();
		e.intersect(getExtent());
		if (!e.valid()) return false;
	}
	return true;
}




SpatRaster SpatRaster::cover(SpatRaster x, double value, SpatOptions &opt) {

	unsigned nl = std::max(nlyr(), x.nlyr());
	SpatRaster out = geometry(nl);
	
	bool rmatch = false;
					 //  lyrs, crs, warncrs, ext, rowcol, res
	if (out.compare_geom(x, false, false, true)) {
		rmatch = true;
	} else {
	//	if (!shared_basegeom(x, 0.1, true)) {
		out.setError("raster dimensions do not match");
		return(out);
	//	} else {
	//		out.msg.has_error = false;
	//		out.msg.error = "";
	//		SpatExtent e = getExtent();
	//		SpatExtent xe = x.getExtent();
	//		double prec = std::min(xres(), yres())/1000;
	//		if (!xe.compare(e, "<=", prec)) {
	//			SpatOptions xopt(opt);
	//			x = x.crop(e, "near", xopt);			
	//		}
	//	}
	}

	if (!x.hasValues()) {
		return *this;
	}
	if (!hasValues()) {
		if (rmatch) {
			return x.deepCopy();
		} else {
			SpatExtent e = getExtent();
			return x.extend(e, opt);
		}
	}
	
	readStart();
	x.readStart();
  	if (!out.writeStart(opt)) { return out; }
	std::vector<double> v, m;
	for (size_t i = 0; i < out.bs.n; i++) {
		v = readValues(out.bs.row[i], out.bs.nrows[i], 0, ncol());
		m = x.readValues(out.bs.row[i], out.bs.nrows[i], 0, ncol());
		recycle(v, m);
		if (std::isnan(value)) {
			for (size_t i=0; i < v.size(); i++) {
				if (std::isnan(v[i])) {
					v[i] = m[i];
				}
			}
		} else {
			for (size_t i=0; i < v.size(); i++) {
				if (v[i] == value) {
					v[i] = m[i];
				}
			}
		}
		if (!out.writeValues(v, out.bs.row[i], out.bs.nrows[i], 0, ncol())) return out;
	}
	out.writeStop();
	readStop();
	x.readStop();
	return(out);
}



SpatRaster SpatRaster::crop(SpatExtent e, std::string snap, SpatOptions &opt) {

	SpatRaster out = geometry();

	if ( !e.valid() ) {
		out.setError("invalid extent");
		return out;
	} 
	e.intersect(out.getExtent());
	if ( !e.valid() ) {
		out.setError("extents do not overlap");
		return out;
	} 

	out.setExtent(e, true, snap);
	if (!hasValues() ) {
		return(out);
	}

	double xr = xres();
	double yr = yres();

	unsigned col1 = colFromX(out.extent.xmin + 0.5 * xr);
	unsigned col2 = colFromX(out.extent.xmax - 0.5 * xr);
	unsigned row1 = rowFromY(out.extent.ymax - 0.5 * yr);
	unsigned row2 = rowFromY(out.extent.ymin + 0.5 * yr);
	if ((row1==0) && (row2==nrow()-1) && (col1==0) && (col2==ncol()-1)) {
		// same extent
		return deepCopy();
	}

	unsigned ncols = out.ncol();
 	if (!out.writeStart(opt)) { return out; }

	readStart();
	std::vector<double> v;
	for (size_t i = 0; i < out.bs.n; i++) {
		v = readValues(row1+out.bs.row[i], out.bs.nrows[i], col1, ncols);
		if (!out.writeValues(v, out.bs.row[i], out.bs.nrows[i], 0, out.ncol())) return out;
	}
	out.writeStop();
	readStop();

	return(out);
}

SpatRaster SpatRaster::flip(bool vertical, SpatOptions &opt) {

	SpatRaster out = geometry();
	if (!hasValues()) return out;
 	if (!out.writeStart(opt)) { return out; }
	readStart();
	std::vector<double> b;
	unsigned nc = ncol();
	unsigned nl = nlyr();

	if (vertical) {
		for (size_t i=0; i < out.bs.n; i++) {
			size_t ii = out.bs.n - 1 - i;
			std::vector<double> a = readBlock(out.bs, ii);
			unsigned lyrrows = nl * out.bs.nrows[ii];
			for (size_t j=0; j < lyrrows; j++) {
				unsigned start = (lyrrows - 1 - j) * nc;
				unsigned end = start + nc;
				b.insert(b.end(), a.begin()+start, a.begin()+end);
			}
			if (!out.writeValues(b, out.bs.row[i], out.bs.nrows[i], 0, ncol())) return out;
			b.resize(0);
		}
	} else {
		for (size_t i=0; i < out.bs.n; i++) {
			std::vector<double> a = readBlock(out.bs, i);
			unsigned lyrrows = nl * out.bs.nrows[i];
			for (size_t j=0; j < lyrrows; j++) {
				unsigned start = j * nc;
				unsigned end = start + nc;
				b.insert(b.end(), a.begin()+start, a.begin()+end);
			}
			if (!out.writeValues(b, out.bs.row[i], out.bs.nrows[i], 0, ncol())) return out;
			b.resize(0);
		}
	}
	out.writeStop();
	readStop();
	return(out);
}



SpatRaster SpatRaster::shift(double x, double y, SpatOptions &opt) {
	SpatRaster out = deepCopy();
	out.extent.xmin = out.extent.xmin + x;
	out.extent.xmax = out.extent.xmax + x;
	out.extent.ymin = out.extent.ymin + y;
	out.extent.ymax = out.extent.ymax + y;
	return out;

}

bool SpatRaster::compare_origin(std::vector<double> x, double tol) {
	std::vector<double> y = origin();
	if (!about_equal(x[0], y[0], xres() * tol)) return false; 
	if (!about_equal(x[1], y[1], yres() * tol)) return false; 
	return true;
}

SpatRaster SpatRasterCollection::merge(SpatOptions &opt) {

	SpatRaster out;
	unsigned n = size();

	if (n == 0) {
		out.setError("empty collection");
		return(out);
	}
	if (n == 1) {
		out = x[0].deepCopy();
		return(out);
	}

	bool any_hasvals = false;
	if (x[0].hasValues()) any_hasvals = true;
	out = x[0].geometry();
	std::vector<double> orig = x[0].origin(); 
	SpatExtent e = x[0].getExtent();
	unsigned nl = x[0].nlyr();
	for (size_t i=1; i<n; i++) {
									 //  lyrs, crs, warncrs, ext, rowcol, res
		if (!out.compare_geom(x[i], false, false, false, false, false, true)) {
			return(out);
		}
		if (!out.compare_origin(x[i].origin(), 0.1)) {
			out.setError("origin of SpatRaster " + std::to_string(i+1) + " does not match the previous SpatRaster(s)");
			return(out);			
		}
		e.unite(x[i].getExtent());
		nl = std::max(nl, x[i].nlyr());
		if (x[i].hasValues()) any_hasvals = true;
	}
	out.setExtent(e, true);
	out = out.geometry(nl);
	if (!any_hasvals) return out;

 //   out.setResolution(xres(), yres());
 	if (!out.writeStart(opt)) { return out; }
	out.fill(NAN);

	for (size_t i=0; i<n; i++) {
		SpatRaster r = x[i];
		if (!r.hasValues()) continue;
		BlockSize bs = r.getBlockSize(4, opt.get_memfrac());
		r.readStart();
		for (size_t j=0; j<bs.n; j++) {
            std::vector<double> v = r.readValues(bs.row[j], bs.nrows[j], 0, r.ncol());
            unsigned row1 = out.rowFromY(r.yFromRow(bs.row[j]));
            unsigned row2 = out.rowFromY(r.yFromRow(bs.row[j]+bs.nrows[j]-1));
            unsigned col1 = out.colFromX(r.xFromCol(0));
            unsigned col2 = out.colFromX(r.xFromCol(r.ncol()-1));
			unsigned ncols = col2-col1+1;
			unsigned nrows = row2-row1+1;
			recycle(v, ncols * nrows * nl);
			
            if (!out.writeValues(v, row1, nrows, col1, ncols)) return out;
		}
		r.readStop();
	}

	out.writeStop();
	return(out);
}




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
	BlockSize bs = getBlockSize(2, 0.5);
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

