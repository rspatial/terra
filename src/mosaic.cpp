/*
#include "spatRasterMultiple.h"

std::vector<double> mean2d(const std::vector<std::vector<double>> &x) {
	size_t n = x[0].size();
	size_t nn = x.size();
	std::vector<double> out(n, NAN);
	size_t d;
	double v;
	for (size_t i=0; i<n; i++) {
		v = 0;
		d = 0;
		for (size_t j=0; j<nn; j++) {
			if (!std::isnan(x[i][j])) {
				v += x[i][j];
				d++;
			}
		}
		if (d > 0) {
			out[i] = v / d;		
		}
	}
	return out;
}	


std::vector<double> SpatRaster::readExtent(SpatExtent e) {
	std::vector<double> out;
	return out;
}


SpatRaster SpatRasterCollection::summary(std::string fun, SpatOptions &opt) {

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
	
	bool anyvals = x[0].hasValues();
	out = x[0].geometry();
	SpatExtent e = x[0].getExtent();
	for (size_t i=n; i>0; i--) {
		// for now, must have same nlyr; but should be easy to recycle.
								 //  lyrs, crs, warncrs, ext, rowcol, res
		if (!x[0].compare_geom(x[i], true, false, false, false, false, true)) {
			out.setError(x[0].msg.error);
			return(out);
		}
		e.unite(x[i].getExtent());
		if (x[i].hasValues()) {
			anyvals = true;
			x[i].readStart();
		} else {
			x.erase(x.begin()+i, x.begin()+i+1);
			n--;
		}
	}
	out.setExtent(e, true);
	e = out.getExtent();
	if (!anyvals) return out;

 	if (!out.writeStart(opt)) { return out; }
//	out.fill(NAN);
    double xmin = e.xmin;
    double xmax = e.xmax;
	
	std::vector<std::vector<double>> v(n);
	for (size_t i=0; i < out.bs.n; i++) {
		double ymin = out.yFromRow(out.bs.row[i]);
		double ymax = out.yFromRow(out.bs.row[i]+out.bs.nrows[i]-1);
		for (size_t j=0; j<n; j++) {
			SpatExtent readext = {xmin, xmax, ymin, ymax};
			v[j] = x[i].readExtent(readext);
		}
		// apply method
	}


	out.writeStop();
	//readStop();
	return(out);
}

*/

/// to be done
/*
SpatRaster SpatRaster::mosaic(SpatRaster x, SpatOptions &opt) {

	SpatRaster out;
	if (!(hasValues() | x.hasValues())) {
		SpatRasterCollection z = { *self, x };
		out = merge(x, opt)
		return out;
	}
	
	out = x[0].geometry();
	SpatExtent e = x[0].getExtent();
	for (size_t i=1; i<n; i++) {
		// for now, must have same nlyr; but should be easy to recycle.
								 //  lyrs, crs, warncrs, ext, rowcol, res
		if (!x[0].compare_geom(x[i], true, false, false, false, false, true)) {
			out.setError(x[0].msg.error);
			return(out);
		}
		e.unite(x[i].getExtent());
		if (x[i].hasValues()) anyvals = true;
	}
	out.setExtent(e, true);
	if (!anyvals) return out;

 //   out.setResolution(xres(), yres());
 	if (!out.writeStart(opt)) { return out; }
	out.fill(NAN);

// to start only consider the case of two rasters. 
    
	e = x[0].getExtent();
	SpatExtent e2 = x[1].getExtent();
	e.intersect(e2);
	if ( e.valid() ) { // there is overlap
		SpatOptions xopt(opt);
		SpatRaster r1 = x[0].crop(e, "near", xopt);
		SpatRaster r2 = x[1].crop(e, "near", xopt);
		// simple, no gradient
		r1.arith(r2, "+", xopt);
	} 
	
	


	for (size_t i=0; i<n; i++) {
		SpatRaster r = x[i];
		if (!r.hasValues()) continue;
		BlockSize bs = r.getBlockSize(opt);
		if (!r.readStart()) {
			out.setError(r.getError());
			return(out);
		}

		for (size_t j=0; j<bs.n; j++) {
            std::vector<double> v = r.readValues(bs.row[j], bs.nrows[j], 0, r.ncol());
            unsigned row1 = out.rowFromY(r.yFromRow(bs.row[j]));
            unsigned row2 = out.rowFromY(r.yFromRow(bs.row[j]+bs.nrows[j]-1));
            unsigned col1 = out.colFromX(r.xFromCol(0));
            unsigned col2 = out.colFromX(r.xFromCol(r.ncol()-1));
            if (!out.writeValues(v, row1, row2-row1+1, col1, col2-col1+1)) return out;
		}
		r.readStop();
	}

	out.writeStop();
	return(out);
}



SpatRaster SpatRasterCollection::mosaic(SpatOptions &opt) {

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

	bool anyvals = false;
	if (x[0].hasValues()) anyvals = true;
	out = x[0].geometry();
	SpatExtent e = x[0].getExtent();
	for (size_t i=1; i<n; i++) {
		// for now, must have same nlyr; but should be easy to recycle.
								 //  lyrs, crs, warncrs, ext, rowcol, res
		if (!x[0].compare_geom(x[i], true, false, false, false, false, true)) {
			out.setError(x[0].msg.error);
			return(out);
		}
		e.unite(x[i].getExtent());
		if (x[i].hasValues()) anyvals = true;
	}
	out.setExtent(e, true);
	if (!anyvals) return out;

 //   out.setResolution(xres(), yres());
 	if (!out.writeStart(opt)) { return out; }
	out.fill(NAN);

// to start only consider the case of two rasters. 
    
	e = x[0].getExtent();
	SpatExtent e2 = x[1].getExtent();
	e.intersect(e2);
	if ( e.valid() ) { // there is overlap
		SpatOptions xopt(opt);
		SpatRaster r1 = x[0].crop(e, "near", xopt);
		SpatRaster r2 = x[1].crop(e, "near", xopt);
		// simple, no gradient
		r1.arith(r2, "+", xopt);
	} 
	
	


	for (size_t i=0; i<n; i++) {
		SpatRaster r = x[i];
		if (!r.hasValues()) continue;
		BlockSize bs = r.getBlockSize(opt);
		if (!r.readStart()) {
			out.setError(r.getError());
			return(out);
		}

		for (size_t j=0; j<bs.n; j++) {
            std::vector<double> v = r.readValues(bs.row[j], bs.nrows[j], 0, r.ncol());
            unsigned row1 = out.rowFromY(r.yFromRow(bs.row[j]));
            unsigned row2 = out.rowFromY(r.yFromRow(bs.row[j]+bs.nrows[j]-1));
            unsigned col1 = out.colFromX(r.xFromCol(0));
            unsigned col2 = out.colFromX(r.xFromCol(r.ncol()-1));
            if (!out.writeValues(v, row1, row2-row1+1, col1, col2-col1+1)) return out;
		}
		r.readStop();
	}

	out.writeStop();
	return(out);
}

*/