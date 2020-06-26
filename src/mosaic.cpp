
#include "spatRasterMultiple.h"

/// to be done


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
		BlockSize bs = r.getBlockSize(4, opt.get_memfrac());
		r.readStart();
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

