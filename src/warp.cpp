#include <vector>
#include "spatRaster.h"
#include "crs.h"
#include "vecmath.h"


SpatRaster SpatRaster::warp(SpatRaster x, std::string method, SpatOptions &opt) {

	unsigned nl = nlyr();
	SpatRaster out = x.geometry(nl);
	out.setNames(getNames());
	std::vector<std::string> f {"bilinear", "neighbor"};
	if (std::find(f.begin(), f.end(), method) == f.end()) {
		out.setError("unknown warp method");
		return out;
	}
	if (!hasValues()) {
		return out;
	}
	
	std::string crsin = getCRS();
	std::string crsout = out.getCRS();
	bool do_prj = true;	
	if ((crsin == "") || (crsout == "")) {
		do_prj = false;	
	} else if (crsin == crsout) {
		do_prj = false;			
	}
	
	if (!do_prj) {
		SpatExtent e = out.extent;
		e.intersect(extent);
		if (!e.valid()) {
			out.setError("No spatial overlap");
			return out;
		}
	}
	
	SpatRaster xx;
	if (do_prj) {
		xx = *this;
	} else {
		unsigned xq = x.xres() / xres();
		unsigned yq = x.yres() / yres();
		if (std::max(xq, yq) > 1) {
			xq = xq == 0 ? 1 : xq;
			yq = yq == 0 ? 1 : yq;
			std::vector<unsigned> agf = {yq, xq, 1};
			SpatOptions agopt;
			if (method == "bilinear") {
				xx = aggregate(agf, "mean", true, agopt);
			} else {
				xx = aggregate(agf, "modal", true, agopt);
			}
		} else {
			xx = *this;
		}
	} 
	unsigned nc = out.ncol();

  	if (!out.writeStart(opt)) { return out; }
	for (size_t i = 0; i < out.bs.n; i++) {
        unsigned firstcell = out.cellFromRowCol(out.bs.row[i], 0);
		unsigned lastcell  = out.cellFromRowCol(out.bs.row[i]+out.bs.nrows[i]-1, nc-1);
		std::vector<double> cells(1+lastcell-firstcell);
		std::iota (std::begin(cells), std::end(cells), firstcell);
        std::vector<std::vector<double>> xy = out.xyFromCell(cells);
		if (do_prj) {
			SpatMessages msg = transform_coordinates(xy[0], xy[1], 
			crsout, crsin);
		}
		std::vector<std::vector<double>> v = xx.extractXY(xy[0], xy[1], method);
		if (!out.writeValues2(v, out.bs.row[i])) return out;
	}
	out.writeStop();

	return(out);
}



SpatRaster SpatRaster::project(std::string crs, std::string method, SpatOptions &opt) {

	SpatRaster temp;
	std::string crsin = getCRS();
	if ((crsin == "") || (crs == "")) {
		temp.setError("insufficient crs info");	
		return temp;
	} else if (crsin == crs) {
		temp.setError("crs are the same");	
		return temp;
	}

	std::vector<std::vector<double>> p = extent.asPoints();
	SpatMessages msg = transform_coordinates(p[0], p[1], crsin, crs);
	if (msg.has_error) {
		temp.msg = msg;
		return temp;
	}
	if (msg.has_warning) {
		msg.setError("cannot do this");
		temp.msg = msg;
		return temp;
	}
	SpatExtent e(vmin(p[0], false), vmax(p[0], false), vmin(p[1], false), vmax(p[1], false));

	temp = SpatRaster(nrow(), ncol(), nlyr(), e, crs);
	return warp(temp, method, opt);
}


