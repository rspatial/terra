/*

#include <vector>
#include "spatRaster.h"
#include "vecmath.h"
#include "file_utils.h"

#ifdef useGDAL
	#include "crs.h"
	#include "warp.h"
#endif


void SpatRaster::project3(SpatRaster &out, std::string method, SpatOptions &opt) {

	std::string src = source[0].filename;
	std::vector<std::string> crsout = out.getSRS();
	std::vector<double> e = out.extent.asVector();
	std::vector<std::string> warpops = {"-t_srs", crsout[1], "-overwrite", 
					"-ts", std::to_string(out.ncol()), std::to_string(out.nrow()), 
					"-te", std::to_string(e[0]), std::to_string(e[2]), std::to_string(e[1]), std::to_string(e[3]), "-r", method, "-ovr", "NONE"};
	std::string dst = opt.filename;
	if (dst=="") {
		dst = tempFile(opt.get_tempdir(), ".tif");
	} else {
		if (opt.overwrite) {
			warpops.push_back("-overwrite");
		}
	}
	std::vector<std::string> empty;
	gdalwarp(src, dst, warpops, empty, empty);
	out.addWarning(dst);
	out = SpatRaster(dst);
}

*/

/*
void SpatRaster::project3(SpatRaster &out, std::string method, SpatOptions &opt) {

	unsigned nc = out.ncol();
	std::string crsin = getCRS();
	std::string crsout = out.getCRS();

  	if (!out.writeStart(opt)) { return; }
	for (size_t i = 0; i < out.bs.n; i++) {
        unsigned firstcell = out.cellFromRowCol(out.bs.row[i], 0);
		unsigned lastcell  = out.cellFromRowCol(out.bs.row[i]+out.bs.nrows[i]-1, nc-1);
		std::vector<double> cells(1+lastcell-firstcell);
		std::iota (std::begin(cells), std::end(cells), firstcell);
        std::vector<std::vector<double>> xy = out.xyFromCell(cells);
		#ifdef useGDAL
		out.msg = transform_coordinates(xy[0], xy[1], crsout, crsin);
		#else
		out.setError("GDAL is needed for crs transformation, but not available");
		return out;
		#endif
		std::vector<std::vector<double>> v = extractXY(xy[0], xy[1], method);
		if (!out.writeValues2(v, out.bs.row[i], out.bs.nrows[i], 0, out.ncol())) return;
	}
	out.writeStop();
}

*/

/*

SpatRaster SpatRaster::project2(SpatRaster &x, std::string method, SpatOptions &opt) {

	unsigned nl = nlyr();
	SpatRaster out = x.geometry(nl);
	out.setNames(getNames());
	std::vector<std::string> f {"bilinear", "ngb"};
	if (std::find(f.begin(), f.end(), method) == f.end()) {
		out.setError("unknown project method");
		return out;
	}
	if (!hasValues()) {
		return out;
	}
	
	if (srs.is_equal(out.srs))  {
		out.setError("rasters have the same crs");
		return out;
	}
	if (srs.is_empty() || out.srs.is_empty()) {
		out.setError("rasters must have a known crs");
		return out;
	}

// Need to check the resolution of the output to see
// is aggregation is needed

	//unsigned xq = x.xres() / xres();
	//unsigned yq = x.yres() / yres();
	//if (std::max(xq, yq) > 1) {
	//  SpatRaster xx;
	//	xq = xq == 0 ? 1 : xq;
	//	yq = yq == 0 ? 1 : yq;
	//	std::vector<unsigned> agf = {yq, xq, 1};
	//	SpatOptions agopt;
	//	if (method == "bilinear") {
	//		xx = aggregate(agf, "mean", true, agopt);
	//	} else {
	//		xx = aggregate(agf, "modal", true, agopt);
	//	}
	//  xx.project3(out, method, opt);
	//} else {
	
	project3(out, method, opt);	
	// }
	return out;
}



SpatRaster SpatRaster::project1(std::string newcrs, std::string method, SpatOptions &opt) {

	SpatRaster out;

	#ifndef useGDAL
	out.setError("GDAL is not available");
	return out;
	#else 
	
	if (srs.is_empty() || (newcrs == "")) {
		out.setError("insufficient crs info");	
		return out;
//	} else if (oldcrs == newcrs) {
//		out.setError("input and output crs are the same");	
//		return out;
	}

	std::vector<std::string> oldcrs = getSRS();
	std::vector<std::vector<double>> p = extent.asPoints();
	out.msg = transform_coordinates(p[0], p[1], oldcrs[1], newcrs);
	if (out.hasError()) {
		out.msg = msg;
		return out;
	}
	if (out.hasWarning()) {
		out.setError("cannot do this");
		return out;
	}
	SpatExtent e(vmin(p[0], false), vmax(p[0], false), vmin(p[1], false), vmax(p[1], false));

	out = SpatRaster(nrow(), ncol(), nlyr(), e, newcrs);
	project3(out, method, opt);
	return out;
	#endif	
}

*/

