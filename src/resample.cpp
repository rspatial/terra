#include <vector>
#include "spatRaster.h"
#include "vecmath.h"



void SpatRaster::resample2(SpatRaster &out, const std::string &method, SpatOptions &opt) {

	unsigned nc = out.ncol();
  	if (!out.writeStart(opt)) { return; }
	for (size_t i = 0; i < out.bs.n; i++) {
        unsigned firstcell = out.cellFromRowCol(out.bs.row[i], 0);
		unsigned lastcell  = out.cellFromRowCol(out.bs.row[i]+out.bs.nrows[i]-1, nc-1);
		std::vector<double> cells(1+lastcell-firstcell);
		std::iota (std::begin(cells), std::end(cells), firstcell);
        std::vector<std::vector<double>> xy = out.xyFromCell(cells);
		std::vector<std::vector<double>> v = extractXY(xy[0], xy[1], method);
		if (!out.writeValues2(v, out.bs.row[i], out.bs.nrows[i], 0, out.ncol())) return;
	}
	out.writeStop();

}


SpatRaster SpatRaster::resample1(SpatRaster &x, const std::string &method, SpatOptions &opt) {

	unsigned nl = nlyr();
	SpatRaster out = x.geometry(nl);
	out.setNames(getNames());
	std::vector<std::string> f {"bilinear", "ngb"};
	if (std::find(f.begin(), f.end(), method) == f.end()) {
		out.setError("unknown resample method");
		return out;
	}
	if (!hasValues()) {
		return out;
	}
	
	if ((!srs.is_empty()) && (!out.srs.is_empty())) {
		if (!srs.is_equal(out.srs)) {
			out.addWarning("Rasters have different crs");
		}
	}
	
	unsigned xq = x.xres() / xres();
	unsigned yq = x.yres() / yres();
	if (std::max(xq, yq) > 1) {
		SpatRaster xx;
		xq = xq == 0 ? 1 : xq;
		yq = yq == 0 ? 1 : yq;
		std::vector<unsigned> agf = {yq, xq, 1};
		SpatOptions agopt;
		if (method == "bilinear") {
			xx = aggregate(agf, "mean", true, agopt);
		} else {
			xx = aggregate(agf, "modal", true, agopt);
		}
		xx.resample2(out, method, opt);
	} else {
		if ((x.xres() == xres()) && (x.yres() == yres()) & (extent.compare(x.extent, "==", std::min(xres(), yres())/1000))) {
			out = *this;
			out.extent = x.extent;
		} else {
			resample2(out, method, opt);
		}
	}
	return out;
}

