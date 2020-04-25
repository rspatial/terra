#include "spatVector.h"
#include "distance.h"


SpatVector SpatVector::point_buffer(double d, unsigned quadsegs) { 

	std::vector<std::vector<double>> xy = coordinates();
	SpatVector out;
	out.lyr.srs = lyr.srs;
	size_t n = quadsegs * 4;	
	std::vector<double> px(n);
	std::vector<double> py(n);
	double step = 360.0 / n;
	SpatGeom g(polygons);
	g.addPart(SpatPart(0, 0));
	size_t npts = size();

	if (is_lonlat()) {
		double a=6378137;
		double f=1/298.257223563;
		std::vector<double> brng(n);
		for (size_t i=0; i<n; i++) {
			brng[i] = i * step;
		}
		for (size_t i=0; i<npts; i++) {
			std::vector<std::vector<double>> dp = destpoint_lonlat(xy[0][i], xy[1][i], brng, d, a, f);
			g.setPart(SpatPart(dp[0], dp[1]), 0);
			out.addGeom(g);
		}

	} else {
		std::vector<double> cosb(n);
		std::vector<double> sinb(n);
		for (size_t i=0; i<n; i++) {
			double brng = i * step;
			brng = toRad(brng);
			cosb[i] = d * cos(brng);
			sinb[i] = d * sin(brng);
		}
		for (size_t i=0; i<npts; i++) {
			for (size_t j=0; j<n; j++) {
				px[j] = xy[0][i] + cosb[j];
				py[j] = xy[1][i] + sinb[j];
			}
			g.setPart(SpatPart(px, py), 0);
			out.addGeom(g);
		}
	}
	return(out);
}




SpatVector SpatVector::buffer(double d, unsigned segments, unsigned capstyle){
	SpatVector out;
	std::string vt = type();
	if (vt != "points") {
		out.setError("must be points");
		return out;
	}
	if ((vt == "points") && (d <= 0)) {
		out.setError("buffer size must be >= 0 with points");
		return out;
	}
	out = point_buffer(d, segments); 
	return out;
}

