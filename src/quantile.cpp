#include <algorithm>
#include <cmath>
#include <vector>
#include <algorithm>
#include "SpatRaster.h"

#include "vecmath.h"
#include "math_utils.h"
#include "string_utils.h"


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
		unsigned x1 = floor(x);
		unsigned x2 = ceil(x);
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

