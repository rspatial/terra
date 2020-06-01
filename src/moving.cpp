#include <vector>
#include <algorithm>
#include <cmath>

double median(const std::vector<double>& v) {
	size_t n = v.size();
	std::vector<double> vv;
	vv.reserve(n);
	for (size_t i=0; i<n; i++) {
        if (!std::isnan(v[i])) {
            vv.push_back(v[i]);
        }
	}
	n = vv.size();
	if (n == 0) {
		return(NAN);
	}
	size_t n2 = n / 2;
	std::nth_element(vv.begin(), vv.begin()+n2, vv.end());
	double med = vv[n2];
	return med;
}



std::vector<double> movingMedian(const std::vector<double> &x, size_t n) {
	std::vector<double> out(x.size());
	std::vector<double> d(n, NAN);
	size_t half = (n/2);
	size_t half1 = half+1;
	// fill left side
	for (size_t i=0; i<half; i++) {
		for (size_t j=0; j< (half1+i); j++) {
			d[j] = x[j];
		}
		out[i] = median(d);
	}
	// middle
	size_t maxn = out.size() - half;
	std::vector<double> v;
	for (size_t i=half; i<maxn; i++) {
		v = std::vector<double>(x.begin()+i-half, x.begin()+i+half1);
 		out[i] = median(v);
	}
	// right side
	int j=0;
	for (size_t i=maxn; i<out.size(); i++) {
		v[j++] = NAN;
		out[i] = median(v);
	}
	return(out);
}

