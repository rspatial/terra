
template <typename T>
void recycle(std::vector<T> &v, unsigned n) {
	size_t s = v.size();
	v.resize(n);
	for (size_t i=s; i<n; i++) {
		 v[i] = v[i % s];
	}
}


template <typename T>
void recycle(std::vector<T> &x, std::vector<T> &y) {
	size_t xsize = x.size();
	size_t ysize = y.size();
	if (xsize != ysize) {
		size_t n = std::max(xsize, ysize);
		if (xsize > ysize) {
			size_t s = y.size();
			y.resize(n);
			for (size_t i=ysize; i < n; i++) {
				y[i] = y[i % s];
			} 				
		} else {
			size_t s = x.size();
			x.resize(n);
			for (size_t i=xsize; i<n; i++) {
				x[i] = x[i % s];
			} 				
		}
	}
}

