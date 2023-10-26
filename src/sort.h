
#ifndef SORT_GUARD
#define SORT_GUARD

#include <numeric>

std::vector<std::size_t> sort_order_nan_a(const std::vector<double> &x);
std::vector<std::size_t> sort_order_nan_d(const std::vector<double> &x);
std::vector<std::size_t> sort_order_nal_a(const std::vector<long> &x);
std::vector<std::size_t> sort_order_nal_d(const std::vector<long> &x);
std::vector<std::size_t> sort_order_nas_a(const std::vector<std::string> &x);
std::vector<std::size_t> sort_order_nas_d(const std::vector<std::string> &x);

template <typename T>
std::vector<std::size_t> sort_order_a(const std::vector<T> &x){
	std::vector<std::size_t> p(x.size());
	std::iota(p.begin(), p.end(), 0);
	std::sort(p.begin(), p.end(),
            [&](std::size_t i, std::size_t j){ return (x[i] < x[j]); });
	return p;
}

template <typename T>
std::vector<std::size_t> sort_order_d(const std::vector<T> &x){
	std::vector<std::size_t> p(x.size());
	std::iota(p.begin(), p.end(), 0);
	std::sort(p.begin(), p.end(),
            [&](std::size_t i, std::size_t j){ return (x[i] > x[j]); });
	return p;
}


template <typename T>
void permute(std::vector<T> &x, const std::vector<std::size_t> &p) {
	std::vector<bool> done(x.size());
	for (std::size_t i = 0; i < x.size(); ++i)  {
		if (done[i]) {
			continue;
		}	
		done[i] = true;
		size_t prev_j = i;
		size_t j = p[i];
		while (i != j) {
			std::swap(x[prev_j], x[j]);
			done[j] = true;
			prev_j = j;
			j = p[j];
		}
	}
}

#endif
