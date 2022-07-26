#include <vector>
#include <string>
#include <numeric>
#include <algorithm>

#include "spatFactor.h"


template <typename T>
std::vector<T> unique_values(std::vector<T> d) {
  std::sort(d.begin(), d.end());
  d.erase(std::unique(d.begin(), d.end()), d.end());
  //d.erase(std::remove(d.begin(), d.end(), na), d.end());
  return d;
}

std::vector<std::string> string_values(std::vector<std::string>& v) {
  return v;
}

template <typename T>
std::vector<std::string> string_values(std::vector<T>& v) {
  std::vector<std::string> result;
  std::transform(std::begin(v),
                 std::end(v), std::back_inserter(result),
                 [](T d) {
                    std::string s = std::to_string(d);
                     s.erase( s.find_last_not_of('0') + 1, std::string::npos );
                     s.erase( s.find_last_not_of('.') + 1, std::string::npos );
                     return s;
                   }
  );
  return result;
}


//void SpatFactor::compute_levels() {
//	levels = unique_values(v);
//	labels = string_values(levels);
//}


SpatFactor::SpatFactor(std::vector<unsigned> _values, std::vector<std::string> _labels) {
	v = _values;
	labels = _labels;
	// is this needed?
	//levels.resize(labels.size());
	//std::iota(levels.begin(), levels.end(), 0);
}

SpatFactor::SpatFactor(std::vector<unsigned> _values) {
	std::vector<unsigned> u = unique_values(_values);
	size_t n = _values.size();
	size_t un = u.size();
		
	labels = string_values(u);
	//levels.resize(un);
	//std::iota(levels.begin(), levels.end(), 0);
	v.resize(n);
	for (size_t i=0; i<n; i++) {
		for (size_t j=0; j<un; j++) {
			if (u[j] == _values[i]) {
				v[i] = j;
			}
		}
	}
}

SpatFactor::SpatFactor(std::vector<std::string> _values) {
	std::vector<std::string> u = unique_values(_values);
	size_t n = _values.size();
	size_t un = u.size();
		
	labels = string_values(u);
	//levels.resize(un);
	//std::iota(levels.begin(), levels.end(), 0);
	v.resize(n);
	for (size_t i=0; i<n; i++) {
		for (size_t j=0; j<un; j++) {
			if (u[j] == _values[i]) {
				v[i] = j;
			}
		}
	}
}

/*
bool SpatFactor::set_levels(std::vector<unsigned> _levels, std::vector<std::string> _labels) {
	if (_levels.size() == _labels.size()) {
		levels = _levels;
		labels = _labels;
	} else {
		return false;
	}
	return true;
}
*/

bool SpatFactor::set_labels(std::vector<std::string> _labels) {
	//if (levels.size() == _labels.size()) {
	labels = _labels;
	//} else {
	//	return false;
	//}
	return true;
}



SpatFactor SpatFactor::subset(std::vector<unsigned> i) {
	SpatFactor out;
	out.labels = labels;
	size_t n = i.size();
	out.reserve(n);
	for (size_t j=0; j<n; j++) {
		out.v.push_back(v[i[j]]);
	}
	return out;
}

std::string SpatFactor::getLabel(size_t i) {
	if (i < v.size()) {
		unsigned j = v[i] - 1;
		if (j < labels.size()) {
			return labels[j];
		}
	}
	return "";
}

std::vector<std::string> SpatFactor::getLabels() {
	std::vector<std::string> out;
	size_t n = v.size();
	size_t m = labels.size() + 1;
	out.reserve(n);

	for (size_t i=0; i<n; i++) {
		if (v[i] < m) {
			out.push_back(labels[v[i]-1]);
		} else {
			out.push_back("");
		}
	}
	return out;
}



