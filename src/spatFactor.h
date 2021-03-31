#include <vector>
#include <string>

//#include "NA.h"
#include <algorithm>

template <typename T>
std::vector<T> unique_values(std::vector<T> d) {
  std::sort(d.begin(), d.end());
  d.erase(std::unique(d.begin(), d.end()), d.end());

  //d.erase(std::remove(d.begin(), d.end(), na), d.end());  
  return d;
}

/*
std::string double_2_string(double x) { 
  std::string s = std::to_string (x);
  s.erase( s.find_last_not_of('0') + 1, std::string::npos );
  s.erase( s.find_last_not_of('.') + 1, std::string::npos );
  return s;
}


std::vector<std::string> double_2_string(const std::vector<double> &x, std::string prep) { 
  std::vector<std::string> out(x.size());
  for (size_t i=0; i<x.size(); i++) {
    out[i] = prep + double_2_string (x[i]);
  }
  return out;
}
*/

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


class SpatFactor {
public:
	std::vector<unsigned> v;
	std::vector<unsigned> levels;
	std::vector<std::string> labels;
  
	size_t size() { return v.size(); }
	
	void compute_levels() {
		levels = unique_values(v);
		labels = string_values(levels);
	}

	template <typename T>
	void set_values(std::vector<T> _values) {
    std::vector<T> u = unique_values(_values);		
	  size_t n = _values.size();
	  size_t un = u.size();
	  
	labels = string_values(u);
    levels.resize(un);
    std::iota(levels.begin(), levels.end(), 0);
    v.resize(n);
    for (size_t i=0; i<n; i++) {
      for (size_t j=0; j<un; j++) {
        if (u[j] == _values[i]) {
          v[i] = j;
        }
      }  
    }
	}

	bool set_levels(std::vector<unsigned> _levels, std::vector<std::string> _labels) {
	  if (_levels.size() == _labels.size()) {
		  levels = _levels;
		  labels = _labels;
	  } else {
	    return false;
	  }
	  return true;
	}
	
	bool set_labels(std::vector<std::string> _labels) {
	  if (levels.size() == _labels.size()) {
	    labels = _labels;
	  } else {
	    return false;
	  }
	  return true;
	}
	
	void resize(size_t n) {
		//v.resize(n, NA<unsigned>::value);
	  v.resize(n);
	}
	void reserve(size_t n) {
	  v.reserve(n);
	}	
	template <typename T>
	  SpatFactor(std::vector<T> _v) {
	   set_values(_v);
	}
	
	
	SpatFactor subset(std::vector<unsigned> i) {
	  SpatFactor out;
	  out.levels = levels;
	  out.labels = labels;
	  size_t n = i.size();
	  out.reserve(n); 
	  for (size_t j=0; j<n; j++) {
	    out.v[j] = v[i[j]]; 
	  }
	  return out;
	}
	
	
	SpatFactor(){} ;
};

