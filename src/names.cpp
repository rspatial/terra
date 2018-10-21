using namespace std;
#include "spat.h"
#include "util.h"


template <typename T>
vector<size_t> sort_indexes(const vector<T> &v) {

  // initialize original index locations
  vector<size_t> idx(v.size());
  iota(idx.begin(), idx.end(), 0);

  // sort indexes based on comparing values in v
  std::sort(idx.begin(), idx.end(),
       [&v](size_t i1, size_t i2) {return v[i1] < v[i2];});

  return idx;
}

void make_valid(std::vector<string> &s) {
    for (size_t i=0; i<s.size(); i++) {
        lrtrim(s[i]);
        if (s[i] == "") s[i] = "X";
        if (isdigit(s[i][0])) s[i] = "X" + s[i];
//        if ((s[i][0] == ".") & (s[i].size() > 1)) {
//			if (isdigit(s[i][1])) s[i] = "X" + s[i];
//		}
        std::replace(s[i].begin(), s[i].end(), ' ', '.');
    }
}

void make_unique(std::vector<string> &s) {
    std::vector<long unsigned> x = sort_indexes(s);
    std::sort(s.begin(), s.end());
    std::vector<string> ss = s;
    unsigned j = 1;
    for (size_t i=1; i<s.size(); i++) {
        if (s[i] == s[i-1]) {
            ss[i-1] = s[i-1] + "_" + to_string(j);
            ss[i] = s[i] + "_" + to_string(j + 1);
            j++;
        } else {
            j = 1;
        }
    }
    for (size_t i=0; i<s.size(); i++) {
        s[x[i]] = ss[i];
    }
}


std::vector<string> SpatRaster::getNames() {
	std::vector<string> x;
	for (size_t i=0; i<source.size(); i++) {
		x.insert(x.end(), source[i].names.begin(), source[i].names.end());
	}
	return(x);
}


bool SpatRaster::setNames(std::vector<string> names) {
	if (names.size() != nlyr()) {
		return false;
	} else {
        make_valid(names);
        make_unique(names);
        size_t begin=0;
        size_t end;
        for (size_t i=0; i<source.size(); i++)	{
            end = begin + source[i].nlyr;
            source[i].names = std::vector<string> (names.begin() + begin, names.begin() + end);
            begin = end;
        }
        return true;
	}
}

