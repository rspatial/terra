// Copyright (c) 2018-2023  Robert J. Hijmans
//
// This file is part of the "spat" library.
//
// spat is free software: you can redistribute it and/or modify it
// under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 2 of the License, or
// (at your option) any later version.
//
// spat is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with spat. If not, see <http://www.gnu.org/licenses/>.

#include <vector>
#include <cmath>
#include <algorithm>
#include <string>
#include <functional>
#include <map>


double median_se_rm(const std::vector<double>& v, size_t s, size_t e) {
	size_t n = (e - s) + 1;
	std::vector<double> vv;
	vv.reserve(n);
	for (size_t i=s; i<e; i++) {
        if (!std::isnan(v[i])) {
            vv.push_back(v[i]);
        }
	}
	n = vv.size();
	if (n == 0) {
		return(NAN);
	}
	if (n == 1) {
		return(vv[0]);
	}
	size_t n2 = n / 2;
	if (n % 2) {
		std::nth_element(vv.begin(), vv.begin()+n2, vv.end());
		return vv[n2];
	} else {
		std::sort(vv.begin(), vv.end());
		return (vv[n2] + vv[n2-1]) / 2;
	}
}


double median_se(const std::vector<double>& v, size_t s, size_t e) {
	size_t n = (e - s) + 1;
	std::vector<double> vv;
	vv.reserve(n);
	for (size_t i=s; i<e; i++) {
        if (std::isnan(v[i])) {
			return(NAN);
		} else {
            vv.push_back(v[i]);
        } 
	}
	n = vv.size();
	if (n == 0) {
		return(NAN);
	}
	if (n == 1) {
		return(vv[0]);
	}
	size_t n2 = n / 2;
	if (n % 2) {
		std::nth_element(vv.begin(), vv.begin()+n2, vv.end());
		return vv[n2];
	} else {
		std::sort(vv.begin(), vv.end());
		return (vv[n2] + vv[n2-1]) / 2;
	}
}



double sum_se_rm(const std::vector<double>& v, size_t s, size_t e) {
	double x = v[s];
	for (size_t i=(s+1); i<e; i++) {
		if (std::isnan(x)) {
			x = v[i];
		} else if (!std::isnan(v[i])) {
			x += v[i];
		}
	}
	return x;
}

double sum_se(const std::vector<double>& v, size_t s, size_t e) {
	double x = v[s];
	if (std::isnan(x)) {
		return(x);
	}
	for (size_t i=(s+1); i<e; i++) {
		if (std::isnan(v[i])) {
			return NAN;
			break;
		} else {
			x += v[i];
		}
	}
	return x;
}


double sum2_se_rm(const std::vector<double>& v, size_t s, size_t e) {
	double x = v[s] * v[s];
	for (size_t i=(s+1); i<e; i++) {
		if (std::isnan(x)) {
			x = v[i] * v[i];
		} else if (!std::isnan(v[i])) {
			x += v[i] * v[i];
		}
	}
	return x;
}

double sum2_se(const std::vector<double>& v, size_t s, size_t e) {
	double x = v[s] * v[s];
	if (std::isnan(x)) return(x);
	for (size_t i=(s+1); i<e; i++) {
		if (std::isnan(v[i])) {
			return(v[i]);
		} else {
			x += v[i] * v[i];
		}
	}
	return x;
}



double prod_se_rm(const std::vector<double>& v, size_t s, size_t e) {
	double x = v[s];
	for (size_t i=(s+1); i<e; i++) {
		if (std::isnan(x)) {
			x = v[i];
		} else if (!std::isnan(v[i])) {
			x *= v[i];
		}
	}
	return x;
}


double prod_se(const std::vector<double>& v, size_t s, size_t e) {
	double x = v[s];
	if (std::isnan(x)) {
		return(NAN);
	}
	for (size_t i=(s+1); i<e; i++) {
		if (std::isnan(v[i])) {
			return(NAN);
		} else {
			x *= v[i];
		}
	}
	return x;
}



double mean_se_rm(const std::vector<double>& v, size_t s, size_t e) {
	double x = 0;
	unsigned d = 0;
	for (size_t i=s; i<e; i++) {
		if (!std::isnan(v[i])) {
			x += v[i];
			d++;
		}
	}
	if (d > 0) {
		x /= (double) d;
	} else {
		x = NAN;
	}
	return x;
}


double mean_se(const std::vector<double>& v, size_t s, size_t e) {
	double x = 0;
	unsigned d = 0;
	for (size_t i=s; i<e; i++) {
		if (std::isnan(v[i])) {
			return(NAN);
		} else {
			x += v[i];
			d++;
		}
	}
	if (d > 0) {
		x /= (double) d;
	} else {
		x = NAN;
	}
	return x;
}


double sd_se_rm(const std::vector<double>& v, size_t s, size_t e) {
	double m = mean_se_rm(v, s, e);
	if (std::isnan(m)) return m;
	double x = 0;
	size_t n = 0;
	for (size_t i=s; i<e; i++) {
		if (!std::isnan(v[i])) {
			double d = (v[i] - m);
			x += (d * d);
			n++;
		}
	}
	n--;
	x = sqrt(x / n);
	return x;
}


double sd_se(const std::vector<double>& v, size_t s, size_t e) {
	double m = mean_se(v, s, e);
	if (std::isnan(m)) return m;
	double x = 0;
	size_t n = 0;
	for (size_t i=s; i<e; i++) {
		double d = (v[i] - m);
		x += (d * d);
		n++;
	}
	n--;
	return sqrt(x / n);
}


double sdpop_se_rm(const std::vector<double>& v, size_t s, size_t e) {
	double m = mean_se_rm(v, s, e);
	if (std::isnan(m)) return m;
	double x = 0;
	size_t n = 0;
	for (size_t i=s; i<e; i++) {
		if (!std::isnan(v[i])) {
			double d = (v[i] - m);
			x += d * d;
			n++;
		}
	}
	return sqrt(x / n);
}



double sdpop_se(const std::vector<double>& v, size_t s, size_t e) {
	double m = mean_se(v, s, e);
	if (std::isnan(m)) return m;
	double x = 0;
	size_t n = 0;
	for (size_t i=s; i<e; i++) {
		double d = (v[i] - m);
		x += d * d;
		n++;
	}
	return sqrt(x / n);
}


double min_se_rm(const std::vector<double>& v, size_t s, size_t e) {
	double x = v[s];
	for (size_t i=(s+1); i<e; i++) {
		if (!std::isnan(v[i])) {
			if (std::isnan(x)) {
				x = v[i];
			} else {
				x = std::min(x, v[i]);
			}
		}
	}
	return x;
}


double min_se(const std::vector<double>& v, size_t s, size_t e) {
	double x = v[s];
	if (std::isnan(x)) return x;
	for (size_t i=(s+1); i<e; i++) {
		if (std::isnan(v[i])) {
			return NAN;
		} else {
			x = std::min(x, v[i]);
		}
	}
	return x;
}


double max_se_rm(const std::vector<double>& v, size_t s, size_t e) {
	double x = v[s];
	for (size_t i=(s+1); i<e; i++) {
		if (!std::isnan(v[i])) {
			if (std::isnan(x)) {
				x = v[i];
			} else {
				x = std::max(x, v[i]);
			}
		}
	}
	return x;
}


double max_se(const std::vector<double>& v, size_t s, size_t e) {
	double x = v[s];
	if (std::isnan(x)) return x;
	for (size_t i=(s+1); i<e; i++) {
		if (std::isnan(v[i])) {
			return NAN;
		} else {
			x = std::max(x, v[i]);
		}
	}
	return x;
}


double first_se_rm(std::vector<double>& v, size_t s, size_t e) {
	for (size_t i=s; i<e; i++) {
		if (!std::isnan(v[i])) {
			return v[i];
		}
	} 
	return NAN;
}

double first_se(std::vector<double>& v, size_t s, size_t e) {
	return v[s];
}


double which_se_rm(const std::vector<double>& v, size_t s, size_t e) {
	for (size_t i=s; i<e; i++) {
		if ((!std::isnan(v[i])) && (v[i] != 0)) {
			return (i+1) - s; // +1 for R
		}
	}
	return NAN;
}

double which_se(const std::vector<double>& v, size_t s, size_t e) {
	return which_se_rm(v, s, e);
}


double whichmin_se_rm(const std::vector<double>& v, size_t s, size_t e) {
	double x = v[s];
	double out = s;
	if (std::isnan(x)) {
		out = NAN;
	} 
	for (size_t i=(s+1); i<e; i++) {
		if (!std::isnan(v[i])) {
			if (std::isnan(out)) {
				x = v[i];
				out = i;
			} else if (v[i] < x) {
				x = v[i];
				out = i;
			}
		}
	} 
	out++;  // +1 for R
	return (out - s);
}

double whichmin_se(const std::vector<double>& v, size_t s, size_t e) {
	return whichmin_se_rm(v, s, e);
}


double whichmax_se_rm(const std::vector<double>& v, size_t s, size_t e) {
	double x = v[s];
	double out = s;
	if (std::isnan(x)) {
		out = NAN;
	}
	for (size_t i=(s+1); i<e; i++) {
		if (!std::isnan(v[i])) {
			if (std::isnan(out)) {
				x = v[i];
				out = i;
			} else if (v[i] > x) {
				x = v[i];
				out = i;
			}
		}
	}
	out++;  // +1 for R
	return (out - s); 
}

double whichmax_se(const std::vector<double>& v, size_t s, size_t e) {
	return whichmax_se_rm(v, s, e);
}


double all_se_rm(const std::vector<double>& v, size_t s, size_t e) {
	double x = 1;
	for (size_t i=s; i<e; i++) {
		if (!std::isnan(v[i])) {
			if (v[i] == 0) {
				x = v[i];
				return x;
			} 
		}
    }
	return x;
}


double all_se(const std::vector<double>& v, size_t s, size_t e) {
	double x = 1;
    for (size_t i=s; i<e; i++) {
        if (std::isnan(v[i]) || (v[i] == 0)) {
            x = v[i];
			return x;
		}
	}
	return x;
}


double any_se_rm(const std::vector<double>& v, size_t s, size_t e) {
	double x = NAN;
	for (size_t i=s; i<e; i++) {
		if (!std::isnan(v[i])) {
			if (v[i] != 0) {
				x = 1;
				break;
			} else {
				x = 0;
			}
		}
	}
	return x;
}

double any_se(const std::vector<double>& v, size_t s, size_t e) {
	double x = 0;
	for (size_t i=s; i<e; i++) {
		if (std::isnan(v[i])) {
			return(NAN);
		}
		if (v[i] != 0) {
			x = 1;
		}
	}
	return x;
}


std::vector<double> range_se_rm(std::vector<double>& v, size_t s, size_t e) {

	std::vector<double> x = { v[s], v[s] };
	for (size_t i=(s+1); i<e; i++) {
		if (!std::isnan(v[i])) {
			if (std::isnan(x[0])) {
				x[0] = v[i];
				x[1] = v[i];
			} else {
				x[0] = std::min(x[0], v[i]);
				x[1] = std::max(x[1], v[i]);
			}
		}
	} 
	return x;
}


std::vector<double> range_se(std::vector<double>& v, size_t s, size_t e) {

	std::vector<double> x = { v[s], v[s] };
	if (!std::isnan(x[0])) { return x; }

	for (size_t i=(s+1); i<e; i++) {
		if (std::isnan(v[i])) {
			x[0] = NAN;
			x[1] = NAN;
			break;
		} else {
			x[0] = std::min(x[0], v[i]);
			x[1] = std::max(x[1], v[i]);
		}
	}
	return x;
	
}



double modal_se_rm(std::vector<double>& v, size_t s, size_t e) {

	std::map<double, size_t> count;
	for_each( v.begin()+s, v.begin()+e, [&count]( double val ){
			if(!std::isnan(val)) count[val]++;
		}
	);

	if (count.size() == 0) return NAN;
	
    std::map<double, size_t>::iterator mode =	
        std::max_element(count.begin(), count.end(),[] (const std::pair<double, size_t>& a, 
		const std::pair<double, size_t>& b)->bool{ return a.second < b.second; } );
		
    return mode->first;
}

double modal_se(std::vector<double>& v, size_t s, size_t e) {
	std::map<double, size_t> count;
	for(size_t i=s; i<e; i++) {
		if (std::isnan(v[i])) {
			return NAN;
		} else {
			count[v[i]]++;
		}
	}

    std::map<double, size_t>::iterator mode =	
        std::max_element(count.begin(), count.end(),[] (const std::pair<double, size_t>& a, 
		const std::pair<double, size_t>& b)->bool{ return a.second < b.second; } );
		
    return mode->first;
}


double isna_se(const std::vector<double>& v, size_t s, size_t e) {
	double x = 0;
	for (size_t i=s; i<e; i++) {
		if (std::isnan(v[i])) {
			x++;
		}
	}
	return x;
}


double isnotna_se(const std::vector<double>& v, size_t s, size_t e) {
	double x = 0;
	for (size_t i=s; i<e; i++) {
		if (!std::isnan(v[i])) {
			x++;
		}
	}
	return x;
}


void cumsum_se_rm(std::vector<double>& v, size_t s, size_t e) {
    for (size_t i=(s+1); i<e; i++) {
        if (std::isnan(v[i])) {
            v[i] = v[i-1];
        } else if (!std::isnan(v[i-1])){
            v[i] += v[i-1];
        }
    }
}

void cumsum_se(std::vector<double>& v, size_t s, size_t e) {
    for (size_t i=(s+1); i<e; i++) {
        if (std::isnan(v[i]) || std::isnan(v[i-1])) {
            v[i] = NAN;
        } else {
            v[i] += v[i-1];
        }
    }
}



void cumprod_se_rm(std::vector<double>& v, size_t s, size_t e) {
    for (size_t i=(s+1); i<e; i++) {
        if (std::isnan(v[i])) {
            v[i] = v[i-1];
        } else if (!std::isnan(v[i-1])){
            v[i] *= v[i-1];
        }
    }
}

void cumprod_se(std::vector<double>& v, size_t s, size_t e) {
    for (size_t i=(s+1); i<e; i++) {
        if (std::isnan(v[i]) || std::isnan(v[i-1])) {
            v[i] = NAN;
        } else {
			v[i] *= v[i-1]; 
        }
    }
}



void cummax_se_rm(std::vector<double>& v, size_t s, size_t e) {
    for (size_t i=(s+1); i<e; i++) {
        if (std::isnan(v[i])) {
            v[i] = v[i-1];
        } else if (!std::isnan(v[i-1])){
           v[i] = std::max(v[i], v[i-1]);
        }
    }
}

void cummax_se(std::vector<double>& v, size_t s, size_t e) {
    for (size_t i=(s+1); i<e; i++) {
        if (std::isnan(v[i]) || std::isnan(v[i-1])) {
			v[i] = NAN;
        } else {
			v[i] = std::max(v[i], v[i-1]);
        }
    }
}


void cummin_se_rm(std::vector<double>& v, size_t s, size_t e) {
    for (size_t i=(s+1); i<e; i++) {
        if (std::isnan(v[i])) {
            v[i] = v[i-1];
        } else if (!std::isnan(v[i-1])){
            v[i] = std::min(v[i], v[i-1]);
        }
    }
}


void cummin_se(std::vector<double>& v, size_t s, size_t e) {
    for (size_t i=(s+1); i<e; i++) {
        if (std::isnan(v[i]) || std::isnan(v[i-1])) {
            v[i] = NAN;
        } else {
            v[i] = std::min(v[i], v[i-1]);
        }
    }
}



bool haveseFun(std::string fun) {
	std::vector<std::string> f {"sum", "mean", "median", "modal", "which", "which.min", "which.max", "min", "max", "prod", "any", "all", "sd", "std", "first"};
	auto it = std::find(f.begin(), f.end(), fun);
	if (it == f.end()) {
		return false;
	}
	return true;
}


bool getseFun(std::function<double(std::vector<double>&, size_t, size_t)> &fun,
				std::string fname, bool narm) {
					
	if (fname == "mean") {
		fun = narm ? mean_se_rm : mean_se;
	} else if (fname == "sum") {
		fun = narm ? sum_se_rm : sum_se;
	} else if (fname == "sum2") {
		fun = narm ? sum2_se_rm : sum2_se;
	} else if (fname == "min") {
		fun = narm ? min_se_rm : min_se;
	} else if (fname == "max") {
		fun = narm ? max_se_rm : max_se;
	} else if (fname == "median") {
		fun = narm ? median_se_rm : median_se;
	} else if (fname == "modal") {
		fun = narm ? modal_se_rm : modal_se;
	} else if (fname == "prod") {
		fun = narm ? prod_se_rm : prod_se;
	} else if (fname == "which") {
		fun = narm ? which_se_rm : which_se;
	} else if (fname == "which.min") {
		fun = narm ? whichmin_se_rm : whichmin_se;
	} else if (fname == "which.max") {
		fun = narm ? whichmax_se_rm : whichmax_se;
	} else if (fname == "any") {
		fun = narm ? any_se_rm : any_se;
	} else if (fname == "all") {
		fun = narm ? all_se_rm : all_se;
	} else if (fname == "sd") {
		fun = narm ? sd_se_rm : sd_se;
	} else if (fname == "std") {
		fun = narm ? sdpop_se_rm : sdpop_se;
	} else if (fname == "first") {
		fun = narm ? first_se_rm : first_se;
	} else {
		return false;
	}
	return true;
}



double wsum_se_rm(const std::vector<double>& v, const std::vector<double>& w, size_t s, size_t e) {
	if (w.size() == 0) return NAN;
	double x = 0;
	bool allna = true;
	for (size_t i=s; i<e; i++) {
		if (!std::isnan(v[i])) {
			x += v[i] * w[i];
			allna = false;
		}
	}
	if (allna) {
		return NAN;
	} else {
		return x;
	}
}


double wsum_se(const std::vector<double>& v, const std::vector<double>& w, size_t s, size_t e) {
	if (w.size() == 0) return NAN;
	double x = 0;
	for (size_t i=s; i<e; i++) {
		x += v[i] * w[i];
	}
	return x;
}


double wmean_se_rm(const std::vector<double>& v, const std::vector<double>& w, size_t s, size_t e) {
	double sv = 0;
	double sw = 0;
	for (size_t i=s; i<e; i++) {
		if ((!std::isnan(v[i])) && (!std::isnan(w[i]))) {
			sv += (v[i] * w[i]);
			sw += w[i];			
		}
	}
	return sv / sw;
}


double wmean_se(const std::vector<double>& v, const std::vector<double>& w, size_t s, size_t e) {
	double sv = 0;
	double sw = 0;
	for (size_t i=s; i<e; i++) {
		if (!std::isnan(w[i])) {
			sv += (v[i] * w[i]);
			sw += w[i];			
		}
	}
	return sv / sw;
}

double wmin_se_rm(const std::vector<double>& v, const std::vector<double>& w, size_t s, size_t e) {
	double x = NAN;
	for (size_t i=(s); i<e; i++) {
		if (!std::isnan(w[i]) && !std::isnan(v[i])) {
			if (std::isnan(x)) {
				x = v[i];
			} else {
				x = std::min(x, v[i]);
			}
		}
	}
	return x;
}

double wmin_se(const std::vector<double>& v, const std::vector<double>& w, size_t s, size_t e) {
	double x = NAN;
	for (size_t i=(s); i<e; i++) {
		if (!std::isnan(w[i])) {
			if (std::isnan(v[i])) {
				return NAN;
			}
			if (std::isnan(x)) {
				x = v[i];
			} else {
				x = std::min(x, v[i]);
			}
		}
	}
	return x;
}

double wmax_se_rm(const std::vector<double>& v, const std::vector<double>& w, size_t s, size_t e) {
	double x = NAN;
	for (size_t i=(s); i<e; i++) {
		if (!std::isnan(w[i]) && !std::isnan(v[i])) {
			if (std::isnan(x)) {
				x = v[i];
			} else {
				x = std::max(x, v[i]);
			}
		}
	}
	return x;
}

double wmax_se(const std::vector<double>& v, const std::vector<double>& w, size_t s, size_t e) {
	double x = NAN;
	for (size_t i=(s); i<e; i++) {
		if (!std::isnan(w[i])) {
			if (std::isnan(v[i])) {
				return NAN;
			}
			if (std::isnan(x)) {
				x = v[i];
			} else {
				x = std::max(x, v[i]);
			}
		}
	}
	return x;
}


bool haveseWFun(std::string fun) {
	std::vector<std::string> f {"sum", "mean", "min", "max"};
	auto it = std::find(f.begin(), f.end(), fun);
	if (it == f.end()) {
		return false;
	}
	return true;
}


bool getseWfun(std::function<double(std::vector<double>&, std::vector<double>&, size_t, size_t)> &fun, 
			std::string fname, bool narm) {
	if (fname == "mean") {
		fun = narm ? wmean_se_rm : wmean_se;
	} else if (fname == "sum") {
		fun = narm ? wsum_se_rm : wsum_se;
	} else if (fname == "min") {
		fun = narm ? wmin_se_rm : wmin_se;
	} else if (fname == "max") {
		fun = narm ? wmax_se_rm : wmax_se;
	} else {
		return false;
	}
	return true;
}

