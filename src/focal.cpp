// Copyright (c) 2018-2021  Robert J. Hijmans
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

#include "spatRaster.h"
#include "vecmath.h"


std::vector<double> rcValue(std::vector<double> &d, const int& nrow, const int& ncol, const unsigned& nlyr, const int& row, const int& col) {
  
	std::vector<double> out(nlyr, NAN);
	if ((row < 0) || (row > (nrow -1)) || (col < 0) || (col > (ncol-1))) {
		return out;
	} else {
		unsigned nc = nrow * ncol;
		unsigned cell = row * ncol + col;
		for (size_t i=0; i<nlyr; i++) {
			unsigned lcell = cell + i * nc;
			out[i] = d[lcell];
		}
	}
	return out;
}

// todo: three dimensional focal




std::vector<double> SpatRaster::focal_values(std::vector<unsigned> w, double fillvalue, int_64 row, int_64 nrows, SpatOptions &ops) {

	if (nlyr() > 1) {
		std::vector<unsigned> lyr = {0};
		SpatRaster s = subset(lyr, ops);
		s.focal_values(w, fillvalue, row, nrows, ops); 
	}

	if ((w[0] % 2 == 0) || (w[1] % 2 == 0)) {
		setError("weights matrix must have uneven sides");
		std::vector<double> d;
		return(d);
	}

	int_64 wr = w[0] / 2;
	int_64 nr = nrow();
	int_64 nc = ncol();
	wr = std::min(wr, nr-1);

	int_64 startrow = row-wr;
	startrow = startrow < 0 ? 0 : startrow;
	int_64 startoff = row-startrow;

	nrows = std::max((int_64)1, nrows);
	int_64 readnrows = nrows+startoff+wr;
	int_64 endoff = wr;
	if ((startrow+readnrows) > nr ) {
		readnrows = nr-startrow;
		endoff = readnrows - (nrows+startoff);
	}
	std::vector<double> d;
	readValues(d, startrow, readnrows, 0, nc);

//	get_focal(f, d, nrows, nc, w[0], w[1], offset, endoff, fillvalue);
//  get_focal(std::vector<double> &out, const std::vector<double> &d, int nrow, int ncol, int wrows, int wcols, int startoff, int endoff,  double fill) {

// ??
	//wr = std::min(wr, std::max((int_64)1, nrows-1));
	int_64 wc = w[1] / 2;
	wc = std::min(wc, nc-1);

	size_t n = nrows * nc * w[0] * w[1];
	std::vector<double> out(n, fillvalue);
	int_64 nrmax = nrows + startoff + endoff - 1;
	//int nrmax = d.size() / ncol - 1;
	size_t f = 0;
	const bool global = is_global_lonlat();

//Rcpp::Rcout << "sr" << startrow << " so " << startoff << " rnr " << readnrows << " wr " << wr << " wc " << wc << " nrows " << nrows << std::endl; 

	for (int_64 r=0; r < nrows; r++) {
		for (int_64 c=0; c < nc; c++) {
			for (int_64 i = -wr; i <= wr; i++) {
				int_64 row = r+startoff+i;
				if ((row < 0) || (row > nrmax)) {
					f += w[1];
				} else {
					size_t bcell = row * nc;
					for (int_64 j = -wc; j <= wc; j++) {
						int_64 col = c + j;
						if ((col >= 0) && (col < nc)) {
							out[f] = d[bcell+col];
						} else if (global) {
							if (col < 0) {
								col = nc + col;
							} else if (col >= nc) {
								col = col - nc;
							} 
							out[f] = d[bcell+col]; 
						}
						f++;
					}
				}
			}
		}
	}
	return out;
}





void focal_win_fun(const std::vector<double> &d, std::vector<double> &out, int nc, int srow, int nr,
                    std::vector<double> window, int wnr, int wnc, double fill, bool narm, bool naonly, bool naomit, bool expand, bool global, std::function<double(std::vector<double>&, bool)> fun) {

	out.resize(nc * nr);
	int hwc = wnc / 2;
	int hwr = wnr / 2;
	std::vector<bool> winNA(window.size(), false);
	for (size_t i=0; i<window.size(); i++) {
		if (std::isnan(window[i])) winNA[i] = true; 
	}

	int wr1 = wnr - 1;
	int wc1 = wnc - 1;
	int nc1 = nc - 1;
	
	bool checkNA = naonly || naomit;

	for (int r=0; r < nr; r++) {
		int rread = r+srow;
		for (int c=0; c < nc; c++) {
			size_t cell = r*nc + c;
			if (checkNA) {
				size_t readcell = rread*nc + c;
				if (naonly) {
					if (!std::isnan(d[readcell])) {
						out[cell] = d[readcell];
						continue;
					}					
				} else if (std::isnan(d[readcell])) {
					out[cell] = d[readcell];						
					continue;
				}
			}
			std::vector<double> v;
			v.reserve(wnr * wnc);
			for (int rr=0; rr<wnr; rr++) {
				int offr = wr1 - rr;
				for (int cc=0; cc < wnc; cc++)  {
					int offc = wc1 - cc; 
					//int wi = wnc * offr + offc;
					int wi = wnc * rr + cc;
					if (winNA[wi]) {
						continue;
					}
					int row = rread + hwr - offr;
					int col = c + hwc - offc;
					if (global) {
						col = col < 0 ? nc + col : col;
						col = col > nc1 ? col - nc : col;
						v.push_back(d[nc*row + col] * window[wi]);
					} else if (expand) {
						col = col < 0 ? 0 : col;
						col = col > nc1 ? nc1 : col;
						v.push_back(d[nc*row + col] * window[wi]);
					} else {
						if (col >= 0 && col < nc) {
							v.push_back(d[nc*row + col] * window[wi]);
						} else {
							v.push_back(fill * window[wi]);
						}
					}
				}
			}
			out[cell] = fun(v, narm);
		}
	}
}



void focal_win_sum(const std::vector<double> &d, std::vector<double> &out, int nc, int srow, int nr,
                    std::vector<double> window, int wnr, int wnc, double fill, bool narm, bool naonly, bool naomit, bool expand, bool global) {

	out.resize(nc*nr, NAN);
	int hwc = wnc / 2;
	int hwr = wnr / 2;
	bool nafill = std::isnan(fill);
	bool dofill = !(narm && nafill);
	std::vector<bool> winNA(window.size(), false);
	for (size_t i=0; i<window.size(); i++) {
		if (std::isnan(window[i])) winNA[i] = true; 
	}

	int wr1 = wnr - 1;
	int wc1 = wnc - 1;
	int nc1 = nc - 1;
	bool checkNA = naonly || naomit;

	for (int r=0; r < nr; r++) {
		int rread = r+srow;
		for (int c=0; c < nc; c++) {
			size_t cell = r*nc + c;
			if (checkNA) {
				size_t readcell = rread*nc + c;
				if (naonly) {
					if (!std::isnan(d[readcell])) {
						out[cell] = d[readcell];
						continue;
					}					
				} else if (std::isnan(d[readcell])) {
					continue;
				}
			}
			double value = 0;
			bool found = false;
			for (int rr=0; rr<wnr; rr++) {
				int offr = wr1 - rr;
				for (int cc=0; cc < wnc; cc++)  {
					int offc = wc1 - cc; 
					//int wi = wnc * offr + offc;
					int wi = wnc * rr + cc;
					if (winNA[wi]) {
						continue;
					}
					int row = rread + hwr - offr;
					int col = c + hwc - offc;
					if (global) {
						col = col < 0 ? nc + col : col;
						col = col > nc1 ? col - nc : col;
						if (narm) {
							if (!std::isnan(d[nc * row + col])) {
								value += d[nc*row + col] * window[wi];
								found = true;
							}
						} else {
							value += d[nc*row + col] * window[wi];
						}
					} else if (expand) {
						col = col < 0 ? 0 : col;
						col = col > nc1 ? nc1 : col;
						if (narm) {
							if (!std::isnan(d[nc * row + col])) {
								value += d[nc*row + col] * window[wi];
								found = true;
							}
						} else {
							value += d[nc*row + col] * window[wi];
						}
					} else {
						if (col >= 0 && col < nc) {
							if (narm) {
								if (!std::isnan(d[nc*row + col])) {
									value += d[nc*row + col] * window[wi];
									found = true;
								}
							} else {
								value += d[nc*row + col] * window[wi]; 
							}
						} else if (dofill) {
							value += fill * window[wi];
						}
					}
				}
			}
			if (narm) {
				if (found) {
					out[cell] = value;
				}
			} else {
				out[cell] = value;
			}
		}
	}
}

void focal_win_mean(const std::vector<double> &d, std::vector<double> &out, int nc, int srow, int nr,
                    std::vector<double> window, int wnr, int wnc, double fill, bool narm, bool naonly, bool naomit, bool expand, bool global) {

	out.resize(nc*nr, NAN);
	int hwc = wnc / 2;
	int hwr = wnr / 2;
	bool nafill = std::isnan(fill);
	bool dofill = !(narm && nafill);
	std::vector<bool> winNA(window.size(), false);
	double winsum = 0;
	std::vector<double> poswin = window;
	for (size_t i=0; i<window.size(); i++) {
		if (std::isnan(window[i])) {
			winNA[i] = true; 
		} else {
			if (window[i] < 0) {
				poswin[i] = -poswin[i];
			} 
			winsum += poswin[i]; 
		}
	}

	int wr1 = wnr - 1;
	int wc1 = wnc - 1;
	int nc1 = nc - 1;

	bool checkNA = naonly || naomit;

	for (int r=0; r<nr; r++) {
		int rread = r+srow;
		for (int c=0; c < nc; c++) {
			size_t cell = r*nc + c;
			if (checkNA) {
				size_t readcell = rread*nc + c;				
				if (naonly) {
					if (!std::isnan(d[readcell])) {
						out[cell] = d[readcell];
						continue;
					}					
				} else if (std::isnan(d[readcell])) {
					continue;
				}
			}
			double value = 0;
			if (narm) winsum = 0;
			for (int rr=0; rr<wnr; rr++) {
				int offr = wr1 - rr;
				for(int cc=0; cc < wnc; cc++)  {
					int offc = wc1 - cc; 
					//int wi = wnc * offr + offc;
					int wi = wnc * rr + cc;
					if (winNA[wi]) {
						continue;
					}
					int row = rread + hwr - offr;
					int col = c + hwc - offc;
					if (global) {
						col = col < 0 ? nc + col : col;
						col = col > nc1 ? col - nc : col;
						if (narm) {
							if (!std::isnan(d[nc * row + col])) {
								value += d[nc * row + col] * window[wi];
								winsum += poswin[wi];
							}
						} else {
							value += d[nc * row + col] * window[wi]; 
						}
					} else if (expand) {
						col = col < 0 ? 0 : col;
						col = col > nc1 ? nc1 : col;
						if (narm) {
							if (!std::isnan(d[nc * row + col])) {
								value += d[nc * row + col] * window[wi];
								winsum += poswin[wi];
							}
						} else {
							value += d[nc * row + col] * window[wi]; 
						}
					} else {
						if (col >= 0 && col < nc) {
							if (narm) {
								if (!std::isnan(d[nc * row + col])) {
									value += d[nc * row + col] * window[wi];
									winsum += poswin[wi];
								}
							} else {
								value += d[nc * row + col] * window[wi]; 
							}
						} else if (dofill) {
							value += fill;
							if (narm) {
								winsum += poswin[wi];
							}
						}
					}
				}
			}
			if (winsum > 0) {
				out[cell] = value / winsum;
			}
		}
	}
}




SpatRaster SpatRaster::focal(std::vector<unsigned> w, std::vector<double> m, double fillvalue, bool narm, bool naonly, bool naomit, std::string fun, bool expand, SpatOptions &opt) {

	SpatRaster out = geometry();
	bool global = is_global_lonlat();
	size_t nl = nlyr();

	if (!source[0].hasValues) { return(out); }

	if (w.size() != 2) {
		out.setError("size of w is not 1 or 2");
		return out;
	}
	if ((w[0] % 2) == 0 || (w[1] % 2) == 0) {
		out.setError("w must be odd sized");
		return out;
	}
	unsigned ww = w[0] * w[1];
	if (ww < 3) {
		out.setError("not a meanigful window");
		return out;
	}

	if (ww != m.size()) {
		out.setError("weights matrix size does not match prod(w)");
		return out;
	}

	size_t nc = ncol();
	size_t nr = nrow();
	if (w[0] > (nr*2)) {
		out.setError("nrow(w) > 2 * nrow(x)");
		return out;
	}
	if (w[1] > (nc*2)) { 
		out.setError("ncol(w) > 2 * ncol(x)");
		return out;
	}

	if (!readStart()) {
		out.setError(getError());
		return(out);
	}
	opt.ncopies += 2;
	opt.minrows = w[0] > nr ? nr : w[0];

 	if (!out.writeStart(opt)) {
		readStop();
		return out;
	}
	size_t hw0 = w[0]/2;
	size_t dhw0 = hw0 * 2;
	//size_t fsz = hw0*nc;
	size_t fsz2 = dhw0*nc;

	bool dofun = false;
	std::function<double(std::vector<double>&, bool)> fFun;
	if ((fun != "mean") && (fun != "sum")) {
		if (!haveFun(fun)) {
			out.setError("unknown function argument");
			return out;
		}
		fFun = getFun(fun);
		dofun = true;
	}

	std::vector<double> fill;
	for (size_t i = 0; i < out.bs.n; i++) {
		unsigned rstart, roff;
		unsigned rnrows = out.bs.nrows[i];
		if (i == 0) {
			rstart = 0;
			roff = dhw0;
			if (i != (out.bs.n-1)) {
				rnrows += hw0;
			}
		} else {
			rstart = out.bs.row[i] + hw0;
			roff = hw0;
			if (i == (out.bs.n-1)) {
				rnrows -= hw0; 
			}
		}

		std::vector<double> vout, voutcomb;
		std::vector<double> vin, vincomb;
		size_t off=0;
		if (nl > 1) {
			readValues(vincomb, rstart, rnrows, 0, nc);
			off = nc * out.bs.nrows[i];
			voutcomb.reserve(vincomb.size());
		} else {
			readValues(vin, rstart, rnrows, 0, nc);
		}
		for (size_t lyr=0; lyr<nl; lyr++) {
			vout.clear();
			if (nl > 1) {
				size_t lyroff = lyr * off;
				vin = {vincomb.begin() + lyroff, vincomb.begin() + lyroff + off}; 
			}
			if (i==0) {
				if (expand) {
					fill = {vin.begin(), vin.begin()+nc};
					for (size_t i=1; i<dhw0; i++) {
						fill.insert(fill.end(), fill.begin(), fill.end());
					}
				} else {
					fill.resize(fsz2, fillvalue);
				}
			}
			vin.insert(vin.begin(), fill.begin(), fill.end());

			if (i == (out.bs.n-1)) {
				if (expand) {
					fill = {vin.end()-nc, vin.end()};
					for (size_t i=1; i<dhw0; i++) {
						fill.insert(fill.end(), fill.begin(), fill.end());
					}
				} else {
					std::fill(fill.begin(), fill.end(), fillvalue);
				}
				vin.insert(vin.end(), fill.begin(), fill.end());
			}

			if (dofun) {
				focal_win_fun(vin, vout, nc, roff, out.bs.nrows[i], m, w[0], w[1], fillvalue, narm, naonly, naomit, expand, global, fFun);
			} else if (fun == "mean") {
				focal_win_mean(vin, vout, nc, roff, out.bs.nrows[i], m, w[0], w[1], fillvalue, narm, naonly, naomit, expand, global);
			} else {
				focal_win_sum(vin, vout, nc, roff, out.bs.nrows[i], m, w[0], w[1], fillvalue, narm, naonly, naomit, expand, global);
			}

			if (i != (out.bs.n-1)) {
				fill = {vin.end() - fsz2, vin.end() };  
			}

			if (naonly) {
				for (size_t j=0; j<vout.size(); j++) {
					size_t k = fsz2 + j;
					if (!std::isnan(vin[k])) {
						vout[j] = vin[k];
					}
				}
			}
			if (nl > 1) {
				voutcomb.insert(voutcomb.end(), vout.begin(), vout.end());
			}
		}
		if (nl > 1) {
			if (!out.writeBlock(voutcomb, i)) return out;
		} else {
			if (!out.writeBlock(vout, i)) return out;
		}
	}

	out.writeStop();
	readStop();
	return(out);
}




/*



SpatRaster SpatRaster::focal1(std::vector<unsigned> w, std::vector<double> m, double fillvalue, bool narm, bool naonly, std::string fun, SpatOptions &opt) {

	SpatRaster out = geometry(1);
	if (nlyr() > 1) {
		SpatOptions ops(opt);
		out.addWarning("focal computations are only done for the first input layer");
		std::vector<unsigned> lyr = {0};
		*this = subset(lyr, ops);
	}

	if (!source[0].hasValues) { return(out); }

	bool wmat = false;
	if (m.size() > 1) {
		wmat = true;
	} else if (w.size() == 1) {
		w.push_back(w[0]);
	}
	if (w.size() != 2) {
		out.setError("size of w is not 1 or 2");
		return out;
	}
	if ((w[0] % 2) == 0 || (w[1] % 2) == 0) {
		out.setError("w must be odd sized");
		return out;
	}
	if (w[0] < 3 && w[1] < 3) {
		out.setError("w must be > 1");
		return out;
	}

	unsigned ww = w[0] * w[1];
	if (ww < 9) {
		out.setError("not a meaningful window");
		return out;
	}
	if (wmat && (ww != m.size())) {
		out.setError("weight matrix error");
		return out;
	}
	if (!readStart()) {
		out.setError(getError());
		return(out);
	}
	opt.ncopies += 2 * ww;
 	if (!out.writeStart(opt)) {
		readStop();
		return out;
	}

	std::function<double(std::vector<double>&, bool)> fFun = getFun(fun);
	std::vector<double> v;
	if (naonly) {
		int mid = ww / 2;
		for (size_t i = 0; i < out.bs.n; i++) {
			std::vector<double> fv = focal_values(w, fillvalue, out.bs.row[i], out.bs.nrows[i]);
			v.resize(out.bs.nrows[i] * ncol(), NAN);
			if (wmat) {
				for (size_t j=0; j<v.size(); j++) {
					double midv = fv[j*ww+mid];
					if (std::isnan(midv)) {
						v[j] = 0;
						for (size_t k=0; k<ww; k++) {
							v[j] += fv[j*ww+k] * m[k];
						}
					} else {
						v[j] = midv;
					}
				}
			} else {
				for (size_t j=0; j<v.size(); j++) {
					unsigned off = j*ww;
					std::vector<double> x(fv.begin()+off, fv.begin()+off+ww);
					double midv = x[mid];
					if (std::isnan(midv)) {
						v[j] = fFun(x, narm);
					} else {
						v[j] = midv;
					}
				}
			}
			if (!out.writeBlock(v, i)) return out;
		}
	} else {
		for (size_t i = 0; i < out.bs.n; i++) {


			std::vector<double> fv = focal_values(w, fillvalue, out.bs.row[i], out.bs.nrows[i]);
			v.resize(out.bs.nrows[i] * ncol());
			if (wmat) {
				if (narm) {
					for (size_t j=0; j<v.size(); j++) {
						v[j] = 0;
						size_t cnt = 0;
						for (size_t k=0; k<ww; k++) {
							double vv = fv[j*ww+k] * m[k];
							if (!std::isnan(vv)) {
								v[j] += vv;
								cnt++;
							}
						}
						if (cnt == 0) v[j] = NAN;
					}
				} else {
					for (size_t j=0; j<v.size(); j++) {
						v[j] = 0;
						for (size_t k=0; k<ww; k++) {
							v[j] += fv[j*ww+k] * m[k];
						}
					}
				}
			} else {
				for (size_t j=0; j<v.size(); j++) {
					unsigned off = j*ww;
					std::vector<double> x(fv.begin()+off, fv.begin()+off+ww);
					v[j] = fFun(x, narm);
				}
			}
			if (!out.writeBlock(v, i)) return out;
		}
	}
	out.writeStop();
	readStop();
	return(out);
}


void focalrow(std::vector<double> &x, std::vector<double> n, 
            const size_t &w1, const size_t &ww, const size_t &nc, 
            const std::vector<double> &fill) {
	x.erase(x.begin(), x.begin() + w1);
	x.resize(x.size() + w1);
	n.insert(n.begin(), fill.begin(), fill.end());
	n.insert(n.end(), fill.begin(), fill.end());
	for (size_t i=0; i<nc; i++) {
		size_t off = (i+1)*ww - w1;
		for (size_t j=0; j<w1; j++) {
			x[off + j] = n[i + j];
		}
	}
}




SpatRaster SpatRaster::focal2(std::vector<unsigned> w, std::vector<double> m, double fillvalue, bool narm, bool naonly, std::string fun, SpatOptions &opt) {

	SpatRaster out = geometry(1);
	if (nlyr() > 1) {
		SpatOptions ops(opt);
		out.addWarning("focal computations are only done for the first input layer");
		std::vector<unsigned> lyr = {0};
		*this = subset(lyr, ops);
	}

	if (!source[0].hasValues) { return(out); }

	bool wmat = false;
	if (m.size() > 1) {
		wmat = true;
	} else if (w.size() == 1) {
		w.push_back(w[0]);
	}
	if (w.size() != 2) {
		out.setError("size of w is not 1 or 2");
		return out;
	}
	if ((w[0] % 2) == 0 || (w[1] % 2) == 0) {
		out.setError("w must be odd sized");
		return out;
	}
	if (w[0] < 3 && w[1] < 3) {
		out.setError("w must be > 1");
		return out;
	}
	size_t nc = ncol();
	size_t nr = nrow();
	if (w[0] > (nr*2)) {
		out.setError("nrow(w) > 2 * nrow(x)");
		return out;
//		w[0] = nr*2;
//		w[0] = std::fmod(w[0], 2) == 0 ? w[0]+1 : w[0];
	}
	if (w[1] > (nc*2)) { 
		out.setError("ncol(w) > 2 * ncol(x)");
		return out;
//		w[1] = nc*2;
//		w[1] = std::fmod(w[1], 2) == 0 ? w[1]+1 : w[1];
	}
	unsigned ww = w[0] * w[1];
	if (ww < 9) {
		out.setError("not a meaningful window");
		return out;
	}
	if (wmat && (ww != m.size())) {
		out.setError("weight matrix error");
		return out;
	}
	if (!readStart()) {
		out.setError(getError());
		return(out);
	}
	opt.ncopies += std::max(1, (int)(nc * ww / ncell()));
	opt.minrows = w[0] > nr ? nr : w[0];

 	if (!out.writeStart(opt)) {
		readStop();
		return out;
	}

	std::function<double(std::vector<double>&, bool)> fFun = getFun(fun);

	size_t hw0 = w[0]/2;
	std::vector<double> fill(hw0, fillvalue);

	std::vector<double> fvals(nc * ww, fillvalue);

	std::vector<double> vout;
	for (size_t i = 0; i < out.bs.n; i++) {
		if (i==0) {
			if (i == (out.bs.n-1)) {
				vout.resize(out.bs.nrows[i]* nc);
			} else {
				vout.resize((out.bs.nrows[i]-hw0) * nc);
			}
		} else if (i == (out.bs.n-1)) {
			vout.resize((out.bs.nrows[i] + 2 * hw0)* nc);
		} else if (i == 1) {
			vout.resize(out.bs.nrows[i]* nc);
		}

		std::vector<double> vin = readValues(out.bs.row[i], out.bs.nrows[i]);
		size_t start= 0;
		if (i == 0) {
			for (size_t r=0; r<hw0; r++) {
				unsigned off = r*nc;
				std::vector<double> row(vin.begin()+off, vin.begin()+off+nc);
				focalrow(fvals, row, w[1], ww, nc, fill);
				start++;
			}
		}

		for (size_t r=start; r< out.bs.nrows[i]; r++) {
			unsigned off1 = r*nc;
			std::vector<double> row(vin.begin()+off1, vin.begin()+off1+nc);
			focalrow(fvals, row, w[1], ww, nc, fill);
			unsigned off2 = (r-start)*nc;
			for (size_t j=0; j<nc; j++) {
				unsigned offset = j*ww;
				std::vector<double> x(fvals.begin()+offset, fvals.begin()+offset+ww);
				if (wmat) {
					vout[off2+j] = 0;
					for (size_t k=0; k<ww; k++) {
						vout[off2+j] += w[k] * x[k];
					}
				} else {
					vout[off2+j] = fFun(x, narm);
				}
			}
		}

		if (i == (out.bs.n - 1)) {
			size_t m = i == 0 ? 1 : 2;
			std::vector<double> row(nc, fillvalue);
			start = out.bs.nrows[i] - hw0;
			for (size_t r=0; r<(m*hw0); r++) {
				unsigned off = (start+r)*nc;
				focalrow(fvals, row, w[1], ww, nc, fill);
				for (size_t j=0; j<nc; j++) {
					unsigned offset = j*ww;
					std::vector<double> x(fvals.begin()+offset, fvals.begin()+offset+ww);
					if (wmat) {
						vout[off+j] = 0;
						for (size_t k=0; k<ww; k++) {
							vout[off+j] += w[k] * x[k];
						}
					} else {
						vout[off+j] = fFun(x, narm);
					}
				}
			}
		}

		int startrow = out.bs.row[i]; 
		int nrows = out.bs.nrows[i];
		if (i==0) {
			if (out.bs.n != 1) {
				nrows -= hw0;
			}
		} else if (i == (out.bs.n-1)) {
			startrow -= hw0;
			nrows += hw0;
		} else {
			startrow -= hw0;
		}
		if (!out.writeValues(vout, startrow, nrows, 0, nc)) return out;
	}
	out.writeStop();
	readStop();
	return(out);
}
*/