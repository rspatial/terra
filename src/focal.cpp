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

	std::vector<double> error;
	if (w.size() < 2) {
		setError("weights matrix must have more than one side");
		return(error);
	}

	for (size_t i=0; i<w.size(); i++) {
		if (w[i] % 2 == 0) {
			setError("weights matrix must have uneven sides");
			return(error);
		}
		if (w[i] < 1) {
			setError("weights matrix must have dimensions larger than 0");
			return(error);
		}
	}

	const bool global = is_global_lonlat();

	int_64 nr = nrow();
	nrows = std::min(nrows, nr - row + 1);

	int_64 nc = ncol();
	int_64 wr = w[0] / 2;
	int_64 wc = w[1] / 2;
	//may be unexpected
	//wr = std::min(wr, nr-1);
	//wc = std::min(wc, nc-1);

	int_64 startrow = row-wr;
	startrow = startrow < 0 ? 0 : startrow;
	int_64 startoff = row-startrow;

	nrows = nrows < 1 ? 1 : nrows;
	int_64 readnrows = nrows+startoff+wr;
	int_64 endoff = wr;
	if ((startrow+readnrows) > nr ) {
		readnrows = nr-startrow;
		endoff = readnrows - (nrows+startoff);
	}

// ??
	//wr = std::min(wr, std::max((int_64)1, nrows-1));

	size_t n = nrows * nc * w[0] * w[1];
	int_64 nrmax = nrows + startoff + endoff - 1;
	//int nrmax = d.size() / ncol - 1;
	size_t f = 0;

	std::vector<double> d;
	readValues(d, startrow, readnrows, 0, nc);
	std::vector<double> out(n, fillvalue);

// << "sr " << startrow << " so " << startoff << " rnr " << readnrows << " wr " << wr << " wc " << wc << " nrows " << nrows << std::endl;


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
							size_t idx = bcell+col;
							out[f] = d[idx];
						} else if (global) {
							if (col < 0) {
								col = nc + col;
							} else if (col >= nc) {
								col = col - nc;
							}
							size_t idx = bcell+col;
							out[f] = d[idx];
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
//	opt.ncopies += 2;
	opt.minrows = w[0] > nr ? nr : w[0];

 	if (!out.writeStart(opt, filenames())) {
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
	if (nl == 1) {
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
			std::vector<double> vout, vin;
			readValues(vin, rstart, rnrows, 0, nc);
			vout.clear();

			if (i==0) {
				if (expand) {
					fill.reserve(dhw0 * nc);
					for (size_t i=0; i<dhw0; i++) {
						fill.insert(fill.end(), vin.begin(), vin.begin()+nc);
					}
				} else {
					fill.resize(fsz2, fillvalue);
				}
			}

			vin.insert(vin.begin(), fill.begin(), fill.end());

			if (i == (out.bs.n-1)) {
				if (expand) {
					fill.resize(0);
					fill.reserve(dhw0 * nc);
					for (size_t i=1; i<dhw0; i++) {
						fill.insert(fill.end(), vin.end()-nc, vin.end());
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
			if (!out.writeBlock(vout, i)) return out;
		}
	} else {
		for (size_t i = 0; i < out.bs.n; i++) {
			unsigned rstart, roff;
			unsigned rnrows = out.bs.nrows[i];
			if (i != (out.bs.n-1)) {
				rnrows += hw0;
			}
			if (i == 0) {
				rstart = 0;
				roff = dhw0;
			} else {
				rstart = out.bs.row[i] - hw0;
				roff = hw0;
				rnrows += hw0;
			}

			std::vector<double> vout, voutcomb, vin, vincomb;
			size_t off=0;
			readValues(vincomb, rstart, rnrows, 0, nc);
			off = nc * rnrows; //out.bs.nrows[i];
			voutcomb.reserve(vincomb.size());
			for (size_t lyr=0; lyr<nl; lyr++) {
				vout.clear();
				size_t lyroff = lyr * off;
				vin = {vincomb.begin() + lyroff, vincomb.begin() + lyroff + off};
				if (i==0) {
					if (expand) {
						fill.resize(0);
						fill.reserve(dhw0 * nc);
						for (size_t i=0; i<dhw0; i++) {
							fill.insert(fill.end(), vin.begin(), vin.begin()+nc);
						}
					} else {
						fill.resize(fsz2, fillvalue);
					}
					vin.insert(vin.begin(), fill.begin(), fill.end());
				}

				if (i == (out.bs.n-1)) {
					if (expand) {
						fill.resize(0);
						fill.reserve(dhw0 * nc);
						for (size_t i=1; i<dhw0; i++) {
							fill.insert(fill.end(), vin.end()-nc, vin.end());
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
				voutcomb.insert(voutcomb.end(), vout.begin(), vout.end());
			}
			if (!out.writeBlock(voutcomb, i)) return out;
		}
	}
	out.writeStop();
	readStop();
	return(out);
}

