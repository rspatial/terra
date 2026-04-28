// Copyright (c) 2018-2026  Robert J. Hijmans
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
#include "sort.h"

#include <algorithm>
#include <cmath>
#include <cstring>
#include <vector>

#include "tbb_helper.h"


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




std::vector<double> SpatRaster::focal_values(std::vector<unsigned> w, double fillvalue, int64_t row, int64_t nrows, SpatOptions &ops) {

	if (nlyr() > 1) {
		std::vector<size_t> lyr = {0};
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

	int64_t nr = nrow();
	nrows = std::min(nrows, nr - row + 1);

	int64_t nc = ncol();
	int64_t wr = w[0] / 2;
	int64_t wc = w[1] / 2;
	//may be unexpected
	//wr = std::min(wr, nr-1);
	//wc = std::min(wc, nc-1);

	int64_t startrow = row-wr;
	startrow = startrow < 0 ? 0 : startrow;
	int64_t startoff = row-startrow;

	nrows = nrows < 1 ? 1 : nrows;
	int64_t readnrows = nrows+startoff+wr;
	int64_t endoff = wr;
	if ((startrow+readnrows) > nr ) {
		readnrows = nr-startrow;
		endoff = readnrows - (nrows+startoff);
	}

// ??
	//wr = std::min(wr, std::max((int64_t)1, nrows-1));

	size_t n = nrows * nc * w[0] * w[1];
	int64_t nrmax = nrows + startoff + endoff - 1;

	std::vector<double> d;
	readValues(d, startrow, readnrows, 0, nc);
	std::vector<double> out(n, fillvalue);

	const size_t per_row_out = (size_t)nc * w[0] * w[1];
	const size_t per_cell_out = (size_t)w[0] * w[1];
	const int64_t wcols = (int64_t)w[1];

	// Per-row kernel. f is computed from (r, c) so disjoint rows write to
	// disjoint output slices and the loop is safe to parallelise.
	//
	// Two fast paths inside the column loop:
	//   * interior columns (c-wc >= 0 && c+wc < nc): the entire window row
	//     fits in the source row, so we can std::copy wnc doubles in one go
	//     and skip per-column branches.
	//   * non-global edge columns: the output is already preinitialised to
	//     fillvalue, so we only need to write the in-bounds slice.
	auto kernel = [&](int64_t r_begin, int64_t r_end) {
		for (int64_t r = r_begin; r < r_end; r++) {
			for (int64_t c = 0; c < nc; c++) {
				size_t f = (size_t)r * per_row_out + (size_t)c * per_cell_out;
				const bool col_interior = (c - wc >= 0) && (c + wc < nc);
				for (int64_t i = -wr; i <= wr; i++) {
					int64_t row = r + startoff + i;
					if ((row < 0) || (row > nrmax)) {
						// out-of-range row: fillvalue is already in `out`
						f += wcols;
						continue;
					}
					size_t bcell = (size_t)row * (size_t)nc;
					if (col_interior) {
						std::copy(d.begin() + bcell + c - wc,
						          d.begin() + bcell + c + wc + 1,
						          out.begin() + f);
						f += wcols;
					} else if (global) {
						for (int64_t j = -wc; j <= wc; j++) {
							int64_t col = c + j;
							if (col < 0)        col += nc;
							else if (col >= nc) col -= nc;
							out[f++] = d[bcell + col];
						}
					} else {
						// non-global edge: write only the in-bounds slice.
						// Out-of-bounds cells keep the preinitialised fill.
						int64_t j_lo = std::max<int64_t>(-wc, -c);
						int64_t j_hi = std::min<int64_t>(wc, nc - 1 - c);
						if (j_lo <= j_hi) {
							size_t off_in   = bcell + c + j_lo;
							size_t off_out  = f + (size_t)(j_lo + wc);
							std::copy(d.begin() + off_in,
							          d.begin() + off_in + (size_t)(j_hi - j_lo + 1),
							          out.begin() + off_out);
						}
						f += wcols;
					}
				}
			}
		}
	};

#if defined(USE_TBB)
	// Parallelise the row loop. Threshold avoids task overhead on tiny
	// blocks (e.g. when the R-side `focal()` calls focalValues on a single
	// row to probe the user's function). The grain size is in source rows;
	// we pick something that keeps each task ~O(64KB-1MB) of output.
	if (ops.parallel && nrows > 4 && per_row_out >= 4096) {
		int64_t grain = std::max<int64_t>(1, 65536 / (int64_t)per_row_out);
		terra_parallel_for(ops, tbb::blocked_range<int64_t>(0, nrows, grain),
			[&](const tbb::blocked_range<int64_t> &r) {
				kernel(r.begin(), r.end());
			});
	} else {
		kernel(0, nrows);
	}
#else
	kernel(0, nrows);
#endif

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



// ── old focal, see focal2 ────────────────────────────────────────

SpatRaster SpatRaster::focal(std::vector<unsigned> w, std::vector<double> m, double fillvalue, bool narm, bool naonly, bool naomit, std::string fun, bool expand, SpatOptions &opt) {

	std::vector<std::string> f {"modal", "min", "max", "first"};
	auto it = std::find(f.begin(), f.end(), fun);
	bool props = (it != f.end());

	SpatRaster out = geometry(nlyr(), props);
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

// ===  focal2 ========================================

//  faster than focal for two common cases:
//
//   1. PLAIN UNIFORM WINDOW with constant weights 
//      per-pixel sliding accumulator (separable column / row pass).
//
//   2. RANK-1 WEIGHTED WINDOW (m == col_kernel ⊗ row_kernel, no NA in m)
//      two-pass convolution.


// ── window-shape detection ---------------------------------------------------

// All weights identical (no NaN). Constant weight stored in *w_const*.
static bool window_is_plain(const std::vector<double> &m, double &w_const) {
	if (m.empty()) return false;
	double w0 = m[0];
	if (std::isnan(w0)) return false;
	for (size_t i = 1; i < m.size(); i++) {
		if (std::isnan(m[i]) || m[i] != w0) return false;
	}
	w_const = w0;
	return true;
}

// Rank-1 outer-product check: m[r*wnc + c] == col_k[r] * row_k[c] for all r,c.
// Uses a tiny tolerance so that floating-point input (e.g. Gaussian kernels)
// is recognised. Fills row_k (size wnc) and col_k (size wnr) on success.
static bool window_is_rank1(const std::vector<double> &m, int wnr, int wnc,
                            std::vector<double> &row_k, std::vector<double> &col_k) {
	for (double x : m) if (std::isnan(x)) return false;

	int piv_r = -1, piv_c = -1;
	double piv_v = 0.0;
	double piv_a = 0.0;
	for (int r = 0; r < wnr; r++) {
		for (int c = 0; c < wnc; c++) {
			double a = std::abs(m[r * wnc + c]);
			if (a > piv_a) { piv_a = a; piv_r = r; piv_c = c; piv_v = m[r * wnc + c]; }
		}
	}
	if (piv_r < 0) {                          // all zero
		row_k.assign(wnc, 0.0);
		col_k.assign(wnr, 0.0);
		return true;
	}

	row_k.assign(wnc, 0.0);
	col_k.assign(wnr, 0.0);
	for (int c = 0; c < wnc; c++) row_k[c] = m[piv_r * wnc + c] / piv_v;
	for (int r = 0; r < wnr; r++) col_k[r] = m[r * wnc + piv_c];

	double scale = std::max(1.0, piv_a);
	double tol = 1e-10 * scale;
	for (int r = 0; r < wnr; r++) {
		for (int c = 0; c < wnc; c++) {
			double e = col_k[r] * row_k[c];
			if (std::abs(m[r * wnc + c] - e) > tol) return false;
		}
	}
	return true;
}


// ── helper: clamp/wrap a column index given the boundary mode ----------------
//
// returns a (column, fill_flag) pair. fill_flag == true means the caller
// must use the fillvalue contribution rather than reading d[col].

static inline int boundary_col(int col, int nc, bool expand, bool global) {
	if (col >= 0 && col < nc) return col;
	if (global) {
		if (col < 0)   return col + nc;
		else           return col - nc;
	}
	if (expand) {
		if (col < 0)   return 0;
		else           return nc - 1;
	}
	return -1;                                // fillvalue territory
}


// ── PLAIN UNIFORM SLIDING BOX  (sum & mean) ---------------------------------
//
// Maintains, for the current output row r:
//   colsum[c] = Σ_{i=-hwr..hwr}  contribution(d[(r+i)*nc + c])
//   colcnt[c] = Σ_{i=-hwr..hwr}  is_valid(d[(r+i)*nc + c])
//
// "contribution" is the value if non-NaN, 0 if NaN. is_valid is the
// indicator of non-NaN.
//
// Within a row a horizontal running window of width wnc is then maintained
// in O(1) per output column. Over the whole raster the total work is
// O(nr * nc) — independent of the window size — so for a 21x21 mean over
// a 1000x1000 raster this is ~441x faster than the naive O(nr*nc*W²) loop.
//
// Boundary modes:
//   * fillvalue: a virtual column made of `fillvalue` rows is used outside.
//                Its colsum / colcnt values are precomputed once.
//   * expand   : column 0 (or nc-1) is replicated. We just clamp the col.
//   * global   : column index wraps. We just wrap the col.
//
// emit() is templated so the same machinery serves both "sum" and "mean".

// Process a contiguous output-row range [r_start, r_end). The function is
// self-contained: it builds its own vertical column accumulator from vin,
// so disjoint row ranges can run on different threads concurrently.

template <typename Emit>
static void box_run_range(const std::vector<double> &d, std::vector<double> &out,
                          int nc, int srow, int r_start, int r_end,
                          int wnr, int wnc, double w_const,
                          double fill, bool narm, bool naonly, bool naomit,
                          bool expand, bool global, Emit emit) {

	const int hwr = wnr / 2;
	const int hwc = wnc / 2;

	const bool checkNA = naonly || naomit;
	const bool fill_is_nan = std::isnan(fill);

	// colsum[c]/colcnt[c] cover rows [srow+r_start-hwr .. srow+r_start+hwr]
	std::vector<double> colsum(nc, 0.0);
	std::vector<int>    colcnt(nc, 0);
	for (int i = -hwr; i <= hwr; i++) {
		int rr = srow + r_start + i;
		const double *row = d.data() + (size_t) rr * nc;
		for (int c = 0; c < nc; c++) {
			double v = row[c];
			if (!std::isnan(v)) { colsum[c] += v; colcnt[c]++; }
		}
	}

	// precomputed colsum/colcnt for an out-of-bounds (fillvalue) column.
	double fill_colsum;
	int    fill_colcnt;
	if (fill_is_nan) {
		fill_colsum = 0.0;
		fill_colcnt = 0;
	} else {
		fill_colsum = (double) wnr * fill;
		fill_colcnt = wnr;
	}

	auto get_colsum = [&](int c) -> double {
		if (c >= 0 && c < nc) return colsum[c];
		if (global)            return colsum[(c + nc) % nc];
		if (expand)            return colsum[c < 0 ? 0 : nc - 1];
		return fill_colsum;
	};
	auto get_colcnt = [&](int c) -> int {
		if (c >= 0 && c < nc) return colcnt[c];
		if (global)            return colcnt[(c + nc) % nc];
		if (expand)            return colcnt[c < 0 ? 0 : nc - 1];
		return fill_colcnt;
	};

	const int wsize = wnr * wnc;

	for (int r = r_start; r < r_end; r++) {
		double winsum = 0.0;
		int    wincnt = 0;
		for (int c = -hwc; c <= hwc; c++) {
			winsum += get_colsum(c);
			wincnt += get_colcnt(c);
		}

		for (int c = 0; c < nc; c++) {
			size_t cell = (size_t) r * nc + c;
			if (checkNA) {
				double iv = d[(size_t) (srow + r) * nc + c];
				if (naonly) {
					if (!std::isnan(iv)) {
						out[cell] = iv;
						goto slide_h;
					}
				} else {
					if (std::isnan(iv)) {
						out[cell] = iv;
						goto slide_h;
					}
				}
			}

			out[cell] = emit(winsum, wincnt, wsize, w_const, narm, fill_is_nan);

slide_h:
			if (c < nc - 1) {
				int cdrop = c - hwc;
				int cadd  = c + hwc + 1;
				winsum += get_colsum(cadd) - get_colsum(cdrop);
				wincnt += get_colcnt(cadd) - get_colcnt(cdrop);
			}
		}

		if (r < r_end - 1) {
			int rdrop = srow + r - hwr;
			int radd  = srow + r + hwr + 1;
			const double *drop_row = d.data() + (size_t) rdrop * nc;
			const double *add_row  = d.data() + (size_t) radd  * nc;
			for (int c = 0; c < nc; c++) {
				double dv = drop_row[c];
				double av = add_row[c];
				if (!std::isnan(dv)) { colsum[c] -= dv; colcnt[c]--; }
				if (!std::isnan(av)) { colsum[c] += av; colcnt[c]++; }
			}
		}
	}
}


// Driver: dispatch the row range in parallel chunks when TBB is available.
// Each chunk pays an O(nc * wnr) setup cost to build its own column
// accumulator, so we only parallelize when each chunk has enough rows to
// amortize that cost.

template <typename Emit>
static void box_run(const std::vector<double> &d, std::vector<double> &out,
                    int nc, int srow, int nr,
                    int wnr, int wnc, double w_const,
                    double fill, bool narm, bool naonly, bool naomit,
                    bool expand, bool global, const SpatOptions &opt, Emit emit) {
#if defined(USE_TBB)
	if (opt.parallel && nr > 1 && (long long) nr * nc > 50000) {
		// each chunk should be at least ~max(wnr, 16) rows to amortize
		// the O(nc * wnr) accumulator-setup overhead.
		size_t grain = (size_t) std::max(wnr, 16);
		terra_parallel_for(opt, tbb::blocked_range<int>(0, nr, grain),
			[&](const tbb::blocked_range<int> &range) {
				box_run_range(d, out, nc, srow, range.begin(), range.end(),
				              wnr, wnc, w_const, fill,
				              narm, naonly, naomit, expand, global, emit);
			});
		return;
	}
#else
	(void) opt;
#endif
	box_run_range(d, out, nc, srow, 0, nr, wnr, wnc, w_const, fill,
	              narm, naonly, naomit, expand, global, emit);
}


// emitters --------------------------------------------------------------------

static inline double emit_sum(double s, int cnt, int wsize, double w, bool narm, bool fill_is_nan) {
	if (narm) {
		// fillvalue NA + narm → no contribution. winsum already excludes it.
		return cnt > 0 ? w * s : NAN;
	}
	// not narm: any NA inside the window → NA.
	// Number of "missing" cells = wsize - cnt minus the fillvalue NA cells
	// which are already excluded. With fillvalue NA, those cells should
	// also propagate NA → so condition simplifies to: cnt < wsize.
	// With finite fillvalue, the column accumulator already counted the
	// fill cells as "valid" (we set fill_colcnt = wnr).
	if (fill_is_nan) {
		return cnt == wsize ? w * s : NAN;
	}
	return cnt == wsize ? w * s : NAN;
}

static inline double emit_mean(double s, int cnt, int wsize, double w,
                               bool narm, bool fill_is_nan) {
	(void) wsize; (void) w;
	if (cnt <= 0) return NAN;
	if (!narm && fill_is_nan && cnt < wsize) return NAN;
	if (!narm && cnt < wsize) return NAN;
	return s / cnt;
}


// ── RANK-1 SEPARABLE CONVOLUTION  (sum & mean) ------------------------------
//
// Two passes:
//   pass 1: convolve every input row by row_k -> tmp     (size: nrows_in × nc)
//   pass 2: convolve every column of tmp by col_k -> out (size: nr × nc)
//
// Each pass touches O(W) per pixel: O((wnc + wnr) * nr * nc) total.
//
// The "valid" mask is propagated alongside the value (NaN-skipping), and
// the final emitter applies the sum/mean rule. The denominator for the
// mean is Σ |effective_weight| over non-NaN cells (matches focal_win_mean).

static void separable_run(const std::vector<double> &d, std::vector<double> &out,
                          int nc, int srow, int nr, int wnr, int wnc,
                          const std::vector<double> &row_k, const std::vector<double> &col_k,
                          double fill, bool narm, bool naonly, bool naomit,
                          bool expand, bool global, bool is_mean, const SpatOptions &opt) {

	const int hwr = wnr / 2;
	const int hwc = wnc / 2;
	const int rows_in = nr + 2 * hwr;
	const bool checkNA = naonly || naomit;
	const bool fill_is_nan = std::isnan(fill);

	std::vector<double> row_kabs(wnc), col_kabs(wnr);
	double row_kabs_sum = 0.0, col_kabs_sum = 0.0;
	for (int c = 0; c < wnc; c++) { row_kabs[c] = std::abs(row_k[c]); row_kabs_sum += row_kabs[c]; }
	for (int r = 0; r < wnr; r++) { col_kabs[r] = std::abs(col_k[r]); col_kabs_sum += col_kabs[r]; }
	const double full_kabs_sum = row_kabs_sum * col_kabs_sum;

	std::vector<double> tmp_val((size_t) rows_in * nc);
	std::vector<double> tmp_w  ((size_t) rows_in * nc);

	// Pass 1: 1-D row convolution with row_k → tmp. Embarrassingly parallel
	// over input rows (each output row reads only its own input row).
	auto pass1_range = [&](int rr_begin, int rr_end) {
		for (int rr = rr_begin; rr < rr_end; rr++) {
			const double *src = d.data() + (size_t) (srow - hwr + rr) * nc;
			double *dv = tmp_val.data() + (size_t) rr * nc;
			double *dw = tmp_w  .data() + (size_t) rr * nc;
			for (int c = 0; c < nc; c++) {
				double s = 0.0;
				double wsum = 0.0;
				for (int k = 0; k < wnc; k++) {
					int sc = c + k - hwc;
					int eff = boundary_col(sc, nc, expand, global);
					double v = (eff < 0) ? fill : src[eff];
					if (std::isnan(v)) {
						if (!narm) s = NAN;
					} else {
						if (!std::isnan(s)) s += row_k[k] * v;
						wsum += row_kabs[k];
					}
				}
				dv[c] = std::isnan(s) ? NAN : s;
				dw[c] = wsum;
			}
		}
	};

	// Pass 2: 1-D column convolution with col_k → out. Embarrassingly
	// parallel over output rows (only reads tmp, only writes its own row).
	auto pass2_range = [&](int r_begin, int r_end) {
		for (int r = r_begin; r < r_end; r++) {
			double *outrow = out.data() + (size_t) r * nc;
			const double *iv_row = checkNA ? d.data() + (size_t) (srow + r) * nc : nullptr;

			for (int c = 0; c < nc; c++) {
				if (checkNA) {
					double iv = iv_row[c];
					if (naonly) {
						if (!std::isnan(iv)) { outrow[c] = iv; continue; }
					} else {
						if (std::isnan(iv))  { outrow[c] = iv; continue; }
					}
				}
				double s = 0.0;
				double wsum = 0.0;
				for (int k = 0; k < wnr; k++) {
					size_t srcrow = (size_t) (r + k);
					double v  = tmp_val[srcrow * nc + c];
					double wv = tmp_w  [srcrow * nc + c];
					if (std::isnan(v)) {
						if (!narm) s = NAN;
					} else {
						if (!std::isnan(s)) s += col_k[k] * v;
						wsum += col_kabs[k] * wv;
					}
				}
				if (std::isnan(s)) {
					outrow[c] = NAN;
				} else if (is_mean) {
					outrow[c] = (wsum > 0.0) ? (s / wsum) : NAN;
				} else {
					if (!narm && (wsum < full_kabs_sum - 1e-12 * full_kabs_sum) && fill_is_nan) {
						outrow[c] = NAN;
					} else {
						outrow[c] = s;
					}
				}
			}
		}
	};

#if defined(USE_TBB)
	if (opt.parallel && (long long) rows_in * nc > 20000) {
		terra_parallel_for(opt, tbb::blocked_range<int>(0, rows_in, 8),
			[&](const tbb::blocked_range<int> &r) { pass1_range(r.begin(), r.end()); });
		terra_parallel_for(opt, tbb::blocked_range<int>(0, nr, 8),
			[&](const tbb::blocked_range<int> &r) { pass2_range(r.begin(), r.end()); });
		return;
	}
#else
	(void) opt;
#endif
	pass1_range(0, rows_in);
	pass2_range(0, nr);
}


// ── focal2 ────────────────────────────────────────

SpatRaster SpatRaster::focal2(std::vector<unsigned> w, std::vector<double> m, double fillvalue, bool narm, bool naonly, bool naomit, std::string fun, bool expand, SpatOptions &opt) {

	std::vector<std::string> props_funs {"modal", "min", "max", "first"};
	bool props = std::find(props_funs.begin(), props_funs.end(), fun) != props_funs.end();

	SpatRaster out = geometry(nlyr(), props);
	bool global = is_global_lonlat();
	size_t nl = nlyr();

	if (!source[0].hasValues) { return out; }

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
	if (w[0] > (nr * 2)) { out.setError("nrow(w) > 2 * nrow(x)"); return out; }
	if (w[1] > (nc * 2)) { out.setError("ncol(w) > 2 * ncol(x)"); return out; }

	if (!readStart()) { out.setError(getError()); return out; }

	opt.minrows = w[0] > nr ? nr : w[0];
	if (!out.writeStart(opt, filenames())) { readStop(); return out; }

	const size_t hw0 = w[0] / 2;
	const size_t dhw0 = hw0 * 2;
	const size_t fsz2 = dhw0 * nc;

	// dispatch decisions (precomputed once)
	const bool is_sum  = (fun == "sum");
	const bool is_mean = (fun == "mean");
	bool dofun = false;
	std::function<double(std::vector<double>&, bool)> fFun;
	if (!is_sum && !is_mean) {
		if (!haveFun(fun)) {
			out.setError("unknown function argument");
			readStop();
			return out;
		}
		fFun = getFun(fun);
		dofun = true;
	}

	double w_const = 0.0;
	bool plain = window_is_plain(m, w_const);

	std::vector<double> row_k, col_k;
	bool rank1 = false;
	if (!plain && (is_sum || is_mean)) {
		rank1 = window_is_rank1(m, (int) w[0], (int) w[1], row_k, col_k);
	}

	// per-block driver. opt is captured by reference so opt.parallel and
	// opt.threads (the global TBB thread cap) are honored consistently.
	auto run_block = [&](std::vector<double> &vin, std::vector<double> &vout, size_t roff, size_t bnr) {
		if (plain && (is_sum || is_mean)) {
			vout.assign(nc * bnr, NAN);
			if (is_sum) {
				box_run(vin, vout, (int) nc, (int) roff, (int) bnr,
				        (int) w[0], (int) w[1], w_const,
				        fillvalue, narm, naonly, naomit, expand, global, opt, emit_sum);
			} else {
				box_run(vin, vout, (int) nc, (int) roff, (int) bnr,
				        (int) w[0], (int) w[1], w_const,
				        fillvalue, narm, naonly, naomit, expand, global, opt, emit_mean);
			}
		} else if (rank1) {
			vout.assign(nc * bnr, NAN);
			separable_run(vin, vout, (int) nc, (int) roff, (int) bnr,
			              (int) w[0], (int) w[1], row_k, col_k,
			              fillvalue, narm, naonly, naomit, expand, global, is_mean, opt);
		} else if (dofun) {
			focal_win_fun(vin, vout, (int) nc, (int) roff, (int) bnr, m,
			              (int) w[0], (int) w[1], fillvalue,
			              narm, naonly, naomit, expand, global, fFun);
		} else if (is_mean) {
			focal_win_mean(vin, vout, (int) nc, (int) roff, (int) bnr, m,
			               (int) w[0], (int) w[1], fillvalue,
			               narm, naonly, naomit, expand, global);
		} else {
			focal_win_sum(vin, vout, (int) nc, (int) roff, (int) bnr, m,
			              (int) w[0], (int) w[1], fillvalue,
			              narm, naonly, naomit, expand, global);
		}
	};

	std::vector<double> fill;
	if (nl == 1) {
		for (size_t i = 0; i < out.bs.n; i++) {
			unsigned rstart, roff;
			unsigned rnrows = out.bs.nrows[i];
			if (i == 0) {
				rstart = 0;
				roff = dhw0;
				if (i != (out.bs.n - 1)) rnrows += hw0;
			} else {
				rstart = out.bs.row[i] + hw0;
				roff = hw0;
				if (i == (out.bs.n - 1)) rnrows -= hw0;
			}
			std::vector<double> vout, vin;
			readValues(vin, rstart, rnrows, 0, nc);
			vout.clear();

			if (i == 0) {
				if (expand) {
					fill.reserve(dhw0 * nc);
					for (size_t k = 0; k < dhw0; k++) {
						fill.insert(fill.end(), vin.begin(), vin.begin() + nc);
					}
				} else {
					fill.resize(fsz2, fillvalue);
				}
			}
			vin.insert(vin.begin(), fill.begin(), fill.end());

			if (i == (out.bs.n - 1)) {
				if (expand) {
					fill.resize(0);
					fill.reserve(dhw0 * nc);
					for (size_t k = 1; k < dhw0; k++) {
						fill.insert(fill.end(), vin.end() - nc, vin.end());
					}
				} else {
					std::fill(fill.begin(), fill.end(), fillvalue);
				}
				vin.insert(vin.end(), fill.begin(), fill.end());
			}

			run_block(vin, vout, roff, out.bs.nrows[i]);

			if (i != (out.bs.n - 1)) {
				fill = {vin.end() - fsz2, vin.end()};
			}
			if (!out.writeBlock(vout, i)) { readStop(); return out; }
		}
	} else {
		for (size_t i = 0; i < out.bs.n; i++) {
			unsigned rstart, roff;
			unsigned rnrows = out.bs.nrows[i];
			if (i != (out.bs.n - 1)) rnrows += hw0;
			if (i == 0) {
				rstart = 0;
				roff = dhw0;
			} else {
				rstart = out.bs.row[i] - hw0;
				roff = hw0;
				rnrows += hw0;
			}

			std::vector<double> vout, voutcomb, vin, vincomb;
			size_t off = 0;
			readValues(vincomb, rstart, rnrows, 0, nc);
			off = nc * rnrows;
			voutcomb.reserve(vincomb.size());
			for (size_t lyr = 0; lyr < nl; lyr++) {
				vout.clear();
				size_t lyroff = lyr * off;
				vin = {vincomb.begin() + lyroff, vincomb.begin() + lyroff + off};
				if (i == 0) {
					if (expand) {
						fill.resize(0);
						fill.reserve(dhw0 * nc);
						for (size_t k = 0; k < dhw0; k++) {
							fill.insert(fill.end(), vin.begin(), vin.begin() + nc);
						}
					} else {
						fill.resize(fsz2, fillvalue);
					}
					vin.insert(vin.begin(), fill.begin(), fill.end());
				}
				if (i == (out.bs.n - 1)) {
					if (expand) {
						fill.resize(0);
						fill.reserve(dhw0 * nc);
						for (size_t k = 1; k < dhw0; k++) {
							fill.insert(fill.end(), vin.end() - nc, vin.end());
						}
					} else {
						std::fill(fill.begin(), fill.end(), fillvalue);
					}
					vin.insert(vin.end(), fill.begin(), fill.end());
				}

				run_block(vin, vout, roff, out.bs.nrows[i]);
				voutcomb.insert(voutcomb.end(), vout.begin(), vout.end());
			}
			if (!out.writeBlock(voutcomb, i)) { readStop(); return out; }
		}
	}
	out.writeStop();
	readStop();
	return out;
}


















// ===  moving window ========================================


inline double rarea(const double &Ax, const double &Ay, const double &Bx, const double &By, const double &Cx, const double &Cy) {
   return std::abs( (Bx*Ay - Ax*By) + (Cx*By - Bx*Cy) + (Ax*Cy - Cx*Ay) ) / 2;
}


void sortvecs(std::vector<double> &X, std::vector<double> &Y, std::vector<double> &Z) {  
	std::vector<std::size_t> p = sort_order_a(X);
	permute(X, p);
	permute(Y, p);
	permute(Z, p);
  
	p = sort_order_a(Y);
	permute(X, p);
	permute(Y, p);
	permute(Z, p);
}


std::vector<std::vector<double>> SpatRaster::win_rect(std::vector<double> x, std::vector<double> y, std::vector<double> z, std::vector<double> win, SpatOptions &opt) {

	sortvecs(x, y, z);

	win[0] = std::abs(win[0]);
	win[1] = std::abs(win[1]);
    const double h = win[0] / 2; 
    const double w = win[1] / 2;

	// multiply for floating point imprecision
	const double rar = win[0] * win[1] * 1.00000001;
	double angle = std::fmod(win[2], 360.0);
	if (angle < 0) angle += 360.0;
    const bool rotated = angle != 0.0;

	double cphi=0, sphi=0, bigw=0, bigh=0, offh;
	double wcphi=0, hcphi=0, wsphi=0, hsphi=0;
	std::vector<double> ox(4);
	std::vector<double> oy(4);

	if (rotated) {
		angle = angle * M_PI / 180.0;
		cphi = cos(angle);
		sphi = sin(angle);
		wcphi= cphi * w;
		hcphi= cphi * h;
		wsphi= sphi * w;
		hsphi= sphi * h;

		ox[0] = -wcphi - hsphi;
		oy[0] = -wsphi + hcphi;
		ox[1] =  wcphi - hsphi;
		oy[1] =  wsphi + hcphi;
		ox[2] =  wcphi + hsphi;
		oy[2] =  wsphi - hcphi;
		ox[3] = -wcphi + hsphi;
		oy[3] = -wsphi - hcphi;

		bigw = (vmax(ox, false) - vmin(ox, false))/2;
		bigh = (vmax(oy, false) - vmin(oy, false))/2;
		offh = bigh * 1.00000001;
	} else {
		offh = h * 1.00000001;
	}

	const size_t nc = ncol();
	const size_t nr = nrow();

	size_t np = x.size() * 2;
	std::vector<std::vector<double>> out(2);
	out[0].reserve(np);
	out[1].reserve(np);

	size_t minpt = win[3] < 2 ? 1 : win[3];
	std::vector<double> rx(4);
	std::vector<double> ry(4);

	std::vector<int64_t> cols(nc);
	std::iota(cols.begin(), cols.end(), 0);
	std::vector<double> xc = xFromCol(cols);

	if (minpt < 2) {

		if (rotated) {
			for (size_t r=0; r<nr; r++) {
				double yrow = yFromRow(r);
				double ytop = yrow + offh;
				double ybot = yrow - offh;
				double rnc = r * nc;
				ry[0] = yrow + oy[0];
				ry[1] = yrow + oy[1];
				ry[2] = yrow + oy[2];
				ry[3] = yrow + oy[3];
				np = y.size();
				for (long i=(np-1); i>=0; i--) {
					if (y[i] > ytop) {
						y.pop_back(); // above current row
					} else if (y[i] >= ybot) {
						for (long j=(nc-1); j>=0; j--) {
							double dist = x[i] - xc[j];
							if (dist > bigw) {
								break;
							} else if (std::abs(dist) <= bigw) {
								rx[0] = xc[j] + ox[0];
								rx[1] = xc[j] + ox[1];
								rx[2] = xc[j] + ox[2];
								rx[3] = xc[j] + ox[3];
								// triangles apd, dpc, cpb, bpa
								double area = rarea(rx[0], ry[0], x[i], y[i], rx[3], ry[3]);
								area += rarea(rx[3], ry[3], x[i], y[i], rx[2], ry[2]);
								area += rarea(rx[2], ry[2], x[i], y[i], rx[1], ry[1]);
								area += rarea(rx[1], ry[1], x[i], y[i], rx[0], ry[0]);
								if (area < rar) {
									out[0].push_back(rnc+j);
									out[1].push_back(z[i]);
								}
							} 
						}
					} else {
						break;
					}
				}
			}
		} else {
			for (size_t r=0; r<nr; r++) {
				double yrow = yFromRow(r);
				double ytop = yrow + offh;
				double ybot = yrow - offh;
				double rnc = r*nc;
				np = y.size();
				for (long i=(np-1); i>=0; i--) {
					if (y[i] > ytop) {
						y.pop_back(); // above current row
					} else if (y[i] >= ybot) {
						for (long j=(nc-1); j>=0; j--) {
							double dist = x[i] - xc[j];
							if (dist > w) {
								break;
							} else if (std::abs(dist) <= w) {
								out[0].push_back(rnc+j);
								out[1].push_back(z[i]);
							}
						}
					} else {
						break;
					}
				} 
			}
		}
	} else {

		if (rotated) {
			for (size_t r=0; r<nr; r++) {
				double yrow = yFromRow(r);
				double ytop = yrow + offh;
				double ybot = yrow - offh;
				double rnc = r*nc;
				ry[0] = yrow + oy[0];
				ry[1] = yrow + oy[1];
				ry[2] = yrow + oy[2];
				ry[3] = yrow + oy[3];
				np = y.size();
				std::vector<double> tmp0, tmp1;
				for (long i=(np-1); i>=0; i--) {
					if (y[i] > ytop) {
						y.pop_back(); // above current row
					} else if (y[i] >= ybot) {
						bool found = false;
						size_t minlim = 0;
						for (long j=(nc-1); j>=0; j--) {
							double dist = x[i] - xc[j];
							if (dist > bigw) {
								break;
							} else if (std::abs(dist) <= bigw) {
								rx[0] = xc[j] + ox[0];
								rx[1] = xc[j] + ox[1];
								rx[2] = xc[j] + ox[2];
								rx[3] = xc[j] + ox[3];
								// triangles apd, dpc, cpb, bpa
								double area = rarea(rx[0], ry[0], x[i], y[i], rx[3], ry[3]);
								area += rarea(rx[3], ry[3], x[i], y[i], rx[2], ry[2]);
								area += rarea(rx[2], ry[2], x[i], y[i], rx[1], ry[1]);
								area += rarea(rx[1], ry[1], x[i], y[i], rx[0], ry[0]);
								if (area < rar) {
									tmp0.push_back(rnc+j);
									tmp1.push_back(z[i]);
									found = true;
									minlim++;
								}
							} 
						}
						if (found) {
							if (minlim >= minpt) {
								out[0].insert(out[0].end(), tmp0.begin(), tmp0.end());
								out[1].insert(out[1].end(), tmp1.begin(), tmp1.end());
							}
							tmp0.resize(0);
							tmp1.resize(0);
							tmp0.reserve(10);
							tmp1.reserve(10);
						}
					} else {
						break;
					}
				}
			}
		} else {
			for (size_t r=0; r<nr; r++) {
				double yrow = yFromRow(r);
				double ytop = yrow + offh;
				double ybot = yrow - offh;
				double rnc = r*nc;
				np = y.size();
				std::vector<double> tmp0, tmp1;
				for (long i=(np-1); i>=0; i--) {
					if (y[i] > ytop) {
						y.pop_back(); // above current row
					} else if (y[i] >= ybot) {
						bool found = false;
						size_t minlim = 0;
						for (long j=(nc-1); j>=0; j--) {
							double dist = x[i] - xc[j];
							if (dist > w) {
								break;
							} else if (std::abs(dist) <= w) {
								tmp0.push_back(rnc+j);
								tmp1.push_back(z[i]);
								found = true;
								minlim++;
							}
						}
						if (found) {
							if (minlim >= minpt) {
								out[0].insert(out[0].end(), tmp0.begin(), tmp0.end());
								out[1].insert(out[1].end(), tmp1.begin(), tmp1.end());
							}
							tmp0.resize(0);
							tmp1.resize(0);
							tmp0.reserve(10);
							tmp1.reserve(10);
						}
					} else {
						break;
					}
				} 
			}
		}
	}
    return out;
}




std::vector<std::vector<double>> SpatRaster::win_circle(std::vector<double> x, std::vector<double> y, std::vector<double> z, std::vector<double> win, SpatOptions &opt) {

// the basic approach is from GDALGRID

	sortvecs(x, y, z);

    const double radius1 = win[0] * win[0]; 
    const double radius2 = win[1] * win[1];
	const double R12 = radius1 * radius2;

	double angle = std::fmod(win[2], 360.0);
	if (angle < 0) angle += 360.0;
    const bool rotated = angle != 0.0;
    angle = angle * M_PI / 180.0;

    //  coefficients for coordinate system rotation.
    const double cf1 = rotated ? cos(angle) : 0.0;
    const double cf2 = rotated ? sin(angle) : 0.0;

	size_t minpt = win[3] < 2 ? 1 : win[3];

	// for now assuming circles
    const double h = std::max(win[0], win[1]); 
//    const double w = h;

	const size_t nc = ncol();
	const size_t nr = nrow();

	size_t np = x.size() * 2;
	std::vector<std::vector<double>> out(2);
	out[0].reserve(np);
	out[1].reserve(np);

	std::vector<double> rx(4);
	std::vector<double> ry(4);

	std::vector<int64_t> cols(nc);
	std::iota(cols.begin(), cols.end(), 0);
	std::vector<double> xc = xFromCol(cols);


	if (minpt < 2) {
		for (size_t r=0; r<nr; r++) {
			double yrow = yFromRow(r);
			double ytop = yrow + h;
			double ybot = yrow - h;
			double rnc = r*nc;
			np = y.size();
			for (long i=(np-1); i>=0; i--) {
				if (y[i] > ytop) {
					y.pop_back(); // above current row
				} else if (y[i] >= ybot) {
					bool found = false;
					for (long j=(nc-1); j>=0; j--) {
						double RX = x[i] - xc[j];
						double RY = y[i] - yrow;
						if (rotated) {
							RX = RX * cf1 + RY * cf2;
							RY = RY * cf1 - RX * cf2;
						}
						if ((radius2 * RX * RX + radius1 * RY * RY) <= R12) {
							out[0].push_back(rnc+j);
							out[1].push_back(z[i]);
							found = true;
						} else {
							if (found) break;
						}
					}
				} else {
					break;
				}
			}
		}
	} else {
		for (size_t r=0; r<nr; r++) {
			double yrow = yFromRow(r);
			double ytop = yrow + h;
			double ybot = yrow - h;
			double rnc = r*nc;
			np = y.size();
			std::vector<double> tmp0, tmp1;
			for (long i=(np-1); i>=0; i--) {
				if (y[i] > ytop) {
					y.pop_back(); // above current row
				} else if (y[i] >= ybot) {
					bool found = false;
					size_t minlim = 0;
					for (long j=(nc-1); j>=0; j--) {
						double RX = x[i] - xc[j];
						double RY = y[i] - yrow;
						if (rotated) {
							RX = RX * cf1 + RY * cf2;
							RY = RY * cf1 - RX * cf2;
						}
						if ((radius2 * RX * RX + radius1 * RY * RY) <= R12) {
							tmp0.push_back(rnc+j);
							tmp1.push_back(z[i]);
							found = true;
							minlim++;
						}
					}
					if (found) {
						if (minlim >= minpt) {
							out[0].insert(out[0].end(), tmp0.begin(), tmp0.end());
							out[1].insert(out[1].end(), tmp1.begin(), tmp1.end());
						}
						tmp0.resize(0);
						tmp1.resize(0);
						tmp0.reserve(10);
						tmp1.reserve(10);
					}
				} else {
					break;
				}
			}
		}
	}
    return out;
}


// ===  moving window (end) ================================
