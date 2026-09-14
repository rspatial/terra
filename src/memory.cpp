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
#include "ram.h"
#include <limits>


namespace {

size_t size_gcd(size_t a, size_t b) {
	while (b != 0) {
		size_t t = a % b;
		a = b;
		b = t;
	}
	return a;
}

size_t size_lcm(size_t a, size_t b) {
	if ((a == 0) || (b == 0)) return 0;
	size_t g = size_gcd(a, b);
	if (a / g > (std::numeric_limits<size_t>::max() / b)) return 0;
	return (a / g) * b;
}

// Common GDAL block height to align read chunks to, or 0 if none / not useful.
// Strip files (height 1) and single-block-in-y files are ignored.
size_t file_align_rows(const std::vector<SpatRasterSource> &source, size_t nr) {
	size_t bh = 0;
	for (size_t i=0; i<source.size(); i++) {
		if (source[i].memory) continue;
		for (size_t j=0; j<source[i].blockrows.size(); j++) {
			int br = source[i].blockrows[j];
			if (br <= 1) continue;
			size_t b = (size_t) br;
			if (b >= nr) continue;
			if (bh == 0) {
				bh = b;
			} else {
				bh = size_lcm(bh, b);
				if ((bh == 0) || (bh >= nr)) return 0;
			}
		}
	}
	return bh;
}

// File row of raster row 0. 0 if no window or sources disagree.
size_t file_align_offrow(const std::vector<SpatRasterSource> &source) {
	size_t off = 0;
	bool seen = false;
	for (size_t i=0; i<source.size(); i++) {
		if (source[i].memory) continue;
		if (!source[i].hasWindow) continue;
		size_t o = source[i].window.off_row;
		if (!seen) {
			off = o;
			seen = true;
		} else if (o != off) {
			return 0;
		}
	}
	return seen ? off : 0;
}

} // namespace



bool SpatRaster::canProcessInMemory(SpatOptions &opt) {
	if (opt.get_todisk()) return false;
	double demand = size() * opt.ncopies;
	if (demand < opt.get_memmin()) {
		return true;
	}
	double supply;
	if (opt.get_memmax() > 0) {
		supply = opt.get_memmax() * opt.get_memfrac();
		supply = std::min(supply, availableRAM());
	} else {
		supply = availableRAM() * opt.get_memfrac();
	}
	std::vector<double> v;
	double maxsup = v.max_size(); //for 32 bit systems
	supply = std::min(supply, maxsup);
	return (demand < supply);
}


size_t SpatRaster::chunkSize(SpatOptions &opt) {
	double n = opt.ncopies;
	double frac = opt.get_memfrac();

	double demand = size() * n;
	if (demand < opt.get_memmin()) {
		return nrow();
	}

	double cells_in_row = ncol() * nlyr() * n;
	double supply;

	if (opt.get_memmax() > 0) {
		supply = opt.get_memmax() * opt.get_memfrac();
		supply = std::min(supply, availableRAM());
	} else {
		supply = availableRAM() * opt.get_memfrac();
	}
	double rows = supply * frac / cells_in_row;
	//double maxrows = 10000;
	//rows = std::min(rows, maxrows);
	size_t urows = floor(rows);
	urows = std::max(urows, (size_t)opt.minrows);
	if (urows < 1) return (1);
	if (urows < nrow()){
		return(urows);
	} else {
		return (nrow());
	}
}


std::vector<double> SpatRaster::mem_needs(SpatOptions &opt) {
	//returning bytes
	unsigned n = opt.ncopies;
	double memneed  = ncell() * (nlyr() * n);
	double memavail;
	if (opt.get_memmax() > 0) {
		memavail = std::min(opt.get_memmax(), availableRAM());
	} else {
		memavail = availableRAM();
	}
	double frac = opt.get_memfrac();
	double csize = chunkSize(opt);
	double inmem = canProcessInMemory(opt);
	std::vector<double> out = {memneed, memavail, frac, csize, inmem} ;
	return out;
}

//BlockSize SpatRaster::getBlockSize(unsigned n, double frac, unsigned steps) {
BlockSize SpatRaster::getBlockSize( SpatOptions &opt) {

	BlockSize bs;
	size_t nr = nrow();
	size_t cs = chunkSize(opt);
	bs.n = std::ceil(nr / double(cs));
	size_t steps = opt.get_steps();

	if (steps > 0) {
		if (steps > nr) {
			steps = nr;
		}
		bs.n = std::max(steps, bs.n);
	}

	cs = nr / bs.n;
	// methods that read beyond block boundaries (e.g. focal) set opt.minrows
	// to assure that all blocks, including the last one, have enough rows (#2138)
	size_t minr = std::min(opt.minrows, nr);
	if (cs < minr) {
		cs = minr;
	}

	// When GDAL reports a tiled/chunked layout, read full-width strips that
	// are a multiple of the native block height and start on a block boundary.
	// That way each file tile is decoded once. If a window does not start on
	// a file-block row, the first chunk is shortened to the next boundary.
	size_t lead = 0;
	size_t bh = file_align_rows(source, nr);
	if ((bh > 1) && (cs < nr)) {
		if (cs >= bh) {
			size_t down = (cs / bh) * bh;
			if (down >= minr) {
				cs = down;
			} else {
				size_t up = ((minr + bh - 1) / bh) * bh;
				if ((up > 0) && (up < nr)) {
					cs = up;
				}
			}
		}
		if (cs >= bh) {
			size_t phase = file_align_offrow(source) % bh;
			if (phase != 0) {
				lead = bh - phase;
				while ((lead < minr) && (lead + bh < nr)) {
					lead += bh;
				}
				if (lead >= nr) {
					lead = 0;
				}
			}
		}
	}

	if (lead == 0) {
		bs.n = std::ceil(nr / double(cs));
		size_t lastrows = nr - (bs.n - 1) * cs;
		if ((lastrows < minr) && (bs.n > 1)) {
			// merge a too-small remainder block with the one before it (#2138)
			bs.n -= 1;
			lastrows += cs;
		}
		bs.row = std::vector<size_t>(bs.n);
		bs.nrows = std::vector<size_t>(bs.n, cs);
		size_t r = 0;
		for (size_t i=0; i<bs.n; i++) {
			bs.row[i] = r;
			r += cs;
		}
		bs.nrows[bs.n-1] = lastrows;
	} else {
		size_t rest = nr - lead;
		size_t nmid = rest / cs;
		size_t rem = rest % cs;
		if ((rem > 0) && (rem < minr)) {
			if (nmid > 0) {
				rem += cs;
				nmid -= 1;
			} else {
				lead += rem;
				rem = 0;
			}
		}
		bs.n = 1 + nmid + (rem > 0 ? 1 : 0);
		bs.row.resize(bs.n);
		bs.nrows.resize(bs.n);
		bs.row[0] = 0;
		bs.nrows[0] = lead;
		size_t r = lead;
		for (size_t i=0; i<nmid; i++) {
			bs.row[i+1] = r;
			bs.nrows[i+1] = cs;
			r += cs;
		}
		if (rem > 0) {
			bs.row[bs.n-1] = r;
			bs.nrows[bs.n-1] = rem;
		}
	}
	return bs;
}

