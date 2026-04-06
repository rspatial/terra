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

#include <sstream>
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <iomanip>

#include "spatRasterMultiple.h"

#ifdef useGDAL
#include "ogr_spatialref.h"
#endif

#include "file_utils.h"
#include "string_utils.h"


static std::string basename_trunc(const std::string &path, size_t maxlen = 150) {
	std::string b = basename(path);
	if (b.size() > maxlen) {
		b = b.substr(0, maxlen) + "~";
	}
	return b;
}

static std::string pad_right(const std::string &s, size_t width) {
	if (s.size() >= width) return s;
	return std::string(width - s.size(), ' ') + s;
}

// Extent coordinates: SpatExtent uses full double precision; SpatVector / SpatRaster use fewer digits.
static std::string format_extent_double(double v, int significant_digits = 17) {
	std::ostringstream o;
	o << std::setprecision(significant_digits) << std::defaultfloat << v;
	return o.str();
}


// C++ equivalent of R's .name_or_proj4
static std::string crs_description(const std::string &wkt, const std::string &proj4) {
	if (wkt.empty() && proj4.empty()) return "";

#ifdef useGDAL
	OGRSpatialReference oSRS;
	OGRErr err = oSRS.SetFromUserInput(wkt.c_str());
	if (err == OGRERR_NONE) {
		std::string node = oSRS.IsGeographic() ? "geogcs" : "projcs";

		std::string name;
		const char *val = oSRS.GetAttrValue(node.c_str());
		if (val) name = val;

		std::string authority;
		val = oSRS.GetAuthorityName(node.c_str());
		if (val) authority = val;

		std::string code;
		val = oSRS.GetAuthorityCode(node.c_str());
		if (val) code = val;

		if (!name.empty() && name != "unknown" && name != "unnamed") {
			std::string r;
			if (proj4.substr(0, 13) == "+proj=longlat") {
				r = "lon/lat " + name;
			} else {
				r = name;
			}
			if (!code.empty()) {
				r += " (" + authority + ":" + code + ")";
			}
			return r;
		}
	}
#endif

	if (!proj4.empty()) return proj4;

	// try to extract the name from the WKT
	size_t p1 = wkt.find('"');
	if (p1 != std::string::npos) {
		size_t p2 = wkt.find('"', p1 + 1);
		if (p2 != std::string::npos) {
			return wkt.substr(p1 + 1, p2 - p1 - 1);
		}
	}
	return "";
}


// Match R show.R printDF(..., first=TRUE): names / types / first n data rows.
static std::string trunc_cell(const std::string &s, size_t mx) {
	if (s.size() <= mx + 2) return s;
	return s.substr(0, mx - 1) + "~";
}

static std::string format_df_double(double v) {
	char buf[64];
	std::snprintf(buf, sizeof(buf), "%g", v);
	return std::string(buf);
}

static std::string df_cell_as_string(const SpatDataFrame &df, size_t row, size_t col) {
	size_t it = df.itype[col];
	switch (it) {
		case 0: {
			double v = df.getDvalue(row, col);
			if (std::isnan(v)) return "NA";
			return format_df_double(v);
		}
		case 1: {
			long v = df.getIvalue(row, col);
			if (v == df.NAL) return "NA";
			return std::to_string(v);
		}
		case 2: {
			std::string v = df.getSvalue(row, col);
			if (v == df.NAS) return "NA";
			return v;
		}
		case 3: {
			int8_t b = df.getBvalue(row, col);
			if (b == 2) return "NA";
			return b ? "TRUE" : "FALSE";
		}
		case 4: {
			SpatTime_t t = df.getTvalue(row, col);
			if (t == df.NAT) return "NA";
			std::string step = df.tv[df.iplace[col]].step;
			std::vector<int> d = get_date(t);
			char buf[32];
			if (step == "days") {
				std::snprintf(buf, sizeof(buf), "%04d-%02d-%02d", d[0], d[1], d[2]);
			} else {
				std::snprintf(buf, sizeof(buf), "%04d-%02d-%02d %02d:%02d:%02d", d[0], d[1], d[2], d[3], d[4], d[5]);
			}
			return std::string(buf);
		}
		case 5: {
			SpatFactor f = df.getF(col);
			return f.getLabel(row);
		}
		default:
			return "?";
	}
}

static void append_spatvector_df_preview(std::ostringstream &s, const SpatDataFrame &df, size_t max_value_rows = 3) {
	size_t nr = df.nrow();
	size_t nc = df.ncol();
	if (nc == 0) return;

	// R printDF: at most 10 columns when there are more than 12.
	size_t nc_show = (nc > 12) ? 10 : nc;
	size_t nr_show = (max_value_rows == 0) ? 0 : std::min(max_value_rows, nr);
	size_t d2 = std::max(nc_show, (size_t)1);
	size_t mx = (size_t)std::max(15.0, 100.0 / (double)d2);

	std::vector<std::string> cnames(nc_show);
	std::vector<std::string> ctypes(nc_show);
	std::vector<std::vector<std::string>> vals(nr_show, std::vector<std::string>(nc_show));

	static const char *type_str[] = {"<num>", "<int>", "<chr>", "<lgl>", "", "<fact>"};

	for (size_t j = 0; j < nc_show; j++) {
		cnames[j] = trunc_cell(df.names[j], mx);
		size_t it = df.itype[j];
		if (it == 4) {
			ctypes[j] = (df.tv[df.iplace[j]].step == "days") ? "<Date>" : "<POSIXt>";
		} else {
			ctypes[j] = (it < 6) ? type_str[it] : "<?>";
		}
	}

	for (size_t r = 0; r < nr_show; r++) {
		for (size_t j = 0; j < nc_show; j++) {
			vals[r][j] = trunc_cell(df_cell_as_string(df, r, j), mx);
		}
	}

	std::vector<size_t> w(nc_show, 0);
	for (size_t j = 0; j < nc_show; j++) {
		w[j] = std::max(cnames[j].size(), ctypes[j].size());
		for (size_t r = 0; r < nr_show; r++) {
			w[j] = std::max(w[j], vals[r][j].size());
		}
	}

	auto join_padded = [&](const std::vector<std::string> &cells) {
		std::string line;
		for (size_t j = 0; j < cells.size(); j++) {
			if (j) line += ' ';
			line += pad_right(cells[j], w[j]);
		}
		return line;
	};

	s << "names       : " << join_padded(cnames);
	if (nc > nc_show) {
		s << "   (and " << (nc - nc_show) << " more)";
	}
	s << "\n";
	s << "type        : " << join_padded(ctypes) << "\n";

	if (nr_show == 0) return;

	static const char *VAL_PREFIX = "values      : ";
	static const char *CONT_PREFIX = "              ";

	s << VAL_PREFIX << join_padded(vals[0]) << "\n";
	for (size_t r = 1; r < nr_show; r++) {
		s << CONT_PREFIX << join_padded(vals[r]) << "\n";
	}
	if (nr > nr_show) {
		s << CONT_PREFIX << "...\n";
	}
}


std::string SpatExtent::show() {
	std::ostringstream s;
	s << "SpatExtent : " << format_extent_double(xmin) << ", " << format_extent_double(xmax) << ", "
	  << format_extent_double(ymin) << ", " << format_extent_double(ymax) << " (xmin, xmax, ymin, ymax)\n";
	return s.str();
}


std::string SpatVector::show() {
	std::ostringstream s;
	SpatExtent e = getExtent();
	size_t nr = nrow();
	size_t nc = ncol();

	s << "class       : SpatVector\n";
	s << "geometry    : " << type() << "\n";
	s << "dimensions  : " << nr << ", " << nc << "  (geometries, attributes)\n";
	s << "extent      : " << format_extent_double(e.xmin, 7) << ", " << format_extent_double(e.xmax, 7) << ", "
	  << format_extent_double(e.ymin, 7) << ", " << format_extent_double(e.ymax, 7) << "  (xmin, xmax, ymin, ymax)\n";

	if (!source.empty()) {
		std::string bn = basename_trunc(source);
		if (!bn.empty()) {
			std::string lyr = source_layer;
			std::string fn_no_ext = bn;
			size_t dot = fn_no_ext.rfind('.');
			if (dot != std::string::npos) fn_no_ext = fn_no_ext.substr(0, dot);
			if (!lyr.empty() && lyr != fn_no_ext) {
				s << "source      : " << bn << " (" << lyr << ")\n";
			} else {
				s << "source      : " << bn << "\n";
			}
		}
	}

	std::string wkt = srs.get("wkt");
	std::string proj4 = srs.get("proj4");
	std::string crs = crs_description(wkt, proj4);
	s << "coord. ref. : " << crs << "\n";

	if (nc > 0) {
		append_spatvector_df_preview(s, df);
	}

	return s.str();
}


std::string SpatRaster::show(bool one_based) {
	std::ostringstream s;
	size_t nr = nrow();
	size_t nc = ncol();
	size_t nl = nlyr();
	size_t nsr = nsrc();

	s << "class       : SpatRaster\n";
	s << "size        : " << nr << ", " << nc << ", " << nl
	  << "  (nrow, ncol, nlyr)\n";

	if ((nsr == 1) && source[0].is_multidim) {
		std::vector<std::vector<std::string>> dnms = dim_names();
		std::vector<std::vector<size_t>> dsz = dim_size();
		if (!dnms.empty() && !dnms[0].empty()) {
			std::string dnames, dsizes;
			for (size_t i = dnms[0].size(); i > 0; i--) {
				if (!dnames.empty()) { dnames += ", "; dsizes += ", "; }
				dnames += dnms[0][i-1];
				dsizes += std::to_string(dsz[0][i-1]);
			}
			s << "dimensions  : " << dnames << " (" << dsizes << "}\n";
		}
	}

	s << "resolution  : " << xres() << ", " << yres() << "  (x, y)\n";

	std::vector<bool> hw = hasWindow();
	bool any_win = false, all_win = true;
	for (bool w : hw) { any_win |= w; all_win &= w; }

	SpatExtent e = getExtent();
	if (any_win) {
		s << (all_win ? "window      : " : "extent (win): ");
	} else {
		s << "extent      : ";
	}
	s << format_extent_double(e.xmin, 7) << ", " << format_extent_double(e.xmax, 7) << ", "
	  << format_extent_double(e.ymin, 7) << ", " << format_extent_double(e.ymax, 7)
	  << "  (xmin, xmax, ymin, ymax)\n";

	std::string wkt = getSRS("wkt");
	std::string proj4 = getSRS("proj4");
	s << "coord. ref. : " << crs_description(wkt, proj4) << "\n";

	if (hasValues()) {
		const size_t mnr = 6;
		std::vector<std::string> ln = getNames();

		if (nl > mnr) {
			ln.resize(mnr);
			ln.push_back("...");
		}
		size_t lnmx = 60 / std::min(mnr, ln.size());
		for (auto &nm : ln) {
			if (nm.size() > lnmx + 2) {
				size_t mid = lnmx / 2;
				nm = nm.substr(0, mid) + "~" + nm.substr(nm.size() - mid);
			}
		}

		// sources
		std::vector<bool> mem = inMemory();
		std::vector<std::string> fnames = filenames();
		std::vector<std::string> srcs(fnames.size());
		for (size_t i = 0; i < fnames.size(); i++) {
			if (mem[i]) {
				srcs[i] = "memory";
			} else {
				std::string f = fnames[i];
				f.erase(std::remove(f.begin(), f.end(), '"'), f.end());
				if (f.substr(0, 5) != "HDF5:") {
					f = basename_trunc(f);
				}
				srcs[i] = f;
			}
		}

		bool all_mem = true;
		for (bool m : mem) all_mem &= m;

		if (all_mem) {
			s << "source(s)   : memory\n";
		} else if (nsr > 1) {
			size_t mxsrc = 3;
			std::vector<size_t> lbs = nlyrBySource();
			s << "sources     : " << srcs[0];
			if (lbs[0] != 1) s << " (" << lbs[0] << " layers)";
			s << "\n";
			for (size_t i = 1; i < std::min(mxsrc, nsr); i++) {
				s << "              " << srcs[i];
				if (lbs[i] != 1) s << " (" << lbs[i] << " layers)";
				s << "\n";
			}
			if (nsr > mxsrc) {
				if (nsr == mxsrc + 1) {
					s << "              " << srcs[mxsrc];
					if (lbs[mxsrc] != 1) s << " (" << lbs[mxsrc] << " layers)";
					s << "\n";
				} else {
					s << "              ... and " << (nsr - mxsrc) << " more sources\n";
				}
			}
		} else {
			s << "source      : " << srcs[0] << "\n";
		}

		if (!rgbtype.empty()) {
			std::vector<int> rgb = getRGB();
			s << "colors " << rgbtype << "  :";
			for (size_t i = 0; i < rgb.size(); i++) {
				if (i) s << ",";
				int v = rgb[i];
				if (one_based) {
					v++;
				}
				s << " " << v;
			}
			s << "\n";
		}

		std::vector<bool> hasct = hasColors();
		bool any_ct = false;
		for (bool c : hasct) any_ct |= c;
		if (any_ct) {
			s << "color table : ";
			bool first = true;
			for (size_t i = 0; i < hasct.size(); i++) {
				if (hasct[i]) {
					if (!first) s << ", ";
					s << (i + 1);
					first = false;
				}
			}
			s << "\n";
		}

		// varnames
		std::vector<std::string> varnms = getSourceNames();
		std::vector<std::string> lngnms = getLongSourceNames();
		bool show_varnames = false;
		for (size_t i = 0; i < varnms.size(); i++) {
			if (!varnms[i].empty()) {
				std::string fn_base = srcs[i];
				size_t dot = fn_base.rfind('.');
				if (dot != std::string::npos) fn_base = fn_base.substr(0, dot);
				if (fn_base != varnms[i]) { show_varnames = true; break; }
			}
		}
		if (show_varnames) {
			for (size_t i = 0; i < varnms.size(); i++) {
				if (i < lngnms.size() && !lngnms[i].empty()) {
					varnms[i] += " (" + lngnms[i] + ")";
				}
			}
			if (nsr == 1) {
				s << "varname     : " << varnms[0] << "\n";
			} else {
				s << "varnames    : " << varnms[0] << "\n";
				for (size_t i = 1; i < std::min(nsr, (size_t)3); i++) {
					s << "              " << varnms[i] << "\n";
				}
				if (nsr > 3) s << "              ...\n";
			}
		}

		// units
		std::vector<std::string> uts = getUnit();
		bool hasunits = false;
		for (auto &u : uts) { if (!u.empty()) { hasunits = true; break; } }

		// min/max
		std::vector<bool> hMM = hasRange();
		bool any_MM = false;
		for (bool h : hMM) any_MM |= h;

		std::vector<bool> isf = hasCategories();
		bool any_f = false;
		for (bool f : isf) any_f |= f;

		if (any_MM || any_f) {
			std::vector<double> rmin = range_min();
			std::vector<double> rmax = range_max();

			std::vector<std::string> minv(rmin.size());
			std::vector<std::string> maxv(rmax.size());
			for (size_t i = 0; i < rmin.size(); i++) {
				if (i < hMM.size() && !hMM[i]) {
					minv[i] = " ? ";
					maxv[i] = " ? ";
				} else {
					// Known range: show NaN like R when all cell values are missing (not "?")
					minv[i] = std::isnan(rmin[i]) ? "NaN" : double_to_string(rmin[i]);
					maxv[i] = std::isnan(rmax[i]) ? "NaN" : double_to_string(rmax[i]);
				}
			}

			if (nl > mnr) {
				minv.resize(mnr); minv.push_back("...");
				maxv.resize(mnr); maxv.push_back("...");
			}

			// column-align names, min, max
			size_t ncols = ln.size();
			std::vector<size_t> w(ncols);
			for (size_t i = 0; i < ncols; i++) {
				w[i] = ln[i].size();
				if (i < minv.size()) w[i] = std::max(w[i], minv[i].size());
				if (i < maxv.size()) w[i] = std::max(w[i], maxv[i].size());
			}
			for (size_t i = 0; i < ncols; i++) {
				ln[i] = pad_right(ln[i], w[i]);
				if (i < minv.size()) minv[i] = pad_right(minv[i], w[i]);
				if (i < maxv.size()) maxv[i] = pad_right(maxv[i], w[i]);
			}

			if (ncols == 1) {
				s << "name        : " << ln[0] << "\n";
				s << "min value   : " << minv[0] << "\n";
				s << "max value   : " << maxv[0] << "\n";
			} else {
				auto join = [](const std::vector<std::string> &v) {
					std::string out;
					for (size_t i = 0; i < v.size(); i++) {
						if (i) out += ", ";
						out += v[i];
					}
					return out;
				};
				s << "names       : " << join(ln)   << "\n";
				s << "min values  : " << join(minv)  << "\n";
				s << "max values  : " << join(maxv)  << "\n";
			}

			if (hasunits) {
				std::vector<std::string> utsu;
				for (auto &u : uts) {
					if (std::find(utsu.begin(), utsu.end(), u) == utsu.end())
						utsu.push_back(u);
				}
				if (utsu.size() == 1) {
					s << "unit        : " << utsu[0] << "\n";
				} else {
					if (nl > mnr) { uts.resize(mnr); uts.push_back("..."); }
					for (size_t i = 0; i < uts.size() && i < ncols; i++)
						uts[i] = pad_right(uts[i], w[i]);
					std::string uj;
					for (size_t i = 0; i < uts.size(); i++) {
						if (i) uj += ", ";
						uj += uts[i];
					}
					s << "unit        : " << uj << "\n";
				}
			}

		} else {
			// no min/max: just names (and possibly units)
			if (nl == 1) {
				s << "name        : " << ln[0] << "\n";
			} else {
				std::string nj;
				for (size_t i = 0; i < ln.size(); i++) {
					if (i) nj += ", ";
					nj += ln[i];
				}
				s << "names       : " << nj << "\n";
			}
			if (hasunits) {
				std::vector<std::string> utsu;
				for (auto &u : uts) {
					if (std::find(utsu.begin(), utsu.end(), u) == utsu.end())
						utsu.push_back(u);
				}
				if (utsu.size() == 1) {
					s << "unit        : " << utsu[0] << "\n";
				} else {
					if (nl > mnr) { uts.resize(mnr); uts.push_back("..."); }
					std::string uj;
					for (size_t i = 0; i < uts.size(); i++) {
						if (i) uj += ", ";
						uj += uts[i];
					}
					s << "unit        : " << uj << "\n";
				}
			}
		}
	}

	// depth
	if (hasDepth()) {
		std::string dname = getDepthName();
		std::string dunit = getDepthUnit();
		std::string label;
		if (dname == "depth") {
			label = (dunit.empty()) ? "" : "[" + dunit + "]: ";
		} else {
			if (dunit.empty() || dunit == "unknown") {
				label = dname + ": ";
			} else {
				label = dname + " [" + dunit + "]: ";
			}
		}
		std::vector<double> dpth = getDepth();
		std::sort(dpth.begin(), dpth.end());
		dpth.erase(std::unique(dpth.begin(), dpth.end()), dpth.end());
		s << "depth       : ";
		if (dpth.size() > 1) {
			s << dpth.front() << " to " << dpth.back()
			  << " (" << label << dpth.size() << " steps)";
		} else if (!dpth.empty()) {
			s << dpth[0];
		}
		s << "\n";
	}

	// time
	if (hasTime()) {
		std::string tims = source[0].timestep;
		std::string label = "time        ";
		if (tims == "yearmonths")      label = "time (ymnts)";
		else if (tims == "months")     label = "time (mnts) ";
		else if (tims == "years")      label = "time (years)";
		else if (tims == "days")       label = "time (days) ";
		else if (tims == "raw")        label = "time (raw)  ";

		std::vector<std::string> ts = getTimeStr(false, " ");
		if (!ts.empty()) {
			std::string first_t = ts.front();
			std::string last_t  = ts.back();
			// unique count
			std::vector<std::string> ut = ts;
			std::sort(ut.begin(), ut.end());
			ut.erase(std::unique(ut.begin(), ut.end()), ut.end());

			if (ut.size() > 1) {
				s << label << ": " << first_t << " to " << last_t;
				s << " (" << ut.size() << " steps)";
			} else {
				s << label << ": " << first_t;
			}
			s << "\n";
		}
	}

	return s.str();
}


// Match R .sources() for SpatRasterDataset (basename, strip quotes, colon paths, empty -> memory, unique).
static std::string stack_source_one(const std::string &raw) {
	std::string f = basename_trunc(raw);
	f.erase(std::remove(f.begin(), f.end(), '"'), f.end());
	if (f.find(':') != std::string::npos) {
		size_t c = f.find(':');
		f = basename_trunc(f.substr(0, c));
		f.erase(std::remove(f.begin(), f.end(), '"'), f.end());
	}
	if (f.empty()) f = "memory";
	return f;
}

static std::vector<std::string> sources_unique_stack(const std::vector<std::string> &filenames) {
	std::vector<std::string> out;
	for (const auto &raw : filenames) {
		std::string s = stack_source_one(raw);
		bool dup = false;
		for (const auto &u : out) {
			if (u == s) {
				dup = true;
				break;
			}
		}
		if (!dup) out.push_back(s);
	}
	return out;
}


std::string SpatRasterStack::show() {
	std::ostringstream s;
	s << "class       : SpatRasterDataset\n";
	size_t ns = nsds();
	s << "subdatasets : " << ns << "\n";
	if (ns == 0) return s.str();

	s << "dimensions  : " << nrow() << ", " << ncol() << " (nrow, ncol)\n";

	std::vector<size_t> nss = nlyr();
	std::vector<std::string> nlyr_txt;
	nlyr_txt.reserve(nss.size());
	for (size_t v : nss) nlyr_txt.push_back(std::to_string(v));
	if (nlyr_txt.size() > 10) {
		std::vector<std::string> t(nlyr_txt.begin(), nlyr_txt.begin() + 9);
		t.push_back("...");
		nlyr_txt = std::move(t);
	}
	s << "nlyr        : ";
	for (size_t i = 0; i < nlyr_txt.size(); i++) {
		if (i) s << ", ";
		s << nlyr_txt[i];
	}
	s << "\n";

	std::vector<double> xy = resolution();
	s << "resolution  : " << xy[0] << ", " << xy[1] << "  (x, y)\n";

	SpatExtent e = getExtent();
	s << "extent      : " << format_extent_double(e.xmin, 7) << ", " << format_extent_double(e.xmax, 7) << ", "
	  << format_extent_double(e.ymin, 7) << ", " << format_extent_double(e.ymax, 7) << "  (xmin, xmax, ymin, ymax)\n";

	std::string wkt = getSRS("wkt");
	std::string proj4 = getSRS("proj4");
	s << "coord. ref. : " << crs_description(wkt, proj4) << "\n";

	std::vector<std::string> srcs = sources_unique_stack(filenames());
	if (srcs.size() > 6) {
		srcs.resize(6);
		srcs.push_back("...");
	}
	s << "source(s)   : ";
	for (size_t i = 0; i < srcs.size(); i++) {
		if (i) s << ", ";
		s << srcs[i];
	}
	s << "\n";

	std::vector<std::string> ln = get_names();
	bool anynm = false;
	for (const auto &nm : ln) {
		if (!nm.empty()) {
			anynm = true;
			break;
		}
	}
	if (anynm) {
		if (ln.size() > 6) {
			ln.resize(6);
			ln.push_back("...");
		}
		s << "names       : ";
		for (size_t i = 0; i < ln.size(); i++) {
			if (i) s << ", ";
			s << ln[i];
		}
		s << "\n";
	}

	return s.str();
}


std::string SpatRasterCollection::show() {
	std::ostringstream s;
	s << "class       : SpatRasterCollection\n";
	size_t nr = size();
	s << "length      : " << nr << "\n";

	if (nr > 0) {
		std::vector<size_t> d = dims();
		size_t ncol = (nr > 14) ? 15 : nr;
		std::vector<std::vector<std::string>> cells(3, std::vector<std::string>(ncol));
		for (size_t j = 0; j < ncol; j++) {
			if (nr > 14 && j == 14) {
				cells[0][j] = cells[1][j] = cells[2][j] = "...";
			} else {
				cells[0][j] = std::to_string(d[j]);
				cells[1][j] = std::to_string(d[j + nr]);
				cells[2][j] = std::to_string(d[j + 2 * nr]);
			}
		}
		std::vector<size_t> colw(ncol, 0);
		for (size_t j = 0; j < ncol; j++) {
			for (size_t r = 0; r < 3; r++) {
				colw[j] = std::max(colw[j], cells[r][j].size());
			}
		}
		auto join_row = [&](size_t r) {
			std::string line;
			for (size_t j = 0; j < ncol; j++) {
				if (j) line += ", ";
				line += pad_right(cells[r][j], colw[j]);
			}
			return line;
		};
		s << "nrow        : " << join_row(0) << "\n";
		s << "ncol        : " << join_row(1) << "\n";
		s << "nlyr        : " << join_row(2) << "\n";
	}

	SpatExtent e = getExtent();
	s << "extent      : " << format_extent_double(e.xmin, 7) << ", " << format_extent_double(e.xmax, 7) << ", "
	  << format_extent_double(e.ymin, 7) << ", " << format_extent_double(e.ymax, 7) << "  (xmin, xmax, ymin, ymax)\n";

	if (!ds.empty()) {
		std::string wkt = ds[0].getSRS("wkt");
		std::string proj4 = ds[0].getSRS("proj4");
		std::string crs = crs_description(wkt, proj4);
		if (!crs.empty()) {
			s << "crs (first) : " << crs << "\n";
		}
	}

	std::vector<std::string> ln = get_names();
	bool anynm = false;
	for (const auto &nm : ln) {
		if (!nm.empty()) {
			anynm = true;
			break;
		}
	}
	if (anynm) {
		if (ln.size() > 6) {
			ln.resize(6);
			ln.push_back("...");
		}
		s << "names       : ";
		for (size_t i = 0; i < ln.size(); i++) {
			if (i) s << ", ";
			s << ln[i];
		}
		s << "\n";
	}

	return s.str();
}


std::string SpatVectorCollection::show() {
	std::ostringstream s;
	s << " class       : SpatVectorCollection\n";
	size_t n = v.size();
	s << " length      : " << n << "\n";
	size_t nn = (n > 15) ? 15 : n;
	if (n > 0) {
		for (size_t i = 0; i < nn; i++) {
			std::string gt = v[i].type();
			size_t nr = v[i].nrow();
			if (i == 0) {
				s << " geometry    : " << gt << " (" << nr << ")\n";
			} else {
				s << "               " << gt << " (" << nr << ")\n";
			}
		}
		if (n > nn) {
			s << "               "
			  << "   and " << (n - nn) << "more\n";
		}
		std::string wkt = v[0].getSRS("wkt");
		std::string proj4 = v[0].getSRS("proj4");
		std::string crs = crs_description(wkt, proj4);
		if (!crs.empty()) {
			s << " crs (first) : " << crs << "\n";
		}
		std::vector<std::string> nms = names;
		if (nms.size() > 10) {
			nms.resize(9);
			nms.push_back("...");
		}
		s << " names       : ";
		for (size_t i = 0; i < nms.size(); i++) {
			if (i) s << ", ";
			s << nms[i];
		}
		s << "\n";
	}
	return s.str();
}


std::string SpatVectorProxy::show() {
	std::ostringstream s;
	SpatExtent e = v.getExtent();
	size_t nr = v.nrow();
	size_t nc = v.ncol();

	s << " class       : SpatVectorProxy\n";
	s << " geometry    : " << v.type() << "\n";
	s << " dimensions  : " << nr << ", " << nc << "  (geometries, attributes)\n";
	s << " extent      : " << format_extent_double(e.xmin, 7) << ", " << format_extent_double(e.xmax, 7) << ", "
	  << format_extent_double(e.ymin, 7) << ", " << format_extent_double(e.ymax, 7) << "  (xmin, xmax, ymin, ymax)\n";

	std::string bn = basename_trunc(v.source);
	std::string lyr = v.source_layer;
	std::string fn_no_ext = bn;
	size_t dot = fn_no_ext.rfind('.');
	if (dot != std::string::npos) fn_no_ext = fn_no_ext.substr(0, dot);
	if (!lyr.empty() && lyr != fn_no_ext) {
		s << " source      : " << bn << " (" << lyr << ")\n";
	} else {
		s << " source      : " << bn << "\n";
	}

	std::string wkt = v.srs.get("wkt");
	std::string proj4 = v.srs.get("proj4");
	s << " coord. ref. : " << crs_description(wkt, proj4) << "\n";

	if (nc > 0) {
		append_spatvector_df_preview(s, v.df, 0);
	}

	return s.str();
}
