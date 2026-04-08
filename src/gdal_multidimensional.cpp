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

#include "spatRasterMultiple.h"
#include "proj.h"
#include "ogr_spatialref.h"
#include "gdal_priv.h"
#include "gdal.h"
#include "gdalio.h"
#include "crs.h"
#include "string_utils.h"
#include "file_utils.h"
#include "vecmath.h"
#include "recycle.h"
#include <stddef.h>
#include <cmath>
#include <algorithm>
#include <list>

//#include <cstdint>

namespace {


static bool md_is_col_dim_name(const std::string &nm) {
	std::string n = lower_case(lrtrim_copy(nm));
	if (n == "longitude" || n == "easting" || n == "eastings") return true;
	if (n == "lon" || n == "long" || n == "x" || n == "proj_x") return true;
	return false;
}

static bool md_is_row_dim_name(const std::string &nm) {
	std::string n = lower_case(lrtrim_copy(nm));
	if (n == "latitude" || n == "northing" || n == "northings") return true;
	if (n == "lat" || n == "y" || n == "proj_y") return true;
	return false;
}

static bool md_is_time_dim_name(const std::string &nm) {
	std::string n = lower_case(lrtrim_copy(nm));
	if (n == "t" || n == "time") return true;
	return in_string(n, "time");
}

static bool md_is_vertical_dim_name(const std::string &nm) {
	std::string n = lower_case(lrtrim_copy(nm));
	if (n == "z" || n == "lev" || n == "sigma") return true;
	return in_string(n, "layer") || in_string(n, "level") || in_string(n, "pressure")
		|| in_string(n, "depth") || in_string(n, "altitude") || in_string(n, "height")
		|| in_string(n, "plev");
}

static int md_find_col_dim(const std::vector<std::string> &dimnames) {
	for (size_t i = 0; i < dimnames.size(); i++) {
		if (md_is_col_dim_name(dimnames[i])) return (int) i;
	}
	return -1;
}

static int md_find_row_dim(const std::vector<std::string> &dimnames) {
	for (size_t i = 0; i < dimnames.size(); i++) {
		if (md_is_row_dim_name(dimnames[i])) return (int) i;
	}
	return -1;
}

static void md_classify_two_extra_dims(int &iz, int &it, size_t e0, size_t e1,
		const std::vector<std::string> &dimnames, const std::vector<std::string> &dimcalendar) {
	bool t0 = md_is_time_dim_name(dimnames[e0]) || !dimcalendar[e0].empty();
	bool t1 = md_is_time_dim_name(dimnames[e1]) || !dimcalendar[e1].empty();
	bool v0 = md_is_vertical_dim_name(dimnames[e0]);
	bool v1 = md_is_vertical_dim_name(dimnames[e1]);
	if (t0 && !t1 && !v0 && v1) {
		it = (int) e0;
		iz = (int) e1;
		return;
	}
	if (t1 && !t0 && !v1 && v0) {
		it = (int) e1;
		iz = (int) e0;
		return;
	}
	if (v0 && !v1 && !t0 && t1) {
		iz = (int) e0;
		it = (int) e1;
		return;
	}
	if (v1 && !v0 && !t1 && t0) {
		iz = (int) e1;
		it = (int) e0;
		return;
	}
	if (v0 && !v1) {
		iz = (int) e0;
		it = (int) e1;
		return;
	}
	if (v1 && !v0) {
		iz = (int) e1;
		it = (int) e0;
		return;
	}
	if (t0 && !t1) {
		it = (int) e0;
		iz = (int) e1;
		return;
	}
	if (t1 && !t0) {
		it = (int) e1;
		iz = (int) e0;
		return;
	}
	if (e0 < e1) {
		iz = (int) e0;
		it = (int) e1;
	} else {
		iz = (int) e1;
		it = (int) e0;
	}
}

// Row-major layer index over extra dimensions: extras[0] slowest, extras.back() fastest.
static void md_layer_to_indices(size_t layer, const std::vector<size_t> &sizes,
		std::vector<size_t> &idx) {
	idx.resize(sizes.size());
	size_t rem = layer;
	for (size_t j = 0; j < sizes.size(); j++) {
		size_t stride = 1;
		for (size_t t = j + 1; t < sizes.size(); t++) {
			stride *= sizes[t];
		}
		idx[j] = rem / stride;
		rem %= stride;
	}
}

// GDAL row-major: last dimension in the array (highest index) varies fastest in the
// read buffer. Terra uses one row = fixed latitude, columns = longitude (lon fastest
// within the row). If lon's index > lat's (ix > iy), the read buffer is already
// row-major in that sense. If lat's index > lon's (iy > ix), GDAL uses lat-fast
// layout buf[c*nrows+r] for (lat r, lon c); convert to buf[r*ncols+c].
static void md_reorder_spatial_gdal_to_terra(std::vector<double> &v, size_t offset, size_t nrows, size_t ncols) {
	if (nrows <= 1 || ncols <= 1) {
		return;
	}
	std::vector<double> tmp(nrows * ncols);
	for (size_t r = 0; r < nrows; r++) {
		for (size_t c = 0; c < ncols; c++) {
			tmp[r * ncols + c] = v[offset + c * nrows + r];
		}
	}
	std::copy(tmp.begin(), tmp.end(), v.begin() + offset);
}

}  // namespace

//#if INTPTR_MAX != INT32_MAX
//    #define IS_64_BIT
//#endif


bool parse_ncdf_time(SpatRasterSource &s, const std::string unit, const std::string calendar, std::vector<double> raw, std::string &msg) {

	std::vector<int64_t> out;
	std::string origin = unit;
	bool years = false;
	bool yearsbp = false;
	bool yearmonths = false;
	bool months = false;
	bool days = false;
	bool hours = false;
	bool minutes = false;
	bool seconds = false;
	bool foundorigin = false;
	std::string step;
	
	lowercase(origin);
	if ((origin.find("seconds")) != std::string::npos) {
		seconds = true;
	} else if ((origin.find("minutes")) != std::string::npos) {
		minutes = true;
	} else if ((origin.find("hours")) != std::string::npos) {
		hours = true;
	} else if ((origin.find("days")) != std::string::npos) {
		days = true;
	} else if ((origin.find("months since")) != std::string::npos) {
		yearmonths = true;
		foundorigin = true;
	} else if ((origin.find("months")) != std::string::npos) {
		months = true;
		foundorigin = true;
	} else if ((origin.find("years before present")) != std::string::npos) {
		yearsbp = true;
		foundorigin = true;
	} else if ((origin.find("years")) != std::string::npos) {
		years = true;
	}

	if (!foundorigin) {
		size_t pos;		
		if ((pos = origin.find("from")) != std::string::npos) {
			origin.erase(0, pos + 5);
			foundorigin = true;
		} else if ((pos = origin.find("since")) != std::string::npos) {
			origin.erase(0, pos + 6);
			foundorigin = true;
		}
	}

	SpatTime_t offset = 0;
	if (foundorigin) {
		step = "seconds";
		out.reserve(raw.size());
		std::string cal = "366";
		if (calendar == "360_day" || calendar == "360 day") {
			cal = "360";
		} else if (calendar == "noleap" || calendar == "365_day" || calendar == "365 day") {
			cal = "365";		
		} else if (calendar =="gregorian" || calendar =="proleptic_gregorian" || calendar=="standard" || calendar == "julian") {
			cal = "366";
		} else if (!(months || years || yearmonths || yearsbp)) {
			//cal = "366";
			msg = "unknown calendar (assuming standard): " + calendar;			
		}

		// this shortcut means that 360/noleap calendars loose only have dates, no time
		// to be refined
		if ((hours || minutes || seconds) && (cal == "360")) {
			int div = 24;
			double add = 0;
			std::vector<int> ymd = getymd(origin);
			if (hours) {
				hours = false;
				add = ymd[3] + ymd[4] / 60 + ymd[5] / 3600; 
			} else if (minutes) {
				div = 1440; // 24 * 60
				add = ymd[3] * 60 + ymd[4] + ymd[5] / 60; 
				minutes = false;
			} else if (seconds) {
				div = 86400; // 24 * 3600
				add = ymd[3] * 3600 + ymd[4] * 60 + ymd[5]; 
				seconds = false;
			}
			for (size_t i=0; i<raw.size(); i++) {
				raw[i] = (raw[i]+add) / div; 
			}
			days = true;
		} 
	
		if (days) {
			step = "days";
			std::vector<int> ymd = getymd(origin);
			if (cal == "365") {
				for (size_t i=0; i<raw.size(); i++) out.push_back(
					get_time_noleap(ymd[0], ymd[1], ymd[2], 0, 0, 0, raw[i], "days"));
			} else if (cal == "360") {
				for (size_t i=0; i<raw.size(); i++) out.push_back(
					time_from_day_360(ymd[0], ymd[1], ymd[2], raw[i]));
			} else {
				for (size_t i=0; i<raw.size(); i++) out.push_back(
					time_from_day(ymd[0], ymd[1], ymd[2], raw[i]));
			}
		} else if (hours) {
			if (cal == "365") {
				std::vector<int> ymd = getymd(origin);
				for (size_t i=0; i<raw.size(); i++) out.push_back(
						get_time_noleap(ymd[0], ymd[1], ymd[2], ymd[3], 0, 0, raw[i], "hours")
					);
			} else {
				offset = get_time_string(origin);
				for (size_t i=0; i<raw.size(); i++) out.push_back(raw[i]*3600+offset);
			}
		} else if (minutes) {
			if (cal == "365") {
				std::vector<int> ymd = getymd(origin);
				for (size_t i=0; i<raw.size(); i++) out.push_back(
						get_time_noleap(ymd[0], ymd[1], ymd[2], ymd[3], ymd[4], 0, raw[i], "minutes")
					);
			} else {
				offset = get_time_string(origin);
				for (size_t i=0; i<raw.size(); i++) out.push_back(60*raw[i]+offset);
			}
		} else if (seconds) {
			if (cal == "365") {
				std::vector<int> ymd = getymd(origin);
				for (size_t i=0; i<raw.size(); i++) out.push_back(
						get_time_noleap(ymd[0], ymd[1], ymd[2], ymd[3], ymd[4], 0, raw[i], "minutes")
					);
			} else {
				offset = get_time_string(origin);
				for (size_t i=0; i<raw.size(); i++) out.push_back(raw[i]+offset);
			}
		} else if (years) {
			step = "years";
			int syear = getyear(origin);
			for (size_t i=0; i<raw.size(); i++) out.push_back(get_time(syear+raw[i], 6, 30, 0, 0, 0));
		} else if (yearsbp) {
			step = "years";
			int syear = 1950;
			for (size_t i=0; i<raw.size(); i++) out.push_back(get_time(syear+raw[i], 6, 30, 0, 0, 0));
		} else if (yearmonths) {
			step = "yearmonths";
			int syear = getyear(origin);
			for (size_t i=0; i<raw.size(); i++) {
				long year = std::floor(raw[i] / 12.0);
				int month = raw[i] - (12 * year) + 1;
				out.push_back(get_time(syear+year, month, 15, 0, 0, 0));
			}
		} else if (months) {
			step = "months";
			// check for 0..11 range
			int zero = vmin(raw, true) == 0.0;
			for (size_t i=0; i<raw.size(); i++) {
				unsigned m = std::ceil(raw[i] + zero);
				out.push_back(get_time(1970, m, 15, 0, 0, 0));
			}
		} else {
			step = "raw";
			for (size_t i=0; i<raw.size(); i++) out.push_back(raw[i]);
		}

		if (step == "days" && raw.size() == out.size()) {
			for (size_t i = 0; i < raw.size(); i++) {
				double fr = raw[i] - std::floor(raw[i]);
				if (fr > 1e-6 && fr < 1.0 - 1e-6) {
					step = "seconds";
					break;
				}
			}
		}
	}

	s.time = out;
	s.timestep = step;
	s.hasTime = true;
	return true;
}




#if (GDAL_VERSION_MAJOR > 3) || (GDAL_VERSION_MAJOR == 3 && GDAL_VERSION_MINOR >= 4) // && defined(IS_64_BIT)


std::vector<std::string> GetArrayNames(std::shared_ptr<GDALGroup> x, bool filter) {
// FROM GDAL 3.11 (while not widely available).
// * Author:   Even Rouault <even.rouault at spatialys.com>
// * Copyright (c) 2019, Even Rouault <even.rouault at spatialys.com>
    std::vector<std::string> ret;
    std::list<std::shared_ptr<GDALGroup>> stackGroups;
    stackGroups.push_back(nullptr);  // nullptr means this
    while (!stackGroups.empty()) {
        std::shared_ptr<GDALGroup> groupPtr = std::move(stackGroups.front());
        stackGroups.erase(stackGroups.begin());
        const GDALGroup *poCurGroup = groupPtr ? groupPtr.get() : x.get();
        for (const std::string &arrayName :  poCurGroup->GetMDArrayNames(nullptr)) {
            std::string osFullName = poCurGroup->GetFullName();
            if (!osFullName.empty() && osFullName.back() != '/') 
                osFullName += '/';
            osFullName += arrayName;
            ret.push_back(std::move(osFullName));
        }
        auto insertionPoint = stackGroups.begin();
        for (const auto &osSubGroup : poCurGroup->GetGroupNames(nullptr)) {
            auto poSubGroup = poCurGroup->OpenGroup(osSubGroup);
            if (poSubGroup)
                stackGroups.insert(insertionPoint, std::move(poSubGroup));
        }
    }
	if (filter) {
		return ncdf_filternames(ret);
	}
    return ret;
}

// Arrays usable as SpatRaster md sources: at least 2 dimensions; higher dimension count first.
static std::vector<std::string> md_arrays_usable_for_raster(
		std::shared_ptr<GDALGroup> poRootGroup,
		const std::vector<std::string> &candidates) {
	struct NameDim {
		std::string name;
		size_t ndim;
	};
	std::vector<NameDim> tmp;
	tmp.reserve(candidates.size());
	std::string startgroup = "";
	for (const std::string &nm : candidates) {
		auto poVar = poRootGroup->ResolveMDArray(nm.c_str(), startgroup, nullptr);
		if (!poVar) {
			continue;
		}
		size_t nd = poVar->GetDimensions().size();
		if (nd >= 2) {
			tmp.push_back({poVar->GetFullName(), nd});
		}
	}
	std::stable_sort(tmp.begin(), tmp.end(), [](const NameDim &a, const NameDim &b) {
		if (a.ndim != b.ndim) {
			return a.ndim > b.ndim;
		}
		return a.name < b.name;
	});
	std::vector<std::string> out;
	out.reserve(tmp.size());
	for (const auto &p : tmp) {
		out.push_back(p.name);
	}
	return out;
}


static bool md_fill_source_from_marray(
	SpatRaster &parent,
	const std::string &fname,
	const std::string &array_request_name,
	std::shared_ptr<GDALMDArray> poVar,
	std::vector<std::string> options,
	bool noflip,
	bool guessCRS,
	bool errors_are_fatal,
	SpatRasterSource &s) {

	auto fail = [&](const std::string &msg) -> bool {
		if (errors_are_fatal) {
			parent.setError(msg);
		} else {
			parent.addWarning(std::string("skipped multidimensional array: ") + array_request_name + " (" + msg + ")");
		}
		return false;
	};

	s.m_arrayname = poVar->GetFullName();

// dimensions
	std::vector<size_t> dimcount;
	std::vector<std::string> dimnames, dimunits, dimcalendar;
	std::vector<std::vector<double>> dimvals;
	std::vector<std::shared_ptr<GDALDimension>> dimData = poVar->GetDimensions();
	size_t ndim = dimData.size();
	dimvals.reserve(ndim);
	dimcalendar.reserve(ndim);
	
	
    for (size_t i=0; i<ndim; i++) {
		size_t n = dimData[i]->GetSize();
        dimcount.push_back(n);
        std::string name = static_cast<std::string>(dimData[i]->GetName());
		dimnames.push_back(name);
		std::vector<GUInt64> start(1, 0);
		std::vector<size_t> count = {n};
		dimvals.push_back(std::vector<double>(n));

		const auto indvar = dimData[i]->GetIndexingVariable();
		
		if (indvar == NULL) {
			dimvals[i].resize(n);
			std::iota(dimvals[i].begin(), dimvals[i].end(), 1);			
			dimunits.push_back("");
			dimcalendar.push_back("");
		} else {
			dimunits.push_back(static_cast<std::string>(indvar->GetUnit()));
			std::string cal = "";
			auto pcal = indvar->GetAttribute("calendar");
			if (pcal) cal = pcal->ReadAsString();
			dimcalendar.push_back(cal);
			indvar->Read(start.data(), count.data(), nullptr, nullptr, GDALExtendedDataType::Create(GDT_Float64), &dimvals[i][0]);
		}
	}

	s.m_ndims = dimcount.size();
	if (s.m_ndims < 2) {
		return fail("insufficient number of dimensions");
	}
	s.source_name = s.m_arrayname;
	if (!s.source_name.empty() && (s.source_name.front() == '/')) {
		s.source_name.erase(0, 1); 
    }
	
	auto lname = poVar->GetAttribute("long_name");
	if (lname) s.source_name_long = lname->ReadAsString();

	auto mval = poVar->GetAttribute("missing_value");
	if (mval) {
		s.m_missing_value = mval->ReadAsDouble();
		s.m_hasNA = true;
	} else {
		s.m_missing_value = poVar->GetNoDataValueAsDouble(&s.m_hasNA);
	}

	int ix = md_find_col_dim(dimnames);
	int iy = md_find_row_dim(dimnames);
	if (ix < 0 || iy < 0 || ix == iy) {
		ix = (int) ndim - 1;
		iy = (int) ndim - 2;
	}
	int it = -1;
	int iz = -1;
	std::vector<size_t> dimmap_extras;

	if (ndim == 2) {
		s.m_order = {(size_t) ix, (size_t) iy};
	} else {
		std::vector<size_t> extra;
		for (size_t k = 0; k < ndim; k++) {
			if ((int) k != ix && (int) k != iy) {
				extra.push_back(k);
			}
		}
		if (extra.size() == 1) {
			size_t e = extra[0];
			if (md_is_time_dim_name(dimnames[e]) || !dimcalendar[e].empty()) {
				it = (int) e;
			} else {
				iz = (int) e;
			}
			dimmap_extras.push_back(e);
		} else if (extra.size() == 2) {
			md_classify_two_extra_dims(iz, it, extra[0], extra[1], dimnames, dimcalendar);
			if (iz >= 0) {
				dimmap_extras.push_back(iz);
			}
			if (it >= 0) {
				dimmap_extras.push_back(it);
			}
		} else {
			dimmap_extras = extra;
			for (size_t e : extra) {
				if (it < 0 && (md_is_time_dim_name(dimnames[e]) || !dimcalendar[e].empty())) {
					it = (int) e;
				}
			}
			for (size_t e : extra) {
				if ((int) e != it && iz < 0 && md_is_vertical_dim_name(dimnames[e])) {
					iz = (int) e;
					break;
				}
			}
		}
		s.m_order.clear();
		s.m_order.push_back((size_t) ix);
		s.m_order.push_back((size_t) iy);
		for (size_t e : dimmap_extras) {
			s.m_order.push_back(e);
		}
	}

	std::vector<int> dimmap;
	dimmap.reserve(2 + dimmap_extras.size());
	dimmap.push_back(ix);
	dimmap.push_back(iy);
	for (size_t e : dimmap_extras) {
		dimmap.push_back((int) e);
	}

	for (size_t ii = 0; ii < ndim; ii++) {
		if ((int) ii != ix && (int) ii != iy) {
			continue;
		}
		if (dimvals[ii].size() <= 2) {
			continue;
		}
		const auto indvar2 = dimData[ii]->GetIndexingVariable();
		if (indvar2 == NULL) {
			continue;
		}
		double res = dimvals[ii][1] - dimvals[ii][0];
		if (!indvar2->IsRegularlySpaced(dimvals[ii][0], res)) {
			return fail(dimnames[ii] + " is not regularly spaced");
		}
	}

	SpatExtent e;
 	s.ncol = dimcount[ix];
	s.nrow = dimcount[iy];

	double start = dimvals[ix][0];
	double end = dimvals[ix][dimvals[ix].size()-1];
	double res = (end - start) / (s.ncol-1);
	e.xmin = start - 0.5 * res;
	e.xmax = end + 0.5 * res;

	start = dimvals[iy][0];
	end = dimvals[iy][dimvals[iy].size()-1];
	res = (end - start) / (s.nrow-1);
	e.ymax = end + 0.5 * res;
	e.ymin = start - 0.5 * res;
	
	s.flipped = false;
	if ((!noflip) && (e.ymin > e.ymax)) {
		std::swap(e.ymin, e.ymax);
		s.flipped = true;
	}
//	s.m_names.push_back(dimnames[ix]);
//	s.m_names.push_back(dimnames[iy]);


	std::vector<size_t> extra_sizes;
	for (size_t e : dimmap_extras) {
		extra_sizes.push_back(dimcount[e]);
	}
	s.nlyr = 1;
	for (size_t sz : extra_sizes) {
		s.nlyr *= sz;
	}

	std::vector<int64_t> time_coord;
	if (it >= 0) {
		std::string msg;
		std::string cal = "standard";
		if (static_cast<size_t>(it) < dimcalendar.size() && (!dimcalendar[it].empty())) {
			cal = dimcalendar[it];
		}
		parse_ncdf_time(s, dimunits[it], cal, dimvals[it], msg);
		time_coord = s.time;
	}

	size_t pos_it = (size_t) -1;
	size_t pos_iz = (size_t) -1;
	for (size_t j = 0; j < dimmap_extras.size(); j++) {
		if ((int) dimmap_extras[j] == it) {
			pos_it = j;
		}
		if ((int) dimmap_extras[j] == iz) {
			pos_iz = j;
		}
	}

	s.nlyrfile = s.nlyr;
	s.resize(s.nlyr);
	s.layers.resize(s.nlyr);
    std::iota(s.layers.begin(), s.layers.end(), 0);

	std::vector<size_t> idx;
	if (it >= 0 && pos_it != (size_t) -1 && !time_coord.empty()) {
		for (size_t L = 0; L < s.nlyr; L++) {
			md_layer_to_indices(L, extra_sizes, idx);
			s.time[L] = time_coord[idx[pos_it]];
		}
	}
	if (iz >= 0 && pos_iz != (size_t) -1) {
		s.depthname = dimnames[iz];
		s.hasDepth = true;
		const std::vector<double> &dvz = dimvals[iz];
		for (size_t L = 0; L < s.nlyr; L++) {
			md_layer_to_indices(L, extra_sizes, idx);
			s.depth[L] = dvz[idx[pos_iz]];
		}
	}
	
	for (size_t i = 0; i < dimmap.size(); i++) {
		s.m_dims.push_back((size_t) dimmap[i]);
	}
	s.extent = e;


	std::string wkt = "";
	auto srs = poVar->GetSpatialRef();

	if (srs != NULL) {
		char *cp;
		const char *options[3] = { "MULTILINE=YES", "FORMAT=WKT2", NULL };
		OGRErr err = srs->exportToWkt(&cp, options);
		if (err == OGRERR_NONE) {
			wkt = std::string(cp);
		}
		CPLFree(cp);
	} 
	if (wkt.empty()) {
		// temporary work-around for https://github.com/rspatial/terra/issues/2068
		std::vector<std::string> ops;
		try {
			std::vector<std::string> empty_dom;
			SpatRasterStack rstack(fname, {0}, true, ops, true, true, empty_dom);
			wkt = rstack.getSRS("wkt");
		} catch(...) {}
	} 
	
	if (guessCRS && wkt.empty()) {
		
		if (s.extent.xmin >= -181 && s.extent.xmax <= 361 && s.extent.ymin >= -91 && s.extent.ymax <= 91) {
			wkt = "OGC:CRS84";
			s.parameters_changed = true;
			parent.addWarning("guessed crs");
		}
	}
	std::string msg = "";
	if (!s.srs.set({wkt}, msg)) {
		parent.addWarning(msg);
	}	

	bool app_so = true;
	size_t opsz = options.size();
	if (opsz > 0) {
		if (options[opsz-1] == "so=false") {
			app_so = false;
			options.resize(opsz-1); 
		}
	}

	if (app_so) {
		bool hasScale=false;
		double scale = poVar->GetScale(&hasScale, nullptr);
		if (scale == 1) hasScale = false;
		bool hasOffset=false;
		double offset = poVar->GetOffset(&hasOffset, nullptr);
		if (offset == 0) hasOffset = false;
		if (hasScale || hasOffset) {
			s.has_scale_offset = std::vector<bool>(s.nlyr, true);
			s.offset = std::vector<double>(s.nlyr, offset);
			s.scale = std::vector<double>(s.nlyr, scale);
		}
	}
	
	s.rotated = false;
	s.memory = false;
	s.filename = fname;
	s.hasValues = true;
	s.unit = std::vector<std::string>(s.nlyr, poVar->GetUnit());
	s.is_multidim = true;

// layer names
	std::vector<std::string> nms;
	std::vector<std::string> arn = strsplit_last(s.m_arrayname, "/");
	std::string arname = arn[arn.size()-1];
	nms.resize(s.nlyr);
	if (dimmap_extras.empty()) {
		if (s.nlyr > 1) {
			for (size_t i = 0; i < s.nlyr; i++) {
				nms[i] = arname + "-" + std::to_string(i + 1);
			}
		} else {
			nms[0] = arname;			
		}
	} else {
		for (size_t L = 0; L < s.nlyr; L++) {
			md_layer_to_indices(L, extra_sizes, idx);
			std::string nm = arname;
			bool name_has_dim = false;
			bool skipped_time_for_name = false;
			for (size_t j = 0; j < dimmap_extras.size(); j++) {
				if (it >= 0 && (int) dimmap_extras[j] == it) {
					skipped_time_for_name = true;
					continue;
				}
				name_has_dim = true;
				nm += "_" + dimnames[dimmap_extras[j]] + "="
					+ double_to_string(dimvals[dimmap_extras[j]][idx[j]]);
			}
			if (!name_has_dim && (s.nlyr > 1)) {
				// Only non-spatial dim is time (omitted from label): number layers
				nm += "_" + std::to_string(L + 1);
			} else if (skipped_time_for_name && (s.nlyr > 1)) {
				// Time is in metadata, not in the label; add 1-based time step so
				// names match rast(, md=FALSE), e.g. t2m_expver=1_1 .. _24
				size_t tidx = (pos_it != (size_t) -1) ? (idx[pos_it] + 1) : (L + 1);
				nm += "_" + std::to_string(tidx);
			}
			nms[L] = nm;
		}
	}
	s.names = nms;

// time
	s.m_size = dimcount;
	s.m_names = dimnames;

	return true;
}


bool SpatRaster::constructFromFileMulti(std::string fname, std::vector<int> subds, std::vector<std::string> subname, std::vector<std::string> drivers, std::vector<std::string> options, std::vector<int> dims, bool noflip, bool guessCRS, std::vector<std::string> domains) {
//	(void) dims;
//	(void) domains;


	char ** drvs = NULL;
	for (size_t i=0; i<drivers.size(); i++) {
		drvs = CSLAddString(drvs, drivers[i].c_str());
	}

    auto poDataset = std::unique_ptr<GDALDataset>(GDALDataset::Open(fname.c_str(), GDAL_OF_MULTIDIM_RASTER, drvs));
	CSLDestroy(drvs);
	drvs = NULL;
    if( !poDataset ) {
		if (!file_exists(fname)) {
			setError("file does not exist: " + fname);
		} else if (drivers.size() > 0) {
			setError("cannot read multidim from this file or with this driver");
		} else {
			setError("cannot read multidim from this file");
        }
		return false;
    }

	std::shared_ptr<GDALGroup> poRootGroup = poDataset->GetRootGroup();
    if( !poRootGroup ) {
		setError("dataset has no root group");
		return false;
    }

	std::vector<std::string> anms_all = GetArrayNames(poRootGroup, true);
	if (anms_all.empty()) {
		setError("no MD arrays found in file");
		return false;
	}

	std::vector<std::string> anms = md_arrays_usable_for_raster(poRootGroup, anms_all);
	if (anms.empty()) {
		setError("file has no array with at least 2 dimensions");
		return false;
	}

	// Variables to load (like rast(, md=FALSE) combining all SDS)
	std::vector<std::string> arrays_to_use;
	if ((subname.size() > 0) && (!subname[0].empty())) {
		arrays_to_use.push_back(subname[0]);
	} else if (subds.size() > 0 && subds[0] >= 0) {
		if ((size_t) subds[0] >= anms.size()) {
			setError("array index is out of range (there are " + std::to_string(anms.size()) + " arrays)");
			return false;
		} else {
			arrays_to_use.push_back(anms[subds[0]]);
		}
	} else {
		arrays_to_use = anms;
	}

	const bool single_var = (arrays_to_use.size() == 1);
	SpatOptions opt;
	size_t nvar_ok = 0;
	size_t max_nlyr_var = 0;
	size_t min_nlyr_var = (size_t) -1;

	for (size_t ai = 0; ai < arrays_to_use.size(); ai++) {
		std::string startgroup = "";
		auto poVar = poRootGroup->ResolveMDArray(arrays_to_use[ai].c_str(), startgroup, nullptr);
		if (!poVar) {
			if (single_var) {
				setError(std::string("cannot find array: \"") + arrays_to_use[ai] + "\".\nAvailable arrays: " + concatenate(anms, ", "));
				return false;
			}
			addWarning(std::string("skipped multidimensional array (not found): ") + arrays_to_use[ai]);
			continue;
		}
		{
			size_t nd = poVar->GetDimensions().size();
			if (nd < 2) {
				if (single_var) {
					setError("array \"" + arrays_to_use[ai] + "\" has " + std::to_string(nd) +
						" dimension(s); rast(, md=TRUE) requires at least 2 dimensions");
					return false;
				}
				addWarning(std::string("skipped multidimensional array (<2 dims): ") + arrays_to_use[ai]);
				continue;
			}
		}

		SpatRasterSource s;
		s.open_drivers = drivers;
		std::vector<std::string> opts_copy = options;
		if (!md_fill_source_from_marray(*this, fname, arrays_to_use[ai], poVar, std::move(opts_copy), noflip, guessCRS, single_var, s)) {
			if (single_var) {
				return false;
			}
			continue;
		}

		if (nvar_ok == 0) {
			setSource(s);
		} else {
			SpatRaster chunk;
			chunk.setSource(s);
			if (!chunk.compare_geom(*this, false, false, 0.1)) {
				addWarning(std::string("skipped multidimensional array (different geometry): ") + arrays_to_use[ai]);
				continue;
			}
			addSource(chunk, false, opt);
		}
		max_nlyr_var = std::max(max_nlyr_var, (size_t) s.nlyr);
		if (min_nlyr_var == (size_t) -1) {
			min_nlyr_var = s.nlyr;
		} else {
			min_nlyr_var = std::min(min_nlyr_var, (size_t) s.nlyr);
		}
		nvar_ok++;
	}

	if (nvar_ok == 0) {
		setError(std::string("could not load multidimensional data. Arrays: ") + concatenate(anms, ", "));
		return false;
	}


//	if (verbose) {
	if (arrays_to_use.size() > 1 && max_nlyr_var > 1) {
		std::string w = "combined " + std::to_string(nvar_ok) + " variables";
		w += " (";
		size_t nshow = std::min(source.size(), (size_t) 6);
		for (size_t i = 0; i < nshow; i++) {
			if (i > 0) {
				w += ", ";
			}
			w += source[i].source_name.empty() ? source[i].m_arrayname : source[i].source_name;
		}
		if (source.size() > 5) {
			w += ", ...)";
		} else {
			w += ")";
		}
		addWarning(w);
	}
//	}
	return true;
}


bool SpatRaster::readStartMulti(size_t src) {

	char ** drvs = NULL;
	for (size_t i=0; i<source[src].open_drivers.size(); i++) {
		drvs = CSLAddString(drvs, source[src].open_drivers[i].c_str());
	}

    auto poDataset = std::unique_ptr<GDALDataset>(GDALDataset::Open(source[src].filename.c_str(), GDAL_OF_MULTIDIM_RASTER, drvs));
    if( !poDataset ) {
		setError("not a good dataset");
        return false;
    }


	std::shared_ptr<GDALGroup> poRootGroup = poDataset->GetRootGroup();
    if( !poRootGroup ) {
		setError("no roots");
		return false;
    }

	std::string startgroup="";
	auto poVar = poRootGroup->ResolveMDArray(source[src].m_arrayname.c_str(), startgroup, nullptr);

//    auto poVar = poRootGroup->OpenMDArray(source[src].m_arrayname.c_str());
    if( !poVar )   {
		setError("cannot find: " + source[src].m_arrayname);
		return false;
    }


	if (source[src].has_scale_offset[0]) {
		source[src].m_array = poVar->GetUnscaled();
	} else {
		source[src].m_array = poVar;
	}
	source[src].open_read = true;
	return true;
}


bool SpatRaster::readStopMulti(size_t src) {
//	Rcpp::Rcout << "readStopMulti\n";
	source[src].open_read = false;
	source[0].m_array.reset();
	return true;
}



bool SpatRaster::readChunkMulti(std::vector<double> &data, size_t src, size_t row, size_t nrows, size_t col, size_t ncols) {

	std::vector<GUInt64> offset(source[src].m_ndims, 0);
	std::vector<size_t> dims = source[src].m_dims;

/*
	std::vector<size_t> sizes = source[src].m_size;
	Rcpp::Rcout << "dims: ";
	for (size_t i=0; i<dims.size(); i++) {Rcpp::Rcout << dims[i] << " ";}
	Rcpp::Rcout << "\nsize: ";
	for (size_t i=0; i<sizes.size(); i++) {Rcpp::Rcout << sizes[i] << " ";}
	Rcpp::Rcout << "\nrc: ";
	Rcpp::Rcout << col << " " << ncols << " " << row << " " << nrows << "\n";
*/

	offset[source[src].m_dims[0]] = col;
	offset[source[src].m_dims[1]] = row;
	size_t ndim = source[src].m_dims.size();
	std::vector<size_t> count(source[src].m_ndims, 1);
	count[source[src].m_dims[0]] = ncols;
	count[source[src].m_dims[1]] = nrows;

	const size_t rowdim = source[src].m_dims[1];
	std::vector<long long int> stride;
	const long long int *stride_arg = nullptr;
	if (!source[src].flipped) {
		stride.resize(source[src].m_ndims, 1);
		stride[rowdim] = -1;
		offset[rowdim] = nrow() - row - 1;
		stride_arg = stride.data();
	}

	size_t insize = data.size();

	auto dt = GDALExtendedDataType::Create(GDT_Float64);

	const bool md_lat_fast =
		source[src].m_dims.size() >= 2 && source[src].m_dims[1] > source[src].m_dims[0];

	// Two-dimensional variable (lon x lat only): one Read matches terra layout.
	if (ndim == 2) {
		size_t n = vprod(count, false);
		data.resize(insize + n);
		source[src].m_array->Read(&offset[0], &count[0], stride_arg, NULL, dt, &data[insize], NULL, 0);
		if (md_lat_fast) {
			md_reorder_spatial_gdal_to_terra(data, insize, nrows, ncols);
		}
	} else {
		// ndim >= 3: always read one terra-layer at a time. A single Read over all
		// extra dimensions fills the buffer in GDAL's native dimension order; terra
		// stores values layer-major (ncell contiguous per layer), same as
		// readChunkGDAL / RasterIO.
		std::vector<size_t> extra_sizes;
		for (size_t j = 2; j < ndim; j++) {
			size_t gd = source[src].m_dims[j];
			count[gd] = 1;
			extra_sizes.push_back(source[src].m_size[gd]);
		}
		const size_t block = ncols * nrows;
		data.resize(insize + block * source[src].layers.size());
		std::vector<size_t> idx;
		for (size_t i = 0; i < source[src].layers.size(); i++) {
			md_layer_to_indices(source[src].layers[i], extra_sizes, idx);
			for (size_t j = 0; j < extra_sizes.size(); j++) {
				offset[source[src].m_dims[2 + j]] = idx[j];
			}
			double *dest = &data[insize + i * block];
			source[src].m_array->Read(&offset[0], &count[0], stride_arg, NULL, dt, dest, NULL, 0);
			if (md_lat_fast) {
				md_reorder_spatial_gdal_to_terra(data, insize + i * block, nrows, ncols);
			}
		}
	}

	if (source[src].m_hasNA) {
//		Rcpp::Rcout << source[src].m_missing_value << std::endl;
		std::replace (data.begin()+insize, data.end(), source[src].m_missing_value, (double)NAN);
	}

//	data.insert(data.end(), out.begin(), out.end());
	return true;
}

bool SpatRaster::readRowColMulti(size_t src, std::vector<std::vector<double>> &out, size_t outstart, std::vector<int64_t> &rows, const std::vector<int64_t> &cols) {
	
//	Rcpp::Rcout << "readRowColMulti " << src << "\n";
	if (!readStartMulti(src)) {
		return false;
	}
	size_t n = rows.size();
	size_t nl = source[src].layers.size();

	out.resize(outstart + nl);
	for (size_t i = outstart; i < outstart + nl; i++) {
		out[i].clear();
		out[i].reserve(n);
	}

	std::vector<GUInt64> offset(source[src].m_ndims, 0);

	size_t ndim = source[src].m_dims.size();
	std::vector<size_t> count(source[src].m_ndims, 1);

	const size_t rowdim = source[src].m_dims[1];
	std::vector<size_t> extra_sizes;
	for (size_t j = 2; j < ndim; j++) {
		extra_sizes.push_back(source[src].m_size[source[src].m_dims[j]]);
	}
	std::vector<size_t> idx;

	auto dt = GDALExtendedDataType::Create(GDT_Float64);

	// For ndim > 2, read one terra-layer at a time (same as readChunkMulti).
	// The former in_order "fast path" read the full extra-dimensional block into v[0..nl),
	// which assumes GDAL's flattened dimension order matches md_layer_to_indices;
	if (ndim > 2) {
		for (size_t j = 2; j < ndim; j++) {
			count[source[src].m_dims[j]] = 1;
		}
		for (size_t p = 0; p < n; p++) {
			if (std::isnan(cols[p]) || std::isnan(rows[p])) {
				for (size_t j = 0; j < nl; j++) {
					out[outstart + j].push_back(NAN);
				}
				continue;
			}
			offset[source[src].m_dims[0]] = cols[p];
			if (!source[src].flipped) {
				offset[rowdim] = nrow() - rows[p] - 1;
			} else {
				offset[rowdim] = rows[p];
			}
			for (size_t j = 0; j < nl; j++) {
				md_layer_to_indices(source[src].layers[j], extra_sizes, idx);
				for (size_t e = 0; e < extra_sizes.size(); e++) {
					offset[source[src].m_dims[2 + e]] = idx[e];
				}
				double val = NAN;
				source[src].m_array->Read(&offset[0], &count[0], nullptr, NULL, dt, &val, NULL, 0);
				if (source[src].m_hasNA) {
					if (val == source[src].m_missing_value) {
						val = NAN;
					}
				}
				out[outstart + j].push_back(val);
			}
		}
	} else {
		if (source[src].in_order(true)) {
			for (size_t j = 2; j < ndim; j++) {
				size_t gd = source[src].m_dims[j];
				count[gd] = source[src].m_size[gd];
			}
		} else {
			for (size_t j = 2; j < ndim; j++) {
				count[source[src].m_dims[j]] = 1;
			}
		}
		std::vector<double> v(nl, NAN);
		for (size_t p = 0; p < n; p++) {
			if (std::isnan(cols[p]) || std::isnan(rows[p])) {
				for (size_t j = 0; j < nl; j++) {
					out[outstart + j].push_back(NAN);
				}
				continue;
			}
			offset[source[src].m_dims[0]] = cols[p];
			if (!source[src].flipped) {
				offset[rowdim] = nrow() - rows[p] - 1;
			} else {
				offset[rowdim] = rows[p];
			}

			if (source[src].in_order(true)) {
				source[src].m_array->Read(&offset[0], &count[0], nullptr, NULL, dt, &v[0], NULL, 0);
			} else {
				for (size_t j = 0; j < nl; j++) {
					md_layer_to_indices(source[src].layers[j], extra_sizes, idx);
					for (size_t e = 0; e < extra_sizes.size(); e++) {
						offset[source[src].m_dims[2 + e]] = idx[e];
					}
					source[src].m_array->Read(&offset[0], &count[0], NULL, NULL, dt, &v[j], NULL, 0);
				}
			}
			if (source[src].m_hasNA) {
				std::replace(v.begin(), v.end(), source[src].m_missing_value, (double)NAN);
			}
			for (size_t j = 0; j < nl; j++) {
				out[outstart + j].push_back(v[j]);
			}
		}
	}

	readStopMulti(src);	
	return true;
}




bool SpatRaster::writeStartMulti(SpatOptions &opt, const std::vector<std::string> &srcnames) {

	if (!hasValues()) {
		setError("there are no cell values");
		return false;
	}

	std::string filename = opt.get_filename();	
	if (filename.empty()) {
		setError("empty filename");
		return(false);
	} 
	// assure filename won't be used again
	opt.set_filenames({""});
	std::string driver = "netCDF";
/*
	std::string driver = opt.get_filetype();
	getGDALdriver(filename, driver);
	if (driver.empty()) {
		setError("cannot guess file type from filename");
		return(false);
	}
	if (driver != "netCDF") {
		if (driver.empty()) {
			setError("multi-dim only implemented for netCDF");
			return(false);
		}
	}
*/
	
    GDALDriver *poDriver;
    poDriver = GetGDALDriverManager()->GetDriverByName(driver.c_str());

	GDALDataset *poDS;
    poDS = poDriver->CreateMultiDimensional(filename.c_str(), NULL, NULL);
	if (poDS == NULL) {
		setError("failed writing "+ driver + " file");
		GDALClose( (GDALDatasetH) poDS );
		return false;
	}

    auto rg = poDS->GetRootGroup();
	
	std::vector<std::shared_ptr<GDALDimension>> dim_ptrs;
	auto dt = GDALExtendedDataType::Create(GDT_Float64);


	size_t nz = nlyr();
	size_t ny = nrow();
	size_t nx = ncol();
	dim_ptrs.push_back(rg->CreateDimension("Z", "", "", nz));
	dim_ptrs.push_back(rg->CreateDimension("Y", "", "", ny));
	dim_ptrs.push_back(rg->CreateDimension("X", "", "", nx));


    auto var = rg->CreateMDArray("Z", {dim_ptrs[0]}, dt);
    var = rg->OpenMDArray("Z");

	std::vector<double> dvals(nz);
	std::iota(dvals.begin(), dvals.end(), 0);

	std::vector<size_t> count = {nz};
	std::vector<GUInt64> start = {0};

	var->Write(start.data(), count.data(), nullptr, nullptr, dt, &dvals[0]); 

    var = rg->CreateMDArray("Y", {dim_ptrs[1]}, dt);
	yFromRow(dvals);
	count = {ny};
	var->Write(start.data(), count.data(), nullptr, nullptr, dt, &dvals[0]); 

    var = rg->CreateMDArray("X", {dim_ptrs[2]}, dt);
	xFromCol(dvals);
	count = {nx};
	var->Write(start.data(), count.data(), nullptr, nullptr, dt, &dvals[0]); 
	
	std::string vname = source[0].source_name.empty() ? "array" : source[0].source_name;
	
    var = rg->CreateMDArray(vname, dim_ptrs, GDALExtendedDataType::Create(GDT_Float64));
	

	std::string wkt = source[0].srs.wkt;
	if (!wkt.empty()) {
		OGRSpatialReference *srs = NULL;
		srs = new OGRSpatialReference;
		const char *cp = wkt.c_str();
		srs->importFromWkt(cp);
		if (srs != NULL) {
			if (!var->SetSpatialRef(srs)) {
				addWarning("failed to assign CRS to array");
			}
			CPLFree(srs);
		}
	}

	source[0].m_array = var;
	source[0].gdalconnection = poDS;

	return true;
}

bool SpatRaster::writeValuesMulti(std::vector<double> &vals, size_t startrow, size_t nrows, size_t startcol, size_t ncols){
	SpatOptions opt;
	std::vector<GUInt64> start = {0, 0, 0};
	std::vector<size_t> count = {nlyr(), nrow(), ncol()};
	source[0].m_array->Write(start.data(), count.data(), nullptr, nullptr, GDALExtendedDataType::Create(GDT_Float64), &vals[0]); 
	return true;
}

bool SpatRaster::writeStopMulti() {
	//GDALMDArrayRelease(
	source[0].m_array.reset();
	GDALClose( source[0].gdalconnection );
	return true;
}

#else



bool SpatRaster::constructFromFileMulti(std::string fname, std::vector<int> subds, std::vector<std::string> subname, std::vector<std::string> drivers, std::vector<std::string> options, std::vector<int> dims, bool noflip, bool guessCRS, std::vector<std::string> domains) {
	setError("multidim is not supported with GDAL < 3.4 or on 32-bit systems");
	return false;
}

bool SpatRaster::readStartMulti(size_t src) {
	return false;
}

bool SpatRaster::readStopMulti(size_t src) {
	return false;
}

bool SpatRaster::readChunkMulti(std::vector<double> &data, size_t src, size_t row, size_t nrows, size_t col, size_t ncols) {
	return false;
}

bool SpatRaster::readRowColMulti(size_t src, std::vector<std::vector<double>> &out, size_t outstart, std::vector<int64_t> &rows, const std::vector<int64_t> &cols) {
	return false;
}

bool SpatRaster::writeStartMulti(SpatOptions &opt, const std::vector<std::string> &srcnames) {
	return false;
}

bool SpatRaster::writeValuesMulti(std::vector<double> &vals, size_t startrow, size_t nrows, size_t startcol, size_t ncols){
	return false;
}

bool SpatRaster::writeStopMulti() {
	return false;
}

#endif


std::vector<double> SpatRaster::readValuesMulti(size_t src, size_t row, size_t nrows, size_t col, size_t ncols, int lyr) {

	std::vector<double> out;
	if (lyr < 0) {
		if (!readStartMulti(src)) {
			return out;
		}
		readChunkMulti(out, src, row, nrows, col, ncols);
		readStopMulti(src);
		return out;
	}
	return out;
}



void getSampleRowCol2(std::vector<int64_t> &oldrow, std::vector<int64_t> &oldcol, size_t nrows, size_t ncols, size_t snrow, size_t sncol) {

	double rf = nrows / (double)(snrow);
	double cf = ncols / (double)(sncol);
	//double rstart = std::floor(0.5 * rf);
	//double cstart = std::floor(0.5 * cf);
	double rstart = 0.5 * rf;
	double cstart = 0.5 * cf;
	
	std::vector<int64_t> xcol, xrow;
	xcol.reserve(sncol);
	for (size_t i =0; i<sncol; i++) {
        xcol.push_back(i * cf + cstart);
	}
	xrow.reserve(snrow);
	for (size_t i =0; i<snrow; i++) {
        xrow.push_back(i * rf + rstart);
	}
	oldrow.reserve(sncol * snrow);
	oldcol.reserve(sncol * snrow);
	for (size_t i =0; i<snrow; i++) {
		for (size_t j=0; j<sncol; j++) {
			oldrow.push_back(xrow[i]);
			oldcol.push_back(xcol[j]);
		}
	}
}


std::vector<double> SpatRaster::readSampleMulti(size_t src, size_t srows, size_t scols, bool overview) {
	(void) overview;
	std::vector<int64_t> colnr, rownr;
	getSampleRowCol2(rownr, colnr, nrow(), ncol(), srows, scols);
	const size_t n = rownr.size();
	const size_t nl = source[src].layers.size();
	std::vector<std::vector<double>> out(nl);
	if (!readRowColMulti(src, out, 0, rownr, colnr)) {
		return std::vector<double>();
	}
	if (hasError()) {
		return std::vector<double>();
	}
	// Same band layout as readGDALsample / readChunkGDAL: layer-major (cells, then next layer).
	std::vector<double> ret(n * nl);
	for (size_t lyr = 0; lyr < nl; lyr++) {
		if (out[lyr].size() != n) {
			setError("internal error in readSampleMulti: unexpected sample size");
			return std::vector<double>();
		}
		double *dest = ret.data() + lyr * n;
		std::copy(out[lyr].begin(), out[lyr].end(), dest);
	}
	return ret;
}


SpatRaster SpatRaster::writeRasterM(SpatOptions &opt) {
	SpatRaster out;
	
	std::vector<std::string> fnames = opt.get_filenames();

	if (!writeStartMulti(opt, {""})) {
		out.setError(getError());
		return out;
	}
	std::vector<double> vals = getValues(-1, opt);
	if (!writeValuesMulti(vals, 0, nrow(), 0, ncol())) {
		out.setError(getError());
		return out;		
	}
	writeStopMulti();

	std::vector<std::string> empty;
	std::vector<int> dims = {-1};
	
	out.constructFromFileMulti(fnames[0], {0}, empty, empty, empty, dims, false, false, {""});
	return out;
}


std::vector<bool> SpatRaster::is_multidim() {
	std::vector<bool> out;
	out.reserve(source.size());
	for (size_t i=0; i<source.size(); i++) {
		out.push_back(source[i].is_multidim);
	}
	return(out);
}

std::vector<std::vector<std::string>> SpatRaster::dim_names() {
	std::vector<std::vector<std::string>> out(source.size());
	for (size_t i=0; i<source.size(); i++) {
		if (source[i].is_multidim) {
			out[i] = source[i].m_names;
		}
	}
	return(out);
	
}

std::vector<std::vector<size_t>> SpatRaster::dim_order() {
	std::vector<std::vector<size_t>> out(source.size());
	for (size_t i=0; i<source.size(); i++) {
		if (source[i].is_multidim) {
			out[i] = source[i].m_order;
		}
	}
	return out;
}

std::vector<std::vector<size_t>> SpatRaster::dim_size() {
	std::vector<std::vector<size_t>> out(source.size());
	for (size_t i=0; i<source.size(); i++) {
		if (source[i].is_multidim) {
			out[i] = source[i].m_size;
		}
	}
	return out;
}
