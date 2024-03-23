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
#include "string_utils.h"
#include "file_utils.h"
#include "spatTime.h"
#include "recycle.h"
#include "vecmath.h"

#include <set>

#ifdef useGDAL
#include "crs.h"
#endif


SpatRaster::SpatRaster(std::string fname, std::vector<int> subds, std::vector<std::string> subdsname, std::vector<std::string> drivers, std::vector<std::string> options) {
#ifdef useGDAL
	constructFromFile(fname, subds, subdsname, drivers, options);
#endif
}


SpatRaster::SpatRaster(std::vector<std::string> fname, std::vector<int> subds, std::vector<std::string> subdsname, bool multi, std::vector<std::string> drivers, std::vector<std::string> options,std::vector<size_t> xyz) {
// argument "x" is ignored. It is only there to have four arguments such that the  module
// can distinguish this constructor from another with three arguments.
	if (fname.empty()) {
		setError("no filename");
		return;
	}

#ifdef useGDAL
	if (multi) {
		constructFromFileMulti(fname[0], subds, subdsname, drivers, options, xyz);
		return;
	}

	if (!constructFromFile(fname[0], subds, subdsname, drivers, options)) {
		//setError("cannot open file: " + fname[0]);
		return;
	}
	SpatOptions opt;
	for (size_t i=1; i<fname.size(); i++) {
		SpatRaster r;
		bool ok = r.constructFromFile(fname[i], subds, subdsname, drivers, options);
		if (r.msg.has_warning) {
			addWarning(r.msg.warnings[0]);
		}

		if (ok) {
			addSource(r, false, opt);
			if (r.msg.has_error) {
				setError(r.msg.error);
				return;
			}
		} else {
			if (r.msg.has_error) {
				setError(r.msg.error);
			}
			return;
		}
	}
#endif
}


void SpatRaster::setSources(std::vector<SpatRasterSource> &s) {
	source = s;
}


void SpatRaster::setSource(SpatRasterSource &s) {
	s.resize(s.nlyr); // appears to be necessary!
	source = {s};
}


SpatRaster::SpatRaster(SpatRasterSource &s) {
	source = {s};
}


SpatRaster::SpatRaster() {

	SpatRasterSource s;
	s.nrow = 10;
	s.ncol = 10;
	s.extent = SpatExtent();
	s.memory = true;
	s.filename = "";
	//s.driver = "";
	s.nlyr = 1; // or 0?
	s.resize(1);

	s.hasRange = { false };
	s.hasValues = false;
	s.valueType = { 0 };
	s.layers.resize(1, 0);
	s.dtype = "";
	s.names = {"lyr.1"};
	s.srs.proj4 = "+proj=longlat +datum=WGS84";
	s.srs.wkt = "GEOGCRS[\"WGS 84\", DATUM[\"World Geodetic System 1984\", ELLIPSOID[\"WGS 84\",6378137,298.257223563, LENGTHUNIT[\"metre\",1]]], PRIMEM[\"Greenwich\",0, ANGLEUNIT[\"degree\",0.0174532925199433]], CS[ellipsoidal,2], AXIS[\"geodetic latitude (Lat)\",north, ORDER[1], ANGLEUNIT[\"degree\",0.0174532925199433]], AXIS[\"geodetic longitude (Lon)\",east, ORDER[2], ANGLEUNIT[\"degree\",0.0174532925199433]], USAGE[ SCOPE[\"Horizontal component of 3D system.\"], AREA[\"World.\"], BBOX[-90,-180,90,180]], ID[\"EPSG\",4326]]";
	setSource(s);
}

/*
SpatRaster SpatRaster::dropSource() {
	SpatRaster out = geometry();
	out.source.resize(0);
	return out;
}
*/


SpatRaster SpatRaster::subsetSource(size_t snr) {
	if (snr >= source.size()) {
		SpatRaster out;
		out.setError("invalid source number");
		return out;
	}
	SpatRaster out(source[snr]);
	return out;
}

bool SpatRaster::hasValues() {
//	if (source.size() == 0) {
//		return false;
//	} else {
	return source[0].hasValues ;
//	}
}


SpatRaster::SpatRaster(std::vector<unsigned> rcl, std::vector<double> ext, std::string crs) {

	SpatRasterSource s;
	rcl.resize(3, 1);
	s.nrow=rcl[0];
	s.ncol=rcl[1];
	s.extent.xmin = ext[0];
	s.extent.xmax = ext[1];
	s.extent.ymin = ext[2];
	s.extent.ymax = ext[3];
	s.hasValues = false;
	s.hasRange = {false};
	s.valueType = { 0 };

	s.memory = true;
	s.filename = "";
	//s.driver = "";
	s.nlyr = rcl[2];
	s.layers.resize(1, 0);
	//s.layers.resize(1, s.nlyr);
	//std::iota(s.layers.begin(), s.layers.end(), 0);

	s.dtype = "";

#ifdef useGDAL
	std::string msg;
	if (!s.srs.set( crs, msg )) {
		setError(msg);
		return;
	} else if (!msg.empty()) {
		addWarning(msg);
	}

#else
	s.srs.proj4 = lrtrim_copy(crs);
#endif

	for (unsigned i=0; i < rcl[2]; i++) {
		s.names.push_back("lyr." + std::to_string(i+1)) ;
	}

	setSource(s);
}



SpatRaster::SpatRaster(unsigned nr, unsigned nc, unsigned nl, SpatExtent ext, std::string crs) {

	SpatRasterSource s;
	s.nrow = nr;
	s.ncol = nc;
	s.extent = ext;
	s.hasValues = false;
	s.memory = true;
	s.filename = "";
	//s.driver = "";
	s.nlyr = nl;
	s.hasRange = { false };
	s.valueType = { 0 };

	s.layers.resize(1, 0);
	//s.layers.resize(1, _nlyr);
	//std::iota(s.layers.begin(), s.layers.end(), 0);
	s.dtype = "";
#ifdef useGDAL
	std::string msg;
	if (!s.srs.set(crs, msg )) {
		setError(msg);
		return;
	} else if (!msg.empty()) {
		addWarning(msg);
	}
#else
	s.srs.proj4 = lrtrim_copy(crs);
#endif
	for (unsigned i=0; i < nl; i++) {
		s.names.push_back("lyr." + std::to_string(i+1)) ;
	}
	setSource(s);
}


/*
SpatRaster::SpatRaster(const SpatRaster &r) {
	source.nrow = r.nrow;
	source.ncol = r.ncol;
	source.extent = r.extent;
	source.crs = r.crs;
	source.memory = true;
	nlyrs = (nlyrs < 1) ? nlyr(): nlyrs;
	source.resize(nlyrs);
	source.values.resize(0);

	std::vector<std::string> nms(s.nlyr);
	for (size_t i=0; i < s.nlyr; i++) { nms[i] = "lyr" + std::to_string(i+1); }
	source.names = nms;
	// would still need "setSource" to set
}
*/



SpatRaster SpatRaster::geometry(long nlyrs, bool properties, bool time, bool units, bool keeptags) {
	SpatRasterSource s;
	//s.values.resize(0);
	s.nrow = nrow();
	s.ncol = ncol();
	s.extent = getExtent();
	s.srs = source[0].srs;
	//s.prj = prj;
	s.memory = true;
	s.hasValues = false;
	long nl = nlyr();
	bool keepnlyr = ((nlyrs == nl) || (nlyrs < 1));
	nlyrs = (keepnlyr) ? nlyr(): nlyrs;
	
	if (properties) {
		s.hasColors = hasColors();
		s.cols = getColors();
		s.hasCategories = hasCategories();
		s.cats = getCategories();
	}
	s.resize(nlyrs);
	std::vector<std::string> nms;
	if (keepnlyr) {
		nms = getNames();
		if (time && hasTime()) {
			s.hasTime = true;
			s.timestep = getTimeStep();
			s.timezone = getTimeZone();
			s.time = getTime();
		}
		if (units && hasUnit()) {
			s.hasUnit = true;
			s.unit = getUnit();
		}
		
		std::vector<std::string> un = getSourceNames();
		std::sort(un.begin(), un.end() );
		un.erase(std::unique(un.begin(), un.end()), un.end());
		if (un.size() == 1) {
			s.source_name = un[0];
		}
		un = getLongSourceNames();
		std::sort(un.begin(), un.end() );
		un.erase(std::unique(un.begin(), un.end()), un.end());
		if (un.size() == 1) {
			s.source_name_long = un[0];
		}
	} else {
		for (size_t i=0; i < s.nlyr; i++) {
			nms.push_back("lyr" + std::to_string(i+1));
		}
	}
	s.names = nms;
	SpatRaster out(s);
	if (properties) {
		out.rgb = rgb;
		out.rgbtype = rgbtype;
		out.rgblyrs = rgblyrs;
	}	
	if (keeptags) {
		out.tags = tags;
		out.lyrTags = lyrTags;
	}
	return out;
}


SpatRaster SpatRaster::geometry_opt(long nlyrs, bool properties, bool time, bool units, bool tags, bool datatype, SpatOptions &opt) {

	if (datatype && hasValues() && (!opt.datatype_set)) {
		std::vector<std::string> dt = getDataType(true, true);
		if ((dt.size() == 1) && !dt[0].empty()) {
			if (!hasScaleOffset()) {
				opt.set_datatype(dt[0]);
			}
		}
	}	
	
	return geometry(nlyrs, properties, time, units, tags);
}

SpatRaster SpatRaster::deepCopy() {
	return *this;
}



std::vector<double> SpatRaster::resolution() {
	SpatExtent extent = getExtent();
	return std::vector<double> { (extent.xmax - extent.xmin) / ncol(), (extent.ymax - extent.ymin) / nrow() };
}


SpatRaster SpatRaster::setResolution(double xres, double yres) {
	SpatRaster out;

	if ((xres <= 0) | (yres <= 0)) {
		out.setError("resolution must be larger than 0");
		return(out);
	}
	SpatExtent e = getExtent();
	unsigned nc = std::max(1., round((e.xmax-e.xmin) / xres));
	unsigned nr = std::max(1., round((e.ymax-e.ymin) / yres));

	double xmax = e.xmin + nc * xres;
	double ymax = e.ymin + nr * yres;
	unsigned nl = nlyr();
	std::vector<unsigned> rcl = {nr, nc, nl};
	std::vector<double> ext = {e.xmin, xmax, e.ymin, ymax};

	out = SpatRaster(rcl, ext, {""});
	out.source[0].srs = source[0].srs;
	return out;
}


size_t SpatRaster::ncol() {
	if (source.empty()) {
		return 0;
	} else {
		return source[0].ncol;
	}
}

size_t SpatRaster::nrow() {
	if (source.empty()) {
		return source[0].ncol;
	} else {
		return source[0].nrow;
	}
}


unsigned SpatRaster::nlyr() {
	unsigned x = 0;
	for (size_t i=0; i<source.size(); i++) {
		x += source[i].nlyr;
	}
	return(x);
}

std::vector<std::string> SpatRaster::filenames() {
	std::vector<std::string> x(source.size());
	for (size_t i=0; i<x.size(); i++) { x[i] = source[i].filename; }
	return(x);
}

std::vector<bool> SpatRaster::inMemory() {
	std::vector<bool> m(source.size());
	for (size_t i=0; i<m.size(); i++) { m[i] = source[i].memory; }
	return(m);
}

std::vector<bool> SpatRaster::hasRange() {
	std::vector<bool> x;
	for (size_t i=0; i<source.size(); i++) {
		x.insert(x.end(), source[i].hasRange.begin(), source[i].hasRange.end());
	}
	return(x);
}

std::vector<int> SpatRaster::getValueType(bool unique) {
	std::vector<int> d;
	d.reserve(nlyr());
	for (size_t i=0; i<source.size(); i++) {
		d.insert(d.end(), source[i].valueType.begin(), source[i].valueType.end());
	}
	if (unique) {
		std::sort(d.begin(), d.end());
		d.erase(std::unique(d.begin(), d.end()), d.end());	
	}
	return(d);
}


bool SpatRaster::setValueType(unsigned char d) {
	if (d > 3) {
		return false;
	}
	for (size_t i=0; i<source.size();i++) {
		source[i].valueType = std::vector<unsigned char>(source[i].nlyr, d);
	}
	return true;
}


std::vector<double> SpatRaster::range_min() {
	std::vector<double> x;
	for (size_t i=0; i<source.size(); i++) {
		x.insert(x.end(), source[i].range_min.begin(), source[i].range_min.end());
	}
	return(x);
}

std::vector<double> SpatRaster::range_max() {
	std::vector<double> x;
	for (size_t i=0; i<source.size(); i++) {
		x.insert(x.end(), source[i].range_max.begin(), source[i].range_max.end());
	}
	return(x);
}


bool SpatRaster::is_lonlat() {
	if (source[0].srs.is_lonlat()) {
		SpatExtent e = getExtent();
		if ((e.xmin < -181) || (e.xmax > 361) || (e.ymin < -90.001) || (e.ymax > 90.001)) {
			addWarning("coordinates are out of range for lon/lat");
		}
		return true;
	}
	return false;
}

bool SpatRaster::could_be_lonlat() {
	if (is_lonlat()) return true;
	SpatExtent e = getExtent();
	return source[0].srs.could_be_lonlat(e);
}


bool SpatRaster::is_global_lonlat() {
	SpatExtent e = getExtent();
	return source[0].srs.is_global_lonlat(e);
}

int SpatRaster::ns_polar() {
	int polar = 0;
	if (!is_lonlat()) {
		return polar;
	}
	SpatExtent e = getExtent();
	if ((e.ymax - 90) > -0.00001) {
		polar = 1;
	}
	if ((e.ymin + 90) < 0.00001) {
		polar = polar == 1 ? 2 : -1;
	}
	return polar;
}



bool SpatRaster::sources_from_file() {
	for (size_t i=0; i<source.size(); i++) {
		if (!source[i].memory) {
			return true;
		}
	}
	return false;
}


SpatRaster SpatRaster::sources_to_disk(std::vector<std::string> &tmpfs, bool unique, SpatOptions &opt) {
// if a tool needs to read from disk, perhaps from unique filenames
// use writeRaster to write to a single file.
	SpatRaster out;
	size_t nsrc = source.size();
	std::set<std::string> ufs;
	size_t ufsize = ufs.size();

	std::string tmpbasename = tempFile(opt.get_tempdir(), opt.tmpfile, "_temp_");

	SpatOptions ops(opt);
	for (size_t i=0; i<nsrc; i++) {
		bool write = false;
		if (!source[i].in_order() || source[i].memory) {
			write = true;
		} else if (unique) {
			ufs.insert(source[i].filename);
			if (ufs.size() == ufsize) {
				write = true;
			} else {
				ufsize++;
			}
		}
		SpatRaster rs(source[i]);
		if (write) {
			std::string fname = tmpbasename + std::to_string(i) + ".tif";
			opt.set_filenames({fname});
			tmpfs.push_back(fname);
			rs = rs.writeRaster(opt);
		}
		if (i == 0) {
			out.setSource(rs.source[0]);
		} else {
			out.addSource(rs, false, ops);
		}
	}
	return out;
}

bool SpatRaster::setSRS(std::string crs) {
	std::string msg;
	SpatSRS srs;
	if (!srs.set(crs, msg )) {
		addWarning("Cannot set raster SRS: "+ msg);
		return false;
	} else if (!msg.empty()) {
		addWarning(msg);
	}

	for (size_t i = 0; i < nsrc(); i++) {
		source[i].srs = srs;
		if (!source[i].memory) {
			source[i].parameters_changed = true;
		}
	}
	return true;
}


/*
#ifdef useGDAL
bool SpatRaster::setSRS(OGRSpatialReference *poSRS, std::string &msg) {
	if (!srs.set(poSRS, msg)){
		addWarning("Cannot set raster SRS: "+ msg);
		return false;
	}
	for (size_t i = 0; i < nsrc(); i++) {
		source[i].srs = srs;
	}
	return true;
}
#endif
*/

std::string  SpatRaster::getSRS(std::string x) {
	return source[0].srs.get(x);
}



std::vector<std::string> SpatRaster::getNames() {
	std::vector<std::string> x;
	for (size_t i=0; i<source.size(); i++) {
		x.insert(x.end(), source[i].names.begin(), source[i].names.end());
	}
	return(x);
}


bool SpatRaster::setNames(std::vector<std::string> names, bool make_valid) {
	if (names.size() == 1) {
		recycle(names, nlyr());
	}

	if (names.size() != nlyr()) {
		return false;
	} else {
		if (make_valid) {
			make_valid_names(names);
			make_unique_names(names);
        }
		size_t begin=0;
        size_t end;
        for (size_t i=0; i<source.size(); i++)	{
            end = begin + source[i].nlyr;
            source[i].names = std::vector<std::string> (names.begin() + begin, names.begin() + end);
            begin = end;
        }
        return true;
	}
}


std::vector<std::string> SpatRaster::getLongSourceNames() {
	std::vector<std::string> x;
	x.reserve(source.size());
	for (size_t i=0; i<source.size(); i++) {
		x.push_back(source[i].source_name_long);
	}
	return(x);
}


bool SpatRaster::setLongSourceNames(std::vector<std::string> names) {
	if (names.size() == 1) {
		for (size_t i=0; i<source.size(); i++)	{
			source[i].source_name_long = names[0];
		}
	} else if (names.size() == nsrc()) {
		for (size_t i=0; i<source.size(); i++)	{
			source[i].source_name_long = names[i];
		}
	} else {
		return false;
	}
	return true;
}



std::vector<std::string> SpatRaster::getSourceNames() {
	std::vector<std::string> x;
	x.reserve(source.size());
	for (size_t i=0; i<source.size(); i++) {
		x.push_back(source[i].source_name);
	}
	return(x);
}


bool SpatRaster::setSourceNames(std::vector<std::string> names) {
	if (names.size() == 1) {
		for (size_t i=0; i<source.size(); i++)	{
			source[i].source_name = names[0];
		}
	} else if (names.size() == nsrc()) {
		for (size_t i=0; i<source.size(); i++)	{
			source[i].source_name = names[i];
		}
	} else {
		return false;
	}
	return true;
}


bool SpatRaster::setNAflag(std::vector<double> flag) {
	size_t sz = source.size();
	if (flag.size() == 1) recycle(flag, sz);
	if (flag.size() != sz) {
		return false;
	}
	double na = NAN;
	for (size_t i=0; i<sz; i++)	{
		if (std::isnan(flag[i])) {
			source[i].hasNAflag = false;
			source[i].NAflag = NAN;
		} else {

			if (source[i].memory) {
				source[i].hasNAflag = false;
				std::replace(source[i].values.begin(), source[i].values.end(), flag[i], na);
				source[i].setRange();
			} else {
				source[i].hasNAflag = true;
				source[i].NAflag = flag[i];
			}
		}
	}
	return true;
}


std::vector<double> SpatRaster::getNAflag() {
	std::vector<double> out(source.size(), NAN);
	for (size_t i=0; i<source.size(); i++)	{
		if (source[i].hasNAflag) {
			out[i] = source[i].NAflag;
		}
	}
	return out;
}


bool SpatRaster::hasTime() {
	bool test = source[0].hasTime;
	for (size_t i=1; i<source.size(); i++) {
		test = test && source[i].hasTime;
	}
	return test;
}

/*
std::vector<double> SpatRaster::getTimeDbl() {
	std::vector<int_64> t64 = getTime();
	std::vector<double> out(t64.size());
	for (size_t i=0; i < out.size(); i++) {
		out[i] = t64[i];
	}
	return out;
}
*/

std::vector<std::string> SpatRaster::getTimeStr(bool addstep) {
	std::vector<std::string> out;
	std::vector<int_64> time = getTime();
	out.reserve(time.size()+addstep);
	if (addstep) out.push_back(source[0].timestep);
	if (source[0].timestep == "seconds") {
		for (size_t i=0; i < time.size(); i++) {
			std::vector<int> x = get_date(time[i]);
			if (x.size() > 2) {
				out.push_back( std::to_string(x[0]) + "-"
						  + std::to_string(x[1]) + "-"
						  + std::to_string(x[2]) + " "
						  + std::to_string(x[3]) + ":"
						  + std::to_string(x[4]) + ":"
						  + std::to_string(x[5]) );

			} else {
				out.push_back("");
			}
		}
	} else if (source[0].timestep == "days") {
		for (size_t i=0; i < time.size(); i++) {
			std::vector<int> x = get_date(time[i]);
			if (x.size() > 2) {
				out.push_back( std::to_string(x[0]) + "-"
						  + std::to_string(x[1]) + "-"
						  + std::to_string(x[2]) );

			} else {
				out.push_back("");
			}
		}
	} else {
		for (size_t i=0; i < time.size(); i++) {
			out.push_back( std::to_string(time[i]));
		}
	}
	return out;
}


std::vector<int_64> SpatRaster::getTime() {
	std::vector<int_64> x;
	for (size_t i=0; i<source.size(); i++) {
		if (source[i].time.size() != source[i].nlyr) {
			std::vector<double> nas(source[i].nlyr, 0);
			x.insert(x.end(), nas.begin(), nas.end());
		} else {
			x.insert(x.end(), source[i].time.begin(), source[i].time.end());
		}
	}
	return(x);
}

std::string SpatRaster::getTimeStep() {
	return source[0].timestep;
}

std::string SpatRaster::getTimeZone() {
	return source[0].timezone;
}

bool SpatRaster::setTime(std::vector<int_64> time, std::string step, std::string zone) {

	if (time.empty() || step == "remove") {
		for (size_t i=0; i<source.size(); i++)	{
			source[i].time = std::vector<int_64> (source[i].nlyr);
			source[i].timestep = "";
			source[i].timezone = "";
			source[i].hasTime = false;
		}
		return true;
	}

	if (time.size() != nlyr()) {
		return false;
	}
	
	std::vector<std::string> steps = {"seconds", "raw", "days", "yearmonths", "years", "months"};
	if (!is_in_vector(step, steps)) {
		return false;
	}
	size_t begin=0;
	for (size_t i=0; i<source.size(); i++)	{
		size_t end = begin + source[i].nlyr;
        source[i].time = std::vector<int_64> (time.begin() + begin, time.begin() + end);
		source[i].timestep = step;
		source[i].timezone = zone;
		source[i].hasTime = true;
        begin = end;
    }

    return true;
}


std::vector<double> SpatRaster::getDepth() {
	std::vector<double> x;
	for (size_t i=0; i<source.size(); i++) {
		if (source[i].depth.size() != source[i].nlyr) {
			std::vector<double> nas(source[i].nlyr, NAN);
			x.insert(x.end(), nas.begin(), nas.end());
		} else {
			x.insert(x.end(), source[i].depth.begin(), source[i].depth.end());
		}
	}
	return(x);
}



bool SpatRaster::setDepth(std::vector<double> depths) {

	if (depths.empty()) {
		for (size_t i=0; i<source.size(); i++)	{
			source[i].depth = std::vector<double>(source[i].nlyr);
		}
		return true;
	}

	if (depths.size() == 1) {
        for (size_t i=0; i<source.size(); i++)	{
            source[i].depth = std::vector<double> (source[i].nlyr, depths[0]);
        }
        return true;
	} else if (depths.size() != nlyr()) {
		return false;
	} else {
        size_t begin=0;
        for (size_t i=0; i<source.size(); i++)	{
            size_t end = begin + source[i].nlyr;
            source[i].depth = std::vector<double> (depths.begin() + begin, depths.begin() + end);
            begin = end;
        }
        return true;
	}
}



bool SpatRaster::setUnit(std::vector<std::string> units) {
	if (units.size() == 1) {
		bool hu = true;
		if (units[0].empty()) {
			hu = false;
		}
        for (size_t i=0; i<source.size(); i++)	{
            source[i].unit = std::vector<std::string> (source[i].nlyr, units[0]);
			source[i].hasUnit = hu;
        }
        return true;
	} else if (units.size() != nlyr()) {
		return false;
	} else {
        size_t begin=0;
        for (size_t i=0; i<source.size(); i++)	{
            size_t end = begin + source[i].nlyr;
            source[i].unit = std::vector<std::string> (units.begin() + begin, units.begin() + end);
            source[i].hasUnit = true;
			begin = end;
        }
        return true;
	}
}

bool SpatRaster::hasUnit() {
	bool test = source[0].hasUnit;
	for (size_t i=1; i<source.size(); i++) {
		test = test && source[i].hasUnit;
	}
	return test;
}


std::vector<std::string> SpatRaster::getUnit() {
	std::vector<std::string> x;
	for (size_t i=0; i<source.size(); i++) {
		if (source[i].unit.size() != source[i].nlyr) {
			std::vector<std::string> nas(source[i].nlyr, "");
			x.insert(x.end(), nas.begin(), nas.end());
		} else {
			x.insert(x.end(), source[i].unit.begin(), source[i].unit.end());
		}
	}
	return(x);
}




double SpatRaster::xres() {
	SpatExtent extent = getExtent();
	return (extent.xmax - extent.xmin) / ncol() ;
}

double SpatRaster::yres() {
	SpatExtent extent = getExtent();
	return (extent.ymax - extent.ymin) / nrow() ;
}


std::vector<bool> SpatRaster::is_rotated() {
	std::vector<bool> b(source.size(), false);
	for (size_t i=0; i<source.size(); i++) {
		if (source[i].rotated) {
			b[i] = true;
		}
	}
	return b;
}


bool SpatRaster::valid_sources(bool files, bool rotated) {
	std::vector<std::string> ff;
	for (size_t i=0; i<source.size(); i++) {
		std::string f = source[i].filename;
		if (f.empty()) continue;
		if (files) {
			std::size_t found = f.find(':'); // perhaps http: or PG:xxx
			if ((found == 1) || (found == std::string::npos)) {
				if (!file_exists(f)) {
					setError("missing source: " + f);
					return false;
				}
			}
		}
		if (rotated) {
			if (source[i].rotated) {
				setError(f + " is rotated");
				return false;
			}
		}
	}
	return true;
}

std::vector<bool> SpatRaster::hasWindow() {
	std::vector<bool> out;
	out.reserve(nlyr());
	for (size_t i=0; i<nsrc(); i++) {
		for (size_t j=0; j<source[i].nlyr; j++) {
			out.push_back(source[i].hasWindow);
		}
	}
	return out;
}


bool SpatRaster::removeWindow() {
	for (size_t i=0; i<nsrc(); i++) {
		if (source[i].hasWindow) {
			SpatExtent e = source[0].window.full_extent;
			setExtent(e, true, true, "");
			for (size_t i=0; i<source.size(); i++) {
				source[i].hasWindow = false;
				source[i].nrow = source[0].window.full_nrow;
				source[i].ncol = source[0].window.full_ncol;
			}
		}
	}
	return true;
}


bool SpatRaster::setWindow(SpatExtent x) {

	if ( !x.valid() ) {
		setError("invalid extent");
		return false;
	}

	removeWindow();
	x = align(x, "near");
	SpatExtent e = getExtent();
	if (x.compare(e, "==", 0.1 * xres())) {
		return true;
	}

	e = e.intersect(x);
	if ( !e.valid() ) {
		setError("extents do not overlap");
		return false;
	}

// get read-window
	double xr = xres();
	double yr = yres();

	bool expand = false;
	std::vector<size_t> rc(2);
	std::vector<size_t> exp(4, 0);

	int_64 r = rowFromY(x.ymax - 0.5 * yr);
	if (r < 0) {
		rc[0] = 0;
		expand = true;
		exp[0] = trunc(abs(e.ymax - x.ymax) / yr);
	} else {
		rc[0] = r;
	}
	r = rowFromY(x.ymin + 0.5 * yr);
	if (r < 0) {
		expand = true;
		exp[1] = trunc((e.ymax - x.ymin) / yr);
	}

	r = colFromX(x.xmin + 0.5 * xr);
	if (r < 0) {
		rc[1] = 0;
		expand = true;
		exp[2] = trunc((x.xmin - e.xmin) / xres());
	} else {
		rc[1] = r;
	}
	r = colFromX(x.xmax - 0.5 * xr);
	if (r < 0) {
		expand = true;
		exp[3] = trunc(abs(x.xmin - e.xmin) / xres());
	}

	if (expand) {
		setError("expansion is not yet allowed");
		return false;
	}

	for (size_t i=0; i<source.size(); i++) {
		source[i].window.off_row = rc[0];
		source[i].window.off_col = rc[1];
		source[i].window.expand = exp;
		source[i].window.expanded  = expand;
		source[i].window.full_extent = getExtent();
		source[i].window.full_nrow = source[i].nrow;
		source[i].window.full_ncol = source[i].ncol;
		source[i].hasWindow = true;
	}
	setExtent(x, true, true, "");

	return true;
}

SpatRaster SpatRaster::replace(SpatRaster x, unsigned layer, SpatOptions &opt) {

	SpatRaster out = geometry();
	if (!out.compare_geom(x, false, true, opt.get_tolerance())) {
		return(out);
	}
	SpatOptions fopt(opt);

	size_t n = nlyr();
	if (n == 1) {
		return x;
	}
	std::vector<unsigned> lyrs;
	if (layer == 0) {
		out = x;
		lyrs.resize(n-1);
		std::iota(lyrs.begin(), lyrs.end(), 1);
		SpatRaster r = subset(lyrs, fopt);
		out.addSource(r, false, fopt);
	} else if (layer == n-1) {
		lyrs.resize(n-1);
		std::iota(lyrs.begin(), lyrs.end(), 0);
		out = subset(lyrs, fopt);
		out.addSource(x, false, fopt);
	} else {
		lyrs.resize(layer);
		std::iota(lyrs.begin(), lyrs.end(), 0);
		out = subset(lyrs, fopt);
		out.addSource(x, false, fopt);
		lyrs.resize(n-layer-1);
		std::iota(lyrs.begin(), lyrs.end(), layer+1);
		SpatRaster r = subset(lyrs, fopt);
		out.addSource(r, false, fopt);
	}
	return out;
}


SpatRaster SpatRaster::makeCategorical(long layer, SpatOptions &opt) {

	SpatRaster out;
	if (!hasValues()) {
		out.setError("cannot make categories if the raster has no values");
		return out;
	}

	SpatRaster r;
	SpatOptions fopt(opt);
	if (layer >= 0) {
		if (layer > (long) nlyr()) {
			out.setError("layer number is too high");
			return out;
		}
		std::vector<unsigned> lyrs = {(unsigned) layer};
		r = subset(lyrs, fopt);
	} else {
		r = *this;
	}	
	r.math2("round", 0, fopt);
	std::vector<std::vector<double>> u = r.unique(true, NAN, true, fopt);
	std::vector<std::string> names = r.getNames();
	
	for (size_t i=0; i<r.nlyr(); i++) { 
		std::vector<long> uu(u[i].size());
		std::vector<std::string> s(u[i].size());
		for (size_t j=0; j<s.size(); j++) {
			uu[j] = (long)u[i][j];
			s[j] = std::to_string(uu[j]);
		}
		r.setLabels(i, uu, s, names[i]);
	}
	
	if (nlyr() == r.nlyr()) {
		return r;
	} else {
		return replace(r, layer, opt);
	}
}



bool SpatRaster::createCategories(unsigned layer, SpatOptions &opt) {
	if (layer > (nlyr()-1)) {
		setError("invalid layer number");
		return(false);
	}
	std::vector<unsigned> lyrs(1, layer);
	SpatRaster r = subset(lyrs, opt);
	std::vector<std::vector<double>> u = r.unique(false, NAN, true, opt);
    std::vector<unsigned> sl = findLyr(layer);

	std::vector<std::string> s(u[0].size());
	for (size_t i=0; i<s.size(); i++) {
		s[i] = std::to_string(i+1);
	}
	s.resize(256);
	//std::transform(u[0].begin(), u[0].end(), s.begin(), [](const double& d) {
	//	return std::to_string(d);
	//});
	SpatCategories cat;
	cat.d.add_column(s, "category");
	cat.index = 0;
	source[sl[0]].cats[sl[1]] = cat;
	return true;
}


std::vector<bool> SpatRaster::hasCategories() {
	std::vector<bool> b;
	b.reserve(nlyr());
	std::vector<unsigned> ns = nlyrBySource();
	for (size_t i=0; i<ns.size(); i++) {
		for (size_t j=0; j<ns[i]; j++) {
			b.push_back(source[i].hasCategories[j]);
		}
	}
	return b;
}

std::vector<std::string> SpatRaster::getDataType(bool unique, bool memtype) {
	std::vector<std::string> d;
	size_t n = nsrc();
	d.reserve(n);
	for (size_t i=0; i<n; i++) {
		if (memtype && source[i].memory) {			
			std::vector<unsigned char> v = source[i].valueType;
			std::sort(v.begin(), v.end());
			v.erase(std::unique(v.begin(), v.end()), v.end());	
			if (v.size() == 1) {
				if (v[0] == 1) {
					if (vmax(source[i].range_min, false) > 0) {
						d.push_back("INT4U");						
					} else {
						d.push_back("INT4S");	
					}
				} else if (v[0] == 3) {
					d.push_back("INT1U");	
				}
			} else {
				d.push_back("FLT4S");
			}
		} else {
			d.push_back(source[i].dtype);
		}
	}
	if (unique) {
		std::sort(d.begin(), d.end());
		d.erase(std::unique(d.begin(), d.end()), d.end());
	}
	return d;
}

std::vector<std::string> SpatRaster::dataType() {
	std::vector<std::string> d;
	size_t n = nsrc();
	d.reserve(n);
	for (size_t i=0; i<n; i++) {
		d.push_back(source[i].dtype);
	}
	return d;
}


std::vector<std::vector<std::string>> SpatRaster::getMetadata(bool layers) {
	std::vector<std::vector<std::string>> d;
	size_t n = nsrc();
	if (layers) {
		d.reserve(nlyr());
		for (size_t i=0; i<n; i++) {
			if (source[i].bmdata.empty()) {
				d.resize(d.size() + source[i].nlyr);
			} else {
				d.insert(d.end(), source[i].bmdata.begin(), source[i].bmdata.end());
			}
		}
	} else {
		d.resize(n);
		for (size_t i=0; i<n; i++) {
			if (!source[i].smdata.empty()) {
				d[i] = source[i].smdata;
			}
		}
	}
	return d;
}



bool SpatRaster::setLabels(unsigned layer, std::vector<long> values, std::vector<std::string> labels, std::string name) {

	if (layer > (nlyr()-1)) {
		setError("invalid layer number");
		return(false);
	}
	if (values.size() != labels.size()) {
		setError("size of values is not the same as the size of labels");
		return(false);
	}
	if (values.empty()) {
		addWarning("no labels");
		return(true);
	}

    std::vector<unsigned> sl = findLyr(layer);

	SpatCategories cats;
	cats.d.add_column(values, "ID");
	cats.d.add_column(labels, name);
	cats.index = 1;

	if (source[sl[0]].cats.size() <= sl[1]) {
		source[sl[0]].cats.resize(sl[1]+1);
		source[sl[0]].hasCategories.resize(sl[1]+1);
	}
	source[sl[0]].cats[sl[1]] = cats;
	source[sl[0]].hasCategories[sl[1]] = true;
	return true;
}




bool SpatRaster::setCategories(unsigned layer, SpatDataFrame d, unsigned index) {

	if (layer >= nlyr()) {
		setError("invalid layer number");
		return(false);
	}
    std::vector<unsigned> sl = findLyr(layer);

	SpatCategories cats;
	cats.d = d;
	cats.index = index;

	if (source[sl[0]].cats.size() < sl[1]) {
		source[sl[0]].cats.resize(sl[1]);
	}
	source[sl[0]].cats[sl[1]] = cats;
	source[sl[0]].hasCategories[sl[1]] = true;
	return true;
}


bool SpatRaster::removeCategories(long layer) {
	if (layer > (((long)nlyr())-1)) {
		setError("invalid layer number");
		return(false);
	}
	SpatCategories s;
	if (layer < 0) {
		for (size_t i=0; i<source.size(); i++) {
			for (size_t j=0; j<source[i].cats.size(); j++) {
				source[i].cats[j] = s;
				source[i].hasCategories[j] = false;
			}
		}
	} else {
		std::vector<unsigned> sl = findLyr(layer);
		source[sl[0]].cats[sl[1]] = s;
		source[sl[0]].hasCategories[sl[1]] = false;
	}
	return true;
}

SpatCategories SpatRaster::getLayerCategories(unsigned layer) {
    std::vector<unsigned> sl = findLyr(layer);
	SpatCategories cat = source[sl[0]].cats[sl[1]];
	return cat;
}

std::vector<SpatCategories> SpatRaster::getCategories() {
	std::vector<SpatCategories> cats;
	cats.reserve(nlyr());
	for (size_t i=0; i<source.size(); i++) {
		cats.insert(cats.end(), source[i].cats.begin(), source[i].cats.end());
	}
	return cats;
}


std::vector<std::vector<double>> SpatRaster::getScaleOffset() {
	std::vector<std::vector<double>> so(2);
	so[0].reserve(nlyr());
	so[1].reserve(nlyr());
	for (size_t i=0; i<source.size(); i++) {
		so[0].insert(so[0].end(), source[i].scale.begin(), source[i].scale.end());
		so[1].insert(so[1].end(), source[i].offset.begin(), source[i].offset.end());
	}
	return so;
}

bool SpatRaster::hasScaleOffset() {
	for (size_t i=0; i<source.size(); i++) {
		for (size_t j=0; j<source[i].has_scale_offset.size(); j++) {
			if (source[i].has_scale_offset[j]) return true;
		}
	}
	return false;
}

bool SpatRaster::setScaleOffset(std::vector<double> sc, std::vector<double> of) {
	size_t n = sc.size();
	size_t nl = nlyr();
	if (n != of.size()) {
		setError("length of scale and offset should be the same");
		return false;
	}
	if (n > nl) {
		setError("length of scale and offset cannot exceed the number of layers");
		return false;
	}
	if (n < nl) {
		recycle(sc, nl);
		recycle(of, nl);
		if (n > 1) {
			addWarning("recycling scale and offset to the number of layers");
		}
	}
	size_t k=0;
	size_t nc=ncell();
	for (size_t i=0; i<source.size(); i++)	{
		if (source[i].memory) {
			for (size_t j=0; j<source[i].nlyr; j++) {
				size_t loff = j * nc;
				if ((sc[k] != 1) || (of[k] != 0)) {
					for (size_t p=loff; p<(loff+nc); p++) {
						source[i].values[p] = source[i].values[p] * sc[k] + of[k];
					}
					source[i].range_min[j] = source[i].range_min[j] * sc[k] + of[k];
					source[i].range_max[j] = source[i].range_max[j] * sc[k] + of[k];
				}
				k++;
			}
		} else {
			for (size_t j=0; j<source[i].nlyr; j++) {
				if (source[i].has_scale_offset[j]) {
					source[i].range_min[j] = (source[i].range_min[j] - source[i].offset[j]) / source[i].scale[j];
					source[i].range_max[j] = (source[i].range_max[j] - source[i].offset[j]) / source[i].scale[j];
				}
				source[i].scale[j] = sc[k];
				source[i].offset[j] = of[k];
				if ((sc[k] != 1) || (of[k] != 0)) {
					source[i].has_scale_offset[j] = true;
					source[i].range_min[j] = source[i].range_min[j] * sc[k] + of[k];
					source[i].range_max[j] = source[i].range_max[j] * sc[k] + of[k];
				} else {
					source[i].has_scale_offset[j] = false;
				}
				k++;
			}
		}
	}
	return true;
}


std::vector<std::string> SpatRaster::getLabels(unsigned layer) {
	std::vector<std::string> out;
	if (layer >= nlyr()) return out;

	std::vector<bool> hascat = hasCategories();
	if (!hascat[layer]) return out;

	std::vector<SpatCategories> cats = getCategories();
	SpatCategories cat = cats[layer];

	int nc = cat.d.ncol();
	if (nc <= 0) return out;

	cat.index = cat.index > (nc-1) ? (nc-1) : cat.index;
	out = cat.d.as_string(cat.index);
	return out;
}

bool SpatRaster::setCatIndex(unsigned layer, int index) {
	if (layer > (nlyr()-1)) {
		return(false);
	}
    std::vector<unsigned> sl = findLyr(layer);
	int nc = source[sl[0]].cats[sl[1]].d.ncol();
	if (index < nc) {
		source[sl[0]].cats[sl[1]].index = index;
		if (index >= 0) {
			source[sl[0]].names[sl[1]] = source[sl[0]].cats[sl[1]].d.names[index];
		}
		return true;
	} else {
		return false;
	}
}

int SpatRaster::getCatIndex(unsigned layer) {
	if (layer > (nlyr()-1)) {
		return( -1 );
	}
    std::vector<unsigned> sl = findLyr(layer);
	return source[sl[0]].cats[sl[1]].index;
}

SpatRaster SpatRaster::dropLevels() {
	std::vector<bool> hascats = hasCategories();
	bool bany = false;
	for (size_t i=0; i<hascats.size(); i++) {
		if (hascats[i]) {
			bany = true;
			break;
		}
	}
	if (!bany) return *this;

	std::vector<SpatCategories> cats = getCategories();
	SpatOptions opt;
	SpatRaster out = *this;
	std::vector<std::vector<double>> uvv = unique(true, NAN, true, opt);
	for (size_t i=0; i<hascats.size(); i++) {
		if (hascats[i]) {
			SpatCategories lyrcats = cats[i];
			size_t n = lyrcats.d.nrow();
			std::vector<double> uv = uvv[i];
			std::vector<long> uvi(uv.size());
			for (size_t j=0; j<uv.size(); j++) {
				uvi[j] = uv[j];
			}
			std::vector<long> isin;
			isin.reserve(n);
			for (size_t j=0; j<n; j++) {
				for (size_t k=0; k<uvi.size(); k++) {
					if (lyrcats.d.iv[0][j] == uvi[k]) {
						isin.push_back(j);
						continue;
					}
				}
			}
			lyrcats.d = lyrcats.d.subset_rows(isin);
			if (lyrcats.d.nrow() == 0) {
				out.removeCategories(i);
			} else {
				out.setCategories(i, lyrcats.d, lyrcats.index);
			}
		}
	}
	return out;
}



std::vector<SpatDataFrame> SpatRaster::getColors() {
	std::vector<SpatDataFrame> cols;
	for (size_t i=0; i<source.size(); i++) {
		cols.insert(cols.end(), source[i].cols.begin(), source[i].cols.end());
	}
	return cols;
}



bool SpatRaster::setColors(size_t layer, SpatDataFrame cols) {
	if (cols.ncol() < 4 || cols.ncol() > 5) {
		setError("n columns should be 4 or 5");
		return false;
	}
	if (layer >= nlyr()) {
		setError("layer > nlyr");
		return false;
	}
	if (cols.ncol() == 4) {
		std::vector<long> a(cols.nrow(), 255);
		cols.add_column(a, "alpha");
	}

    std::vector<unsigned> sl = findLyr(layer);
	if (source[sl[0]].cols.size() < (sl[1]+1)) {
		source[sl[0]].cols.resize(sl[1]+1);
	}
	if (source[sl[0]].hasColors.size() < (sl[1]+1)) {
		source[sl[0]].hasColors.resize(sl[1]+1);
	}

	source[sl[0]].cols[sl[1]] = cols;
	source[sl[0]].hasColors[sl[1]] = (cols.nrow() > 0);
	return true;
}


bool SpatRaster::removeColors(size_t layer) {
	if (layer >= nlyr()) {
		return false;
	}
    std::vector<unsigned> sl = findLyr(layer);
	if (source[sl[0]].hasColors[sl[1]]) {
		SpatDataFrame d;
		source[sl[0]].cols[sl[1]] = d;
		source[sl[0]].hasColors[sl[1]] = false;
	}
	return true;
}



std::vector<bool> SpatRaster::hasColors() {
	std::vector<bool> b(nlyr());
	std::vector<unsigned> ns = nlyrBySource();
	unsigned k = 0;
	for (size_t i=0; i<source.size(); i++) {
		for (size_t j=0; j<ns[i]; j++) {
			b[k] = source[i].hasColors[j];
			k++;
		}
	}
	return b;
}


std::vector<double> SpatRaster::cellFromXY (std::vector<double> x, std::vector<double> y, double missing) {
// size of x and y should be the same

	size_t size = x.size();
	std::vector<double> cells(size);

	SpatExtent extent = getExtent();
	double yr_inv = nrow() / (extent.ymax - extent.ymin);
	double xr_inv = ncol() / (extent.xmax - extent.xmin);

	for (size_t i = 0; i < size; i++) {
		// cannot use trunc here because trunc(-0.1) == 0
		long row = std::floor((extent.ymax - y[i]) * yr_inv);
		// points in between rows go to the row below
		// except for the last row, when they must go up
		if (y[i] == extent.ymin) {
			row = nrow()-1 ;
		}

		long col = std::floor((x[i] - extent.xmin) * xr_inv);
		// as for rows above. Go right, except for last column
		if (x[i] == extent.xmax) {
			col = ncol() - 1 ;
		}
		long nr = nrow();
		long nc = ncol();
		if (row < 0 || row >= nr || col < 0 || col >= nc) {
			cells[i] = missing;
		} else {
			cells[i] = row * ncol() + col;
		}
	}

	return cells;
}


double SpatRaster::cellFromXY (double x, double y, double missing) {
	std::vector<double> X = {x};
	std::vector<double> Y = {y};
	std::vector<double> cell = cellFromXY(X, Y, missing);
	return  cell[0];
}


std::vector<double> SpatRaster::cellFromRowCol(std::vector<int_64> row, std::vector<int_64> col) {
	recycle(row, col);
	size_t n = row.size();
	std::vector<double> result(n);
	int_64 nr = nrow();
	int_64 nc = ncol();
	for (size_t i=0; i<n; i++) {
		result[i] = (row[i]<0 || row[i] >= nr || col[i]<0 || col[i] >= nc) ? NAN : (double)row[i] * nc + col[i];
	}
	return result;
}


double SpatRaster::cellFromRowCol (int_64 row, int_64 col) {
	std::vector<int_64> rows = {row};
	std::vector<int_64> cols = {col};
	std::vector<double> cell = cellFromRowCol(rows, cols);
	return  cell[0];
}

std::vector<double> SpatRaster::cellFromRowColCombine(std::vector<int_64> row, std::vector<int_64> col) {
	size_t n = row.size();
	size_t m = col.size();
	double nc = ncol();
	double nr = nrow();
	std::vector<double> x;
	x.reserve(n * m);
	for (size_t i=0; i<n; i++) {
		for (size_t j=0; j<m; j++) {
			if (row[i] < 0 || row[i] >= nr || col[j]<0 || col[j] >= nc) {
				x.push_back(NAN);
			} else {
				x.push_back(row[i] * nc + col[j]);
			}
		}
	}
	return x;
}


double SpatRaster::cellFromRowColCombine(int_64 row, int_64 col) {
	return cellFromRowCol(row, col);
}


std::vector<double> SpatRaster::yFromRow(const std::vector<int_64> &row) {
	size_t size = row.size();
	std::vector<double> result( size );
	SpatExtent extent = getExtent();
	double ymax = extent.ymax;
	double yr = yres();
	int_64 nr = nrow();

	for (size_t i = 0; i < size; i++) {
		result[i] = (row[i] < 0 || row[i] >= nr ) ? NAN : ymax - ((row[i]+0.5) * yr);
	}
	return result;
}

double SpatRaster::yFromRow (int_64 row) {
	std::vector<int_64> rows = {row};
	std::vector<double> y = yFromRow(rows);
	return y[0];
}



std::vector<double> SpatRaster::xFromCol(const std::vector<int_64> &col) {
	size_t size = col.size();
	std::vector<double> result( size );
	SpatExtent extent = getExtent();
	double xmin = extent.xmin;
	double xr = xres();
	int_64 nc = ncol();
	for (size_t i = 0; i < size; i++) {
		result[i] = (col[i] < 0 || col[i] >= nc ) ? NAN : xmin + ((col[i]+0.5) * xr);
	}
	return result;
}

double SpatRaster::xFromCol(int_64 col) {
	std::vector<int_64> cols = {col};
	std::vector<double> x = xFromCol(cols);
	return x[0];
}

std::vector<int_64> SpatRaster::colFromX(const std::vector<double> &x) {

	SpatExtent extent = getExtent();

	double xmin = extent.xmin;
	double xmax = extent.xmax;
	double xr = xres();
	size_t xs = x.size();
	std::vector<int_64> result(xs, -1);

	for (size_t i = 0; i < xs; i++) {
		if (x[i] >= xmin && x[i] < xmax ) {
			result[i] =  trunc((x[i] - xmin) / xr);
		} else if (x[i] == xmax) {
			result[i] = ncol()-1;
		}
	}
	return result;
}


int_64 SpatRaster::colFromX(double x) {
	std::vector<double> xv = {x};
	return colFromX(xv)[0];
}


std::vector<int_64> SpatRaster::rowFromY(const std::vector<double> &y) {

	SpatExtent extent = getExtent();
	double ymin = extent.ymin;
	double ymax = extent.ymax;
	double yr = yres();
	size_t ys = y.size();
	std::vector<int_64> result(ys, -1);

	for (size_t i = 0; i < ys; i++) {
		if (y[i] > ymin && y[i] <= ymax) {
			result[i] = trunc((ymax - y[i]) / yr);
		} else if (y[i] == ymin) {
			result[i] = nrow() - 1;
		}
	}
	return result;
}


int_64 SpatRaster::rowFromY(double y) {
	std::vector<double> Y = {y};
	return rowFromY(Y)[0];
}


void SpatRaster::xyFromCell( std::vector<std::vector<double>> &xy ) {
	
	SpatExtent extent = getExtent();
	double xmin = extent.xmin;
	double ymax = extent.ymax;
	double yr = yres();
	double xr = xres();
	size_t nr = nrow();
	size_t nc = ncol();

	xy[0].reserve(ncell()+2); 
	xy[1].reserve(ncell()+2); 
	for (size_t i = 0; i<nr; i++) {
		for (size_t j = 0; j<nc; j++) {
			xy[0].push_back( xmin + (j + 0.5) * xr );
			xy[1].push_back( ymax - (i + 0.5) * yr );
		}
	}
}


std::vector<std::vector<double>> SpatRaster::xyFromCell( std::vector<double> &cell) {
	size_t n = cell.size();
	SpatExtent extent = getExtent();

	double xmin = extent.xmin;
	double ymax = extent.ymax;
	double yr = yres();
	double xr = xres();
    double ncells = ncell();
    size_t nc = ncol();
	std::vector< std::vector<double> > out(2, std::vector<double> (n, NAN) );
	for (size_t i = 0; i<n; i++) {
		if (std::isnan(cell[i]) || (cell[i] < 0) || (cell[i] >= ncells)) continue;
        size_t row = cell[i] / nc;
        size_t col = cell[i] - (row * nc);
        out[0][i] = xmin + (col + 0.5) * xr;
        out[1][i] = ymax - (row + 0.5) * yr;
	}
	return out;
}


std::vector< std::vector<double>> SpatRaster::xyFromCell( double cell) {
	std::vector<double> vcell = {cell};
	return xyFromCell(vcell);
}


std::vector<std::vector<int_64>> SpatRaster::rowColFromCell(std::vector<double> &cell) {
	size_t cs = cell.size();
	std::vector<std::vector<int_64>> result(2, std::vector<int_64> (cs, -1) );
	double nc = ncell();
	for (size_t i = 0; i < cs; i++) {
		if ((cell[i] >= 0) && (cell[i] < nc )) {
			result[0][i] = trunc(cell[i]/ ncol());
			result[1][i] = (cell[i] - ((result[0][i]) * ncol()));
		}
	}
	return result;
}


std::vector<std::vector<int_64>>  SpatRaster::rowColFromExtent(SpatExtent e) {
	std::vector<std::vector<double>> xy = e.asPoints();
	std::vector<int_64> col = colFromX(xy[0]);
	std::vector<int_64> row = rowFromY(xy[1]);
	std::vector<std::vector<int_64>> out = { row, col };
	return out;
}


std::vector<double> SpatRaster::adjacentMat(std::vector<double> cells, std::vector<bool> mat, std::vector<unsigned> dim, bool include) {
	std::vector<double> out;
	if ((dim.size() != 2) || (dim[0] % 2 == 0) || (dim[1] %2 == 0)) {
		setError("invalid matrix dimensions (must be odd sized)");
		return out;
	}
	if ((dim[0] == 1) && (dim[1] == 1)) {
		setError("invalid matrix dimensions (too small)");
		return out;
	}

	int dy = dim[0] / 2;
	int dx = dim[1] / 2;

	unsigned n = cells.size();
	int nngb = std::accumulate(mat.begin(), mat.end(), 0);
	out.reserve(n * (nngb + include));

    std::vector<int> offcols(nngb);
    std::vector<int> offrows(nngb);

	size_t i = 0;
	size_t j = 0;
	for (int r = -dy; r<=dy; r++) {
		for (int c = -dx; c<=dx; c++) {
			if (mat[i]) {
				offrows[j] = r;
				offcols[j] = c;
				j++;
			}
			i++;
		}
	}

	bool globlatlon = is_global_lonlat();

	std::vector<std::vector<int_64>> rc = rowColFromCell(cells);
	std::vector<int_64> r = rc[0];
	std::vector<int_64> c = rc[1];
    std::vector<int_64> cols(nngb);
    std::vector<int_64> rows(nngb);
    int_64 nc = ncol();
    int_64 lc = nc-1;

	for (size_t i=0; i<n; i++) {
		for (int j = 0; j<nngb; j++) {
			rows[j] = r[i] + offrows[j];
			cols[j] = c[i] + offcols[j];
		}

		if (globlatlon) {
			for (int j = 0; j<nngb; j++) {
				if (cols[j] < 0) cols[j] = nc + cols[j];
				if (cols[j] > lc) cols[j] = cols[j] - nc;
			}
		}
		std::vector<double> adjcells = cellFromRowCol(rows, cols);
        if (include) {
			out.push_back(cells[i]);
        }
		out.insert(out.end(), adjcells.begin(), adjcells.end());
	}
	return out;
}

std::vector<double> SpatRaster::adjacent(std::vector<double> cells, std::string directions, bool include) {

	std::vector<double> out;

	std::vector<std::string> f {"rook", "queen", "bishop", "4", "8", "16"};
	if (std::find(f.begin(), f.end(), directions) == f.end()) {
        setError("argument directions is not valid");
        return(out);
	}
	unsigned n = cells.size();

	unsigned nngb = (directions=="queen" || directions=="8") ? 8 : (directions=="16" ? 16 : 4);
	nngb += include;
	out.reserve(n * nngb);

	std::vector<std::vector<int_64>> rc = rowColFromCell(cells);
	std::vector<int_64> r = rc[0];
	std::vector<int_64> c = rc[1];
	bool globlatlon = is_global_lonlat();
    int_64 nc = ncol();
    int_64 lc = nc-1;
    std::vector<int_64> cols, rows;
	if (directions == "rook" || directions == "4") {
		for (size_t i=0; i<n; i++) {
			rows = {r[i]-1, r[i]   , r[i]  , r[i]+1};
            cols = {c[i]  , c[i]-1 , c[i]+1, c[i]};
            if (globlatlon) {
                if (c[i]==0) {
                    cols[1] = lc;
                } else if (c[i]==lc) {
                    cols[2] = 0;
                }
            }
            if (include) {
				out.push_back(cells[i]);
            }
			std::vector<double> adjcells = cellFromRowCol(rows, cols);
			out.insert(out.end(), adjcells.begin(), adjcells.end());
		}
	} else if (directions == "queen" || directions == "8") {
		for (size_t i=0; i<n; i++) {
            rows = {r[i]-1, r[i]-1, r[i]-1, r[i], r[i], r[i]+1, r[i]+1, r[i]+1};
            cols = {c[i]-1, c[i], c[i]+1, c[i]-1, c[i]+1, c[i]-1, c[i], c[i]+1};
            if (globlatlon) {
                if (c[i]==0) {
                    cols = {lc, c[i], c[i]+1, lc, c[i]+1, lc, c[i], c[i]+1};
                } else if (c[i]==lc) {
                    cols = {c[i]-1, c[i], 0, c[i]-1, 0, c[i]-1, c[i], 0};
                }
            }
            if (include) {
				out.push_back(cells[i]);
            }
			std::vector<double> adjcells = cellFromRowCol(rows, cols);
			out.insert(out.end(), adjcells.begin(), adjcells.end());
		}
	} else if (directions == "bishop") {
		for (size_t i=0; i<n; i++) {
            rows = {r[i]-1, r[i]-1, r[i]+1, r[i]+1};
            cols = {c[i]-1, c[i]+1, c[i]-1, c[i]+1};
            if (globlatlon) {
                if (c[i]==0) {
                    cols = {lc, c[i]+1, lc, c[i]+1};
                } else if (c[i]==lc) {
                    cols = {c[i]-1, 0, c[i]-1, 0};
                }
            }
            if (include) {
				out.push_back(cells[i]);
            }
			std::vector<double> adjcells = cellFromRowCol(rows, cols);
			out.insert(out.end(), adjcells.begin(), adjcells.end());
		}
	} else if (directions == "16") {
		for (size_t i=0; i<n; i++) {
            rows = {r[i]-2, r[i]-2, r[i]-1, r[i]-1, r[i]-1, r[i]-1, r[i]-1, r[i]  , r[i]  , r[i]+1, r[i]+1, r[i]+1, r[i]+1, r[i]+1, r[i]+2, r[i]+2};
            cols = {c[i]-1, c[i]+1, c[i]-2, c[i]-1, c[i],   c[i]+1, c[i]+2, c[i]-1, c[i]+1, c[i]-2, c[i]-1, c[i]  , c[i]+1, c[i]+2, c[i]-1, c[i]+1};
            if (globlatlon) {
                if ((c[i]==0) || (c[i]==1)) {
                    for (size_t j=0; j<16; j++) {
                        cols[j] = (cols[j] < 0) ? nc-cols[j] : cols[j];
                    }
                } else if (c[i]==nc) {
                    for (size_t j=0; j<16; j++) {
                        cols[j] = (cols[j] > lc) ? cols[j]-nc : cols[j];
                    }
                }
            }
            if (include) {
				out.push_back(cells[i]);
            }
			std::vector<double> adjcells = cellFromRowCol(rows, cols);
			out.insert(out.end(), adjcells.begin(), adjcells.end());
		}
	}
	return(out);
}


SpatVector SpatRaster::as_multipoints(bool narm, bool nall, SpatOptions &opt) {

	BlockSize bs = getBlockSize(opt);
    size_t ncl = ncell();
	SpatVector pv;
	pv.reserve(1);

    std::vector<std::vector<double>> xy;
	if (!narm) {
        for (size_t i=0; i<ncl; i++) {
            xy = xyFromCell(i);
			SpatPart p(xy[0], xy[1]);
			SpatGeom g(p, points);
			pv.addGeom(g);
			g.parts.resize(0);
        }
		return pv;
	}

	if (!readStart()) {
		pv.setError(getError());
		return(pv);
	}

	size_t nc = ncol();
	unsigned nl = nlyr();
	std::vector<double> v, x, y;
	for (size_t i = 0; i < bs.n; i++) {
		readValues(v, bs.row[i], bs.nrows[i], 0, nc);
        size_t off1 = (bs.row[i] * nc);
 		size_t vnc = bs.nrows[i] * nc;
		for (size_t j=0; j<vnc; j++) {
			if (nall) {
				bool allna = true;
				for (size_t lyr=0; lyr<nl; lyr++) {
					size_t off2 = lyr*vnc;
					if (!std::isnan(v[off2+j])) {
						allna = false;
						continue;
					}
				}
				if (allna) continue;
			} else {
				bool foundna = false;
				for (size_t lyr=0; lyr<nl; lyr++) {
					size_t off2 = lyr*vnc;
					if (std::isnan(v[off2+j])) {
						foundna = true;
						continue;
					}
				}
				if (foundna) continue;
			}
			xy = xyFromCell( off1+j );
			x.push_back(xy[0][0]);
			y.push_back(xy[1][0]);
		}
	}
	SpatPart p(x, y);
	SpatGeom g(p, points);
	pv.addGeom(g);

	readStop();
	pv.srs = source[0].srs;
	return(pv);
}


SpatVector SpatRaster::as_points(bool values, bool narm, bool nall, SpatOptions &opt) {

	BlockSize bs = getBlockSize(opt);
    size_t ncl = ncell();
	SpatVector pv;
	pv.reserve(ncl);
	pv.srs = source[0].srs;

	if (!hasValues()) {
		if (values) {
			pv.addWarning("raster has no values");
		}
		values = false;
		narm = false;
	}

    std::vector<std::vector<double>> xy;
	if ((!values) && (!narm)) {
        for (size_t i=0; i<ncl; i++) {
            xy = xyFromCell(i);
			SpatPart p(xy[0], xy[1]);
			SpatGeom g(p, points);
			pv.addGeom(g);
			//g.parts.resize(0);
        }
		return pv;
	}

	if (values) {
        std::vector<std::string> nms = getNames();
        for (size_t i=0; i<nlyr(); i++) {
            pv.df.add_column(0, nms[i]);
        }
	}
	if (!readStart()) {
		pv.setError(getError());
		return(pv);
	}

	size_t nc = ncol();
	unsigned nl = nlyr();
	std::vector<double> v;
	for (size_t i = 0; i < bs.n; i++) {
		readValues(v, bs.row[i], bs.nrows[i], 0, nc);
        size_t off1 = (bs.row[i] * nc);
 		size_t vnc = bs.nrows[i] * nc;
		if (narm) {
			if (values) {
				pv.df.reserve(ncl);
			}
			for (size_t j=0; j<vnc; j++) {
				if (nall) {
					bool allna = true;
					for (size_t lyr=0; lyr<nl; lyr++) {
						size_t off2 = lyr*vnc;
						if (!std::isnan(v[off2+j])) {
							allna = false;
							continue;
						}
					}
					if (allna) continue;
				} else {
					bool foundna = false;
					for (size_t lyr=0; lyr<nl; lyr++) {
						size_t off2 = lyr*vnc;
						if (std::isnan(v[off2+j])) {
							foundna = true;
							continue;
						}
					}
					if (foundna) continue;
				}


                xy = xyFromCell( off1+j );
                SpatPart p(xy[0], xy[1]);
                SpatGeom g(p, points);
                pv.addGeom(g);
                if (values) {
                    for (size_t lyr=0; lyr<nl; lyr++) {
                        unsigned off2 = lyr*vnc;
                        pv.df.dv[lyr].push_back(v[off2+j]);
                    }
                }
			}
		} else { // if (values) {
			for (size_t j=0; j<vnc; j++) {
                xy = xyFromCell(off1+j);
                SpatPart p(xy[0], xy[1]);
                SpatGeom g(p, points);
                pv.addGeom(g);
			}
            for (size_t lyr=0; lyr<nl; lyr++) {
				size_t off2 = lyr*vnc;
				pv.df.dv[lyr] = std::vector<double>(v.begin()+off2, v.begin()+off2+vnc);
			}
		}
	}
	readStop();
//	pv.srs = source[0].srs;
	return(pv);
}

std::vector<std::vector<double>> SpatRaster::as_points_value(const double& target, SpatOptions &opt) {

	std::vector<std::vector<double>> xy(2);
	if (nlyr() > 1) {
		setError("can only process one layer");
		return xy;
	}

	BlockSize bs = getBlockSize(opt);

	if (!readStart()) {
		return(xy);
	}

	size_t nc = ncol();
    size_t ncl = ncell();
	std::vector<double> cells;
	cells.reserve(std::min(ncl/10, (size_t)10000));

	std::vector<double> v;
	if (std::isnan(target)) {
		for (size_t i = 0; i < bs.n; i++) {
			readValues(v, bs.row[i], bs.nrows[i], 0, nc);
			size_t base = (bs.row[i] * nc);
			size_t szv = v.size();
			for (size_t j=0; j<szv; j++) {
				if (std::isnan(v[j])) {
					cells.push_back(base+j);
				}
			}
		}
	} else {
		for (size_t i = 0; i < bs.n; i++) {
			readValues(v, bs.row[i], bs.nrows[i], 0, nc);
			size_t base = (bs.row[i] * nc);
			size_t szv = v.size();
			for (size_t j=0; j<szv; j++) {
				if (v[j] == target) {
					cells.push_back(base+j);
				}
			}
		}
	}
	readStop();
	return xyFromCell(cells);
}



std::vector<std::vector<double>> SpatRaster::coordinates(bool narm, bool nall, SpatOptions &opt) {

    std::vector<std::vector<double>> xy(2);

	if ( !(narm) || (!hasValues()) ) {
        xyFromCell(xy);
		return xy;
	}

	BlockSize bs = getBlockSize(opt);

	if (!readStart()) {
		return(xy);
	}
	size_t nc = ncol();
	unsigned nl = nlyr();
	std::vector<double> v;
	for (size_t i = 0; i < bs.n; i++) {
		readValues(v, bs.row[i], bs.nrows[i], 0, nc);
        size_t off1 = (bs.row[i] * nc);
 		size_t vnc = bs.nrows[i] * nc;
		for (size_t j=0; j<vnc; j++) {
			if (nall) {
				bool allna = true;
				size_t off2 = 0;
				for (size_t lyr=0; lyr<nl; lyr++) {
					if (!std::isnan(v[off2+j])) {
						allna = false;
						continue;
					}
					off2 += vnc;
				}
				if (allna) continue;
			} else {
				bool foundna = false;
				size_t off2 = 0;
				for (size_t lyr=0; lyr<nl; lyr++) {
					if (std::isnan(v[off2+j])) {
						foundna = true;
						continue;
					}
					off2 += vnc;
				}
				if (foundna) continue;
			}
			std::vector<std::vector<double>> xyc = xyFromCell( off1+j );
			xy[0].push_back(xyc[0][0]);
			xy[1].push_back(xyc[1][0]);
		}
	}
	readStop();
	return(xy);
}


std::vector<std::vector<double>> SpatRaster::cells_notna(SpatOptions &opt) {

	std::vector<std::vector<double>> out(2);
	if (nlyr() > 1) {
		setError("can only process one layer");
		return out;
	}

	BlockSize bs = getBlockSize(opt);

	if (!readStart()) {
		return(out);
	}

	size_t nc = ncol();
    size_t ncl = ncell();
	size_t rs = std::max(ncl/50, (size_t)10000);
	out[0].reserve(rs);
	out[1].reserve(rs);

	for (size_t i = 0; i < bs.n; i++) {
		std::vector<double> v;
		readValues(v, bs.row[i], bs.nrows[i], 0, nc);
		size_t base = (bs.row[i] * nc);
		size_t szv = v.size();
		for (size_t j=0; j<szv; j++) {
			if (!std::isnan(v[j])) {
				out[0].push_back(base+j); // cell
				out[1].push_back(v[j]); // value
			}
		}
	}
	readStop();
	return out;
}

std::vector<double> SpatRaster::cells_notna_novalues(SpatOptions &opt) {


	if (nlyr() > 1) {
		SpatOptions topt(opt);
		SpatRaster x = nonan(true, topt);
		return x.cells_notna_novalues(opt);
	}
	
	std::vector<double> out;
	BlockSize bs = getBlockSize(opt);

	if (!readStart()) {
		return(out);
	}

	size_t nc = ncol();
    size_t ncl = ncell();
	size_t rs = std::max(ncl/500, (size_t)10000);
	out.reserve(rs);

	for (size_t i = 0; i < bs.n; i++) {
		std::vector<double> v;
		readValues(v, bs.row[i], bs.nrows[i], 0, nc);
		size_t base = (bs.row[i] * nc);
		size_t szv = v.size();
		for (size_t j=0; j<szv; j++) {
			if (!std::isnan(v[j])) {
				out.push_back(base+j); // cell
			}
		}
	}
	readStop();
	return out;
}


void getCorners(std::vector<double> &x,  std::vector<double> &y, const double &X, const double &Y, const double &xr, const double &yr) {
	x[0] = X - xr;
	y[0] = Y - yr;
	x[1] = X - xr;
	y[1] = Y + yr;
	x[2] = X + xr;
	y[2] = Y + yr;
	x[3] = X + xr;
	y[3] = Y - yr;
	x[4] = x[0];
	y[4] = y[0];
}


SpatVector SpatRaster::as_polygons(bool round, bool dissolve, bool values, bool narm, bool nall, int digits, SpatOptions &opt) {

	if (!hasValues()) {
		values = false;
		narm = false;
		dissolve=false;
	}

	if (dissolve) {
		return polygonize(round, values, narm, dissolve, digits, opt);
	}

	SpatVector vect;
	opt.ncopies = 12;
	if (!canProcessInMemory(opt)) {
		if (ncell() > 1000000) { // for testing with canPIM=false
			vect.setError("the raster is too large");
			return vect;
		}
	}

	bool remove_values = false;
	if (narm) {
		if (!values) remove_values = true;
		values=true;
	}

	unsigned nl = nlyr();
	unsigned nc = ncell();
	if (values) {
		std::vector<double> v = getValues(-1, opt);
		std::vector<std::string> nms = getNames();
		make_unique_names(nms);
		for (size_t i=0; i<nl; i++) {
			size_t offset = i * nc;
			std::vector<double> vv(v.begin()+offset, v.begin()+offset+nc);
			vect.add_column(vv, nms[i]);
		}
	}

	SpatGeom g;
	g.gtype = polygons;
	double xr = xres()/2;
	double yr = yres()/2;
	std::vector<double> x(5);
	std::vector<double> y(5);

	std::vector<double> cells(ncell()) ;
	std::iota (std::begin(cells), std::end(cells), 0);
	std::vector< std::vector<double> > xy = xyFromCell(cells);
	vect.reserve(cells.size());
	for (int i=nc-1; i>=0; i--) {
		if (narm) {
			bool erase;
			if (nall) {
				erase = true;
				for (size_t j=0; j<nl; j++) {
					if (!std::isnan(vect.df.dv[j][i])) {
						erase=false;
						break;
					}
				}
			} else {
				erase = false;
				for (size_t j=0; j<nl; j++) {
					if (std::isnan(vect.df.dv[j][i])) {
						erase=true;
						break;
					}
				}
			}
			if (erase) {
				for (size_t j=0; j<nl; j++) {
					vect.df.dv[j].erase (vect.df.dv[j].begin()+i);
				}
				continue; // skip the geom
			}
		}
		getCorners(x, y, xy[0][i], xy[1][i], xr, yr);
		SpatPart p(x, y);
		g.addPart(p);
		vect.addGeom(g);
		g.parts.resize(0);
	}

	std::reverse(std::begin(vect.geoms), std::end(vect.geoms));

//	if (dissolve) {
//		vect = vect.aggregate(vect.get_names()[0], true);
//	}

	if (remove_values) {
		vect.df = SpatDataFrame();
	}
	vect.srs = source[0].srs;
	return(vect);
}


SpatVector SpatRaster::as_lines(SpatOptions &opt) {

	SpatVector vect;
	opt.ncopies = 12;
	if (!canProcessInMemory(opt)) {
		if (ncell() > 1000000) { // for testing with canPIM=false
			vect.setError("the raster is too large");
			return vect;
		}
	}

	SpatGeom g;
	g.gtype = lines;

	std::vector<int_64> cols(ncol());
	std::vector<int_64> rows(nrow());
	std::iota(std::begin(rows), std::end(rows), 0);
	std::iota(std::begin(cols), std::end(cols), 0);
	std::vector<double> x = xFromCol(cols);
	std::vector<double> y = yFromRow(rows);

	double xr = xres()/2;
	double yr = yres()/2;
	for (double &d : x) d = d - xr;
	for (double &d : y) d = d + yr;
	x.push_back(x[x.size()-1] + xres());
	y.push_back(y[y.size()-1] - yres());

	SpatExtent e = getExtent();
	for (size_t i=0; i<x.size(); i++) {

		std::vector<double> xc = {x[i], x[i]};
		std::vector<double> yc = {e.ymin, e.ymax};
		SpatPart p(xc, yc);
		g.addPart(p);
		vect.addGeom(g);
		g.parts.resize(0);
	}
	for (size_t i=0; i<y.size(); i++) {
		std::vector<double> xc = {e.xmin, e.xmax};
		std::vector<double> yc = {y[i], y[i]};
		SpatPart p(xc, yc);
		g.addPart(p);
		vect.addGeom(g);
		g.parts.resize(0);
	}

	vect.srs = source[0].srs;
	return(vect);
}



bool SpatRaster::setRGB(int r, int g, int b, int alpha, std::string type) {
	std::vector<int> channels;
	if (alpha >= 0) {
		channels = {r, g, b, alpha};
	} else {
		channels = {r, g, b};
	}
	size_t mxlyr = vmax( channels, false );
	if (nlyr() <= mxlyr) {
		//addWarning("layer number for R, G, B, cannot exceed the number of layers");
		return false;
	} else {
		size_t mnlyr =  vmin( channels, false );;
		if (mnlyr >= 0) {
			rgblyrs = channels;
			std::vector<std::string> f = {"rgb", "hsv", "hsi", "hsl"};
			std::transform(type.begin(), type.end(), type.begin(), ::tolower);
			if (std::find(f.begin(), f.end(), type) == f.end()) {
				addWarning("color type must be one of: 'rgb', 'hsv', 'hsi', 'hsl'");
				type = "rgb";
			}
			rgbtype = type; 
			rgb = true;
		} else {
			rgb = false;
			return false;
		}
	}
	return true;
}

std::vector<int> SpatRaster::getRGB(){
	return rgblyrs;
}

void SpatRaster::removeRGB(){
	rgblyrs = std::vector<int>(0);
	rgbtype = "";
	rgb = false;
}


bool SpatRaster::to_memory(SpatOptions &opt) {
	if ((nsrc() == 1) && (source[0].memory)) {
		return true;
	}
	SpatRaster g = geometry();
	SpatRasterSource s = g.source[0];
	s.hasValues = true;
	s.memory = true;
	s.names = getNames();
	s.driver = "memory";
	source[0].values = getValues(-1, opt);
	return true;
}


SpatRaster SpatRaster::to_memory_copy(SpatOptions &opt) {
	SpatRaster m = geometry();
	std::vector<double> v = getValues(-1, opt);
	m.setValues(v, opt);
	return m;
}


std::vector<int> SpatRaster::getFileBlocksize() {
	std::vector<int> b;
	b.reserve(2 * nlyr());
	for (size_t i=0; i<source.size(); i++) {
		b.insert(b.end(), source[i].blockrows.begin(), source[i].blockrows.end());
	}
	for (size_t i=0; i<source.size(); i++) {
		b.insert(b.end(), source[i].blockcols.begin(), source[i].blockcols.end());
	}
	return b;

}


bool SpatRaster::addTag(std::string name, std::string value) {
	lrtrim(name);
	lrtrim(value);
	if (value == "") {
		return removeTag(name);
	} else if (name != "") {
		tags[name] = value;
		return true;
	} 
	return false;
}

bool SpatRaster::removeTag(std::string name) {
	std::map<std::string, std::string>::iterator it = tags.find(name);
	if (it == tags.end()) return false;
	tags.erase(it);
	return true;
}

std::string SpatRaster::getTag(std::string name) {
	std::map<std::string, std::string>::iterator it = tags.find(name);
	if (it != tags.end()) return it->second;
	return "";
}

std::vector<std::string> SpatRaster::getTags() {
	std::vector<std::string> out;
	out.reserve(2 * tags.size());
	for(auto e : tags) {
		out.push_back(e.first);
		out.push_back(e.second);
	}
	return out;
}



void SpatRaster::addLyrTags(std::vector<size_t> lyrs, std::vector<std::string> names, std::vector<std::string> values) {

	size_t n = std::max(std::max(lyrs.size(), names.size()), values.size());
	if (n == 0) return;
	
	recycle(lyrs, n);
	recycle(names, n);
	recycle(values, n);
	
	size_t nl = nlyr();
	for (size_t i=0; i<lyrs.size(); i++) {
		if (lyrs[i] >= nl) continue;
		lrtrim(names[i]);
		lrtrim(values[i]);
		if (values[i] == "") {
			removeLyrTag(lyrs[i], names[i]);
		} else {
			if (lyrs[i] >= lyrTags.size()) lyrTags.resize(lyrs[i]+1);
			if (names[i] != "") {
				lyrTags[lyrs[i]][names[i]] = values[i];
			} 
		}
	}
}

bool SpatRaster::removeLyrTag(size_t lyr, std::string name) {
	if (lyr >= lyrTags.size()) return false;
	std::map<std::string, std::string>::iterator it = lyrTags[lyr].find(name);
	if (it == lyrTags[lyr].end()) return false;
	lyrTags[lyr].erase(it);
	return true;
}

bool SpatRaster::removeLyrTags() {
	lyrTags.resize(0);
	return true;
}


std::string SpatRaster::getLyrTag(size_t lyr, std::string name) {
	if (lyr >= lyrTags.size()) return "";
	std::map<std::string, std::string>::iterator it = lyrTags[lyr].find(name);
	if (it != lyrTags[lyr].end()) return it->second;
	return "";
}

std::vector<std::string> SpatRaster::getLyrTags(std::vector<size_t> lyrs) {
	std::vector<std::string> out;
	out.reserve(lyrs.size());
	for (size_t i=0; i<lyrs.size(); i++) {
		if (lyrs[i] < lyrTags.size()) {
			for(auto e : lyrTags[lyrs[i]]) {
				out.push_back(std::to_string(lyrs[i]));
				out.push_back(e.first);
				out.push_back(e.second);
			}
		}
	}
	return out;
}

