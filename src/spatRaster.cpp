// Copyright (c) 2018-2020  Robert J. Hijmans
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
#include <set>

#ifdef useGDAL
#include "crs.h"
#endif


SpatRaster::SpatRaster(std::string fname, int subds, std::string subdsname) {
#ifdef useGDAL
	constructFromFile(fname, subds, subdsname);
#endif
}


SpatRaster::SpatRaster(std::vector<std::string> fname, int subds, std::string subdsname, std::string x) {
// argument "x" is ignored. It is only there to have four arguments such that the Rcpp module
// can distinguish this constructor from another with three arguments. 	
#ifdef useGDAL
	constructFromFile(fname[0], subds, subdsname);
	for (size_t i=1; i<fname.size(); i++) {
		SpatRaster r;
		bool ok = r.constructFromFile(fname[i], subds, subdsname);
		if (ok) {
			addSource(r);
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


void SpatRaster::setSources(std::vector<RasterSource> s) {
	source = s;
//	extent = s[0].extent;
//	srs = s[0].srs;
}


void SpatRaster::setSource(RasterSource s) {
	s.resize(s.nlyr);
	std::vector<RasterSource> vs = {s};
	setSources(vs);
}


SpatRaster::SpatRaster(RasterSource s) {
	std::vector<RasterSource> vs = {s};
	setSources(vs);
}


SpatRaster::SpatRaster() {

	RasterSource s;
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
	s.layers.resize(1, 0);
	s.datatype = "";
	s.names = {"lyr.1"};
	s.srs.proj4 = "+proj=longlat +datum=WGS84";
	s.srs.wkt = "GEOGCS[\"WGS 84\", DATUM[\"WGS_1984\", SPHEROID[\"WGS 84\",6378137,298.257223563]], PRIMEM[\"Greenwich\",0], UNIT[\"degree\",0.0174532925199433]]";
	setSource(s);
}



SpatRaster::SpatRaster(std::vector<unsigned> rcl, std::vector<double> ext, std::string crs) {

	RasterSource s;
	s.nrow=rcl[0];
	s.ncol=rcl[1];
	s.extent.xmin = ext[0];
	s.extent.xmax = ext[1];
	s.extent.ymin = ext[2];
	s.extent.ymax = ext[3];
	s.hasValues = false;
	s.hasRange = {false};

	s.memory = true;
	s.filename = "";
	//s.driver = "";
	s.nlyr = rcl[2];
	s.layers.resize(1, 0);
	//s.layers.resize(1, s.nlyr);
	//std::iota(s.layers.begin(), s.layers.end(), 0);

	s.datatype = "";

#ifdef useGDAL
	std::string msg;
	if (!s.srs.set( crs, msg )) {
		setError(msg);
		return;
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

	RasterSource s;
	s.nrow = nr;
	s.ncol = nc;
	s.extent = ext;
	s.hasValues = false;
	s.memory = true;
	s.filename = "";
	//s.driver = "";
	s.nlyr = nl;
	s.hasRange = { false };
	s.layers.resize(1, 0);
	//s.layers.resize(1, _nlyr);
	//std::iota(s.layers.begin(), s.layers.end(), 0);
	s.datatype = "";
#ifdef useGDAL
	std::string msg;
	if (!s.srs.set(crs, msg )) {
		setError(msg);
		return;
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



SpatRaster SpatRaster::geometry(long nlyrs) {
	RasterSource s;
	s.values.resize(0);
	s.nrow = nrow();
	s.ncol = ncol();
	s.extent = getExtent();
	s.srs = source[0].srs;
	//s.prj = prj;
	s.memory = true;
	s.hasValues = false;
	long nl = nlyr();
	bool keepnlyr = ((nlyrs == nl) | (nlyrs < 1));
	nlyrs = (keepnlyr) ? nlyr(): nlyrs;
	s.resize(nlyrs);
	std::vector<std::string> nms;
	if (keepnlyr) {
		nms = getNames();
	} else {
		for (size_t i=0; i < s.nlyr; i++) {
			nms.push_back("lyr" + std::to_string(i+1));
		}
	}
	s.names = nms;
	SpatRaster out(s);
	return out;
}


SpatRaster SpatRaster::deepCopy() {
	SpatRaster out = *this;
	return out;
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
	unsigned nc = round((e.xmax-e.xmin) / xres);
	unsigned nr = round((e.ymax-e.ymin) / yres);
	double xmax = e.xmin + nc * xres;
	double ymax = e.ymin + nr * yres;
	unsigned nl = nlyr();
	std::vector<unsigned> rcl = {nr, nc, nl};
	std::vector<double> ext = {e.xmin, xmax, e.ymin, ymax};

	out = SpatRaster(rcl, ext, {""});
	out.source[0].srs = source[0].srs;
	return out;
}


uint_64 SpatRaster::ncol() {
	if (source.size() > 0) {
		return source[0].ncol;
	} else {
		return 0;
	}
}

uint_64 SpatRaster::nrow() {
	if (source.size() > 0) {
		return source[0].nrow;
	} else {
		return 0;
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
	return source[0].srs.is_lonlat();
}

bool SpatRaster::could_be_lonlat() {
	SpatExtent e = getExtent();
	return source[0].srs.could_be_lonlat(e);
}


bool SpatRaster::is_global_lonlat() {
	SpatExtent e = getExtent();
	return source[0].srs.is_global_lonlat(e);
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

	std::string tmpbasename = tempFile(opt.get_tempdir(), "_temp_");


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
			out.addSource(rs);
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
	}
	for (size_t i = 0; i < nsrc(); i++) { 
		source[i].srs = srs; 
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


bool SpatRaster::setNames(std::vector<std::string> names) {
	if (names.size() != nlyr()) {
		return false;
	} else {
        make_valid_names(names);
        make_unique_names(names);
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


std::vector<bool> SpatRaster::hasTime() {
	std::vector<bool> x(source.size());
	for (size_t i=0; i<source.size(); i++) {
		x[i] = source[i].hasTime; 
	}
	return(x);
}


std::vector<double> SpatRaster::getTime() {
	std::vector<double> x;
	for (size_t i=0; i<source.size(); i++) {
		if (source[i].time.size() != source[i].nlyr) {
			std::vector<double> nas(source[i].nlyr, NAN);
			x.insert(x.end(), nas.begin(), nas.end());			
		} else {
			x.insert(x.end(), source[i].time.begin(), source[i].time.end());
		}
	}
	return(x);
}


bool SpatRaster::setTime(std::vector<double> times) {
	if (times.size() != nlyr()) {
		return false;
	} else {
        size_t begin=0;
        size_t end;
        for (size_t i=0; i<source.size(); i++)	{
            end = begin + source[i].nlyr;
            source[i].time = std::vector<double> (times.begin() + begin, times.begin() + end);
            begin = end;
        }
        return true;
	}
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
	if (depths.size() != nlyr()) {
		return false;
	} else {
        size_t begin=0;
        size_t end;
        for (size_t i=0; i<source.size(); i++)	{
            end = begin + source[i].nlyr;
            source[i].depth = std::vector<double> (depths.begin() + begin, depths.begin() + end);
            begin = end;
        }
        return true;
	}
}



bool SpatRaster::setUnit(std::vector<std::string> units) {
	if (units.size() == 1) {
        size_t begin=0;
        for (size_t i=0; i<source.size(); i++)	{
            size_t end = begin + source[i].nlyr;
			size_t sz =  end - begin + 1;
            source[i].unit = std::vector<std::string> (sz, units[0]);
            begin = end;
        }
        return true;
	} else if (units.size() != nlyr()) {
		return false;
	} else {
        size_t begin=0;
        size_t end;
        for (size_t i=0; i<source.size(); i++)	{
            end = begin + source[i].nlyr;
            source[i].unit = std::vector<std::string> (units.begin() + begin, units.begin() + end);
            begin = end;
        }
        return true;
	}
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

