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
#include "time.h"
#include "recycle.h"

#include <set>

#ifdef useGDAL
#include "crs.h"
#endif


SpatRaster::SpatRaster(std::string fname, std::vector<int> subds, std::vector<std::string> subdsname) {
#ifdef useGDAL
	constructFromFile(fname, subds, subdsname);
#endif
}


SpatRaster::SpatRaster(std::vector<std::string> fname, std::vector<int> subds, std::vector<std::string> subdsname, std::string x) {
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



SpatRaster SpatRaster::geometry(long nlyrs, bool properties) {
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


std::vector<std::string> SpatRaster::getLongNames() {
	std::vector<std::string> x;
	for (size_t i=0; i<source.size(); i++) {
		x.insert(x.end(), source[i].long_names.begin(), source[i].long_names.end());
	}
	return(x);
}


bool SpatRaster::setLongNames(std::vector<std::string> nms) {
	if (nms.size() == 1) {
        size_t begin=0;
        for (size_t i=0; i<source.size(); i++)	{
            size_t end = begin + source[i].nlyr;
			size_t sz =  end - begin + 1;
            source[i].long_names = std::vector<std::string> (sz, nms[0]);
            begin = end;
        }
        return true;
	} else if (nms.size() != nlyr()) {
		return false;
	} else {
        size_t begin=0;
        size_t end;
        for (size_t i=0; i<source.size(); i++)	{
            end = begin + source[i].nlyr;
            source[i].long_names = std::vector<std::string> (nms.begin() + begin, nms.begin() + end);
            begin = end;
        }
        return true;
	}
}



bool SpatRaster::hasTime() {
	bool test = true;
	for (size_t i=0; i<source.size(); i++) {
		test = test & source[i].hasTime; 
	}
	return(test);
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

std::vector<std::string> SpatRaster::getTimeStr() {
	std::vector<std::string> out;
	if (source[0].timestep == "seconds") {
		std::vector<int_64> time = getTime();
		out.reserve(time.size());
		for (size_t i=0; i < out.size(); i++) {
			std::vector<int> x = get_date(time[i]);
			if (x.size() > 2) {
				out.push_back( std::to_string(x[0]) + "-" 
						  + std::to_string(x[1]) + "-"
						  + std::to_string(x[2]) );
						  
			} else {
				out.push_back("");
			}
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

bool SpatRaster::setTime(std::vector<int_64> time, std::string step) {
	if (time.size() != nlyr()) {
		return false;
	} 
	if (!(step == "seconds") || (step == "raw")) {  // "days", "months", "years"
		return false;
	} 
	size_t begin=0;
	for (size_t i=0; i<source.size(); i++)	{
		size_t end = begin + source[i].nlyr;
        source[i].time = std::vector<int_64> (time.begin() + begin, time.begin() + end);
		source[i].timestep = step;
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
        for (size_t i=0; i<source.size(); i++)	{
            source[i].unit = std::vector<std::string> (source[i].nlyr, units[0]);
        }
        return true;
	} else if (units.size() != nlyr()) {
		return false;
	} else {
        size_t begin=0;
        for (size_t i=0; i<source.size(); i++)	{
            size_t end = begin + source[i].nlyr;
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


bool SpatRaster::valid_sources(bool files, bool rotated) {
	std::vector<std::string> ff;
	for (size_t i=0; i<source.size(); i++) { 
		std::string f = source[i].filename; 
		if (f == "") continue;
		if (files) {
			std::size_t found = f.find(":"); // perhaps http: or PG:xxx
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

bool SpatRaster::hasWindow() {
	return source[0].hasWindow;
}

bool SpatRaster::removeWindow() {
	if (hasWindow()) {
		SpatExtent e = source[0].window.full_extent;
		setExtent(e, true, "");
		for (size_t i=0; i<source.size(); i++) {
			source[i].hasWindow = false;
			source[i].nrow = source[0].window.full_nrow;
			source[i].ncol = source[0].window.full_ncol;
		}
	} 
	return true;
}





bool SpatRaster::setWindow(SpatExtent x) {

	if ( !x.valid() ) {
		setError("invalid extent");
		return false;
	} 

	if (hasWindow()) {
		removeWindow();
	}

	x = align(x, "near");
	SpatExtent e = getExtent();
	if (x.compare(e, "==", 0.1 * xres())) {
		return true;
	}
	
	e.intersect(x);
	if ( !e.valid() ) {
		setError("extents do not overlap");
		return false;
	} 

// get read-window
	double xr = xres();
	double yr = yres();

	bool expand = false;
	std::vector<uint_64> rc(2);
	std::vector<uint_64> exp(4, 0);

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

	for (size_t i=0; i<source.size(); i++) {
		source[i].window.off_row = rc[0];
		source[i].window.off_col = rc[1];
		source[i].window.expand = exp;
		source[i].window.expanded  = expand;
		source[i].window.full_extent = getExtent();
		source[i].window.full_nrow   = source[i].nrow;
		source[i].window.full_ncol   = source[i].ncol;
		source[i].hasWindow     = true;
	}
	setExtent(x, true, "");		

	return true;
}


