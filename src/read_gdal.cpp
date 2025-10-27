// Copyright (c) 2018-2025  Robert J. Hijmans
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


#include <stdexcept>
#include <algorithm>
#include <stdint.h>
#include <vector>
//#include <regex>

//#include "spatRaster.h"
#include "spatRasterMultiple.h"

#include "vecmath.h"
#include "file_utils.h"
#include "string_utils.h"
#include "spatTime.h"
#include "recycle.h"
#include "gdalio.h"

//#include "NA.h"

#include "gdal_priv.h"
#include "cpl_conv.h" // for CPLMalloc()
#include "cpl_string.h"
#include "ogr_spatialref.h"

#include "gdal_rat.h"
//#include "hdr.h"

#if GDAL_VERSION_MAJOR >= 3
#include "proj.h"
#endif


void SpatRaster::gdalogrproj_init(std::string path) {
    GDALAllRegister();
    OGRRegisterAll();
	//GDALregistred = true;
#if GDAL_VERSION_MAJOR >= 3
 #ifdef PROJ_6
	if (!path.empty()) {
		const char *cp = path.c_str();
		proj_context_set_search_paths(PJ_DEFAULT_CTX, 1, &cp);
	}
 #endif
 #ifdef PROJ_71
	#ifndef __EMSCRIPTEN__
		proj_context_set_enable_network(PJ_DEFAULT_CTX, 1);
	#endif
 #endif
#endif

}

/*
bool GetTime(std::string filename, std::vector<int64_t> &time, std::string &timestep, size_t nl) {
	filename += ".time";
	if (!file_exists(filename)) {
		return false;
	}
	std::vector<std::string> s = read_text(filename);
	if (nl != (s.size()-1)) return false;
	time.reserve(nl);
	timestep = s[0];
	for (size_t i=0; i<nl; i++) {
		time.push_back( parse_time(s[i+1]) );
	}
	return true;
}

bool GetUnits(std::string filename, std::vector<std::string> &units, size_t nl) {
	filename += ".unit";
	if (!file_exists(filename)) {
		return false;
	}
	units = read_text(filename);
	if (nl != units.size()) return false;
	return true;
}
*/

bool read_aux_json(std::string filename, std::vector<int64_t> &time, std::string &timestep, std::vector<std::string> &units, size_t nlyr) {
	filename += ".aux.json";
	if (!file_exists(filename)) return false;
	std::vector<std::string> s = read_text(filename);
	int itime=-1, istep=-1, iunit=-1;
	for (size_t i=0; i<s.size(); i++) {
		std::vector<std::string> x = strsplit_first(s[i], ":");
		if (x.size() != 2) continue;
		x[0].erase(std::remove(x[0].begin(), x[0].end(), '\"'), x[0].end());
		if (x[0] == "time") itime = i;
		if (x[0] == "timestep") istep = i;
		if (x[0] == "unit") iunit = i;
	}
	if (itime >= 0) {
		std::vector<std::string> x = strsplit_first(s[itime], "[");
		if (x.size() == 2) {
			x = strsplit(x[1], "]");
			x = strsplit(x[0], ",");
			std::vector<int64_t> tm;
			for (size_t i=0; i<x.size(); i++) {
				unquote(x[i]);
				tm.push_back( parse_time(x[i]) );
			}
			if (tm.size() == nlyr) {
				time = tm;
			}
		}
		if ((istep >= 0) && !time.empty()) {
			std::vector<std::string> x = strsplit_first(s[istep], ":");
			if (x.size() == 2) {
				x = strsplit(x[1], ",");
				unquote(x[0]);
				timestep = x[0];
			}
		}
	}
	if (iunit >= 0) {
		std::vector<std::string> x = strsplit_first(s[iunit], "[");
		if (x.size() == 2) {
			x = strsplit(x[1], "]");
			x = strsplit(x[0], ",");
			if (x.size() == nlyr) {
				for (size_t i=0; i< x.size(); i++) {
					unquote(x[i]);
				}
				units = x;
			}
		}
	}
	return false;
}


void prints(std::vector<std::string> x) {
	for (size_t i=0; i<x.size(); i++) {Rcpp::Rcout << x[i] << " ";}
	Rcpp::Rcout << "\n";	
}


bool GetRAT(GDALRasterAttributeTable *pRAT, SpatCategories &cats, const std::string &driver) {

/*
	const char *GFU_type_string[] = {"GFT_Integer", "GFT_Real","GFT_String"};
	const char *GFU_usage_string[] = {"GFU_Generic", "GFU_PixelCount", "GFU_Name", "GFU_Min",
		"GFU_Max", "GFU_MinMax", "GFU_Red", "GFU_Green", "GFU_Blue", "GFU_Alpha", "GFU_RedMin",
		"GFU_GreenMin", "GFU_BlueMin", "GFU_AlphaMin", "GFU_RedMax", "GFU_GreenMax", "GFU_BlueMax",
		"GFU_AlphaMax", "GFU_MaxCount"};
	std::vector<std::string> GFT_type;
	std::vector<std::string> GFT_usage;
*/

//	auto tabtype = pRAT->GetTableType(); perhaps check for "GRTT_ATHEMATIC" #1845
	
	size_t nc = (int) pRAT->GetColumnCount();
	size_t nr = (int) pRAT->GetRowCount();

	std::vector<std::string> ss = {"histogram", "count", "red", "green", "blue", "alpha", "opacity", "r", "g", "b", "a"};

	std::vector<std::string> ratnms;
	std::vector<int> id, id2;

	bool hasvalue=false;
	for (size_t i=0; i<nc; i++) {
		std::string name = pRAT->GetNameOfCol(i);
		ratnms.push_back(name);
		lowercase(name);
		if (!hasvalue && ((name == "value") || (name == "id") || (name == "ids"))) {
			id.insert(id.begin(), i);
			hasvalue = true;
		} else {
			int k = where_in_vector(name, ss, false);
			if (k >= 0) {
				id2.push_back(i);
			} else {
				id.push_back(i);
			}
		}
	}
	bool good_rat = true;
	size_t sid = id.size();
//	Rcpp::Rcout << hasvalue << " " << sid << std::endl;
	if ((hasvalue && sid == 1) || ((!hasvalue) && sid == 0)) {
// #790 avoid having just "count" or "histogram". return false for #1845
		good_rat = false;
		return false;
	}
	id.insert(id.end(), id2.begin(), id2.end());

	if (driver == "AIG") {
		std::vector<std::string> compnms = {"ID", "VALUE", "COUNT"};
		if ((id.size() == 3) && (ratnms == compnms)) {
			cats.index = -1;
			return false;
		}
	}

	if (!hasvalue) {
		std::vector<long> vid(nr);
		std::iota(vid.begin(), vid.end(), 0);
		cats.d.add_column(vid, "value");
	}

	int first_string = -1;

	for (size_t k=0; k<id.size(); k++) {
		size_t i = id[k];
		std::string name = pRAT->GetNameOfCol(i);
		GDALRATFieldType nc_type = pRAT->GetTypeOfCol(i);
//		GFT_type.push_back(GFU_type_string[nc_types[i]]);
//		GDALRATFieldUsage nc_usage = pRAT->GetUsageOfCol(i);
//		GFT_usage.push_back(GFU_usage_string[nc_usages[i]]);

		if (nc_type == GFT_Integer) {
			std::vector<long> d(nr);
			for (size_t j=0; j<nr; j++) {
				d[j] = (int) pRAT->GetValueAsInt(j, i);
			}
			cats.d.add_column(d, name);
		} else if (nc_type == GFT_Real) {
			std::vector<double> d(nr);
			for (size_t j=0; j<nr; j++) {
				d[j] = (double) pRAT->GetValueAsDouble(j, i);
			}
			cats.d.add_column(d, name);
		} else if (nc_type == GFT_String) {
			std::vector<std::string> d(nr);
			for (size_t j=0; j<nr; j++) {
				d[j] = (std::string) pRAT->GetValueAsString(j, i);
			}
			if (first_string < 0) first_string = cats.d.ncol();
			cats.d.add_column(d, name);
		}
	}	
	if (cats.d.nrow() == 0) {
		return false;
	}
	cats.index = good_rat ? (first_string >= 0 ? first_string :(cats.d.ncol() > 1 ? 1 : 0)) : -1;
	return true;
}


bool GetVAT(std::string filename, SpatCategories &vat) {

	filename += ".vat.dbf";
	if (!file_exists(filename)) {
		return false;
	}

	SpatVector v, fvct;
	std::vector<double> fext;

	v.read(filename, "", "", fext, fvct, false, "", "", {}); 
	if (v.df.nrow() == 0) return false;


	std::vector<std::string> nms = v.df.get_names();
	std::vector<std::string> ss = {"count", "histogram"};

	std::vector<size_t> rng;
	rng.reserve(nms.size());

	for (size_t i=0; i<nms.size(); i++) {
		int j = where_in_vector(nms[i], ss, true);
		if (j < 0) rng.push_back(i);
	}

	if (rng.size() > 1) {
		vat.d = v.df.subset_cols(rng);
//		vat.d.names[0] = "ID";
		vat.index = 1;
		std::string sc = vat.d.names[1];
		lowercase(sc);
		if (sc == "count") {
			if (rng.size() == 2) {
				return false;
			} else {
				vat.index = 2;
			}
		}
		return true;
	}
	return false;
}


SpatDataFrame GetCOLdf(GDALColorTable *pCT) {

	SpatDataFrame out;
	size_t nc = (int) pCT->GetColorEntryCount();

	out.add_column(1, "value");
	out.add_column(1, "red");
	out.add_column(1, "green");
	out.add_column(1, "blue");
	out.add_column(1, "alpha");
	out.reserve(nc);

	for (size_t i=0; i<nc; i++) {
		const GDALColorEntry * col = pCT->GetColorEntry(i);
		out.iv[0].push_back(i);
		out.iv[1].push_back(col->c1);
		out.iv[2].push_back(col->c2);
		out.iv[3].push_back(col->c3);
		out.iv[4].push_back(col->c4);
	}
	return(out);
}

bool getIntFromDoubleCol(std::vector<double> & dv, std::vector<long> &iv) {
	double dmn = vmin(dv, true);
	if (dmn < 0) return false;
	double dmx = vmax(dv, true);
	if (dmx > 255) {
		return false;
	}
	iv.resize(0);
	iv.reserve(dv.size());
	if (dmx <= 1) {
		for (size_t i=0; i<dv.size(); i++) {
			iv.push_back( dv[i] * 255 );
		}
	} else {
		for (size_t i=0; i<dv.size(); i++) {
			iv.push_back( dv[i] );
		}
	}
	return true;
}

bool setIntCol(SpatDataFrame &d, SpatDataFrame &out, int k, std::string name) {
	if (d.itype[k] == 0) {
		std::vector<long> iv;
		size_t j = d.iplace[k];
		if (getIntFromDoubleCol(d.dv[j], iv)) {
			out.add_column(iv, name);
		} else {
			return false;
		}
	} else if (d.itype[k] == 1) {
		size_t j = d.iplace[k];
		long dmn = vmin(d.iv[j], true);
		if (dmn < 0) return false;
		long dmx = vmax(d.iv[j], true);
		if (dmx > 255) return false;
		out.add_column(d.iv[j], name);
	} else {
		return false;
	}
	return true;
}


bool colsFromRat(SpatDataFrame &d, SpatDataFrame &out) {

	if ((d.nrow() == 0) || (d.ncol() == 0)) {
		return false;
	}

	std::vector<std::string> ss = d.get_names();
	for (size_t i=0; i<ss.size(); i++) {
		lowercase(ss[i]);
	}
//	int k = where_in_vector("value", ss, true);
//	if (k >= 0) {
	int k = 0;  
	size_t j = d.iplace[k];
		
	if (d.itype[k] == 1) {
		out.add_column(d.iv[j], "value");
	} else if (d.itype[k] == 0) {
		std::vector<long> x;
		x.reserve(d.nrow());
		for (size_t i=0; i<d.nrow(); i++) {
			x.push_back(d.dv[j][i]);
		}
		out.add_column(x, "value");
	} else {
		return false;
	}

	std::vector<std::string> cols1 = {"red", "green", "blue"};
	std::vector<std::string> cols2 = {"r", "g", "b"};
	for (size_t i=0; i<3; i++) {
		int k = where_in_vector(cols1[i], ss, true);
		if (k >= 0) {
			if (!setIntCol(d, out, k, cols1[i])) return false;
		} else {
			int k = where_in_vector(cols2[i], ss, true);
			if (k >= 0) {
				if (!setIntCol(d, out, k, cols1[i])) return false;
			} else {
				return false;
			}
		}
	}
	k = where_in_vector("alpha", ss, true);
	if (k >= 0) {
		setIntCol(d, out, k, "alpha");
	} else {
		int k = where_in_vector("transparency", ss, true);
		if (k >= 0) {
			setIntCol(d, out, k, "alpha");
		} else {
			int k = where_in_vector("opacity", ss, true);
			if (k >= 0) {
				setIntCol(d, out, k, "alpha");
			} else {
				std::vector<long> a(out.nrow(), 255);
				out.add_column(a, "alpha");
			}
		}
	}
	return true;
}


/*
SpatDataFrame GetColFromRAT(SpatDataFrame &rat) {

	SpatDataFrame out;
	size_t nr = rat.nrow();
	if (nr > 256) return out;
	std::vector<std::string> nms = rat.get_names();
	int red = where_in_vector("red", nms, false);
	int green = where_in_vector("green", nms, false);
	int blue = where_in_vector("blue", nms, false);
	int alpha = where_in_vector("alpha", nms, true);
	std::vector<unsigned> r {(unsigned)red, (unsigned)green, (unsigned)blue};
	if (alpha >= 0) {
		r.push_back(alpha);
	}
	out = rat.subset_cols(r);
	if (alpha < 0) {
		std::vector<long> a(nr, 255);
		out.add_column(a, "alpha");
	}
	out.names = {"red", "green", "blue", "alpha"};
	return out;
}
*/

SpatCategories GetCategories(char **pCat, std::string name) {
	long n = CSLCount(pCat);

	SpatCategories scat;

	std::vector<long> id;
	std::vector<std::string> nms;
	id.reserve(n);
	nms.reserve(n);

	for (long i = 0; i<n; i++) {
		const char *field = CSLGetField(pCat, i);
		std::string s = field;
		if (!s.empty()) {
			id.push_back(i);
			nms.push_back(field);
		}
	}

	scat.d.add_column(id, "value");
	name = name.empty() ? "category" : name;
	scat.d.add_column(nms, name);
	scat.index = 1;
	return(scat);
}


std::string strend(std::string f, size_t n) {
	n = std::min(n, f.length());
	std::string end = f.substr(f.length() - n);
	return end;
}

std::string basename_sds(std::string f) {
	const size_t i = f.find_last_of("\\/");
	if (std::string::npos != i) {
		f.erase(0, i + 1);
	}
	// this may be incorrect of the variable name includes a ":" ?
	const size_t j = f.find_last_of(':');
	if (std::string::npos != j) {
		f.erase(0, j + 1);
	}

	std::string end = strend(f, 3);
	if ((end == ".h5") || (end == ".nc")) {
		f.erase( f.end()-3, f.end() );
	} else if (strend(f, 4) == ".hdf")  {
		f.erase( f.end()-4, f.end() );
	}
	f.erase(std::remove(f.begin(), f.end(), '"'), f.end());

/*
	f = std::regex_replace(f, std::regex("\\.h5$"), "");
	f = std::regex_replace(f, std::regex("\\.hdf$"), "");
	f = std::regex_replace(f, std::regex("\\.nc$"), "");
	f = std::regex_replace(f, std::regex("\""), "");
*/
	return f;
}



std::string getDsWKT(GDALDataset *poDataset) {
	std::string wkt = "";
#if GDAL_VERSION_MAJOR >= 3
	const OGRSpatialReference *srs = poDataset->GetSpatialRef();
	if (srs == NULL) return wkt;
	char *cp;
	const char *options[3] = { "MULTILINE=YES", "FORMAT=WKT2", NULL };
	OGRErr err = srs->exportToWkt(&cp, options);
	if (err == OGRERR_NONE) {
		wkt = std::string(cp);
	}
	CPLFree(cp);

#else
	if (poDataset->GetProjectionRef() != NULL) {
		char *cp;
		OGRSpatialReference oSRS(poDataset->GetProjectionRef());

#if GDAL_VERSION_MAJOR >= 3
		const char *options[3] = { "MULTILINE=NO", "FORMAT=WKT2", NULL };
		OGRErr err = oSRS.exportToWkt(&cp, options);
#else
		OGRErr err = oSRS.exportToWkt(&cp);
#endif
		if (err == OGRERR_NONE) {
			wkt = std::string(cp);
		}
        CPLFree(cp);
	}
#endif
	return wkt;
}

std::string getDsPRJ(GDALDataset *poDataset) {
	std::string prj = "";
#if GDAL_VERSION_MAJOR >= 3
	const OGRSpatialReference *srs = poDataset->GetSpatialRef();
	if (srs == NULL) return prj;
	char *cp;
	OGRErr err = srs->exportToProj4(&cp);
	if (err == OGRERR_NONE) {
		prj = std::string(cp);
	}
        CPLFree(cp);
#else
	if( poDataset->GetProjectionRef() != NULL ) {
		OGRSpatialReference oSRS(poDataset->GetProjectionRef());
		char *pszPRJ = NULL;
		oSRS.exportToProj4(&pszPRJ);
		prj = pszPRJ;
	}
#endif
	return prj;
}


inline std::string dtypename(const std::string &d) {
	if (d == "Float64") return "FLT8S";
	if (d == "Float32") return "FLT4S";
	if (d == "Int64") return "INT8S";
	if (d == "Int32") return "INT4S";
	if (d == "Int16") return "INT2S";
	if (d == "Int8") return "INT1S";
	if (d == "UInt64") return "INT8U";
	if (d == "UInt32") return "INT4U";
	if (d == "UInt16") return "INT2U";
	if (d == "Byte") return "INT1U";
	return "FLT4S";
}


void get_tags(std::vector<std::string> meta, std::string prefix,  std::vector<std::string> &name, std::vector<std::string> &value) {
	if (!meta.empty()) {
		for (size_t i=0; i<meta.size(); i++) {
			size_t tagpos = meta[i].find(prefix);
			if (tagpos != std::string::npos) {		
				size_t pos = meta[i].find("=");
				if (pos != std::string::npos) {
					std::string mn = meta[i].substr(prefix.size(), pos-tagpos-prefix.size());
					if (!((mn == "_FillValue") || (mn == "grid_mapping") || (mn == "Conventions") || (mn == "created_by") || (mn == "created_date"))) {
						name.push_back(mn);
						value.push_back(meta[i].substr(pos+1, meta[i].length()));
					}
				}
			}
		}
	}
	
}

std::vector<std::string> get_metadata(std::string filename, std::vector<std::string> options) {
	std::vector<std::string> metadata;
    GDALDataset *poDataset = openGDAL(filename, GDAL_OF_RASTER | GDAL_OF_READONLY | GDAL_OF_VERBOSE_ERROR, {}, options);
    if( poDataset == NULL )  {
		return metadata;
	}
	char **m = poDataset->GetMetadata();
	if (m) {
		while (*m != nullptr) {
			metadata.push_back(*m++);
		}
	}
	GDALClose( (GDALDatasetH) poDataset );
	return metadata;
}


SpatRasterStack::SpatRasterStack(std::string fname, std::vector<int> ids, bool useids, std::vector<std::string> options, bool noflip, bool guessCRS, std::vector<std::string> domains) {

    GDALDataset *poDataset = openGDAL(fname, GDAL_OF_RASTER | GDAL_OF_READONLY | GDAL_OF_VERBOSE_ERROR, {}, {});
    if( poDataset == NULL )  {
		if (!file_exists(fname)) {
			setError("file does not exist: " + fname);
		} else {
			setError("cannot read from " + fname );
		}
		return;
	}

	std::string delim = "NAME=";
	char **metadata = poDataset->GetMetadata("SUBDATASETS");

	if (metadata == NULL) {
		GDALClose( (GDALDatasetH) poDataset );
		SpatRaster sub;
		if (sub.constructFromFile(fname, {-1}, {""}, {}, options, false, guessCRS, domains)) {
			std::string sname = sub.source[0].source_name;
			push_back(sub, sname, sub.source[0].source_name_long, sub.source[0].unit[0], true);
		}
		return;
	}

	std::vector<std::string> meta;
    for (size_t i=0; metadata[i] != NULL; i++) {
		meta.push_back(metadata[i]);
	}

	if (!useids) {
		ids.resize(meta.size());
		std::iota(ids.begin(), ids.end(), 0);
	}
	int idssz = ids.size();
	int metsz = meta.size();

	if (metsz == 0) {
		setError("file does not consist of subdatasets");
	} else {
		for (int i=0; i<idssz; i++) {
			if ((ids[i] < 0) || ((2*ids[i]) >= metsz)) {
				continue;
			}
			std::string s = meta[ids[i]*2];
			size_t pos = s.find(delim);
			if (pos != std::string::npos) {
				s.erase(0, pos + delim.length());
				SpatRaster sub;
				if (sub.constructFromFile(s, {-1}, {""}, {}, options, false, guessCRS, domains)) {
					std::string sname = sub.source[0].source_name.empty() ? basename_sds(s) : sub.source[0].source_name;
					if (!push_back(sub, sname, sub.source[0].source_name_long, sub.source[0].unit[0], true)) {
						addWarning("skipped (different geometry): " + s);
					}
				} else {
					addWarning("skipped (fail): " + s);
				}
			}
		}
	}
	meta.resize(0);
	char **m = poDataset->GetMetadata();
	if (m) {
		while (*m != nullptr) {
			meta.push_back(*m++);
		}
	}
	GDALClose( (GDALDatasetH) poDataset );

	std::vector<std::string> tagnames, tagvalues;
//	get_tags(meta, "NC_GLOBAL#TAG_", tagnames, tagvalues);
	get_tags(meta, "NC_GLOBAL#", tagnames, tagvalues);
	for (size_t i=0; i<tagnames.size(); i++) addTag(tagnames[i], tagvalues[i], "GLOBAL");

}


SpatRasterCollection::SpatRasterCollection(std::string fname, std::vector<int> ids, bool useids, std::vector<std::string> options, bool noflip, bool guessCRS, std::vector<std::string> domains) {

//	std::vector<std::string> ops;
    GDALDataset *poDataset = openGDAL(fname, GDAL_OF_RASTER | GDAL_OF_READONLY | GDAL_OF_VERBOSE_ERROR, {}, {});
    if( poDataset == NULL )  {
		if (!file_exists(fname)) {
			setError("file does not exist: " + fname);
		} else {
			setError("cannot read from " + fname );
		}
		return;
	}

	std::string delim = "NAME=";
	char **metadata = poDataset->GetMetadata("SUBDATASETS");

	if (metadata == NULL) {
		GDALClose( (GDALDatasetH) poDataset );
		SpatRaster sub;
		if (sub.constructFromFile(fname, {-1}, {""}, {}, options, false, guessCRS, domains)) {
			std::string sname = sub.source[0].source_name;
			push_back(sub, sname);
		}
		return;
	}

	std::vector<std::string> meta;
    for (size_t i=0; metadata[i] != NULL; i++) {
		meta.push_back(metadata[i]);
	}

	if (!useids) {
		ids.resize(meta.size());
		std::iota(ids.begin(), ids.end(), 0);
	}
	int idssz = ids.size();
	int metsz = meta.size();

	if (metsz == 0) {
		setError("file does not consist of subdatasets");
	} else {
		for (int i=0; i<idssz; i++) {
			if ((ids[i] < 0) || ((2*ids[i]) >= metsz)) {
				continue;
			}
			std::string s = meta[ids[i]*2];
			size_t pos = s.find(delim);
			if (pos != std::string::npos) {
				s.erase(0, pos + delim.length());
				SpatRaster sub;
				if (sub.constructFromFile(s, {-1}, {""}, {}, options, false, guessCRS, domains)) {
					push_back(sub, basename_sds(s));
				} else {
					addWarning("skipped (fail): " + s);
				}
			}
		}
	}
	meta.resize(0);
	char **m = poDataset->GetMetadata();
	if (m) {
		while (*m != nullptr) {
			meta.push_back(*m++);
		}
	}
	GDALClose( (GDALDatasetH) poDataset );

	std::vector<std::string> tagnames, tagvalues;
//	get_tags(meta, "NC_GLOBAL#TAG_", tagnames, tagvalues);
	get_tags(meta, "NC_GLOBAL#", tagnames, tagvalues);
	for (size_t i=0; i<tagnames.size(); i++) addTag(tagnames[i], tagvalues[i], "GLOBAL");

}


/*
SpatRaster SpatRaster::fromFiles(std::vector<std::string> fname, std::vector<int> subds, std::vector<std::string> subdsname, std::vector<std::string> drivers, std::vector<std::string> options) {
	SpatRaster out;
	out.constructFromFile(fname[0], subds, subdsname, options);
	if (out.hasError()) return out;
	SpatOptions opt;
	for (size_t i=1; i<fname.size(); i++) {
		SpatRaster r;
		bool ok = r.constructFromFile(fname[i], subds, subdsname, drivers, options);
		if (r.msg.has_warning) {
			out.addWarning(r.msg.warnings[0]);
		}
		if (ok) {
			out.addSource(r, false, opt);
			if (r.msg.has_error) {
				out.setError(r.msg.error);
				return out;
			}
		} else {
			if (r.msg.has_error) {
				out.setError(r.msg.error);
			}
			return out;
		}
	}
	return out;
}

*/

bool getGCPs(GDALDataset *poDataset, SpatRasterSource &s) {
	int n = poDataset->GetGCPCount();
//	Rcpp::Rcout << "n GCP " << n << std::endl;
	if (n == 0) return false;
	const GDAL_GCP *gcp;
	gcp	= poDataset->GetGCPs();
	
	double adfGeoTransform[6];
	if (GDALGCPsToGeoTransform(n, gcp, adfGeoTransform, true)) {
		//for (size_t i=0; i<6; i++) {
		//	Rcpp::Rcout << adfGeoTransform[i] << " ";
		//}
		//Rcpp::Rcout << std::endl;
		double xmin = adfGeoTransform[0]; /* left x */
		double xmax = xmin + adfGeoTransform[1] * s.ncol; /* w-e resolution */
		if (xmin > xmax) {
			std::swap(xmin, xmax);			
		}
		double ymax = adfGeoTransform[3]; // top y
		double ymin = ymax + s.nrow * adfGeoTransform[5];
		if (adfGeoTransform[5] > 0) {
			s.flipped = true;
			std::swap(ymin, ymax);
		}
		SpatExtent e(xmin, xmax, ymin, ymax);
		s.extent = e;
		if (adfGeoTransform[2] != 0 || adfGeoTransform[4] != 0) {
			s.rotated = true;
		}
		return true;
	}
	return false;
}



bool SpatRaster::constructFromFile(std::string fname, std::vector<int> subds, std::vector<std::string> subdsname, std::vector<std::string> drivers, std::vector<std::string> options, bool noflip, bool guessCRS, std::vector<std::string> domains) {



	if (fname == "WCS:") {
		// for https://github.com/rspatial/terra/issues/1505
		setError("no raster data in WCS:");
		return false;
	}
	
	bool apply_so = true;
	std::vector<std::string> clean_ops = options;
	size_t opsz = options.size();
	if (opsz > 0) {
		if (options[opsz-1] == "so=false") {
			apply_so = false;
			clean_ops.resize(opsz-1);
		}
	}

    GDALDataset *poDataset = openGDAL(fname, GDAL_OF_RASTER | GDAL_OF_READONLY | GDAL_OF_VERBOSE_ERROR, drivers, clean_ops);

    if( poDataset == NULL )  {
		if (!file_exists(fname)) {
			setError("file does not exist: " + fname);
		} else {
			setError("cannot open this file as a SpatRaster: " + fname);
		}
		return false;
	}

	int nl = poDataset->GetRasterCount();
	std::string gdrv = poDataset->GetDriver()->GetDescription();

	char **metasds = poDataset->GetMetadata("SUBDATASETS");
	
	if (metasds != NULL) {
		std::vector<std::string> meta;
		for (size_t i=0; metasds[i] != NULL; i++) {
			meta.push_back(metasds[i]);
		}
		GDALClose( (GDALDatasetH) poDataset );
		return constructFromSDS(fname, meta, subds, subdsname, options, gdrv, noflip, guessCRS, domains);

	} else if (nl==0) {
		setError("no raster data in " + fname);
		return false;
	}

	for (size_t i=0; i<domains.size(); i++) {
		char **meterra = poDataset->GetMetadata(domains[i].c_str());
		if (meterra != NULL) {
			std::vector<std::string> meta;
			for (size_t j=0; meterra[j] != NULL; j++) {
				std::string s = meterra[j];
				size_t pos = s.find("=");
				if (pos != std::string::npos) {
					std::string name = s.substr(0, pos);
					std::string value = s.substr(pos+1); 
					addTag(name, value, domains[i]);
				}
			}
		}
	}
	
	SpatRasterSource s;

	char **metasrc = poDataset->GetMetadata();
	while (metasrc != nullptr && *metasrc != nullptr) {
		s.smdata.push_back(*metasrc++);
	}
	s.ncol = poDataset->GetRasterXSize();
	s.nrow = poDataset->GetRasterYSize();
	s.nlyr = nl;
	s.nlyrfile = nl;
	s.resize(nl);

	s.flipped = false;
	s.rotated = false;
	double adfGeoTransform[6];

	bool hasExtent = true;
	if (poDataset->GetGeoTransform( adfGeoTransform ) == CE_None) {

		double xmin = adfGeoTransform[0]; /* left x */
		double xmax = xmin + adfGeoTransform[1] * s.ncol; /* w-e resolution */
		//xmax = roundn(xmax, 9);
		double ymax = adfGeoTransform[3]; // top y
		double ymin = ymax + s.nrow * adfGeoTransform[5];
		//ymin = roundn(ymin, 9);

		if (adfGeoTransform[5] > 0) {
			std::swap(ymin, ymax);
			s.extset = true;
			s.flipped = true;
		}

		SpatExtent e(xmin, xmax, ymin, ymax);
		s.extent = e;

		if (adfGeoTransform[2] != 0 || adfGeoTransform[4] != 0) {
			s.rotated = true;
			addWarning("the data in this file are rotated. Use 'rectify' to fix that");
		}
	} else if (getGCPs(poDataset, s)) {
		if (s.rotated) {
			addWarning("the data in this file are rotated. Use 'rectify' to fix that");
		}
	} else {
		bool warn=true;
		hasExtent = false;
		if (adfGeoTransform[5] > 0) {
			if (noflip) {
				warn = false;
				s.extset = true;
			} else {
				s.flipped = true;
			}
		}
		SpatExtent e(0, s.ncol, 0, s.nrow);
		s.extent = e;

		if ((gdrv=="netCDF") || (gdrv == "HDF5")) {
			#ifndef standalone
			setMessage("ncdf extent");
			#else
			addWarning("unknown extent. Cells not equally spaced?");
			#endif
		} else if (warn) {
			addWarning("unknown extent");
		}

		// seems to cause more harm then benefit #1627
		//try {
		//	s.flipped = adfGeoTransform[5] > 0;
		//} catch(...) {}
	}

	s.memory = false;
	s.filename = fname;
	s.open_ops = clean_ops;
	//s.open_drivers = {gdrv}; // failed for some hdf
	s.open_drivers = drivers;
	//s.driver = "gdal";

/*
	if( poDataset->GetProjectionRef() != NULL ) {
		OGRSpatialReference oSRS(poDataset->GetProjectionRef());
		char *pszPRJ = NULL;
		oSRS.exportToProj4(&pszPRJ);
		s.crs = pszPRJ;
	} else {
		s.crs = "";
	}
*/
	std::string crs = getDsWKT(poDataset);

	if (crs.empty()) {
		if (guessCRS && hasExtent && s.extent.xmin >= -180.1 && s.extent.xmax <= 360.1 && s.extent.ymin >= -90.1 && s.extent.ymax <= 90.1) {
			crs = "OGC:CRS84";
			s.parameters_changed = true;
		} else {
			crs = "";
		}
	}

	std::string msg;
	if (!s.srs.set(crs, msg)) {
		addWarning(msg);
	}

	GDALRasterBand  *poBand;
	//int nBlockXSize, nBlockYSize;
	double adfMinMax[2];
	int bGotMin, bGotMax;

//	s.layers.resize(1);
//	std::string unit = "";

	s.source_name = basename_noext(fname);
	std::vector<std::vector<std::string>> bandmeta(s.nlyr);
	bool getCols = s.nlyr == 3;
	std::vector<unsigned> rgb_lyrs(3, -99);

	s.hasUnit = true;
	s.hasTime = true;
	std::vector<std::string> datm, unts;
	datm.reserve(s.nlyr);
	unts.reserve(s.nlyr);
	
	int bs1, bs2;
	for (size_t i = 0; i < s.nlyr; i++) {
		
		poBand = poDataset->GetRasterBand(i+1);

		if (s.hasTime) {
			const char* dtm = poBand->GetMetadataItem("DATE_TIME");
			if (dtm != NULL) {
				datm.push_back(dtm);
			} else {
				s.hasTime = false;
			}
		}
		if (s.hasUnit) {
			const char* ut = poBand->GetMetadataItem("UNIT");
			if (ut != NULL) {
				unts.push_back(ut);
			} else {
				s.hasUnit = false;
			}
		}

		// if ((gdrv=="netCDF") || (gdrv == "HDF5") || (gdrv == "GRIB") || (gdrv == "GTiff")) {
			char **m = poBand->GetMetadata();
			while (m != nullptr && *m != nullptr) {
				bandmeta[i].push_back(*m++);
			}
			//for (size_t j = 0; j<bandmeta[i].size(); j++) {
			//	Rcpp::Rcout << bandmeta[i][j] << std::endl;
			//}
			
			char **meterra = poBand->GetMetadata("USER_TAGS");
			if (meterra != NULL) {
//				std::vector<std::string> meta;
				for (size_t j=0; meterra[j] != NULL; j++) {
					std::string ms = meterra[j];
					size_t pos = ms.find("=");
					if (pos != std::string::npos) {
						std::string name = ms.substr(0, pos);
						std::string value = ms.substr(pos+1); 
						s.addLyrTag(i, name, value);
					}
				}
			}
		// }

		int success;
	//	double naflag = poBand->GetNoDataValue(&success);
	//	if (success) {
	//		s.NAflag = naflag;
	//	} else {
	//		s.NAflag = NAN;
	//	}

		adfMinMax[0] = poBand->GetMinimum( &bGotMin );
		adfMinMax[1] = poBand->GetMaximum( &bGotMax );

		s.has_scale_offset[i] = false;

		if (apply_so) {
			double offset = poBand->GetOffset(&success);
			if (success) {
				if (offset != 0) {
					s.offset[i] = offset;
					s.has_scale_offset[i] = true;
				}
			}
			double scale = poBand->GetScale(&success);
			if (success) {
				if (scale != 1) {
					s.scale[i] = scale;
					s.has_scale_offset[i] = true;
				}
			}
			if (s.has_scale_offset[i]) {
				adfMinMax[0] = adfMinMax[0] * s.scale[i] + s.offset[i];
				adfMinMax[1] = adfMinMax[1] * s.scale[i] + s.offset[i];
			}
		}

		if( (bGotMin && bGotMax) ) {
			s.hasRange[i] = true;
			s.range_min[i] = adfMinMax[0];
			s.range_max[i] = adfMinMax[1];
		}

		poBand->GetBlockSize(&bs1, &bs2);
		s.blockcols[i] = bs1;
		s.blockrows[i] = bs2;
		s.dtype = dtypename(GDALGetDataTypeName(poBand->GetRasterDataType()));


		//if( poBand->GetOverviewCount() > 0 ) printf( "Band has %d overviews.\n", poBand->GetOverviewCount() );

		if (getCols) {
			if (poBand->GetColorInterpretation() == GCI_RedBand) {
				rgb_lyrs[0] = i;
			} else if (poBand->GetColorInterpretation() == GCI_GreenBand) {
				rgb_lyrs[1] = i;
			} else if (poBand->GetColorInterpretation() == GCI_BlueBand) {
				rgb_lyrs[2] = i;
			}
		}
		GDALColorTable *ct = poBand->GetColorTable();
		if( ct != NULL ) {
			s.hasColors[i] = true;
			s.cols[i] = GetCOLdf(ct);
		}

		std::string bandname = poBand->GetDescription();

		char **cat = poBand->GetCategoryNames();
		if (cat != NULL)	{
			SpatCategories scat = GetCategories(cat, bandname);

			s.cats[i] = scat;
			s.hasCategories[i] = true;
		}

		SpatCategories crat;
		//bool found_rat = false;
		
		if (!s.hasCategories[i]) {
			GDALRasterAttributeTable *rat = poBand->GetDefaultRAT();
			if (rat != NULL) {
				if (GetRAT(rat, crat, gdrv)) {
					s.cats[i] = crat;
					s.hasCategories[i] = true;
//				} else {
//					found_rat = false;
				}
			}
		}
		//	} else {
		//		s.cats[i].d.cbind(crat.d); // needs more checking.
		//	} else {

		if (!s.hasCategories[i]) {
			if (GetVAT(fname, crat)) {
				s.cats[i] = crat;
				s.hasCategories[i] = true;
				//found_rat = true;
			}
		}

		if ((!s.hasColors[i]) && s.hasCategories[i]) {
			SpatDataFrame ratcols;
			if (colsFromRat(crat.d, ratcols)) {
				s.hasColors[i] = true;
				s.cols[i] = ratcols;
			}
		}

		std::string nm = "";
		if (s.hasCategories[i]) {
			if ((s.cats[i].index >= 0) && (s.cats[i].index < (int)s.cats[i].d.ncol())) {
				std::vector<std::string> nms = s.cats[i].d.get_names();
				nm = nms[s.cats[i].index];
			}
		}
		
		if (nm.empty()) {
			if (!bandname.empty()) {
				nm = bandname;
			} else if (s.nlyr > 1) {
				nm = s.source_name + "_" + std::to_string(i+1);
			} else {
				nm = basename_noext(fname) ;
			}
		}

		std::string dtype = GDALGetDataTypeName(poBand->GetRasterDataType());
		if ((!s.has_scale_offset[i]) && (in_string(dtype, "Int") || (dtype == "Byte"))) {
			s.valueType[i] = 1;
		}
		s.names[i] = nm;
	}

	if (s.hasTime) {
			
		if (datm[0].find('T') != std::string::npos) {
			s.timestep = "seconds";
		} else {
			// backwards compatibility
			if (datm[0].length() == 4) {
				s.timestep = "years";
			} else if (datm[0].length() == 7) {
				if (datm[0].substr(0, 5) == "0000-") {
					s.timestep = "months";					
				} else {
					s.timestep = "yearmonths";
				}
			} else if (datm[0].length() == 10) {
			// current formats, always xxxx-xx-xx where possible
				if (datm[0].substr(7, 3) == "-00") {
					if (datm[0].substr(4, 3) == "-00") {
						s.timestep = "years";
					} else if (datm[0].substr(0, 4) == "0000") {
						s.timestep = "months";		
					} else {
						s.timestep = "yearmonths";
					}
				} else {
					s.timestep = "days";
				}
			}
		}
		for (size_t i=0; i<datm.size(); i++) {
			s.time[i] = parse_time(datm[i]);
		}
		

// try units from json
		std::vector<int64_t> timestamps;
		std::string timestep="raw";
		//std::vector<std::string> units;
		try {
			read_aux_json(fname, timestamps, timestep, unts, s.nlyr);
		} catch(...) {
			unts.resize(0);
			addWarning("could not parse aux.json");
		}
		if (!unts.empty()) {
			s.hasUnit = true;
		}
		
		
	} else {

		std::vector<int64_t> timestamps;
		std::string timestep="raw";
//		std::vector<std::string> units;
		if (unts.empty()) {
			try {
				read_aux_json(fname, timestamps, timestep, unts, s.nlyr);
			} catch(...) {
				timestamps.resize(0);
				unts.resize(0);
				addWarning("could not parse aux.json");
			}
			if (!timestamps.empty()) {
				s.time = timestamps;
				s.timestep = timestep;
				s.hasTime = true;
			}
			if (!unts.empty()) {
				s.hasUnit = true;
			}
		}
	}
	if (s.hasUnit) {
		s.unit = unts;
	}


	msg = "";
	std::vector<std::string> metadata;

	if ((gdrv=="netCDF") || (gdrv == "HDF5"))  {
		
		char **m = poDataset->GetMetadata();
		if (m) {
			while (*m != nullptr) {
				metadata.push_back(*m++);
			}
		}
		s.set_names_time_ncdf(metadata, bandmeta, msg);

		if (s.srs.is_empty() && guessCRS) {

			bool lat = false;
			bool lon = false;
			for (size_t i=0; i<metadata.size(); i++) {				
				if (!lat) lat = metadata[i].find("long_name=latitude") != std::string::npos;
				if (!lon) lon = metadata[i].find("long_name=longitude") != std::string::npos;
			}
			if (lon && lat && s.extent.ymin > -91 && s.extent.ymax < 91 && s.extent.xmin > -361  && s.extent.xmax < 361) {
				if (s.srs.set("OGC:CRS84", msg)) {
					s.parameters_changed = true;
				}
			}
		}
	} else if (gdrv == "GRIB") {	
		s.set_names_time_grib(bandmeta, msg);
	} else if (gdrv == "GTiff") {	
	// needs to get its own generic one 
	//	s.set_names_time_tif(bandmeta, msg);
	}
	s.bmdata = bandmeta;
	if (msg.size() > 1) {
		addWarning(msg);
	}

	GDALClose( (GDALDatasetH) poDataset );
	s.hasValues = true;
	setSource(s);

	if ((!metadata.empty())) {
		std::vector<std::string> tagnames, tagvalues;
//		std::string stag = s.source_name + "#TAG_";
		std::string stag = s.source_name + "#";
		get_tags(metadata, stag, tagnames, tagvalues);
		for (size_t i=0; i<tagnames.size(); i++) addTag(tagnames[i], tagvalues[i], gdrv);

	}
	if (getCols) {
		setRGB(rgb_lyrs[0], rgb_lyrs[1], rgb_lyrs[2], -99, "rgb");
	}

	return true;
}


bool SpatRaster::readStartGDAL(size_t src) {

    GDALDataset *poDataset = openGDAL(source[src].filename, GDAL_OF_RASTER | GDAL_OF_READONLY, source[src].open_drivers, source[src].open_ops);


    if( poDataset == NULL )  {
		size_t ncolon = std::count(source[src].filename.begin(), source[src].filename.end(), ':');
		if ((ncolon < 2) && (!file_exists(source[src].filename ))) {
			setError("file does not exist: " + source[src].filename);
		} else {
			if (source[src].filename.substr(0, 4) == "HDF4") {
				setError("cannot read from " + source[src].filename + "\n(Only 32 open datasets allowed with HDF4)");
			} else {
				setError("cannot read from " + source[src].filename);
			}
		}
		return false;
	}

    source[src].gdalconnection = poDataset;
	source[src].open_read = true;
	return(true);
}

bool SpatRaster::readStopGDAL(size_t src) {
	if (source[src].gdalconnection != NULL) {
		GDALClose( (GDALDatasetH) source[src].gdalconnection);
	}
	source[src].open_read = false;
	return true;
}



void NAso(std::vector<double> &d, size_t n, const std::vector<double> &flags, const std::vector<double> &scale, const std::vector<double>  &offset, const std::vector<bool> &haveso, const bool haveUserNAflag, const double userNAflag){
	size_t nl = flags.size();
	double na = NAN;

	for (size_t i=0; i<nl; i++) {
		size_t start = i*n;
		if (!std::isnan(flags[i])) {
			double flag = flags[i];
			// a hack to avoid problems with double derived from float - double comparison
			if (flag < -3.4e+37) {
				flag = -3.4e+37;
				for (size_t j=start; j<(start+n); j++) {
					if (d[j] < flag) {
						d[j] = NAN;
					}
				}
			} else {
				std::replace(d.begin()+start, d.begin()+start+n, flag, na);
			}
		}
		
		if (haveso[i]) {
			for (size_t j=start; j<(start+n); j++) {
				d[j] = d[j] * scale[i] + offset[i];
			}
		}
	}
	if (haveUserNAflag) {
		std::replace(d.begin(), d.end(), userNAflag, na);
	}
}


void vflip(std::vector<double> &v, const size_t &ncell, const size_t &nrows, const size_t &ncols, const size_t &nl) {
	for (size_t i=0; i<nl; i++) {
		size_t off = i*ncell;
		size_t nr = nrows/2;
		for (size_t j=0; j<nr; j++) {
			size_t d1 = off + j * ncols;
			size_t d2 = off + (nrows-j-1) * ncols;
			std::vector<double> r(v.begin()+d1, v.begin()+d1+ncols);
			std::copy(v.begin()+d2, v.begin()+d2+ncols, v.begin()+d1);
			std::copy(r.begin(), r.end(), v.begin()+d2);
		}
	}
}


void SpatRaster::readChunkGDAL(std::vector<double> &data, size_t src, size_t row, size_t nrows, size_t col, size_t ncols) {

	if (source[src].is_multidim) {
		readChunkMulti(data, src, row, nrows, col, ncols);
		return;
	}
	std::vector<double> errout;
	if (source[src].rotated) {
		setError("cannot read from rotated files. First use 'rectify'");
		return;
	}

	if (!(source[src].open_read || source[src].open_write)) {
		setError("the file is not open for reading");
		return;
	}


	if (source[src].flipped) {
		row = nrow() - row - nrows;
	}


	if (source[src].hasWindow) { // ignoring the expanded case.
		row = row + source[src].window.off_row;
		col = col + source[src].window.off_col;
	}

	size_t ncell = ncols * nrows;
	size_t nl = source[src].nlyr;
	std::vector<double> out(ncell * nl);
	int hasNA;
	std::vector<double> naflags(nl, NAN);
	CPLErr err = CE_None;

	std::vector<int> panBandMap;
	if (!source[src].in_order(true)) {
		panBandMap.reserve(nl);
		for (size_t i=0; i < nl; i++) {
			panBandMap.push_back(source[src].layers[i]+1);
		}
	}

	if (panBandMap.empty()) {
		err = source[src].gdalconnection->RasterIO(GF_Read, col, row, ncols, nrows, &out[0], ncols, nrows, GDT_Float64, nl, NULL, 0, 0, 0, NULL);
	} else {
		err = source[src].gdalconnection->RasterIO(GF_Read, col, row, ncols, nrows, &out[0], ncols, nrows, GDT_Float64, nl, &panBandMap[0], 0, 0, 0, NULL);
	}

	GDALRasterBand  *poBand;
	if (err == CE_None ) {
		for (size_t i=0; i<nl; i++) {
			poBand = source[src].gdalconnection->GetRasterBand(source[src].layers[i]+1);
			double naflag = poBand->GetNoDataValue(&hasNA);
			if (hasNA)  naflags[i] = naflag;
		}
		NAso(out, ncell, naflags, source[src].scale, source[src].offset, source[src].has_scale_offset, source[src].hasNAflag, source[src].NAflag);
	}

/*
	for (size_t i=0; i < nl; i++) {
		cell = ncell * i;
		poBand = source[src].gdalconnection->GetRasterBand(source[src].layers[i] + 1);
		double naflag = poBand->GetNoDataValue(&hasNA);
		if (!hasNA) { naflag = NAN; }
		GDALDataType gdtype = poBand->GetRasterDataType();
		if (gdtype == GDT_Float64) {
			err = poBand->RasterIO(GF_Read, col, row, ncols, nrows, &out[cell], ncols, nrows, gdtype, 0, 0);
			if (err != CE_None ) { break; }
			set_NA(out, naflag);
		}
	}
*/
	if (err != CE_None ) {
		setError("cannot read values");
		return;
	}

	if (source[src].flipped) {
		vflip(out, ncell, nrows, ncols, nl);
	}
	data.insert(data.end(), out.begin(), out.end());
}




std::vector<double> SpatRaster::readValuesGDAL(size_t src, size_t row, size_t nrows, size_t col, size_t ncols, int lyr) {

	if (source[src].is_multidim) {
		return readValuesMulti(src, row, nrows, col, ncols, lyr);
	}

	std::vector<double> errout;
	if (source[src].rotated) {
		setError("cannot read from rotated files. First use 'rectify'");
		return errout;
	}
	if (source[src].flipped) {
		row = nrow() - row - nrows;
	}

	if (source[src].hasWindow) { // ignoring the expanded case.
		row = row + source[src].window.off_row;
		col = col + source[src].window.off_col;
	}

    GDALDataset *poDataset = openGDAL(source[src].filename, GDAL_OF_RASTER | GDAL_OF_READONLY, source[src].open_drivers, source[src].open_ops);
	
    if( poDataset == NULL )  {
		if (!file_exists(source[src].filename )) {
			setError("file does not exist: " + source[src].filename);
		} else {
			setError("cannot read from " + source[src].filename  );
		}
		return errout;
	}

	GDALRasterBand *poBand;

	unsigned ncell = ncols * nrows;
	unsigned nl;
	std::vector<int> panBandMap;
	if (lyr < 0) {
		nl = source[src].nlyr;
		if (!source[src].in_order(true)) {
			panBandMap.reserve(nl);
			for (size_t i=0; i < nl; i++) {
				panBandMap.push_back(source[src].layers[i]+1);
			}
		}
	} else {
		nl = 1;
		panBandMap.push_back(source[src].layers[lyr]+1);
	}

	std::vector<double> out(ncell*nl);
	int hasNA;
	std::vector<double> naflags(nl, NAN);
	CPLErr err = CE_None;
	if (panBandMap.empty()) {
		err = poDataset->RasterIO(GF_Read, col, row, ncols, nrows, &out[0], ncols, nrows, GDT_Float64, nl, NULL, 0, 0, 0, NULL);
	} else {
		err = poDataset->RasterIO(GF_Read, col, row, ncols, nrows, &out[0], ncols, nrows, GDT_Float64, nl, &panBandMap[0], 0, 0, 0, NULL);
	}

	if (err == CE_None ) {
		for (size_t i=0; i<nl; i++) {
			poBand = poDataset->GetRasterBand(source[src].layers[i]+1);
			double naf = poBand->GetNoDataValue(&hasNA);
			if (hasNA)  naflags[i] = naf;
		}
		NAso(out, ncell, naflags, source[src].scale, source[src].offset, source[src].has_scale_offset, source[src].hasNAflag, source[src].NAflag);
	}

	GDALClose((GDALDatasetH) poDataset);
	if (err != CE_None ) {
		setError("cannot read values");
		return errout;
	}

	if (source[src].flipped) {
		vflip(out, ncell, nrows, ncols, nl);
	}
	return out;
}



std::vector<double> SpatRaster::readGDALsample(size_t src, size_t srows, size_t scols, bool overview) {

	if (source[src].is_multidim) {
		return readSampleMulti(src, srows, scols, overview);
	}

	std::vector<double> errout;
	if (source[src].rotated) {
		setError("cannot read from rotated files. First use 'rectify'");
		return errout;
	}

	size_t row =0, col=0, nrows=nrow(), ncols=ncol();
	if (source[src].hasWindow) {
		row = row + source[0].window.off_row;
		col = col + source[0].window.off_col;
		srows = std::min(srows, nrows);
		scols = std::min(scols, ncols);
	}

	std::vector<std::string> openops = source[src].open_ops;
	
	#if GDAL_VERSION_MAJOR <= 3 && GDAL_VERSION_MINOR < 3
	// do nothing
	#else 
	if (!overview) {
		openops.push_back("OVERVIEW_LEVEL=NONE");
	}
	#endif
	
    GDALDataset *poDataset = openGDAL(source[src].filename, GDAL_OF_RASTER | GDAL_OF_READONLY, source[src].open_drivers, openops);

    if( poDataset == NULL )  {
		if (!file_exists(source[src].filename )) {
			setError("file does not exist: " + source[src].filename);
		} else {
			setError("cannot read from " + source[src].filename  );
		}
		return errout;
	}

	size_t ncell = scols * srows;
	size_t nl = source[src].nlyr;
	std::vector<double> out(ncell*nl);
	int hasNA;
	CPLErr err = CE_None;

	std::vector<double> naflags(nl, NAN);

	std::vector<int> panBandMap;
	if (!source[src].in_order(true)) {
		panBandMap.reserve(nl);
		for (size_t i=0; i < nl; i++) {
			panBandMap.push_back(source[src].layers[i]+1);
		}
	}
/*
	if (panBandMap.size() > 0) {
		err = poDataset->RasterIO(GF_Read, col, row, ncols, nrows, &out[0], scols, srows, GDT_Float64, nl, &panBandMap[0], 0, 0, 0, NULL);
	} else {
		err = poDataset->RasterIO(GF_Read, col, row, ncols, nrows, &out[0], scols, srows, GDT_Float64, nl, NULL, 0, 0, 0, NULL);
	}
*/

	if (panBandMap.empty()) {
		err = poDataset->RasterIO(GF_Read, col, row, ncols, nrows, &out[0], scols, srows, GDT_Float64, nl, NULL, 0, 0, 0, NULL);
	} else {
		err = poDataset->RasterIO(GF_Read, col, row, ncols, nrows, &out[0], scols, srows, GDT_Float64, nl, &panBandMap[0], 0, 0, 0, NULL);
	}

	if (err == CE_None ) {
		GDALRasterBand *poBand;
		for (size_t i=0; i<nl; i++) {
			poBand = poDataset->GetRasterBand(source[src].layers[i]+1);
			double naflag = poBand->GetNoDataValue(&hasNA);
			if (hasNA)  naflags[i] = naflag;
		}
		NAso(out, ncell, naflags, source[src].scale, source[src].offset, source[src].has_scale_offset, source[src].hasNAflag, source[src].NAflag);
	}


/*
	for (size_t i=0; i < nl; i++) {
		poBand = poDataset->GetRasterBand(source[src].layers[i] + 1);
		size_t off = i * ncell;
		err = poBand->RasterIO(GF_Read, 0, 0, ncol(), nrow(), &out[off], scols, srows, GDT_Float64, 0, 0);
		if (err != CE_None ) { break; }
		double naflag = poBand->GetNoDataValue(&hasNA);
		if (!hasNA) { naflag = NAN; }
		setNAso(out, off, ncell, naflag, source[src].scale[i], source[src].offset[i], source[src].has_scale_offset[i]);
	}
*/

	GDALClose((GDALDatasetH) poDataset);
	if (err != CE_None ) {
		setError("cannot read values");
		return errout;
	}

	if (source[src].flipped) {
		vflip(out, ncell, srows, scols, nl);
	}

	return out;
}


void SpatRaster::readRowColGDAL(size_t src, std::vector<std::vector<double>> &out, size_t outstart, std::vector<int64_t> &rows, const std::vector<int64_t> &cols) {


	if (source[src].is_multidim) {
		readRowColMulti(src, out, outstart, rows, cols);
		return;
	}


	if (source[src].rotated) {
		setError("cannot read from rotated files. First use 'rectify'");
		return;
	}
	size_t n = rows.size();
	if (n < 1) {
		addWarning("nothing to extract");
		return;
	}

    GDALDataset *poDataset = openGDAL(source[src].filename, GDAL_OF_RASTER | GDAL_OF_READONLY, source[src].open_drivers, source[src].open_ops);

    if( poDataset == NULL )  {
		if (!file_exists(source[src].filename )) {
			setError("file does not exist: " + source[src].filename);
		} else {
			setError("cannot read from " + source[src].filename  );
		}
		return;
	}

	std::vector<size_t> lyrs = source[src].layers;
	size_t nl = lyrs.size();
	size_t outend = outstart + nl;

	size_t fnr = nrow() - 1;
	if (source[src].flipped) {
		for (size_t i=0; i<n; i++) {
			rows[i] = fnr - rows[i];
		}
	}

	std::vector<int> panBandMap;
	if (!source[src].in_order(true)) {
		panBandMap.reserve(nl);
		for (size_t i=0; i < nl; i++) {
			panBandMap.push_back(lyrs[i]+1);
		}
	}

	for (size_t i=outstart; i<outend; i++) {
		out[i] = std::vector<double> (n, NAN);
	}


	int64_t nr1 = nrow()-1;
	int64_t nc1 = ncol()-1;
	if (source[src].hasWindow) {
		nr1 = source[src].window.full_nrow - 1;
		nc1 = source[src].window.full_ncol - 1;
	}

	CPLErr err = CE_None;
	std::vector<double> value(nl);

	if (panBandMap.empty()) {
		for (size_t i=0; i < n; i++) {
			if ((cols[i] < 0) || (cols[i] > nc1) || (rows[i] < 0) || (rows[i] > nr1) ) continue;
			err = poDataset->RasterIO(GF_Read, cols[i], rows[i], 1, 1, &value[0], 1, 1, GDT_Float64, nl, NULL, 0, 0, 0, NULL);
			if (err == CE_None) {
				for (size_t j=0; j<nl; j++) {
					out[outstart+j][i] = value[j];
				}
			} else {
				break;
			}
		}
	} else {
		for (size_t i=0; i < n; i++) {
			if ((cols[i] < 0) || (rows[i] < 0)) continue;
			err = poDataset->RasterIO(GF_Read, cols[i], rows[i], 1, 1, &value[0], 1, 1, GDT_Float64, nl, &panBandMap[0], 0, 0, 0, NULL);
			if (err == CE_None) {
				for (size_t j=0; j<nl; j++) {
					out[outstart+j][i] = value[j];
				}
			} else {
				break;
			}
		}
	}
		//std::vector<double> naflags(nl, NAN);

	if (err == CE_None) {
		int hasNA;
		GDALRasterBand *poBand;
		for (size_t i=0; i<nl; i++) {
			poBand = poDataset->GetRasterBand(lyrs[i]+1);
			double naflag = poBand->GetNoDataValue(&hasNA);
			if (!hasNA) naflag = NAN;
			NAso(out[outstart+i], n, {naflag}, source[src].scale, source[src].offset, source[src].has_scale_offset, source[src].hasNAflag, source[src].NAflag);
		}
	}
	GDALClose((GDALDatasetH) poDataset);
	if (err != CE_None ) {
		setError("cannot read values");
		return;
	}
	

}


/*
std::vector<std::vector<double>> SpatRaster::readRowColGDAL(size_t src, std::vector<int64_t> &rows, const std::vector<int64_t> &cols) {

	std::vector<std::vector<double>> errout;
	if (source[src].rotated) {
		setError("cannot read from rotated files. First use 'rectify'");
		return errout;
	}

    GDALDataset *poDataset = openGDAL(source[src].filename, GDAL_OF_RASTER | GDAL_OF_READONLY, source[src].open_drivers, source[src].open_ops);

    if( poDataset == NULL )  {
		if (!file_exists(source[src].filename )) {
			setError("file does not exist: " + source[src].filename);
		} else {
			setError("cannot read from " + source[src].filename  );
		}
		return errout;
	}

	GDALRasterBand *poBand;

	std::vector<size_t> lyrs = source[src].layers;
	size_t nl = lyrs.size();
	size_t n = rows.size();

	size_t fnr = nrow() - 1;
	if (source[src].flipped) {
		for (size_t i=0; i<n; i++) {
			rows[i] = fnr - rows[i];
		}
	}

	std::vector<int> panBandMap;
	if (!source[src].in_order(true)) {
		panBandMap.reserve(nl);
		for (size_t i=0; i < nl; i++) {
			panBandMap.push_back(lyrs[i]+1);
		}
	}

	std::vector<double> out(n * nl, NAN);
	CPLErr err = CE_None;
	for (size_t i=0; i < n; i++) {
		if ((cols[i] < 0) || (rows[i] < 0)) continue;
		if (panBandMap.empty()) {
			err = poDataset->RasterIO(GF_Read, cols[i], rows[i], 1, 1, &out[i*nl], 1, 1, GDT_Float64, nl, NULL, 0, 0, 0, NULL);
		} else {
			err = poDataset->RasterIO(GF_Read, cols[i], rows[i], 1, 1, &out[i*nl], 1, 1, GDT_Float64, nl, &panBandMap[0], 0, 0, 0, NULL);
		}
		if (err != CE_None ) {
			break;x`
		}
	}

	if (err == CE_None ) {
		std::vector<double> naflags(nl, NAN);
		int hasNA;
		for (size_t i=0; i<nl; i++) {
			poBand = poDataset->GetRasterBand(lyrs[i]+1);
			double naflag = poBand->GetNoDataValue(&hasNA);
			if (hasNA)  naflags[i] = naflag;
		}
		NAso(out, n, naflags, source[src].scale, source[src].offset, source[src].has_scale_offset, source[src].hasNAflag, source[src].NAflag);
	}

	GDALClose((GDALDatasetH) poDataset);
	if (err != CE_None ) {
		setError("cannot read values");
		return errout;
	}


	size_t nr = rows.size();
	std::vector<std::vector<double>> r(nl, std::vector<double> (nr));
	for (size_t i=0; i<nr; i++) {
		for (size_t j=0; j<nl; j++) {
			size_t k = (i*nl) + j;
			r[j][i] = out[k];
		}
	}
	return r;
}
*/

/*
std::vector<double> SpatRaster::readRowColGDALFlat(size_t src, std::vector<int64_t> &rows, const std::vector<int64_t> &cols) {

	std::vector<double> errout;
	if (source[src].rotated) {
		setError("cannot read from rotated files. First use 'rectify'");
		return errout;
	}

    GDALDataset *poDataset = openGDAL(source[src].filename, GDAL_OF_RASTER | GDAL_OF_READONLY, source[src].open_drivers, source[src].open_ops);

    if( poDataset == NULL )  {
		if (!file_exists(source[src].filename )) {
			setError("file does not exist: " + source[src].filename);
		} else {
			setError("cannot read from " + source[src].filename  );
		}
		return errout;
	}

	GDALRasterBand *poBand;

	std::vector<size_t> lyrs = source[src].layers;
	size_t nl = lyrs.size();
	size_t n = rows.size();

	size_t fnr = nrow() - 1;
	if (source[src].flipped) {
		for (size_t i=0; i<n; i++) {
			rows[i] = fnr - rows[i];
		}
	}

	std::vector<int> panBandMap;
	if (!source[src].in_order(true)) {
		panBandMap.reserve(nl);
		for (size_t i=0; i < nl; i++) {
			panBandMap.push_back(lyrs[i]+1);
		}
	}

	std::vector<double> out(n * nl, NAN);
	CPLErr err = CE_None;
	for (size_t j=0; j < n; j++) {
		if ((cols[j] < 0) || (rows[j] < 0)) continue;
		if (panBandMap.empty()) {
			err = poDataset->RasterIO(GF_Read, cols[j], rows[j], 1, 1, &out[j*nl], 1, 1, GDT_Float64, nl, NULL, 0, 0, 0, NULL);
		} else {
			err = poDataset->RasterIO(GF_Read, cols[j], rows[j], 1, 1, &out[j*nl], 1, 1, GDT_Float64, nl, &panBandMap[0], 0, 0, 0, NULL);
		}
		if (err != CE_None ) {
			break;
		}
	}

	if (err == CE_None ) {
		std::vector<double> naflags(nl, NAN);
		int hasNA;
		for (size_t i=0; i<nl; i++) {
			poBand = poDataset->GetRasterBand(lyrs[i]+1);
			double naflag = poBand->GetNoDataValue(&hasNA);
			if (hasNA)  naflags[i] = naflag;
		}
		NAso(out, n, naflags, source[src].scale, source[src].offset, source[src].has_scale_offset, source[src].hasNAflag, source[src].NAflag);
	}

	GDALClose((GDALDatasetH) poDataset);
	if (err != CE_None ) {
		setError("cannot read values");
		return errout;
	}

	return out;
}
*/


// ncdf

bool ncdf_good_ends(std::string const &s) {
	std::vector<std::string> end = {"_bnds", "_bounds", "lat", "lon", "longitude", "latitude"};
	for (size_t i=0; i<end.size(); i++) {
		if (s.length() >= end[i].length()) {
			if (s.compare(s.length() - end[i].length(), s.length(), end[i]) == 0) {
				return false;
			}
		}
	}
	if (s == "x" || s == "y" || s == "northing" || s == "easting") {
		return false;
	}
	return true;
}

void ncdf_pick_most(std::vector<std::string> &sd, std::vector<std::string> &varname, std::vector<std::string> &longname, std::vector<int> &dim1, std::vector<int> &dim2) {
	if (sd.size() < 2) return;
	std::vector<int> ud = dim1;
	std::sort(ud.begin(), ud.end());
	ud.erase(std::unique(ud.begin(), ud.end()), ud.end());
	if (ud.size() > 1) {
		std::vector<std::string> tmpsd, tmpvarname, tmplongname;
		std::vector<int> tmpdim1, tmpdim2;
		int mx = ud[ud.size()-1];
		for (size_t i=0; i<sd.size(); i++) {
			if (dim1[i] == mx) {
				tmpsd.push_back(sd[i]);
				tmpvarname.push_back(varname[i]);
				tmplongname.push_back(longname[i]);
				tmpdim1.push_back(dim1[i]);
				tmpdim2.push_back(dim2[i]);
			}
		}
		sd = tmpsd;
		varname = tmpvarname;
		longname = tmplongname;
		dim1 = tmpdim1;
		dim2 = tmpdim2;
	}
}



bool SpatRaster::constructFromSDS(std::string filename, std::vector<std::string> meta, std::vector<int> subds, std::vector<std::string> subdsname, std::vector<std::string> options, std::string driver, bool noflip, bool guessCRS, std::vector<std::string> domains) {

	bool ncdf = driver =="netCDF";
	bool gtiff = driver == "GTiff";

	std::vector<std::vector<std::string>> info = parse_metadata_sds(meta);
	int n = info[0].size();

	if (gtiff && (subds[0] < 0) && subdsname[0].empty()) {
		subds.resize(n);
		std::iota(subds.begin(), subds.end(), 0);
	}
	std::vector<std::string> sd, varname, srcname;

// std::vector<unsigned> varnl;
// for selection based on nlyr

	if (info[0].empty()) {
		return false;
	}
	// select sds by index
	if ((!subds.empty() && (subds[0] >= 0))) {
		for (size_t i=0; i<subds.size(); i++) {
			if (subds[i] >=0 && subds[i] < n) {
				sd.push_back(info[0][subds[i]]);
				varname.push_back(info[1][i]);
			} else {
				std::string emsg = std::to_string(subds[i]+1) + " is not valid. There are " + std::to_string(info[0].size()) + " subdatasets\n";
				setError(emsg);
				return false;
			}
		}
	// select by name
	} else if (!subdsname.empty() && !subdsname[0].empty()) {
		for (size_t i=0; i<subdsname.size(); i++) {
			int w = where_in_vector(subdsname[i], info[1], false);
			if (w >= 0) {
				sd.push_back(info[0][w]);
				varname.push_back(info[1][w]);
			} else {
				std::string emsg = concatenate(info[1], ", ");
				emsg = subdsname[i] + " not found. Choose one of:\n" + emsg;
				setError(emsg);
				return false;
			}
		}

	// select all
	} else {
		// eliminate sources based on names like "*_bnds" and "lat"
		std::vector<int> rows, cols;
		for (size_t i=0; i<info[1].size(); i++) {			
			if (ncdf_good_ends(info[1][i])) {
				sd.push_back(info[0][i]);
				varname.push_back(info[1][i]);
				srcname.push_back(info[2][i]);
				try {
					rows.push_back(std::stol(info[3][i]));
				} catch(...) {
					rows.push_back(0);					
				}
				try {
					cols.push_back(std::stol(info[4][i]));
				} catch(...) {
					cols.push_back(0);					
				}
			}
		}
		if (sd.empty()) { // all were removed
			std::vector<size_t> nl(n);
			for (size_t i=0; i<nl.size(); i++) {
				nl[i] = stol(info[5][i]);
			}
			size_t mxnl = *max_element(nl.begin(), nl.end());
			for (size_t i=0; i<nl.size(); i++) {
				if (nl[i] == mxnl) {
					sd.push_back(info[0][i]);
					varname.push_back(info[1][i]);
					srcname.push_back(info[2][i]);
					rows.push_back(std::stol(info[3][i]));
					cols.push_back(std::stol(info[4][i]));
				}
			}
		}
		// pick the ones with most rows and then cols
		// to avoid picking the 1 or 2 "row" datasets
		ncdf_pick_most(sd, varname, srcname, rows, cols);
		ncdf_pick_most(sd, varname, srcname, cols, rows);
	}

	std::vector<size_t> srcnl;
	size_t cnt;

    for (cnt=0; cnt < sd.size(); cnt++) {
		if (constructFromFile(sd[cnt], {-1}, {""}, {}, options, noflip, guessCRS, domains)) break;
	}
//	source[0].source_name = srcname[cnt];

	std::vector<std::string> skipped, used;
	srcnl.push_back(nlyr());
	used.push_back(varname[0]);
	SpatRaster out;
	SpatOptions opt;
    for (size_t i=(cnt+1); i < sd.size(); i++) {
//		printf( "%s\n", sd[i].c_str() );
		bool success = out.constructFromFile(sd[i], {-1}, {""}, {}, options, noflip, guessCRS, domains);
		if (success) {
			if (out.compare_geom(*this, false, false, 0.1)) {
//				out.source	[0].source_name = srcname[i];
				addSource(out, false, opt);
				srcnl.push_back(out.nlyr());
				used.push_back(varname[i]);
			} else {
				skipped.push_back(varname[i]);
			}
		} else {
			skipped.push_back(varname[i]);
		}
	}

	if (!skipped.empty()) {
		std::string s="skipped sub-datasets (see 'describe(sds=TRUE)'):\n" + skipped[0];
		for (size_t i=1; i<skipped.size(); i++) {
			s += ", " + skipped[i];
			if ((i%3) == 0) s += "\n";
		}
		addWarning(s);
	}

	if (!(ncdf || gtiff)) {
		std::vector<std::string> lyrnames;
		for (size_t i=0; i<used.size(); i++) {
			std::vector<std::string> nms = { basename(used[i]) };
			recycle(nms, srcnl[i]);
			make_unique_names(nms);
			lyrnames.insert(lyrnames.end(), nms.begin(), nms.end());
		}
		if (!lyrnames.empty()) {
			setNames(lyrnames, false);
		}
	}

	size_t opsz = options.size();
	if (opsz > 0) {
		if (options[opsz-1] == "so=false") {
			options.resize(opsz-1);
		}
	}

	std::vector<std::string> metadata = get_metadata(filename, options);
	if (!metadata.empty()) {
		std::vector<std::string> tagnames, tagvalues;
//		get_tags(metadata, "NC_GLOBAL#TAG_", tagnames, tagvalues);
		get_tags(metadata, "NC_GLOBAL#", tagnames, tagvalues);
		for (size_t i=0; i<tagnames.size(); i++) {
			std::string mn = tagnames[i];
			if (!((mn == "_FillValue") || (mn == "grid_mapping") || (mn == "Conventions") || (mn == "created_by") || (mn == "created_date"))) {
				addTag(tagnames[i], tagvalues[i], driver);
			}
		}
	}

	return true;
}



std::vector<int64_t> ncdf_str2int64v(std::string s, std::string delim) {
	std::vector<int64_t> out;
	size_t pos = 0;
	while ((pos = s.find(delim)) != std::string::npos) {
		std::string v = s.substr(0, pos);
		s.erase(0, pos + 1);
		out.push_back(std::stoll(v));
	}
	out.push_back(std::stoll(s));
	return out;
}


bool get_long(std::string input, long &output) {
    try  {
		output = std::stol(input);
		return true;
    } catch (std::invalid_argument &e)  {
		return false;
    }
}


bool get_double(std::string input, double &output) {
    try {
		output = std::stod(input);
		return true;
    } catch (std::invalid_argument &e)  {
		return false;
    }
}


std::vector<int64_t> ncdf_time(const std::vector<std::string> &metadata, std::vector<std::string> vals, std::string &step, std::string &msg) {

	std::vector<int64_t> out, bad;
	if (vals.empty()) {
		step = "";
		return out;
	}

	std::vector<double> raw;
	raw.reserve(vals.size());
	for (size_t i=0; i<vals.size(); i++) {
		
		double dval;
		if (get_double(vals[i], dval)) {
			raw.push_back(dval);
		} else {
			step = "";
			return out;
		}
	}

	bool fu=false;
	bool fc=false;
	std::string origin;
	std::string calendar = "standard";
	for (size_t i=0; i<metadata.size(); i++) {

		if (!fc) {
			std::string pattern = "time#calendar=";
			std::size_t found = metadata[i].find(pattern);
			if (found != std::string::npos) {
				calendar = metadata[i];
				calendar.erase(calendar.begin(), calendar.begin()+pattern.size());
				fc = true;
			}
		}
		if (!fu) {
			std::string pattern = "time#units=";
			std::size_t found = metadata[i].find(pattern);
			if (found != std::string::npos) {
				origin = metadata[i];
				origin.erase(origin.begin(), origin.begin()+pattern.size());
				fu = true;
			}
		}
		if (fc & fu) break;
	}

	bool years = false;
	bool yearsbp = false;
	bool yearmonths = false;
	bool months = false;
	bool days = false;
	bool hours = false;
	bool minutes = false;
	bool seconds = false;
	bool foundorigin = false;

	if (fu) {
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
	}

	return out;
}




//NETCDF_DIM_k=0
//NETCDF_DIM_tile=0
//NETCDF_DIM_time=0
//NETCDF_VARNAME=NVEL



void ncdf_names(const std::vector<std::vector<std::string>> &m, std::vector<std::vector<std::string>> &out, std::vector<double> &depth, bool &has_depth, std::string &depth_name) {

	out.resize(3);
	if (m.empty()) return;

	std::string vname, lname, units = "";
	std::vector<std::string> b = m[0];
	for (size_t j=0; j<b.size(); j++) {
		size_t pos = b[j].find("NETCDF_VARNAME");
		if (pos != std::string::npos) {
			vname = b[j].erase(0, pos+15);
			continue;
		}
		pos = b[j].find("units=");
		if (pos != std::string::npos) {
			units = b[j].erase(0, pos+6);
			continue;
		}
		pos = b[j].find("long_name=");
		if (pos != std::string::npos) {
			lname = b[j].erase(0, pos+10);
			continue;
		}
		pos = b[j].find("standard_name=");
		if (pos != std::string::npos) {
			if (lname.empty()) {
				lname = b[j].erase(0, pos+14);
			}
		}
	}
	out[2] = {vname, lname, units};

// the below could be found analytically, but this is easy and safe
	for (size_t i=0; i<m.size(); i++) {
		std::string dim;
		b = m[i];
		for (size_t j=0; j<b.size(); j++) {
			size_t pos = b[j].find("NETCDF_DIM_");
			if (pos != std::string::npos) {
				size_t pos = b[j].find("NETCDF_DIM_time");
				if (pos != std::string::npos) {
					out[0].push_back( b[j].erase(0, pos+16) );
				} else {
					dim += b[j].erase(0, pos+11);
				}
			}
		}
		size_t pos = dim.find("=");
		double v;
		if (pos != std::string::npos) {
			if (i == 0) {
				depth_name = dim.substr(1, pos-1);
			}
			std::string dim2 = dim;
			dim2.erase(0, pos+1);
			try {
				v = std::stod(dim2);
			} catch(...) {
				v = NAN;
				has_depth = false;
			}	
			depth.push_back(v);
		} else {
			has_depth = false;
		}
		out[1].push_back(vname + dim);
	}

	return;
}


void SpatRasterSource::set_names_time_ncdf(std::vector<std::string> metadata, std::vector<std::vector<std::string>> bandmeta, std::string &msg) {

	if (bandmeta.empty()) return;

/*
	for (size_t i=0; i<bandmeta.size(); i++) {
		for (size_t j=0; j < bandmeta[i].size(); j++) {
			Rcpp::Rcout << bandmeta[i][j] << " ";
		}
		Rcpp::Rcout << std::endl;
	}
*/

	std::vector<std::vector<std::string>> nms;
	std::vector<double> mdepth;
	bool hasdepth = true;
	ncdf_names(bandmeta, nms, mdepth, hasdepth, depthname);
	

/*
	for (size_t i=0; i<nms.size(); i++) {
		for (size_t j=0; j<nms[i].size(); j++) {
			Rcpp::Rcout << nms[i][j] << " ";
		}
		Rcpp::Rcout << std::endl;
	}
*/
	
	if (hasdepth) {
		depth = mdepth;
		hasDepth = true;
		std::string pattern = depthname+"#units=";
		for (size_t i=0; i<metadata.size(); i++) {
			std::size_t found = metadata[i].find(pattern);
			if (found != std::string::npos) {
				depthunit = metadata[i];
				depthunit.erase(depthunit.begin(), depthunit.begin()+pattern.size());
				break;
			}
		}
	}
	
	if (!nms[1].empty()) {
		names = nms[1];
		make_unique_names(names);
	}
	source_name = nms[2][0];
	source_name_long = nms[2][1];

	if (!hasUnit) {
		if (nms[2][2].empty()) {
			unit = {""};
			hasUnit = false;
		} else {
			unit = {nms[2][2]};
			hasUnit = true;
		}
		recycle(unit, nlyr);
	}
	
	if (!nms[0].empty()) {
		std::string step;
		std::vector<int64_t> x;
		try {
			x = ncdf_time(metadata, nms[0], step, msg);
			if (x.size() == nlyr) {
				time = x;
				timestep = step;
				hasTime = true;
			}
		} catch(...) {
			msg = "could not extract time scale";
		}
	}
}




std::vector<std::vector<std::string>> grib_names(const std::vector<std::vector<std::string>> &m) {

	std::vector<std::vector<std::string>> out(4);
	if (m.empty()) return out;
	
	bool ft1 = false;
	bool ft2 = false;

	for (size_t i=0; i<m.size(); i++) {

		std::string comm, time1, time2, units = "";

		for (size_t j=0; j<m[i].size(); j++) {
			
			size_t pos = m[i][j].find("GRIB_COMMENT=");
			if (pos != std::string::npos) {
				comm = m[i][j];
				comm.erase(0, pos+13);
				lrtrim(comm);
				continue;
			}
			pos = m[i][j].find("GRIB_UNIT=");
			if (pos != std::string::npos) {
				units = m[i][j];
				units.erase(0, pos+10);
				str_replace(units, "[", "");
				str_replace(units, "]", "");
				lrtrim(units);
				continue;
			}
			pos = m[i][j].find("GRIB_VALID_TIME=");
			if (pos != std::string::npos) {
				std::string tmp = m[i][j];
				tmp.erase(0, pos+16);
				pos = tmp.find("sec");
				if (pos != std::string::npos) {
					tmp.erase(tmp.begin()+pos, tmp.end());
				}
				time1 = tmp;
				ft1 = true;
				continue;
			}
			pos = m[i][j].find("TIME=");
			if (pos != std::string::npos) {
				std::string tmp = m[i][j];
				tmp.erase(0, pos+5);
				pos = tmp.find("sec");
				if (pos != std::string::npos) {
					tmp.erase(tmp.begin()+pos, tmp.end());
				}
				time2 = tmp;
				ft2 = true;
			}
		}
		out[0].push_back(comm);
		out[1].push_back(units);
		out[2].push_back(time1);
		out[3].push_back(time2);
	}
	
	if (!ft1) {
		if (ft2) {
			out[2] = out[3];
		}
	}
	out.resize(3);
	return out;
}



void SpatRasterSource::set_names_time_grib(std::vector<std::vector<std::string>> bandmeta, std::string &msg) {

	if (bandmeta.empty()) return;
	
	std::vector<std::vector<std::string>> nms = grib_names(bandmeta);

	if (nms[0].size() != names.size()) return;

	for (size_t i=0; i<names.size(); i++) {
		names[i] += "; " + nms[0][i];
		str_replace(names[i], "0[-] ", "");
		str_replace_all(names[i], "\"", "");
	}

	if (nms[1].size() == nms[0].size()) {
		unit = {nms[1]};
	}

	bool hastime = false;
	std::vector<int64_t> tm;
	if (nms[2].size() == nms[0].size()) {
		hastime = true;
		int64_t tim;
		for (size_t i=0; i<nms[2].size(); i++) {
			if (nms[2][i].empty()) {
				hastime = false;
				break;
			}
			try {
				tim = stol(nms[2][i]);
			} catch(...) {
				hastime = false;
				break;
			}
			tm.push_back(tim);
		}
	}

	if (hastime) {
		time = tm;
		timestep = "seconds";
		hasTime = true;
	}
	
}



std::vector<std::vector<std::string>> tiff_names(const std::vector<std::vector<std::string>> &m) {

	std::vector<std::vector<std::string>> out(4);
	if (m.empty()) return out;
	
	for (size_t i=0; i<m.size(); i++) {

		std::string time, units = "";
		bool keepgoing = false;

		for (size_t j=0; j<m[i].size(); j++) {			
			size_t pos = m[i][j].find("UNIT=");
			if (pos == std::string::npos) {
				pos = m[i][j].find("unit=");
			}	
			if (pos != std::string::npos) {
				units = m[i][j];
				units.erase(0, pos+5);
				str_replace(units, "[", "");
				str_replace(units, "]", "");
				lrtrim(units);
				keepgoing = true;
				continue;
			}
			pos = m[i][j].find("TIME=");
			if (pos == std::string::npos) {
				pos = m[i][j].find("time=");
			}	
			if (pos != std::string::npos) {
				std::string tmp = m[i][j];
				tmp.erase(0, pos+5);
				pos = tmp.find("sec");
				if (pos != std::string::npos) {
					tmp.erase(tmp.begin()+pos, tmp.end());
				}
				time = tmp;
				keepgoing = true;
			}
			if (!keepgoing) return out;
		}
		out[1].push_back(units);
		out[2].push_back(time);
	}
	
	return out;
}


void SpatRasterSource::set_names_time_tif(std::vector<std::vector<std::string>> bandmeta, std::string &msg) {

	if (bandmeta.empty()) return;
	
	std::vector<std::vector<std::string>> nms = tiff_names(bandmeta);

	if (nms[1].size() == nlyr) {
		unit = {nms[0]};
	}

	bool hastime = false;
	std::vector<int64_t> tm;
	if (nms[1].size() == nlyr) {
		hastime = true;
		int64_t tim;
		for (size_t i=0; i<nms[1].size(); i++) {
			if (nms[1][i].empty()) {
				hastime = false;
				break;
			}
			try {
				tim = stol(nms[1][i]);
			} catch(...) {
				hastime = false;
				break;
			}
			tm.push_back(tim);
		}
	}

	if (hastime) {
		time = tm;
		timestep = "seconds";
		hasTime = true;
	}
	
}

