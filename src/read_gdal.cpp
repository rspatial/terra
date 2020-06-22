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


#include <algorithm>
#include <stdint.h>
#include <regex>

//#include "spatRaster.h"
#include "spatRasterMultiple.h"

#include "file_utils.h"
#include "string_utils.h"
#include "NA.h"
#include "date.h"


#include "gdal_priv.h"
#include "cpl_conv.h" // for CPLMalloc()
#include "cpl_string.h"
#include "ogr_spatialref.h"

#include "gdal_rat.h"
//#include "hdr.h"

#include "gdal_errors.h"

#if GDAL_VERSION_MAJOR >= 3
#include "proj.h"
#endif


void SpatRaster::gdalogrproj_init(std::string path) {
    GDALAllRegister();
    OGRRegisterAll(); 
	//GDALregistred = true;
#if GDAL_VERSION_MAJOR >= 3
	if (path != "") {
		const char *cp = path.c_str();
		proj_context_set_search_paths(PJ_DEFAULT_CTX, 1, &cp);
	}
#endif
}


SpatDataFrame GetRATdf(GDALRasterAttributeTable *pRAT) {

	SpatDataFrame out;
/*
	const char *GFU_type_string[] = {"GFT_Integer", "GFT_Real","GFT_String"};
	const char *GFU_usage_string[] = {"GFU_Generic", "GFU_PixelCount", "GFU_Name", "GFU_Min",
		"GFU_Max", "GFU_MinMax", "GFU_Red", "GFU_Green", "GFU_Blue", "GFU_Alpha", "GFU_RedMin",
		"GFU_GreenMin", "GFU_BlueMin", "GFU_AlphaMin", "GFU_RedMax", "GFU_GreenMax", "GFU_BlueMax",
		"GFU_AlphaMax", "GFU_MaxCount"};
	std::vector<std::string> GFT_type;
	std::vector<std::string> GFT_usage;
*/
	size_t nc = (int) pRAT->GetColumnCount();
	size_t nr = (int) pRAT->GetRowCount();

	for (size_t i=0; i<nc; i++) {
		GDALRATFieldType nc_type = pRAT->GetTypeOfCol(i);
//		GFT_type.push_back(GFU_type_string[nc_types[i]]);
//		GDALRATFieldUsage nc_usage = pRAT->GetUsageOfCol(i);
//		GFT_usage.push_back(GFU_usage_string[nc_usages[i]]);
		std::string name = pRAT->GetNameOfCol(i);
		if (nc_type == GFT_Integer) {
			std::vector<long> d(nr);
			for (size_t j=0; j<nr; j++) {
				d[j] = (int) pRAT->GetValueAsInt(j, i);
			}
			out.add_column(d, name);
		} else if (nc_type == GFT_Real) {
			std::vector<double> d(nr);
			for (size_t j=0; j<nr; j++) {
				d[j] = (double) pRAT->GetValueAsDouble(j, i);
			}
			out.add_column(d, name);
		} else if (nc_type == GFT_String) {
			std::vector<std::string> d(nr);
			for (size_t j=0; j<nr; j++) {
				d[j] = (std::string) pRAT->GetValueAsString(j, i);
			}
			out.add_column(d, name);
		}
	}
	return(out);
}



SpatCategories GetCategories(char **pCat) {
	size_t n = CSLCount(pCat);
	std::vector<std::string> nms(n);
	for (size_t i = 0; i<n; i++) {
		const char *field = CSLGetField(pCat, i);
		nms[i] = field;
	}
	SpatCategories scat;
	scat.labels = nms;
	scat.levels.resize(nms.size());
	std::iota(scat.levels.begin(), scat.levels.end(), 0);
	return(scat);
}


std::string basename_sds(std::string f) {
	const size_t i = f.find_last_of("\\/");
	if (std::string::npos != i) {
		f.erase(0, i + 1);
	}
	const size_t j = f.find_last_of(":");
	if (std::string::npos != j) {
		f.erase(0, j + 1);
	}
	f = std::regex_replace(f, std::regex(".hdf$"), "");
	f = std::regex_replace(f, std::regex(".nc$"), "");
	f = std::regex_replace(f, std::regex("\""), "");

	return f;
}

std::vector<std::vector<std::string>> metatime(std::vector<std::string> meta) {
	std::vector<std::vector<std::string>> out(meta.size());
	std::string delim = "=";
	for (size_t i=0; i<meta.size(); i++) {
		std::string s = meta[i];
		size_t pos = s.find(delim);
		if (pos != std::string::npos) {
			out[i].push_back(s.erase(pos+1, std::string::npos));
			out[i].push_back(s.erase(0, pos+1));
		} else {
			out[i].push_back(s);
		}
	}
	return out;
}


bool SpatRaster::constructFromSubDataSets(std::string filename, std::vector<std::string> meta, int subds, std::string subdsname) {

	std::vector<std::string> sd; //, nms;
	std::vector<std::string> dc; //, nms;
	std::string ndelim = "NAME=";
	std::string ddelim = "DESC=";
	for (size_t i=0; i<meta.size(); i++) {
		std::string s = meta[i];
		size_t pos = s.find(ndelim);
		if (pos != std::string::npos) {
			s.erase(0, pos + ndelim.length());
			sd.push_back(s);
		} else {
			size_t pos = s.find(ddelim);
			if (pos != std::string::npos) {
				s.erase(0, pos + ddelim.length());
				dc.push_back(s);
			}
		}
	}
	if (sd.size() == 0) {
		return false;
	}
	bool useDC = (dc.size() == sd.size());
	int sdsize = sd.size();
	if (subds >=0) {
		if (subds < sdsize) {
			sd = {sd[subds]};
			if (useDC) {
				dc = {dc[subds]};
			}
		} else {
			std::string emsg = std::to_string(subds) + " is not valid. There are " + std::to_string(sd.size()) + " subdatasets\n";
			setError(emsg);
			return false;
		}
	} else if (subdsname != "") {
		std::vector<std::string> shortnames = getlastpart(sd, ":");
		int w = where_in_vector(subdsname, shortnames);
		if (w >= 0) {
			sd = {sd[w]};
			if (useDC) {
				dc = {dc[w]};
			}			
		} else {
			std::string emsg = concatenate(shortnames, ", ");
			emsg = subdsname + " not found. Choose one of:\n" + emsg;
			setError(emsg);
			return false;
		}
	}
	
	bool success = constructFromFile(sd[0], -1, "");
	if (!success) {
		return false;
	}
	SpatRaster out;
    for (size_t i=1; i < sd.size(); i++) {
//		printf( "%s\n", sd[i].c_str() );
		success = out.constructFromFile(sd[i], -1, "");
		if (success) {
//			out.source[0].subdataset = true;
			addSource(out);
			if (out.msg.has_error) {
				//setError(out.msg.error);
				//return false;
				addWarning("skipped (different geometry): " + sd[i]);
			}
		} else {
			if (out.msg.has_error) {
				setError(out.msg.error);
			}
			return false;
		}
	}

	for (std::string& s : sd) s = basename_sds(s);
	success = setNames(sd);

	return true;
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
		OGRErr err = oSRS.exportToPrettyWkt(&cp);
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


SpatRasterStack::SpatRasterStack(std::string fname, std::vector<int> ids, bool useids) {

    GDALDataset *poDataset;
    poDataset = (GDALDataset *) GDALOpen( fname.c_str(), GA_ReadOnly );
    if( poDataset == NULL )  {
		if (!file_exists(fname)) {
			setError("file does not exist");
		} else {
			setError("cannot read from " + fname );
		}
		return;
	}
	
	std::string delim = "NAME=";
	char **metadata = poDataset->GetMetadata("SUBDATASETS");

	if (metadata == NULL) {
		setError("file has no subdatasets");
		GDALClose( (GDALDatasetH) poDataset );	
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
	SpatRaster sub;
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
				if (sub.constructFromFile(s, -1, "")) {
					if (!push_back(sub, basename_sds(s))) {
						addWarning("skipped (different geometry): " + s);
					}
				} else {
					addWarning("skipped (fail): " + s);
				}
			} 
		} 
	}
	GDALClose( (GDALDatasetH) poDataset );	
}


/*
bool ncdf_time(std::string filename, int &startdate, std::string &calendar) {
    GDALDataset *poDataset;
	std::string ftime = "NETCDF:\"" + filename + "\":time_bnds" ;
    poDataset = (GDALDataset *) GDALOpen( ftime.c_str(), GA_ReadOnly );
    if( poDataset == NULL )  {
		return false;
	}
	GDALRasterBand *poBand;
	poBand = poDataset->GetRasterBand(1);
	const char *pszv = nullptr;
	if (( pszv = poBand->GetMetadataItem("units")) != nullptr ) {			
		std::string s = pszv;
		std::string delim = "days since ";
		size_t pos = s.find(delim);
		if (pos == std::string::npos) {
			return false;
		}
		s.erase(0, delim.length());
		if (s.size() > 9) {
			try {
				int y = std::stoi(s.substr(0,4));
				int m = std::stoi(s.substr(5,2));
				int d = std::stoi(s.substr(8,2));
				std::vector<int> ymd = {y, m, d};
				startdate = date_from_ymd(ymd);
			} catch (...) {
				return false;
			}
		}
		const char *calpt = nullptr;
		if (( calpt = poBand->GetMetadataItem("calendar")) != nullptr ) {			
			calendar = calpt;
		}
	}
	GDALClose( (GDALDatasetH) poDataset );

	return true;
}

bool fixTime(std::vector<double> &time, int &startdate, std::string &calendar) {
	int nday = 0;
	if (calendar =="gregorian" || calendar =="proleptic_gregorian" || calendar=="standard") {
		nday = 366;
	} else if ((calendar == "365 day") || (calendar == "365_day")) {
		nday = 365;
	} else if (calendar == "noleap") {
		nday = 360;
	} else {
		nday = 360;
	}
	if (nday < 366) {
		return false;
	//	startyear = as.numeric( format(startDate, "%Y") )
	//	startmonth = as.numeric( format(startDate, "%m") )
	//	startday = as.numeric( format(startDate, "%d") )
	//	year = trunc( as.numeric(time)/nday )
	//	doy = (time - (year * nday))
	//	origin = paste(year+startyear, "-", startmonth, "-", startday, sep='')
	//	time = as.Date(doy, origin=origin)		
	} else {
		for (double& d : time) d += startdate;
		return true;
	}
}


//#include <iostream>
//#include "Rcpp.h"
*/

bool SpatRaster::constructFromFile(std::string fname, int subds, std::string subdsname) {

    GDALDataset *poDataset;
    poDataset = (GDALDataset *) GDALOpen(fname.c_str(), GA_ReadOnly );

    if( poDataset == NULL )  {
		if (!file_exists(fname)) {
			setError("file does not exist");
		} else {
			setError("cannot read from " + fname );
		}
		return false;
	}

	unsigned nl = poDataset->GetRasterCount();

/*
	char **metadata = poDataset->GetMetadataDomainList();
    for (size_t i=0; metadata[i] != NULL; i++) {
		Rcpp::Rcout << metadata[i] << std::endl;
	}
*/	
	
	if (nl == 0) {
		std::vector<std::string> meta;
		char **metadata = poDataset->GetMetadata("SUBDATASETS");
		if (metadata != NULL) {
			for (size_t i=0; metadata[i] != NULL; i++) {
				meta.push_back(metadata[i]);
			}
			return constructFromSubDataSets(fname, meta, subds, subdsname);
		} else {
			setError("no data detected in " + fname);
			return false;
		}
	}

	
	RasterSource s;
	s.ncol = poDataset->GetRasterXSize();
	s.nrow = poDataset->GetRasterYSize();
	s.nlyr = nl;
	s.nlyrfile = nl;
	s.layers.resize(nl);
    std::iota(s.layers.begin(), s.layers.end(), 0);

	s.flipped = false;
	s.rotated = false;
	double adfGeoTransform[6];	
	if( poDataset->GetGeoTransform( adfGeoTransform ) == CE_None ) {
		double xmin = adfGeoTransform[0]; /* left x */
		double xmax = xmin + adfGeoTransform[1] * s.ncol; /* w-e pixel resolution */
		//xmax = roundn(xmax, 9);
		double ymax = adfGeoTransform[3]; // top y 
		double ymin = ymax + s.nrow * adfGeoTransform[5]; 
		//ymin = roundn(ymin, 9);

		if (adfGeoTransform[5] > 0) {
			s.flipped = true;
			std::swap(ymin, ymax);
		}

		SpatExtent e(xmin, xmax, ymin, ymax);
		s.extent = e;
		
		if (adfGeoTransform[2] != 0 || adfGeoTransform[4] != 0) {
			s.rotated = true;
			addWarning("the data in this file are rotated. Use 'rectify' to fix that");
		}
	}

	s.memory = false;
	s.filename = fname;
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
	if (crs == "") {
		if (s.extent.xmin >= -180 && s.extent.xmax <= 360 && s.extent.ymin >= -90 && s.extent.ymax <= 90) {
			crs = "+proj=longlat +datum=WGS84";
		}
	}
	std::string msg;
	if (!s.srs.set({crs}, msg)) {
		setError(msg);
		return false;
	}

	GDALRasterBand  *poBand;
	//int nBlockXSize, nBlockYSize;
	double adfMinMax[2];
	int bGotMin, bGotMax;

//	s.layers.resize(1);


	/*
	std::string gdrv = poDataset->GetDriver()->GetDescription();
	Rcpp::Rcout << "driver: " << gdrv << std::endl;
	int startdate=0;
	std::string calendar = "";
	std::string unit = "";
	if (gdrv == "netCDF") {
		poBand = poDataset->GetRasterBand(1);

		const char *punit = nullptr;
		if (( punit = poBand->GetMetadataItem("units")) != nullptr ) {			
			unit = punit;								
		}		

		const char *pmeta = nullptr;
		if ((pmeta = poBand->GetMetadataItem("NETCDF_DIM_time")) != nullptr ) {
			s.time.resize(s.nlyr, NAN);
			s.time[0] = CPLAtofM(pmeta);
			s.hasTime = true;

			ncdf_time(fname, startdate, calendar);
		}
	}
	*/

	for (size_t i = 0; i < s.nlyr; i++) {
		poBand = poDataset->GetRasterBand(i+1);

		/*
		if (s.hasTime) {
			const char *pszValue = nullptr;
			if( (pszValue = poBand->GetMetadataItem("NETCDF_DIM_time")) != nullptr ) {
				s.time[i] = CPLAtofM(pszValue);
			}
		}
		*/

		int success;
	//	double naflag = poBand->GetNoDataValue(&success);
	//	if (success) {
	//		s.NAflag = naflag;
	//	} else {
	//		s.NAflag = NAN;
	//	}
		double offset = poBand->GetOffset(&success);
		if (success) {
			s.offset.push_back(offset);
			s.has_scale_offset.push_back(true);
		} else {
			s.offset.push_back(0);
			s.has_scale_offset.push_back(false);
		}
		double scale = poBand->GetScale(&success);
		if (success) {
			s.scale.push_back(scale);
			s.has_scale_offset[i] = true;
		} else {
			s.scale.push_back(1);
		}


		std::string dtype = GDALGetDataTypeName(poBand->GetRasterDataType());

		adfMinMax[0] = poBand->GetMinimum( &bGotMin );
		adfMinMax[1] = poBand->GetMaximum( &bGotMax );
		if( (bGotMin && bGotMax) ) {
			s.hasRange.push_back(true);
			s.range_min.push_back( adfMinMax[0] );
			s.range_max.push_back( adfMinMax[1] );
		} else {
			s.hasRange.push_back(false);
			s.range_min.push_back( NAN );
			s.range_max.push_back( NAN );
		}

		//if( poBand->GetOverviewCount() > 0 ) printf( "Band has %d overviews.\n", poBand->GetOverviewCount() );

		//GDALGetColorInterpretationName( poBand->GetColorInterpretation()) );

		GDALColorTable *ct;
		ct = poBand->GetColorTable();
		if( ct != NULL )	{
			s.hasColors.push_back(true);
		} else {
			s.hasColors.push_back(false);
		}

		GDALRasterAttributeTable *rat = poBand->GetDefaultRAT();
		if( rat != NULL )	{
			s.hasAttributes.push_back(true);
			SpatDataFrame df = GetRATdf(rat);
			s.atts.resize(i+1);
			s.atts[i] = df;
		} else {
			s.hasAttributes.push_back(false);
		}

		char **cat = poBand->GetCategoryNames();
		if( cat != NULL )	{
			s.hasCategories.push_back(true);
			SpatCategories scat = GetCategories(cat);
			s.cats.resize(i+1);
			s.cats[i] = scat;
		} else {
			s.hasCategories.push_back(false);
		}

		std::string bandname = poBand->GetDescription();
		if (bandname != "") {
			s.names.push_back(bandname);
		} else {
			if (s.nlyr > 1) {
				s.names.push_back(basename_noext(fname).substr(0,3) + "_" + std::to_string(i+1) ) ;
			} else {
				s.names.push_back(basename_noext(fname)) ;
			}
		}
	}
	GDALClose( (GDALDatasetH) poDataset );

	s.hasValues = true;
	setSource(s);

	//if (unit != "") setUnit({unit});
	//if (s.hasTime) fixTime(s.time, startdate, calendar);

	return true;
}


bool SpatRaster::readStartGDAL(unsigned src) {
    GDALDataset *poDataset;
    //GDALAllRegister();
	const char* pszFilename = source[src].filename.c_str();
	poDataset = (GDALDataset *) GDALOpen( pszFilename, GA_ReadOnly );
    if( poDataset == NULL )  {
		setError("cannot read from " + source[src].filename );
		return false;
	}
    source[src].gdalconnection = poDataset;
	source[src].open_read = true;
	return(true);
}

bool SpatRaster::readStopGDAL(unsigned src) {
	if (source[src].gdalconnection != NULL) {
		GDALClose( (GDALDatasetH) source[src].gdalconnection);
	}
	source[src].open_read = false;
	return true;
}


//#include <iostream>
//#include "Rcpp.h"

void NAso(std::vector<double> &d, size_t n, const std::vector<double> &flags, const std::vector<double> &scale, const std::vector<double>  &offset, const std::vector<bool> &haveso){
	size_t nl = flags.size();

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
				double na = NAN;
				std::replace(d.begin()+start, d.begin()+start+n, flag, na); 
			}
		} 
		if (haveso[i]) {
			for (size_t j=start; j<(start+n); j++) {
				d[j] = d[j] * scale[i] + offset[i];
			}
		}
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


std::vector<double> SpatRaster::readChunkGDAL(unsigned src, unsigned row, unsigned nrows, unsigned col, unsigned ncols) {

	std::vector<double> errout;
	if (source[src].rotated) {
		setError("cannot read from rotated files. First use 'rectify'");
		return errout;
	}

	GDALRasterBand  *poBand;
	unsigned ncell = ncols * nrows;
	unsigned nl = source[src].nlyr;
	std::vector<double> out(ncell * nl);
	int hasNA;
	std::vector<double> naflags(nl, NAN);
	CPLErr err = CE_None;

	std::vector<int> panBandMap;
	if (!source[src].in_order()) {
		panBandMap.reserve(nl);
		for (size_t i=0; i < nl; i++) {
			panBandMap.push_back(source[src].layers[i]+1);
		}
	}

		
	err = source[src].gdalconnection->RasterIO(GF_Read, col, row, ncols, nrows, &out[0], ncols, nrows, GDT_Float64, nl, &panBandMap[0], 0, 0, 0, NULL);
	if (err == CE_None ) { 
		for (size_t i=0; i<nl; i++) {
			poBand = source[src].gdalconnection->GetRasterBand(source[src].layers[i]+1);
			double naflag = poBand->GetNoDataValue(&hasNA);
			if (hasNA)  naflags[i] = naflag;
		}
		NAso(out, ncell, naflags, source[src].scale, source[src].offset, source[src].has_scale_offset);
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
		return errout;
	}

	if (source[src].flipped) {
		vflip(out, ncell, nrows, ncols, nl);
	}

	return(out);
}





std::vector<double> SpatRaster::readValuesGDAL(unsigned src, unsigned row, unsigned nrows, unsigned col, unsigned ncols) {

	std::vector<double> errout;
	if (source[src].rotated) {
		setError("cannot read from rotated files. First use 'rectify'");
		return errout;
	}

    GDALDataset *poDataset;
	GDALRasterBand *poBand;
	const char* pszFilename = source[src].filename.c_str();
    poDataset = (GDALDataset *) GDALOpen(pszFilename, GA_ReadOnly);
    if( poDataset == NULL )  {
		setError("cannot read values. Does the file still exist?");
		return errout;
	}
	unsigned ncell = ncols * nrows;
	unsigned nl = source[src].nlyr;
	std::vector<double> out(ncell*nl);

/*	int* panBandMap = NULL;
	if (!source[src].in_order()) {
		panBandMap = (int *) CPLMalloc(sizeof(int) * nl);
		for (size_t i=0; i < nl; i++) {
			panBandMap[i] = source[src].layers[i] + 1;
		}
	}
*/
	std::vector<int> panBandMap;
	if (!source[src].in_order()) {
		panBandMap.reserve(nl);
		for (size_t i=0; i < nl; i++) {
			panBandMap.push_back(source[src].layers[i]+1);
		}
	}
		
	int hasNA;
	std::vector<double> naflags(nl, NAN);
	CPLErr err = CE_None;
	err = poDataset->RasterIO(GF_Read, col, row, ncols, nrows, &out[0], ncols, nrows, GDT_Float64, nl, &panBandMap[0], 0, 0, 0, NULL);

	if (err == CE_None ) { 
		for (size_t i=0; i<nl; i++) {
			poBand = poDataset->GetRasterBand(source[src].layers[i]+1);
			double naf = poBand->GetNoDataValue(&hasNA);
			if (hasNA)  naflags[i] = naf;
		}
		NAso(out, ncell, naflags, source[src].scale, source[src].offset, source[src].has_scale_offset);
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



std::vector<double> SpatRaster::readGDALsample(unsigned src, unsigned srows, unsigned scols) {

	std::vector<double> errout;
	if (source[src].rotated) {
		setError("cannot read from rotated files. First use 'rectify'");
		return errout;
	}

    GDALDataset *poDataset;
	GDALRasterBand *poBand;
    //GDALAllRegister();
	const char* pszFilename = source[src].filename.c_str();
    poDataset = (GDALDataset *) GDALOpen(pszFilename, GA_ReadOnly);
    if( poDataset == NULL )  {
		setError("no data");
		return errout;
	}
	unsigned ncell = scols * srows;
	unsigned nl = source[src].nlyr;
	std::vector<double> out(ncell*nl);
	int hasNA;
	CPLErr err = CE_None;

	std::vector<double> naflags(nl, NAN);

	std::vector<int> panBandMap;
	if (!source[src].in_order()) {
		panBandMap.reserve(nl);
		for (size_t i=0; i < nl; i++) {
			panBandMap.push_back(source[src].layers[i]+1);
		}
	}
		
	err = poDataset->RasterIO(GF_Read, 0, 0, ncol(), nrow(), &out[0], scols, srows, GDT_Float64, nl, &panBandMap[0], 0, 0, 0, NULL);
	if (err == CE_None ) { 
		for (size_t i=0; i<nl; i++) {
			poBand = poDataset->GetRasterBand(source[src].layers[i]+1);
			double naflag = poBand->GetNoDataValue(&hasNA);
			if (hasNA)  naflags[i] = naflag;
		}
		NAso(out, ncell, naflags, source[src].scale, source[src].offset, source[src].has_scale_offset);
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



std::vector<std::vector<double>> SpatRaster::readRowColGDAL(unsigned src, std::vector<unsigned> &rows, const std::vector<unsigned> &cols) {

	std::vector<std::vector<double>> errout;
	if (source[src].rotated) {
		setError("cannot read from rotated files. First use 'rectify'");
		return errout;
	}

    GDALDataset *poDataset;
	GDALRasterBand *poBand;
    //GDALAllRegister();
	const char* pszFilename = source[src].filename.c_str();
    poDataset = (GDALDataset *) GDALOpen(pszFilename, GA_ReadOnly);
    if( poDataset == NULL )  {
		return errout;
	}

	std::vector<unsigned> lyrs = source[src].layers;
	unsigned nl = lyrs.size();
	unsigned n = rows.size();

	size_t fnr = nrow() - 1;
	if (source[src].flipped) {
		for (size_t i=0; i<n; i++) {
			rows[i] = fnr - rows[i];
		}
	}

	std::vector<int> panBandMap;
	if (!source[src].in_order()) {
		panBandMap.reserve(nl);
		for (size_t i=0; i < nl; i++) {
			panBandMap.push_back(lyrs[i]+1);
		}
	}

	std::vector<double> out(n * nl, NAN);
	CPLErr err = CE_None;
	for (size_t j=0; j < n; j++) {
		if (is_NA(cols[j]) | is_NA(rows[j])) continue;
		err = poDataset->RasterIO(GF_Read, cols[j], rows[j], 1, 1, &out[j*nl], 1, 1, GDT_Float64, nl, &panBandMap[0], 0, 0, 0, NULL);
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
		NAso(out, n, naflags, source[src].scale, source[src].offset, source[src].has_scale_offset);
	}

	GDALClose((GDALDatasetH) poDataset);
	if (err != CE_None ) {
		setError("cannot read values");
		return errout;
	}

	
	std::vector<std::vector<double>> r(nl);
	for (size_t i=0; i<nl; i++) {
		size_t off = n * i;
		r[i] = std::vector<double>(out.begin()+off, out.begin()+off+n);  
	}
	return r;
}




/*
void setNAso(std::vector<double> &d, size_t start, size_t n, double flag, double scale, double offset, bool haveso){
	n += start;
	if (!std::isnan(flag)) {
		std::replace(d.begin()+start, d.begin()+n, flag, NAN); 
	}
	if (haveso) {
		for (size_t j=start; j<n; j++) {
			d[j] = d[j] * scale + offset;
		}
	}
}

*/


/*

template <typename T>
void set_NA_so(const std::vector<T> &lyr, double naflag, std::vector<double> &out, const size_t &cell, double scale, double offset, bool haveso){
	size_t n = lyr.size();
	if (haveso) {
		if (!std::isnan(naflag)) {
			T flag = naflag;
			for (size_t j=0; j<n; j++) {
				if (lyr[j] == flag) {
					out[cell+j] = NAN;
				} else {
					out[cell+j] = lyr[j] * scale + offset;
				}
			}
		} else {
			for (size_t j=0; j<n; j++) {
				out[cell+j] = lyr[j] * scale + offset;
			}
		}
	} else {
		if (!std::isnan(naflag)) {
			T flag = naflag;
			for (size_t j=0; j<n; j++) {
				if (lyr[j] == flag) {
					out[cell+j] = NAN;
				} else {
					out[cell+j] = lyr[j];
				}
			}
		} else {
			for (size_t j=0; j<n; j++) {
				out[cell+j] = lyr[j];
			}
		}
	}
}


template <typename T>
void set_NA_so2(const std::vector<T> &lyr, double naflag, std::vector<double> &out, double scale, double offset, bool haveso){
	size_t n = lyr.size();
	if (haveso) {
		if (!std::isnan(naflag)) {
			T flag = naflag;
			for (size_t j=0; j<n; j++) {
				if (lyr[j] == flag) {
					out[j] = NAN;
				} else {
					out[j] = lyr[j] * scale + offset;
				}
			}
		} else {
			for (size_t j=0; j<n; j++) {
				out[j] = lyr[j] * scale + offset;
			}
		}
	} else {
		if (!std::isnan(naflag)) {
			T flag = naflag;
			for (size_t j=0; j<n; j++) {
				if (lyr[j] == flag) {
					out[j] = NAN;
				} else {
					out[j] = lyr[j];
				}
			}
		} else {
			for (size_t j=0; j<n; j++) {
				out[j] = lyr[j];
			}
		}
	}
}

*/




/*
bool SpatRaster::constructFromFiles(std::vector<std::string> fnames) {

	SpatRaster r = SpatRaster(fnames[0], -1);
	setSource(r.source[0]);
	for (size_t i=1; i<fnames.size(); i++) {
		r = SpatRaster(fnames[i], -1);
		if (!compare_geom(r, false, true, true)) {
			setError("geometry of " + fnames[i] + " does not match previous sources");
			return false;
		} else {
			addSource(r);
		}
	}
	return true;
}
*/



/*
void applyScaleOffset(std::vector<double> &d, double scale, double offset, bool haveso) {
	if (haveso) {
		for (size_t i=0; i<d.size(); i++) {
			d[i] = d[i] * scale + offset;
		}
	}
}
*/

/*
std::vector<double> SpatRaster::readValuesGDAL(unsigned src, unsigned row, unsigned nrows, unsigned col, unsigned ncols) {

	std::vector<double> errout;
    GDALDataset *poDataset;
	GDALRasterBand *poBand;
    //GDALAllRegister();
	const char* pszFilename = source[src].filename.c_str();
    poDataset = (GDALDataset *) GDALOpen(pszFilename, GA_ReadOnly);
    if( poDataset == NULL )  {
		return errout;
	}
	unsigned ncell = ncols * nrows;
	unsigned nl = source[src].nlyr;
	std::vector<double> out(ncell*nl);
	int hasNA;
	CPLErr err = CE_None;
	for (size_t i=0; i < nl; i++) {
		unsigned cell = ncell * i;
		unsigned thislayer = source[src].layers[i];
		poBand = poDataset->GetRasterBand(thislayer + 1);
		GDALDataType gdtype = poBand->GetRasterDataType();
		double naflag = poBand->GetNoDataValue(&hasNA);
		if (!hasNA) { naflag = NAN; }
		if (gdtype == GDT_Float64) {
			std::vector<double> lyrout(ncell);
			err = poBand->RasterIO(GF_Read, col, row, ncols, nrows, &lyrout[0], ncols, nrows, gdtype, 0, 0);
			if (err != CE_None ) { break; }
			set_NA_so(lyrout, naflag, out, cell, source[src].scale[i], source[src].offset[i], source[src].has_scale_offset[i]);
		} else if (gdtype == GDT_Float32) {
			std::vector<float> lyrout(ncell);
			err = poBand->RasterIO(GF_Read, col, row, ncols, nrows, &lyrout[0], ncols, nrows, gdtype, 0, 0);
			if (err != CE_None ) { break; }
			set_NA_so(lyrout, naflag, out, cell, source[src].scale[i], source[src].offset[i], source[src].has_scale_offset[i]);
		} else if (gdtype == GDT_Byte) {
			std::vector<uint8_t> lyrout(ncell);
			err = poBand->RasterIO(GF_Read, col, row, ncols, nrows, &lyrout[0], ncols, nrows, gdtype, 0, 0);
			if (err != CE_None ) { break; }
			set_NA_so(lyrout, naflag, out, cell, source[src].scale[i], source[src].offset[i], source[src].has_scale_offset[i]);
		} else if (gdtype == GDT_Int16) {
			std::vector<int16_t> lyrout(ncell);
			err = poBand->RasterIO(GF_Read, col, row, ncols, nrows, &lyrout[0], ncols, nrows, gdtype, 0, 0);
			if (err != CE_None ) { break; }
			set_NA_so(lyrout, naflag, out, cell, source[src].scale[i], source[src].offset[i], source[src].has_scale_offset[i]);
		} else if (gdtype == GDT_UInt16) {
			std::vector<uint16_t> lyrout(ncell);
			err = poBand->RasterIO(GF_Read, col, row, ncols, nrows, &lyrout[0], ncols, nrows, gdtype, 0, 0);
			if (err != CE_None ) { break; }
			set_NA_so(lyrout, naflag, out, cell, source[src].scale[i], source[src].offset[i], source[src].has_scale_offset[i]);
		} else if (gdtype == GDT_Int32) {
			std::vector<int32_t> lyrout(ncell);
			err = poBand->RasterIO(GF_Read, col, row, ncols, nrows, &lyrout[0], ncols, nrows, gdtype, 0, 0);
			if (err != CE_None ) { break ;}
			set_NA_so(lyrout, naflag, out, cell, source[src].scale[i], source[src].offset[i], source[src].has_scale_offset[i]);
		} else if (gdtype == GDT_UInt32) {
			std::vector<uint32_t> lyrout(ncell);
			err = poBand->RasterIO(GF_Read, col, row, ncols, nrows, &lyrout[0], ncols, nrows, gdtype, 0, 0);
			if (err != CE_None ) { break ;}
			set_NA_so(lyrout, naflag, out, cell, source[src].scale[i], source[src].offset[i], source[src].has_scale_offset[i]);
		} else {
			//int tbd
		}
	}

	GDALClose((GDALDatasetH) poDataset);
	if (err != CE_None ) {
		setError("cannot read values");
		return errout;
	}
	return out;
}
*/




/*
std::vector<std::vector<double>> SpatRaster::readRowColGDAL(unsigned src, const std::vector<unsigned> &rows, const std::vector<unsigned> &cols) {
    GDALDataset *poDataset;
	GDALRasterBand *poBand;
    //GDALAllRegister();
	const char* pszFilename = source[src].filename.c_str();
	std::vector<std::vector<double>> errout;
    poDataset = (GDALDataset *) GDALOpen(pszFilename, GA_ReadOnly);
    if( poDataset == NULL )  {
		return errout;
	}

	std::vector<unsigned> lyrs = source[src].layers;
	unsigned nl = lyrs.size();
	unsigned n = rows.size();
	std::vector<std::vector<double>> out(nl, std::vector<double>(n, NAN));

	CPLErr err = CE_None;
	int hasNA;
//	unsigned offset;

	for (size_t i=0; i<nl; i++) {
		//offset = n * i;
		poBand = poDataset->GetRasterBand(lyrs[i] + 1);
		GDALDataType gdtype = poBand->GetRasterDataType();
		double naflag = poBand->GetNoDataValue(&hasNA);
		if (!hasNA) { naflag = NAN; }

		if (gdtype == GDT_Float64) {
			std::vector<double> lyrout(n, NAN);
			for (size_t j=0; j < n; j++) {
				if (is_NA(cols[j]) | is_NA(rows[j])) continue;
				err = poBand->RasterIO(GF_Read, cols[j], rows[j], 1, 1, &lyrout[j], 1, 1, gdtype, 0, 0);
				if (err != CE_None ) { break ;}
			}
			set_NA_so2(lyrout, naflag, out[i], source[src].scale[i], source[src].offset[i], source[src].has_scale_offset[i]);
		} else if (gdtype == GDT_Float32) {
			std::vector<float> lyrout(n, NAN);
			for (size_t j=0; j < n; j++) {
				if (is_NA(cols[j]) | is_NA(rows[j])) continue;
				err = poBand->RasterIO(GF_Read, cols[j], rows[j], 1, 1, &lyrout[j], 1, 1, gdtype, 0, 0);
				if (err != CE_None ) { break ;}
			}
			set_NA_so2(lyrout, naflag, out[i], source[src].scale[i], source[src].offset[i], source[src].has_scale_offset[i]);
		} else if (gdtype == GDT_Byte) {
			int8_t navalue = NA<int8_t>::value;
			std::vector<int8_t> lyrout(n, navalue);
			for (size_t j=0; j < n; j++) {
				if (is_NA(cols[j]) | is_NA(rows[j])) continue;
				err = poBand->RasterIO(GF_Read, cols[j], rows[j], 1, 1, &lyrout[j], 1, 1, gdtype, 0, 0);
				if (err != CE_None ) { break ;}
			}
			set_NA_so2(lyrout, naflag, out[i], source[src].scale[i], source[src].offset[i], source[src].has_scale_offset[i]);
		} else if (gdtype == GDT_Int16) {
			int16_t navalue = NA<int16_t>::value;
			std::vector<int16_t> lyrout(n, navalue);
			for (size_t j=0; j < n; j++) {
				if (is_NA(cols[j]) | is_NA(rows[j])) continue;
				err = poBand->RasterIO(GF_Read, cols[j], rows[j], 1, 1, &lyrout[j], 1, 1, gdtype, 0, 0);
				if (err != CE_None ) { break ;}
			}
			set_NA_so2(lyrout, naflag, out[i], source[src].scale[i], source[src].offset[i], source[src].has_scale_offset[i]);
		} else if (gdtype == GDT_UInt16) {
			uint16_t navalue = NA<uint16_t>::value;
			std::vector<uint16_t> lyrout(n, navalue);
			for (size_t j=0; j < n; j++) {
				if (is_NA(cols[j]) | is_NA(rows[j])) continue;
				err = poBand->RasterIO(GF_Read, cols[j], rows[j], 1, 1, &lyrout[j], 1, 1, gdtype, 0, 0);
				if (err != CE_None ) { break ;}
			}
			set_NA_so2(lyrout, naflag, out[i], source[src].scale[i], source[src].offset[i], source[src].has_scale_offset[i]);
		} else if (gdtype == GDT_Int32) {
			int32_t navalue = NA<int32_t>::value;
			std::vector<int32_t> lyrout(n, navalue);
			for (size_t j=0; j < n; j++) {
				if (is_NA(cols[j]) | is_NA(rows[j])) continue;
				err = poBand->RasterIO(GF_Read, cols[j], rows[j], 1, 1, &lyrout[j], 1, 1, gdtype, 0, 0);
				if (err != CE_None ) { break ;}
			}
			set_NA_so2(lyrout, naflag, out[i], source[src].scale[i], source[src].offset[i], source[src].has_scale_offset[i]);
		} else if (gdtype == GDT_UInt32) {
			uint32_t navalue = NA<uint32_t>::value;
			std::vector<uint32_t> lyrout(n, navalue);
			for (size_t j=0; j < n; j++) {
				if (is_NA(cols[j]) | is_NA(rows[j])) continue;
				err = poBand->RasterIO(GF_Read, cols[j], rows[j], 1, 1, &lyrout[j], 1, 1, gdtype, 0, 0);
				if (err != CE_None ) { break ;}
			}
			set_NA_so2(lyrout, naflag, out[i], source[src].scale[i], source[src].offset[i], source[src].has_scale_offset[i]);
		} else {
			setError("unknown data type");
		}
	}
	GDALClose((GDALDatasetH) poDataset);
	if (err != CE_None ) {
		setError("cannot read values");
		return errout;
	}
	return out;
}
*/


/*
std::vector<double> SpatRaster::readGDALsample(unsigned src, unsigned srows, unsigned scols) {

    GDALDataset *poDataset;
	GDALRasterBand *poBand;
    //GDALAllRegister();
	const char* pszFilename = source[src].filename.c_str();
    poDataset = (GDALDataset *) GDALOpen(pszFilename, GA_ReadOnly);
	std::vector<double> errout;
    if( poDataset == NULL )  {
		return errout;
	}
	unsigned ncell = scols * srows;
	unsigned cell;
	unsigned nl = source[src].nlyr;
	std::vector<double> out(ncell*nl);
	int hasNA;
	CPLErr err = CE_None;
	for (size_t i=0; i < nl; i++) {
		cell = ncell * i;
		poBand = poDataset->GetRasterBand(source[src].layers[i] + 1);
		GDALDataType gdtype = poBand->GetRasterDataType();
		double naflag = poBand->GetNoDataValue(&hasNA);
		if (!hasNA) { naflag = NAN; }
		if (gdtype == GDT_Float64) {
			std::vector<double> lyrout(ncell);
			err = poBand->RasterIO(GF_Read, 0, 0, ncol(), nrow(), &lyrout[0], scols, srows, gdtype, 0, 0);
			if (err != CE_None ) { break; }
			set_NA_so(lyrout, naflag, out, cell, source[src].scale[i], source[src].offset[i], source[src].has_scale_offset[i]);
		} else if (gdtype == GDT_Float32) {
			std::vector<float> lyrout(ncell);
			err = poBand->RasterIO(GF_Read, 0, 0, ncol(), nrow(), &lyrout[0], scols, srows, gdtype, 0, 0);
			if (err != CE_None ) { break; }
			set_NA_so(lyrout, naflag, out, cell, source[src].scale[i], source[src].offset[i], source[src].has_scale_offset[i]);
		} else if (gdtype == GDT_Byte) {
			std::vector<int8_t> lyrout(ncell);
			err = poBand->RasterIO(GF_Read, 0, 0, ncol(), nrow(), &lyrout[0], scols, srows, gdtype, 0, 0);
			if (err != CE_None ) { break; }
			set_NA_so(lyrout, naflag, out, cell, source[src].scale[i], source[src].offset[i], source[src].has_scale_offset[i]);
		} else if (gdtype == GDT_Int16) {
			std::vector<int16_t> lyrout(ncell);
			err = poBand->RasterIO(GF_Read, 0, 0, ncol(), nrow(), &lyrout[0], scols, srows, gdtype, 0, 0);
			if (err != CE_None ) { break; }
			set_NA_so(lyrout, naflag, out, cell, source[src].scale[i], source[src].offset[i], source[src].has_scale_offset[i]);
		} else if (gdtype == GDT_UInt16) {
			std::vector<uint16_t> lyrout(ncell);
			err = poBand->RasterIO(GF_Read, 0, 0, ncol(), nrow(), &lyrout[0], scols, srows, gdtype, 0, 0);
			if (err != CE_None ) { break; }
			set_NA_so(lyrout, naflag, out, cell, source[src].scale[i], source[src].offset[i], source[src].has_scale_offset[i]);
		} else if (gdtype == GDT_Int32) {
			std::vector<int32_t> lyrout(ncell);
			err = poBand->RasterIO(GF_Read, 0, 0, ncol(), nrow(), &lyrout[0], scols, srows, gdtype, 0, 0);
			if (err != CE_None ) { break ;}
			set_NA_so(lyrout, naflag, out, cell, source[src].scale[i], source[src].offset[i], source[src].has_scale_offset[i]);
		} else if (gdtype == GDT_UInt32) {
			std::vector<uint32_t> lyrout(ncell);
			err = poBand->RasterIO(GF_Read, 0, 0, ncol(), nrow(), &lyrout[0], scols, srows, gdtype, 0, 0);
			if (err != CE_None ) { break ;}
			set_NA_so(lyrout, naflag, out, cell, source[src].scale[i], source[src].offset[i], source[src].has_scale_offset[i]);
		} else {
			//int tbd
		}
	}
	GDALClose((GDALDatasetH) poDataset);
	if (err != CE_None ) {
		setError("cannot read values");
		return errout;
	}
	return out;
}
*/
