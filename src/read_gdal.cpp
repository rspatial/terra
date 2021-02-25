// Copyright (c) 2018-2021  Robert J. Hijmans
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

//#include "spatRaster.h"
#include "spatRasterMultiple.h"

#include "file_utils.h"
#include "string_utils.h"
#include <regex>
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



SpatDataFrame GetCOLdf(GDALColorTable *pCT) {

	SpatDataFrame out;
	size_t nc = (int) pCT->GetColorEntryCount();

	out.add_column(1, "red");
	out.add_column(1, "green");
	out.add_column(1, "blue");
	out.add_column(1, "alpha");
	out.reserve(nc);

	for (size_t i=0; i<nc; i++) {	
		const GDALColorEntry * col = pCT->GetColorEntry(i);
		out.iv[0].push_back(col->c1);
		out.iv[1].push_back(col->c2);
		out.iv[2].push_back(col->c3);
		out.iv[3].push_back(col->c4);
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
				if (sub.constructFromFile(s, {-1}, {""})) {
					if (!push_back(sub, basename_sds(s), "", "", true)) {
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



bool SpatRaster::constructFromFile(std::string fname, std::vector<int> subds, std::vector<std::string> subdsname) {

    GDALDataset *poDataset;
    poDataset = (GDALDataset *) GDALOpen(fname.c_str(), GA_ReadOnly );

    if( poDataset == NULL )  {
		if (!file_exists(fname)) {
			setError("file does not exist: " + fname);
		} else {
			setError("cannot read from " + fname );
		}
		return false;
	}

	int nl = poDataset->GetRasterCount();
	std::string gdrv = poDataset->GetDriver()->GetDescription();

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
			return constructFromSDS(fname, meta, subds, subdsname, gdrv=="netCDF"); 
		} else {
			setError("no data detected in " + fname);
			return false;
		}
	}


	SpatRasterSource s;
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
	} else {
		SpatExtent e(0, 1, 0, 1);
		s.extent = e;
		if (gdrv=="netCDF") {
			#ifndef standalone
			setMessage("ncdf extent");
			#else 
			addWarning("unknown extent. Cells not equally spaced?");
			#endif
		} else {
			addWarning("unknown extent");
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
			crs = "EPSG:4326";
			s.parameters_changed = true;
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

	//Rcpp::Rcout << "driver: " << gdrv << std::endl;

//	std::string unit = "";

	std::string varname = basename_noext(fname).substr(0,3);
	std::vector<std::vector<std::string>> bandmeta(s.nlyr);


	for (size_t i = 0; i < s.nlyr; i++) {
		poBand = poDataset->GetRasterBand(i+1);

		if (gdrv == "netCDF") {
			char **m = poBand->GetMetadata();
			while (*m != nullptr) {
				bandmeta[i].push_back(*m++);
			}
		}

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

		GDALColorTable *ct = poBand->GetColorTable();
		if( ct != NULL ) {
			s.hasColors.push_back(true);
			s.cols.resize(i+1);
			s.cols[i] = GetCOLdf(ct);
		} else {
			s.hasColors.push_back(false);
		}


		GDALRasterAttributeTable *rat = poBand->GetDefaultRAT();
		if( rat != NULL )	{
			s.hasAttributes.push_back(true);
			s.atts.resize(i+1);
			s.atts[i] = GetRATdf(rat);;
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
				s.names.push_back(varname + "_" + std::to_string(i+1) ) ;
			} else {
				s.names.push_back(basename_noext(fname)) ;
			}
		}
	}

	if (gdrv == "netCDF") {
		std::vector<std::string> metadata;
		char **m = poDataset->GetMetadata();
		while (*m != nullptr) {
			metadata.push_back(*m++);
		}
		std::string msg;
		s.set_names_time_ncdf(metadata, bandmeta, msg);
		if (msg.size() > 1) {
			addWarning(msg);
		}
	}

	GDALClose( (GDALDatasetH) poDataset );
	s.hasValues = true;
	setSource(s);

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


void SpatRaster::readChunkGDAL(std::vector<double> &data, unsigned src, size_t row, unsigned nrows, size_t col, unsigned ncols) {

	if (source[src].hasWindow) { // ignoring the expanded case.
		row = row + source[src].window.off_row;
		col = col + source[src].window.off_col;
	}

	std::vector<double> errout;
	if (source[src].rotated) {
		setError("cannot read from rotated files. First use 'rectify'");
		return;
	}

	if (!source[src].open_read) {
		setError("the file is not open for reading");
		return;
	}

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


	if (panBandMap.size() > 0) {
		err = source[src].gdalconnection->RasterIO(GF_Read, col, row, ncols, nrows, &out[0], ncols, nrows, GDT_Float64, nl, &panBandMap[0], 0, 0, 0, NULL);
	} else {
		err = source[src].gdalconnection->RasterIO(GF_Read, col, row, ncols, nrows, &out[0], ncols, nrows, GDT_Float64, nl, NULL, 0, 0, 0, NULL);	
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





std::vector<double> SpatRaster::readValuesGDAL(unsigned src, size_t row, size_t nrows, size_t col, size_t ncols, int lyr) {


	std::vector<double> errout;
	if (source[src].rotated) {
		setError("cannot read from rotated files. First use 'rectify'");
		return errout;
	}

	if (source[src].hasWindow) { // ignoring the expanded case.
		row = row + source[src].window.off_row;
		col = col + source[src].window.off_col;
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
	unsigned nl;
	std::vector<int> panBandMap;
	if (lyr < 0) {
		nl = source[src].nlyr;
		if (!source[src].in_order()) {
			panBandMap.reserve(nl);
			for (size_t i=0; i < nl; i++) {
				panBandMap.push_back(source[src].layers[i]+1);
			}
		}
	} else {
		nl = 1;
		panBandMap.push_back(lyr+1);
	}

	std::vector<double> out(ncell*nl);
	int hasNA;
	std::vector<double> naflags(nl, NAN);
	CPLErr err = CE_None;
	if (panBandMap.size() > 0) {
		err = poDataset->RasterIO(GF_Read, col, row, ncols, nrows, &out[0], ncols, nrows, GDT_Float64, nl, &panBandMap[0], 0, 0, 0, NULL);
	} else {
		err = poDataset->RasterIO(GF_Read, col, row, ncols, nrows, &out[0], ncols, nrows, GDT_Float64, nl, NULL, 0, 0, 0, NULL);
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



std::vector<double> SpatRaster::readGDALsample(unsigned src, size_t srows, size_t scols) {

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
	
	if (panBandMap.size() > 0) {
		err = poDataset->RasterIO(GF_Read, col, row, ncols, nrows, &out[0], scols, srows, GDT_Float64, nl, &panBandMap[0], 0, 0, 0, NULL);
	} else {
		err = poDataset->RasterIO(GF_Read, col, row, ncols, nrows, &out[0], scols, srows, GDT_Float64, nl, NULL, 0, 0, 0, NULL);
	}
	if (err == CE_None ) { 
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



std::vector<std::vector<double>> SpatRaster::readRowColGDAL(unsigned src, std::vector<int_64> &rows, const std::vector<int_64> &cols) {

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
		if ((cols[j] < 0) || (rows[j] < 0)) continue;
		if (panBandMap.size() > 0) {
			err = poDataset->RasterIO(GF_Read, cols[j], rows[j], 1, 1, &out[j*nl], 1, 1, GDT_Float64, nl, &panBandMap[0], 0, 0, 0, NULL);
		} else {
			err = poDataset->RasterIO(GF_Read, cols[j], rows[j], 1, 1, &out[j*nl], 1, 1, GDT_Float64, nl, NULL, 0, 0, 0, NULL);
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

