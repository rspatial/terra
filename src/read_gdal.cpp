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

#include "spatRaster.h"
#include "file_utils.h"
#include "string_utils.h"
#include "NA.h"

#include "gdal_priv.h"
#include "cpl_conv.h" // for CPLMalloc()
#include "cpl_string.h"
#include "ogr_spatialref.h"

#include "gdal_rat.h"
//#include "hdr.h"


void SpatRaster::spatinit() {
    GDALAllRegister();
    OGRRegisterAll(); // should go to SpatVector
	//GDALregistred = true;
}





bool SpatRaster::constructFromFiles(std::vector<std::string> fnames) {

	SpatRaster r = SpatRaster(fnames[0]);
	setSource(r.source[0]);
	for (size_t i=1; i<fnames.size(); i++) {
		r = SpatRaster(fnames[i]);
		if (!compare_geom(r, false, true, true)) {
			setError("geometry of " + fnames[i] + " does not match previous sources");
			return false;
		} else {
			addSource(r);
		}
	}
	return true;
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

bool SpatRaster::constructFromSubDataSets(std::string filename, std::vector<std::string> sds) {

	std::vector<std::string> sd; //, nms;
	std::string delim = "NAME=";
	for (size_t i=0; i<sds.size(); i++) {
		std::string s = sds[i];
		size_t pos = s.find(delim);
		if (pos != std::string::npos) {
			s.erase(0, pos + delim.length());
			sd.push_back(s);
		} 
		//else {
			// _DESC=
		//	s.erase(0, pos + delim.length());
			//nms.push_back(s);
			//printf( "%s\n", s.c_str() );

		//}
	}

	constructFromFile(sd[0]);
	SpatRaster out;
	bool success;
    for (size_t i=1; i < sd.size(); i++) {
//		printf( "%s\n", sd[i].c_str() );
		success = out.constructFromFile(sd[i]);
		if (success) {
//			out.source[0].subdataset = true;
			addSource(out);
			if (out.msg.has_error) {
				setError(out.msg.error);
				return false;
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


bool SpatRaster::constructFromFile(std::string fname) {

    GDALDataset *poDataset;
    
	//if (!GDALregistred) spatinit(); //
	GDALAllRegister();

	const char* pszFilename = fname.c_str();
    poDataset = (GDALDataset *) GDALOpen( pszFilename, GA_ReadOnly );

    if( poDataset == NULL )  {
		if (!file_exists(fname)) {
			setError("file does not exist");
		} else {
			setError("cannot read from " + fname );
		}
		return false;
	}

	unsigned nl = poDataset->GetRasterCount();

	if (nl == 0) {
		std::vector<std::string> meta;
		char **metadata = poDataset->GetMetadata("SUBDATASETS");
	    for (size_t i=0; metadata[i] != NULL; i++) {
			meta.push_back(metadata[i]);
		}
		if (meta.size() > 0) {
			return constructFromSubDataSets(fname, meta);
		}
	}

	RasterSource s;
	s.ncol = poDataset->GetRasterXSize();
	s.nrow = poDataset->GetRasterYSize();
	s.nlyr = nl;
	s.nlyrfile = nl;
	s.layers.resize(nl);
    std::iota(s.layers.begin(), s.layers.end(), 0);

	double adfGeoTransform[6];
	if( poDataset->GetGeoTransform( adfGeoTransform ) == CE_None ) {
		// the rounding below is to address a design flaw in GDAL
		// GDAL provides the coordinates of one corner and the resolution,
		// instead of the coordinates of all (two opposite) corners.
		// This makes computation of the opposite corner coordinates only
		// approximate for large rasters with a high resolution.
		double xmin = adfGeoTransform[0]; /* left x */
		double xmax = xmin + adfGeoTransform[1] * s.ncol; /* w-e pixel resolution */
		//xmax = roundn(xmax, 9);
		double ymax = adfGeoTransform[3]; /* top y */
		double ymin = ymax + s.nrow * adfGeoTransform[5]; /* n-s pixel resolution (negative value) */
		//ymin = roundn(ymin, 9);
		SpatExtent e(xmin, xmax, ymin, ymax);
		s.extent = e;
	}

	s.memory = false;
	s.filename = fname;
	s.driver = "gdal";


	if( poDataset->GetProjectionRef() != NULL ) {
		OGRSpatialReference oSRS(poDataset->GetProjectionRef());
		char *pszPRJ = NULL;
		oSRS.exportToProj4(&pszPRJ);
		s.crs = pszPRJ;
	} else {
		s.crs = "";
	}


	GDALRasterBand  *poBand;
	//int nBlockXSize, nBlockYSize;
	double adfMinMax[2];
	int bGotMin, bGotMax;

//	s.layers.resize(1);
	for (size_t i = 0; i < s.nlyr; i++) {
		poBand = poDataset->GetRasterBand(i+1);

//		source.layers[0].push_back(i+1);
		//poBand->GetBlockSize( &nBlockXSize, &nBlockYSize );


//	double GDALRasterBand::GetOffset 	( 	int *  	pbSuccess = nullptr	)
//	double GDALRasterBand::GetScale 	( 	int *  	pbSuccess = nullptr	)


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
				s.names.push_back(basename_noext(fname) + std::to_string(i+1) ) ;
			} else {
				s.names.push_back(basename_noext(fname)) ;
			}
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
	GDALClose( (GDALDatasetH) source[src].gdalconnection);
	source[src].open_read = false;
	return true;
}

std::vector<double> SpatRaster::readChunkGDAL(unsigned src, unsigned row, unsigned nrows, unsigned col, unsigned ncols) {
	GDALRasterBand  *poBand;
	unsigned ncell = ncols * nrows;
	unsigned cell;
	std::vector<double> errout;
	unsigned nl = source[src].nlyr;
	std::vector<double> out(ncell * nl);
	int hasNA;
	CPLErr err = CE_None;
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
		} else {

		}
	}
	if (err != CE_None ) {
		setError("cannot read values");
		return errout;
	}
	return(out);
}




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



/*
void applyScaleOffset(std::vector<double> &d, double scale, double offset, bool haveso) {
	if (haveso) {
		for (size_t i=0; i<d.size(); i++) {
			d[i] = d[i] * scale + offset;
		}
	}
}
*/


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

