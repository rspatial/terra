// Copyright (c) 2018-2019  Robert J. Hijmans
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
#include "math_utils.h"
#include "string_utils.h"
#include "file_utils.h"

#include <unordered_map>

#include "gdal_priv.h"
#include "cpl_conv.h" // for CPLMalloc()
#include "cpl_string.h"
#include "ogr_spatialref.h"



void GDALformat(std::string &filename, std::string &format) {

	lrtrim(format);
	lrtrim(filename);
	std::string ext = getFileExt(filename);
    lowercase(ext);
	
	std::unordered_map<std::string, std::string>
	drivers = {
		{".tif","GTiff"}, {".tiff","GTiff"},
		{".nc","netCDF"}, {".cdf","netCDF"}, {".ncdf","netCDF"},
		{".img","HFA"},
		{".flt","EHdr"},
		{".grd","RRASTER"},
		{".sgrd","SAGA"}, {".sdat","SAGA"},
		{".bil","BIL"},
		{".bsq","BSQ"},
		{".bip","BIP"},
		{".rst","RST"},
		{".envi","ENVI"},
		{".asc","AAIGrid"}
	};

    auto i = drivers.find(ext);
    if (i == drivers.end()) {
		format = "GTiff";
	} else {
		format = i->second;
	}
}


bool SpatRaster::writeStartGDAL(std::string filename, std::string format, std::string datatype, bool overwrite) {

	SpatMessages m = can_write(filename, overwrite);
	if (m.has_error) {
		msg = m;
		return(false);
	}
	
	GDALformat(filename, format);
	const char *pszFormat = format.c_str();
	const char *pszDstFilename = filename.c_str();
    GDALDriver *poDriver;
    char **papszMetadata;
	GDALAllRegister();
	
    poDriver = GetGDALDriverManager()->GetDriverByName(pszFormat);
    if(poDriver == NULL) {
		setError("driver failure");
		return (false);
	}
    papszMetadata = poDriver->GetMetadata();
    if(! CSLFetchBoolean( papszMetadata, GDAL_DCAP_CREATE, FALSE)) return (false);
 
	GDALDataset *poDstDS;
	char **papszOptions = NULL;
	
	poDstDS = poDriver->Create( pszDstFilename, ncol(), nrow(), nlyr(), GDT_Float64, papszOptions);

	std::vector<double> rs = resolution();
	double adfGeoTransform[6] = { extent.xmin, rs[0], 0, extent.ymax, 0, -1 * rs[1] };
	poDstDS->SetGeoTransform(adfGeoTransform);
	
	std::string prj = getCRS();
	OGRSpatialReference oSRS;
	OGRErr erro = oSRS.importFromProj4(&prj[0]); 
	if (erro == 4) { 
		setError("CRS failure");
		return false ;
	}
	char *pszSRS_WKT = NULL;	
	oSRS.exportToWkt(&pszSRS_WKT);
	poDstDS->SetProjection(pszSRS_WKT);
	CPLFree(pszSRS_WKT);

	source[0].resize(nlyr());
	source[0].gdalconnection = poDstDS;
	source[0].datatype = datatype;
	for (size_t i =0; i<nlyr(); i++) {
		source[0].range_min[i] = std::numeric_limits<double>::max();
		source[0].range_max[i] = std::numeric_limits<double>::lowest();
	}
	source[0].driver = "gdal" ;

	return true;
}

/*
bool SpatRaster::writeValuesGDAL(std::vector<double> vals, unsigned row){
	unsigned nrows = vals.size() / (nlyr() * ncol());
	unsigned start;
	CPLErr err = CE_None;
	GDALRasterBand *poBand;
	double vmin, vmax;
	unsigned nc = nrows * ncol();
	GDALDataType gdtype;
	for (size_t i=0; i < nlyr(); i++) {
		start = nc * i;

		std::string datatype = source[0].datatype;
		if (datatype == "FLT8S") {
			gdtype = GDT_Float64;
			std::vector<double> vv(vals.begin(), vals.end());
			poBand = source[0].gdalconnection->GetRasterBand(i+1);
			err = poBand->RasterIO(GF_Write, 0, row, ncol(), nrows, &vv[start], ncol(), nrows, gdtype, 0, 0 );
			if (err == 4) break;
		} else if (datatype == "FLT4S") {
			gdtype = GDT_Float32;
			std::vector<float> vv(vals.begin(), vals.end());
			poBand = source[0].gdalconnection->GetRasterBand(i+1);
			err = poBand->RasterIO(GF_Write, 0, row, ncol(), nrows, &vv[start], ncol(), nrows, gdtype, 0, 0 );
			if (err == 4) break;
			//std::cout <<  "\n" << vv[0] << "\n";
		} else if (datatype == "INT4S") {
			gdtype = GDT_Int32;			
			std::vector<int32_t> vv(vals.begin(), vals.end());
			poBand = source[0].gdalconnection->GetRasterBand(i+1);
			err = poBand->RasterIO(GF_Write, 0, row, ncol(), nrows, &vv[start], ncol(), nrows, gdtype, 0, 0 );
			if (err == 4) break;
		} else if (datatype == "INT2S") {
			gdtype = GDT_Int16;
			std::vector<int16_t> vv(vals.begin(), vals.end());
			poBand = source[0].gdalconnection->GetRasterBand(i+1);
			err = poBand->RasterIO(GF_Write, 0, row, ncol(), nrows, &vv[start], ncol(), nrows, gdtype, 0, 0 );
			if (err == 4) break;
		} else if (datatype == "INT4U") {
			gdtype = GDT_UInt32;			
			std::vector<uint32_t> vv(vals.begin(), vals.end());
			poBand = source[0].gdalconnection->GetRasterBand(i+1);
			err = poBand->RasterIO(GF_Write, 0, row, ncol(), nrows, &vv[start], ncol(), nrows, gdtype, 0, 0 );
			if (err == 4) break;
		} else if (datatype == "INT2U") {
			gdtype = GDT_UInt16;			
			std::vector<uint16_t> vv(vals.begin(), vals.end());
			poBand = source[0].gdalconnection->GetRasterBand(i+1);
			err = poBand->RasterIO(GF_Write, 0, row, ncol(), nrows, &vv[start], ncol(), nrows, gdtype, 0, 0 );
			if (err == 4) break;
		} else if (datatype == "INT1U") {
			gdtype = GDT_Byte;		
			std::vector<int8_t> vv(vals.begin(), vals.end());
			poBand = source[0].gdalconnection->GetRasterBand(i+1);
			err = poBand->RasterIO(GF_Write, 0, row, ncol(), nrows, &vv[start], ncol(), nrows, gdtype, 0, 0 );
			if (err == 4) break;
		}	

		minmax(vals.begin()+start, vals.begin()+start+nc, vmin, vmax);
		source[0].range_min[i] = std::min(source[0].range_min[i], vmin);
		source[0].range_max[i] = std::max(source[0].range_max[i], vmax); 
		
	}
	
	if (err != CE_None ) {
		setError("cannot write values");
		return false;
	}	
	return true;
}
*/


bool SpatRaster::fillValuesGDAL(double fillvalue) {
	CPLErr err = CE_None;
	GDALRasterBand *poBand;
	for (size_t i=0; i < nlyr(); i++) {
		poBand = source[0].gdalconnection->GetRasterBand(i+1);
		err = poBand->Fill(fillvalue);
	}
	if (err != CE_None ) {
		setError("cannot fill values");
		return false;
	}	
	return true;
}



bool SpatRaster::writeValuesGDAL(std::vector<double> &vals, unsigned startrow, unsigned nrows, unsigned startcol, unsigned ncols){
	CPLErr err = CE_None;
	GDALRasterBand *poBand;
	double vmin, vmax;
	unsigned nc = nrows * ncols;
	for (size_t i=0; i < nlyr(); i++) {
		unsigned start = nc * i;
		std::string datatype = source[0].datatype;
		poBand = source[0].gdalconnection->GetRasterBand(i+1);
		if (datatype == "FLT8S") {
			//std::vector<double> vv(vals.begin(), vals.end());
			err = poBand->RasterIO(GF_Write, startcol, startrow, ncols, nrows, &vals[start], ncols, nrows, GDT_Float64, 0, 0 );
		} else if (datatype == "FLT4S") {
			std::vector<float> vv(vals.begin(), vals.end());
			err = poBand->RasterIO(GF_Write, startcol, startrow, ncols, nrows, &vv[start], ncols, nrows, GDT_Float32, 0, 0 );
		} else if (datatype == "INT4S") {
			std::vector<int32_t> vv(vals.begin(), vals.end());
			err = poBand->RasterIO(GF_Write, startcol, startrow, ncols, nrows, &vv[start], ncols, nrows, GDT_Int32, 0, 0 );
		} else if (datatype == "INT2S") {
			std::vector<int16_t> vv(vals.begin(), vals.end());
			err = poBand->RasterIO(GF_Write, startcol, startrow, ncols, nrows, &vv[start], ncols, nrows, GDT_Int16, 0, 0 );
		} else if (datatype == "INT4U") {
			std::vector<uint32_t> vv(vals.begin(), vals.end());
			err = poBand->RasterIO(GF_Write, startcol, startrow, ncols, nrows, &vv[start], ncols, nrows, GDT_UInt32, 0, 0 );
		} else if (datatype == "INT2U") {
			std::vector<uint16_t> vv(vals.begin(), vals.end());
			err = poBand->RasterIO(GF_Write, startcol, startrow, ncols, nrows, &vv[start], ncols, nrows, GDT_UInt16, 0, 0 );
		} else if (datatype == "INT1U") {
			std::vector<int8_t> vv(vals.begin(), vals.end());
			err = poBand->RasterIO(GF_Write, startcol, startrow, ncols, nrows, &vv[start], ncols, nrows, GDT_Byte, 0, 0 );
		}	
		if (err == 4) break;

		minmax(vals.begin()+start, vals.begin()+start+nc, vmin, vmax);
		source[0].range_min[i] = std::min(source[0].range_min[i], vmin);
		source[0].range_max[i] = std::max(source[0].range_max[i], vmax); 
		
	}
	
	if (err != CE_None ) {
		setError("cannot write values");
		return false;
	}	
	return true;
}


bool SpatRaster::writeStopGDAL() {
	GDALRasterBand *poBand;
	source[0].hasRange.resize(nlyr());
	for (size_t i=0; i < nlyr(); i++) {
		poBand = source[0].gdalconnection->GetRasterBand(i+1);
		poBand->SetStatistics(source[0].range_min[i], source[0].range_max[i], -9999., -9999.);
		source[0].hasRange[i] = true;
	}
	GDALClose( (GDALDatasetH) source[0].gdalconnection );
	return true;
}




bool SpatRaster::writeRasterGDAL(std::string filename, std::string format, std::string datatype, bool overwrite) {
	SpatRaster r = geometry();
	bool values = true;
	if (!hasValues()) {
		addWarning("there are no cell values");
		values = false;
	}
	if (!r.writeStartGDAL(filename, format, datatype, overwrite)) {
		return false;
	}
	if (values) {
		std::vector<double> v = getValues();
		if (!r.writeValuesGDAL(v, 0, nrow(), 0, ncol())) {
			return false;
		}
	}
	if (!r.writeStopGDAL()) {
		setError("cannot close file");
		return false;
	}
	return true;
}


