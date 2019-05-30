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


bool SpatRaster::writeStartGDAL(std::string filename, std::string format, std::string datatype) {

	GDALformat(filename, format);
	const char *pszFormat = format.c_str();
	const char *pszDstFilename = filename.c_str();
    GDALDriver *poDriver;
    char **papszMetadata;
	GDALAllRegister();
	
    poDriver = GetGDALDriverManager()->GetDriverByName(pszFormat);
    if(poDriver == NULL) return (false);
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
	if (erro == 4) { return false ; }	// ??
	
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

	return true;
}


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
	bool success;
	SpatRaster r = geometry();
	
	if (!hasValues()) {
		addWarning("none of the cells have values");
	}
	success = r.writeStartGDAL(filename, format, datatype);
	if (!success) {
		setError("cannot open file");
		return false;
	}
	success = r.writeValuesGDAL(getValues(), 0);
	if (!success) {
		setError("cannot write values to file");
		return false;
	}
	success = r.writeStopGDAL();
	if (!success) {
		setError("cannot close file");
		return false;
	}
	return true;
}


