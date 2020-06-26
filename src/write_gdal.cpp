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
#include "math_utils.h"
#include "string_utils.h"
#include "file_utils.h"

#include <unordered_map>

#include "gdal_priv.h"
#include "cpl_conv.h" // for CPLMalloc()
#include "cpl_string.h"
#include "ogr_spatialref.h"



void getGDALdriver(std::string &filename, std::string &driver) {

	lrtrim(driver);
	if (driver != "") return;
	
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
    if (i != drivers.end()) {
		driver = i->second;
	}
}


bool getGDALDataType(std::string datatype, GDALDataType &gdt) {
	if (datatype=="FLT4S") {
		gdt = GDT_Float32;
	} else if (datatype == "INT4S") {
		gdt = GDT_Int32;
	} else if (datatype == "FLT8S") {
		gdt = GDT_Float64;
	} else if (datatype == "INT2S") {
		gdt = GDT_Int16;
	} else if (datatype == "INT4U") {
		gdt = GDT_UInt32;
	} else if (datatype == "INT2U") {
		gdt = GDT_UInt16;
	} else if (datatype == "INT1U") {
		gdt = GDT_Byte;
	} else {
		gdt = GDT_Float32;
		return false;
	}
	return true;
}



CPLErr setBandCategories(GDALRasterBand *poBand, std::vector<std::string> cats) {
	char **names = NULL;
	for (size_t i = 0; i < cats.size(); i++) {
		names = CSLAddString(names, cats[i].c_str());
	}
	CPLErr err = poBand->SetCategoryNames(names);
	return err;
}


//#include <iostream>
//#include "Rcpp.h"

bool SpatRaster::writeStartGDAL(std::string filename, std::string driver, std::string datatype, bool overwrite, SpatOptions &opt) {

	std::string errmsg;
	if (!can_write(filename, overwrite, errmsg)) {
		setError(errmsg);
		return(false);
	}

	getGDALdriver(filename, driver);
	if (driver == "") {
		setError("cannot guess file type from filename");
		return(false);	
	}
	
	const char *pszFormat = driver.c_str();
	const char *pszDstFilename = filename.c_str();
    GDALDriver *poDriver;
    char **papszMetadata;
	//GDALAllRegister();

    poDriver = GetGDALDriverManager()->GetDriverByName(pszFormat);
    if(poDriver == NULL) {
		setError("driver failure");
		return (false);
	}
    papszMetadata = poDriver->GetMetadata();
    if(! CSLFetchBoolean( papszMetadata, GDAL_DCAP_CREATE, FALSE)) return (false);

	GDALDataset *poDstDS;
	char **papszOptions = NULL;

	for (size_t i=0; i<opt.gdal_options.size(); i++) {
		std::vector<std::string> gopt = strsplit(opt.gdal_options[i], "=");
		if (gopt.size() == 2) {
			papszOptions = CSLSetNameValue( papszOptions, gopt[0].c_str(), gopt[1].c_str() );
		}
	}

//	Rcpp::Rcout << datatype << std::endl;
	
	GDALDataType gdt;
	if (!getGDALDataType(datatype, gdt)) {
		addWarning("unknown datatype = " + datatype);
	}

	poDstDS = poDriver->Create( pszDstFilename, ncol(), nrow(), nlyr(), gdt, papszOptions);

	CSLDestroy( papszOptions );

	if (opt.names.size() == nlyr()) {
		setNames(opt.names);
	}
	GDALRasterBand *poBand;
	std::vector<std::string> nms = getNames();
	for (size_t i=0; i < nlyr(); i++) {
		poBand = poDstDS->GetRasterBand(i+1);
		poBand->SetDescription(nms[i].c_str());
		if ((i==0) || (driver != "GTiff")) {
			// to avoid "Setting nodata to nan on band 2, but band 1 has nodata at nan." 
			if (datatype == "INT4S") {
				poBand->SetNoDataValue(INT32_MIN); //-2147483648; 
			} else if (datatype == "INT2S") {
				poBand->SetNoDataValue(INT16_MIN); 
			} else {
				poBand->SetNoDataValue(NAN); 
			}
		}
	}
	std::vector<double> rs = resolution();
	double adfGeoTransform[6] = { extent.xmin, rs[0], 0, extent.ymax, 0, -1 * rs[1] };
	poDstDS->SetGeoTransform(adfGeoTransform);

	std::string crs = srs.wkt;
	OGRSpatialReference oSRS;
	OGRErr erro = oSRS.SetFromUserInput(&crs[0]);
	if (erro == 4) {
		setError("CRS failure");
		return false ;
	}
	char *pszSRS_WKT = NULL;
	oSRS.exportToWkt(&pszSRS_WKT);
	poDstDS->SetProjection(pszSRS_WKT);
	CPLFree(pszSRS_WKT);

	source[0].resize(nlyr());
	source[0].nlyrfile = nlyr();
	source[0].gdalconnection = poDstDS;
	source[0].datatype = datatype;
	for (size_t i =0; i<nlyr(); i++) {
		source[0].range_min[i] = std::numeric_limits<double>::max();
		source[0].range_max[i] = std::numeric_limits<double>::lowest();
	}
	source[0].driver = "gdal" ;
	source[0].filename = filename;
	source[0].memory = false;
	return true;
}

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
	//GDALRasterBand *poBand;
	double vmin, vmax;
	unsigned nc = nrows * ncols;
	size_t nl = nlyr();
	//for (size_t i=0; i < nl; i++) {
	//unsigned start = nc * i;
	std::string datatype = source[0].datatype;
	//poBand = source[0].gdalconnection->GetRasterBand(i+1);

	for (size_t i=0; i < nl; i++) {
		unsigned start = nc * i;
		minmax(vals.begin()+start, vals.begin()+start+nc, vmin, vmax);
		source[0].range_min[i] = std::min(source[0].range_min[i], vmin);
		source[0].range_max[i] = std::max(source[0].range_max[i], vmax);
	}

	if (datatype == "FLT8S") {
			//std::vector<double> vv(vals.begin(), vals.end());
		err = source[0].gdalconnection->RasterIO(GF_Write, startcol, startrow, ncols, nrows, &vals[0], ncols, nrows, GDT_Float64, nl, NULL, 0, 0, 0, NULL );
	} else if (datatype == "FLT4S") {
		std::vector<float> vv(vals.begin(), vals.end());
		err = source[0].gdalconnection->RasterIO(GF_Write, startcol, startrow, ncols, nrows, &vv[0], ncols, nrows, GDT_Float32, nl, NULL, 0, 0, 0, NULL );
	} else if (datatype == "INT4S") {
		std::transform(vals.begin(), vals.end(), vals.begin(),
			[](double v) { return (std::isnan(v) ? double(INT32_MIN) : v); } );
		// std::replace(vals.begin(), vals.end(), NAN, (double) -2147483648); //works
		std::vector<int32_t> vv(vals.begin(), vals.end());
		err = source[0].gdalconnection->RasterIO(GF_Write, startcol, startrow, ncols, nrows, &vv[0], ncols, nrows, GDT_Int32, nl, NULL, 0, 0, 0, NULL );
	} else if (datatype == "INT2S") {
		std::transform(vals.begin(), vals.end(), vals.begin(),
			[](double v) { return (std::isnan(v) ? double(INT16_MIN) : v); } );
		//std::replace(vals.begin(), vals.end(), NAN, -32768); //works not
		std::vector<int16_t> vv(vals.begin(), vals.end());
		err = source[0].gdalconnection->RasterIO(GF_Write, startcol, startrow, ncols, nrows, &vv[0], ncols, nrows, GDT_Int16, nl, NULL, 0, 0, 0, NULL );
	} else if (datatype == "INT4U") {
		std::vector<uint32_t> vv(vals.begin(), vals.end());
		err = source[0].gdalconnection->RasterIO(GF_Write, startcol, startrow, ncols, nrows, &vv[0], ncols, nrows, GDT_UInt32, nl, NULL, 0, 0, 0, NULL );
	} else if (datatype == "INT2U") {
		std::vector<uint16_t> vv(vals.begin(), vals.end());
		err = source[0].gdalconnection->RasterIO(GF_Write, startcol, startrow, ncols, nrows, &vv[0], ncols, nrows, GDT_UInt16, nl, NULL, 0, 0, 0, NULL );
	} else if (datatype == "INT1U") {
		std::vector<int8_t> vv(vals.begin(), vals.end());
		err = source[0].gdalconnection->RasterIO(GF_Write, startcol, startrow, ncols, nrows, &vv[0], ncols, nrows, GDT_Byte, nl, NULL, 0, 0, 0, NULL );
	}
	//if (err == 4) break;

	if (err != CE_None ) {
		setError("cannot write values");
		return false;
	}
	return true;
}


bool SpatRaster::writeStopGDAL() {
	GDALRasterBand *poBand;
	source[0].hasRange.resize(nlyr());
	std::string datatype = source[0].datatype;

	for (size_t i=0; i < nlyr(); i++) {
		poBand = source[0].gdalconnection->GetRasterBand(i+1);
		if (datatype == "INT4S") {
			source[0].range_min[i] = (long) source[0].range_min[i]; 
			source[0].range_max[i] = (long) source[0].range_max[i]; 
		} else if (datatype == "INT2S") {
			source[0].range_min[i] = (int) source[0].range_min[i]; 
			source[0].range_max[i] = (int) source[0].range_max[i]; 
		}
		poBand->SetStatistics(source[0].range_min[i], source[0].range_max[i], -9999., -9999.);
		source[0].hasRange[i] = true;
	}
	GDALClose( (GDALDatasetH) source[0].gdalconnection );
	source[0].hasValues = true;
	return true;
}




SpatRaster SpatRaster::writeRasterGDAL(std::string filename, std::string format, std::string datatype, bool overwrite, SpatOptions &opt) {
	SpatRaster out = geometry();
	bool values = true;
	if (!hasValues()) {
		addWarning("there are no cell values");
		values = false;
	}
	if (!out.writeStartGDAL(filename, format, datatype, overwrite, opt)) {
		return out;
	}
	if (values) {
		std::vector<double> v = getValues();
		if (!out.writeValuesGDAL(v, 0, nrow(), 0, ncol())) {
			return out;
		}
	}
	if (!out.writeStopGDAL()) {
		setError("cannot close file");
		return out;
	}
	return out;
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

