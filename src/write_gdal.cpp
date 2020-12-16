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

#include "gdalio.h"


void getGDALdriver(std::string &filename, std::string &driver) {

	lrtrim(filename);
	lrtrim(driver);

	if (driver != "") {
		if (driver == "RST") {
			filename = noext(filename) + ".rst";
		}
	
		return;
	}
	
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
		{".rst","RST"},
		{".envi","ENVI"},
		{".asc","AAIGrid"}
	};

    auto i = drivers.find(ext);
    if (i != drivers.end()) {
		driver = i->second;
	}
}



bool setCats(GDALRasterBand *poBand, SpatCategories &cats) {
	char **names = NULL;
	for (size_t i = 0; i < cats.labels.size(); i++) {
		names = CSLAddString(names, cats.labels[i].c_str());
	}
	CPLErr err = poBand->SetCategoryNames(names);
	return (err == CE_None);
}


bool setCT(GDALRasterBand *poBand, SpatDataFrame &d) {
	CPLErr err = poBand->SetColorInterpretation(GCI_PaletteIndex);
	GDALColorTable *poCT = new GDALColorTable(GPI_RGB);
	GDALColorEntry col;
	for (size_t j=0; j< d.nrow(); j++) {
		col.c1 = (short)d.iv[0][j];
		col.c2 = (short)d.iv[1][j];
		col.c3 = (short)d.iv[2][j];
		col.c4 = (short)d.iv[3][j];
		poCT->SetColorEntry(j, &col);
	}
	err = poBand->SetColorTable(poCT);
	delete poCT;
	return (err == CE_None);
}



bool SpatRaster::writeStartGDAL(SpatOptions &opt) {


	std::string filename = opt.get_filename();
	if (filename == "") {
		setError("empty filename");
		return(false);
	} else {
		// make sure filename won't be used again
		opt.set_filenames({""});
	}

	std::string driver = opt.get_filetype();
	getGDALdriver(filename, driver);
	if (driver == "") {
		setError("cannot guess file type from filename");
		return(false);	
	}
	if (driver == "AAIGrid" && nlyr() > 1) {
		setError("AAIGrid can only have one layer");
		return false;
	}	
	//if (driver == "netCDF") {
	//	setError("netCDF writing is only supported through 'writeCDF'");
	//	return false;
	//}	

	std::string errmsg;
	if (!can_write(filename, opt.get_overwrite(), errmsg)) {
		setError(errmsg);
		return(false);
	}

		
	//std::string ext = getFileExt(filename);
	//lowercase(ext);
	std::string datatype = opt.get_datatype();
	std::vector<bool> hasCats = hasCategories();
	std::vector<bool> hasCT = hasColors();
	std::vector<SpatDataFrame> ct = getColors();
	if (hasCT[0] || hasCats[0]) { 
		// must be INT1U for color table with gtiff
		// perhaps do not do this if datatype was explicitly set by user
		datatype = "INT1U";
	} else if (datatype != "INT1U") {
		std::fill(hasCT.begin(), hasCT.end(), false);
	}

	GDALDataType gdt;
	if (!getGDALDataType(datatype, gdt)) {
		setError("invalid datatype");
		return false;
		//addWarning("unknown datatype = " + datatype + ". Set to FLT4S");
		// datatype = "FLT4S"
	}
	source[0].datatype = datatype;
	
	int dsize = std::stoi(datatype.substr(3,1));
	GIntBig diskNeeded = ncell() * nlyr() * dsize;
	std::string dname = dirname(filename);
	GIntBig diskAvailable = VSIGetDiskFreeSpace(dname.c_str());
	if ((diskAvailable > -1) && (diskAvailable < diskNeeded)) {
		setError("insufficient disk space (perhaps from temporary file)");
		return(false);			
	}

    #ifdef useRcpp
	if (opt.verbose) {
		double gb = 1073741824 / 8;
		Rcpp::Rcout<< "filename      : " << filename << std::endl;
		Rcpp::Rcout<< "driver        : " << driver   << std::endl;
		if (diskAvailable > 0) {
			Rcpp::Rcout<< "disk available: " << roundn(diskAvailable / gb, 1) << " GB" << std::endl;
		}
		Rcpp::Rcout<< "disk needed   : " << roundn(diskNeeded / gb, 1) << " GB" << std::endl;
	}
	#endif

    GDALDriver *poDriver;
    poDriver = GetGDALDriverManager()->GetDriverByName(driver.c_str());
    if(poDriver == NULL) {
		setError("invalid driver");
		return (false);
	}

	GDALDataset *poDS;
	char **papszOptions = NULL;
	for (size_t i=0; i<opt.gdal_options.size(); i++) {
		std::vector<std::string> gopt = strsplit(opt.gdal_options[i], "=");
		if (gopt.size() == 2) {
			papszOptions = CSLSetNameValue( papszOptions, gopt[0].c_str(), gopt[1].c_str() );
		}
	}

    char **papszMetadata;
    papszMetadata = poDriver->GetMetadata();

	//bool isncdf = ((driver == "netCDF" && opt.get_ncdfcopy()));

    if (CSLFetchBoolean( papszMetadata, GDAL_DCAP_CREATE, FALSE)) {
		poDS = poDriver->Create(filename.c_str(), ncol(), nrow(), nlyr(), gdt, papszOptions);
	} else if (CSLFetchBoolean( papszMetadata, GDAL_DCAP_CREATECOPY, FALSE)) {
		copy_driver = driver;
		if (canProcessInMemory(opt)) {
			poDriver = GetGDALDriverManager()->GetDriverByName("MEM");
			poDS = poDriver->Create("", ncol(), nrow(), nlyr(), gdt, papszOptions);
		} else {
			std::string f = tempFile(opt.get_tempdir(), ".tif");
			copy_filename = f;
			poDriver = GetGDALDriverManager()->GetDriverByName("GTiff");
			poDS = poDriver->Create(f.c_str(), ncol(), nrow(), nlyr(), gdt, papszOptions);
		}
	} else {
		setError("cannot create this format: "+ driver);
		return false;
	}
	CSLDestroy( papszOptions );

	if (opt.names.size() == nlyr()) {
		setNames(opt.names);
	}
	GDALRasterBand *poBand;
	std::vector<std::string> nms = getNames();
	double naflag=NAN; 
	bool hasNAflag = opt.has_NAflag(naflag);

	for (size_t i=0; i < nlyr(); i++) {

		poBand = poDS->GetRasterBand(i+1);

		if (hasCT[i]) {
			if (!setCT(poBand, ct[i])) {
				addWarning("could not write the color table");
			}
		}
		if (hasCats[i]) {
			SpatCategories cats = getLayerCategories(i);
			if (!setCats(poBand, cats)) {
				addWarning("could not write categories");
			}
		}
		/*
		if (isncdf) {
			std::string opt = "NETCDF_VARNAME";
			char ** papszMetadata; 
			papszMetadata = CSLSetNameValue( papszOptions, opt.c_str(), nms[i].c_str() );
			poBand->SetMetadata(papszMetadata);

		} else {
		*/	
		poBand->SetDescription(nms[i].c_str());
		
		if ((i==0) || (driver != "GTiff")) {
			// to avoid "Setting nodata to nan on band 2, but band 1 has nodata at nan." 
			if (hasNAflag) {
				poBand->SetNoDataValue(naflag); 
			} else if (datatype == "INT4S") {
				poBand->SetNoDataValue(INT32_MIN); //-2147483648; 
			} else if (datatype == "INT2S") {
				poBand->SetNoDataValue(INT16_MIN); 
			} else if (datatype == "INT4U") {
				double na = (double)INT32_MAX * 2 - 1;
				poBand->SetNoDataValue(na); 
			} else if (datatype == "INT2U") {
				double na = (double)INT16_MAX * 2 - 1;
				poBand->SetNoDataValue(na); 
			} else if (datatype == "INT1U") {
				poBand->SetNoDataValue(255); 
			} else {
				poBand->SetNoDataValue(NAN); 
			}
		}
	}

	std::vector<double> rs = resolution();
	SpatExtent extent = getExtent();
	double adfGeoTransform[6] = { extent.xmin, rs[0], 0, extent.ymax, 0, -1 * rs[1] };
	poDS->SetGeoTransform(adfGeoTransform);
	std::string crs = source[0].srs.wkt;
	OGRSpatialReference oSRS;
	OGRErr erro = oSRS.SetFromUserInput(&crs[0]);
	if (erro == 4) {
		setError("CRS failure");
		GDALClose( (GDALDatasetH) poDS );
		return false ;
	}
	char *pszSRS_WKT = NULL;
	oSRS.exportToWkt(&pszSRS_WKT);
	poDS->SetProjection(pszSRS_WKT);
	CPLFree(pszSRS_WKT);
	// destroySRS(oSRS) ?
	
	source[0].gdalconnection = poDS;

	source[0].resize(nlyr());
	source[0].nlyrfile = nlyr();
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


/*
void min_max_na(std::vector<double> &vals, const double &na, const double &mn, const double &mx) {
	for (double &v : vals) { 
		v = std::isnan(v) ? na : (v < mn ? na : (v > mx ? na : v)); 
	}
}
*/

template <typename T>
void tmp_min_max_na(std::vector<T> &out, const std::vector<double> &v, const double &na, const double &mn, const double &mx) {
	size_t n = v.size();
	out.reserve(n);
	for (size_t i=0; i<n; i++) { 
		out.push_back(std::isnan(v[i]) ? na : (v[i] < mn ? na : (v[i] > mx ? na : v[i]))); 
	}
}


bool SpatRaster::writeValuesGDAL(std::vector<double> &vals, uint_64 startrow, uint_64 nrows, uint_64 startcol, uint_64 ncols){

	CPLErr err = CE_None;
	double vmin, vmax;
	uint_64 nc = nrows * ncols;
	size_t nl = nlyr();
	std::string datatype = source[0].datatype;

	
	for (size_t i=0; i < nl; i++) {
		uint_64 start = nc * i;
		minmax(vals.begin()+start, vals.begin()+start+nc, vmin, vmax);
		if (std::isnan(source[0].range_min[i])) {
			source[0].range_min[i] = vmin;
			source[0].range_max[i] = vmax;			
		} else if (!std::isnan(vmin)) {
			source[0].range_min[i] = std::min(source[0].range_min[i], vmin);
			source[0].range_max[i] = std::max(source[0].range_max[i], vmax);
		}
	}

	if ((datatype == "FLT8S") || (datatype == "FLT4S")) {
		err = source[0].gdalconnection->RasterIO(GF_Write, startcol, startrow, ncols, nrows, &vals[0], ncols, nrows, GDT_Float64, nl, NULL, 0, 0, 0, NULL );
	} else {
		int hasNA=0;
		double na = source[0].gdalconnection->GetRasterBand(1)->GetNoDataValue(&hasNA);
		if (!hasNA) {
			na = NAN;
		}
		if (datatype == "INT4S") {
			//min_max_na(vals, na, (double)INT32_MIN, (double)INT32_MAX);
			//std::vector<int32_t> vv(vals.begin(), vals.end());
			std::vector<int32_t> vv;
			tmp_min_max_na(vv, vals, na, (double)INT32_MIN, (double)INT32_MAX);
			err = source[0].gdalconnection->RasterIO(GF_Write, startcol, startrow, ncols, nrows, &vv[0], ncols, nrows, GDT_Int32, nl, NULL, 0, 0, 0, NULL );
		} else if (datatype == "INT2S") {					
			//min_max_na(vals, na, (double)INT16_MIN, (double)INT16_MAX); 
			//std::vector<int16_t> vv(vals.begin(), vals.end());
			std::vector<int16_t> vv;
			tmp_min_max_na(vv, vals, na, (double)INT16_MIN, (double)INT16_MAX);
			err = source[0].gdalconnection->RasterIO(GF_Write, startcol, startrow, ncols, nrows, &vv[0], ncols, nrows, GDT_Int16, nl, NULL, 0, 0, 0, NULL );
		} else if (datatype == "INT4U") {
			//min_max_na(vals, na, 0, (double)INT32_MAX * 2 - 1);
			//std::vector<uint32_t> vv(vals.begin(), vals.end());
			std::vector<uint32_t> vv;
			tmp_min_max_na(vv, vals, na, 0, (double)INT32_MAX * 2 - 1);
			err = source[0].gdalconnection->RasterIO(GF_Write, startcol, startrow, ncols, nrows, &vv[0], ncols, nrows, GDT_UInt32, nl, NULL, 0, 0, 0, NULL );
		} else if (datatype == "INT2U") {
			//min_max_na(vals, na, 0, (double)INT16_MAX * 2 - 1); 
			//std::vector<uint16_t> vv(vals.begin(), vals.end());
			std::vector<uint16_t> vv;
			tmp_min_max_na(vv, vals, na, 0, (double)INT16_MAX * 2 - 1); 
			err = source[0].gdalconnection->RasterIO(GF_Write, startcol, startrow, ncols, nrows, &vv[0], ncols, nrows, GDT_UInt16, nl, NULL, 0, 0, 0, NULL );
		} else if (datatype == "INT1U") {
			//min_max_na(vals, na, 0, 255);
			//std::vector<int8_t> vv(vals.begin(), vals.end());
			std::vector<int8_t> vv;
			tmp_min_max_na(vv, vals, na, 0, 255);
			err = source[0].gdalconnection->RasterIO(GF_Write, startcol, startrow, ncols, nrows, &vv[0], ncols, nrows, GDT_Byte, nl, NULL, 0, 0, 0, NULL );
		} else {
			setError("bad datatype");
			GDALClose( source[0].gdalconnection );
			return false;
		}
	}

	if (err != CE_None ) {
		setError("cannot write values (err: " + std::to_string(err) +")");
		GDALClose( source[0].gdalconnection );
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
	if (copy_driver != "") {
		GDALDataset *newDS;
		GDALDriver *poDriver;
		poDriver = GetGDALDriverManager()->GetDriverByName(copy_driver.c_str());
		copy_driver = "";		
		if (copy_filename == "") {
			newDS = poDriver->CreateCopy(source[0].filename.c_str(),
				source[0].gdalconnection, FALSE, NULL, NULL, NULL);
			GDALClose( (GDALDatasetH) newDS );
			GDALClose( (GDALDatasetH) source[0].gdalconnection );
		} else {
			GDALClose( (GDALDatasetH) source[0].gdalconnection );
			GDALDataset *oldDS;
			oldDS = (GDALDataset *) GDALOpen(copy_filename.c_str(), GA_ReadOnly );
			if( oldDS == NULL )  {
				setError("something went terribly wrong");
				return false;
			}
			newDS = poDriver->CreateCopy(source[0].filename.c_str(),
				oldDS, FALSE, NULL, NULL, NULL);
			GDALClose( (GDALDatasetH) oldDS );
			GDALClose( (GDALDatasetH) newDS );
			copy_filename = "";
		}
	} else {
		GDALClose( (GDALDatasetH) source[0].gdalconnection );
	}
	source[0].hasValues = true;
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
