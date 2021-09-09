// Copyright (c) 2018-2021  Robert J. Hijmansf
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
#include "gdal_rat.h"

#include "gdalio.h"


bool setCats(GDALRasterBand *poBand, std::vector<std::string> &labels) {
	char **labs = NULL;
	for (size_t i = 0; i < labels.size(); i++) {
		labs = CSLAddString(labs, labels[i].c_str());
	}
	CPLErr err = poBand->SetCategoryNames(labs);
	return (err == CE_None);
}



bool setRat(GDALRasterBand *poBand, SpatDataFrame &d) {

	GDALRasterAttributeTable *pRat = poBand->GetDefaultRAT();
	if (pRat == NULL) {
		return false;
	}
	size_t nr = d.nrow();
	
	for (size_t i=0; i<d.ncol(); i++) {
		const char *fn = d.names[i].c_str();
		if (d.itype[i] == 0) {
			if (pRat->CreateColumn(fn, GFT_Real, GFU_Generic) != CE_None) {
				delete pRat;
				return false;
			};
		} else if (d.itype[i] == 1) {
			if (pRat->CreateColumn(fn, GFT_Integer, GFU_Generic) != CE_None) {
				delete pRat;
				return false;
			}
		} else {
			if (pRat->CreateColumn(fn, GFT_String, GFU_Generic) != CE_None) {
				delete pRat;
				return false;
			}
		}
	}

	pRat->SetRowCount(nr);
	for (size_t i=0; i<d.ncol(); i++) {
		if (d.itype[i] == 0) {
			std::vector<double> v = d.dv[d.iplace[i]];
			if( pRat->ValuesIO(GF_Write, i, 0, nr, &v[0]) != CE_None ) {
				return false;				
			}
		} else if (d.itype[i] == 1) {
			std::vector<long> v = d.iv[d.iplace[i]];
			for (size_t j=0; j<v.size(); j++) {
				pRat->SetValue(j, i, (int)v[j]);
			}
		} else {
			std::vector<std::string> v = d.sv[d.iplace[i]];
			for (size_t j=0; j<v.size(); j++) {
				pRat->SetValue(j, i, v[j].c_str());
			}
		}
	}

	CPLErr err = poBand->SetDefaultRAT(pRat);
	delete pRat;
	return (err == CE_None);
}



bool setCT(GDALRasterBand *poBand, SpatDataFrame &d) {
	CPLErr err = poBand->SetColorInterpretation(GCI_PaletteIndex);
	if (err != CE_None) {
		return false;
	}
	GDALColorTable *poCT = new GDALColorTable(GPI_RGB);
	GDALColorEntry col;
	for (size_t j=0; j< d.nrow(); j++) {
	if (d.iv[3][j] == 0) { // maintain transparaency in gtiff
			col.c1 = 255;
			col.c2 = 255;
			col.c3 = 255;
			col.c4 = 0;
		} else {
			col.c1 = (short)d.iv[0][j];
			col.c2 = (short)d.iv[1][j];
			col.c3 = (short)d.iv[2][j];
			col.c4 = (short)d.iv[3][j];
		}
		poCT->SetColorEntry(j, &col);
	}
	err = poBand->SetColorTable(poCT);
	delete poCT;
	return (err == CE_None);
}


SpatDataFrame grayColorTable() {
	SpatDataFrame coltab;
	std::vector<long> col(256);
	std::iota(col.begin(), col.end(), 0);
	coltab.add_column(col, "red");
	coltab.add_column(col, "green");
	coltab.add_column(col, "blue");
	std::fill(col.begin(), col.end(), 255);
	coltab.add_column(col, "alpha");
	return coltab;
}


bool checkFormatRequirements(const std::string &driver, std::string &filename, std::string &msg) {

	if (driver == "SAGA") {
		std::string ext = getFileExt(filename);
		if (ext != ".sdat") {
			msg = "SAGA filenames must end on '.sdat'";
			return false;
		}
	}

	return true;
}


void stat_options(int sstat, bool &compute_stats, bool &gdal_stats, bool &gdal_minmax, bool &gdal_approx) {
	compute_stats = true;
	gdal_stats  = true;
	gdal_minmax = false;
	if (sstat == 1) {
		gdal_stats = false;
	} else if (sstat == 2) {
		gdal_stats = true;
		gdal_approx = true;
	} else if (sstat == 3) {
		gdal_stats = true;
		gdal_approx = false;
	} else if (sstat == 4) {
		gdal_minmax = true;
		gdal_approx = true;
	} else if (sstat == 5) {
		gdal_minmax = true;
		gdal_approx = false;
	} else {
		compute_stats = false;
	}
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
	std::string datatype = opt.get_datatype();
	
	bool writeRGB = (rgb && nlyr() == 3 && rgblyrs.size() == 3);
	if (writeRGB) {
		datatype = "INT1U";
	}
	std::string errmsg;
	if (!checkFormatRequirements(driver, filename, errmsg)) {
		setError(errmsg);
		return false;
	}
	
	bool append = std::find(opt.gdal_options.begin(), opt.gdal_options.end(), "APPEND_SUBDATASET=YES") != opt.gdal_options.end();
	if (append) {
		if (!file_exists(filename)) {
			setError("cannot append to a file that does not exist");
			return(false);
		} 
	} else if (file_exists(filename) & (!opt.get_overwrite())) {
		setError("file exists. You can use 'overwrite=TRUE' to overwrite it");
		return(false);
	}
	
//	if (!can_write(filename, opt.get_overwrite(), errmsg)) {
//		setError(errmsg);
//		return(false);
//	}
	
	std::vector<bool> hasCT = hasColors();
	std::vector<bool> hasCats = hasCategories();
	std::vector<SpatDataFrame> ct = getColors();
	if (hasCT[0] || hasCats[0]) { 
		datatype = "INT1U";
	} else if (datatype != "INT1U") {
		std::fill(hasCT.begin(), hasCT.end(), false);
	}

	if (opt.datatype_set) {
		if (datatype != opt.get_datatype()) {
			addWarning("changed datatype to " + datatype);
		}
	}

	GDALDataType gdt;
	if (!getGDALDataType(datatype, gdt)) {
		setError("invalid datatype");
		return false;
	}
	source[0].datatype = datatype;

	int dsize = std::stoi(datatype.substr(3,1));
	GIntBig diskNeeded = ncell() * nlyr() * dsize;
	std::string dname = dirname(filename);
	GIntBig diskAvailable = VSIGetDiskFreeSpace(dname.c_str());
	if ((diskAvailable > -1) && (diskAvailable < diskNeeded)) {
		setError("insufficient disk space (perhaps from temporary files?)");
		return(false);		
	}


	stat_options(opt.get_statistics(), compute_stats, gdal_stats, gdal_minmax, gdal_approx);

    GDALDriver *poDriver;
    poDriver = GetGDALDriverManager()->GetDriverByName(driver.c_str());
    if(poDriver == NULL) {
		setError("invalid driver");
		return (false);
	}

	char **papszOptions = set_GDAL_options(driver, diskNeeded, writeRGB, opt.gdal_options);

/*	if (driver == "GTiff") {
		GDAL_tiff_options(diskNeeded > 4194304000, writeRGB, opt);
	}
	
	for (size_t i=0; i<opt.gdal_options.size(); i++) {
		std::vector<std::string> gopt = strsplit(opt.gdal_options[i], "=");
		if (gopt.size() == 2) {
			papszOptions = CSLSetNameValue( papszOptions, gopt[0].c_str(), gopt[1].c_str() );
		}
	}
	if (writeRGB) {
		papszOptions = CSLSetNameValue( papszOptions, "PHOTOMETRIC", "RGB");
	}
*/

    char **papszMetadata;
    papszMetadata = poDriver->GetMetadata();
    if (!CSLFetchBoolean( papszMetadata, GDAL_DCAP_RASTER, FALSE)) {
		setError(driver + " is not a raster format");
		return false;
	}
	//bool isncdf = ((driver == "netCDF" && opt.get_ncdfcopy()));

	GDALDataset *poDS;
	if (CSLFetchBoolean( papszMetadata, GDAL_DCAP_CREATE, FALSE)) {
		poDS = poDriver->Create(filename.c_str(), ncol(), nrow(), nlyr(), gdt, papszOptions);
	} else if (CSLFetchBoolean( papszMetadata, GDAL_DCAP_CREATECOPY, FALSE)) {
		copy_driver = driver;
		gdal_options = opt.gdal_options;
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
		setError("cannot write this format: "+ driver);
		CSLDestroy( papszOptions );
		return false;
	}

	CSLDestroy( papszOptions );
	if (poDS == NULL) {
		if (!filepath_exists(filename)) {
			setError("failed writing "+ driver + " file. Path does not exist");
		} else {
			setError("failed writing "+ driver + " file");
		}
		GDALClose( (GDALDatasetH) poDS );
		return false;	
	}

    #ifdef useRcpp
	if (opt.verbose) {
		double gb = 1073741824;
		char **filelist = poDS->GetFileList();
		std::vector <std::string> files;
		if (filelist != NULL) {
			for (size_t i=0; filelist[i] != NULL; i++) {
				files.push_back(filelist[i]);
				std::replace( files[i].begin(), files[i].end(), '\\', '/'); 
			}
		}
		CSLDestroy( filelist );
		for (size_t i=0; i<files.size(); i++) {
			Rcpp::Rcout<< "filename      : " << files[i] << std::endl;
		}
		Rcpp::Rcout<< "compute stats : " << compute_stats;
		if (compute_stats) {
			Rcpp::Rcout << ", GDAL: "   << gdal_stats << ", minmax: " 
			<< gdal_minmax << ", approx: " << gdal_approx;
		} 
		Rcpp::Rcout << std::endl;

		Rcpp::Rcout<< "driver        : " << driver   << std::endl;
		if (diskAvailable > 0) {
			Rcpp::Rcout<< "disk available: " << roundn(diskAvailable / gb, 1) << " GB" << std::endl;
		}
		Rcpp::Rcout<< "disk needed   : " << roundn(diskNeeded / gb, 1) << " GB" << std::endl;
	}
	#endif


	if (opt.names.size() == nlyr()) {
		setNames(opt.names);
	}
	GDALRasterBand *poBand;
	std::vector<std::string> nms = getNames();
	double naflag=NAN; 
	bool hasNAflag = opt.has_NAflag(naflag);

	if (writeRGB) nms = {"red", "green", "blue"};

	for (size_t i=0; i < nlyr(); i++) {

		poBand = poDS->GetRasterBand(i+1);

		if (hasCT[i]) {
			if (!setCT(poBand, ct[i])) {
				addWarning("could not write the color table");
			}
		}
		if (hasCats[i]) {
			//bool catsdone = false;
			//SpatCategories cats = getLayerCategories(i);
			//if (cats.d.ncol() > 2) {
			//	if (!setRat(poBand, cats.d)) {
			//		catsdone = true;
			//	}
			//}
			//if (!catsdone) {
				std::vector<std::string> labs = getLabels(i);
				if (!setCats(poBand, labs)) {
					addWarning("could not write categories");
				}
			//}
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
				//double na = (double)UINT32_MAX;
				poBand->SetNoDataValue(UINT32_MAX); 
			} else if (datatype == "INT2U") {
				//double na = (double)INT16_MAX * 2 - 1;
				poBand->SetNoDataValue(UINT16_MAX); 
			} else if (datatype == "INT1U") {
				poBand->SetNoDataValue(255); 
			} else {
				poBand->SetNoDataValue(NAN); 
			}
		}

		if (writeRGB) {
			if (rgblyrs[i]==0) {
				poBand->SetColorInterpretation(GCI_RedBand);
			} else if (rgblyrs[i]==1) {
				poBand->SetColorInterpretation(GCI_GreenBand);
			} else if (rgblyrs[i]==2) {
				poBand->SetColorInterpretation(GCI_BlueBand);
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
		source[0].range_min[i] = NAN; //std::numeric_limits<double>::max();
		source[0].range_max[i] = NAN; //std::numeric_limits<double>::lowest();
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


bool SpatRaster::writeValuesGDAL(std::vector<double> &vals, size_t startrow, size_t nrows, size_t startcol, size_t ncols){

	CPLErr err = CE_None;
	double vmin, vmax;
	size_t nc = nrows * ncols;
	size_t nl = nlyr();
	std::string datatype = source[0].datatype;

	if ((compute_stats) && (!gdal_stats)) {
		for (size_t i=0; i < nl; i++) {
			size_t start = nc * i;
			if (datatype == "INT4S") {
				minmaxlim(vals.begin()+start, vals.begin()+start+nc, vmin, vmax, (double)INT32_MIN, (double)INT32_MAX);
			} else if (datatype == "INT2S") {
				minmaxlim(vals.begin()+start, vals.begin()+start+nc, vmin, vmax, (double)INT16_MIN, (double)INT16_MAX);
			} else if (datatype == "INT4U") {
				minmaxlim(vals.begin()+start, vals.begin()+start+nc, vmin, vmax, 0.0, (double)UINT32_MAX);
			} else if (datatype == "INT2U") {
				minmaxlim(vals.begin()+start, vals.begin()+start+nc, vmin, vmax, 0.0, (double)UINT16_MAX);
			} else if (datatype == "INT1U") {
				minmaxlim(vals.begin()+start, vals.begin()+start+nc, vmin, vmax, 0.0, 255.0);
			} else {
				minmax(vals.begin()+start, vals.begin()+start+nc, vmin, vmax);
			}
			if (!std::isnan(vmin)) {
				if (std::isnan(source[0].range_min[i])) {
					source[0].range_min[i] = vmin;
					source[0].range_max[i] = vmax;		
				} else {
					source[0].range_min[i] = std::min(source[0].range_min[i], vmin);
					source[0].range_max[i] = std::max(source[0].range_max[i], vmax);
				}
			}
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
			tmp_min_max_na(vv, vals, na, 0, (double)UINT32_MAX);
			err = source[0].gdalconnection->RasterIO(GF_Write, startcol, startrow, ncols, nrows, &vv[0], ncols, nrows, GDT_UInt32, nl, NULL, 0, 0, 0, NULL );
		} else if (datatype == "INT2U") {
			//min_max_na(vals, na, 0, (double)INT16_MAX * 2 - 1); 
			//std::vector<uint16_t> vv(vals.begin(), vals.end());
			std::vector<uint16_t> vv;
			tmp_min_max_na(vv, vals, na, 0, (double)UINT16_MAX); 
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

		if (compute_stats) {
			if (gdal_stats) {
				double mn, mx, av=-9999, sd=-9999;
				//int approx = gdal_approx;
				if (gdal_minmax) {
					double adfMinMax[2];
					poBand->ComputeRasterMinMax(gdal_approx, adfMinMax);
					mn = adfMinMax[0];
					mx = adfMinMax[1];
				} else {
					poBand->ComputeStatistics(gdal_approx, &mn, &mx, &av, &sd, NULL, NULL);				
				}
				poBand->SetStatistics(mn, mx, av, sd);
			} else {
				if (datatype.substr(0,3) == "INT") {
					source[0].range_min[i] = trunc(source[0].range_min[i]); 
					source[0].range_max[i] = trunc(source[0].range_max[i]); 		
				}
				poBand->SetStatistics(source[0].range_min[i], source[0].range_max[i], -9999., -9999.);
			}
			source[0].hasRange[i] = true;
		} else {
			source[0].hasRange[i] = false;
		}
	}
	
	if (copy_driver != "") {
		GDALDataset *newDS;
		GDALDriver *poDriver;
		char **papszOptions = set_GDAL_options(copy_driver, 0.0, false, gdal_options);
		poDriver = GetGDALDriverManager()->GetDriverByName(copy_driver.c_str());
		if (copy_filename == "") {
			newDS = poDriver->CreateCopy(source[0].filename.c_str(),
				source[0].gdalconnection, FALSE, papszOptions, NULL, NULL);
			if( newDS == NULL )  {
				setError("mem copy create failed for "+ copy_driver);
				copy_driver = "";
				GDALClose( (GDALDatasetH) newDS );
				GDALClose( (GDALDatasetH) source[0].gdalconnection );
				return false;	
			}
			copy_driver = "";
			GDALClose( (GDALDatasetH) newDS );
			GDALClose( (GDALDatasetH) source[0].gdalconnection );
		} else {
			GDALClose( (GDALDatasetH) source[0].gdalconnection );
			GDALDataset *oldDS = openGDAL(copy_filename.c_str(), GDAL_OF_RASTER | GDAL_OF_READONLY, source[0].open_ops);

			if( oldDS == NULL )  {
				setError("file copy create failed for "+ copy_driver);
				copy_driver = "";
				copy_filename = "";
				GDALClose( (GDALDatasetH) oldDS );
				return false;
			}

			newDS = poDriver->CreateCopy(source[0].filename.c_str(),
				oldDS, FALSE, papszOptions, NULL, NULL);
			if( newDS == NULL )  {
				setError("copy create failed for "+ copy_driver);
				copy_driver = "";
				copy_filename = "";
				GDALClose( (GDALDatasetH) oldDS );
				GDALClose( (GDALDatasetH) newDS );
				return false;	
			}
			copy_driver = "";
			copy_filename = "";
			GDALClose( (GDALDatasetH) oldDS );
			GDALClose( (GDALDatasetH) newDS );
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
