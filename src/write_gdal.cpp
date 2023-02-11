// Copyright (c) 2018-2023  Robert J. Hijmans
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
#include "vecmath.h"
#include "recycle.h"

#include <unordered_map>
#include <vector>

#include "gdal_priv.h"
#include "cpl_conv.h" // for CPLMalloc()
#include "cpl_string.h"
#include "ogr_spatialref.h"
#include "gdal_rat.h"

#include "gdalio.h"
/*
void add_quotes(std::vector<std::string> &s) {
	for (size_t i=0; i< s.size(); i++) {
		s[i] = "\"" + s[i] + "\"";
	}
}
*/


std::string quoted_csv(const std::vector<std::string> &s) {
	std::string ss;
	if (s.empty()) {
		ss = "";
		return ss;
	}
	ss = "\"" + s[0] + "\"";
	for (size_t i=1; i< s.size(); i++) {
		ss += ",\"" + s[i] + "\"";
	}
	return ss;
}

bool SpatRaster::write_aux_json(std::string filename) {
	filename += ".aux.json";
	std::ofstream f;
	bool wunits = hasUnit();
	bool wtime = hasTime();
	if (wunits || wtime) {
		f.open(filename);
		if (f.is_open()) {
			f << "{" << std::endl;
			if (wtime) {
				std::vector<std::string> tstr = getTimeStr(false);
				std::string ss = quoted_csv(tstr);
				f << "\"time\":[" << ss << "]," << std::endl;
				f << "\"timestep\":\"" << source[0].timestep << "\"";
				if (wunits) f << ",";
				f << std::endl;
			}
			if (wunits) {
				std::vector<std::string> units = getUnit();
				std::string ss = quoted_csv(units);
				f << "\"unit\":[" << ss << "]" << std::endl;
			}
			f << "}" << std::endl;
		} else {
			f.close();
			return false;
		}
		f.close();
		return true;
	}
	return true;
}




bool setRat(GDALRasterBand *poBand, SpatDataFrame &d) {

	size_t nr = d.nrow();
	if (nr == 0) return true;

//	GDALRasterAttributeTable *pRat = poBand->GetDefaultRAT();
	GDALDefaultRasterAttributeTable *pRat = new GDALDefaultRasterAttributeTable();

	for (size_t i=0; i<d.ncol(); i++) {
		const char *fn = d.names[i].c_str();
		if (d.itype[i] == 0) {
			if (pRat->CreateColumn(fn, GFT_Real, GFU_Generic) != CE_None) {
				return false;
			};
		} else if (d.itype[i] == 1) {
			if (pRat->CreateColumn(fn, GFT_Integer, GFU_Generic) != CE_None) {
				return false;
			}
		} else {
			if (pRat->CreateColumn(fn, GFT_String, GFU_Generic) != CE_None) {
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


bool is_rat(SpatDataFrame &d) {
	if (d.nrow() == 0) return false;
	if (d.ncol() > 2) return true;
	if (d.itype[0] == 1) {
		long dmin = vmin(d.iv[0], true);
		long dmax = vmax(d.iv[0], true);
		if (dmin >= 0 && dmax <= 255) {
			return false;
		}
	} else if (d.itype[0] == 0) {
		double dmin = vmin(d.dv[0], true);
		double dmax = vmax(d.dv[0], true);
		if (dmin >= 0 && dmax <= 255) {
			return false;
		}
	}
	return true;
}

/*
bool setCats(GDALRasterBand *poBand, std::vector<std::string> &labels) {
	char **labs = NULL;
	for (size_t i = 0; i < labels.size(); i++) {
		labs = CSLAddString(labs, labels[i].c_str());
	}
	CPLErr err = poBand->SetCategoryNames(labs);
	return (err == CE_None);
}
*/

bool setBandCategories(GDALRasterBand *poBand, std::vector<long> value, std::vector<std::string> labs) {

	if (labs.size() != value.size()) return false;
	if (vmin(value, false) < 0) return false;
	if (vmax(value, false) > 255) return false;
	std::vector<std::string> s(256, "");

	for (size_t i=0; i<value.size(); i++) {
		s[value[i]] = labs[i];
	}
	char **slabs = NULL;
	for (size_t i = 0; i < s.size(); i++) {
		slabs = CSLAddString(slabs, s[i].c_str());
	}
	CPLErr err = poBand->SetCategoryNames(slabs);
	return (err == CE_None);
}



bool setCT(GDALRasterBand *poBand, SpatDataFrame &d) {

	if (d.ncol() < 5) return false;
	if (d.itype[0] != 1) return false;
	if (d.itype[1] != 1) return false;
	if (d.itype[2] != 1) return false;
	if (d.itype[3] != 1) return false;
	if (d.itype[4] != 1) return false;

	long dmin = vmin(d.iv[0], true);
	long dmax = vmax(d.iv[0], true);
	if (dmin < 0 || dmax > 255) {
		return false;
	}

	SpatDataFrame s;
	s.add_column(1, "red");
	s.add_column(1, "green");
	s.add_column(1, "blue");
	s.add_column(1, "alpha");
	s.resize_rows(256);
	for (size_t i=0; i<d.nrow(); i++) {
		s.iv[0][d.iv[0][i]] = d.iv[1][i];
		s.iv[1][d.iv[0][i]] = d.iv[2][i];
		s.iv[2][d.iv[0][i]] = d.iv[3][i];
		s.iv[3][d.iv[0][i]] = d.iv[4][i];
	}

	CPLErr err = poBand->SetColorInterpretation(GCI_PaletteIndex);
	if (err != CE_None) {
		return false;
	}
	GDALColorTable *poCT = new GDALColorTable(GPI_RGB);
	GDALColorEntry col;
	for (size_t j=0; j< s.nrow(); j++) {
		if (s.iv[3][j] == 0) { // maintain transparency in gtiff
			col.c1 = 255;
			col.c2 = 255;
			col.c3 = 255;
			col.c4 = 0;
		} else {
			col.c1 = (short)s.iv[0][j];
			col.c2 = (short)s.iv[1][j];
			col.c3 = (short)s.iv[2][j];
			col.c4 = (short)s.iv[3][j];
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


void removeVatJson(std::string filename) {
	std::vector<std::string> exts = {".vat.dbf", ".vat.cpg", ".json"};
	for (size_t i=0; i<exts.size(); i++) {
		std::string f = filename + exts[i];
		if (file_exists(f)) {
			remove(f.c_str());
		}
	}
}


bool SpatRaster::writeStartGDAL(SpatOptions &opt, const std::vector<std::string> &srcnames) {

	std::string filename = opt.get_filename();
	if (filename.empty()) {
		setError("empty filename");
		return(false);
	} else {
		// make sure filename won't be used again
		opt.set_filenames({""});
	}

	std::string driver = opt.get_filetype();
	getGDALdriver(filename, driver);
	if (driver.empty()) {
		setError("cannot guess file type from filename");
		return(false);
	}
    GDALDriver *poDriver;
    poDriver = GetGDALDriverManager()->GetDriverByName(driver.c_str());
    if(poDriver == NULL) {
		setError("invalid driver");
		return (false);
	}
    char **papszMetadata;
    papszMetadata = poDriver->GetMetadata();
    if (!CSLFetchBoolean( papszMetadata, GDAL_DCAP_RASTER, FALSE)) {
		setError(driver + " is not a raster format");
		return false;
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

	std::string appstr = "APPEND_SUBDATASET=YES";
	bool append = std::find(opt.gdal_options.begin(), opt.gdal_options.end(), appstr) != opt.gdal_options.end();
	if (append && (!CSLFetchBoolean( papszMetadata, GDAL_DMD_SUBDATASETS, FALSE))) {
		setError("cannot append datasets with this file format");
		return false;
	}
	if (append && opt.get_overwrite()) {
		setError("cannot append and overwrite at the same time");
		return false;
	}
	if (!append) {
		std::string msg;
		if (!can_write({filename}, srcnames, opt.get_overwrite(), msg)) {
			setError(msg);
			return false;
		}
	}
	removeVatJson(filename);

// what if append=true?
	std::string auxf = filename + ".aux.xml";
	remove(auxf.c_str());
	auxf = filename + ".aux.json";
	remove(auxf.c_str());

	std::vector<bool> hasCT = hasColors();
	std::vector<bool> hasCats = hasCategories();
	std::vector<SpatDataFrame> ct = getColors();
	bool cat = hasCats[0];
	bool warnCT = true;

	bool rat = cat ? is_rat(source[0].cats[0].d) : false;
	if (rat) {
		// needs redesign. Is CT also part of RAT? 
		// other layers affected? etc.
		warnCT = false;
		if (hasCT[0]) {
			if (ct[0].nrow() < 256) {
				if (opt.datatype_set && (datatype != opt.get_datatype())) {
					addWarning("change datatype to INT1U to write the color-table");
				} else {
					datatype = "INT1U";
				}
			}
		} else {
			//if (opt.datatype_set) {
			//	std::string sdt = opt.get_datatype().substr(0, 3); 
			//	if (sdt != "INT") {
			//		addWarning("change datatype to an INT type to write the categories");
			//	}
			//} else {
			//	datatype = "INT4S";
			//}
			if (!opt.datatype_set && (driver != "GPKG")) {
				datatype = "INT4S";
			}
		}
	} else if (hasCT[0] || cat) {
		if (opt.datatype_set && (datatype != "INT1U")) {
			addWarning("change datatype to INT1U to write the color-table");
		} else {
			datatype = "INT1U";
		}
	} else if (datatype != "INT1U") {
		std::fill(hasCT.begin(), hasCT.end(), false);
	}
	//if (opt.datatype_set) {
	//	if (datatype != opt.get_datatype()) {
	//		addWarning("changed datatype to " + datatype);
	//	}
	//}

	GDALDataType gdt;
	if (!getGDALDataType(datatype, gdt)) {
		setError("invalid datatype");
		return false;
	}

	int dsize = std::stoi(datatype.substr(3,1));
	GIntBig diskNeeded = ncell() * nlyr() * dsize;
	std::string dname = dirname(filename);
	GIntBig diskAvailable = VSIGetDiskFreeSpace(dname.c_str());
	if ((diskAvailable > -1) && (diskAvailable < diskNeeded)) {
		long gb = 1073741824;
		std::string msg = "Estimated disk space needed without compression: " + std::to_string(diskNeeded/gb) + "GB. Available: " + std::to_string(diskAvailable/gb) + " GB.";
		// was an error, but actual file size is not known
		addWarning(msg);
	}

	stat_options(opt.get_statistics(), compute_stats, gdal_stats, gdal_minmax, gdal_approx);
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
			//std::string driver = opt.get_filetype();
			//std::string f = tempFile(opt.get_tempdir(), opt.pid, "");
			//getGDALdriver(f, driver);
			//if (driver == "") {
			//	setError("invalid default temp filetype");
			//	return(false);
			//}

			std::string f, driver;
			if (!getTempFile(f, driver, opt)) {
				return false;
			}
			//std::string f = tempFile(opt.get_tempdir(), opt.pid, ".tif");
			copy_filename = f;

			GDALDriver *poDriver;
			//poDriver = GetGDALDriverManager()->GetDriverByName("GTiff");
			poDriver = GetGDALDriverManager()->GetDriverByName(driver.c_str());
			if(poDriver == NULL) {
				setError("invalid driver");
				return false;
			}
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
			setError("failed writing "+ driver + " file. Path does not exist:\n   " + filename);
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

	std::vector<double> scale = opt.get_scale();
	std::vector<double> offset = opt.get_offset();
	size_t nl = nlyr();
	if (((scale.size() > 1) || (offset.size())) || 
		((scale[0] != 1) || (offset[0] != 0))) {	
		recycle(scale, nl);
		recycle(offset, nl);
	}
	bool scoff = false;
	for (size_t i=0; i<scale.size(); i++) {
		//if (scale[i] == 0) scale[i] = 1;
		if ((scale[i] != 1) || (offset[i] != 0)) {
			if (!scoff) {
				source[0].has_scale_offset = std::vector<bool>(nl, false);
				scoff = true;
			}
			source[0].has_scale_offset[i] = true;
		}
	}
	if (scoff) {
		source[0].scale  = scale;
		source[0].offset = offset;
	}

	bool scoffwarning = false;
	
	for (size_t i=0; i < nlyr(); i++) {

		poBand = poDS->GetRasterBand(i+1);
		if ((i==0) && hasCT[i]) {
			if (!setCT(poBand, ct[i])) {
				if (warnCT) {
					addWarning("could not write the color table");
				}
			}
		}
		if (hasCats[i]) {
			if (is_rat(source[0].cats[i].d)) {
				if (!setRat(poBand, source[0].cats[i].d)) {
					addWarning("could not write attribute table");
				}
			} else {
				SpatCategories lyrcats = getLayerCategories(i);
				if (lyrcats.d.ncol() == 2) {
					std::vector<std::string> labs = getLabels(i);
					std::vector<long> ind = lyrcats.d.as_long(0);
					if (!setBandCategories(poBand, ind, labs)) {
						addWarning("could not write categories");
					}
				}
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
				poBand->SetNoDataValue(UINT32_MAX);
			} else if (datatype == "INT2U") {
				//double na = (double)INT16_MAX * 2 - 1;
				poBand->SetNoDataValue(UINT16_MAX);
			} else if (datatype == "INT1U") {
				poBand->SetNoDataValue(255);
			} else if (datatype == "INT1S") {
				poBand->SetNoDataValue(-128); //GDT_Int8
#if GDAL_VERSION_MAJOR <= 3 && GDAL_VERSION_MINOR < 5
// no Int64
#else 
			} else if (datatype == "INT8S") {
				//INT64_MIN == -9223372036854775808;
				if (poBand->SetNoDataValueAsInt64(INT64_MIN) != CE_None) {
					addWarning("no data problem");
				}			
			} else if (datatype == "INT8U") {
				if (poBand->SetNoDataValueAsUInt64(UINT64_MAX-1101) != CE_None) {
					addWarning("no data problem");
				}			
#endif
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
		
		if (scoff) {
			if (source[0].has_scale_offset[i]) {
				bool failed = (poBand->SetScale(scale[i])) != CE_None;
				if (!failed) {
					failed = ((poBand->SetOffset(offset[i])) != CE_None);
				}
				if (failed) {
					source[0].has_scale_offset[i] = false;
					source[0].scale[i]  = 1;
					source[0].offset[i] = 0;
					scoffwarning = true;
				}
			}
		}
	}

	if (scoffwarning) {
		addWarning("could not set offset");
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

#if GDAL_VERSION_MAJOR >= 3
	const char *options[3] = { "MULTILINE=YES", "FORMAT=WKT2", NULL };
	oSRS.exportToWkt(&pszSRS_WKT, options);
#else
	oSRS.exportToWkt(&pszSRS_WKT);
#endif

	poDS->SetProjection(pszSRS_WKT);
	CPLFree(pszSRS_WKT);
	// destroySRS(oSRS) ?

	source[0].resize(nlyr());
	source[0].nlyrfile = nlyr();
	source[0].dtype = datatype;
	for (size_t i =0; i<nlyr(); i++) {
		source[0].range_min[i] = NAN; //std::numeric_limits<double>::max();
		source[0].range_max[i] = NAN; //std::numeric_limits<double>::lowest();
	}
	source[0].driver = "gdal" ;
	source[0].filename = filename;
	source[0].memory = false;
	write_aux_json(filename);

/*
	if (append) {
		GDALClose( (GDALDatasetH) poDS );
		std::vector<std::string> ops;
		poDS = openGDAL(filename, GDAL_OF_RASTER | GDAL_OF_UPDATE, ops);
		std::vector<std::string> subds;
		char **metadata = poDS->GetMetadata("SUBDATASETS");
		if (metadata != NULL) {
			for (size_t i=0; metadata[i] != NULL; i++) {
				subds.push_back(metadata[i]);
			}
			std::vector<std::vector<std::string>> s = parse_metadata_sds(subds);
			GDALClose( (GDALDatasetH) poDS );
			filename = s[0].back();
			poDS = openGDAL(filename, GDAL_OF_RASTER | GDAL_OF_UPDATE, ops);
		}
	}
*/

	source[0].gdalconnection = poDS;
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


template <typename Iterator>
void minmaxlim(Iterator start, Iterator end, double &vmin, double &vmax, const double &lmin, const double &lmax, bool& outrange) {
    vmin = std::numeric_limits<double>::max();
    vmax = std::numeric_limits<double>::lowest();
    bool none = true;
	for (Iterator v = start; v !=end; ++v) {
		if (!std::isnan(*v)) {
			if (*v >= lmin && *v <= lmax) {
				if (*v > vmax) {
					vmax = *v;
					none = false;
				}
				if (*v < vmin) {
					vmin = *v;
				}
			} else {
				outrange = true;
			}
		}
    }
    if (none) {
        vmin = NAN;
        vmax = NAN;
    }
	vmin = std::trunc(vmin);	
	vmax = std::trunc(vmax);	
}



bool SpatRaster::writeValuesGDAL(std::vector<double> &vals, size_t startrow, size_t nrows, size_t startcol, size_t ncols){

	CPLErr err = CE_None;
	double vmin, vmax;
	size_t nc = nrows * ncols;
	size_t nl = nlyr();
	std::string datatype = source[0].dtype;

	size_t n = vals.size() / nl;
	for (size_t i=0; i<nl; i++) {
		if (source[0].has_scale_offset[i]) {
			size_t start = i*n;
			for (size_t j=start; j<(start+n); j++) {
				vals[j] = (vals[j] - source[0].offset[i]) / source[0].scale[i];
			}
		}
	}

	if ((compute_stats) && (!gdal_stats)) {
		bool invalid = false;
		for (size_t i=0; i < nl; i++) {
			size_t start = nc * i;
			if (datatype == "INT8S") {
				minmaxlim(vals.begin()+start, vals.begin()+start+nc, vmin, vmax, (double)INT64_MIN, (double)INT64_MAX, invalid);
			} else if (datatype == "INT4S") {
				minmaxlim(vals.begin()+start, vals.begin()+start+nc, vmin, vmax, (double)INT32_MIN, (double)INT32_MAX, invalid);
			} else if (datatype == "INT2S") {
				minmaxlim(vals.begin()+start, vals.begin()+start+nc, vmin, vmax, (double)INT16_MIN, (double)INT16_MAX, invalid);
			} else if (datatype == "INT8U") {
				minmaxlim(vals.begin()+start, vals.begin()+start+nc, vmin, vmax, 0.0, (double)UINT64_MAX, invalid);
			} else if (datatype == "INT4U") {
				minmaxlim(vals.begin()+start, vals.begin()+start+nc, vmin, vmax, 0.0, (double)UINT32_MAX, invalid);
			} else if (datatype == "INT2U") {
				minmaxlim(vals.begin()+start, vals.begin()+start+nc, vmin, vmax, 0.0, (double)UINT16_MAX, invalid);
			} else if (datatype == "INT1U") {
				minmaxlim(vals.begin()+start, vals.begin()+start+nc, vmin, vmax, 0.0, 255.0, invalid);
			} else if (datatype == "INT1S") {
				minmaxlim(vals.begin()+start, vals.begin()+start+nc, vmin, vmax, -128.0, 127.0, invalid);
			} else {
				minmax(vals.begin()+start, vals.begin()+start+nc, vmin, vmax);
			}
			if (source[0].has_scale_offset[i]) {
				vmin = vmin * source[0].scale[i] + source[0].offset[i];
				vmax = vmax * source[0].scale[i] + source[0].offset[i];
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
		if (invalid) {
			addWarning("detected values outside of the limits of datatype " + datatype);
		}
	}

	int hasNA = 0;
	double na = source[0].gdalconnection->GetRasterBand(1)->GetNoDataValue(&hasNA);
	if ((datatype == "FLT8S") || (datatype == "FLT4S")) {
		if (hasNA) {
			size_t n = vals.size();
			for (size_t i=0; i<n; i++) {
				if (std::isnan(vals[i])) vals[i] = na;
			}
		}
		err = source[0].gdalconnection->RasterIO(GF_Write, startcol, startrow, ncols, nrows, &vals[0], ncols, nrows, GDT_Float64, nl, NULL, 0, 0, 0, NULL );
	} else {
		if (datatype == "INT8S") {
#if GDAL_VERSION_MAJOR <= 3 && GDAL_VERSION_MINOR < 5
			setError("cannot write INT8S values with GDAL < 3.5");
			GDALClose( source[0].gdalconnection );
			return false;	
#else 			
			std::vector<int64_t> vv;
			tmp_min_max_na(vv, vals, na, (double)INT64_MIN, (double)INT64_MAX);
			err = source[0].gdalconnection->RasterIO(GF_Write, startcol, startrow, ncols, nrows, &vv[0], ncols, nrows, GDT_Int64, nl, NULL, 0, 0, 0, NULL );
#endif
		} else if (datatype == "INT4S") {
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
		} else if (datatype == "INT1S") {
#if GDAL_VERSION_MAJOR <= 3 && GDAL_VERSION_MINOR < 7
			setError("cannot write INT1S values with GDAL < 3.7");
			GDALClose( source[0].gdalconnection );
			return false;	
#else 			
			std::vector<int8_t> vv;
			tmp_min_max_na(vv, vals, na, -127.0, 128.0);
			err = source[0].gdalconnection->RasterIO(GF_Write, startcol, startrow, ncols, nrows, &vv[0], ncols, nrows, GDT_Int8, nl, NULL, 0, 0, 0, NULL );
#endif

		} else if (datatype == "INT8U") {
#if GDAL_VERSION_MAJOR <= 3 && GDAL_VERSION_MINOR < 5
			setError("cannot write INT8U values with GDAL < 3.5");
			GDALClose( source[0].gdalconnection );
			return false;	
#else 			
			std::vector<uint64_t> vv;
			tmp_min_max_na(vv, vals, na, 0, (double)UINT64_MAX);
			err = source[0].gdalconnection->RasterIO(GF_Write, startcol, startrow, ncols, nrows, &vv[0], ncols, nrows, GDT_UInt64, nl, NULL, 0, 0, 0, NULL );
#endif			
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
			std::vector<uint8_t> vv;
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
	std::string datatype = source[0].dtype;

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
				} else if (datatype == "FLT4S") { // match precision
					source[0].range_min[i] = (float) source[0].range_min[i]; 
					source[0].range_max[i] = (float) source[0].range_max[i]; 
				}
				poBand->SetStatistics(source[0].range_min[i], source[0].range_max[i], -9999., -9999.);
			}
			source[0].hasRange[i] = true;
		} else {
			source[0].hasRange[i] = false;
		}
	}

	if (copy_driver.empty()) {
		GDALClose( (GDALDatasetH) source[0].gdalconnection );
	} else {
		GDALDataset *newDS;
		GDALDriver *poDriver;
		char **papszOptions = set_GDAL_options(copy_driver, 0.0, false, gdal_options);
		poDriver = GetGDALDriverManager()->GetDriverByName(copy_driver.c_str());
		if (copy_filename.empty()) {
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

			GDALDataset *oldDS = openGDAL(copy_filename.c_str(), GDAL_OF_RASTER | GDAL_OF_READONLY, source[0].open_drivers, source[0].open_ops);

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
		CSLDestroy(papszOptions);
	}
	source[0].hasValues = true;

	return true;
}



bool SpatRaster::fillValuesGDAL(double fillvalue) {
	CPLErr err = CE_None;
	GDALRasterBand *poBand;
	int hasNA;
	for (size_t i=0; i < nlyr(); i++) {
		poBand = source[0].gdalconnection->GetRasterBand(i+1);
		if (std::isnan(fillvalue)) {
			double naflag = poBand->GetNoDataValue(&hasNA);
			if (hasNA) {
				err = poBand->Fill(naflag);
			} else {
				err = poBand->Fill(fillvalue);
			}
		} else {
			err = poBand->Fill(fillvalue);
		}
	}
	if (err != CE_None ) {
		setError("cannot fill values");
		return false;
	}
	return true;
}


bool SpatRaster::update_meta(bool names, bool crs, bool ext, SpatOptions &opt) { 
	if ((!names) & (!crs) & (!ext)) {
		addWarning("nothing to do");
		return false;
	}
	GDALDatasetH hDS;
	GDALRasterBandH poBand;
	size_t n=0;
	for (size_t i=0; i<nsrc(); i++) {
		if (source[i].memory) continue;
		n++;
		if (!open_gdal(hDS, i, true, opt)) {
			setError("cannot open source " + std::to_string(i+1));
			return false;
		}
		if (names) {
			for (size_t b=0; b < source[i].nlyr; b++) {
				poBand = GDALGetRasterBand(hDS, b+1);
				GDALSetDescription(poBand, source[i].names[b].c_str());
			}
		} 
		if (crs) {
			std::string crs = source[i].srs.wkt;
			OGRSpatialReference oSRS;
			OGRErr erro = oSRS.SetFromUserInput(&crs[0]);
			if (erro == 4) {
				setError("CRS failure");
				GDALClose( hDS );
				return false ;
			}
			char *pszSRS_WKT = NULL;
		#if GDAL_VERSION_MAJOR >= 3
			const char *options[3] = { "MULTILINE=YES", "FORMAT=WKT2", NULL };
			oSRS.exportToWkt(&pszSRS_WKT, options);
		#else
			oSRS.exportToWkt(&pszSRS_WKT);
		#endif
			GDALSetProjection(hDS, pszSRS_WKT);
			CPLFree(pszSRS_WKT);
		}
		if (ext) {
			std::vector<double> rs = resolution();
			SpatExtent extent = getExtent();
			double adfGeoTransform[6] = { extent.xmin, rs[0], 0, extent.ymax, 0, -1 * rs[1] };
			GDALSetGeoTransform(hDS, adfGeoTransform);
		}
		GDALClose(hDS);
	}
	if (n == 0) {
		addWarning("no sources on disk");
		return false;
	}
	return true;
}
