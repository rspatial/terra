#include <unordered_map>

#include "ogr_spatialref.h"

#include "spatRaster.h"
#include "string_utils.h"
#include "file_utils.h"
#include "crs.h"
//#include <vector>
//#include <string>

#include "cpl_port.h"
#include "cpl_conv.h" // CPLFree()
//#include "gdal_version.h"




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
		{".nc","netCDF"}, {".cdf","netCDF"},
		{".gpkg","GPKG"},
		{".img","HFA"}, {".ige","HFA"},
		{".bmp","BMP"},
		{".flt","EHdr"},
		{".grd","RRASTER"}, {".gri","RRASTER"},
		{".sgrd","SAGA"}, {".sdat","SAGA"},
		{".rst","RST"}, {".rdc","RST"},
		{".envi","ENVI"},
		{".asc","AAIGrid"},
		{".bmp","BMP"},
//		{".jpg","JPEG"}, or JPEG2000?
		{".png","PNG"},
		{".gif","GIF"},
	};

    auto i = drivers.find(ext);
    if (i != drivers.end()) {
		driver = i->second;
	}
}

bool SpatRaster::getTempFile(std::string &filename, std::string &driver, SpatOptions& opt) {

	driver = opt.get_def_filetype();
	if ((driver == "") || (driver == "GTiff")) {
		driver = "GTiff";
		filename = tempFile(opt.get_tempdir(), opt.pid, ".tif");
		return true;
	}
	filename = tempFile(opt.get_tempdir(), opt.pid, "");
	std::unordered_map<std::string, std::string>
	exts = {
		{"GTiff", ".tif"},
		{"NetCDF", ".nc"},
		{"GPKG", ".gpkg"},
		{"HFA", ".img"},
		{"RRASTER", ".grd"},
		{"SAGA", ".sgrd"},
		{"RST", ".rst"},
		{"ENVI", ".envi"},
		{"AAIGrid", ".asc"},
	};
    auto i = exts.find(driver);
    if (i != exts.end()) {
		filename += i->second;
	}
	return true;
}


/*
std::string sectostr(int x) {
	char buffer[20];
	time_t now = x;
	tm *utc = gmtime(&now);
	strftime (buffer, 20, "%Y-%m-%d %H:%M:%S", utc);
	std::string s = buffer;
	return s;
}
*/

GDALDataset* openGDAL(std::string filename, unsigned OpenFlag, std::vector<std::string> allowed_drivers, std::vector<std::string> open_options) {

	char ** openops = NULL;
	// for ncdf
//	openops = CSLSetNameValue(openops, "@HONOUR_VALID_RANGE", "NO");

	for (size_t i=0; i<open_options.size(); i++) {
		std::vector<std::string> opt = strsplit(open_options[i], "=");
		if (opt.size() == 2) {
//			Rcpp::Rcout << opt[0] << "=" << opt[1] << std::endl;
			openops = CSLSetNameValue(openops, opt[0].c_str(), opt[1].c_str());
		}
	}
	char ** drivers = NULL;
	for (size_t i=0; i<allowed_drivers.size(); i++) {
		drivers = CSLAddString(drivers, allowed_drivers[i].c_str());
	}

	GDALDataset *poDataset = static_cast<GDALDataset*>(GDALOpenEx( filename.c_str(), OpenFlag, drivers, openops, NULL));
	CSLDestroy(openops);
	CSLDestroy(drivers);
	return poDataset;
}


std::vector<std::string> get_metadata(std::string filename) {

	std::vector<std::string> out;

	std::vector<std::string> ops;
    GDALDataset *poDataset = openGDAL(filename, GDAL_OF_RASTER | GDAL_OF_READONLY, ops, ops);
    if( poDataset == NULL )  {
		return out;
	}

	char **m = poDataset->GetMetadata();
	if (m != NULL) { // needed
		while (*m != nullptr) {
			out.push_back(*m++);
		}
	}
	GDALClose( (GDALDatasetH) poDataset );
	return out;
}


std::vector<std::string> get_metadata_sds(std::string filename) {
	std::vector<std::string> meta;
    GDALDataset *poDataset = openGDAL(filename, GDAL_OF_RASTER | GDAL_OF_READONLY, meta, meta);
    if( poDataset == NULL )  {
		return meta;
	}
	char **metadata = poDataset->GetMetadata("SUBDATASETS");
	if (metadata != NULL) {
		for (size_t i=0; metadata[i] != NULL; i++) {
			meta.push_back(metadata[i]);
		}
	}
	GDALClose( (GDALDatasetH) poDataset );
	return meta;
}

std::vector<std::vector<std::string>> parse_metadata_sds(std::vector<std::string> meta) {

	std::vector<std::string> name, var, desc, nr, nc, nl;
	std::string ndelim = "NAME=";
	std::string ddelim = "DESC=";

    for (size_t i=0; i<meta.size(); i++) {
		std::string s = meta[i];
		size_t pos = s.find(ndelim);
		if (pos != std::string::npos) {
			s.erase(0, pos + ndelim.length());
			name.push_back(s);
			std::string vdelim = ":";
			size_t pos = s.find_last_of(vdelim);
			if (pos != std::string::npos) {
				s.erase(0, pos + vdelim.length());
				var.push_back(s);
			} else {
				std::string v = "v" + std::to_string(i);
				var.push_back(v);
			}
		} else {
			size_t pos = s.find(ddelim);
			if (pos != std::string::npos) {
				s.erase(0, pos + ddelim.length());
				pos = s.find("]");
				std::string dims = s.substr(1, pos-1);

				std::vector<std::string> d = strsplit(dims, "x");
				if (d.size() < 2) {
					nl.push_back("0");
					nr.push_back("0");
					nc.push_back("0");
				} else if (d.size() == 2) {
					nl.push_back("1");
					nr.push_back(d[0]);
					nc.push_back(d[1]);
				} else {
					size_t ds = d.size()-1;
					size_t nls = 0;
					try {
						nls = stoi(d[ds-2]);
						for (size_t i=0; i<(ds-2); i++) {
							nls *= stoi(d[i]);
						}
					} catch(...) {}
					nl.push_back(std::to_string(nls));
					nr.push_back(d[ds-1]);
					nc.push_back(d[ds]);
				}
				//desc.push_back(std::string(pos, s.size()));
				s = s.substr(pos+2, s.size());
				pos = s.find(" ");
				s = s.substr(0, pos);
				desc.push_back(s);


			//	nr.push_back( std::to_string(sub.nrow()));
			//	nc.push_back(std::to_string(sub.ncol()));
			//	nl.push_back(std::to_string(sub.nlyr()));

			} else {
				desc.push_back("");
			}
		}
	}
	std::vector<std::vector<std::string>> out(6);
	out[0] = name;
	out[1] = var;
	out[2] = desc;
	out[3] = nr;
	out[4] = nc;
	out[5] = nl;
	return out;
}



std::vector<std::vector<std::string>> sdinfo(std::string fname) {

	std::vector<std::vector<std::string>> out(6);
	std::vector<std::string> ops;
    GDALDataset *poDataset = openGDAL(fname, GDAL_OF_RASTER | GDAL_OF_READONLY, ops, ops);

    if( poDataset == NULL ) {
		if (!file_exists(fname)) {
			out[0] = std::vector<std::string> {"no such file"};
		} else {
			out[0] = std::vector<std::string> {"cannot open file"};
		}
		return out;
	}
	char **metadata = poDataset->GetMetadata("SUBDATASETS");
	if (metadata == NULL) {
		out[0] = std::vector<std::string> {"no subdatasets"};
		GDALClose( (GDALDatasetH) poDataset );
		return out;
	}
	std::vector<std::string> meta;
	for (size_t i=0; metadata[i] != NULL; i++) {
		meta.push_back(metadata[i]);
	}
	if (meta.size() == 0) {
		GDALClose( (GDALDatasetH) poDataset );
		out[0] = std::vector<std::string> {"no subdatasets"};
		return out;
	}

	SpatRaster sub;
	std::vector<std::string> name, var, desc, nr, nc, nl;
	std::string ndelim = "NAME=";
	std::string ddelim = "DESC=";

    for (size_t i=0; i<meta.size(); i++) {
		std::string s = meta[i];
		//Rcpp::Rcout << s << std::endl;
		size_t pos = s.find(ndelim);
		if (pos != std::string::npos) {
			s.erase(0, pos + ndelim.length());
			name.push_back(s);
			std::string vdelim = ":";
			size_t pos = s.find_last_of(vdelim);
			if (sub.constructFromFile(s, {-1}, {""}, {}, {})) {
				nr.push_back( std::to_string(sub.nrow()));
				nc.push_back(std::to_string(sub.ncol()));
				nl.push_back(std::to_string(sub.nlyr()));
			}
			if (pos != std::string::npos) {
				s.erase(0, pos + vdelim.length());
				var.push_back(s);
			} else {
				var.push_back("v" + std::to_string(i));
			}
		} else {
			size_t pos = s.find(ddelim);
			if (pos != std::string::npos) {
				s.erase(0, pos + ddelim.length());
				desc.push_back(s);
			} else {
				desc.push_back("");
			}
		}
	}
	GDALClose( (GDALDatasetH) poDataset );
	//var.resize(name.size());
	out[0] = name;
	out[1] = var;
	out[2] = desc;
	out[3] = nr;
	out[4] = nc;
	out[5] = nl;
	return out;
}




#if GDAL_VERSION_MAJOR <= 2 && GDAL_VERSION_MINOR < 1

SpatRaster SpatRaster::make_vrt(std::vector<std::string> filenames, std::vector<std::string> options, SpatOptions &opt) {
	SpatRaster out;
	out.setError( "GDAL version >= 2.1 required for vrt");
	return out;
}


std::string gdalinfo(std::string filename, std::vector<std::string> options, std::vector<std::string> oo) {
	std::string out = "GDAL version >= 2.1 required for gdalinfo";
	return out;
}


#else

# include "gdal_utils.h" // requires >= 2.1

SpatRaster SpatRaster::make_vrt(std::vector<std::string> filenames, std::vector<std::string> options, SpatOptions &opt) {

	SpatRaster out;
	std::string outfile = opt.get_filename();
	if (outfile == "") {
		outfile = tempFile(opt.get_tempdir(), opt.pid, ".vrt");
	} else if (file_exists(outfile) && (!opt.get_overwrite())) {
		out.setError("output file exists. You can use 'overwrite=TRUE' to overwrite it");
		return(out);
	}

/*
	std::vector<GDALDataset *> tiles;
	std::vector<std::string> ops;

	for (std::string& f : filenames) {
		GDALDataset *poDataset = openGDAL(f, GDAL_OF_RASTER | GDAL_OF_READONLY, ops, ops);
		if( poDataset == NULL )  {
			for (size_t j= 0; j<tiles.size(); j++) GDALClose(tiles[j]);
			out.setError("cannot open " + f);
			return out;
		}
		tiles.push_back(poDataset);
	}
*/

	char **names = NULL;
	for (std::string& f : filenames) {
		names = CSLAddString(names, f.c_str());
	}

//	psOptions * vrtops;
/*gdalbuildvrt [-tileindex field_name]
            [-resolution {highest|lowest|average|user}]
            [-te xmin ymin xmax ymax] [-tr xres yres] [-tap]
            [-separate] [-b band]* [-sd subdataset]
            [-allow_projection_difference] [-q]
            [-optim {[AUTO]/VECTOR/RASTER}]
            [-addalpha] [-hidenodata]
            [-srcnodata "value [value...]"] [-vrtnodata "value [value...]"]
            [-ignore_srcmaskband]
            [-a_srs srs_def]
            [-r {nearest,bilinear,cubic,cubicspline,lanczos,average,mode}]
            [-oo NAME=VALUE]*
*/

	std::vector <char *> vops = string_to_charpnt(options);
	GDALBuildVRTOptions* vrtops = GDALBuildVRTOptionsNew(vops.data(), NULL);
	if (vrtops == NULL) {
		out.setError("options error");
		CSLDestroy( names );
		return(out);
	}
	int pbUsageError;
//	GDALDataset *ds = (GDALDataset *) GDALBuildVRT(outfile.c_str(), tiles.size(), (GDALDatasetH *) tiles.data(), nullptr, vrtops, &pbUsageError);

	GDALDataset *ds = (GDALDataset *) GDALBuildVRT(outfile.c_str(), filenames.size(), nullptr, names, vrtops, &pbUsageError);
	
	GDALBuildVRTOptionsFree(vrtops);
	CSLDestroy( names );

	//for (size_t i= 0; i<tiles.size(); i++) GDALClose(tiles[i]);
	if(ds == NULL )  {
		out.setError("cannot create vrt. Error #"+ std::to_string(pbUsageError));
		return out;
	}
	GDALClose(ds);
	if (!out.constructFromFile(outfile, {-1}, {""}, {}, {})) {
		out.setError("cannot open created vrt");
		return out;
	}
	return out;
}



std::string gdalinfo(std::string filename, std::vector<std::string> options, std::vector<std::string> openopts) {
// adapted from the 'sf' package by Edzer Pebesma et al

	std::string out = "";
	char ** opops = NULL;
	for (size_t i=0; i<openopts.size(); i++) {
		std::vector<std::string> opt = strsplit(openopts[i], "=");
		if (opt.size() == 2) {
			opops = CSLSetNameValue(opops, opt[0].c_str(), opt[1].c_str());
		}
	}
	GDALDatasetH ds = GDALOpenEx(filename.c_str(), GA_ReadOnly, NULL, opops, NULL);
	//if (opops != NULL) CSLDestroy(opops);
	if (ds == NULL) return out;

	std::vector <char *> options_char = string_to_charpnt(options);
	GDALInfoOptions* opt = GDALInfoOptionsNew(options_char.data(), NULL);
	char *val = GDALInfo(ds, opt);
	out = val;
	CPLFree(val);
	GDALClose(ds);
	GDALInfoOptionsFree(opt);
	return out;
}


#endif

bool getNAvalue(GDALDataType gdt, double & naval) {
	if (gdt == GDT_Float32) {
		naval = NAN;
	} else if (gdt == GDT_Int32) {
		naval = INT32_MIN;
	} else if (gdt == GDT_Float64) {
		naval = NAN;
	} else if (gdt == GDT_Int16) {
		naval = INT16_MIN;
	} else if (gdt == GDT_UInt32) {
		naval = UINT32_MAX;
	} else if (gdt == GDT_UInt16) {
		naval = UINT16_MAX;
	} else if (gdt == GDT_Byte) {
		naval = 255;
	} else {
		naval = NAN;
		return false;
	}
	return true;
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





bool GDALsetSRS(GDALDatasetH &hDS, const std::string &crs) {
	OGRSpatialReferenceH hSRS = OSRNewSpatialReference( NULL );
	OGRErr erro = OSRSetFromUserInput(hSRS, crs.c_str());
	if (erro == 4) {
		return false ;
	}
	char *pszSRS_WKT = NULL;
	OSRExportToWkt( hSRS, &pszSRS_WKT );
	OSRDestroySpatialReference( hSRS );
	GDALSetProjection( hDS, pszSRS_WKT );
	CPLFree( pszSRS_WKT );
	return true;
}


bool SpatRaster::as_gdalvrt(GDALDatasetH &hVRT, SpatOptions &opt) {
// all sources should be on disk
	GDALDriverH hDrv = GDALGetDriverByName("MEM");
	hVRT = GDALCreate(hDrv, "", ncol(), nrow(), nlyr(), GDT_Float64, NULL);

	std::vector<double> rs = resolution();
	SpatExtent extent = getExtent();

	double adfGeoTransform[6] = { extent.xmin, rs[0], 0, extent.ymax, 0, -1 * rs[1] };
	GDALSetGeoTransform(hVRT, adfGeoTransform);

	if (!GDALsetSRS(hVRT, source[0].srs.wkt)) {
		addWarning("cannot set SRS");
	}

	char** papszOptions = NULL;
	SpatRaster RS;
	GDALDatasetH DS;
	for (size_t i=0; i<nlyr(); i++) {
		RS = SpatRaster(source[i]);
		std::string filename = source[i].filename;
		if (!SpatRaster::open_gdal(DS, i, false, opt)) {
			setError("cannot open datasource");
			return false;
		}
		papszOptions = CSLSetNameValue(papszOptions, "SourceFilename", filename.c_str());
		size_t n = source[i].layers.size();
		for (size_t j=0; j<n; j++) {
			std::string sband = std::to_string(source[i].layers[j] + 1);
			papszOptions = CSLSetNameValue(papszOptions, "SourceBand", sband.c_str());
			GDALAddBand(hVRT, GDT_Float64, papszOptions);
		}
	}
	CSLDestroy(papszOptions);
	return true;
}



bool SpatRaster::open_gdal(GDALDatasetH &hDS, int src, bool update, SpatOptions &opt) {

	size_t isrc = src < 0 ? 0 : src;

	bool fromfile = !source[isrc].memory;

	if (fromfile & (nsrc() > 1) & (src < 0)) {
		if (canProcessInMemory(opt)) {
			fromfile = false;
		} else {
			// make VRT
			setError("right now this method can only handle one file source at a time");
			return false;
		}
	}

	if (fromfile) {
		std::string f;
		//if (source[src].parameters_changed) {
			// make a copy to get the write the new crs or extent
			// can we use a VRT instead?
		//	f = tempFile(opt.get_tempdir(), ".tif");
		//	SpatRaster tmp(source[src]);
		//	SpatOptions topt(opt);
		//	topt.set_filenames({f});
		//	tmp.writeRaster(topt);
		//} else {
		f = source[src].filename;
		//}
		//hDS = GDALOpenShared(f.c_str(), GA_ReadOnly);

		if (update) {
			hDS = openGDAL(f, GDAL_OF_RASTER | GDAL_OF_UPDATE | GDAL_OF_SHARED, source[src].open_drivers, source[src].open_ops);
		/*
			if (hDS != NULL) { // for user-set extents
				std::vector<double> rs = resolution();
				SpatExtent extent = getExtent();
				double adfGeoTransform[6] = { extent.xmin, rs[0], 0, extent.ymax, 0, -1 * rs[1] };
				GDALSetGeoTransform(hDS, adfGeoTransform);
			}
		*/
		} else {
			hDS = openGDAL(f, GDAL_OF_RASTER | GDAL_OF_READONLY | GDAL_OF_SHARED, source[src].open_drivers, source[src].open_ops);
		}
		return (hDS != NULL);

	} else { // in memory


		bool hasval = source[isrc].hasValues;

		size_t nl;
		if (src < 0) {
			nl = nlyr();
		} else {
			nl = source[src].layers.size();
		}
		size_t ncls = nrow() * ncol();
		GDALDriverH hDrv = GDALGetDriverByName("MEM");

/*https://gis.stackexchange.com/questions/196048/how-to-reuse-memory-pointer-of-gdal-memory-driver
		char **papszOptions = NULL;
		hDS = GDALCreate(hDrv, "", ncol(), nrow(), 0, GDT_Float64, papszOptions);
		if (hDS == NULL) return false;
		std::vector<double> vals;
		for(size_t i=0; i<nl; i++)	{
			size_t off = ncls * i;
			vals = std::vector<double>(source[0].values.begin() +off, source[0].values.begin() +off+ncls);
			char szPtrValue[128] = { '\0' };
			int nRet = CPLPrintPointer( szPtrValue, reinterpret_cast<void*>(&vals[0]), sizeof(szPtrValue) );
			szPtrValue[nRet] = 0;
			papszOptions = CSLSetNameValue(papszOptions, "DATAPOINTER", szPtrValue);
			GDALAddBand(hDS, GDT_Float64, papszOptions);
		}
		CSLDestroy(papszOptions);
*/

		size_t nr = nrow();
		size_t nc = ncol();
		hDS = GDALCreate(hDrv, "", nc, nr, nl, GDT_Float64, NULL);
		if (hDS == NULL) return false;

		std::vector<double> rs = resolution();
		SpatExtent extent = getExtent();
		double adfGeoTransform[6] = { extent.xmin, rs[0], 0, extent.ymax, 0, -1 * rs[1] };
		GDALSetGeoTransform(hDS, adfGeoTransform);

		if (!GDALsetSRS(hDS, source[0].srs.wkt)) {
			setError("cannot set SRS");
			return false;
		}

		CPLErr err = CE_None;

		if (hasval) {
			std::vector<std::string> nms;
			if (src < 0) {
				nms = getNames();
			} else {
				nms = source[src].names;
			}

			std::vector<double> vv, vals;
			if (src < 0) {
				vv = getValues(-1, opt);
			} else {
				if (!getValuesSource(src, vv)) {
					setError("cannot read from source");
					return false;
				}
			}

			for (size_t i=0; i < nl; i++) {
				GDALRasterBandH hBand = GDALGetRasterBand(hDS, i+1);
				GDALSetRasterNoDataValue(hBand, NAN);
				GDALSetDescription(hBand, nms[i].c_str());

				size_t offset = ncls * i;
				vals = std::vector<double>(vv.begin() + offset, vv.begin() + offset+ncls);
				err = GDALRasterIO(hBand, GF_Write, 0, 0, nc, nr, &vals[0], nc, nr, GDT_Float64, 0, 0);
				if (err != CE_None) {
					return false;
				}
			}
		}
	}

	return true;
}


bool SpatRaster::from_gdalMEM(GDALDatasetH hDS, bool set_geometry, bool get_values) {

	if (set_geometry) {
		SpatRasterSource s;
		s.ncol = GDALGetRasterXSize( hDS );
		s.nrow = GDALGetRasterYSize( hDS );
		s.nlyr = GDALGetRasterCount( hDS );

		double adfGeoTransform[6];
		if( GDALGetGeoTransform( hDS, adfGeoTransform ) != CE_None ) {
			setError("Cannot get geotransform");
			return false;
		}
		double xmin = adfGeoTransform[0];
		double xmax = xmin + adfGeoTransform[1] * s.ncol;
		double ymax = adfGeoTransform[3];
		double ymin = ymax + s.nrow * adfGeoTransform[5];
		s.extent = SpatExtent(xmin, xmax, ymin, ymax);

		s.memory = true;
		s.names = source[0].names;
		std::string wkt;

#if GDAL_VERSION_MAJOR >= 3

		std::string errmsg;
		OGRSpatialReferenceH srs = GDALGetSpatialRef( hDS );
		if (srs == NULL) {
			return false;
		}
		const char *options[3] = { "MULTILINE=YES", "FORMAT=WKT2", NULL };
		char *cp;
		OGRErr err = OSRExportToWktEx(srs, &cp, options);
		if (is_ogr_error(err, errmsg)) {
			CPLFree(cp);
			return false;
		}
		wkt = std::string(cp);
		CPLFree(cp);
#else
		const char *pszSrc = GDALGetProjectionRef( hDS );
		if (pszSrc != NULL) {
			wkt = std::string(pszSrc);
		} else {
			return false;
		}

		//OGRSpatialReferenceH srs = GDALGetProjectionRef( hDS );
		//OGRSpatialReference oSRS(poDataset->GetProjectionRef());
		//OGRErr err = oSRS.exportToPrettyWkt(&cp);
#endif
		std::string msg;
		if (!s.srs.set({wkt}, msg)) {
			setError(msg);
			return false;
		} else if (msg != "") {
			addWarning(msg);
		}

		setSource(s);

	}

	if (get_values) {
		source[0].values.resize(0);
		source[0].values.reserve(ncell() * nlyr());
		CPLErr err = CE_None;
		int hasNA;
		size_t nl = nlyr();
		for (size_t i=0; i < nl; i++) {
			GDALRasterBandH hBand = GDALGetRasterBand(hDS, i+1);
			std::vector<double> lyrout( ncell() );
			err = GDALRasterIO(hBand, GF_Read, 0, 0, ncol(), nrow(), &lyrout[0], ncol(), nrow(), GDT_Float64, 0, 0);
			if (err != CE_None ) {
				setError("CE_None");
				return false;
			}
			//double naflag = -3.4e+38;
			double naflag = GDALGetRasterNoDataValue(hBand, &hasNA);
			if (hasNA && (!std::isnan(naflag))) {
				if (naflag < -3.4e+37) {
					naflag = -3.4e+37;
					for (size_t i=0; i<lyrout.size(); i++) {
						if (lyrout[i] < naflag) {
							lyrout[i] = NAN;
						}
					}
				} else {
					std::replace(lyrout.begin(), lyrout.end(), naflag, (double) NAN);
				}
			}


			bool has_so = false;;
			int success;
			double moffset = GDALGetRasterOffset(hBand, &success);
			if (success) {
				if (moffset != 0) {
					has_so = true;
				}
			} else {
				moffset = 0;
			}

			double mscale = GDALGetRasterScale(hBand, &success);
			if (success) {
				if (mscale != 1) {
					has_so = true;
				}
			} else {
				mscale = 1;
			}
			if (has_so) {
				for (double &d : lyrout) { d = d * mscale + moffset;}
			}
			source[0].values.insert(source[0].values.end(), lyrout.begin(), lyrout.end());
		}

		source[0].hasValues = true;
		source[0].memory = true;
		source[0].driver = "memory";
		source[0].setRange();
	}

	return true;
}



char ** set_GDAL_options(std::string driver, double diskNeeded, bool writeRGB, std::vector<std::string> gdal_options) {

	char ** gdalops = NULL;
	if (driver == "GTiff") {
		bool lzw = true;
		bool compressed = true;
		for (size_t i=0; i<gdal_options.size(); i++) {
			if (gdal_options[i].substr(0, 8) == "COMPRESS") {
				lzw = false;
				if (gdal_options[i].substr(9, 4) == "NONE") {
					compressed = false;
				}
				break;
			}
		}
		if (lzw) {
			gdalops = CSLSetNameValue( gdalops, "COMPRESS", "LZW");
		}
		if (compressed & (diskNeeded > 4194304000)) {
			bool big = true;
			for (size_t i=0; i<gdal_options.size(); i++) {
				if (gdal_options[i].substr(0, 7) == "BIGTIFF") {
					big = false;
					break;
				}
			}
			if (big) {
				gdalops = CSLSetNameValue( gdalops, "BIGTIFF", "YES");
			}
		}
		if (writeRGB) {
			gdalops = CSLSetNameValue( gdalops, "PROFILE", "GeoTIFF");
		}
	} else if (writeRGB) {
		gdalops = CSLSetNameValue( gdalops, "PHOTOMETRIC", "RGB");
	}

	for (size_t i=0; i<gdal_options.size(); i++) {
		std::vector<std::string> gopt = strsplit(gdal_options[i], "=");
		if (gopt.size() == 2) {
			gdalops = CSLSetNameValue(gdalops, gopt[0].c_str(), gopt[1].c_str() );
		}
	}
	return gdalops;
}



bool SpatRaster::create_gdalDS(GDALDatasetH &hDS, std::string filename, std::string driver, bool fill, double fillvalue, std::vector<bool> has_so, std::vector<double> scale, std::vector<double> offset, SpatOptions& opt) {

	has_so.resize(nlyr(), false);

	const char *pszFormat = driver.c_str();
	GDALDriverH hDrv = GDALGetDriverByName(pszFormat);

	double naflag = opt.NAflag; //NAN;
	GDALDataType gdt;
	char **papszOptions = NULL;
	if (driver != "MEM") {
		std::string datatype = opt.get_datatype();
		if (!getGDALDataType(datatype, gdt)) {
			addWarning("unknown datatype = " + datatype);
			getGDALDataType("FLT4S", gdt);
		}
		int dsize = std::stoi(datatype.substr(3,1));
		GIntBig diskNeeded = ncell() * nlyr() * dsize;
		std::string dname = dirname(filename);
		GIntBig diskAvailable = VSIGetDiskFreeSpace(dname.c_str());
		if ((diskAvailable > -1) && (diskAvailable < diskNeeded)) {
			setError("insufficient disk space (perhaps from temporary files?)");
			return(false);
		}

		papszOptions = set_GDAL_options(driver, diskNeeded, false, opt.gdal_options);

		if (datatype == "INT4S") {
			naflag = INT32_MIN; //-2147483648;
		} else if (datatype == "INT2S") {
			naflag = INT16_MIN;
		} else if (datatype == "INT4U") {
			naflag = UINT32_MAX;
		} else if (datatype == "INT2U") {
			naflag = UINT16_MAX;
		} else if (datatype == "INT1U") {
			naflag = 255; // ?;
		}
	} else {
		getGDALDataType(opt.get_datatype(), gdt);
	}
	const char *pszFilename = filename.c_str();
	hDS = GDALCreate(hDrv, pszFilename, ncol(), nrow(), nlyr(), gdt, papszOptions );
	CSLDestroy( papszOptions );

	std::vector<std::string> nms = getNames();
	std::vector<bool> hasCats = hasCategories();

	for (size_t i=0; i < nlyr(); i++) {
		GDALRasterBandH hBand = GDALGetRasterBand(hDS, i+1);
		GDALSetDescription(hBand, nms[i].c_str());
		GDALSetRasterNoDataValue(hBand, naflag);
		//GDALSetRasterNoDataValue(hBand, -3.4e+38);
		if (fill) GDALFillRaster(hBand, fillvalue, 0);

		if (has_so[i]) {
			GDALSetRasterOffset(hBand, offset[i]);
			GDALSetRasterScale(hBand, scale[i]);
		}

		if (hasCats[i]) {
			std::vector<std::string> cats = getLabels(i);
			char **names = NULL;
			for (size_t j = 0; j < cats.size(); j++) {
				names = CSLAddString(names, cats[j].c_str());
			}
			CPLErr err = GDALSetRasterCategoryNames(hBand, names);
			if (err != CE_None) {
				addWarning("could not write categories");
			}
			CSLDestroy( names );
		}
	}

	std::vector<double> rs = resolution();
	SpatExtent e = getExtent();
	double adfGeoTransform[6] = { e.xmin, rs[0], 0, e.ymax, 0, -1 * rs[1] };
	GDALSetGeoTransform( hDS, adfGeoTransform);

	std::string wkt = getSRS("wkt");
	if (wkt != "") {
		OGRSpatialReferenceH hSRS = OSRNewSpatialReference( NULL );
		OGRErr erro = OSRSetFromUserInput(hSRS, wkt.c_str());
		if (erro == 4) {
			setError("CRS failure");
			OSRDestroySpatialReference( hSRS );
			return false;
		}
		char *pszSRS_WKT = NULL;
		OSRExportToWkt( hSRS, &pszSRS_WKT );
		GDALSetProjection( hDS, pszSRS_WKT );
		CPLFree(pszSRS_WKT);
		OSRDestroySpatialReference( hSRS );
	}
	return true;
}




/*

bool SpatRaster::setValuesMEM(GDALDatasetH hDS, bool set_geometry) {

	if (set_geometry) {
		SpatRasterSource s;
		s.ncol = GDALGetRasterXSize( hDS );
		s.nrow = GDALGetRasterYSize( hDS );
		s.nlyr = GDALGetRasterCount( hDS );

		double adfGeoTransform[6];
		if( GDALGetGeoTransform( hDS, adfGeoTransform ) != CE_None ) {
			setError("Cannot get geotransform");
			return false;
		}
		double xmin = adfGeoTransform[0];
		double xmax = xmin + adfGeoTransform[1] * s.ncol;
		double ymax = adfGeoTransform[3];
		double ymin = ymax + s.nrow * adfGeoTransform[5];
		s.extent = SpatExtent(xmin, xmax, ymin, ymax);

		s.driver = "memory";
		s.names = source[0].names;
		setSource(s);

		OGRSpatialReferenceH srs = GDALGetSpatialRef( hDS );
		char *cp;
		const char *options[3] = { "MULTILINE=YES", "FORMAT=WKT2", NULL };
		OGRErr err = OSRExportToWktEx(srs, &cp, options);
		std::string errmsg;
		if (is_ogr_error(err, errmsg)) {
			CPLFree(cp);
			return false;
		}
		std::string wkt = std::string(cp);
		CPLFree(cp);
		std::string msg;
		if (!s.srs.set({wkt}, msg)) {
			setError(msg);
			return false;
		}
	}


	source[0].values.reserve(ncell() * nlyr());
	CPLErr err = CE_None;
	int hasNA;
	for (size_t i=0; i < nlyr(); i++) {
		GDALRasterBandH hBand = GDALGetRasterBand(hDS, i+1);
		std::vector<double> lyrout( ncell() );
		err = GDALRasterIO(hBand, GF_Read, 0, 0, ncol(), nrow(), &lyrout[0], ncol(), nrow(), GDT_Float64, 0, 0);
		if (err != CE_None ) {
			setError("CE_None");
			return false;
		}

		//double naflag = -3.4e+38;
		double naflag = GDALGetRasterNoDataValue(hBand, &hasNA);
		if (hasNA) std::replace(lyrout.begin(), lyrout.end(), naflag, (double) NAN);
		source[0].values.insert(source[0].values.end(), lyrout.begin(), lyrout.end());

	}
	source[0].hasValues = TRUE;
	source[0].memory = TRUE;
	source[0].driver = "memory";
	source[0].setRange();
	return true;
}



bool gdal_ds_create(SpatRaster x, GDALDatasetH &hDS, std::string filename, std::string driver, std::vector<std::string> foptions, bool fill, std::string &msg) {

	msg = "";

	char **papszOptions = NULL;
	for (size_t i=0; i < foptions.size(); i++) {
		std::vector<std::string> wopt = strsplit(foptions[i], "=");
		if (wopt.size() == 2) {
			papszOptions = CSLSetNameValue( papszOptions, wopt[0].c_str(), wopt[1].c_str() );
		}
	}

	const char *pszFormat = driver.c_str();
	GDALDriverH hDrv = GDALGetDriverByName(pszFormat);

	const char *pszFilename = filename.c_str();
	hDS = GDALCreate(hDrv, pszFilename, x.ncol(), x.nrow(), x.nlyr(), GDT_Float64, papszOptions );
	CSLDestroy( papszOptions );

	GDALRasterBandH hBand;
	std::vector<std::string> nms = x.getNames();
	for (size_t i=0; i < x.nlyr(); i++) {
		hBand = GDALGetRasterBand(hDS, i+1);
		GDALSetDescription(hBand, nms[i].c_str());
		GDALSetRasterNoDataValue(hBand, NAN);
		//GDALSetRasterNoDataValue(hBand, -3.4e+38);
		if (fill) GDALFillRaster(hBand, NAN, 0);
	}

	std::vector<double> rs = x.resolution();
	SpatExtent e = x.getExtent();
	double adfGeoTransform[6] = { e.xmin, rs[0], 0, e.ymax, 0, -1 * rs[1] };
	GDALSetGeoTransform( hDS, adfGeoTransform);

	OGRSpatialReferenceH hSRS = OSRNewSpatialReference( NULL );
	std::vector<std::string> srs = x.getSRS();
	std::string wkt = srs[1];
	OGRErr erro = OSRSetFromUserInput(hSRS, wkt.c_str());
	if (erro == 4) {
		msg = "CRS failure";
		OSRDestroySpatialReference( hSRS );
		return false;
	}

	char *pszSRS_WKT = NULL;
	OSRExportToWkt( hSRS, &pszSRS_WKT );
	GDALSetProjection( hDS, pszSRS_WKT );
	CPLFree(pszSRS_WKT);
	OSRDestroySpatialReference( hSRS );
	return true;
}



bool SpatRaster::open_gdal(GDALDatasetH &hDS) {
	// needs to loop over sources. thus should vector of GDALDatasetH
	// for now just doing the first
	if (!source[0].hasValues) {
		return false;

	} else if (source[0].driver == "gdal") {
		std::string f = source[0].filename;
		hDS = GDALOpen(f.c_str(), GA_ReadOnly);
		return(hDS != NULL);

	} else { // in memory

		size_t nl = nlyr();
		size_t ncls = nrow() * ncol();
		GDALDriverH hDrv = GDALGetDriverByName("MEM");

//https://gis.stackexchange.com/questions/196048/how-to-reuse-memory-pointer-of-gdal-memory-driver
//		char **papszOptions = NULL;
//		hDS = GDALCreate(hDrv, "", ncol(), nrow(), 0, GDT_Float64, papszOptions);
//		if (hDS == NULL) return false;
//		std::vector<double> vals;
//		for(size_t i=0; i<nl; i++)	{
//			size_t off = ncls * i;
//			vals = std::vector<double>(source[0].values.begin() +off, source[0].values.begin() +off+ncls);
//			char szPtrValue[128] = { '\0' };
//			int nRet = CPLPrintPointer( szPtrValue, reinterpret_cast<void*>(&vals[0]), sizeof(szPtrValue) );
//			szPtrValue[nRet] = 0;
//			papszOptions = CSLSetNameValue(papszOptions, "DATAPOINTER", szPtrValue);
//			GDALAddBand(hDS, GDT_Float64, papszOptions);
//		}
//		CSLDestroy(papszOptions);


		size_t nr = nrow();
		size_t nc = ncol();
		hDS = GDALCreate(hDrv, "", nc, nr, nl, GDT_Float64, NULL);
		if (hDS == NULL) return false;

		std::vector<double> rs = resolution();
		double adfGeoTransform[6] = { extent.xmin, rs[0], 0, extent.ymax, 0, -1 * rs[1] };
		GDALSetGeoTransform(hDS, adfGeoTransform);

		OGRSpatialReferenceH hSRS = OSRNewSpatialReference( NULL );
		std::string crs = srs.wkt;
		OGRErr erro = OSRSetFromUserInput(hSRS, crs.c_str());
		if (erro == 4) {
			setError("CRS failure");
			return false ;
		}
		char *pszSRS_WKT = NULL;
		OSRExportToWkt( hSRS, &pszSRS_WKT );
		OSRDestroySpatialReference( hSRS );
		GDALSetProjection( hDS, pszSRS_WKT );
		CPLFree( pszSRS_WKT );

		CPLErr err = CE_None;
		std::vector<double> vals;

		std::vector<std::string> nms = getNames();

		for (size_t i=0; i < nl; i++) {
			GDALRasterBandH hBand = GDALGetRasterBand(hDS, i+1);
			GDALSetRasterNoDataValue(hBand, NAN);
			GDALSetDescription(hBand, nms[i].c_str());

			size_t offset = ncls * i;
			vals = std::vector<double>(source[0].values.begin() + offset, source[0].values.begin() + offset + ncls);
			err = GDALRasterIO(hBand, GF_Write, 0, 0, nc, nr, &vals[0], nc, nr, GDT_Float64, 0, 0 );
			if (err != CE_None) {
				return false;
			}
		}
	}

	return true;
}




*/

