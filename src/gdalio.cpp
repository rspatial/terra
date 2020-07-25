#include "ogr_spatialref.h"

#include "spatRaster.h"
#include "string_utils.h"
#include "crs.h"


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
		setError("cannot set SRS");
		return false;
	}
	
	char** papszOptions = NULL;
	SpatRaster RS;
	GDALDatasetH DS;
	for (size_t i=0; i<nlyr(); i++) {
		RS = SpatRaster(source[i]);
		std::string filename = source[i].filename;
		if (!SpatRaster::open_gdal(DS, i, opt)) {
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


bool SpatRaster::to_memory() {
	setValues(getValues());
	return true;
}

SpatRaster SpatRaster::to_memory_copy() {
	SpatRaster m = geometry();
	m.setValues(getValues());	
	return m;
}



bool SpatRaster::open_gdal(GDALDatasetH &hDS, int src, SpatOptions &opt) {
	// needs to loop over sources. 
	// Should be a vector of GDALDatasetH
	// Or can we combine them here into a VRT?
	// for now just doing the first
	
	size_t isrc = src < 0 ? 0 : src;
		
	bool hasval = source[isrc].hasValues;
	bool fromfile = !source[isrc].memory;

	if (fromfile & (nsrc() > 1) & (src < 0)) {
		if (canProcessInMemory(4, opt)) {
			fromfile = false;
		} else {
			// make VRT
			setError("right now this method can only handle one file source at a time");
			return false;
		}
	}
	
	if (fromfile) {

		std::string f = source[src].filename;
		hDS = GDALOpenShared(f.c_str(), GA_ReadOnly);
		return(hDS != NULL);
		
	} else { // in memory
				
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
				vv = getValues();
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
		RasterSource s;
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
		}
		
		setSource(s);
		
	}
	
	if (get_values) {
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
			if (hasNA && (!std::isnan(naflag))) {
				std::replace(lyrout.begin(), lyrout.end(), naflag, (double) NAN);
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


bool SpatRaster::create_gdalDS(GDALDatasetH &hDS, std::string filename, std::string driver, bool fill, double fillvalue, SpatOptions& opt) {

	std::vector<std::string> foptions = opt.gdal_options;
	char **papszOptions = NULL;
	for (size_t i=0; i < foptions.size(); i++) {
		std::vector<std::string> wopt = strsplit(foptions[i], "=");
		if (wopt.size() == 2) {
			papszOptions = CSLSetNameValue( papszOptions, wopt[0].c_str(), wopt[1].c_str() );
		}
	}

	const char *pszFormat = driver.c_str();
	GDALDriverH hDrv = GDALGetDriverByName(pszFormat);

	std::string datatype = opt.get_datatype();
	GDALDataType gdt;
	if (!getGDALDataType(datatype, gdt)) {
		addWarning("unknown datatype = " + datatype);
	}

	const char *pszFilename = filename.c_str();
	hDS = GDALCreate(hDrv, pszFilename, ncol(), nrow(), nlyr(), gdt, papszOptions );
	CSLDestroy( papszOptions );

	GDALRasterBandH hBand;
	std::vector<std::string> nms = getNames();
	for (size_t i=0; i < nlyr(); i++) {
		hBand = GDALGetRasterBand(hDS, i+1);
		GDALSetDescription(hBand, nms[i].c_str());
		GDALSetRasterNoDataValue(hBand, NAN);
		//GDALSetRasterNoDataValue(hBand, -3.4e+38);
		if (fill) GDALFillRaster(hBand, fillvalue, 0);
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
		RasterSource s;
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

