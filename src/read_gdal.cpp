#include "spatraster.h"
#include "util.h"

#include "gdal_priv.h"
#include "cpl_conv.h" // for CPLMalloc()
#include "cpl_string.h"
#include "ogr_spatialref.h"


bool SpatRaster::constructFromFileGDAL(std::string fname) {

    GDALDataset  *poDataset;
    GDALAllRegister();
	const char* pszFilename = fname.c_str();
    poDataset = (GDALDataset *) GDALOpen( pszFilename, GA_ReadOnly );
	
    if( poDataset == NULL )  {	return false;}	
		
	RasterSource s;
	s.ncol = poDataset->GetRasterXSize();
	s.nrow = poDataset->GetRasterYSize();
	s.nlyr = poDataset->GetRasterCount();
	
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

		int hasNA;
		double naflag = poBand->GetNoDataValue(&hasNA);
		if (hasNA) { 
			s.NAflag = naflag; 
		} else {
			s.NAflag = NAN;
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
			s.hasCT.push_back(true);
		} else {
			s.hasCT.push_back(false);
		}
		
		GDALRasterAttributeTable *rat = poBand->GetDefaultRAT(); 	
		if( rat != NULL )	{  // does not appear to work
			s.hasRAT.push_back(true);
		} else {
			s.hasRAT.push_back(false);
		}
		
		s.names.push_back( "lyr" + std::to_string(i+1) ) ;
	}
	GDALClose( (GDALDatasetH) poDataset );

	s.hasValues = true;
	setSource(s);
	return true; 
}


bool SpatRaster::readStartGDAL() {
    GDALDataset  *poDataset;
    GDALAllRegister();	
	const char* pszFilename = source[0].filename.c_str();
	poDataset = (GDALDataset *) GDALOpen( pszFilename, GA_ReadOnly );
    gdalconnection = poDataset; 
	open_read = true;
	return(true);
}


std::vector<double> SpatRaster::readChunkGDAL(unsigned row, unsigned nrows, unsigned col, unsigned ncols, unsigned lyr, unsigned nlyrs) {	
	GDALRasterBand  *poBand;		
	unsigned ncell = ncols * nrows;
	unsigned cell;
	std::vector<double> errout;
	int hasNA;
	std::vector<double> out(ncell * nlyrs);
	for (size_t i=0; i < nlyrs; i++) {		
		cell = ncell * i;
		poBand = gdalconnection->GetRasterBand(lyr + i + 1);
		GDALDataType gdtype = poBand->GetRasterDataType();
		if (gdtype == GDT_Float64) {
			double naflag = poBand->GetNoDataValue(&hasNA);
			CPLErr err = poBand->RasterIO(GF_Read, col, row, ncols, nrows, &out[cell], ncols, nrows, GDT_Float64, 0, 0);
			if (err == 4) {
				setError("cannot read values");
				return errout;
			}
			if (hasNA) {
				double navalue = NA<double>::value;
				std::replace(out.begin()+cell, out.end()+cell+ncell, naflag, navalue);
			}
		} else {
			
		}
	}
	return(out);
}


bool SpatRaster::readStopGDAL() {
	GDALClose( (GDALDatasetH) gdalconnection);
	open_read = false;	
	return true;
}


std::vector<double> SpatRaster::readValuesGDAL(unsigned row, unsigned nrows, unsigned col, unsigned ncols, unsigned lyr, unsigned nlyrs) {
    GDALDataset *poDataset;
	GDALRasterBand *poBand;
    GDALAllRegister();
	const char* pszFilename = source[0].filename.c_str();
    poDataset = (GDALDataset *) GDALOpen(pszFilename, GA_ReadOnly);
	unsigned ncell = ncols * nrows;
	unsigned cell;
	std::vector<double> errout;
	std::vector<double> out(ncell*nlyrs);
	int hasNA;
	
	for (size_t i=0; i < nlyrs; i++) {		
		cell = ncell * i;
		poBand = poDataset->GetRasterBand(lyr + i + 1);
		GDALDataType gdtype = poBand->GetRasterDataType();
		double naflag = poBand->GetNoDataValue(&hasNA);
		if (!hasNA) { naflag = NAN; }
		if (gdtype == GDT_Float64) {
			CPLErr err = poBand->RasterIO(GF_Read, col, row, ncols, nrows, &out[cell], ncols, nrows, gdtype, 0, 0);
			if (err == 4) {
				setError("cannot read values");
				GDALClose((GDALDatasetH) poDataset);	
				return errout;
			} 
			setNAN(out, naflag);
		} else if (gdtype == GDT_Float32) {
			std::vector<float> lyrout(ncell);
			CPLErr err = poBand->RasterIO(GF_Read, col, row, ncols, nrows, &lyrout[cell], ncols, nrows, gdtype, 0, 0);
			if (err == 4) {
				setError("cannot read values");
				GDALClose((GDALDatasetH) poDataset);	
				return errout;
			}
			setNAN(lyrout, naflag);
			for (size_t j=0; j<ncell; j++) {
				out[cell+j] = lyrout[j];
			}
		} else {
			//int tbd
		}
	}
	GDALClose((GDALDatasetH) poDataset);	
	return out;
}



std::vector<double> SpatRaster::readPointsGDAL(const std::vector<double> &x, const std::vector<double> &y) {
    GDALDataset *poDataset;
	GDALRasterBand *poBand;
    GDALAllRegister();
	const char* pszFilename = source[0].filename.c_str();
    poDataset = (GDALDataset *) GDALOpen(pszFilename, GA_ReadOnly);

	std::vector<double> cols = colFromX(x);
	std::vector<double> rows = rowFromY(y);

	unsigned n = rows.size();
	std::vector<double> errout;
	std::vector<double> out(n*nlyr());
	int hasNA;
	unsigned offset;
	for (size_t i=0; i < nlyr(); i++) {		
		offset = n * i;
		poBand = poDataset->GetRasterBand(i + 1);
		GDALDataType gdtype = poBand->GetRasterDataType();
		double naflag = poBand->GetNoDataValue(&hasNA);
		if (!hasNA) { naflag = NAN; }
		if (gdtype == GDT_Float64) {
			for (size_t j=0; j < n; j++) {		
				CPLErr err = poBand->RasterIO(GF_Read, cols[j], rows[j], 1, 1, &out[offset+j], 1, 1, gdtype, 0, 0);
				if (err == 4) {
					setError("cannot read values");
					GDALClose((GDALDatasetH) poDataset);	
					return errout;
				} 
			}
			setNAN(out, naflag);
		} else if (gdtype == GDT_Float32) {
			std::vector<float> lyrout(n);
			for (size_t j=0; j < n; j++) {		
				CPLErr err = poBand->RasterIO(GF_Read, cols[j], rows[j], 1, 1, &lyrout[offset+j], 1, 1, gdtype, 0, 0);
				if (err == 4) {
					setError("cannot read values");
					GDALClose((GDALDatasetH) poDataset);	
					return errout;
				}
			}
			setNAN(lyrout, naflag);
			for (size_t j=0; j<out.size(); j++) {
				out[offset+j] = lyrout[j];
			}
		} else {
			//int tbd
		}
	}
	GDALClose((GDALDatasetH) poDataset);	
	return out;
}

