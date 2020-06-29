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
#include "string_utils.h"

#include "gdal_alg.h"
#include "ogrsf_frmts.h"

#if GDAL_VERSION_MAJOR >= 3

SpatVector SpatRaster::polygonize(bool trunc, SpatOptions &opt) {

	SpatVector out;
	SpatOptions topt(opt);
	SpatRaster tmp = subset({0}, topt);

	// to vectorize all values that are not NAN (or Inf)
	// we could skip this if we know that min(tmp) > 0
	bool usemask;
	std::vector<double> rmin = tmp.range_min();
	SpatRaster mask;
	if (std::isnan(rmin[0]) || rmin[0] > 0) {
		usemask = false;
	} else {
		usemask = true;
		mask = tmp.isfinite(opt);		
	}
	GDALDatasetH rstDS;
	if (! tmp.sources_from_file() ) {
		if (!tmp.open_gdal(rstDS, 0)) {
			out.setError("cannot open dataset");
			return out;
		}
	} else {
		std::string filename = tmp.source[0].filename;
		rstDS = GDALOpen( filename.c_str(), GA_ReadOnly);
		if (rstDS == NULL) {
			out.setError("cannot open dataset from file");
			return out;			
		}
	}
    GDALDataset *srcDS=NULL;
	srcDS = srcDS->FromHandle(rstDS);


	GDALDatasetH rstMask;
	GDALDataset *maskDS=NULL;
	if (usemask) {
		if (! mask.sources_from_file() ) {
			if (!mask.open_gdal(rstMask, 0)) {
				out.setError("cannot open dataset");
				return out;
			}
		} else {
			std::string filename = mask.source[0].filename;
			rstMask = GDALOpen( filename.c_str(), GA_ReadOnly);
			if (rstMask == NULL) {
				out.setError("cannot open dataset from file");
				return out;			
			}
		}
		maskDS = srcDS->FromHandle(rstMask);
	}
	
    GDALDataset *poDS = NULL;
    GDALDriver *poDriver = GetGDALDriverManager()->GetDriverByName( "Memory" );
    if( poDriver == NULL )  {
        out.setError( "cannot create output dataset");
        return out;
    }
    poDS = poDriver->Create("", 0, 0, 0, GDT_Unknown, NULL );
    if( poDS == NULL ) {
        out.setError("Creation of dataset failed" );
        return out;
    }
	std::vector<std::string> nms = getNames();
	std::string name = nms[0];

	OGRSpatialReference *SRS = NULL;
	std::string s = srs.wkt;
	if (s != "") {
		SRS = new OGRSpatialReference;
		OGRErr err = SRS->SetFromUserInput(s.c_str()); 
		if (err != OGRERR_NONE) {
			out.setError("crs error");
			delete SRS;
			return out;
		}
	}

    OGRLayer *poLayer;	
    poLayer = poDS->CreateLayer(name.c_str(), SRS, wkbPolygon, NULL );
    if( poLayer == NULL ) {
        out.setError( "Layer creation failed" );
        return out;
    }
	if (SRS != NULL) SRS->Release();

	OGRFieldDefn oField(name.c_str(), trunc ?  OFTInteger : OFTReal);
	if( poLayer->CreateField( &oField ) != OGRERR_NONE ) {
		out.setError( "Creating field failed");
		return out;
	}

	GDALRasterBand  *poBand;
	poBand = srcDS->GetRasterBand(1);
	//int hasNA=1;
	//poBand->GetNoDataValue(&hasNA);

	CPLErr err;	
	if (usemask) {
		GDALRasterBand  *maskBand;
		maskBand = maskDS->GetRasterBand(1);
		if (trunc) {
			err = GDALPolygonize(poBand, maskBand, poLayer, 0, NULL, NULL, NULL);
		} else {
			err = GDALFPolygonize(poBand, maskBand, poLayer, 0, NULL, NULL, NULL);
		}
		GDALClose(maskDS);
	} else {
		if (trunc) {
			err = GDALPolygonize(poBand, poBand, poLayer, 0, NULL, NULL, NULL);
		} else {
			err = GDALFPolygonize(poBand, poBand, poLayer, 0, NULL, NULL, NULL);
		}
	}
	if (err == 4) {
		out.setError("polygonize error");
		return out;
	}
	GDALClose(srcDS);

	out.read_ogr(poDS);
	GDALClose(poDS);

	out = out.aggregate(name, false);
	
	return out;
}

#else
	
SpatVector SpatRaster::polygonize(bool trunc, SpatOptions &opt) {
	SpatVector out;
	out.setError("not supported with your version of GDAL");
	return out;
}

#endif
