
#include "spatRaster.h"
#include "string_utils.h"

#include "gdal_alg.h"
#include "ogrsf_frmts.h"


SpatVector SpatRaster::polygonize(bool trunc) {

	SpatVector out;
	SpatOptions opt;

	SpatRaster tmp = subset({0}, opt);

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
    GDALDataset *srcDS;
	srcDS = srcDS->FromHandle(rstDS);

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
		out.setError( "Creating Name field failed");
		return out;
	}

	GDALRasterBand  *poBand;
	poBand = srcDS->GetRasterBand(1);

	//char **papszOptions = NULL;
	//if (queen) papszOptions = CSLSetNameValue(papszOptions, "8CONNECTED", "-8");
	
	CPLErr err;	
	if (trunc) {
		err = GDALPolygonize(poBand, poBand, poLayer, 0, NULL, NULL, NULL);
	} else {
		err = GDALFPolygonize(poBand, poBand, poLayer, 0, NULL, NULL, NULL);
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

