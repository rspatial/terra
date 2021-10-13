// Copyright (c) 2018-2021  Robert J. Hijmans
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

#include "spatVector.h"

#ifdef useGDAL

#include "file_utils.h"
#include "ogrsf_frmts.h"

/*
bool SpatVector::ogr_geoms(std::vector<OGRGeometryH> &ogrgeoms, std::string &	message) {

	const char *pszDriverName = "memory";
	GDALDriverH hDriver;
	hDriver = GDALGetDriverByName( pszDriverName );
	GDALDatasetH hDS;
	hDS = GDALCreate( hDriver, "point_out.shp", 0, 0, 0, GDT_Unknown, NULL );
	if (hDS == NULL) {
		return false;
	}

//    char **papszMetadata;

	OGRwkbGeometryType wkb;
	SpatGeomType geomtype = geoms[0].gtype;
	if (geomtype == points) {
		wkb = wkbPoint;
	} else if (geomtype == lines) {
		wkb = wkbMultiLineString;
	} else if (geomtype == polygons) {
		wkb = wkbMultiPolygon;
	} else {
        setError("this geometry type is not supported");
        return false;		
	}

	std::string s = srs.wkt;

//	OGRSpatialReferenceH hSRS = OSRNewSpatialReference( NULL );
//	OGRErr erro = OSRSetFromUserInput(hSRS, s.c_str());
//	if (erro != 4) {
//		return false ;
//		char *pszSRS_WKT = NULL;
//		OSRExportToWkt( hSRS, &pszSRS_WKT );
//		OSRDestroySpatialReference( hSRS );
//		GDALSetProjection( hDS, pszSRS_WKT );
//		CPLFree( pszSRS_WKT );
//	}


	OGRSpatialReference *SRS = NULL;
	if (s != "") {
		SRS = new OGRSpatialReference;
		OGRErr err = SRS->SetFromUserInput(s.c_str()); 
		if (err != OGRERR_NONE) {
			setError("crs error");
			delete SRS;
			return false;
		}
	}

	OGRLayerH hLayer;
	hLayer = GDALDatasetCreateLayer( hDS, "", SRS, wkb, NULL );
    if( hLayer == NULL ) {
        message = "Layer creation failed" ;
        return false;
    }

	if (SRS != NULL) SRS->Release();

	for (size_t i=0; i<ngeoms; i++) {
	
	    OGRFeatureH hFeature;
        OGRGeometryH hPt;

        hFeature = OGR_F_Create( OGR_L_GetLayerDefn( hLayer ) );
		
	
// points -- also need to do multipoints
		OGRGeometryH hPt;
		if (wkb == wkbPoint) {
			SpatGeom g = getGeom(i);
			OGR_G_SetPoint_2D(hPt, 0, g.parts[0].x[0], g.parts[0].y[0]);
			ogrgeoms.push_back(hPt);
// lines		
		} else if (wkb == wkbMultiLineString) {
			SpatGeom g = getGeom(i);
			OGRMultiLineString poGeom;
			for (size_t j=0; j<g.size(); j++) {
				OGRLineString poLine = OGRLineString();
				SpatPart p = g.getPart(j);
				for (size_t k=0; k<p.size(); k++) {
					pt.setX(p.x[k]);
					pt.setY(p.y[k]);
					poLine.setPoint(k, &pt);
				}
				if (poGeom.addGeometry(&poLine) != OGRERR_NONE ) {
					message = "cannot add line";
					return false;
				}
			}
			if (poFeature->SetGeometry( &poGeom ) != OGRERR_NONE) {
				message = "cannot set geometry";
				return false;
			}
		
// polygons		
		} else if (wkb == wkbMultiPolygon) {
			SpatGeom g = getGeom(i);
			OGRPolygon poGeom;
			for (size_t j=0; j<g.size(); j++) {
				OGRLinearRing poRing;
				SpatPart p = g.getPart(j);
				for (size_t k=0; k<p.size(); k++) {
					pt.setX(p.x[k]);
					pt.setY(p.y[k]);
					poRing.setPoint(k, &pt);
				}
				if (poGeom.addRing(&poRing) != OGRERR_NONE ) {
					message = "cannot add ring";
					return false;
				}
			
				if (p.hasHoles()) {
					for (size_t h=0; h < p.nHoles(); h++) {
						SpatHole hole = p.getHole(h);
						OGRLinearRing poHole;
						for (size_t k=0; k<hole.size(); k++) {
							pt.setX(hole.x[k]);
							pt.setY(hole.y[k]);
							poHole.setPoint(k, &pt);
						}					
						if (poGeom.addRing(&poHole) != OGRERR_NONE ) {
							message = "cannot add hole";
							return false;
						}
					}
				}
				//closeRings
			}
			if (poFeature->SetGeometry( &poGeom ) != OGRERR_NONE) {
				message = "cannot set geometry";
				return false;
			}
		} else {
			message = "Only points, lines and polygons are currently supported";
			return false;
		}
	
		if( poLayer->CreateFeature( poFeature ) != OGRERR_NONE ) {
			message = "Failed to create feature";
			return false;
        }

        OGRFeature::DestroyFeature( poFeature );
    }
    //GDALClose( poDS );
	//return true;
	return poDS;
}
*/



GDALDataset* SpatVector::write_ogr(std::string filename, std::string lyrname, std::string driver, bool overwrite) {


    GDALDataset *poDS = NULL;

	if (filename != "") {
		if (file_exists(filename) & (!overwrite)) {
			setError("file exists. Use 'overwrite=TRUE' to overwrite it");
			return(poDS);
		}
		if (nrow() == 0) {
			setError("no geometries to write");
			return(poDS);		
		}
	}
    GDALDriver *poDriver = GetGDALDriverManager()->GetDriverByName( driver.c_str() );
    if( poDriver == NULL )  {
        setError( driver + " driver not available");
        return poDS;
    }
    char **papszMetadata;
    papszMetadata = poDriver->GetMetadata();
    if (!CSLFetchBoolean( papszMetadata, GDAL_DCAP_VECTOR, FALSE)) {
		setError(driver + " is not a vector format");
        return poDS;
	}
    if (!CSLFetchBoolean( papszMetadata, GDAL_DCAP_CREATE, FALSE)) {
		setError("cannot create a "+ driver + " dataset");
        return poDS;
	}
    poDS = poDriver->Create(filename.c_str(), 0, 0, 0, GDT_Unknown, NULL );
    if( poDS == NULL ) {
        setError("Creation of output dataset failed" );
        return poDS;
    }

	OGRwkbGeometryType wkb;
	SpatGeomType geomtype = geoms[0].gtype;
	if (geomtype == points) {
		wkb = wkbPoint;
	} else if (geomtype == lines) {
		wkb = wkbMultiLineString;
	} else if (geomtype == polygons) {
		wkb = wkbMultiPolygon;
	} else {
        setError("this geometry type is not supported: " + type());
        return poDS;		
	}

	std::string s = srs.wkt;

	OGRSpatialReference *SRS = NULL;
	if (s != "") {
		SRS = new OGRSpatialReference;
		OGRErr err = SRS->SetFromUserInput(s.c_str()); 
		if (err != OGRERR_NONE) {
			setError("crs error");
			delete SRS;
			return poDS;
		}
	}

    OGRLayer *poLayer;
    poLayer = poDS->CreateLayer(lyrname.c_str(), SRS, wkb, NULL );
    if( poLayer == NULL ) {
        setError( "Layer creation failed" );
        return poDS;
    }
//	if (SRS != NULL) SRS->Release();
	if (SRS != NULL) OSRDestroySpatialReference(SRS);
	
	std::vector<std::string> nms = get_names();
	std::vector<std::string> tps = df.get_datatypes();
	OGRFieldType otype;
	int nfields = nms.size();
	size_t ngeoms = size();

	for (int i=0; i<nfields; i++) {
		if (tps[i] == "double") {
			otype = OFTReal;
		} else if (tps[i] == "long") {
			otype = OFTInteger64;
		} else {
			otype = OFTString;
		}

		OGRFieldDefn oField(nms[i].c_str(), otype);
		if (otype == OFTString) {
			oField.SetWidth(32); // needs to be computed
		}
		if( poLayer->CreateField( &oField ) != OGRERR_NONE ) {
			setError( "Field creation failed for: " + nms[i]);
			return poDS;
		}
	}

	//unsigned r = 0;

	for (size_t i=0; i<ngeoms; i++) {
	
		OGRFeature *poFeature;
        poFeature = OGRFeature::CreateFeature( poLayer->GetLayerDefn() );
		for (int j=0; j<nfields; j++) {
			if (tps[j] == "double") {
				poFeature->SetField(j, df.getDvalue(i, j));
			} else if (tps[j] == "long") {
				poFeature->SetField(j, (GIntBig)df.getIvalue(i, j));
			} else {
				poFeature->SetField(j, df.getSvalue(i, j).c_str());
			}
		}
		//r++;

// points -- also need to do multi-points
		OGRPoint pt;
		if (wkb == wkbPoint) {
			SpatGeom g = getGeom(i);
			if (!std::isnan(g.parts[0].x[0])) {
				pt.setX( g.parts[0].x[0] );
				pt.setY( g.parts[0].y[0] );
			}
			poFeature->SetGeometry( &pt );
		
// lines		
		} else if (wkb == wkbMultiLineString) {
			SpatGeom g = getGeom(i);
			OGRMultiLineString poGeom;
			for (size_t j=0; j<g.size(); j++) {
				OGRLineString poLine = OGRLineString();
				SpatPart p = g.getPart(j);
				for (size_t k=0; k<p.size(); k++) {
					if (!std::isnan(p.x[k])) {
						pt.setX(p.x[k]);
						pt.setY(p.y[k]);
						poLine.setPoint(k, &pt);
					}
				}
				if (poGeom.addGeometry(&poLine) != OGRERR_NONE ) {
					setError("cannot add line");
					return poDS;
				}
			}
			if (poFeature->SetGeometry( &poGeom ) != OGRERR_NONE) {
				setError("cannot set geometry");
				return poDS;
			}
		
// polygons		
		} else if (wkb == wkbMultiPolygon) {
			SpatGeom g = getGeom(i);
			OGRMultiPolygon poGeom;
			for (size_t j=0; j<g.size(); j++) {
				OGRLinearRing poRing;
				SpatPart p = g.getPart(j);
				for (size_t k=0; k<p.size(); k++) {
					if (!std::isnan(p.x[k])) {
						pt.setX(p.x[k]);
						pt.setY(p.y[k]);
						poRing.setPoint(k, &pt);
					}
				}
				OGRPolygon polyGeom;
				if (polyGeom.addRing(&poRing) != OGRERR_NONE ) {
					setError("cannot add ring");
					return poDS;
				}
			
				if (p.hasHoles()) {
					for (size_t h=0; h < p.nHoles(); h++) {
						SpatHole hole = p.getHole(h);
						OGRLinearRing poHole;
						for (size_t k=0; k<hole.size(); k++) {
							pt.setX(hole.x[k]);
							pt.setY(hole.y[k]);
							poHole.setPoint(k, &pt);
						}					
						if (polyGeom.addRing(&poHole) != OGRERR_NONE ) {
							setError("cannot add hole");
							return poDS;
						}
					}
				}
				poGeom.addGeometry( &polyGeom);
				//closeRings
			}
			
			//OGRMultiPolygon* mGeom = poGeom.toMultiPolygon();	
			if (poFeature->SetGeometry( &poGeom ) != OGRERR_NONE) {
				setError("cannot set geometry");
				return poDS;
			}
		} else {
			setError("Only points, lines and polygons are currently supported");
			return poDS;
		}
	
		if( poLayer->CreateFeature( poFeature ) != OGRERR_NONE ) {
			setError("Failed to create feature");
			return poDS;
        }

        OGRFeature::DestroyFeature( poFeature );
    }
    //GDALClose( poDS );
	//return true;
	return poDS;
}



bool SpatVector::write(std::string filename, std::string lyrname, std::string driver, bool overwrite) {

	GDALDataset *poDS = write_ogr(filename, lyrname, driver, overwrite);
    if (poDS != NULL) GDALClose( poDS );
	if (hasError()) {
		return false;
	} 
	return true;

}

GDALDataset* SpatVector::GDAL_ds() {
	return write_ogr("", "layer", "Memory", true);
}


#include <fstream>

bool SpatDataFrame::write_dbf(std::string filename, bool overwrite, SpatOptions &opt) {
// filename is here "raster.tif"
// to write "raster.tif.vat.dbf"

	if (filename != "") {
		if (file_exists(filename) & (!overwrite)) {
			setError("file exists. Use 'overwrite=TRUE' to overwrite it");
			return(false);
		}
		if (nrow() == 0) {
			setError("nothing to write");
			return(false);		
		}
	}

	std::string fbase = tempFile(opt.get_tempdir(), "");
	std::string f = fbase + ".shp";
    GDALDriver *poDriver = GetGDALDriverManager()->GetDriverByName( "ESRI Shapefile" );
    GDALDataset *poDS = NULL;
    poDS = poDriver->Create(f.c_str(), 0, 0, 0, GDT_Unknown, NULL );
    if( poDS == NULL ) {
        setError("Creation of output dataset failed" );
        return false;
    }

	OGRwkbGeometryType wkb = wkbPoint;
	OGRSpatialReference *SRS = NULL;
    OGRLayer *poLayer;
    poLayer = poDS->CreateLayer("dbf", SRS, wkb, NULL );
    if( poLayer == NULL ) {
        setError( "Layer creation failed" );
        return false;
    }
	std::vector<std::string> nms = get_names();
	std::vector<std::string> tps = get_datatypes();
	OGRFieldType otype;
	int nfields = nms.size();

	for (int i=0; i<nfields; i++) {
		if (tps[i] == "double") {
			otype = OFTReal;
		} else if (tps[i] == "long") {
			otype = OFTInteger64;
		} else {
			otype = OFTString;
		}

		OGRFieldDefn oField(nms[i].c_str(), otype);
		if (otype == OFTString) {
			oField.SetWidth(32); // needs to be computed
		}
		if( poLayer->CreateField( &oField ) != OGRERR_NONE ) {
			setError( "Field creation failed for: " + nms[i]);
			return false;
		}
	}

	for (size_t i=0; i<nrow(); i++) {
	
		OGRFeature *poFeature;
        poFeature = OGRFeature::CreateFeature( poLayer->GetLayerDefn() );
		for (int j=0; j<nfields; j++) {
			if (tps[j] == "double") {
				poFeature->SetField(j, getDvalue(i, j));
			} else if (tps[j] == "long") {
				poFeature->SetField(j, (GIntBig) getIvalue(i, j));
			} else {
				poFeature->SetField(j, getSvalue(i, j).c_str());
			}
		}

		OGRPoint pt;
		pt.setX( 0.0 );
		pt.setY( 0.0 );
		poFeature->SetGeometry( &pt );	
	
		if( poLayer->CreateFeature( poFeature ) != OGRERR_NONE ) {
			setError("Failed to create feature");
			return false;
        }

        OGRFeature::DestroyFeature( poFeature );
    }
    GDALClose( poDS );
	f = fbase + ".dbf";
	filename += ".vat.dbf";
	// c++17 has file_copy
    std::ifstream  src(f.c_str(), std::ios::binary);
    std::ofstream  dst(filename.c_str(),  std::ios::binary);
    dst << src.rdbuf();
	
	filename.erase(filename.length()-3);
	filename += "cpg";
	std::ofstream cpg;
	cpg.open (filename.c_str());
	cpg << "UTF-8";
	cpg.close();
	
	return true;
}


#endif