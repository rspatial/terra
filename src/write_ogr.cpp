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

#include "spatVector.h"

#ifdef useGDAL

#include "file_utils.h"
#include "ogrsf_frmts.h"

bool SpatVector::write_ogr(std::string filename, std::string lyrname, std::string driver, std::string &msg) {
	
    const char *pszDriverName = driver.c_str(); //"ESRI Shapefile";
    GDALDriverH hDriver;
    GDALDatasetH hDS;
    OGRLayerH hLayer;
    OGRFieldDefnH hFieldDefn;
    double x, y;
    char szName[33];

    GDALAllRegister();

    hDriver = GDALGetDriverByName( pszDriverName );
    if( hDriver == NULL ) {
        setError((std::string)pszDriverName + " driver not available");
        return false;
    }

    hDS = GDALCreate( hDriver, filename.c_str(), 0, 0, 0, GDT_Unknown, NULL );
    if( hDS == NULL ) {
        setError("Creation of output file failed" );
        return false;
    }

	OGRwkbGeometryType wkb;
	SpatGeomType geomtype = lyr.geoms[0].gtype;
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
	
	hLayer = GDALDatasetCreateLayer( hDS, lyrname.c_str(), NULL, wkb, NULL );
    if( hLayer == NULL ) {
        setError( "Layer creation failed" );
        return false;
    }

	std::vector<std::string> nms = get_names();
	std::vector<std::string> tps = lyr.df.get_datatypes();
	OGRFieldType otype;
	
	for (size_t i=0; i<nms.size(); i++) {
		
		if (tps[i] == "double") {
			otype = OFTReal;
		} else if (tps[i] == "long") {
			otype = OFTInteger64;
		} else {
			otype = OFTString;
		}
			
		hFieldDefn = OGR_Fld_Create( nms[i].c_str(), otype);
		if (otype == OFTString) {
			OGR_Fld_SetWidth( hFieldDefn, 32); // needs to be computed
		}

		if( OGR_L_CreateField( hLayer, hFieldDefn, TRUE ) != OGRERR_NONE ) {
			setError( "Creating Name field failed" );
			OGR_Fld_Destroy(hFieldDefn);
			GDALClose( hDS );
			return false;
		}

		OGR_Fld_Destroy(hFieldDefn);
	}
	
    while( !feof(stdin) && fscanf( stdin, "%lf,%lf,%32s", &x, &y, szName ) == 3 ) {
        OGRFeatureH hFeature;
        OGRGeometryH hPt;

        hFeature = OGR_F_Create( OGR_L_GetLayerDefn( hLayer ) );

/*		for (size_t i=0; i<nms.size(); i++) {
			if (tps[i] == "double") {
				OGR_F_SetFieldString( hFeature, OGR_F_GetFieldIndex(hFeature, nms[i]), szName );
			} else if (tps[i] == "long") {
				OGR_F_SetFieldString( hFeature, OGR_F_GetFieldIndex(hFeature, nms[i]), szName );
			} else {
				OGR_F_SetFieldString( hFeature, OGR_F_GetFieldIndex(hFeature, nms[i]), szName );
			}
		}
*/		
        hPt = OGR_G_CreateGeometry(wkb);
        OGR_G_SetPoint_2D(hPt, 0, x, y);

        OGR_F_SetGeometry( hFeature, hPt );
        OGR_G_DestroyGeometry(hPt);

        if( OGR_L_CreateFeature( hLayer, hFeature ) != OGRERR_NONE ) {
			setError("Failed to create feature");
			return false;
        }
        OGR_F_Destroy( hFeature );
    }
    GDALClose( hDS );
	return true;
}

/*
bool SpatVector::write(std::string filename, std::string format, bool overwrite) {

	msg.success = true;

    const char *pszDriverName = "ESRI Shapefile";
    GDALDriver *poDriver;
    //GDALAllRegister();
    poDriver = GetGDALDriverManager()->GetDriverByName(pszDriverName );
    if( poDriver == NULL ) {
        setError("driver not available");
		return false;
    }
    GDALDataset *poDS;
    poDS = poDriver->Create(filename.c_str(), 0, 0, 0, GDT_Unknown, NULL );
    if( poDS == NULL ) {
        setError("cannot write file");
		return false;
    }
    OGRLayer *poLayer;
    poLayer = poDS->CreateLayer(basename_noext(filename).c_str(), NULL, wkbPoint, NULL );
    if( poLayer == NULL ) {
        setError("Layer creation failed");
		return false;
    }
    OGRFieldDefn oField( "Name", OFTString );
    oField.SetWidth(32);
    if( poLayer->CreateField( &oField ) != OGRERR_NONE ) {
        setError("Field creation failed");
		return false;
    }
    double x, y;
    char szName[33];
    while( !feof(stdin) && fscanf( stdin, "%lf,%lf,%32s", &x, &y, szName ) == 3 )  {
        OGRFeature *poFeature;
        poFeature = OGRFeature::CreateFeature( poLayer->GetLayerDefn() );
        poFeature->SetField( "Name", szName );
        OGRPoint pt;
        pt.setX( x );
        pt.setY( y );
        poFeature->SetGeometry( &pt );
        if( poLayer->CreateFeature( poFeature ) != OGRERR_NONE ) {
			setError("Failed to create feature");
			return false;
        }
        OGRFeature::DestroyFeature( poFeature );
    }
    GDALClose( poDS );
	return true;
}
*/

#endif
