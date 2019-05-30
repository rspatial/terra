// Copyright (c) 2018-2019  Robert J. Hijmans
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

#ifdef useGDAL


#include "spatRaster.h"
#include "string_utils.h"
#include "ogrsf_frmts.h"


bool SpatVector::write(std::string filename, std::string format, bool overwrite) {

	msg.success = true;

    const char *pszDriverName = "ESRI Shapefile";
    GDALDriver *poDriver;
    GDALAllRegister();
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
    poLayer = poDS->CreateLayer(basename(filename).c_str(), NULL, wkbPoint, NULL );
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


#endif