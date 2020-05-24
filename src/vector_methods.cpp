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
#include "string_utils.h"

#include "gdal_alg.h"
#include "ogrsf_frmts.h"


std::vector<bool> SpatVector::is_valid() {
	std::vector<bool> out;
	out.reserve(nrow());
	GDALDataset* src = write_ogr("", "layer", "Memory", true);
	OGRLayer *inLayer = src->GetLayer(0);
	inLayer->ResetReading();
	OGRFeature *inFeature;
	while( (inFeature = inLayer->GetNextFeature()) != NULL ) {
		OGRGeometry *poGeometry = inFeature->GetGeometryRef();
		out.push_back(poGeometry->IsValid());
		OGRFeature::DestroyFeature( inFeature );
	}
	return out;
}


SpatVector SpatVector::make_valid() {
	SpatVector out;
	GDALDataset* src = write_ogr("", "layer", "Memory", true);
	OGRLayer *inLayer = src->GetLayer(0);
	inLayer->ResetReading();
	OGRFeature *inFeature;
	while( (inFeature = inLayer->GetNextFeature()) != NULL ) {
		OGRGeometry *poGeometry = inFeature->GetGeometryRef();
		//OGRGeometry *poGeom = poGeometry->MakeValid();
		if (inFeature->SetGeometry( poGeometry ) != OGRERR_NONE) {
			out.setError("cannot set geometry");
			return out;
		}
		if (inLayer->SetFeature( inFeature ) != OGRERR_NONE) {
			out.setError("cannot set feature");
			return out;
		}
		OGRFeature::DestroyFeature( inFeature );
	}
	out.read_ogr(src);
	GDALClose(src);
	return out;
}




SpatVector SpatVector::disaggregate() {
	SpatVector out;
	out.srs = srs;
	out.df = df.skeleton();

	if (nrow() == 0) {
		return out;
	}
	
	for (size_t i=0; i<nrow(); i++) {
		SpatGeom g = getGeom(i);
		SpatDataFrame row = df.subset_rows(i);
		for (size_t j=0; j<g.parts.size(); j++) {
			SpatGeom gg = SpatGeom(g.parts[j]);
			gg.gtype = g.gtype;
			out.addGeom(gg);
			if (!out.df.rbind(row)) { 
				out.setError("cannot add row");
				return out;
			}
		}
	}
	
	return out;

}



SpatVector SpatVector::aggregate(std::string field, bool dissolve) {

	SpatVector out;
	
	int i = where_in_vector(field, get_names());
	if (i < 0) {
		out.setError("cannot find field");
		return out;		
	}
	SpatDataFrame uv;
	std::vector<int> idx = df.getIndex(i, uv);

	out.srs = srs;
	out.df  = uv; 

	if (!dissolve) {
		for (size_t i=0; i<uv.nrow(); i++) {
			SpatGeom g;
			g.gtype = geoms[0].gtype;
			for (size_t j=0; j<idx.size(); j++) {
				if (i == (size_t)idx[j]) {
					g.unite( getGeom(j) );
				}	
			}
			out.addGeom(g);
		}
	} else {

		GDALDataset* src = out.write_ogr("", "layer", "Memory", true);
		OGRLayer *inLayer = src->GetLayer(0);
		inLayer->ResetReading();
		OGRFeature *inFeature;
/*
		const OGRSpatialReference *srs = src->GetSpatialRef();
		GDALDataset *dst = NULL;
		GDALDriver *poDriver = GetGDALDriverManager()->GetDriverByName( "Memory" );
		dst = poDriver->Create("", 0, 0, 0, GDT_Unknown, NULL );
		OGRLayer *outLayer;
		outLayer = dst->CreateLayer("lyrname", (OGRSpatialReference *)srs, wkbMultiPolygon, NULL );
		OGRFieldDefn oField("uid", OFTInteger64);
		if( outLayer->CreateField( &oField ) != OGRERR_NONE ) {
			out.setError( "Creating union field failed");
			return out;
		}
		OGRFeature *outFeature;
		i = 0;
*/

		while( (inFeature = inLayer->GetNextFeature()) != NULL ) {
			OGRGeometry *poGeometry = inFeature->GetGeometryRef();
			//OGRMultiPolygon *poGeom = ( OGRMultiPolygon * )poGeometry;
			if (!poGeometry->IsValid()) {
				out.setError("invalid geom");
				return out;				
			}
			OGRGeometry *poGeom = poGeometry->UnionCascaded();	
			if (poGeom == NULL) {
				out.setError("union failed");
				return out;
			}
			//outFeature = OGRFeature::CreateFeature( outLayer->GetLayerDefn() );
			//outFeature->SetField(0, (GIntBig)i); i++;
			if (inFeature->SetGeometry( poGeom ) != OGRERR_NONE) {
				out.setError("cannot set geometry");
				return out;
			}
			if (inLayer->SetFeature( inFeature ) != OGRERR_NONE) {
				out.setError("cannot set feature");
				return out;
			}
			//OGRFeature::DestroyFeature( outFeature );
			OGRFeature::DestroyFeature( inFeature );
		}
		out.read_ogr(src);
		GDALClose(src);
//		GDALClose(dst);
	}
	
	return out;
}


/*
std::vector<OGRGeometry *> geoms_from_ds(GDALDataset* src, int field, int value) {
	std::vector<OGRGeometry *> g;
	OGRLayer *poLayer = src->GetLayer(0);
	poLayer->ResetReading();
	OGRFeature *poFeature;
	
	while( (poFeature = poLayer->GetNextFeature()) != NULL ) {
		OGRGeometry *poGeometry = poFeature->GetGeometryRef();
		g.push_back(poGeometry);
	}
	return g;
}
// create output dataset 
	GDALDataset* dst;
// get unique values in field
// loop over unique values
	// for value in uvalues
	std::vector<OGRGeometry *> gvec = geoms_from_ds(src, field, value);
	OGRGeometry *geom;
	geom = (OGRGeometry *) gvec.data();
	OGRGeometry *gout;
	gout = geom->UnionCascaded();
// set geometry to output	
   return dst;
*/

