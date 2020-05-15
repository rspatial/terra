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

union_cascated
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



SpatVector SpatVector::aggregate(std::string field, bool dissolve) {

	SpatVector out;
	
	int i = where_in_vector(field, get_names());
	if (i < 0) {
		out.setError("cannot find field");
		return out;		
	}
	SpatDataFrame uv;
	std::vector<int> idx = lyr.df.getIndex(i, uv);

	for (size_t i=0; i<uv.nrow(); i++) {
		SpatGeom g;
		g.gtype = lyr.geoms[0].gtype;
		for (size_t j=0; j<idx.size(); j++) {
			if (i == (size_t)idx[j]) {
				g.unite( getGeom(j) );
			}	
		}
		out.addGeom(g);
	}
	
	out.lyr.srs = lyr.srs;
	out.lyr.df  = uv; 

	if (dissolve) {
		out.addWarning("cannot dissolve yet");
	}
	
	return out;
}

