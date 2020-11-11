// Copyright (c) 2018-2020  Robert J. Hijmans
//
// This file is part of the "spat" library
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
#include "vecmath.h"





SpatRaster SpatRaster::rectify(std::string method, SpatRaster aoi, unsigned useaoi, bool snap, SpatOptions &opt) {
	SpatRaster out = geometry(0);

	if (nsrc() > 1) {
		out.setError("you can transform only one data source at a time");
		return(out);
	}
	if (!source[0].rotated) {
		out.setError("this source is not rotated");
		return(out);
	} 
	GDALDataset *poDataset;
	std::string fname = source[0].filename;
	poDataset = (GDALDataset *) GDALOpen(fname.c_str(), GA_ReadOnly );
	if( poDataset == NULL )  {
		setError("cannot read from " + fname );
		return out;
	}
	double gt[6];	
	if( poDataset->GetGeoTransform(gt) != CE_None ) {
		out.setError("can't get geotransform");
		GDALClose( (GDALDatasetH) poDataset );
		return out;
	}
	GDALClose( (GDALDatasetH) poDataset );
	//SpatExtent e = getExtent();
	//std::vector<double> x = {e.xmin, e.xmin, e.xmax, e.xmax };
	//std::vector<double> y = {e.ymin, e.ymax, e.ymin, e.ymax };
	double nc = ncol();
	double nr = nrow();
	std::vector<double> x = {0, 0, nc, nc};
	std::vector<double> y = {0, nr, 0, nr};
	std::vector<double> xx(4);
	std::vector<double> yy(4);
	for (size_t i=0; i<4; i++) {
		xx[i] = gt[0] + x[i]*gt[1] + y[i]*gt[2];
		yy[i] = gt[3] + x[i]*gt[4] + y[i]*gt[5];
	}
	double xmin = vmin(xx, TRUE);
	double xmax = vmax(xx, TRUE);
	double ymin = vmin(yy, TRUE);
	double ymax = vmax(yy, TRUE);
	SpatExtent en(xmin, xmax, ymin, ymax);
	out = out.setResolution(gt[1], -gt[5]);
	out.setExtent(en, false, "out");

	if (useaoi == 1) { // use extent
		en = aoi.getExtent();
		if (snap) {
			en = out.align(en, "near");
			out.setExtent(en, false, "near");
		} else {
			out.setExtent(en, false, "");
		}
	} else if (useaoi == 2){  // extent and resolution
		out = aoi.geometry(0);
	} // else { // if (useaoi == 0) // no aoi
	
		
	out = warper(out, "", method, false, opt);

	return(out);
}
	
	
