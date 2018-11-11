// Author: Robert J. Hijmans
// Date : November 2009
// Version 0.9
// Licence GPL v3

#include <vector>
#include "spatraster.h"

std::vector<double> SpatRaster::readSample(unsigned src, unsigned srows, unsigned scols) {

	double oldrow, oldcol;
	size_t oldcell, newcell, oldnc, newnc;
	std::vector<double>	out(srows * scols);
	double rf = srows / nrow;
	double cf = scols / ncol;
	oldnc = ncell();
	newnc = srows * scols;
	size_t nl = source[src].nlyr;
	
	for (size_t r=0; r<srows; r++) {
		oldrow = r + rf;
		for (size_t c=0; c<scols; c++) {
			oldcol = c + cf;	
			oldcell = oldrow * nrow + oldcol;
			newcell = r * scols + c;
			for (size_t lyr=0; lyr<nl; lyr++) {
				size_t old_offset = lyr * oldnc;
				size_t new_offset = lyr * newnc;
				out[new_offset + newcell] = source[src].values[old_offset + oldcell];
			}
		}
	}
	return out;	
}


SpatRaster SpatRaster::sampleRegular(unsigned size) {

	if (size >= ncell()) return( *this );
	
	double f = size / ncell();
	unsigned nr = ceil(nrow * f); 
	unsigned nc = ceil(ncol * f); 
	if ((nc == ncol) && (nr == nrow)) return( *this );

	SpatRaster out = geometry(); 
	out.source[0].nrow=nr;
	out.source[0].ncol=nc;

	if (!source[0].hasValues) return (out);

	std::vector<double> v;	
	for (size_t src=0; src<nsrc(); src++) {
		if (source[src].memory) {
			v = readSample(src, nr, nc);			
		} else {
			v = readGDALsample(src, nr, nc);
		}
		out.source[0].values.insert(out.source[0].values.end(), v.begin(), v.end());
	}
	return out;
}

