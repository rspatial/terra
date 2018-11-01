// Author: Robert J. Hijmans
// Date : November 2009
// Version 0.9
// Licence GPL v3

#include <vector>
#include "spatraster.h"
using namespace std;


// std::vector<double> SpatRaster::sampleRegular(unsigned size, bool cells, bool asRaster) {


//	stopifnot(hasValues(x) | isTRUE(xy))
//  stopifnot(size > 0)

/*
	bool allx = false;
	if (size >= rcut.ncell()) {
		if (asRaster) {
			SpatRaster r = *this;
			return(r);
		} 
		unsigned n = ncell();
		std::vector<unsigned> cells(n);
		for (size_t i=0; i<n; i++){ cell[i] = i; }
		
	} else {
		double f = nrow/(nrow+ncol) 
		double nrd = f * size;
		double ncd = (1-f) * size;
		f = sqrt(size/(nr*nc));
		nrd = nrd * f;
		ncd = ncd * f;
		unsigned nr = round(nrd);
		unsigned nc = round(ncd);
		std::vector<unsigned> rows(nr);
		std::vector<unsigned> cols(nc);
		nrd = nrd + 1;
		ncd = ncd + 1;
		for (size_t i=0; i<nr; i++){ rows[i] = round((i+1) * nrow/(nrd)); }
		for (size_t i=0; i<nc; i++){ cols[i] = round((i+1) * ncol/(ncd)); }		
		std::vector<unsigned> cells = cellFromRowCol(rows, cols);		
	}
	

	
	if (asRaster) {			
			if (allx) {
				if (!is.null(ext)) {
					return(crop(x, ext))
				} else {
					return(x)
				}
			} 
			
			
			cell = cellFromRowCol(x, rep(rows, each=nc), rep(cols, times=nr))
			if (hv) {
				m = .cellValues(x, cell)
			} else {
				m = NA
			}

			if (is.null(ext))  {
				outras = raster(x)
			} else {
				outras = raster(ext) 
				crs(outras) = crs(x)
			}
			nrow(outras) = nr
			ncol(outras) = nc
			
//		outras = brick(outras, nlyr=nlyr)
		
		outras = setValues(outras, m)
		names(outras) = names(x)
		if (any(is.factor(x))) {
			levels(outras) = levels(x)
		}
		return(outras)
		
	} else {
		
		if (allx) {
			cell <= 1:ncell(rcut)
		} else {
			cell = cellFromRowCol(x, rep(rows, each=nc), rep(cols, times=nr))
		}
		m = NULL
		nstart = 1
		if (xy) {
			m = xyFromCell(x, cell)
			nstart = 3
		}
		if (cells) {
			m = cbind(m, cell=cell)
			nstart = nstart + 1
		} 
		if (hv) {
			m = cbind(m, .cellValues(x, cell))
			colnames(m)[nstart:(nstart+nlyr-1)] = names(x)
		} 
			
		return(m)
	}	
}

*/