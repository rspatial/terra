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

#include <vector>
#include "spatRaster.h"


void getSampleRowCol(std::vector<unsigned> &oldrow, std::vector<unsigned> &oldcol, unsigned nrows, unsigned ncols, unsigned snrow, unsigned sncol) {
	double rf = nrows / (double)(snrow);
	double cf = ncols / (double)(sncol);
	//double rstart = std::floor(0.5 * rf);
	//double cstart = std::floor(0.5 * cf);
	double rstart = 0.5 * rf;
	double cstart = 0.5 * cf;
	oldcol.reserve(sncol);
	for (size_t i =0; i<sncol; i++) {
        oldcol.push_back(i * cf + cstart);
	}
	oldrow.reserve(snrow);
	for (size_t i =0; i<snrow; i++) {
        oldrow.push_back(i * rf + rstart);
	}
}


std::vector<double> SpatRaster::readSample(unsigned src, unsigned srows, unsigned scols) {
	unsigned oldnc = ncell();
	unsigned nl = source[src].nlyr;
	std::vector<unsigned> oldcol, oldrow;
	std::vector<double>	out; 
	getSampleRowCol(oldrow, oldcol, nrow(), ncol(), srows, scols);
	
	out.reserve(srows*scols);
    for (size_t lyr=0; lyr<nl; lyr++) {
        unsigned old_offset = lyr * oldnc;
        for (size_t r=0; r<srows; r++) {
 			unsigned oldc = old_offset + oldrow[r] * ncol();
            for (size_t c=0; c<scols; c++) {
				unsigned oldcell = oldc + oldcol[c];
                out.push_back(source[src].values[oldcell]);
            }
        }
	}
	return out;
}


SpatRaster SpatRaster::sampleRegularRaster(unsigned size) {

	double f = std::min(1.0, sqrt(size / ncell()));
	unsigned nr = ceil(nrow() * f);
	unsigned nc = ceil(ncol() * f);
	if ((size >= ncell()) || ((nc == ncol()) && (nr == nrow()))) {
		return( *this );
	}
	SpatRaster out = geometry(nlyr(), true);
	out.source[0].nrow = nr;
	out.source[0].ncol = nc;
	
	if (!source[0].hasValues) return (out);

	std::vector<double> v;
	for (size_t src=0; src<nsrc(); src++) {
		if (source[src].memory) {
			v = readSample(src, nr, nc);
		//} else if (source[src].driver == "raster") {
		//	v = readSampleBinary(src, nr, nc);
		} else {
		    #ifdef useGDAL
			v = readGDALsample(src, nr, nc);
			#endif
		}
		if (hasError()) return out;
		out.source[0].values.insert(out.source[0].values.end(), v.begin(), v.end());
	}
	out.source[0].memory = true;
	out.source[0].hasValues = true;
	out.source[0].setRange();

	return out;
}

std::vector<std::vector<double>> SpatRaster::sampleRegularValues(unsigned size) {

	std::vector<std::vector<double>> out;
	if (!source[0].hasValues) return (out);
	
	unsigned nsize;
	unsigned nr = nrow();
	unsigned nc = ncol();
	if (size < ncell()) {
		double f = sqrt(size / ncell());
		nr = std::ceil(nrow() * f);
		nc = std::ceil(ncol() * f);
	}
	nsize = nc * nr;
	std::vector<double> v;
	if ((size >= ncell()) || ((nc == ncol()) && (nr == nrow()))) {
		v = getValues() ;
		if (hasError()) return out;
		for (size_t i=0; i<nlyr(); i++) {
			size_t offset = i * nsize;
			std::vector<double> vv(v.begin()+offset, v.begin()+offset+nsize);
			out.push_back(vv);
		}
		return out;
	}

	for (size_t src=0; src<nsrc(); src++) {
		if (source[src].memory) {
			v = readSample(src, nr, nc);
		//} else if (source[src].driver == "raster") {
		//	v = readSampleBinary(src, nr, nc);
		} else {
		    #ifdef useGDAL
			v = readGDALsample(src, nr, nc);
			#endif
		}
		if (hasError()) return out;
		for (size_t i=0; i<source[src].nlyr; i++) {
			size_t offset = i * nsize;
			std::vector<double> vv(v.begin()+offset, v.begin()+offset+nsize);
			out.push_back(vv);
		}
	}	
	return out;
}


/*
std::vector<std::vector<double>> SpatRaster::sampleRandom(unsigned size, unsigned seed) {
    srand(seed); 
	unsigned nc = ncell();
	std::vector<double> r(size);
    for(size_t i = 0; i<size; i++) {
        r[i] = rand() % nc; 
	}
	std::vector<std::vector<double>> d = extractCell(r);
    return d; 
}
*/

#include <random>
#include <unordered_set>


std::vector<double> sample_without_replacement(unsigned size, unsigned N, unsigned seed){
    // Sample "size" elements from [1, N] without replacement
    // see https://stackoverflow.com/questions/4461446/stl-way-to-add-a-constant-value-to-a-stdvector    
	
	size = std::max(unsigned(1), std::min(size, N)); // k <= N
	
    std::default_random_engine gen(seed);   
	std::uniform_int_distribution<> distribution(1, N);
    std::unordered_set<unsigned> samples;
    
    for (size_t r = N - size; r < N; ++r) {
        int v = distribution(gen) - 1;
        if (!samples.insert(v).second) samples.insert(r);
    } 
    std::vector<double> result(samples.begin(), samples.end());
    std::shuffle(result.begin(), result.end(), gen);    
	
    return result;
}


std::vector<double> sample_with_replacement(unsigned size, unsigned N, unsigned seed){
    // Sample "size" elements from [1, N] without replacement
    // see https://stackoverflow.com/questions/4461446/stl-way-to-add-a-constant-value-to-a-stdvector    
	
	size = std::max((unsigned)1, std::min(size, N)); // k <= N
	
    std::default_random_engine gen(seed);   
    std::vector<double> samples;
    std::vector<double> result;
	result.reserve(size);
    
	std::uniform_int_distribution<> distribution(0, N-1);
    for (size_t i=0; i<size; i++) {
        result.push_back( distribution(gen) );
    } 
	
    return result;
}


// need to consider the possibility that ncell() > max_int

// if size is large, use (shuffle(values))[1:size] instead

std::vector<std::vector<double>> SpatRaster::sampleRandomValues(unsigned size, bool replace, unsigned seed) {

	double nc = ncell();
	std::vector<double> dcells;

	if (replace) {
		if (size >= .6 * nc) {
			dcells.resize(nc);
			std::iota(std::begin(dcells), std::end(dcells), 0);
			std::default_random_engine gen(seed);  
			std::shuffle(dcells.begin(), dcells.end(), gen);    
			if (size < nc) {
				dcells.erase(dcells.begin()+size, dcells.end());	
			}
			
		} else {
			dcells = sample_without_replacement(size, nc, seed);
		}
	} else {
		dcells = sample_with_replacement(size, nc, seed);
	}
	
	std::vector<std::vector<double>> d = extractCell(dcells);
	return d; 
}


SpatRaster SpatRaster::sampleRandomRaster(unsigned size, bool replace, unsigned seed) {

	unsigned nsize;
	unsigned nr = nrow();
	unsigned nc = ncol();
	if (size < ncell()) {
		double f = sqrt(size / ncell());
		nr = std::ceil(nrow() * f);
		nc = std::ceil(ncol() * f);
	}
	SpatRaster out = geometry(nlyr(), true);
	out.source[0].nrow = nr;
	out.source[0].ncol = nc;
	if (!source[0].hasValues) return (out);
	
	nsize = nr * nc;
	std::vector<std::vector<double>> vv = sampleRandomValues(nsize, replace, seed);

	for (size_t i=0; i<vv.size(); i++) {
		out.source[0].values.insert(out.source[0].values.end(), vv[i].begin(), vv[i].end());
	}
	out.source[0].memory = true;
	out.source[0].hasValues = true;
	out.source[0].setRange();

	return out;
}



/*
SpatDataFrame sampleCells(unsigned size, std::string type, bool replace) {

	SpatDataFrame out;
	if (size >= ncell() & !replace) {
		out.setError("size >= ncell() & !replace");
		return out;
	}
	if (!source[0].hasValues) {
		out.setError("Raster has no values");
		return (out);
	}
	if (type == "Random") {
		
	} else if (type == "Regular") {
		
	} else { //type == "Stratified" 
	
	} // else "Cluster"
	return out;
}
*/

