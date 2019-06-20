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

#include <vector>
#include "spatRaster.h"

/*
#include "string_utils.h"

void RasterSource::fsopen(std::string filename) {
    std::string grifile = setFileExt(filename, ".gri");
 	std::ofstream fstr(grifile, std::ios::out | std::ios::binary);
    *ofs = &fstr;
}

bool RasterSource::fswrite(std::vector<double> &v) {
	unsigned sz = v.size() * sizeof(double);
	bool result = (*ofs).write(reinterpret_cast<const char*>(&v[0]), sz);
 	return result;
}

void RasterSource::fsclose() {
	(*ofs).close();
}
*/


RasterSource::RasterSource() {
	open_write = false;
	open_read = false;
}


SpatRaster SpatRaster::combineSources(SpatRaster x) {
	SpatRaster out = deepCopy();
	if (compare_geom(x, false, false)) {
        if (!hasValues()) {
            out.source = x.source;
        } else {
            out.source.insert(out.source.end(), x.source.begin(), x.source.end());
        }
	} else {
		out.setError("dimensions and/or extent do not match");
	}
	return(out);
}

void SpatRaster::addSource(SpatRaster x) {
	if (compare_geom(x, false, false)) {
        if (!hasValues()) {  //or if n src == 0?
            source = x.source;
        } else {
            source.insert(source.end(), x.source.begin(), x.source.end());
        }
	} else {
		setError("dimensions and/or extent do not match");
	}
}


unsigned SpatRaster::nsrc() {
    return source.size();
}


int SpatRaster::sourceFromLyr(unsigned lyr) {
    if (lyr >= nlyr()) {
        return(-1);
    }
    unsigned nsrc = 0;
    unsigned nlyrs = -1;
    for (size_t i=0; i<source.size(); i++) {
        nlyrs += source[i].nlyr;
        if (nlyrs >= lyr) break;
        nsrc++;
    }
    return nsrc;
}


std::vector<unsigned> SpatRaster::nlyrBySource() {
    std::vector<unsigned> lyrs(source.size());
    for (size_t i=0; i<source.size(); i++) {
        lyrs[i] = source[i].nlyr;
    }
    return lyrs;
}

std::vector<unsigned> SpatRaster::lyrsBySource() {
    std::vector<unsigned> lyrs(nlyr());
    unsigned start = 0;
    for (size_t i=0; i<source.size(); i++) {
        unsigned n = source[i].nlyr;
        unsigned stop = start + n;
        for (size_t j=start; j<stop; j++) {
            lyrs[j] = i;
        }
        start = stop;
    }
    return lyrs;
}


std::vector<unsigned> SpatRaster::findLyr(unsigned lyr) {
    std::vector<unsigned> sl(2);
    unsigned nlyrs = 0;
    unsigned start = 0;
	bool done = false;
    for (size_t i=0; i<source.size(); i++) {
		if ((nlyrs + source[i].nlyr) >= lyr) {	
			sl[0] = i;
			for (size_t j=start; j<source[i].nlyr; j++) {
				if ((nlyrs + j) == lyr) {
					sl[1] = j;
					done = true;
					break;
				}
			}
		}
		if (done) break;
        nlyrs += source[i].nlyr;
    }
    return sl;
}




std::vector<unsigned> SpatRaster::sourcesFromLyrs(std::vector<unsigned> lyrs) {
    std::vector<unsigned> s(lyrs.size());
    std::vector<unsigned> slyrs = lyrsBySource();
    for (size_t i=0; i<lyrs.size(); i++) {
        s[i] = slyrs[ lyrs[i] ];
    }
    return s;
}



std::vector<double> RasterSource::getValues(unsigned lyr) {
    unsigned nc = nrow * ncol;
    unsigned start = lyr * nc;
    unsigned stop = start + nc;
    std::vector<double> out (values.begin()+start, values.begin()+stop);
    return out;
}


void RasterSource::resize(unsigned n) {
	names.resize(n);
    hasRange.resize(n);
    range_min.resize(n);
    range_max.resize(n);
	has_scale_offset.resize(n);
	scale.resize(n);
	offset.resize(n);
    hasColors.resize(n);
	cols.resize(n);
    hasCategories.resize(n);
	cats.resize(n);
	nlyr = n;
	layers.resize(n);
}


std::vector<RasterSource> RasterSource::subset(std::vector<unsigned> lyrs) {
	std::vector<RasterSource> out;

    unsigned nl = lyrs.size();
    bool all = true;
    if (lyrs.size() == nlyr) {
        for (size_t i=0; i<nl; i++) {
            if (lyrs[i] != i) {
                all = false;
                break;
            }
        }
    } else {
        all = false;
    }
    if (all) {
        out = { *this };
    } else {
        RasterSource rs = *this;
        rs.resize(0);

        if (memory) {
            if (hasValues) {
                rs.values.resize(0);
                for (size_t i=0; i<nl; i++) {
                    unsigned j = lyrs[i];
                    std::vector<double> x = getValues(j);
                    rs.values.insert(rs.values.end(), x.begin(), x.end());
                    rs.names.push_back(names[j]);
                    rs.hasRange.push_back(hasRange[j]);
                    rs.range_min.push_back(range_min[j]);
                    rs.range_max.push_back(range_max[j]);
                    rs.hasColors.push_back(hasColors[j]);
                    rs.hasCategories.push_back(hasCategories[j]);
					//rs.RAT.push_back(RAT[j]);
					rs.has_scale_offset.push_back(has_scale_offset[j]);
					rs.scale.push_back(scale[j]);
					rs.offset.push_back(offset[j]);
                }
            }
        } else {
            rs.layers = lyrs;
            for (size_t i=0; i<nl; i++) {
                unsigned j = lyrs[i];
                rs.names.push_back(names[j]);
                rs.hasRange.push_back(hasRange[j]);
                rs.range_min.push_back(range_min[j]);
                rs.range_max.push_back(range_max[j]);
                rs.hasColors.push_back(hasColors[j]);
                rs.hasCategories.push_back(hasCategories[j]);
				rs.has_scale_offset.push_back(has_scale_offset[j]);
				rs.scale.push_back(scale[j]);
				rs.offset.push_back(offset[j]);
            }
        }
        rs.nlyr = nl;
        out = { rs };
    }
    return out;
}

std::vector<unsigned> validLayers( std::vector<unsigned> lyrs , unsigned nl) {
    unsigned s = lyrs.size();
    for (size_t i=0; i<s; i++) {
        unsigned j = s - i - 1; // start from the back
        if ((lyrs[j] < 0) | (lyrs[j] >= nl)) {
			lyrs.erase(lyrs.begin() + j);
		}
	}

	/* or
    unsigned s = lyrs.size() - 1;
    for (long i=s; i>=0; i--) {
        if ((lyrs[i] < 0) | (lyrs[i] >= nl)) {
			lyrs.erase(lyrs.begin() + i);
		}
	}
	*/

	return lyrs;


}


SpatRaster SpatRaster::subset(std::vector<unsigned> lyrs, SpatOptions &opt) {

    SpatRaster out = geometry();
    out.source.resize(0);

    unsigned oldsize = lyrs.size();
    lyrs = validLayers(lyrs, nlyr());

	if (lyrs.size() == 0) {
		out.setError("no (valid) layer references");
		return(out);
	} else if (lyrs.size() != oldsize) {
        out.addWarning("ignored " + std::to_string(oldsize - lyrs.size()) + " invalid layer reference(s)");
	}

    std::vector<unsigned> srcs = sourcesFromLyrs(lyrs);

    unsigned ss = srcs[0];
    std::vector<unsigned> slyr;
    std::vector<unsigned> lyrbys = nlyrBySource();
    RasterSource rs;
    std::vector<RasterSource> vrs;
    unsigned offset = 0;
    for (size_t i=0; i<ss; i++) { offset += lyrbys[i]; }

    for (size_t i=0; i<lyrs.size(); i++) {
        if (srcs[i] == ss) {
            slyr.push_back( (lyrs[i] - offset) );
        } else {
            vrs = source[ss].subset(slyr);
            out.source.insert(out.source.end(), vrs.begin(), vrs.end());
            ss = srcs[i];
            offset = 0;
            for (size_t i=0; i<ss; i++) { offset += lyrbys[i]; }
            slyr = { lyrs[i] - offset } ;
       }
    }

    vrs = source[ss].subset(slyr);
    out.source.insert(out.source.end(), vrs.begin(), vrs.end());

    if (opt.get_filename() != "") {
        out.writeRaster(opt);
    }
    return out;
}



