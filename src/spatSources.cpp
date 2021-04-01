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

#include <vector>
#include "spatRaster.h"

/*
#include "string_utils.h"

void SpatRasterSource::fsopen(std::string filename) {
    std::string grifile = setFileExt(filename, ".gri");
 	std::ofstream fstr(grifile, std::ios::out | std::ios::binary);
    *ofs = &fstr;
}

bool SpatRasterSource::fswrite(std::vector<double> &v) {
	unsigned sz = v.size() * sizeof(double);
	bool result = (*ofs).write(reinterpret_cast<const char*>(&v[0]), sz);
 	return result;
}

void SpatRasterSource::fsclose() {
	(*ofs).close();
}
*/


SpatRasterSource::SpatRasterSource() {
	open_write = false;
	open_read = false;
}



SpatRaster SpatRaster::combineSources(SpatRaster x) {

	SpatRaster out = geometry();

	if (!out.compare_geom(x, false, false)) {
		return out;
	}

	bool hv = hasValues();
	if (hv != x.hasValues()) {
		out.setError("combined sources must all have values; or none should have values");
		return(out);
	}

	out = deepCopy();
//    if (!hv) {
//       out.source = x.source;
//    } else {
	out.source.insert(out.source.end(), x.source.begin(), x.source.end());
//	}
    // to make names unique
	out.setNames(out.getNames());
	return(out);
}


void SpatRaster::addSource(SpatRaster x) {

	if (compare_geom(x, false, false)) {
        if (!hasValues()) {  //or if n src == 0?
            source = x.source;
        } else {
            source.insert(source.end(), x.source.begin(), x.source.end());
        }
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



std::vector<double> SpatRasterSource::getValues(unsigned lyr) {
	size_t nc ;
	if (hasWindow) {
		nc = window.full_ncol * window.full_nrow;
	} else {
		nc = nrow * ncol;
	}
	size_t start = lyr * nc;
	std::vector<double> out(values.begin()+start, values.begin()+start+nc);
	return out;

}

bool SpatRasterSource::in_order() {
	if (memory) return true;
	if (nlyr != nlyrfile) return false;
	for (size_t i=0; i<layers.size(); i++) {
		if (layers[i] != i) {
			return false;
		}
	} 
	return true;
}


void SpatRasterSource::resize(unsigned n) {
	names.resize(n, "");
	time.resize(n);
	unit.resize(n);
	depth.resize(n);
    hasRange.resize(n, false);
    range_min.resize(n, NAN);
    range_max.resize(n, NAN);
	has_scale_offset.resize(n, false);
	scale.resize(n, 1);
	offset.resize(n, 0);
    hasColors.resize(n, false);
	cols.resize(n);
    hasCategories.resize(n, false);
	cats.resize(n);
    //hasAttributes.resize(n);
	//atts.resize(n);
    //attsIndex.resize(n);
	nlyr = n;
	layers.resize(n);
	std::iota(layers.begin(), layers.end(), 0);
}


//std::vector<SpatRasterSource> SpatRasterSource::subset(std::vector<unsigned> lyrs) {
SpatRasterSource SpatRasterSource::subset(std::vector<unsigned> lyrs) {

	SpatRasterSource out;

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
	out = *this ;
    if (!all) {
        out.resize(0);
        if (memory) {
            if (hasValues) {
                out.values.resize(0);
                for (size_t i=0; i<nl; i++) {
                    unsigned j = lyrs[i];
                    std::vector<double> x = getValues(j);
                    out.values.insert(out.values.end(), x.begin(), x.end());
                    out.names.push_back(names[j]);
					out.time.push_back(time[j]);
					out.depth.push_back(depth[j]);
					out.unit.push_back(unit[j]);
					out.layers.push_back(i);
                    out.hasRange.push_back(hasRange[j]);
                    out.range_min.push_back(range_min[j]);
                    out.range_max.push_back(range_max[j]);
                    out.hasColors.push_back(hasColors[j]);
                    out.cols.push_back(cols[j]);
                    out.hasCategories.push_back(hasCategories[j]);
                    out.cats.push_back(cats[j]);
					//out.hasAttributes.push_back(hasAttributes[j]);
					//out.atts.push_back(atts[j]);
					//out.attsIndex.push_back(attsIndex[j]);
					out.has_scale_offset.push_back(has_scale_offset[j]);
					out.scale.push_back(scale[j]);
					out.offset.push_back(offset[j]);
                }
            }
        } else {
            //out.layers = lyrs;
            for (size_t i=0; i<nl; i++) {
                unsigned j = lyrs[i];
                out.layers.push_back(layers[j]);
                out.names.push_back(names[j]);
				out.time.push_back(time[j]);
				out.unit.push_back(unit[j]);
				out.depth.push_back(depth[j]);
                out.hasRange.push_back(hasRange[j]);
                out.range_min.push_back(range_min[j]);
                out.range_max.push_back(range_max[j]);
                out.hasColors.push_back(hasColors[j]);
                out.cols.push_back(cols[j]);
                out.hasCategories.push_back(hasCategories[j]);
                out.cats.push_back(cats[j]);
				//out.hasAttributes.push_back(hasAttributes[j]);
				//out.atts.push_back(atts[j]);
				//out.attsIndex.push_back(attsIndex[j]);
				out.has_scale_offset.push_back(has_scale_offset[j]);
				out.scale.push_back(scale[j]);
				out.offset.push_back(offset[j]);
            }
        }
        out.nlyr = nl;
    }
    return out;
}

std::vector<unsigned> validLayers( std::vector<unsigned> lyrs , unsigned nl) {
    unsigned s = lyrs.size();
    for (size_t i=0; i<s; i++) {
        unsigned j = s - i - 1; // start from the back
        if (lyrs[j] >= nl) {
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
    SpatRasterSource rs;
    unsigned offset = 0;
    for (size_t i=0; i<ss; i++) { offset += lyrbys[i]; }

    for (size_t i=0; i<lyrs.size(); i++) {
        if (srcs[i] == ss) {
            slyr.push_back( (lyrs[i] - offset) );
        } else {
            rs = source[ss].subset(slyr);
            out.source.push_back(rs);
            ss = srcs[i];
            offset = 0;
            for (size_t i=0; i<ss; i++) { offset += lyrbys[i]; }
            slyr = { lyrs[i] - offset } ;
		}
    }

    rs = source[ss].subset(slyr);
    out.source.push_back(rs);

    if (opt.get_filename() != "") {
        out.writeRaster(opt);
    } else {
		out = out.collapse_sources();
	}
	
    return out;
}




bool SpatRasterSource::combine_sources(const SpatRasterSource &x) {
	if (memory & x.memory) {
		if ((values.size() + x.values.size()) < (values.max_size()/8) ) {
			values.insert(values.end(), x.values.begin(), x.values.end());
			layers.resize(nlyr + x.nlyr);
			std::iota(layers.begin(), layers.end(), 0);
		} else {
			return false;
		}
	} else if (filename == x.filename) {
		layers.insert(layers.end(), x.layers.begin(), x.layers.end());
	} else {
		return false;
	}
	nlyr += x.nlyr;
	names.insert(names.end(), x.names.begin(), x.names.end());
	time.insert(time.end(), x.time.begin(), x.time.end());
	if (!(hasTime & x.hasTime)) {
		hasTime = false;
	}
	unit.insert(unit.end(), x.unit.begin(), x.unit.end());

	depth.insert(depth.end(), x.depth.begin(), x.depth.end());
	hasRange.insert(hasRange.end(), x.hasRange.begin(), x.hasRange.end());
	range_min.insert(range_min.end(), x.range_min.begin(), x.range_min.end());
	range_max.insert(range_max.end(), x.range_max.begin(), x.range_max.end());
	//hasAttributes.insert(hasAttributes.end(), x.hasAttributes.begin(), x.hasAttributes.end());
	//atts.insert(atts.end(), x.atts.begin(), x.atts.end());
	//attsIndex.insert(attsIndex.end(), x.attsIndex.begin(), x.attsIndex.end());
	hasCategories.insert(hasCategories.end(), x.hasCategories.begin(), x.hasCategories.end());
	cats.insert(cats.end(), x.cats.begin(), x.cats.end());
	hasColors.insert(hasColors.end(), x.hasColors.begin(), x.hasColors.end());
	cols.insert(cols.end(), x.cols.begin(), x.cols.end());
	datatype.insert(datatype.end(), x.datatype.begin(), x.datatype.end());
	has_scale_offset.insert(has_scale_offset.end(), x.has_scale_offset.begin(), x.has_scale_offset.end());
	scale.insert(scale.end(), x.scale.begin(), x.scale.end());
	offset.insert(offset.end(), x.offset.begin(), x.offset.end());
	return true;
}



SpatRaster SpatRaster::collapse_sources() {
	SpatRaster out;
	std::vector<SpatRasterSource> src;
	SpatRasterSource s = source[0];
	for (size_t i=1; i<nsrc(); i++) {
		if (! s.combine_sources(source[i])) {
			src.push_back(s);
			s = source[i];
		}
	}
	src.push_back(s);
	out.setSources(src);
	return out;
}


