// Copyright (c) 2018-2022  Robert J. Hijmans
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



SpatRaster SpatRaster::combineSources(SpatRaster x, bool warn) {

	SpatRaster out = geometry();
	if (!hasValues()) {
		if (!x.hasValues()) {
			if (out.compare_geom(x, false, false, 0.1)) {
				out.source.insert(out.source.end(), x.source.begin(), x.source.end());
				out.setNames(out.getNames());
			} else {
				out = x.deepCopy();
				if (warn) {
					out.addWarning("both rasters were empty, but had different geometries. The first one was ignored");
				}
			}
		} else {
			out = x.deepCopy();
			if (warn) {
				out.addWarning("the first raster was empty and ignored");
			}
		}
		return out;
	}

	if (!out.compare_geom(x, false, false, 0.1)) {
		return out;
	}

	out = deepCopy();
	if (!x.hasValues()) {
		out.addWarning("you cannot add SpatRaster with no values to one that has values");
		return(out);
	}
	out.checkTime(x);
	out.source.insert(out.source.end(), x.source.begin(), x.source.end());
    // to make names unique (not great if called several times
	//out.setNames(out.getNames());
	return(out);
}


void SpatRaster::combine(SpatRaster x) {

	if (!compare_geom(x, false, false, 0.1)) {
		return;
	}

	bool hv = hasValues();
	if (hv != x.hasValues()) {
		setError("combined sources must all have values; or none should have values");
		return;
	}

	checkTime(x);
	source.insert(source.end(), x.source.begin(), x.source.end());
	//setNames(getNames());
	return;
}


void SpatRaster::checkTime(SpatRaster &x) {
	if (!hasTime()) {
		std::vector<int_64> time;
		x.setTime(time, "remove", "");
		return;
	}
	if (!x.hasTime()) {
		std::vector<int_64> time;
		setTime(time, "remove", "");
		return;
	}
	std::string s = source[0].timestep;
	std::string xs = x.source[0].timestep;
	if (s == xs) return;
	if ((s == "days") && (xs == "seconds")) {
		x.source[0].timestep = "days";
	} else if ((s == "seconds") && (xs == "days")) {
		for (size_t i=0; i<source.size(); i++) {
			source[i].timestep = "days";
		}
	} else {
		std::vector<int_64> time;
		setTime(time, "remove", "");
		x.setTime(time, "remove", "");
	}
}

void SpatRaster::addSource(SpatRaster x, bool warn, SpatOptions &opt) {

	if (!hasValues()) {
		if (!x.hasValues()) {
			if (compare_geom(x, false, false, 0.1)) {
				source.insert(source.end(), x.source.begin(), x.source.end());
			} else {
				source = x.source;
				if (warn) {
					addWarning("both rasters were empty, but had different geometries. The first one was ignored");
				}
			}
		} else {
			source = x.source;
			if (warn) {
				addWarning("the first raster was empty and was ignored");
			}
		}
		return;
	}
	if (compare_geom(x, false, false, 0.1)) {
		if (!x.hasValues()) {
			x = x.init({NAN}, opt);
		}
		checkTime(x);
        source.insert(source.end(), x.source.begin(), x.source.end());
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

std::vector<unsigned> SpatRaster::getBands() {
	std::vector<unsigned> out;
    for (size_t i=0; i<source.size(); i++) {
		out.insert(out.end(), source[i].layers.begin(), source[i].layers.end());
	}
	return out;
}






std::vector<unsigned> SpatRaster::sourcesFromLyrs(std::vector<unsigned> lyrs) {
    std::vector<unsigned> s(lyrs.size());
    std::vector<unsigned> slyrs = lyrsBySource();
    for (size_t i=0; i<lyrs.size(); i++) {
        s[i] = slyrs[ lyrs[i] ];
    }
    return s;
}


/*
void SpatRasterSource::getValues(std::vector<double> &v, unsigned lyr) {
	size_t nc ;
	if (hasWindow) {
		nc = window.full_ncol * window.full_nrow;
	} else {
		nc = nrow * ncol;
	}
	size_t start = lyr * nc;
	v = std::vector<double>(values.begin()+start, values.begin()+start+nc);
}
*/


void SpatRasterSource::appendValues(std::vector<double> &v, unsigned lyr) {
	size_t nc ;
	if (hasWindow) {
		nc = window.full_ncol * window.full_nrow;
	} else {
		nc = nrow * ncol;
	}
	size_t start = lyr * nc;
	v.insert(v.end(), values.begin()+start, values.begin()+start+nc);
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
	valueType.resize(n, 0);
	dataType.resize(n, "");
    hasRange.resize(n, false);
    range_min.resize(n, NAN);
    range_max.resize(n, NAN);
    blockcols.resize(n);
    blockrows.resize(n);

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


void SpatRasterSource::reserve(unsigned n) {
	names.reserve(n);
	time.reserve(n);
	unit.reserve(n);
	depth.reserve(n);
	valueType.reserve(n);
	dataType.reserve(n);
    hasRange.reserve(n);
    range_min.reserve(n);
    range_max.reserve(n);
    blockcols.reserve(n);
    blockrows.reserve(n);

	has_scale_offset.reserve(n);
	scale.reserve(n);
	offset.reserve(n);
    hasColors.reserve(n);
	cols.reserve(n);
    hasCategories.reserve(n);
	cats.reserve(n);
    //hasAttributes.reserve(n);
	//atts.reserve(n);
    //attsIndex.reserve(n);
	nlyr = n;
	layers.reserve(n);
}


//std::vector<SpatRasterSource> SpatRasterSource::subset(std::vector<unsigned> lyrs) {
SpatRasterSource SpatRasterSource::subset(std::vector<unsigned> lyrs) {

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
		return *this ;
	}

	SpatRasterSource out;
	if (memory) {
		out.srs = srs;
		out.ncol = ncol;
		out.nrow = nrow;
		out.nlyrfile = nlyrfile;
		out.extent = extent;
		out.rotated = rotated;
		out.flipped = flipped;
		out.hasWindow = hasWindow;
		out.window = window;
		out.source_name = source_name;
		out.source_name_long = source_name_long;
		out.timestep = timestep;
		out.hasTime = hasTime;
		out.hasUnit = hasUnit;
		out.memory = memory;
		out.hasValues = hasValues;
		//out.filename = filename;
		//out.driver = driver;
		//out.valueType = valueType;
		//out.hasNAflag = hasNAflag;
		//out.NAflag = NAflag;
	} else {
		//no values, deep copy is cheap
		out = *this;
		out.resize(0);
	}

	out.reserve(nl);
    for (size_t i=0; i<nl; i++) {
        unsigned j = lyrs[i];
		out.names.push_back(names[j]);
		out.time.push_back(time[j]);
		out.depth.push_back(depth[j]);
		out.unit.push_back(unit[j]);
		out.valueType.push_back(valueType[j]);
		out.dataType.push_back(dataType[j]);
        out.hasRange.push_back(hasRange[j]);
        out.range_min.push_back(range_min[j]);
        out.range_max.push_back(range_max[j]);
        out.blockrows.push_back(blockrows[j]);
        out.blockcols.push_back(blockcols[j]);
        out.hasColors.push_back(hasColors[j]);
        out.cols.push_back(cols[j]);
        out.hasCategories.push_back(hasCategories[j]);
        out.cats.push_back(cats[j]);
		out.has_scale_offset.push_back(has_scale_offset[j]);
		out.scale.push_back(scale[j]);
		out.offset.push_back(offset[j]);

		if (memory) {
			out.layers.push_back(i);
			if (hasValues) {
				appendValues(out.values, j);
			}
		} else {
			out.layers.push_back(layers[j]);
		}
    }
    out.nlyr = nl;
	out.hasValues = hasValues;
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

    SpatRaster out = geometry(1);
    out.source.resize(0);

    unsigned oldsize = lyrs.size();
    lyrs = validLayers(lyrs, nlyr());

	if (lyrs.size() == 0) {
		out.setError("no (valid) layer selected");
		return(out);
	} else if (lyrs.size() != oldsize) {
        out.addWarning("ignored " + std::to_string(oldsize - lyrs.size()) + " invalid layer reference(s)");
	}

    std::vector<unsigned> srcs = sourcesFromLyrs(lyrs);
    unsigned ss = srcs[0];
    std::vector<unsigned> slyr;
    std::vector<unsigned> lyrbys = nlyrBySource();
//    SpatRasterSource rs;
    unsigned offset = 0;
    for (size_t i=0; i<ss; i++) { offset += lyrbys[i]; }

    for (size_t i=0; i<lyrs.size(); i++) {
        if (srcs[i] == ss) {
            slyr.push_back( (lyrs[i] - offset) );
        } else {
            out.source.push_back( source[ss].subset(slyr) );
            ss = srcs[i];
            offset = 0;
            for (size_t j=0; j<ss; j++) { offset += lyrbys[j]; }
            slyr = { lyrs[i] - offset } ;
		}
    }

    out.source.push_back( source[ss].subset(slyr) );
    if (opt.get_filename() != "") {
        out = out.writeRaster(opt);
    } //else {
	//	out.collapse();
	//}

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
	valueType.insert(valueType.end(), x.valueType.begin(), x.valueType.end());
	dataType.insert(dataType.end(), x.dataType.begin(), x.dataType.end());
	hasRange.insert(hasRange.end(), x.hasRange.begin(), x.hasRange.end());
	range_min.insert(range_min.end(), x.range_min.begin(), x.range_min.end());
	range_max.insert(range_max.end(), x.range_max.begin(), x.range_max.end());

	blockrows.insert(blockrows.end(), x.blockrows.begin(), x.blockrows.end());
	blockcols.insert(blockcols.end(), x.blockcols.begin(), x.blockcols.end());
	//hasAttributes.insert(hasAttributes.end(), x.hasAttributes.begin(), x.hasAttributes.end());
	//atts.insert(atts.end(), x.atts.begin(), x.atts.end());
	//attsIndex.insert(attsIndex.end(), x.attsIndex.begin(), x.attsIndex.end());
	hasCategories.insert(hasCategories.end(), x.hasCategories.begin(), x.hasCategories.end());
	cats.insert(cats.end(), x.cats.begin(), x.cats.end());
	hasColors.insert(hasColors.end(), x.hasColors.begin(), x.hasColors.end());
	cols.insert(cols.end(), x.cols.begin(), x.cols.end());
	valueType.insert(valueType.end(), x.valueType.begin(), x.valueType.end());
	dataType.insert(dataType.end(), x.dataType.begin(), x.dataType.end());
	has_scale_offset.insert(has_scale_offset.end(), x.has_scale_offset.begin(), x.has_scale_offset.end());
	scale.insert(scale.end(), x.scale.begin(), x.scale.end());
	offset.insert(offset.end(), x.offset.begin(), x.offset.end());
	return true;
}



bool SpatRasterSource::combine(SpatRasterSource &x) {
	if (memory & x.memory) {
		if ((values.size() + x.values.size()) < (values.max_size()/8) ) {
			values.insert(values.end(), x.values.begin(), x.values.end());
			layers.resize(nlyr + x.nlyr);
			std::iota(layers.begin(), layers.end(), 0);
			x.values.resize(0);
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
	dataType.insert(dataType.end(), x.dataType.begin(), x.dataType.end());
	valueType.insert(valueType.end(), x.valueType.begin(), x.valueType.end());
	hasRange.insert(hasRange.end(), x.hasRange.begin(), x.hasRange.end());
	range_min.insert(range_min.end(), x.range_min.begin(), x.range_min.end());
	range_max.insert(range_max.end(), x.range_max.begin(), x.range_max.end());
	blockrows.insert(blockrows.end(), x.blockrows.begin(), x.blockrows.end());
	blockcols.insert(blockcols.end(), x.blockcols.begin(), x.blockcols.end());
	//hasAttributes.insert(hasAttributes.end(), x.hasAttributes.begin(), x.hasAttributes.end());
	//atts.insert(atts.end(), x.atts.begin(), x.atts.end());
	//attsIndex.insert(attsIndex.end(), x.attsIndex.begin(), x.attsIndex.end());
	hasCategories.insert(hasCategories.end(), x.hasCategories.begin(), x.hasCategories.end());
	cats.insert(cats.end(), x.cats.begin(), x.cats.end());
	hasColors.insert(hasColors.end(), x.hasColors.begin(), x.hasColors.end());
	cols.insert(cols.end(), x.cols.begin(), x.cols.end());
	valueType.insert(valueType.end(), x.valueType.begin(), x.valueType.end());
	dataType.insert(dataType.end(), x.dataType.begin(), x.dataType.end());
	has_scale_offset.insert(has_scale_offset.end(), x.has_scale_offset.begin(), x.has_scale_offset.end());
	scale.insert(scale.end(), x.scale.begin(), x.scale.end());
	offset.insert(offset.end(), x.offset.begin(), x.offset.end());
	return true;
}



void SpatRaster::collapse() {
	size_t n = nsrc();
	if (n < 2) return;
	std::vector<size_t> rem;
	for (size_t i=1; i<n; i++) {
		if (source[0].combine(source[i])) {
			rem.push_back(i);
		}
	}
	for (int i=rem.size(); i>= 0; i--) {
		source.erase(source.begin()+i);
	}
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


