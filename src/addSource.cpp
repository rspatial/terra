using namespace std;
#include <vector>
#include "spat.h"

SpatRaster SpatRaster::addSources(SpatRaster x) {
	SpatRaster out = deepCopy();
	out.error = false;
	out.warning = false;
	if (compare_geom(x, false, false)) {
        if (!hasValues()) {
            out.source = x.source;
        } else {
            out.source.insert(out.source.end(), x.source.begin(), x.source.end());
        }
	} else {
		out.error = true;
		out.error_message = "dimensions and/or extent do not match";
	}
	return(out);
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

std::vector<RasterSource> RasterSource::subset(std::vector<unsigned> lyrs) {
    std::vector<RasterSource> out;
    unsigned nl = lyrs.size();
    bool all = true;
    if (lyrs.size() == nl) {
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
        rs.names.resize(0);
        rs.hasRange.resize(0);
        rs.range_min.resize(0);
        rs.range_max.resize(0);

        if (memory) {
            rs.values.resize(0);
            if (hasValues) {
                for (size_t i=0; i<nl; i++) {
                    unsigned j = lyrs[i];
                    std::vector<double> x = getValues(j);
                    rs.values.insert(rs.values.end(), x.begin(), x.end());
                    rs.names.push_back(names[j]);
                    rs.hasRange.push_back(hasRange[j]);
                    rs.range_min.push_back(range_min[j]);
                    rs.range_max.push_back(range_max[j]);
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
            }
        }
        rs.nlyr = nl;
        out = { rs };
    }
    return out;
}


SpatRaster SpatRaster::subset(std::vector<unsigned> lyrs, string filename, bool overwrite) {

    SpatRaster out = geometry();
    out.source.resize(0);

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

    if (filename != "") {
        out.writeRaster(filename, overwrite);
    }
    return out;
}



