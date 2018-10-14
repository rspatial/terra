using namespace std;
#include <vector>
#include "spat.h"


SpatRaster SpatRaster::addSource(SpatRaster x) {
	SpatRaster out = deepCopy();
	out.error = false;
	out.warning = false;
	if (compare_geom(x, false, false)) {
		out.source.insert(out.source.end(), x.source.begin(), x.source.end());
	} else {
		out.error = true;
		out.error_message = "dimensions and/or extent do not match";
	}
	return(out);
}


unsigned SpatRaster::sourceFromLyr(unsigned lyr) {
    if (lyr < 0 || lyr >= nlyr()) {
        return(NAN);
    }
    unsigned nsrc = 0;
    unsigned nlyrs = 0;
    for (size_t i=0; i<source.size(); i++) {
        if (nlyrs >= lyr) break;
        nlyrs += source[i].nlyr;
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




/*
RasterSource SpatRaster::subsetSource(RasterSource s, std::vector<unsigned> lyrs) {
    RasterSource out;

}

std::vector<unsigned> SpatRaster::sourcesFromLyrs(std::vector<unsigned> lyrs) {
    std::vector<unsigned> s(lyrs.size());
    std::vector<unsigned> slyrs = lyrsBySource();
    unsigned nl = nlyr();

    for (size_t i=0; i<lyrs.size(); i++) {
        if (lyrs[i] < 0 || lyrs[i] >= nl) {
            s[i] = NAN;
        } else {
            s[i] = slyrs[i];
        }
    }
    return s;
}


SpatRaster SpatRaster::subset(std::vector<unsigned> i) {
    SpatRaster out;
    unsigned nl = nlyr();
    bool first = true;
    RasterSource rs;
    unsigned lyrcnt = 0;
    unsigned lastsrc = 0;

    for (size_t j=0; j<i.size(); j++) {
        unsigned lyr = i[j];


            size_t s = sourceFromLyr(lyr);




        if (lyr == lyrcnt) {
            lyrnct++;
        } else {

            size_t s = sourceFromLyr(lyr);
            if (s > lastsrc) {
                for (size_t k=lastsrc; k<s; k++) {
                    rs = source[k];
                    if (first) {
                        out.setSource(rs)
                        first = false;
                    } else {
                        out.source.push_back( source[lyr] );
                    }
                }
            } else {



            if (lyr >= 0 && lyr < nl) {
                if (source[lyr].nlyr == 1) {
                    rs = source[lyr];
                } else {
                    if ()
                    rs = ;
                }

                if (first) {
                    out.setSource(rs)
                    first = false;
                } else {
                    out.source.push_back( source[lyr] );
                }
            } else {
                out.warning = true;
                out.warning_message.push_back("index " + std::to_string(lyr) + " is out of range");
            }
        }
    }
    return(out);
}

*/


