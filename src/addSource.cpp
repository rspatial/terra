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

