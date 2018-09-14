using namespace std;
#include "spat.h"
#include "util.h"

bool SpatRaster::compare(unsigned nrows, unsigned ncols, SpatExtent e ) {
	double xrange = extent.xmax - extent.xmin;
	double yrange = extent.ymax - extent.ymin;
	bool e1 = is_equal_range(e.xmax, extent.xmax, xrange, 0.00001);
	bool e2 = is_equal_range(e.xmin, extent.xmin, xrange, 0.00001);
	bool e3 = is_equal_range(e.ymax, extent.ymax, yrange, 0.00001);
	bool e4 = is_equal_range(e.ymin, extent.ymin, yrange, 0.00001);
	bool eOK = (e1 && e2 && e3 && e4);
	bool rcOK = (nrow == nrows) && (ncol == ncols);
	return (rcOK && eOK);
}


