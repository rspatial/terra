// Based on  public-domain code by Darel Rex Finley, 2007
// http://alienryderflex.com/polygon_fill/

#include <vector>
#include "spatraster.h"

std::vector<double> rasterize_polygon(std::vector<double> r, double value, const std::vector<double> &pX, const std::vector<double> &pY, const unsigned nrows, const unsigned ncols, const double xmin, const double ymax, const double rx, const double ry) {

	unsigned n = pX.size();
	std::vector<unsigned> nCol(n);

	for (size_t row=0; row<nrows; row++) {

		double y = ymax - (row+0.5) * ry;

		//  Build a list of nodes.
		unsigned nodes = 0;
		size_t j = n-1;
		for (size_t i=0; i<n; i++) {
			if (((pY[i] < y) && (pY[j] >= y)) || ((pY[j] < y) && (pY[i] >= y))) {
			//	nCol[nodes++]=(int)  (((pX[i] - xmin + (y-pY[i])/(pY[j]-pY[i]) * (pX[j]-pX[i])) + 0.5 * rx ) / rx);
				double nds = ((pX[i] - xmin + (y-pY[i])/(pY[j]-pY[i]) * (pX[j]-pX[i])) + 0.5 * rx ) / rx;
				nds = nds < 0 ? 0 : nds;
		        nds = nds > ncols ? ncols : nds;
				nCol[nodes] = (unsigned) nds;
				nodes++;
			}
			j = i;
		}

		std::sort(nCol.begin(), nCol.begin()+nodes);
		unsigned ncell = ncols * row;

		//  Fill the cells between node pairs.
		for (size_t i=0; i < nodes; i+=2) {
			if (nCol[i+1] > 0 && nCol[i] < ncols) {
				for (size_t col = nCol[i]; col < nCol[i+1]; col++) {
					r[col + ncell] = value;
				}
			}
		}
	}
	return(r);
}



SpatRaster SpatRaster::rasterizePolygons(SpatLayer p, double background, std::string filename, bool overwrite) {

	SpatRaster out = geometry();
  	out.writeStart(filename, overwrite);
	double value = 1;
	double resx = xres();
	double resy = yres();

	SpatGeom poly;
	SpatPart part;
	SpatHole hole;

	unsigned n = p.size();


	for (size_t i = 0; i < out.bs.n; i++) {
		std::vector<double> v(out.bs.nrows[i] * ncol, background);

		for (size_t j = 0; j < n; j++) {

			poly = p.getGeom(j);
			//value = p.getAtt(j);
			//if (std::isnan(value)) { value = j;}
			value = j+1;

			unsigned np = poly.size();

			for (size_t k = 0; k < np; k++) {
				part = poly.getPart(k);
				if (part.hasHoles()) {
					std::vector<double> vv = rasterize_polygon(v, value, part.x, part.y, nrow, ncol, extent.xmin, extent.ymax, resx, resy);
					for (size_t h=0; h < part.nHoles(); h++) {
						hole = part.getHole(h);
						vv = rasterize_polygon(vv, background, hole.x, hole.y, nrow, ncol, extent.xmin, extent.ymax, resx, resy);
					}
					for (size_t q=0; q < vv.size(); q++) {
						if (vv[q] != background) {
							v[q] = vv[q];
						}
					}
				} else {
					v = rasterize_polygon(v, value, part.x, part.y, nrow, ncol, extent.xmin, extent.ymax, resx, resy);
				}
			}
		}
		out.writeValues(v, out.bs.row[i]);
	}
	out.writeStop();

	return(out);

}
