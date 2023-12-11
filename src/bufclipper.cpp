
#include "Rcpp.h"
#include "spatVector.h"
#include "recycle.h"
#include "string_utils.h"
#include "clipper.h"
#include "clipper.core.h"
//#include "../../Utils/clipper.svg.utils.h"
//#include "../../Utils/CommonUtils.h"


bool bufferPolygon(SpatGeom &g, double &d, const Clipper2Lib::JoinType &jtype, const Clipper2Lib::EndType &etype, const double &miter_limit, const int &precision, const double &arc_tolerance) {
	size_t np = g.size();
	Clipper2Lib::PathsD input;
	input.reserve(np);
	Clipper2Lib::PointD pnt;
	for (size_t j=0; j < np; j++) {
		SpatPart p = g.getPart(j);
		Clipper2Lib::PathD part;
		part.reserve(p.x.size());
		for (size_t k=0; k<p.x.size(); k++) {
			pnt.x = p.x[k];
			pnt.y = p.y[k];
			part.push_back(pnt);	
		}
		input.push_back(part);
	}
	Clipper2Lib::PathsD solution = InflatePaths(input, d, jtype, etype, miter_limit, precision, arc_tolerance);
			
	np = solution.size();
	SpatGeom gg;
	gg.gtype = polygons;
	for (size_t j=0; j<np; j++) {
		size_t nn = solution[j].size();
		if (nn > 0) {
			std::vector<double> X, Y;
			X.reserve(nn+1);
			Y.reserve(nn+1);
			for (size_t k=0; k<nn; k++) {
				X.push_back(solution[j][k].x);
				Y.push_back(solution[j][k].y);
			}
			X.push_back(X[0]);
			Y.push_back(Y[0]);
			SpatPart pp(X, Y);
			gg.addPart(pp);
		}
	}
	g = gg;
	return (g.parts.size() > 0);
}

SpatVector SpatVector::buffer4(std::vector<double> d, std::string jointype, double miter_limit, int precision, double arc_tolerance) {
//double miter_limit = 2.0, int precision = 2, double arc_tolerance = 0.0) {

	SpatVector out;
	size_t n = size();
	recycle(d, n);
	std::string vt = type();
	SpatExtent e = getExtent();
	if (precision < 0) {
		double m = std::max(e.xmax - e.xmin, e.ymax - e.ymin);
		m = std::round(sqrt(m));
		precision = std::max(2, (int) (8 - m));
	}
	lowercase(jointype);
	Clipper2Lib::JoinType jtype;
	if (jointype == "square") {
		jtype = Clipper2Lib::JoinType::Square;
	} else if (jointype == "bevel") {
		jtype = Clipper2Lib::JoinType::Bevel;
	} else if (jointype == "miter") {
		jtype = Clipper2Lib::JoinType::Miter;
	} else { // round
		jtype = Clipper2Lib::JoinType::Round;
	}

	if (vt == "points") {
//		Clipper2Lib::EndType etype = Clipper2Lib::EndType::Round;
		return out;
	} else if (vt == "lines") {
//		Clipper2Lib::EndType etype = Clipper2Lib::EndType::Joined;
		return out;
	} else { // polygons
		Clipper2Lib::EndType etype = Clipper2Lib::EndType::Polygon;
		for (size_t i=0; i<n; i++) {
			SpatGeom g = getGeom(i);
			if (bufferPolygon(g, d[i], jtype, etype, miter_limit, precision, arc_tolerance)) {
				out.addGeom(g);
			}
		}
	}
   
	out.srs = srs;
	return out;
	
}


	  