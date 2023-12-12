#include "Rcpp.h"
#include "spatVector.h"
#include "recycle.h"
#include "string_utils.h"
#include "clipper.h"
#include "clipper.core.h"


Clipper2Lib::PathD makePath(std::vector<double> &x, std::vector<double> &y) {
	Clipper2Lib::PathD part;
	part.reserve(x.size());
	Clipper2Lib::PointD pnt;
	for (size_t i=0; i<x.size(); i++) {
		pnt.x = x[i];
		pnt.y = y[i];
		part.push_back(pnt);	
	}
	return part;
}


SpatGeom getPath(Clipper2Lib::PathsD solution, SpatGeomType gtype) {
	size_t np = solution.size();
	SpatGeom gg;
	gg.gtype = gtype;
	for (size_t i=0; i<np; i++) {
		size_t nn = solution[i].size();
		if (nn > 0) {
			std::vector<double> X, Y;
			X.reserve(nn+1);
			Y.reserve(nn+1);
			for (size_t j=0; j<nn; j++) {
				X.push_back(solution[i][j].x);
				Y.push_back(solution[i][j].y);
			}
			X.push_back(X[0]);
			Y.push_back(Y[0]);
			SpatPart pp(X, Y);
			gg.addPart(pp);
		}
	}
	return gg;
}


bool bufferLines(SpatGeom &g, double &d, const Clipper2Lib::JoinType &jtype, const double &miter_limit, const int &precision, const double &arc_tolerance) {
	
	Clipper2Lib::EndType etype = Clipper2Lib::EndType::Joined;
	size_t np = g.size();
	Clipper2Lib::PathsD input;
	input.reserve(np);
	for (size_t j=0; j < np; j++) {
		SpatPart p = g.getPart(j);
		Clipper2Lib::PathD part = makePath(p.x, p.y);
		input.push_back(part);
	}
	Clipper2Lib::PathsD solution = InflatePaths(input, d, jtype, etype, miter_limit, precision, arc_tolerance);
	g = getPath(solution, lines);
	return (g.parts.size() > 0);
}


SpatGeom bufferPolygon(SpatGeom &g, double &d, const Clipper2Lib::JoinType &jtype, const double &miter_limit, const int &precision, const double &arc_tolerance) {

	Clipper2Lib::EndType etype = Clipper2Lib::EndType::Polygon;
	size_t np = g.size();
	Clipper2Lib::PathsD input;
	
	input.reserve(np);
	bool has_holes = false;
	for (size_t i=0; i<np; i++) {
		SpatPart p = g.getPart(i);
		has_holes = has_holes || (p.holes.size() > 0);
		Clipper2Lib::PathD part = makePath(p.x, p.y);
		input.push_back(part);
	}
	SpatGeom out;
	if (has_holes) {
		SpatGeom gg;
		for (size_t i=0; i<input.size(); i++) {
			Clipper2Lib::PathsD inp {input[i]};
			Clipper2Lib::PathsD sol = InflatePaths(inp, d, jtype, etype, miter_limit, precision, arc_tolerance);
			gg = getPath(sol, polygons);
			if (gg.parts.size() == 0) continue;
	
			size_t nh = g.parts[i].holes.size();
			if (nh > 0) {
				Clipper2Lib::PathsD hinp;
				for (size_t j=0; j<nh; j++) {
					Clipper2Lib::PathD part = makePath(g.parts[i].holes[j].x, 
														g.parts[i].holes[j].y);
					hinp.push_back(part);
				}
				Clipper2Lib::PathsD hsol=InflatePaths(hinp, -d, jtype, etype, miter_limit, precision, arc_tolerance);
				SpatGeom gh = getPath(hsol, polygons);
				for (size_t j=0; j<gh.parts.size(); j++) {
					if (d < 0) {
						SpatVector v(gg);
						SpatVector h(gh);
						v = v.erase(h);
						if (v.empty()) continue;
						h = h.crop(v);
						if (h.empty()) continue;
						gg = v.geoms[0];
						gh = h.geoms[0];
					}
					gg.parts[0].addHole(gh.parts[j].x, 
										gh.parts[j].y);	
				}
				SpatVector a(gg);
				a.aggregate(true);
				out.addPart(a.geoms[0].parts[0]);
			} else {
				out.addPart(gg.parts[0]);
			}
		}
	} else {
		Clipper2Lib::PathsD solution = InflatePaths(input, d, jtype, etype, miter_limit, precision, arc_tolerance);
		out = getPath(solution, polygons);
	}
	return out;
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
		for (size_t i=0; i<n; i++) {
			SpatGeom g = getGeom(i);
			g = bufferPolygon(g, d[i], jtype, miter_limit, precision, arc_tolerance);	
			if (g.parts.size() > 0) {
				out.addGeom(g);
			}
		}
	}
   
	out.srs = srs;
	return out;
	
}


	  