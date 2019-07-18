#include "spatVector.h"

#ifdef useGEOS

#include <geos/geom/PrecisionModel.h>
#include <geos/geom/GeometryFactory.h>
#include <geos/geom/Geometry.h>
#include <geos/geom/Point.h>
#include <geos/geom/LinearRing.h>
#include <geos/geom/LineString.h>
#include <geos/geom/MultiLineString.h>
#include <geos/geom/Polygon.h>
#include <geos/geom/MultiPolygon.h>
#include <geos/geom/GeometryCollection.h>
#include <geos/geom/Coordinate.h>
#include <geos/geom/CoordinateSequence.h>
#include <geos/geom/CoordinateArraySequence.h>
#include <geos/geom/IntersectionMatrix.h>
//#include <geos/io/WKBReader.h>
//#include <geos/io/WKBWriter.h>
//#include <geos/io/WKTWriter.h>
#include <geos/util/GeometricShapeFactory.h>
#include <geos/geom/util/SineStarFactory.h>
#include <geos/util/GEOSException.h>
#include <geos/util/IllegalArgumentException.h>
#include <geos/opLinemerge.h>
#include <geos/opPolygonize.h>
//#include <geos/constants.h>
//#include <vector>
//#include <sstream>
//#include <iomanip>
//#include <cstdlib> // exit()


using namespace geos;
using namespace geos::geom;
using namespace geos::operation::polygonize;
using namespace geos::operation::linemerge;
using geos::util::GEOSException;
using geos::util::IllegalArgumentException;


////////////////////////////////////////////////////////////////////////
// UNARY OPERATIONS
////////////////////////////////////////////////////////////////////////

std::vector<Geometry*>* gCentroid(std::vector<Geometry*>* geoms) {
	std::vector<Geometry*>* newgeoms = new std::vector<Geometry*>;
    for(unsigned int i = 0; i < geoms->size(); i++) {
        Geometry* g = (*geoms)[i];
        newgeoms->push_back(g->getCentroid()); //.release());
    }
	return newgeoms;
}

std::vector<Geometry*>* gBuffer(std::vector<Geometry*>* geoms, double d, int segments, int capstyle, bool& success) {
	std::vector<Geometry*>* newgeoms = new std::vector<Geometry*>;
    for(unsigned int i = 0; i < geoms->size(); i++) {
        Geometry* g = (*geoms)[i];
        try {
            Geometry* g2 = g->buffer(d, segments, capstyle);//.release();
            newgeoms->push_back(g2);
        }
        catch(const GEOSException& exc) {
			success = false;
			return newgeoms;
        }
    }
	return newgeoms;
}


std::vector<Geometry*>* gConvexHull(std::vector<Geometry*>* geoms) {
	std::vector<Geometry*>* newgeoms = new std::vector<Geometry*>;
    for(unsigned int i = 0; i < geoms->size(); i++) {
        Geometry* g = (*geoms)[i];
        newgeoms->push_back(g->convexHull()); //.release());
    }
	return newgeoms;
}


void kill_geoms(std::vector<Geometry*>* geoms) {
    for(unsigned int i = 0; i < geoms->size(); i++) {
        delete(*geoms)[i];
    }
    delete geoms;	
}


////////////////////////////////////////////////////////////////////////
// RELATIONAL OPERATORS
////////////////////////////////////////////////////////////////////////

std::vector<bool> disjoint(std::vector<Geometry*>* geoms, bool& success) {
	size_t n = geoms->size();
	std::vector<bool> out;
    for(size_t i = 0; i < n; i++) {
        Geometry* g1 = (*geoms)[i];
        for(size_t j = 0; j < n; j++) {
            Geometry* g2 = (*geoms)[j];
            try {
                if (g1->disjoint(g2)) {
                    out.push_back(true);
                } else {
                    out.push_back(false);
                }
            }
            // Geometry Collection is not a valid argument
            catch(const IllegalArgumentException& exc) {
				success = false;
				return out;
            }
            catch(const std::exception& exc) {
				success = false;
				return out;
            }
        }
    }
	return out;
}




/*  /////////////////////////////////////////////
    // TOUCHES
    /////////////////////////////////////////////

    std::cout << std::endl;
    std::cout << "    TOUCHES   ";
    for(unsigned int i = 0; i < geoms->size(); i++) {
        std::cout << "\t[" << i << "]";
    }
    std::cout << std::endl;
    for(unsigned int i = 0; i < geoms->size(); i++) {
        Geometry* g1 = (*geoms)[i];
        std::cout << "      [" << i << "]\t";
        for(unsigned int j = 0; j < geoms->size(); j++) {
            Geometry* g2 = (*geoms)[j];
            try {
                if(g1->touches(g2)) {
                    std::cout << " 1\t";
                }
                else {
                    std::cout << " 0\t";
                }
            }
            // Geometry Collection is not a valid argument
            catch(const IllegalArgumentException& exc) {
                std::cout << " X\t";
            }
            catch(const std::exception& exc) {
                std::cerr << exc.what() << std::endl;
            }
        }
        std::cout << std::endl;
    }

    /////////////////////////////////////////////
    // INTERSECTS
    /////////////////////////////////////////////

    std::cout << std::endl;
    std::cout << " INTERSECTS   ";
    for(unsigned int i = 0; i < geoms->size(); i++) {
        std::cout << "\t[" << i << "]";
    }
    std::cout << std::endl;
    for(unsigned int i = 0; i < geoms->size(); i++) {
        Geometry* g1 = (*geoms)[i];
        std::cout << "      [" << i << "]\t";
        for(unsigned int j = 0; j < geoms->size(); j++) {
            Geometry* g2 = (*geoms)[j];
            try {
                if(g1->intersects(g2)) {
                    std::cout << " 1\t";
                }
                else {
                    std::cout << " 0\t";
                }
            }
            // Geometry Collection is not a valid argument
            catch(const IllegalArgumentException& exc) {
                std::cout << " X\t";
            }
            catch(const std::exception& exc) {
                std::cerr << exc.what() << std::endl;
            }
        }
        std::cout << std::endl;
    }

    /////////////////////////////////////////////
    // CROSSES
    /////////////////////////////////////////////

    std::cout << std::endl;
    std::cout << "    CROSSES   ";
    for(unsigned int i = 0; i < geoms->size(); i++) {
        std::cout << "\t[" << i << "]";
    }
    std::cout << std::endl;
    for(unsigned int i = 0; i < geoms->size(); i++) {
        Geometry* g1 = (*geoms)[i];
        std::cout << "      [" << i << "]\t";
        for(unsigned int j = 0; j < geoms->size(); j++) {
            Geometry* g2 = (*geoms)[j];
            try {
                if(g1->crosses(g2)) {
                    std::cout << " 1\t";
                }
                else {
                    std::cout << " 0\t";
                }
            }
            // Geometry Collection is not a valid argument
            catch(const IllegalArgumentException& exc) {
                std::cout << " X\t";
            }
            catch(const std::exception& exc) {
                std::cerr << exc.what() << std::endl;
            }
        }
        std::cout << std::endl;
    }

    /////////////////////////////////////////////
    // WITHIN
    /////////////////////////////////////////////

    std::cout << std::endl;
    std::cout << "     WITHIN   ";
    for(unsigned int i = 0; i < geoms->size(); i++) {
        std::cout << "\t[" << i << "]";
    }
    std::cout << std::endl;
    for(unsigned int i = 0; i < geoms->size(); i++) {
        Geometry* g1 = (*geoms)[i];
        std::cout << "      [" << i << "]\t";
        for(unsigned int j = 0; j < geoms->size(); j++) {
            Geometry* g2 = (*geoms)[j];
            try {
                if(g1->within(g2)) {
                    std::cout << " 1\t";
                }
                else {
                    std::cout << " 0\t";
                }
            }
            // Geometry Collection is not a valid argument
            catch(const IllegalArgumentException& exc) {
                std::cout << " X\t";
            }
            catch(const std::exception& exc) {
                std::cerr << exc.what() << std::endl;
            }
        }
        std::cout << std::endl;
    }

    /////////////////////////////////////////////
    // CONTAINS
    /////////////////////////////////////////////

    std::cout << std::endl;
    std::cout << "   CONTAINS   ";
    for(unsigned int i = 0; i < geoms->size(); i++) {
        std::cout << "\t[" << i << "]";
    }
    std::cout << std::endl;
    for(unsigned int i = 0; i < geoms->size(); i++) {
        Geometry* g1 = (*geoms)[i];
        std::cout << "      [" << i << "]\t";
        for(unsigned int j = 0; j < geoms->size(); j++) {
            Geometry* g2 = (*geoms)[j];
            try {
                if(g1->contains(g2)) {
                    std::cout << " 1\t";
                }
                else {
                    std::cout << " 0\t";
                }
            }
            // Geometry Collection is not a valid argument
            catch(const IllegalArgumentException& exc) {
                std::cout << " X\t";
            }
            catch(const std::exception& exc) {
                std::cerr << exc.what() << std::endl;
            }
        }
        std::cout << std::endl;
    }

    /////////////////////////////////////////////
    // OVERLAPS
    /////////////////////////////////////////////

    std::cout << std::endl;
    std::cout << "   OVERLAPS   ";
    for(unsigned int i = 0; i < geoms->size(); i++) {
        std::cout << "\t[" << i << "]";
    }
    std::cout << std::endl;
    for(unsigned int i = 0; i < geoms->size(); i++) {
        Geometry* g1 = (*geoms)[i];
        std::cout << "      [" << i << "]\t";
        for(unsigned int j = 0; j < geoms->size(); j++) {
            Geometry* g2 = (*geoms)[j];
            try {
                if(g1->overlaps(g2)) {
                    std::cout << " 1\t";
                }
                else {
                    std::cout << " 0\t";
                }
            }
            // Geometry Collection is not a valid argument
            catch(const IllegalArgumentException& exc) {
                std::cout << " X\t";
            }
            catch(const std::exception& exc) {
                std::cerr << exc.what() << std::endl;
            }
        }
        std::cout << std::endl;
    }

    /////////////////////////////////////////////
    // RELATE
    /////////////////////////////////////////////

    std::cout << std::endl;
    std::cout << "     RELATE   ";
    for(unsigned int i = 0; i < geoms->size(); i++) {
        std::cout << "\t[" << i << "]";
    }
    std::cout << std::endl;
    for(unsigned int i = 0; i < geoms->size(); i++) {
        Geometry* g1 = (*geoms)[i];
        std::cout << "      [" << i << "]\t";
        for(unsigned int j = 0; j < geoms->size(); j++) {
            Geometry* g2 = (*geoms)[j];
            try {
                // second argument is intersectionPattern
                std::string pattern = "212101212";
                if(g1->relate(g2, pattern)) {
                    std::cout << " 1\t";
                }
                else {
                    std::cout << " 0\t";
                }

                // get the intersectionMatrix itself
                auto im = g1->relate(g2);
            }
            // Geometry Collection is not a valid argument
            catch(const IllegalArgumentException& exc) {
                std::cout << " X\t";
            }
            catch(const std::exception& exc) {
                std::cerr << exc.what() << std::endl;
            }
        }
        std::cout << std::endl;
    }

    /////////////////////////////////////////////
    // EQUALS
    /////////////////////////////////////////////

    std::cout << std::endl;
    std::cout << "     EQUALS   ";
    for(unsigned int i = 0; i < geoms->size(); i++) {
        std::cout << "\t[" << i << "]";
    }
    std::cout << std::endl;
    for(unsigned int i = 0; i < geoms->size(); i++) {
        Geometry* g1 = (*geoms)[i];
        std::cout << "      [" << i << "]\t";
        for(unsigned int j = 0; j < geoms->size(); j++) {
            Geometry* g2 = (*geoms)[j];
            try {
                if(g1->equals(g2)) {
                    std::cout << " 1\t";
                }
                else {
                    std::cout << " 0\t";
                }
            }
            // Geometry Collection is not a valid argument
            catch(const IllegalArgumentException& exc) {
                std::cout << " X\t";
            }
            catch(const std::exception& exc) {
                std::cerr << exc.what() << std::endl;
            }
        }
        std::cout << std::endl;
    }

    /////////////////////////////////////////////
    // EQUALS_EXACT
    /////////////////////////////////////////////

    std::cout << std::endl;
    std::cout << "EQUALS_EXACT  ";
    for(unsigned int i = 0; i < geoms->size(); i++) {
        std::cout << "\t[" << i << "]";
    }
    std::cout << std::endl;
    for(unsigned int i = 0; i < geoms->size(); i++) {
        Geometry* g1 = (*geoms)[i];
        std::cout << "      [" << i << "]\t";
        for(unsigned int j = 0; j < geoms->size(); j++) {
            Geometry* g2 = (*geoms)[j];
            try {
                // second argument is a tolerance
                if(g1->equalsExact(g2, 0.5)) {
                    std::cout << " 1\t";
                }
                else {
                    std::cout << " 0\t";
                }
            }
            // Geometry Collection is not a valid argument
            catch(const IllegalArgumentException& exc) {
                std::cout << " X\t";
            }
            catch(const std::exception& exc) {
                std::cerr << exc.what() << std::endl;
            }
        }
        std::cout << std::endl;
    }

    /////////////////////////////////////////////
    // IS_WITHIN_DISTANCE
    /////////////////////////////////////////////

    std::cout << std::endl;
    std::cout << "IS_WITHIN_DIST";
    for(unsigned int i = 0; i < geoms->size(); i++) {
        std::cout << "\t[" << i << "]";
    }
    std::cout << std::endl;
    for(unsigned int i = 0; i < geoms->size(); i++) {
        Geometry* g1 = (*geoms)[i];
        std::cout << "      [" << i << "]\t";
        for(unsigned int j = 0; j < geoms->size(); j++) {
            Geometry* g2 = (*geoms)[j];
            try {
                // second argument is the distance
                if(g1->isWithinDistance(g2, 2)) {
                    std::cout << " 1\t";
                }
                else {
                    std::cout << " 0\t";
                }
            }
            // Geometry Collection is not a valid argument
            catch(const IllegalArgumentException& exc) {
                std::cout << " X\t";
            }
            catch(const std::exception& exc) {
                std::cerr << exc.what() << std::endl;
            }
        }
        std::cout << std::endl;
    }


////////////////////////////////////////////////////////////////////////
// COMBINATIONS
////////////////////////////////////////////////////////////////////////

    std::cout << std::endl;
    std::cout << "-------------------------------------------------------------------------------" << std::endl;
    std::cout << "COMBINATIONS" << std::endl;
    std::cout << "-------------------------------------------------------------------------------" << std::endl;

    /////////////////////////////////////////////
    // UNION
    /////////////////////////////////////////////

    // Make unions of all geoms
    newgeoms = new std::vector<Geometry*>;
    for(unsigned int i = 0; i < geoms->size() - 1; i++) {
        Geometry* g1 = (*geoms)[i];
        for(unsigned int j = i + 1; j < geoms->size(); j++) {
            Geometry* g2 = (*geoms)[j];
            try {
                Geometry* g3 = g1->Union(g2).release();
                newgeoms->push_back(g3);
            }
            // It's illegal to union a collection ...
            catch(const IllegalArgumentException& ill) {
                //std::cerr <<ill.toString()<<"\n";
            }
            catch(const std::exception& exc) {
                std::cerr << exc.what() << std::endl;
            }
        }
    }

    // Print all unions
    std::cout << std::endl << "----- AND HERE ARE SOME UNION COMBINATIONS ------" << std::endl;
    wkt_print_geoms(newgeoms);

    // Delete the resulting geoms
    for(unsigned int i = 0; i < newgeoms->size(); i++) {
        delete(*newgeoms)[i];
    }
    delete newgeoms;


    /////////////////////////////////////////////
    // INTERSECTION
    /////////////////////////////////////////////

    // Compute intersection of adhiacent geometries
    newgeoms = new std::vector<Geometry*>;
    for(unsigned int i = 0; i < geoms->size() - 1; i++) {
        Geometry* g1 = (*geoms)[i];
        for(unsigned int j = i + 1; j < geoms->size(); j++) {
            Geometry* g2 = (*geoms)[j];
            try {
                Geometry* g3 = g1->intersection(g2).release();
                newgeoms->push_back(g3);
            }
            // Collection are illegal as intersection argument
            catch(const IllegalArgumentException& ill) {
                //std::cerr <<ill.toString()<<"\n";
            }
            catch(const std::exception& exc) {
                std::cerr << exc.what() << std::endl;
            }
        }
    }

    std::cout << std::endl << "----- HERE ARE SOME INTERSECTIONS COMBINATIONS ------" << std::endl;
    wkt_print_geoms(newgeoms);

    // Delete the resulting geoms
    for(unsigned int i = 0; i < newgeoms->size(); i++) {
        delete(*newgeoms)[i];
    }
    delete newgeoms;

    /////////////////////////////////////////////
    // DIFFERENCE
    /////////////////////////////////////////////

    // Compute difference of adhiacent geometries
    newgeoms = new std::vector<Geometry*>;
    for(unsigned int i = 0; i < geoms->size() - 1; i++) {
        Geometry* g1 = (*geoms)[i];
        for(unsigned int j = i + 1; j < geoms->size(); j++) {
            Geometry* g2 = (*geoms)[j];
            try {
                Geometry* g3 = g1->difference(g2).release();
                newgeoms->push_back(g3);
            }
            // Collection are illegal as difference argument
            catch(const IllegalArgumentException& ill) {
                //std::cerr <<ill.toString()<<"\n";
            }
            catch(const std::exception& exc) {
                std::cerr << exc.what() << std::endl;
            }
        }
    }

    std::cout << std::endl << "----- HERE ARE SOME DIFFERENCE COMBINATIONS ------" << std::endl;
    wkt_print_geoms(newgeoms);

    // Delete the resulting geoms
    for(unsigned int i = 0; i < newgeoms->size(); i++) {
        delete(*newgeoms)[i];
    }
    delete newgeoms;

    /////////////////////////////////////////////
    // SYMMETRIC DIFFERENCE
    /////////////////////////////////////////////

    // Compute symmetric difference of adhiacent geometries
    newgeoms = new std::vector<Geometry*>;
    for(unsigned int i = 0; i < geoms->size() - 1; i++) {
        Geometry* g1 = (*geoms)[i];
        for(unsigned int j = i + 1; j < geoms->size(); j++) {
            Geometry* g2 = (*geoms)[j];
            try {
                Geometry* g3 = g1->symDifference(g2).release();
                newgeoms->push_back(g3);
            }
            // Collection are illegal as symdifference argument
            catch(const IllegalArgumentException& ill) {
                //std::cerr <<ill.toString()<<"\n";
            }
            catch(const std::exception& exc) {
                std::cerr << exc.what() << std::endl;
            }
        }
    }

    std::cout << std::endl << "----- HERE ARE SYMMETRIC DIFFERENCES ------" << std::endl;
    wkt_print_geoms(newgeoms);

    // Delete the resulting geoms
    for(unsigned int i = 0; i < newgeoms->size(); i++) {
        delete(*newgeoms)[i];
    }
    delete newgeoms;


    /////////////////////////////////////////////
    // LINEMERGE
    /////////////////////////////////////////////
    LineMerger lm;
    lm.add(geoms);
    std::vector<LineString*>* mls = lm.getMergedLineStrings();
    newgeoms = new std::vector<Geometry*>;
    for(unsigned int i = 0; i < mls->size(); i++) {
        newgeoms->push_back((*mls)[i]);
    }
    delete mls;

    std::cout << std::endl << "----- HERE IS THE LINEMERGE OUTPUT ------" << std::endl;
    wkt_print_geoms(newgeoms);

    // Delete the resulting geoms
    for(unsigned int i = 0; i < newgeoms->size(); i++) {
        delete(*newgeoms)[i];
    }
    delete newgeoms;


    /////////////////////////////////////////////
    // POLYGONIZE
    /////////////////////////////////////////////
    Polygonizer plgnzr;
    plgnzr.add(geoms);
    auto polys = plgnzr.getPolygons();
    newgeoms = new std::vector<Geometry*>;
    for(unsigned int i = 0; i < polys->size(); i++) {
        newgeoms->push_back((*polys)[i].release());
    }

    std::cout << std::endl << "----- HERE IS POLYGONIZE OUTPUT ------" << std::endl;
    wkt_print_geoms(newgeoms);

    // Delete the resulting geoms
    for(unsigned int i = 0; i < newgeoms->size(); i++) {
        delete(*newgeoms)[i];
    }
    delete newgeoms;


    /////////////////////////////////////////////
    // CLEANUP
    /////////////////////////////////////////////

    // Delete base geometries
    for(unsigned int i = 0; i < geoms->size(); i++) {
        delete(*geoms)[i];
    }
    delete geoms;
}


int
main()
{
    std::cout << "GEOS " << geosversion() << " ported from JTS " << jtsport() << std::endl;
    try {
        do_all();
    }
    // All exception thrown by GEOS are subclasses of this
    // one, so this is a catch-all
    catch(const GEOSException& exc) {
        std::cerr << "GEOS Exception: " << exc.what() << "\n";
        exit(1);
    }
    catch(const exception& e) {
        std::cerr << "Standard exception thrown: " << e.what() << std::endl;
        exit(1);
    }
    // and this is a catch-all non standard ;)
    catch(...) {
        std::cerr << "unknown exception trown!\n";
        exit(1);
    }

    // Unload is no more necessary
    //io::Unload::Release();

    exit(0);
}

*/


Point* create_point(double x, double y, GeometryFactory::unique_ptr gf) {
    Coordinate crd(x, y);
    Point* p = gf->createPoint(crd);
    return p;
}


Polygon* create_circle(double centerX, double centerY, double radius, GeometryFactory::unique_ptr gf) {
    geos::util::GeometricShapeFactory shapefactory(gf.get());
    shapefactory.setCentre(Coordinate(centerX, centerY));
    shapefactory.setSize(radius);
    // same as:
    //	shapefactory.setHeight(radius);
    //	shapefactory.setWidth(radius);
    return shapefactory.createCircle();
}


Polygon* create_ellipse(double centerX, double centerY, double width, double height, GeometryFactory::unique_ptr gf) {
    geos::util::GeometricShapeFactory shapefactory(gf.get());
    shapefactory.setCentre(Coordinate(centerX, centerY));
    shapefactory.setHeight(height);
    shapefactory.setWidth(width);
    return shapefactory.createCircle();
}

Polygon* create_rectangle(double llX, double llY, double width, double height, GeometryFactory::unique_ptr gf) {
    geos::util::GeometricShapeFactory shapefactory(gf.get());
    shapefactory.setBase(Coordinate(llX, llY));
    shapefactory.setHeight(height);
    shapefactory.setWidth(width);
    shapefactory.setNumPoints(4); // we don't need more then 4 points for a rectangle...
    // can use setSize for a square
    return shapefactory.createRectangle();
}



std::vector<Geometry*>* create_points(std::vector<double> &x, std::vector<double> &y) {	
	PrecisionModel* pm = new PrecisionModel();
	GeometryFactory::unique_ptr gf = GeometryFactory::create(pm, -1);
	delete pm;
	std::vector<Geometry*>* g = new std::vector<Geometry*>;
    for(unsigned int i = 0; i < x.size(); i++) {
		Coordinate crd(x[i], y[i]);
		Point* p = gf->createPoint(crd);
		g->push_back(p);
	}
	return g;
}


LineString* create_linestring(std::vector<double> x, std::vector<double> y, GeometryFactory::unique_ptr gf) {
    CoordinateArraySequence* cl = new CoordinateArraySequence();
	for (size_t i=0; i<x.size(); i++) {
		cl->add(Coordinate(x[i], y[i]));
	}	
    LineString* ls = gf->createLineString(cl);
    return ls;
}


MultiLineString* create_multilinestring(std::vector<Geometry*>* lines, GeometryFactory::unique_ptr gf) {
 	MultiLineString* ml = gf->createMultiLineString(lines);
	return ml;
}


LinearRing* create_linearring(std::vector<double> x, std::vector<double> y, GeometryFactory::unique_ptr gf) {
    CoordinateArraySequence* cl = new CoordinateArraySequence();
	for (size_t i=0; i<x.size(); i++) {
		cl->add(Coordinate(x[i], y[i]));
	}	
    LinearRing* lr = gf->createLinearRing(cl);
    return lr; 
}


MultiPolygon* create_multipolygon(std::vector<Geometry*>* pols, GeometryFactory::unique_ptr gf) {
 	MultiPolygon* mp = gf->createMultiPolygon(pols);
	return mp;
}


Polygon* create_polygon(std::vector<double> x, std::vector<double> y, 
				std::vector<std::vector<double>> hx, std::vector<std::vector<double>> hy,
				GeometryFactory::unique_ptr gf) {

    LinearRing* outer = create_linearring(x, y, gf);
//	vector<LinearRing*>* holes = new vector<LinearRing*>;
	std::vector<Geometry*>* holes = new std::vector<Geometry*>;
	for (size_t i=0; i<hx.size(); i++) {
		LinearRing* inner = create_linearring(hx[i], hy[i], gf);
		holes->push_back(inner);
	}
    Polygon* poly = gf->createPolygon(outer, holes);
    return poly;
}


// create a GeometryCollection containing copies of all Geometries in given vector.
GeometryCollection* create_simple_collection(std::vector<Geometry*>* geoms, GeometryFactory::unique_ptr gf) {
    return gf->createGeometryCollection(*geoms);
    // to transfer ownership of vector and its elements:
    // return gf->createGeometryCollection(geoms);
}



std::vector<Geometry*>* spat2geos(SpatVector* v) {

	std::vector<Geometry*>* geoms = new std::vector<Geometry*>;	
	PrecisionModel* pm = new PrecisionModel();
	GeometryFactory::unique_ptr geomfact = GeometryFactory::create(pm, -1);
	delete pm;

	std::string vt = v->type();
	if (vt == "points") {
		std::vector<std::vector<double>> xy = v->coordinates();
		geoms = create_points(xy[0], xy[1]);
	} else if (vt == "lines") {
		std::vector<double> x, y;
		size_t n = v->size();
		for (size_t i=0; i<n; i++) {
			SpatGeom g = v->getGeom(i);
			std::vector<Geometry*>* lns = new std::vector<Geometry*>;	
			for (size_t j=0; j < g.size(); j++) {
				SpatPart part = g.getPart(j);			
				LineString* ls = create_linestring(part.x, part.y, geomfact);
				lns->push_back(ls);
			}
			MultiLineString* mls = create_multilinestring(lns, geomfact);
			geoms->push_back(mls);			
		}
	} else if (vt == "polygons") {
		std::vector<double> x, y;
		std::vector<std::vector<double>> hx, hy;
		size_t n = v->size();
		for (size_t i=0; i<n; i++) {
			SpatGeom g = v->getGeom(i);
			std::vector<Geometry*>* pols = new std::vector<Geometry*>;	
			for (size_t j=0; j < g.size(); j++) {
				SpatPart part = g.getPart(j);			
				if (part.hasHoles()) {
					for (size_t k=0; k<part.nHoles(); k++) {
						SpatHole h = part.getHole(k);
						hx.push_back(h.x);
						hy.push_back(h.y);
					}
				}
				Polygon* poly = create_polygon(part.x, part.y, hx, hy, geomfact);
				pols->push_back(poly);
			}
			MultiPolygon* mp = create_multipolygon(pols, geomfact);
			geoms->push_back(mp);			
		}
	}	
	return geoms;
}


SpatVector geos2spat(std::vector<Geometry*>* geoms) {
	SpatVector out;
	size_t ng = geoms->size();
	std::vector<unsigned> gid, gp, hole;
	std::vector<double> x, y;
    for(size_t i = 0; i < ng; i++) {
        Geometry* g = (*geoms)[i];
		size_t np = g->getNumGeometries();
		for(size_t j = 0; j<np; j++) {
			const Geometry* part = g->getGeometryN(j);
			size_t npts = part->getNumPoints();
			CoordinateSequence* crds = part->getCoordinates(); 		
			for (size_t p=0; p < npts; p++) {
				x.push_back(crds->getX(p));
				y.push_back(crds->getY(p));
				gid.push_back(i);			
				gp.push_back(j);			
			}
		}
	}	
	hole.resize(x.size());
	out.setGeometry("polygons", gid, gp, x, y, hole);
	return out;
}



SpatVector SpatVector::buffer2(double d, unsigned segments, unsigned capstyle){
	SpatVector out;
	std::string vt = type();
	if ((vt == "points") && (d <= 0)) {
		out.setError("buffer size must be >= 0 with points");
		return out;
	}
	//std::vector<std::vector<double>> xy = coordinates();
	//std::vector<Geometry*>* points = create_Points(xy[0], xy[1]);
	std::vector<Geometry*>* geoms = spat2geos(this); 
	bool success = true;
	std::vector<Geometry*>* buf   = gBuffer(geoms, d, segments, capstyle, success);
	if (!success) {
		out.setError("buffer failed");
	} else {
		out = geos2spat(buf);
	}
	kill_geoms(geoms);
	kill_geoms(buf);
	return out;	
}

#endif

