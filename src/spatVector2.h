using namespace std;
#include "extent.h"
#include "dataframe.h"

enum GeomType { POINTS, LINES, POLYGONS };

class SpatPart {
	public:
		std::vector<double> x, y; 
		SpatExtent extent;
		bool set(std::vector<double> X, std::vector<double> Y) { 
			x = X; y = Y;  
			extent.xmin = *std::min_element(X.begin(), X.end());
			extent.xmax = *std::max_element(X.begin(), X.end());
			extent.ymin = *std::min_element(Y.begin(), Y.end());
			extent.ymax = *std::max_element(Y.begin(), Y.end());
			return true;
		}

		// for POLYGONS only
		std::vector< std::vector<double>> xHole, yHole; 
		bool hasHoles() { return xHole.size() > 0;}
		unsigned nHoles() { return xHole.size();}
		bool setHole(std::vector<double> X, std::vector<double> Y) { 
			// check if inside pol?
			xHole.push_back(X);
			yHole.push_back(Y);
			return true;
		}
		std::vector<double> getHoleX(unsigned i) { return( xHole[i] ) ; }
		std::vector<double> getHoleY(unsigned i) { return( yHole[i] ) ; }	
};


class SpatGeom {
	public:
		std::vector<SpatPart> parts; 
		SpatExtent extent;
		unsigned size() { return parts.size(); };
		
		SpatPart getPart(unsigned i) { return parts[i]; }
		bool addPart(SpatPart p) { 
			parts.push_back(p); 
			if (parts.size() > 1) {
				extent.unite(p.extent);
			} else {
				extent = p.extent;
			}
			return true; 
		}
		bool setPart(SpatPart p, unsigned i) {
			parts[i] = p;
			// now re-compute extent
			return true;
		}

		double area();
		double length();
		
};


class SpatLayer {
	private:
		std::vector<SpatGeom> geoms; 

	public:
		
		SpatExtent extent;		
		std::string crs;
		unsigned size() { return geoms.size(); };
		GeomType gtype;

		SpatDataFrame df;

		SpatGeom getGeom(unsigned i) { return geoms[i]; };
		bool addGeom(SpatGeom p) { 
			geoms.push_back(p); 
			if (Geoms.size() > 1) {
				extent.unite(p.extent);
			} else {
				extent = p.extent;
			}
			return true; 
		}
		
		SpatGeom subset(std::vector<unsigned> range) { 
			SpatGeom out;
			for (size_t i=0; i < range.size(); i++) {
				out.addGeom( geoms[range[i]] ); 
			}
			out.crs = crs;
			return out;	
		};

		std::vector<double> area();
		std::vector<double> lenght();
};


class SpatVector2 {

	public:
		std::vector<SpatLayer> lyrs;
		unsinged layertypes();
	
		bool read(std::string fname);
		bool write(std::string filename, bool overwrite);

		std::vector<string> names();
		unsigned nrow();
		unsigned ncol();
		string getCRS();
		void setCRS(std::string crs);
		SpatExtent extent();

// attributes		
		std::vector<double> getDv(unsigned i);
		std::vector<long> getIv(unsigned i);
		std::vector<string> getSv(unsigned i);
		std::vector<unsigned> getItype();
		std::vector<unsigned> getIplace();

		bool error = false;
		bool warning = false;
		string error_message;
		std::vector<string> warning_message;
	
};

