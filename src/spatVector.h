#include "extent.h"
#include "dataframe.h"

enum SpatGeomType { points, lines, polygons, unknown };

class SpatHole {
	public:
		std::vector<double> x, y;  
		SpatExtent extent;
		SpatHole();
		SpatHole(std::vector<double> X, std::vector<double> Y);
};

class SpatPart {
	public:		
		std::vector<double> x, y; //, z;
		SpatExtent extent;
		SpatPart();
		SpatPart(std::vector<double> X, std::vector<double> Y);
		SpatPart(double X, double Y);

		// for POLYGONS only
		std::vector< SpatHole > holes; 
		bool addHole(std::vector<double> X, std::vector<double> Y);
		bool addHole(SpatHole h);
		SpatHole getHole(unsigned i) { return( holes[i] ) ; }
		bool hasHoles() { return holes.size() > 0;}
		unsigned nHoles() { return holes.size();}
};


class SpatGeom {
	public:
		SpatGeomType gtype;
		std::vector<SpatPart> parts; 
		SpatExtent extent;	
		SpatGeom();
		SpatGeom(SpatPart p);
		bool addPart(SpatPart p);
		bool addHole(SpatHole h);
		bool setPart(SpatPart p, unsigned i);
		SpatPart getPart(unsigned i);
		double area();
		double length();
		unsigned size() { return parts.size(); };
};


class SpatLayer {
	private:
		std::vector<SpatGeom> geoms; 
		SpatExtent extent;		
	public:
		std::string crs;
		std::vector<std::string> names();
		unsigned nrow();
		unsigned ncol();
		unsigned nxy();
		
		SpatExtent getExtent();
		
		std::string type();
		SpatGeomType getGType(std::string &type);

		std::string getCRS();
		void setCRS(std::string CRS);

		SpatDataFrame df;

		SpatGeom getGeom(unsigned i);
		bool addGeom(SpatGeom p);
		SpatDataFrame getGeometryDF();
		
		SpatLayer subset(std::vector<unsigned> range);
		void setGeometry(std::string type, std::vector<unsigned> id, std::vector<unsigned> part, std::vector<double> x, std::vector<double> y, std::vector<bool> hole);
		
		std::vector<double> area();
		std::vector<double> length();
		unsigned size();
	
		bool read(std::string fname);
		bool write(std::string filename, bool overwrite);

// attributes		
		std::vector<double> getDv(unsigned i);
		std::vector<long> getIv(unsigned i);
		std::vector<std::string> getSv(unsigned i);
		std::vector<unsigned> getItype();
		std::vector<unsigned> getIplace();


		SpatMessages msg;
		void setError(std::string s) { msg.setError(s); }
		void addWarning(std::string s) { msg.addWarning(s); }
	
};



class SpatVector {
	public:
		std::vector<SpatLayer> lyrs;
		unsigned layertypes();
		SpatMessages msg;
		void setError(std::string s) { msg.setError(s); }
		void addWarning(std::string s) { msg.addWarning(s); }
		
};

