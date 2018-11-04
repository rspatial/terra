#include "extent.h"
#include "dataframe.h"


class SpatHole {
	public:
		std::vector<double> x, y; 
		SpatExtent extent;
		SpatHole();
		SpatHole(std::vector<double> X, std::vector<double> Y);
};

class SpatPart {
	public:
		std::vector<double> x, y; 
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
		std::vector<SpatPart> parts; 
		SpatExtent extent;	
		SpatGeom();
		SpatGeom(SpatPart p);
		bool addPart(SpatPart p);
		bool setPart(SpatPart p, unsigned i);
		SpatPart getPart(unsigned i) { return parts[i]; }
		double area();
		double length();
		unsigned size() { return parts.size(); };
};


class SpatLayer {
	private:
		std::vector<SpatGeom> geoms; 
		SpatExtent extent;		
	public:
		enum GeomType { POINTS, LINES, POLYGONS };
		std::string crs;
		std::vector<std::string> names();
		unsigned nrow();
		unsigned ncol();
		unsigned nxy();
		
		SpatExtent getExtent();
		GeomType gtype;
		std::string type();
		std::string getCRS();
		void setCRS(std::string CRS);

		SpatDataFrame df;

		SpatGeom getGeom(unsigned i);
		bool addGeom(SpatGeom p);
		SpatDataFrame getGeometryDF();
		
		SpatLayer subset(std::vector<unsigned> range);

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

		bool error = false;
		bool warning = false;
		std::string error_message;
		std::vector<std::string> warning_message;
	
};


/*
class SpatVector {

	public:
		std::vector<SpatLayer> lyrs;
		unsigned layertypes();

		bool error = false;
		bool warning = false;
		std::string error_message;
		std::vector<std::string> warning_message;

};

*/