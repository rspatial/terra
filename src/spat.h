using namespace std;

#include <string>
#include <vector>
#include <algorithm>
#include <fstream>


class SpatExtent {
	public:
		double xmin, xmax, ymin, ymax;
		SpatExtent() {xmin = -180; xmax = 180; ymin = -90; ymax = 90;}
		SpatExtent(double _xmin, double _xmax, double _ymin, double _ymax) {xmin = _xmin; xmax = _xmax; ymin = _ymin; ymax = _ymax;}
		
		void intersect(SpatExtent e) { 
			xmin = std::max(xmin, e.xmin);
			xmax = std::min(xmax, e.xmax);
			ymin = std::max(ymin, e.ymin);
			ymax = std::min(ymax, e.ymax);
		}

		std::vector<double> asVector() { 
			std::vector<double> e(4);
			e[0] = xmin; e[1] = xmax; e[2] = ymin; e[3] = ymax; 
			return(e);
		}
			
		bool valid() {
			return ((xmax > xmin) && (ymax > ymin));
		}
};




class SpatPolyPart {
	public:
		std::vector<double> x, y; 
		std::vector< std::vector<double>> xHole, yHole; 

		SpatExtent extent;
		bool hasHoles() { return xHole.size() > 0;}
		unsigned nHoles() { return xHole.size();}
		bool set(std::vector<double> X, std::vector<double> Y) { 
			x = X; y = Y;  
			extent.xmin = *std::min_element(X.begin(), X.end());
			extent.xmax = *std::max_element(X.begin(), X.end());
			extent.ymin = *std::min_element(Y.begin(), Y.end());
			extent.ymax = *std::max_element(Y.begin(), Y.end());
			return true;
		}
		bool setHole(std::vector<double> X, std::vector<double> Y) { 
			xHole.push_back(X);
			yHole.push_back(Y);
			return true;
		}
		std::vector<double> getHoleX(unsigned i) { return( xHole[i] ) ; }
		std::vector<double> getHoleY(unsigned i) { return( yHole[i] ) ; }
		
};

class SpatPoly {
	public:
		std::vector<SpatPolyPart> parts; 
		SpatExtent extent;

		unsigned size() { return parts.size(); };
		SpatPolyPart getPart(unsigned i) { return parts[i]; }
		bool addPart(SpatPolyPart p) { 
			parts.push_back(p); 
			if (parts.size() > 1) {
				extent.xmin = std::min(extent.xmin, p.extent.xmin);
				extent.xmax = std::max(extent.xmax, p.extent.xmax);
				extent.ymin = std::min(extent.ymin, p.extent.ymin);
				extent.ymax = std::max(extent.ymax, p.extent.ymax);
			} else {
				extent = p.extent;
			}
			return true; 
		}
};


class SpatPolygons {
	public:
		std::vector<SpatPoly> polys; 
		SpatExtent extent;		
		std::string crs;
		std::vector<double> attr;

		unsigned size() { return polys.size(); };
		SpatPoly getPoly(unsigned i) { return polys[i]; };
		bool addPoly(SpatPoly p) { 
			polys.push_back(p); 
			if (polys.size() > 1) {
				extent.xmin = std::min(extent.xmin, p.extent.xmin);
				extent.xmax = std::max(extent.xmax, p.extent.xmax);
				extent.ymin = std::min(extent.ymin, p.extent.ymin);
				extent.ymax = std::max(extent.ymax, p.extent.ymax);
			} else {
				extent = p.extent;
			}
			attr.push_back(NAN);
			return true; 
		}
		double getAtt(unsigned i) {	return attr[i]; };
		bool setAtt(unsigned i, double a) { attr[i] = a; return true; };
		
		SpatPolygons subset(std::vector<unsigned> range) { 
			SpatPolygons out;
			for (size_t i=0; i < range.size(); i++) {
				out.addPoly( polys[range[i]] ); 
				out.attr.push_back(attr[i]);
			}
			out.crs = crs;
			return out;	
		};
		
};




class RasterSource {
	public:
		std::vector<bool> memory;
		std::vector<string> filename;
		std::vector<string> driver;
		std::vector<unsigned> nlayers;		
		std::vector<std::vector<int> > layers;		
		std::vector<string> datatype;
		std::vector<double> NAflag;
};


class BlockSize {
	public:
		std::vector<unsigned> row;
		std::vector<unsigned> nrows;
		unsigned n;
};




class SpatRaster {
	
	private:
		std::string msg;
		//fstream* fs;
		
	protected:
		SpatExtent extent;
		std::string crs ="+proj=longlat +datum=WGS84";
		void setnlyr() { 
			nlyr = std::accumulate(source.nlayers.begin(), source.nlayers.end(), 0); 
		}
		BlockSize getBlockSize();
		
	public:
		//double NA = std::numeric_limits<double>::quiet_NaN();
		RasterSource source;
		SpatExtent window;
	
	    std::vector<unsigned> getnlayers() {
			return source.nlayers;
		}		
		unsigned nrow, ncol, nlyr;
		unsigned size() { return ncol * nrow * nlyr ; }
		bool hasValues;
		
		BlockSize bs;
		
		std::vector<double> values;
		
		std::vector<bool> hasRange;
		std::vector<bool> hasRAT;
		std::vector<bool> hasCT;
		std::vector<double> range_min;
		std::vector<double> range_max;
		std::vector<string> names;
		std::vector<bool> inMemory() { return source.memory; }

		// constructors
		SpatRaster(std::string fname);
		SpatRaster();
		SpatRaster(std::vector<unsigned> rcl, std::vector<double> ext, std::string _crs);
		SpatRaster(unsigned _nrow, unsigned _ncol, unsigned _nlyr, SpatExtent ext, std::string _crs);

		SpatRaster arith(SpatRaster x, std::string oper, std::string filename="", bool overwrite=false);
		SpatRaster arith(double x, std::string oper, std:: string filename="", bool overwrite=false);
		SpatRaster arith_rev(double x, std::string oper, std::string filename, bool overwrite);
		
		SpatRaster operator + (SpatRaster x) { return arith(x, "+", "", false); }

		
		double ncell() { return nrow * ncol; }

//	void setExtent(std::vector<double> e) {	extent.xmin = e[0]; extent.xmax = e[1]; extent.ymin = e[2]; extent.ymax = e[3]; }
		SpatExtent getExtent() { return extent; }
		void setExtent(SpatExtent e) { extent = e ; }

		void setExtent(SpatExtent ext, bool keepRes=false, std::string snap="");

		
		std::string getCRS()	{ return(crs); }
		void setCRS(std::string _crs) { crs = _crs; }
		std::vector<string> getNames()	{ 
			if (names.size() < 1) {
				return std::vector<string> {"layer"}; // rep for each layer
			}
			return(names); 
		}
		void setNames(std::vector<string> _names) { names = _names; }
	
		std::vector<double> resolution() { return std::vector<double> { (extent.xmax - extent.xmin) / ncol, (extent.ymax - extent.ymin) / nrow };}
		double xres() { return (extent.xmax - extent.xmin) / ncol ;}
		double yres() { return (extent.ymax - extent.ymin) / nrow ;}

		std::vector<double> origin();	
		
		//std::vector<string> filenames() { return source.filename; }

		
		bool compare(unsigned nrows, unsigned ncols, SpatExtent e );
	
		std::vector<double> getValues();
		void setValues(std::vector<double> _values);


		bool constructFromFile(std::string fname);
		bool constructFromFileGDAL(std::string fname);

		std::vector<double> cellFromXY (std::vector<double> x, std::vector<double> y);
		double cellFromXY(double x, double y);
		std::vector<double> cellFromRowCol(std::vector<unsigned> rownr, std::vector<unsigned> colnr);
		double cellFromRowCol(unsigned rownr, unsigned colnr);
		std::vector<double> yFromRow(std::vector<unsigned> rownr);
		double yFromRow(unsigned rownr);
		std::vector<double> xFromCol(std::vector<unsigned> colnr);
		double xFromCol(unsigned colnr);
		std::vector<double> colFromX(std::vector<double> x);
		double colFromX(double x);
		std::vector<double> rowFromY(std::vector<double> y);
		double rowFromY(double y);
		std::vector< std::vector<double> > xyFromCell( std::vector<double> cell );
		std::vector< std::vector<double> > xyFromCell( double cell );
		std::vector< std::vector<double> > rowColFromCell(std::vector<double> cell);
		
		
		double valuesCell(double);
		double valuesCell(int, int);
		std::vector<double> valuesCell(std::vector<double>);	
		std::vector<double> valuesRow(int);	

		void setRange();
		
		
		bool readStart();
		std::vector<double> readValues(unsigned row, unsigned nrows, unsigned col, unsigned ncols);
		bool readStop();
		bool writeStart(std::string filename, bool overwrite);
		//bool writeStartFs(std::string filename, bool overwrite, fstream& f);		
		bool writeValues(std::vector<double> vals, unsigned row);
		bool writeStop();
		bool writeHDR();
		
		std::vector<double> readValuesGDAL(unsigned row, unsigned nrows, unsigned col, unsigned ncols);
		bool writeValuesGDAL(std::string filename, std::vector<double> values, std::string format="GTiff");
		
		
		void openFS(string const &filename);

		
		SpatRaster writeRaster(std::string filename, bool overwrite);
		SpatExtent align(SpatExtent e, string snap="near");
		
		SpatRaster test(string filename);
		SpatRaster crop(SpatExtent e, string filename="", string snap="near", bool overwrite=false);
		SpatRaster trim(unsigned padding=0, std::string filename="", bool overwrite=false);
		SpatRaster mask(SpatRaster mask, string filename="", bool overwrite=false);
		SpatRaster focal(std::vector<double> w, double fillvalue, bool narm, unsigned fun, std::string filename, bool overwrite);
		SpatRaster rasterizePolygons(SpatPolygons p, double background, string filename, bool overwrite);
		
		std::vector<double> focal_values(std::vector<unsigned> w, double fillvalue, unsigned row, unsigned nrows);

		SpatRaster aggregate(std::vector<unsigned> fact, string fun, bool narm, string filename="", bool overwrite=false);
		//std::vector<double> aggregate(std::vector<unsigned> fact, bool narm, string fun, string filename="");

		std::vector<unsigned> get_aggregate_dims( std::vector<unsigned> fact );
		std::vector<std::vector<double> > get_aggregates(std::vector<unsigned> dim);
		
		std::vector<double> sampleRegular(unsigned size, bool cells, bool asRaster);
};


/*
SpatRaster SQRT() {
	SpatRaster r = *this;
	std::transform(r.values.begin(), r.values.end(), r.values.begin(), (double(*)(double)) sqrt);
	return r;
}
		
SpatRaster SQRTfree(SpatRaster* g) {
	SpatRaster r = *g;
	std::transform(r.values.begin(), r.values.end(), r.values.begin(), (double(*)(double)) sqrt);
	return r;
}
*/


