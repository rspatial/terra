using namespace std;

#include <fstream>
#include <numeric>
#include "spatvector.h"


class RasterSource {
	public:
		unsigned ncol;
		unsigned nrow;
		SpatExtent extent;
		unsigned nlyr;
		string crs;
		std::vector<unsigned> layers;
		std::vector<string> names;

		bool hasValues;
		std::vector<bool> hasRange;
		std::vector<double> range_min;
		std::vector<double> range_max;
		std::vector<bool> hasCT;
		std::vector<bool> hasRAT;
		
		bool memory;
		string filename;
		unsigned nlyrfile;
		string driver;
		string bandorder;
		string byteorder;
		string datatype;
		double NAflag;
};


class BlockSize {
	public:
		std::vector<unsigned> row;
		std::vector<unsigned> nrows;
		unsigned n;
};

class AttributeTable {
	public:
		std::vector<unsigned> code;
		std::vector<string> value;
};

class ColorTable {
	public:
		std::vector<unsigned> code;
		std::vector<string> value;
};

class SpatRaster {
	
	private:
		//fstream* fs;
		
	protected:
		SpatExtent extent;
		SpatExtent window;
		std::string crs; 
		BlockSize getBlockSize();
		std::vector<double> values;
		
	public:
		
		//double NA = std::numeric_limits<double>::quiet_NaN();

////////////////////////////////////////////////////	
// properties and property-like methods for entire object
////////////////////////////////////////////////////		
		bool error = false;
		bool warning = false;
		string error_message = "";
		string warning_message = "";

		std::vector<RasterSource> source;
		BlockSize bs;
	
		unsigned nrow, ncol;
		unsigned long size() { return ncol * nrow * nlyr() ; }
		SpatExtent getExtent() { return extent; }
		void setExtent(SpatExtent e) { extent = e ; }
		void setExtent(SpatExtent ext, bool keepRes=false, std::string snap="");  // also set it for sources?
		std::string getCRS() { return(crs); }
		void setCRS(std::string _crs) { crs = _crs; } // also set it for sources?
		std::vector<double> resolution() { return std::vector<double> { (extent.xmax - extent.xmin) / ncol, (extent.ymax - extent.ymin) / nrow };}
		double ncell() { return nrow * ncol; }
		double xres() { return (extent.xmax - extent.xmin) / ncol ;}
		double yres() { return (extent.ymax - extent.ymin) / nrow ;}
		std::vector<double> origin();	
		unsigned nlyr() { 
			unsigned x = 0;
			for (size_t i=0; i<source.size(); i++) { x += source[i].nlyr; } 
			return(x);
		}
		
		// only no values allowed with a single RasterSource
		bool hasValues() { return source[0].hasValues ; };
		std::vector<double> getValues();
		void setValues(std::vector<double> _values);

////////////////////////////////////////////////////
// property like methods for RasterSources
////////////////////////////////////////////////////
		std::vector<string> filenames() { 
			std::vector<string> x(source.size()); 
			for (size_t i=0; i<x.size(); i++) { x[i] = source[i].filename; } 
			return(x);
		}
		std::vector<bool> inMemory() { 
			std::vector<bool> m(source.size()); 
			for (size_t i=0; i<m.size(); i++) { m[i] = source[i].memory; } 
			return(m); 
		}
		

////////////////////////////////////////////////////	
// property like methods for Layers
////////////////////////////////////////////////////		

		std::vector<bool> hasRange() {
			std::vector<bool> x; 
			for (size_t i=0; i<source.size(); i++) { x.insert(x.end(), source[i].hasRange.begin(), source[i].hasRange.end()); }
			return(x);
		}

		std::vector<double> range_min() { 
			std::vector<double> x; 
			for (size_t i=0; i<source.size(); i++) { x.insert(x.end(), source[i].range_min.begin(),source[i].range_min.end()); }
			return(x);
		}

		std::vector<double> range_max() { 
			std::vector<double> x; 
			for (size_t i=0; i<source.size(); i++) { x.insert(x.end(), source[i].range_max.begin(), source[i].range_max.end()); }
			return(x);
		}

		std::vector<string> getNames() { 
			std::vector<string> x; 
			for (size_t i=0; i<source.size(); i++) { x.insert(x.end(), source[i].names.begin(), source[i].names.end()); }
			return(x);
		}

		void setNames(std::vector<string> _names) {
			size_t begin=0;
			size_t end;
			for (size_t i=0; i<source.size(); i++)	{
				end = begin + source[i].nlyr-1;	
				source[i].names = std::vector<string> (_names.begin() + begin, _names.end() + end) ;
				begin = end + 1;
			}
		}
		


////////////////////////////////////////////////////		
// constructors 
////////////////////////////////////////////////////		

		SpatRaster();
		SpatRaster(unsigned _nrow, unsigned _ncol, unsigned _nlyr, SpatExtent ext, std::string _crs);
		SpatRaster(std::vector<unsigned> rcl, std::vector<double> ext, std::string _crs);

		SpatRaster(std::string fname);
		SpatRaster(RasterSource s);
		void setSource(RasterSource s);
		SpatRaster(const SpatRaster& x);
				
		bool constructFromFile(std::string fname);
		bool constructFromFileGDAL(std::string fname);
		
		SpatRaster addSource(SpatRaster x);
		//SpatRaster shallowCopy();

////////////////////////////////////////////////////		
// helper methods
////////////////////////////////////////////////////		

		bool compare_geom(SpatRaster x, bool lyrs, bool crs);
		
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
		
////////////////////////////////////////////////////		
// read and write
////////////////////////////////////////////////////		
		
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
	
////////////////////////////////////////////////////		
// main methods
////////////////////////////////////////////////////		

		SpatRaster arith(SpatRaster x, std::string oper, std::string filename="", bool overwrite=false);
		SpatRaster arith(double x, std::string oper, std:: string filename="", bool overwrite=false);
		SpatRaster arith_rev(double x, std::string oper, std::string filename, bool overwrite);
		
		SpatRaster operator + (SpatRaster x) { return arith(x, "+", "", false); }
		SpatRaster test(string filename);

		SpatRaster aggregate(std::vector<unsigned> fact, string fun, bool narm, string filename="", bool overwrite=false);
		std::vector<unsigned> get_aggregate_dims( std::vector<unsigned> fact );
		std::vector<std::vector<double> > get_aggregates(std::vector<unsigned> dim);

		SpatExtent align(SpatExtent e, string snap="near");
		SpatRaster crop(SpatExtent e, string filename="", string snap="near", bool overwrite=false);
		SpatRaster trim(unsigned padding=0, std::string filename="", bool overwrite=false);
		SpatRaster mask(SpatRaster mask, string filename="", bool overwrite=false);
		SpatRaster focal(std::vector<double> w, double fillvalue, bool narm, unsigned fun, std::string filename, bool overwrite);
		std::vector<double> focal_values(std::vector<unsigned> w, double fillvalue, unsigned row, unsigned nrows);
		SpatRaster rasterizePolygons(SpatPolygons p, double background, string filename, bool overwrite);	
		
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


