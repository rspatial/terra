using namespace std;

#include <fstream>
#include <numeric>
#include <cmath>
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
		string endian;
		string datatype;
		double NAflag;
};


class BlockSize {
	public:
		std::vector<unsigned> row;
		std::vector<unsigned> nrows;
		unsigned n;
};


class SpatRaster {
	
	private:
		//fstream* fs;
		
	protected:
		SpatExtent extent;
		SpatExtent window;
		std::string crs ="+proj=longlat +datum=WGS84";

		void setnlyr() { 
			nlyr=0;
			for (size_t i=0; i < source.size(); i++) {
				nlyr += source[i].layers.size();
			}
		}
		BlockSize getBlockSize();
		
	public:
		std::vector<string> msg;
		bool error;
		bool warning;
		
		//double NA = std::numeric_limits<double>::quiet_NaN();
		std::vector<RasterSource> source;
	
		//unsigned getnlayers() { return nlyr; }
		
		unsigned nrow, ncol, nlyr;
		unsigned long size() { return ncol * nrow * nlyr ; }
		bool hasValues;
		
		BlockSize bs;
		
		std::vector<double> values;
		
		std::vector<string> filenames() { 
			std::vector<string> f(source.size()); 
			for (size_t i=0; i<f.size(); i++) {
				f[i] = source[i].filename; 
			} 
			return(f);
		}

		
		std::vector<bool> hasRange() {
			return source[0].hasRange;
		}
		std::vector<double> range_min(){
			return source[0].range_min;	
		}
		std::vector<double> range_max(){
			return source[0].range_max;	
		}
//		std::vector<string> names() {
//			return source[0].names;
//		}
		std::vector<bool> inMemory() { 
			std::vector<bool> m(source.size()); 
			for (size_t i=0; i<m.size(); i++) {
				m[i] = source[i].memory; 
			} 
			return(m); 
		}
		
		std::vector<string> getNames()	{ 
			return source[0].names;
		}
		void setNames(std::vector<string> _names) { source[0].names = _names; }
	
		// constructors
		SpatRaster(std::string fname);
		SpatRaster();
		SpatRaster(RasterSource s);
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
	
		std::vector<double> resolution() { return std::vector<double> { (extent.xmax - extent.xmin) / ncol, (extent.ymax - extent.ymin) / nrow };}
		double xres() { return (extent.xmax - extent.xmin) / ncol ;}
		double yres() { return (extent.ymax - extent.ymin) / nrow ;}

		std::vector<double> origin();	
		
		
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


