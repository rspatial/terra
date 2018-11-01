	//using namespace std;

#include <fstream>
#include <numeric>
#include "spatvector.h"

#define useGDAL

#ifdef useGDAL
#include "gdal_priv.h"
#endif

class RasterAttributeTable {
	public:
		std::vector<unsigned> code;
		std::vector<std::string> value;
};

class ColorTable {
	public:
		std::vector<unsigned> code;
		std::vector<std::string> value;
};

class RasterSource {
	public:
		unsigned ncol, nrow, nlyr;
		SpatExtent extent;
		std::string crs;
		std::vector<unsigned> layers;
		std::vector<std::string> names;

		//std::vector< std::vector<double> values;
        std::vector<double> values;
        std::vector<int64_t> ivalues;
        std::vector<bool> bvalues;

		std::vector<bool> hasRange;
		std::vector<double> range_min;
		std::vector<double> range_max;
		std::vector<double> time;
		std::vector<bool> hasCT;
		std::vector<bool> hasRAT;
		std::vector<RasterAttributeTable> RAT;
		std::vector<ColorTable> CT;

		bool memory;
		bool hasValues;
		std::string filename;
		//unsigned nlyrfile;

		// for native files
		std::string driver;
		std::string bandorder;
		std::string byteorder;
		std::string datatype;
		double NAflag;

		std::vector<RasterSource> subset(std::vector<unsigned> lyrs);
		std::vector<double> getValues(unsigned lyr);
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
#ifdef useGDAL
		GDALDataset* gdalconnection;
#endif
		bool open_read;
		bool open_write;
	protected:
		SpatExtent extent;
		SpatExtent window;
		std::string crs;

	public:

////////////////////////////////////////////////////
// properties and property-like methods for entire object
////////////////////////////////////////////////////
		unsigned nrow, ncol;
		std::vector<RasterSource> source;
		BlockSize bs;
		BlockSize getBlockSize(unsigned n);

		bool error = false;
		bool warning = false;
		std::string error_message;
		std::vector<std::string> warning_message;


		//double NA = std::numeric_limits<double>::quiet_NaN();

		unsigned long size() { return ncol * nrow * nlyr() ; }
		SpatExtent getExtent() { return extent; }
		void setExtent(SpatExtent e) { extent = e ; }
		void setExtent(SpatExtent ext, bool keepRes=false, std::string snap="");  // also set it for sources?
		std::string getCRS() { return(crs); }
		void setCRS(std::string _crs); 

		std::vector<double> resolution();
		double ncell() { return nrow * ncol; }
		double xres() { return (extent.xmax - extent.xmin) / ncol ;}
		double yres() { return (extent.ymax - extent.ymin) / nrow ;}
		std::vector<double> origin();
		unsigned nlyr();

		// only no values allowed with a single RasterSource
		bool hasValues() { return source[0].hasValues ; };
		std::vector<double> getValues();
		bool setValues(std::vector<double> _values);

////////////////////////////////////////////////////
// property like methods for RasterSources
////////////////////////////////////////////////////
		std::vector<std::string> filenames();
		bool isSource(std::string filename);
		std::vector<bool> inMemory();


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

		std::vector<std::string> getNames();
		bool setNames(std::vector<std::string> names);

////////////////////////////////////////////////////
// constructors
////////////////////////////////////////////////////

		SpatRaster();
		SpatRaster(unsigned _nrow, unsigned _ncol, unsigned _nlyr, SpatExtent ext, std::string _crs);
		SpatRaster(std::vector<unsigned> rcl, std::vector<double> ext, std::string _crs);

		SpatRaster(std::string fname);
		SpatRaster(RasterSource s);
		void setSource(RasterSource s);
		void setSources(std::vector<RasterSource> s);
		//SpatRaster(const SpatRaster& x);

        SpatRaster deepCopy();
        SpatRaster geometry(long nlyrs=-1);

		bool constructFromFile(std::string fname);
		bool constructFromFileGDAL(std::string fname);

		SpatRaster addSources(SpatRaster x);
		SpatRaster subset(std::vector<unsigned> lyrs, std::string filename, bool overwrite);

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
        std::vector<unsigned> sourcesFromLyrs(std::vector<unsigned> lyrs);

		int sourceFromLyr(unsigned lyr);
        std::vector<unsigned> nlyrBySource();
        std::vector<unsigned> lyrsBySource();
        unsigned nsrc();


		double valuesCell(double);
		double valuesCell(int, int);
		std::vector<double> valuesCell(std::vector<double>);
		std::vector<double> valuesRow(int);

		void setRange();

////////////////////////////////////////////////////
// read and write
////////////////////////////////////////////////////

		bool readStart();
		std::vector<double> readValues(unsigned row, unsigned nrows, unsigned col, unsigned ncols, unsigned lyr, unsigned nlyrs);
		std::vector<double> readBlock(BlockSize bs, unsigned i);

		bool readStop();
		std::vector<double> readValuesGDAL(unsigned row, unsigned nrows, unsigned col, unsigned ncols, unsigned lyr, unsigned nlyrs);

		bool readStartGDAL();
		bool readStopGDAL();
		std::vector<double> readChunkGDAL(unsigned row, unsigned nrows, unsigned col, unsigned ncols, unsigned lyr, unsigned nlyrs);


		bool writeStart(std::string filename, bool overwrite);
		bool writeValues(std::vector<double> vals, unsigned row);
		bool writeStop();
		bool writeHDR(std::string filename);

		bool writeStopGDAL();
		bool writeValuesGDAL(std::vector<double> vals, unsigned row);
		bool writeStartGDAL(std::string filename, bool overwrite);

		void openFS(std::string const &filename);

		bool writeRaster(std::string filename, bool overwrite);
		bool writeRasterGDAL(std::string filename, bool overwrite);

		bool canProcessInMemory(unsigned n);
		unsigned chunkSize(unsigned n);


////////////////////////////////////////////////////
// main methods
////////////////////////////////////////////////////

		SpatRaster arith(SpatRaster x, std::string oper, std::string filename, bool overwrite=false);
		SpatRaster arith(double x, std::string oper, std::string filename, bool overwrite=false);
		SpatRaster arith_rev(double x, std::string oper, std::string filename, bool overwrite);
		SpatRaster math(std::string fun, std::string filename, bool overwrite);
		SpatRaster trig(std::string fun, std::string filename, bool overwrite);
		SpatRaster logic(SpatRaster x, std::string oper, std::string filename, bool overwrite=false);
		SpatRaster logic(bool x, std::string oper, std::string filename, bool overwrite=false);
		SpatRaster isnot(std::string filename, bool overwrite);
		SpatRaster cum(std::string fun, bool narm, std::string filename, bool overwrite);
		SpatRaster summary(std::string fun, bool narm, std::string filename, bool overwrite);
		SpatRaster summary_numb(std::string fun, std::vector<double> add, bool narm, std::string filename, bool overwrite);

		SpatRaster operator + (SpatRaster x) { return arith(x, "+", "", false); }
		//SpatRaster test(std::string filename);

		SpatRaster aggregate(std::vector<unsigned> fact, std::string fun, bool narm, std::string filename="", bool overwrite=false);
		std::vector<unsigned> get_aggregate_dims( std::vector<unsigned> fact );
		std::vector<std::vector<double> > get_aggregates(std::vector<unsigned> dim);

		SpatExtent align(SpatExtent e, std::string snap="near");
		SpatRaster crop(SpatExtent e, std::string filename="", std::string snap="near", bool overwrite=false);
		SpatRaster trim(unsigned padding=0, std::string filename="", bool overwrite=false);
		SpatRaster mask(SpatRaster x, std::string filename="", bool overwrite=false);
		SpatRaster focal(std::vector<double> w, double fillvalue, bool narm, unsigned fun, std::string filename, bool overwrite);
		std::vector<double> focal_values(std::vector<unsigned> w, double fillvalue, unsigned row, unsigned nrows);
		SpatRaster rasterizePolygons(SpatPolygons p, double background, std::string filename, bool overwrite);

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


