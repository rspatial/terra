#include <Rcpp.h>
#include "spat.h"

using namespace Rcpp;


NumericMatrix getValuesM(SpatRaster* r) {
	NumericMatrix x(r->ncell(), r->nlyr);
	std::vector<double> v;
	v = r->getValues();
	std::copy(v.begin(), v.end(), x.begin());
	return(x);
}

List getBlockSizeR(SpatRaster* r) {              //+1 for R
	List L = List::create(Named("row") = r->bs.row, Named("nrows") = r->bs.nrows, Named("n") = r->bs.n);
	return(L);
}

RCPP_EXPOSED_CLASS(SpatExtent)

RCPP_EXPOSED_CLASS(RasterSource)
RCPP_EXPOSED_CLASS(SpatRaster)

RCPP_EXPOSED_CLASS(SpatPolyPart)
RCPP_EXPOSED_CLASS(SpatPoly)
RCPP_EXPOSED_CLASS(SpatPolygons)


	
RCPP_MODULE(spat){

    using namespace Rcpp;


    class_<SpatPolyPart>("SpatPolyPart")
		.constructor()
		.field_readonly("x", &SpatPolyPart::x )
		.field_readonly("y", &SpatPolyPart::y )
		.field_readonly("extent", &SpatPolyPart::extent )
		.method("set", &SpatPolyPart::set, "set")
		.method("setHole", &SpatPolyPart::setHole, "setHole")
		.method("getHoleX", &SpatPolyPart::getHoleX, "getHoleX")
		.method("getHoleY", &SpatPolyPart::getHoleY, "getHoleY")
		.method("nHoles", &SpatPolyPart::nHoles, "nHoles")
		.method("hasHoles", &SpatPolyPart::hasHoles, "hasHoles")
		
	;	
    class_<SpatPoly>("SpatPoly")
		.constructor()
		.field_readonly("extent", &SpatPoly::extent )
		.method("getPart", &SpatPoly::getPart, "getPart")
		.method("addPart", &SpatPoly::addPart, "addPart")
		.method("size", &SpatPoly::size, "size")
		
	;	
    class_<SpatPolygons>("SpatPolygons")
//		.field("polygons", &SpatPolygons::polys )
		.field_readonly("extent", &SpatPolygons::extent )
		.field("attr", &SpatPolygons::attr )
		.field("crs", &SpatPolygons::crs )
		.constructor()
		.method("getPoly", &SpatPolygons::getPoly, "getPoly")
		.method("addPoly", &SpatPolygons::addPoly, "addPoly")
		.method("size", &SpatPolygons::size, "size")

		.method("getAtt", &SpatPolygons::getAtt, "getAtt")
		.method("setAtt", &SpatPolygons::setAtt, "setAtt")
	;	

	
    class_<SpatExtent>("SpatExtent")
		.constructor()
		.constructor<double, double, double, double>()
		.property("vector", &SpatExtent::asVector)		
		.property("valid", &SpatExtent::valid)		
	;	

    class_<RasterSource>("RasterSource")
	
		.field_readonly("memory", &RasterSource::memory)
		.field_readonly("filename", &RasterSource::filename)
		.field_readonly("driver", &RasterSource::driver)
		.field_readonly("nlayers", &RasterSource::nlayers)
		.field_readonly("datatype", &RasterSource::datatype)
		.field_readonly("NAflag", &RasterSource::NAflag)
		//std::vector<std::vector<int> > layers;		
	;	

    class_<SpatRaster>("SpatRaster")
		.constructor()
	    .constructor<std::string>()
		.constructor<std::vector<unsigned>, std::vector<double>, std::string>()

		.field_readonly("source", &SpatRaster::source )
		
		.method("cellFromXY", ( std::vector<double> (SpatRaster::*)(std::vector<double>,std::vector<double>) )( &SpatRaster::cellFromXY ))
		.method("cellFromRowCol", ( std::vector<double> (SpatRaster::*)(std::vector<unsigned>,std::vector<unsigned>) )( &SpatRaster::cellFromRowCol ))
		.method("yFromRow", ( std::vector<double> (SpatRaster::*)(std::vector<unsigned>) )( &SpatRaster::yFromRow ))
		.method("xFromCol", ( std::vector<double> (SpatRaster::*)(std::vector<unsigned>) )( &SpatRaster::xFromCol ))
		.method("colFromX", ( std::vector<double> (SpatRaster::*)(std::vector<double>) )( &SpatRaster::colFromX ))
		.method("rowFromY", ( std::vector<double> (SpatRaster::*)(std::vector<double>) )( &SpatRaster::rowFromY ))
		.method("xyFromCell", ( std::vector< std::vector<double> > (SpatRaster::*)(std::vector<double>) )( &SpatRaster::xyFromCell ))
		.method("rowColFromCell", ( std::vector< std::vector<double> > (SpatRaster::*)(std::vector<double>) )( &SpatRaster::rowColFromCell ))

		.method("readStart", &SpatRaster::readStart, "readStart") 
		.method("readStop", &SpatRaster::readStop, "readStop") 
		.method("readValues", &SpatRaster::readValues, "readValues")	
		.method("getValues", &SpatRaster::getValues)
		.method("getBlockSize", &getBlockSizeR)
//		.property("values", &SpatRaster::getValues, &SpatRaster::setValues )	
		.method("setValues", &SpatRaster::setValues)

		.method("setRange", &SpatRaster::setRange, "setRange")
		.method("writeStart", &SpatRaster::writeStart, "writeStart") 
		.method("writeStop", &SpatRaster::writeStop, "writeStop") 
		.method("writeValues", &SpatRaster::writeValues, "writeValues") 

		.method("writeRaster", &SpatRaster::writeRaster, "writeRaster")
		
		.property("crs", &SpatRaster::getCRS, &SpatRaster::setCRS )
		.property("extent", &SpatRaster::getExtent, &SpatRaster::setExtent )
		.property("names", &SpatRaster::getNames, &SpatRaster::setNames )
		.property("res", &SpatRaster::resolution)
		.property("origin", &SpatRaster::origin)
		//.property("layers", &SpatRaster::getnlayers)


		.property("inMemory", &SpatRaster::inMemory )
//		.property("filenames", &SpatRaster::filenames )
		
		.field_readonly("nrow", &SpatRaster::nrow )
		.field_readonly("ncol", &SpatRaster::ncol )
		.field_readonly("nlyr", &SpatRaster::nlyr )
	
		.field_readonly("hasValues", &SpatRaster::hasValues )
		.field_readonly("hasRange", &SpatRaster::hasRange )
		.field_readonly("range_min", &SpatRaster::range_min )
		.field_readonly("range_max", &SpatRaster::range_max )

		.method("test", &SpatRaster::test, "test")
		
		
		.method("rasterizePolygons", &SpatRaster::rasterizePolygons, "rasterizePolygons")
		.method("crop", &SpatRaster::crop, "crop")
		.method("focal", &SpatRaster::focal, "focal")
		.method("focalValues", &SpatRaster::focal_values, "focalValues")
		.method("trim", &SpatRaster::trim, "trim")
		.method("mask", &SpatRaster::mask, "mask")
		.method("aggregate", &SpatRaster::aggregate, "aggregate")
		.method("get_aggregates", &SpatRaster::get_aggregates, "get_aggregates")
		.method("get_aggregate_dims", &SpatRaster::get_aggregate_dims, "get_aggregate_dims")
			
	;	
		
}
