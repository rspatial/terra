#include <Rcpp.h>
#include "spatraster.h"

#include "RcppFunctions.h"

RCPP_EXPOSED_CLASS(SpatExtent)
RCPP_EXPOSED_CLASS(SpatMessages)
RCPP_EXPOSED_CLASS(RasterSource)
RCPP_EXPOSED_CLASS(SpatRaster)
RCPP_EXPOSED_CLASS(SpatLayer)

RCPP_MODULE(spat){

    using namespace Rcpp;

    class_<SpatExtent>("SpatExtent")
		.constructor()
		.constructor<double, double, double, double>()
		.property("vector", &SpatExtent::asVector)		
		.property("valid", &SpatExtent::valid)		
	;	

    class_<SpatMessages>("SpatMessages")
		.constructor()
		.field("success", &SpatMessages::success)		
		.field("has_error", &SpatMessages::has_error)		
		.field("has_warning", &SpatMessages::has_warning)		
		.field("error", &SpatMessages::error)		
		.field("warnings", &SpatMessages::warnings)		
	;	


    class_<SpatLayer>("SpatLayer")
		.constructor()	
		.method("names", &SpatLayer::names, "names")		
		.method("nrow", &SpatLayer::nrow, "nrow")		
		.method("ncol", &SpatLayer::ncol, "ncol")		
		.method("size", &SpatLayer::size, "size")		
		.property("crs", &SpatLayer::getCRS, &SpatLayer::setCRS, "crs")		
		.method("type", &SpatLayer::type, "type")		
		.method("extent", &SpatLayer::getExtent, "extent")		
		.method("read", &SpatLayer::read, "read")		
		.method("getAttributes", &getAttributes, "getAttributes")
		.method("getGeometry", &getGeometry, "getGeometry")
		.method("setGeometry", &SpatLayer::setGeometry, "setGeometry")
		.field("messages", &SpatLayer::msg, "messages")
	;

	

    class_<RasterSource>("RasterSource")	
//		.field_readonly("memory", &RasterSource::memory)
		.field_readonly("filename", &RasterSource::filename)
//		.field_readonly("driver", &RasterSource::driver)
//		.field_readonly("nrow", &RasterSource::nrow)
//		.field_readonly("ncol", &RasterSource::ncol)
		.field_readonly("nlyr", &RasterSource::nlyr)
//		.field_readonly("crs", &RasterSource::crs)
//		.field_readonly("extent", &RasterSource::extent)
//		.field_readonly("datatype", &RasterSource::datatype)
		.field_readonly("NAflag", &RasterSource::NAflag)
		//std::vector<std::vector<int> > layers;		
	;	

    class_<SpatRaster>("SpatRaster")
		.constructor()
	    .constructor<std::string>()
		.constructor<std::vector<unsigned>, std::vector<double>, std::string>()
		
		.property("extent", &SpatRaster::getExtent, &SpatRaster::setExtent )
		.property("crs", &SpatRaster::getCRS, &SpatRaster::setCRS )
		.property("names", &SpatRaster::getNames)
		.property("res", &SpatRaster::resolution)
		.property("origin", &SpatRaster::origin)
		.property("inMemory", &SpatRaster::inMemory )
		.property("filenames", &SpatRaster::filenames )
		.property("hasValues", &SpatRaster::hasValues )

			
		//.field_readonly("hasRAT", &SpatRaster::hasRAT )
		//.field_readonly("hasCT", &SpatRaster::hasCT )
		.property("hasRange", &SpatRaster::hasRange )
		.property("range_min", &SpatRaster::range_min )
		.property("range_max", &SpatRaster::range_max )

		.field_readonly("nrow", &SpatRaster::nrow, "nrow")
		.field_readonly("ncol", &SpatRaster::ncol, "ncol")		
		.method("nsrc", &SpatRaster::nsrc, "nsrc" )	
		.field("messages", &SpatRaster::msg, "messages")
		.method("nlyrBySource", &SpatRaster::nlyrBySource, "nlyrBySource" )		
		.method("nlyr", &SpatRaster::nlyr, "nlyr" )
		.method("setNames", &SpatRaster::setNames, "setNames" )
		.field_readonly("source", &SpatRaster::source )

		
		//.method("addSources", ( std::SpatRaster (SpatRaster::*)(std::SpatRaster) )( &SpatRaster::addSources), "addSources")
		.method("addSources", &SpatRaster::addSources, "addSources")
		.method("subset", &SpatRaster::subset, "subset")

		
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
		.method("setValues", &SpatRaster::setValues)

		.method("setRange", &SpatRaster::setRange, "setRange")
		.method("writeStart", &SpatRaster::writeStart, "writeStart") 
		.method("writeStop", &SpatRaster::writeStop, "writeStop") 
		.method("writeValues", &SpatRaster::writeValues, "writeValues") 

		.method("writeRaster", &SpatRaster::writeRaster, "writeRaster")
		.method("canProcessInMemory", &SpatRaster::canProcessInMemory, "canProcessInMemory")
		.method("chunkSize", &SpatRaster::chunkSize, "chunkSize")
		
		.method("extract", &SpatRaster::extract, "extract")
		
		.method("rasterizePolygons", &SpatRaster::rasterizePolygons, "rasterizePolygons")
		.method("crop", &SpatRaster::crop, "crop")
		.method("focal", &SpatRaster::focal, "focal")
		.method("focalValues", &SpatRaster::focal_values, "focalValues")
		.method("trim", &SpatRaster::trim, "trim")
		.method("mask", &SpatRaster::mask, "mask")
		.method("aggregate", &SpatRaster::aggregate, "aggregate")
		.method("get_aggregates", &SpatRaster::get_aggregates, "get_aggregates")
		.method("get_aggregate_dims", &SpatRaster::get_aggregate_dims, "get_aggregate_dims")
		
		.method("arith_rast", ( SpatRaster (SpatRaster::*)(SpatRaster, std::string, std::string, bool) )( &SpatRaster::arith ))
		.method("arith_numb", ( SpatRaster (SpatRaster::*)(double, std::string, std::string, bool) )( &SpatRaster::arith ))
		.method("arith_rev", &SpatRaster::arith_rev, "arith_rev")
		.method("math", &SpatRaster::math, "math")
		.method("trig", &SpatRaster::trig, "trig")
		.method("cum", &SpatRaster::cum, "cum")
		.method("summary", &SpatRaster::summary, "summary")
		.method("summary_numb", &SpatRaster::summary_numb, "summary_numb")
		.method("logic_rast", ( SpatRaster (SpatRaster::*)(SpatRaster, std::string, std::string, bool) )( &SpatRaster::logic ))
		.method("logic_numb", ( SpatRaster (SpatRaster::*)(bool, std::string, std::string, bool) )( &SpatRaster::logic ))
	;
}


