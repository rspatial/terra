#include <Rcpp.h>
#include "spat.h"

using namespace Rcpp;


NumericMatrix getValuesM(SpatRaster* r) {
	NumericMatrix x(r->ncell(), r->nlyr() );
	std::vector<double> v;
	v = r->getValues();
	std::copy(v.begin(), v.end(), x.begin());
	return(x);
}


List getBlockSizeR(SpatRaster* r, unsigned n) {              //+1 for R
    BlockSize bs = r->getBlockSize(n);
	List L = List::create(Named("row") = bs.row, Named("nrows") = bs.nrows, Named("n") = bs.n);
	return(L);
}

RCPP_EXPOSED_CLASS(SpatExtent)

RCPP_EXPOSED_CLASS(RasterSource)
RCPP_EXPOSED_CLASS(SpatRaster)

RCPP_EXPOSED_CLASS(SpatGeomRing)
RCPP_EXPOSED_CLASS(SpatGeomRings)
RCPP_EXPOSED_CLASS(SpatPolygons)


	
RCPP_MODULE(spat){

    using namespace Rcpp;


    class_<SpatGeomRing>("SpatGeomRing")
		.constructor()
		.field_readonly("x", &SpatGeomRing::x )
		.field_readonly("y", &SpatGeomRing::y )
		.field_readonly("extent", &SpatGeomRing::extent )
		.method("set", &SpatGeomRing::set, "set")
		.method("setHole", &SpatGeomRing::setHole, "setHole")
		.method("getHoleX", &SpatGeomRing::getHoleX, "getHoleX")
		.method("getHoleY", &SpatGeomRing::getHoleY, "getHoleY")
		.method("nHoles", &SpatGeomRing::nHoles, "nHoles")
		.method("hasHoles", &SpatGeomRing::hasHoles, "hasHoles")
		
	;	
    class_<SpatGeomRings>("SpatGeomRings")
		.constructor()
		.field_readonly("extent", &SpatGeomRings::extent )
		.method("getPart", &SpatGeomRings::getGeom, "getPart")
		.method("addPart", &SpatGeomRings::addGeom, "addPart")
		.method("size", &SpatGeomRings::size, "size")
		
	;	
	
    class_<SpatPolygons>("SpatPolygons")
//		.field("polygons", &SpatPolygons::polys )
		.field_readonly("extent", &SpatPolygons::extent )
		.field("attr", &SpatPolygons::attr )
		.field("crs", &SpatPolygons::crs )
		.constructor()
		.method("getPoly", &SpatPolygons::getGeometry, "getPoly")
		.method("addPoly", &SpatPolygons::addGeometry, "addPoly")
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
		.field_readonly("nrow", &RasterSource::nrow)
		.field_readonly("ncol", &RasterSource::ncol)
		.field_readonly("nlyr", &RasterSource::nlyr)
		.field_readonly("crs", &RasterSource::crs)
		.field_readonly("extent", &RasterSource::extent)
		.field_readonly("datatype", &RasterSource::datatype)
		.field_readonly("NAflag", &RasterSource::NAflag)
		//std::vector<std::vector<int> > layers;		
	;	

    class_<SpatRaster>("SpatRaster")
		.constructor()
	    .constructor<std::string>()
		.constructor<std::vector<unsigned>, std::vector<double>, std::string>()
		
		//.method("addSources", ( std::SpatRaster (SpatRaster::*)(std::SpatRaster) )( &SpatRaster::addSources), "addSources")
		.method("addSources", &SpatRaster::addSources, "addSources")
		.method("subset", &SpatRaster::subset, "subset")

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
		.method("canProcessInMemory", &SpatRaster::canProcessInMemory, "canProcessInMemory")
		.method("chunkSize", &SpatRaster::chunkSize, "chunkSize")
		
		.field_readonly("nrow", &SpatRaster::nrow, "nrow")
		.field_readonly("ncol", &SpatRaster::ncol, "ncol")		
		.method("nsrc", &SpatRaster::nsrc, "nsrc" )	
		
		.method("nlyrBySource", &SpatRaster::nlyrBySource, "nlyrBySource" )		
				
		.method("nlyr", &SpatRaster::nlyr, "nlyr" )
		.property("extent", &SpatRaster::getExtent, &SpatRaster::setExtent )
		.property("crs", &SpatRaster::getCRS, &SpatRaster::setCRS )

		.property("names", &SpatRaster::getNames)
		.method("setNames", &SpatRaster::setNames, "setNames" )

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

		
		.field("error", &SpatRaster::error )
		.field("warning", &SpatRaster::warning )
		.field("error_message", &SpatRaster::error_message )
		.field("warning_message", &SpatRaster::warning_message )
		
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


