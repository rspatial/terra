#include "spatRaster.h"

#if GDAL_VERSION_MAJOR >= 3 && GDAL_VERSION_MINOR >= 1

#include "proj.h"


//#include "cpl_conv.h" // for CPLMalloc()
//#include "cpl_string.h"
#include "ogr_spatialref.h"

#include "gdal_priv.h"
#include "gdal.h"
#include "crs.h"
#include "string_utils.h"


bool SpatRaster::constructFromFileMulti(std::string fname, std::string sub, std::vector<size_t> xyz) {

//    GDALAllRegister();
	
	size_t ndims = xyz.size();
	if (ndims < 2 || ndims > 3) {
		setError("need two or three dimension variables");
        return false;
	}
	
	
    auto poDataset = std::unique_ptr<GDALDataset>(
        GDALDataset::Open(fname.c_str(), GDAL_OF_MULTIDIM_RASTER ));
    if( !poDataset ) {
		setError("not a good dataset");
        return false;
    }

	auto poRootGroup = poDataset->GetRootGroup();
    if( !poRootGroup ) {
		setError("no roots");
		return false;
    }

	bool warngroup = false;
	std::vector<std::string> gnames;
	if (sub == "") {
		char** papszOptions = NULL;
		gnames = poRootGroup->GetMDArrayNames(papszOptions);
		CSLDestroy(papszOptions);
		sub = gnames[0];
		if (gnames.size() > 1) warngroup = true;
	}

    auto poVar = poRootGroup->OpenMDArray(sub.c_str());
    if( !poVar )   {
		setError("cannot find: " + sub);
		return false;
    }

	SpatRasterSource s;

//	const auto sref:
	std::string wkt = "";
	std::shared_ptr<OGRSpatialReference> srs = poVar->GetSpatialRef();
	if (srs != NULL) {
		char *cp;
		const char *options[3] = { "MULTILINE=YES", "FORMAT=WKT2", NULL };
		OGRErr err = srs->exportToWkt(&cp, options);
		if (err == OGRERR_NONE) {
			wkt = std::string(cp);
		} 
		CPLFree(cp);
	}
	std::string msg;
	if (!s.srs.set({wkt}, msg)) {
		addWarning(msg);
	}

	s.multidim = true;

    for( const auto poDim: poVar->GetDimensions() ) {
        s.m_dims.push_back(static_cast<size_t>(poDim->GetSize()));
        s.m_dimnames.push_back(static_cast<std::string>(poDim->GetName()));
    }

	if (warngroup) {
	//to do: remove dimnames from gnames
		std::string gn = "";
		for (size_t i=0; i<gnames.size(); i++) {
			if (!is_in_vector(gnames[i], s.m_dimnames) &&  (gnames[i] != sub)) { 
				gn += gnames[i] + ", ";
			}
		}
		addWarning("using: " + sub + ". Other groups are: \n" + gn);
	}
	s.source_name = sub;
	s.source_name_long = poVar->GetFullName();

	Rcpp::Rcout << poVar->GetName() << std::endl;
	Rcpp::Rcout << poVar->GetFullName() << std::endl;

	//std::vector<size_t> xyz = {0,1,2};
	s.ncol = s.m_dims[xyz[0]];
	s.nrow = s.m_dims[xyz[1]];
	if (ndims > 2) {
		s.nlyr = s.m_dims[xyz[2]];
	} else {
		s.nlyr = 1;		
	}
	s.nlyrfile = s.nlyr;
	s.layers.resize(s.nlyr);
    std::iota(s.layers.begin(), s.layers.end(), 0);
	s.flipped = false;
	s.rotated = false;
	s.memory = false;
	s.filename = fname;
	s.hasValues = true;
	s.unit = std::vector<std::string>(s.nlyr, poVar->GetUnit());


// layer names
// time 
// extent

	setSource(s);
	
	for (size_t i=0; i<s.m_dims.size(); i++){
		Rcpp::Rcout << s.m_dims[i] << " " << s.m_dimnames[i] << std::endl;
	}
	return true;
}


/*
read
    std::vector<double> values(nValues);
    poVar->Read(std::vector<GUInt64>{0,0,0}.data(),
                anCount.data(),
                nullptr, // step: defaults to 1,1,1 
                nullptr, // stride: default to row-major convention
                GDALExtendedDataType::Create(GDT_Float64),
                &values[0]);

			
	return values;
*/

  
#else  


bool SpatRaster::constructFromFileMulti(std::string fname, std::string sub, std::vector<size_t> xyz) {
	setError("multidim is not supported by GDAL < 3.1");
	return false;
}


#endif

