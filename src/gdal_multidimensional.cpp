#include "spatRaster.h"

#if GDAL_VERSION_MAJOR >= 3 && GDAL_VERSION_MINOR >= 1

#include "proj.h"


//#include "cpl_conv.h" // for CPLMalloc()
//#include "cpl_string.h"
#include "ogr_spatialref.h"

#include "gdal_priv.h"
#include "gdal.h"
#include "crs.h"



bool SpatRaster::constructFromFileMulti(std::string fname, std::string sub, std::vector<size_t> xyz) {

    GDALAllRegister();
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

	//std::string var = "HWAM";
    auto poVar = poRootGroup->OpenMDArray(sub.c_str());
    if( !poVar )   {
		setError("cannot find: " + sub);
		return false;
    }

//	const auto sref:
	std::shared_ptr<OGRSpatialReference> srs = poVar->GetSpatialRef();
	if (srs != NULL) {
		char *cp;
		std::string msg;
		OGRErr err = srs->exportToProj4(&cp);
		if (is_ogr_error(err, msg)) {
			Rcpp::Rcout << msg << std::endl;
		} else {
			std::string prj = std::string(cp);
			Rcpp::Rcout << prj << std::endl;
		}
		CPLFree(cp);
	} else {
		Rcpp::Rcout << "srs is null" << std::endl;		
	}
	
    std::vector<std::string> dimname;
	SpatRasterSource s;
	s.multidim = true;

    for( const auto poDim: poVar->GetDimensions() ) {
        s.mdims.push_back(static_cast<size_t>(poDim->GetSize()));
        dimname.push_back(static_cast<std::string>(poDim->GetName()));
    }
	
	//std::vector<size_t> xyz = {0,1,2};
	s.ncol = s.mdims[xyz[0]];
	s.nrow = s.mdims[xyz[1]];
	s.nlyr = s.mdims[xyz[2]];
	setSource(s);
	
	for (size_t i=0; i<s.mdims.size(); i++){
		Rcpp::Rcout << s.mdims[i] << " " << dimname[i] << std::endl;
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
	setError("multidim is not supported by GDAL < 3.1")
	return false;
}


#endif

