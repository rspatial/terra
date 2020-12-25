/*
#include "gdal_priv.h"
#include "gdal.h"

#include "spatRaster.h"


std::vector<double> SpatRaster::readmulti(std::string filename, std::string var) {

	SpatRaster out;
    auto poDataset = std::unique_ptr<GDALDataset>(
        GDALDataset::Open(filename.c_str(), GDAL_OF_MULTIDIM_RASTER ));

    if (!poDataset ) {
		setError("not a good dataset");
		return();
    }
    auto poRootGroup = poDataset->GetRootGroup();
    if( !poRootGroup ) {
		setError("no roots");
		return();
    }
    auto poVar = poRootGroup->OpenMDArray(var.c_str());
    if( !poVar )   {
		setError("cannot open this var");
		return();
    }

    size_t nValues = 1;
    std::vector<size_t> anCount;

    for( const auto poDim: poVar->GetDimensions() ) {
        anCount.push_back(static_cast<size_t>(poDim->GetSize()));
        nValues *= anCount.back();
    }

    std::vector<double> values(nValues);
    poVar->Read(std::vector<GUInt64>{0,0,0}.data(),
                anCount.data(),
                nullptr, // step: defaults to 1,1,1 
                nullptr, // stride: default to row-major convention
                GDALExtendedDataType::Create(GDT_Float64),
                &values[0]);
			
	return values;
}

*/