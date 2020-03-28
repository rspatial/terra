//g++ gdalinfo.cpp -lgdal -I/usr/include/gdal -o ginfo
#include "gdal_priv.h"
#include "gdal_frmts.h"
#include "cpl_conv.h" // for CPLMalloc()
#include <iostream>

int main(int argc, char *argv[]) {

    if(argc < 2) {
        std::string path = (std::string)argv[0];
        std::string basename = path.substr(path.find_last_of("/\\") + 1);
        std::string usage = "Usage: " + basename + "method args";
        std::cout << usage << std::endl;
        return 1;
    } else {
        GDALDataset  *poDataset;
        GDALAllRegister();
	const char* pszFilename = argv[1];

        poDataset = (GDALDataset *) GDALOpen( pszFilename, GA_ReadOnly );
        if( poDataset == NULL ) {
           std::cout << "fail" << std::endl;
        } else {
          	double adfGeoTransform[6];
		printf( "Driver: %s/%s\n",
        	poDataset->GetDriver()->GetDescription(),
       		 poDataset->GetDriver()->GetMetadataItem( GDAL_DMD_LONGNAME ) );
		printf( "Size is %dx%dx%d\n",
	        poDataset->GetRasterXSize(), poDataset->GetRasterYSize(),
	        poDataset->GetRasterCount() );
		if( poDataset->GetProjectionRef()  != NULL )
		    printf( "Projection is `%s'\n", poDataset->GetProjectionRef() );
		if( poDataset->GetGeoTransform( adfGeoTransform ) == CE_None )
		{
		    printf( "Origin = (%.6f,%.6f)\n",
		            adfGeoTransform[0], adfGeoTransform[3] );
		    printf( "Pixel Size = (%.6f,%.6f)\n",
	            adfGeoTransform[1], adfGeoTransform[5] );
		}


        }

    }
    return 0;
}

