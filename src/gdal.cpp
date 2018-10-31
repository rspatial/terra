/*
struct Format {
  enum Enum {
    NONE = 0,
    RASTER,
    VECTOR,
  };
};

Format::Enum GetFileFormat( string const & path ) {
	Format::Enum format = Format::NONE;
	unsigned int openFlags = GA_ReadOnly; // | GDAL_OF_VECTOR;
	GDALDataset * pDataset = static_cast<GDALDataset *>( GDALOpenEx( path.c_str(), openFlags, NULL, NULL, NULL ) );
	if ( pDataset )  {
		int rasterCount = pDataset->GetRasterCount();
		int layerCount = pDataset->GetLayerCount();
		if ( rasterCount > 0 && layerCount == 0 ) {
			format = Format::RASTER;
		} else if ( rasterCount == 0 && layerCount > 0 ) {
			format = Format::VECTOR;
		}
	}
	GDALClose( pDataset );
	return format;
}
*/