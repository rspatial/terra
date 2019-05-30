// Copyright (c) 2018-2019  Robert J. Hijmans
//
// This file is part of the "spat" library.
//
// spat is free software: you can redistribute it and/or modify it
// under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 2 of the License, or
// (at your option) any later version.
//
// spat is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with spat. If not, see <http://www.gnu.org/licenses/>.

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