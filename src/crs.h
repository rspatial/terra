// Copyright (c) 2018-2022  Robert J. Hijmans
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

//#ifdef useGDAL
#include "ogr_spatialref.h"

bool can_transform(std::string fromCRS, std::string toCRS);
SpatMessages transform_coordinates(std::vector<double> &x, std::vector<double> &y, std::string fromCRS, std::string toCRS);
bool wkt_from_spatial_reference(const OGRSpatialReference *srs, std::string &wkt, std::string &msg);
bool prj_from_spatial_reference(const OGRSpatialReference *srs, std::string &prj, std::string &msg);
//std::vector<std::string> srefs_from_string(std::string input);
bool wkt_from_string(std::string input, std::string& wkt, std::string& msg);
bool is_ogr_error(OGRErr err, std::string &msg);

//#endif
