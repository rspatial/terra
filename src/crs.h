#include "ogr_spatialref.h"

SpatMessages transform_coordinates(std::vector<double> &x, std::vector<double> &y, std::string fromCRS, std::string toCRS);
std::string wkt_from_spatial_reference(const OGRSpatialReference *srs);
std::string prj_from_spatial_reference(const OGRSpatialReference *srs);
std::vector<std::string> srefs_from_string(std::string input);

