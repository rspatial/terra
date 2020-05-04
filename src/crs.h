#ifdef useGDAL
#include "ogr_spatialref.h"

SpatMessages transform_coordinates(std::vector<double> &x, std::vector<double> &y, std::string fromCRS, std::string toCRS);
bool wkt_from_spatial_reference(const OGRSpatialReference *srs, std::string &wkt, std::string &msg);
bool prj_from_spatial_reference(const OGRSpatialReference *srs, std::string &prj, std::string &msg);
//std::vector<std::string> srefs_from_string(std::string input);
bool wkt_from_string(std::string input, std::string& wkt, std::string& msg);
bool is_ogr_error(OGRErr err, std::string &msg);

#endif
