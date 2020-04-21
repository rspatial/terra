#include <gdal_priv.h> // GDALDriver

void set_error_handler(void);
void unset_error_handler(void);
std::string wkt_from_spatial_reference(const OGRSpatialReference *srs);
std::string prj_from_spatial_reference(const OGRSpatialReference *srs);
std::vector<std::string> srefs_from_string(std::string input);

