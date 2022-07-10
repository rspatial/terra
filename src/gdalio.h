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

bool getGDALDataType(std::string datatype, GDALDataType &gdt);
std::string gdalinfo(std::string filename, std::vector<std::string> options, std::vector<std::string> openopts);
std::vector<std::vector<std::string>> sdinfo(std::string fname);
std::vector<std::string> get_metadata(std::string filename);
std::vector<std::string> get_metadata_sds(std::string filename);
std::vector<std::vector<std::string>> parse_metadata_sds(std::vector<std::string> meta);
void getGDALdriver(std::string &filename, std::string &driver);
bool getNAvalue(GDALDataType gdt, double & naval);
GDALDataset* openGDAL(std::string filename, unsigned OpenFlag, std::vector<std::string> allowed_drivers, std::vector<std::string> open_options);
char ** set_GDAL_options(std::string driver, double diskNeeded, bool writeRGB, std::vector<std::string> gdal_options);

