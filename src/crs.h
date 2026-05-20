// Copyright (c) 2018-2023  Robert J. Hijmans
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



SpatDataFrame get_proj_pipelines(std::string source_crs, std::string target_crs,
		std::string authority, std::vector<double> AOI, std::string use,
		std::string grid_availability, double desired_accuracy,
		bool strict_containment, bool axis_order_authority_compliant);
		


bool can_transform(std::string fromCRS, std::string toCRS);
SpatMessages transform_coordinates(std::vector<double> &x, std::vector<double> &y, std::string fromCRS, std::string toCRS);

// Call at the start of any operation that performs coordinate
// transformations, before any GDAL/PROJ work, to clear stale noise flags
// from previous operations.
void proj_noise_reset();

// Drain any "noisy" PROJ messages that the GDAL error handler has collapsed
// (currently: CDN download failures and cache.db lock failures) into the
// given SpatMessages object as a single warning each, and reset the flags.
// Call this at the end of any operation that performs coordinate transformations

void proj_noise_drain(SpatMessages &m);

// Reset PROJ noise flags on construction and drain any
// collected noise into the target SpatMessages on destruction. Useful for
// functions with many return paths.
struct ProjNoiseScope {
	SpatMessages *m_target;
	explicit ProjNoiseScope(SpatMessages &target) : m_target(&target) {
		proj_noise_reset();
	}
	~ProjNoiseScope() {
		if (m_target) proj_noise_drain(*m_target);
	}
	ProjNoiseScope(const ProjNoiseScope&) = delete;
	ProjNoiseScope& operator=(const ProjNoiseScope&) = delete;
};
bool wkt_from_spatial_reference(const OGRSpatialReference *srs, std::string &wkt, std::string &msg);
bool prj_from_spatial_reference(const OGRSpatialReference *srs, std::string &prj, std::string &msg);
//std::vector<std::string> srefs_from_string(std::string input);
bool wkt_from_string(std::string input, std::string& wkt, std::string& msg);
bool is_ogr_error(OGRErr err, std::string &msg);
void geo_ellipsoid_from_wkt(const std::string &wkt, double &a, double &f);

//#endif
