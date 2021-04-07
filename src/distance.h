// Copyright (c) 2018-2020  Robert J. Hijmans
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

// distance
double distance_plane(const double &x1, const double &y1, const double &x2, const double &y2);
std::vector<double> distance_plane(std::vector<double> &x1, std::vector<double> &y1, std::vector<double> &x2, std::vector<double> &y2);
std::vector<double> distance_plane_vd(std::vector<double> &x1, std::vector<double> &y1, double x2, double y2);

double distance_haversine(double lon1, double lat1, double lon2, double lat2);

double distance_lonlat(const double &lon1, const double &lat1, const double &lon2, const double &lat2);
std::vector<double> distance_lonlat(std::vector<double> &lon1, std::vector<double> &lat1, std::vector<double> &lon2, std::vector<double> &lat2) ;
std::vector<double> distance_lonlat_vd(std::vector<double> &lon1, std::vector<double> &lat1, double lon2, double lat2) ;


// direction
double direction_lonlat(double lon1, double lat1, double lon2, double lat2, bool degrees);
std::vector<double> direction_lonlat(std::vector<double> lon1, std::vector<double> lat1, std::vector<double> lon2, std::vector<double> lat2, bool degrees);
std::vector<double> directionToNearest_lonlat(std::vector<double> lon1, std::vector<double> lat1, std::vector<double> lon2, std::vector<double> lat2, bool degrees, bool from);

double direction_plane(double x1, double y1, double x2, double y2, bool degrees);
std::vector<double> direction_plane(std::vector<double> x1, std::vector<double> y1, std::vector<double> x2, std::vector<double> y2, bool degrees);
std::vector<double> directionToNearest_plane(std::vector<double> x1, std::vector<double> y1, std::vector<double> x2, std::vector<double> y2, bool degrees, bool from);

// destination
std::vector<double> destpoint_lonlat(double longitude, double latitude, double  bearing, double distance);

std::vector<std::vector<double> > destpoint_lonlat(std::vector<double>& longitude, std::vector<double>& latitude, std::vector<double>& bearing, std::vector<double>& distance);

std::vector<std::vector<double> > destpoint_lonlat(double& longitude, double& latitude, std::vector<double>& bearing, double& distance);

std::vector<double> destpoint_plane(double x, double y, double bearing, double distance);
std::vector<std::vector<double> > destpoint_plane(std::vector<double>  x, std::vector<double>  y, std::vector<double>  bearing, std::vector<double>  distance);

double toRad(double &deg);
