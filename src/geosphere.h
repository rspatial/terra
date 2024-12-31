// Copyright (c) 2018-2025  Robert J. Hijmans
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


//double distance_cos(const double &lon1, const double &lat1, const double &lon2, const double &lat2, const double &r = 6378137);

inline void deg2rad(std::vector<double> &x) {
	const double f = 0.0174532925199433;
	for (double& d : x) d *= f;
}

inline void deg2rad(double &x) {
	const double f = 0.0174532925199433;
	x *= f;
}


inline double distance_cos(double lon1, double lat1, double lon2, double lat2) {
	const double r = 6378137.;
	return r * acos((sin(lat1) * sin(lat2) + cos(lat1) * cos(lat2) * cos(lon1-lon2)));
}



inline double distance_hav(double lon1, double lat1, double lon2, double lat2) {
	const double r = 6378137.;
	double dLat = lat2-lat1;
	double dLon = lon2-lon1;
	double a = sin(dLat/2.) * sin(dLat/2.) + cos(lat1) * cos(lat2) * sin(dLon/2.) * sin(dLon/2.);
	return 2. * atan2(sqrt(a), sqrt(1. - a)) * r;
}

double distance_geo(double lon1, double lat1, double lon2, double lat2);


double direction_cos(double& lon1, double& lat1, double& lon2, double& lat2);

double dist2segment_cos(double plon, double plat, double lon1, double lat1, double lon2, double lat2, double r=6378137.);
double dist2segment_geo(double plon, double plat, double lon1, double lat1, double lon2, double lat2, double notused=0.);
