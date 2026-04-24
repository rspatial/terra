// Copyright (c) 2018-2026  Robert J. Hijmans
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

#include "spatVector.h"
#include "geodesic.h"
#include "recycle.h"
#include "geosphere.h"

#ifndef M_PI
#define M_PI  (3.1415926535897932384626433)
#endif

#ifndef M_2PI
#define M_2PI (M_PI * 2.0)
#endif

#ifndef M_hPI
#define M_hPI (M_PI / 2.0)
#endif

#ifndef WGS84_a
#define WGS84_a 6378137.0
#endif

#ifndef WGS84_f
#define WGS84_f 1/298.257223563
#endif


inline void normLon(double &lon) {
	lon = fmod(lon + 180, 360.) - 180;
}

inline void normLonRad(double &lon) {
	lon = fmod(lon + M_PI, M_2PI) - M_PI;
}




inline double get_sign(const double &x) {
	return (x > 0.0) ? 1.0 : (x < 0.0) ? -1.0 : 0;
}


double distance_geo(double lon1, double lat1, double lon2, double lat2) {
	double s12, azi1, azi2;
	struct geod_geodesic g;
	geod_init(&g, WGS84_a, WGS84_f);
	geod_inverse(&g, lat1, lon1, lat2, lon2, &s12, &azi1, &azi2);
  	return s12;
}

inline double distance_cos_r(double lon1, double lat1, double lon2, double lat2, double r = 6378137.) {
	return r * acos((sin(lat1) * sin(lat2) + cos(lat1) * cos(lat2) * cos(lon1-lon2)));
}

inline double distance_hav_r(double lon1, double lat1, double lon2, double lat2, const double r = 6378137.) {
	double dLat = lat2-lat1;
	double dLon = lon2-lon1;
	double a = sin(dLat/2.) * sin(dLat/2.) + cos(lat1) * cos(lat2) * sin(dLon/2.) * sin(dLon/2.);
	return 2. * atan2(sqrt(a), sqrt(1. - a)) * 6378137.0;
}


/*
double distance_cosdeg(double lon1, double lat1, double lon2, double lat2, double r = 6378137.) {
	deg2rad(lon1);
	deg2rad(lon2);
	deg2rad(lat1);
	deg2rad(lat2);
	return r * acos((sin(lat1) * sin(lat2) + cos(lat1) * cos(lat2) * cos(lon1-lon2)));
}
*/


void dest_geo(double slon, double slat, double sazi, double dist, double &dlon, double &dlat, double &dazi) {
	struct geod_geodesic g;
	geod_init(&g, WGS84_a, WGS84_f);
	geod_direct(&g, slat, slon, sazi, dist, &dlat, &dlon, &dazi);
}


double direction_geo(double lon1, double lat1, double lon2, double lat2) {
	double s12, azi1, azi2;
	struct geod_geodesic g;
	geod_init(&g, WGS84_a, WGS84_f);
	geod_inverse(&g, lat1, lon1, lat2, lon2, &s12, &azi1, &azi2);
	return( azi1) ;
}


double direction_cos(double& lon1, double& lat1, double& lon2, double& lat2) {
	if ((lon1 == lon2) && (lat1 == lat2)) return 0; // NAN?
	double dLon = lon2 - lon1;
	double y = sin(dLon) * cos(lat2);
	double x = cos(lat1) * sin(lat2) - sin(lat1) * cos(lat2) * cos(dLon);
	double azm = atan2(y, x);
	// normalize to [0, 2*pi), consistent with direction_plane
	if (azm < 0) azm += M_2PI;
	return azm;
}



double dist2track_geo(double lon1, double lat1, double lon2, double lat2, double plon, double plat, bool sign, double r=6378137) {
	double a = 1;
	double f = 0;
	struct geod_geodesic geod;
	geod_init(&geod, a, f);
	double d, b2, b3, azi;
	geod_inverse(&geod, lat1, lon1, lat2, lon2, &d, &b2, &azi);
	geod_inverse(&geod, lat1, lon1, plat, plon, &d, &b3, &azi);
	double toRad = M_PI / 180.;
	b2 *= toRad;
	b3 *= toRad;
	double xtr = asin(sin(b3-b2) * sin(d)) * r;
	return sign ? xtr : fabs(xtr);
}


inline double dist2track_cos(double lon1, double lat1, double lon2, double lat2, double plon, double plat, bool sign, double r=6378137) {
	double b2 = direction_cos(lon1, lat1, lon2, lat2);
	double b3 = direction_cos(lon1, lat1, plon, plat);
	double d = distance_cos_r(lon1, lat1, plon, plat, 1);
	double xtr = asin(sin(b3-b2) * sin(d)) * r;
	return sign ? xtr : fabs(xtr);
}


inline double dist2track_hav(double lon1, double lat1, double lon2, double lat2, double plon, double plat, bool sign, double r=6378137) {
	double b2 = direction_cos(lon1, lat1, lon2, lat2);
	double b3 = direction_cos(lon1, lat1, plon, plat);
	double d = distance_hav_r(lon1, lat1, plon, plat, 1);
	double xtr = asin(sin(b3-b2) * sin(d)) * r;
	return sign ? xtr : fabs(xtr);
}



// signed along-track distance
double alongTrackDistance_geo(double lon1, double lat1, double lon2, double lat2, double plon, double plat, double r=6378137) {
	double a = 1;
	double f = 0;
	struct geod_geodesic geod;
	geod_init(&geod, a, f);
	double d, b2, b3, azi;
	geod_inverse(&geod, lat1, lon1, lat2, lon2, &d, &b2, &azi);
	geod_inverse(&geod, lat1, lon1, plat, plon, &d, &b3, &azi);
	deg2rad(b2);
	deg2rad(b3);
	double xtr = asin(sin(b3-b2) * sin(d));
	double bsign = get_sign(cos(b2-b3));
	double angle = cos(d) / cos(xtr);
	angle = angle > 1 ? 1 : angle < -1 ? -1 : angle;
	return bsign * acos(angle) * r;
}


double alongTrackDistance_cos(double lon1, double lat1, double lon2, double lat2, double plon, double plat, double r=6378137) {

	double tc = direction_cos(lon1, lat1, lon2, lat2); // * toRad
	double tcp = direction_cos(lon1, lat1, plon, plat); // * toRad
    double dp = distance_cos_r(lon1, lat1, plon, plat, 1);
	double xtr = asin(sin(tcp-tc) * sin(dp));

// +1/-1 for ahead/behind [lat1,lon1]
	double bearing = get_sign(cos(tc - tcp));
	double angle = cos(dp) / cos(xtr);

// Fixing limits for the angle between [-1, 1] to avoid NaNs from acos
	angle = angle > 1 ? 1 : angle < -1 ? -1 : angle;
	return bearing * acos(angle) * r;
}



double alongTrackDistance_hav(double lon1, double lat1, double lon2, double lat2, double plon, double plat, double r=6378137) {

	double tc = direction_cos(lon1, lat1, lon2, lat2); // * toRad
	double tcp = direction_cos(lon1, lat1, plon, plat); // * toRad
    double dp = distance_hav_r(lon1, lat1, plon, plat, 1);
	double xtr = asin(sin(tcp-tc) * sin(dp));

// +1/-1 for ahead/behind [lat1,lon1]
	double bearing = get_sign(cos(tc - tcp));
	double angle = cos(dp) / cos(xtr);

// Fixing limits for the angle between [-1, 1] to avoid NaNs from acos
	angle = angle > 1 ? 1 : angle < -1 ? -1 : angle;
	return bearing * acos(angle) * r;
}


// the alongTrackDistance is the length of the path along the great circle to the point of intersection
// there are two, depending on which node you start
// we want to use the min, but the max needs to be < segment length

double dist2segment_geo(double plon, double plat, double lon1, double lat1, double lon2, double lat2, double notused) {

	double seglength = distance_geo(lon1, lat1, lon2, lat2);
	double trackdist1 = alongTrackDistance_geo(lon1, lat1, lon2, lat2, plon, plat);
	double trackdist2 = alongTrackDistance_geo(lon2, lat2, lon1, lat1, plon, plat);
	if ((trackdist1 < 0) || (trackdist2 < 0) ||
		(trackdist1 > seglength) || (trackdist2 > seglength)) {
		double d1 = distance_geo(lon1, lat1, plon, plat);
		double d2 = distance_geo(lon2, lat2, plon, plat);
		return d1 < d2 ? d1 : d2;
	}
	return dist2track_geo(lon1, lat1, lon2, lat2, plon, plat, false);
}


double dist2segment_cos(double plon, double plat, double lon1, double lat1, double lon2, double lat2, double r) {
	double seglength = distance_cos_r(lon1, lat1, lon2, lat2, r);
	double trackdist1 = alongTrackDistance_cos(lon1, lat1, lon2, lat2, plon, plat, r);
	double trackdist2 = alongTrackDistance_cos(lon2, lat2, lon1, lat1, plon, plat, r);
	if ((trackdist1 < 0) || (trackdist2 < 0) ||
		(trackdist1 > seglength) || (trackdist2 > seglength)) {
		double d1 = distance_cos_r(lon1, lat1, plon, plat, r);
		double d2 = distance_cos_r(lon2, lat2, plon, plat, r);
		return d1 < d2 ? d1 : d2;
	}
	return dist2track_cos(lon1, lat1, lon2, lat2, plon, plat, false, r);
}

double dist2segment_hav(double plon, double plat, double lon1, double lat1, double lon2, double lat2, double r) {
	double seglength = distance_hav_r(lon1, lat1, lon2, lat2, r);
	double trackdist1 = alongTrackDistance_hav(lon1, lat1, lon2, lat2, plon, plat, r);
	double trackdist2 = alongTrackDistance_hav(lon2, lat2, lon1, lat1, plon, plat, r);
	if ((trackdist1 < 0) || (trackdist2 < 0) ||
		(trackdist1 > seglength) || (trackdist2 > seglength)) {
		double d1 = distance_hav_r(lon1, lat1, plon, plat, r);
		double d2 = distance_hav_r(lon2, lat2, plon, plat, r);
		return d1 < d2 ? d1 : d2;
	}
	return dist2track_hav(lon1, lat1, lon2, lat2, plon, plat, false, r);
}


// [[Rcpp::export]]
double dist2segmentPoint_geo(double plon, double plat, double lon1, double lat1, double lon2, double lat2, double &ilon, double &ilat) {

	double seglength = distance_geo(lon1, lat1, lon2, lat2);
	double trackdist1 = alongTrackDistance_geo(lon1, lat1, lon2, lat2, plon, plat);
	double trackdist2 = alongTrackDistance_geo(lon2, lat2, lon1, lat1, plon, plat);
	if ((trackdist1 < 0) || (trackdist2 < 0) ||
		(trackdist1 > seglength) || (trackdist2 > seglength)) {
		double d1 = distance_geo(lon1, lat1, plon, plat);
		double d2 = distance_geo(lon2, lat2, plon, plat);
		if (d1 < d2) {
			ilon = lon1;
			ilat = lat1;
			return d1;
		} else {
			ilon = lon2;
			ilat = lat2;
			return d2;
		}
	}
	double azi;
	double crossd = dist2track_geo(lon1, lat1, lon2, lat2, plon, plat, false);
	if (trackdist1 < trackdist2) {
		double bear = direction_geo(lon1, lat1, lon2, lat2);
		dest_geo(lon1, lat1, bear, trackdist1, ilon, ilat, azi);
	} else {
		double bear = direction_geo(lon2, lat2, lon1, lat1);
		dest_geo(lon2, lat2, bear, trackdist2, ilon, ilat, azi);
	}
	return crossd;
}



// [[Rcpp::export(name = "intermediate")]]
std::vector<std::vector<double>> intermediate(double lon1, double lat1, double lon2, double lat2, int n, double distance) {
	double a = 6378137.0;
	double f = 1/298.257223563;
	struct geod_geodesic geod;
	geod_init(&geod, a, f);
	double d, azi1, azi2;

	std::vector<std::vector<double>> out(2);
	if (n <= 0) {
		if (distance <= 0) {
			out[0] = {lon1, lon2};
			out[1] = {lon1, lon2};
			return out;
		} else {
			geod_inverse(&geod, lat1, lon1, lat2, lon2, &d, &azi1, &azi2);
			n = std::round(d / distance);
			if (n < 2) {
				out[0] = {lon1, lon2};
				out[1] = {lon1, lon2};
				return out;
			}
			distance = d / n;
		}
	} else if (n == 1) {
		out[0] = {lon1, lon2};
		out[1] = {lon1, lon2};
		return out;
	} else {
		geod_inverse(&geod, lat1, lon1, lat2, lon2, &d, &azi1, &azi2);
		//distance = d / n;
	}
	out[0].resize(n+1);
	out[1].resize(n+1);
	out[0][0] = lon1;
	out[1][0] = lat1;
	for (int i=1; i<n; i++) {
		geod_direct(&geod, lat1, lon1, azi1, distance*i, &out[1][i], &out[0][i], &azi2);
	}
	out[0][n] = lon2;
	out[1][n] = lat2;
	return out;
}




std::vector<bool> antipodal(std::vector<double> lon1, std::vector<double> lat1, std::vector<double> lon2, std::vector<double> lat2, double tol=1e-9) {
	recycle(lon1, lon2);
	recycle(lat1, lat2);
	std::vector<bool> out;
	out.reserve(lon1.size());
	double Pi180 = M_PI / 180.;
	for (size_t i=0; i<lon1.size(); i++){
		normLon(lon1[i]);
		normLon(lon2[i]);
		double diflon = fabs(lon1[i] - lon2[i]);
		double diflat = fabs(lat1[i] + lat2[i]);
		out.push_back(
			(diflat < tol) && ((cos(lat2[i] * Pi180) * fabs(fmod(diflon, 360.) - 180)) < tol)
		);
	}
	return out;
}


void antipodes(std::vector<double> &lon, std::vector<double> &lat) {
	size_t n=lon.size();
	for (size_t i=0; i<n; i++) {
		lon[i] = lon[i] + 180;
		normLon(lon[i]);
		lat[i] = -lat[i];
	}
}



// Spherical point-in-polygon test by ray casting along each test point's meridian, 
// counting how many polygon edges the ray (from the point to the north
// pole) crosses. Odd = inside. Holes invert the result.
//
// Implementation notes:
//  - For each edge V1->V2, compute the longitude offsets d1, d2 of V1
//    and V2 relative to the test point's longitude, normalized to
//    (-pi, pi]. Anchor d2 within pi of d1 to follow the shorter arc.
//  - The arc crosses the meridian iff (d1 < 0) XOR (e2 < 0), where
//    e2 = d1 + delta_short. This half-open convention prevents
//    double-counting at vertices that lie exactly on the meridian.
//  - The crossing latitude is found from the great-circle equation
//        N . P = 0   where N = V1 x V2 (Cartesian on unit sphere)
//    giving tan(lat) = -(Nx*cos(plon) + Ny*sin(plon)) / Nz.
//  - Edges whose great-circle goes through both poles (Nz ~ 0) are
//    skipped: their intersection with the meridian is 0-dimensional
//    measure and the topological count is unaffected.
//  - Edges of zero length are skipped.

inline bool pip_geo_edge_crosses_north(double lon1, double lat1, double lon2, double lat2, double plon, double plat) {
	static const double D2R = M_PI / 180.0;
	static const double TWOPI = 2.0 * M_PI;
	if (lon1 == lon2 && lat1 == lat2) return false;

	double rl1 = lon1 * D2R, rl2 = lon2 * D2R;
	double rp1 = lat1 * D2R, rp2 = lat2 * D2R;
	double rpl = plon * D2R, rpt = plat * D2R;

	double d1 = rl1 - rpl;
	while (d1 >  M_PI) d1 -= TWOPI;
	while (d1 <= -M_PI) d1 += TWOPI;
	double d2 = rl2 - rpl;
	while (d2 >  M_PI) d2 -= TWOPI;
	while (d2 <= -M_PI) d2 += TWOPI;

	double delta = d2 - d1;
	while (delta >  M_PI) delta -= TWOPI;
	while (delta <= -M_PI) delta += TWOPI;
	double e2 = d1 + delta;

	bool d1neg = (d1 < 0.0);
	bool e2neg = (e2 < 0.0);
	if (d1neg == e2neg) return false;

	double cl1 = cos(rp1), sl1 = sin(rp1);
	double cl2 = cos(rp2), sl2 = sin(rp2);
	double x1 = cl1 * cos(rl1), y1 = cl1 * sin(rl1), z1 = sl1;
	double x2 = cl2 * cos(rl2), y2 = cl2 * sin(rl2), z2 = sl2;
	double Nx = y1 * z2 - z1 * y2;
	double Ny = z1 * x2 - x1 * z2;
	double Nz = x1 * y2 - y1 * x2;
	if (std::fabs(Nz) < 1e-15) return false;
	double cp = cos(rpl), sp = sin(rpl);
	double lat_cross = atan(-(Nx * cp + Ny * sp) / Nz);
	return lat_cross > rpt;
}


bool pip_geo_in_ring(double plon, double plat, const std::vector<double> &rx, const std::vector<double> &ry) {
	int crossings = 0;
	size_t n = rx.size();
	if (n < 3) return false;
	// Treat ring as closed even if not explicitly so.
	bool closed = (rx[0] == rx[n-1]) && (ry[0] == ry[n-1]);
	size_t last = closed ? n - 1 : n;
	for (size_t i = 0; i < last; i++) {
		size_t j = (i + 1) % n;
		if (!closed && j == 0) break;
		if (pip_geo_edge_crosses_north(rx[i], ry[i], rx[j], ry[j], plon, plat)) {
			crossings++;
		}
	}
	return (crossings & 1) == 1;
}


// Prepare per-edge invariants for fast point-in-polygon queries.
//
// For each edge V1->V2 we precompute the unit-sphere cross product
// N = V1 x V2; the great-circle through V1 and V2 has equation N.P = 0,
// so the meridian-crossing latitude of the arc at longitude plon is
//
//     tan(lat_cross) = -(Nx*cos(plon) + Ny*sin(plon)) / Nz
//
// which is now a couple of multiplies + an atan per edge per point, with
// no trig of the edge endpoints required at query time.
//
// We also cache:
//   * a safe upper bound on the latitude reached by any arc
//     (max_arc_lat), used to short-circuit rings whose entire arc lies
//     south of the test point's latitude;
//   * whether the ring topologically encircles a pole (wraps_pole) and
//     whether the north pole sits inside the ring (flip_north). The
//     ray-cast goes from the test point north to the pole; if that ray
//     terminates inside the ring, the parity of edge crossings has the
//     opposite meaning from the usual "odd = inside" convention.
void prepare_pip_geo_ring(const std::vector<double> &rx,
                          const std::vector<double> &ry,
                          PipGeoRing &out) {
	static const double D2R = M_PI / 180.0;
	static const double R2D = 180.0 / M_PI;

	out.lon.clear();
	out.lat.clear();
	out.Nx.clear();
	out.Ny.clear();
	out.Nz.clear();
	out.max_arc_lat = -90.0;
	out.wraps_pole = false;
	out.flip_north = false;

	size_t n = rx.size();
	if (n < 3) return;

	// Close the ring if the caller didn't.
	bool closed = (rx[0] == rx[n-1]) && (ry[0] == ry[n-1]);
	out.lon = rx;
	out.lat = ry;
	if (!closed) {
		out.lon.push_back(rx[0]);
		out.lat.push_back(ry[0]);
	}
	n = out.lon.size();
	size_t nedges = n - 1;
	out.Nx.assign(nedges, 0.0);
	out.Ny.assign(nedges, 0.0);
	out.Nz.assign(nedges, 0.0);

	double total_dlon = 0.0;
	double sum_lat = 0.0;
	double max_lat = -90.0;

	for (size_t i = 0; i < nedges; i++) {
		double lon1 = out.lon[i];
		double lat1 = out.lat[i];
		double lon2 = out.lon[i+1];
		double lat2 = out.lat[i+1];

		sum_lat += lat1;

		// Signed shortest-arc dlon in degrees, summed for pole-wrap detection.
		double dlon_deg = lon2 - lon1;
		while (dlon_deg >  180.0) dlon_deg -= 360.0;
		while (dlon_deg <= -180.0) dlon_deg += 360.0;
		total_dlon += dlon_deg;

		if (lon1 == lon2 && lat1 == lat2) {
			// Degenerate edge; leave N = 0 so it is skipped at query time.
			if (lat1 > max_lat) max_lat = lat1;
			continue;
		}

		double rl1 = lon1 * D2R, rl2 = lon2 * D2R;
		double rp1 = lat1 * D2R, rp2 = lat2 * D2R;
		double cl1 = cos(rp1), sl1 = sin(rp1);
		double cl2 = cos(rp2), sl2 = sin(rp2);
		double x1 = cl1 * cos(rl1), y1 = cl1 * sin(rl1), z1 = sl1;
		double x2 = cl2 * cos(rl2), y2 = cl2 * sin(rl2), z2 = sl2;
		double Nx = y1 * z2 - z1 * y2;
		double Ny = z1 * x2 - x1 * z2;
		double Nz = x1 * y2 - y1 * x2;
		out.Nx[i] = Nx;
		out.Ny[i] = Ny;
		out.Nz[i] = Nz;

		// Upper bound on the latitude reached by this arc.
		// The full great circle's max |sin(lat)| is sqrt(1 - (Nz/|N|)^2);
		// for short arcs the actual extremum lies between the endpoints,
		// so endpoint max is a safe (and usually tight) bound. For longer
		// arcs we need to use the great-circle bulge upper bound.
		double endpoint_max = std::max(lat1, lat2);
		double Nmag2 = Nx*Nx + Ny*Ny + Nz*Nz;
		double arc_max = endpoint_max;
		if (Nmag2 > 1e-30) {
			double s = 1.0 - (Nz*Nz) / Nmag2;
			if (s < 0.0) s = 0.0;
			double gc_max = asin(sqrt(s)) * R2D;
			// Use the gc bulge bound only when the arc is long enough
			// for the bulge to actually exceed the endpoints. A safe
			// proxy: if endpoint span in lon is small, no meaningful
			// bulge above endpoints. Otherwise fall back to gc_max.
			double dlon_abs = std::fabs(dlon_deg);
			if (dlon_abs > 1.0 && gc_max > endpoint_max) {
				arc_max = gc_max;
			}
		}
		if (arc_max > max_lat) max_lat = arc_max;
	}
	out.max_arc_lat = max_lat;

	// A ring whose signed dlons sum to ~+/-360 either encircles a pole
	// or wraps all the way around in longitude. In practice (for our
	// callers) such rings encircle a pole; pick which pole by the sign
	// of the average vertex latitude. This is a robust heuristic for
	// caps/skirts that lie predominantly in one hemisphere.
	double avg_lat = sum_lat / static_cast<double>(nedges);
	out.wraps_pole = std::fabs(std::fabs(total_dlon) - 360.0) < 1.0;
	out.flip_north = out.wraps_pole && (avg_lat > 0.0);
	if (out.wraps_pole) {
		// A pole-enclosing ring has no useful lat upper bound.
		out.max_arc_lat = 90.0;
	}
}


bool pip_geo_in_ring_prepared(double plon_deg, double plat_deg,
                              double cos_plon_rad, double sin_plon_rad,
                              double plat_rad,
                              const PipGeoRing &ring) {
	size_t nedges = ring.Nx.size();
	if (nedges == 0) return false;

	// Cheap ring-level prune: if the ray (going north from plat to the
	// north pole) starts above every edge's max latitude, no edge can be
	// crossed. Skip the entire ring. Disabled when the ring encloses a
	// pole (no useful lat bound).
	if (plat_deg > ring.max_arc_lat) {
		return ring.flip_north;
	}

	int crossings = 0;
	for (size_t i = 0; i < nedges; i++) {
		double Nz = ring.Nz[i];
		// Skip degenerate / pole-to-pole great circles.
		if (Nz == 0.0) continue;

		// Straddle test in degrees: does the shorter arc from V1 to V2
		// cross the meridian of the test point?
		double d1 = ring.lon[i] - plon_deg;
		while (d1 >  180.0) d1 -= 360.0;
		while (d1 <= -180.0) d1 += 360.0;
		double d2 = ring.lon[i+1] - plon_deg;
		while (d2 >  180.0) d2 -= 360.0;
		while (d2 <= -180.0) d2 += 360.0;
		double delta = d2 - d1;
		while (delta >  180.0) delta -= 360.0;
		while (delta <= -180.0) delta += 360.0;
		double e2 = d1 + delta;
		if ((d1 < 0.0) == (e2 < 0.0)) continue;

		if (std::fabs(Nz) < 1e-15) continue;
		double lat_cross = atan(-(ring.Nx[i] * cos_plon_rad + ring.Ny[i] * sin_plon_rad) / Nz);
		if (lat_cross > plat_rad) crossings++;
	}
	bool odd = (crossings & 1) == 1;
	return odd ^ ring.flip_north;
}

