// author Robert Hijmans
// October 2009
// version 0.1
// license GPL

/*

#ifndef M_PI
#define M_PI (3.14159265358979323846)
#endif

#ifndef toRAD
#define toRAD (M_PI / 180)
#endif

#ifndef geod_a
#define geod_a (6378137)
#endif

#ifndef geod_f
#define geod_f (0.0033528106647474805)
#endif


#include <vector>
#include <math.h>
#include "ggeodesic.h"


void DegtoRad(double &deg) {
	deg *= 0.0174532925199433;
}

void normalizeLonDeg(double &x) {
	x = fmod((x + 180), 360 - 180);
}

void normalizeLonRad(double &x) {
	x = fmod((x + M_PI), (2 * M_PI) - M_PI);
}





bool antipodal(double x1, double y1, double x2, const double y2, const double &tol=1e-9) {
	normalizeLonDeg(x1);
	normalizeLonDeg(x2);
	double diflon = fabs(x1 - x2);
	double diflat = fabs(y1 - y2);
	DegtoRad(y1);
	return ((diflat < tol) && (cos(y1) * fabs(fmod(diflon, 360) - 180) < tol));
}


void antipode(const double &x, const double &y, double &xa, double &ya) {
	xa = x + 180;
	normalizeLonDeg(xa);
	ya = -y;
}

double bearing(double lon1, double lat1, double lon2, double lat2 ) {
	double s12, azi1, azi2;
	struct geod_geodesic g;
	geod_init(&g, geod_a, geod_f);
	geod_inverse(&g, lat1, lon1, lat2, lon2, &s12, &azi1, &azi2);
  	return azi1;
}


void destPoint(double lon1, double lat1, double lon2, double lat2, double d, double &dlon, double &dlat) {
	double s12, azi1, azi2;
	struct geod_geodesic g;
	geod_init(&g, geod_a, geod_f);
	geod_inverse(&g, lat1, lon1, lat2, lon2, &s12, &azi1, &azi2);
	geod_direct(&g, lat1, lon1, azi1, d, &dlat, &dlon, &azi2);
}


void intermediate_points(double lon1, double lat1, double lon2, double lat2, size_t n, std::vector<double> &lon, std::vector<double> &lat) {

	lon.resize(0);
	lat.resize(0);

	if ((lon1 == lon2) && (lat1 == lat2)) {
		return;
	}

	if (antipodal(lon1, lat1, lon2, lat2)) {
		lon.resize(1, NAN);
		lat.resize(1, NAN);
		return;
	}

	lon.reserve(n);
	lat.reserve(n);
	lon.push_back(lon1);
	lat.push_back(lon2);
	
	double azi1, azi2, s12, dlon, dlat;
	struct geod_geodesic g;
	geod_init(&g, geod_a, geod_f);
	geod_inverse(&g, lat1, lon1, lat2, lon2, &s12, &azi1, &azi2);

	double step = s12 / (n+1);
	for (size_t i=0; i<n; i++) {
		geod_direct(&g, lat1, lon1, azi1, step*(i+1), &dlat, &dlon, &azi2);
		lon.push_back(dlon);
		lat.push_back(dlat);
	}
	lon.push_back(lon2);
	lat.push_back(lat2);

	//if (breakAtDateLine) {
	//	res[[i]] <- .breakAtDateLine(x)
	//}

}


double dist2segment (double lon1, double lat1, double lon2, double lat2, double lon3, double lat3) {

	bool sign = false;

	double s12, azi1, azi2;
	struct geod_geodesic g;
	//geod_init(&g, r, 0);
	geod_init(&g, geod_a, geod_f);

	geod_inverse(&g, lat1, lon1, lat2, lon2, &s12, &azi1, &azi2);
	double tc = azi1 * toRAD;
	geod_inverse(&g, lat1, lon1, lat3, lon3, &s12, &azi1, &azi2);
	double tcp = azi1 * toRAD;
	//geod_init(&g, 1, 0);
	geod_init(&g, 1, geod_f);
	geod_inverse(&g, lat1, lon1, lat3, lon3, &s12, &azi1, &azi2);
	double xtr = (asin(sin(tcp-tc) * sin(s12)) * geod_a);
	xtr = sign ? xtr : fabs(xtr);
	return xtr;
}

double alongTrackDistance(double lon1, double lat1, double lon2, double lat2, double lon3, double lat3) {

	double s12, azi1, azi2;
	struct geod_geodesic g;
	//geod_init(&g, r, 0);
	geod_init(&g, geod_a, geod_f);

	geod_inverse(&g, lat1, lon1, lat2, lon2, &s12, &azi1, &azi2);
	double tc = azi1 * toRAD;
	geod_inverse(&g, lat1, lon1, lat3, lon3, &s12, &azi1, &azi2);
	double tcp = azi1 * toRAD;
	//geod_init(&g, 1, 0);
	geod_init(&g, 1, geod_f);
	geod_inverse(&g, lat1, lon1, lat3, lon3, &s12, &azi1, &azi2);

	double xtr = asin(sin(tcp-tc) * sin(s12));

// +1/-1 for ahead/behind [lat1,lon1]
// is std::abs enough?? or do we need 0 as well?
//	bearing = sign(cos(tc - tcp)) ;
//	double dist = bearing * acos(cos(s12) / cos(xtr)) * geod_a;
	double dist = acos(cos(s12) / cos(xtr)) * geod_a;
	return fabs(dist);
}
*/

/*

dist2Line <- function(p, line, distfun=distGeo) {

	line <- .pointsToMatrix(line) 
	line1 <- line[-nrow(line), ,drop=FALSE]
	line2 <- line[-1, ,drop=FALSE]
	seglength  <- distfun(line1, line2)
	
	res <- matrix(nrow=nrow(p), ncol=3)
	colnames(res) <- c("distance","lon","lat")
	
	for (i in 1:nrow(p)) {
		xy <- p[i,]
# the shortest distance of a point to a great circle
		crossdist <- abs(dist2gc(line1, line2, xy))
		
# the alongTrackDistance is the length of the path along the great circle to the point of intersection
# there are two, depending on which node you start
# we want to use the min, but the max needs to be < segment length
		trackdist1 <- alongTrackDistance(line1, line2, xy)
		trackdist2 <- alongTrackDistance(line2, line1, xy)
		mintrackdist <- pmin(trackdist1, trackdist2)
		maxtrackdist <- pmax(trackdist1, trackdist2)
		crossdist[maxtrackdist >= seglength] <- NA 
		
# if the crossdist is NA, we use the distance to the nodes
		nodedist <- distfun(xy, line)
		
		warnopt = getOption('warn')
	 	options('warn'=-1) 		
		distmin1 <- min(nodedist, na.rm=TRUE)
		distmin2 <- min(crossdist, na.rm=TRUE)
		options('warn'= warnopt) 
		
		if (distmin1 <= distmin2) {
			j <- which.min(nodedist)
			res[i,] <- c(distmin1, line[j,])
		} else {
			j <- which.min(crossdist)
			# if else to determine from which node to start
			if (trackdist1[j] < trackdist2[j]) {
				bear <- bearing(line1[j,], line2[j,])
				pt <- destPoint(line1[j,], bear, mintrackdist[j])
				res[i,] <- c(crossdist[j], pt)
			} else {
				bear <- bearing(line2[j,], line1[j,])
				pt <- destPoint(line2[j,], bear, mintrackdist[j])
				res[i,] <- c(crossdist[j], pt)	
			}
		}
	}
	return(res)
}

 */


/*
double distance_haversine(double lon1, double lat1, double lon2, double lat2) {
// Haversine formula to calculate distance between two points specified by 
// from: Haversine formula - R.W. Sinnott, "Virtues of the Haversine",
//  Sky and Telescope, vol 68, no 2, 1984
//  http:#//www.census.gov/cgi-bin/geo/gisfaq?Q5.1

	double r=6378137;
	lon1 = lon1 * toRAD;
	lat1 = lat1 * toRAD;
	lon2 = lon2 * toRAD;
	lat2 = lat2 * toRAD;
		
	double dLat = lat2-lat1;
	double dLon = lon2-lon1;
	double a = pow(sin(dLat/2), 2) + cos(lat1) * cos(lat2) * pow(sin(dLon/2), 2);
	// to avoid values of 'a' that are a sliver above 1, which may occur at antipodes
	// https://stackoverflow.com/q/45889616/635245 
	a = a > 1 ? 1 : a;
	return 2 * atan2(sqrt(a), sqrt(1-a)) * r;
}

*/
