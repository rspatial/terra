#include "spatVector.h"
#include <numeric>
#include "vecmath.h"
#include <array>
#include <map>
#include <unordered_map>

//#include "gdal_alg.h"
//#include "ogrsf_frmts.h"

//#ifdef useGDAL
	#include "crs.h"
//#endif



// Build a hexagonal tessellation covering `e`.
// "size" is the across-flats distance:
//   pointy-top: hex width  (horizontal, between vertical edges)
//   flat-top:   hex height (vertical, between horizontal edges)
//
// If anchor_x or anchor_y is NaN the grid is anchored at the lower-left of
// the extent: the first hex center is at (xmin + size/2, ymin + s) for
// pointy-top (or (xmin + s, ymin + size/2) for flat-top). Otherwise the grid
// is anchored globally so that a hex is centered at (anchor_x, anchor_y) and
// the same hex layout is produced regardless of the extent (different extents
// produce mutually consistent, aligned cells).

SpatVector SpatVector::hexagons(SpatExtent e, double size, std::string crs, bool flat_top, double anchor_x, double anchor_y) {

	SpatVector out;
	if (size <= 0) {
		out.setError("size must be a positive number");
		return out;
	}
	if (!e.valid_notempty()) {
		out.setError("invalid extent");
		return out;
	}

	const double sqrt3 = std::sqrt(3.0);
	const double s = size / sqrt3;     // hex side length

	std::vector<double> vx_off(6), vy_off(6);
	double dx, dy;                     // center-to-center spacing in x and y
	double half_w, half_h;             // half hex bounding-box width / height

	if (flat_top) {
		// vertices at angles 0, 60, 120, 180, 240, 300
		for (int k = 0; k < 6; k++) {
			double a = M_PI / 180.0 * (60.0 * k);
			vx_off[k] = s * std::cos(a);
			vy_off[k] = s * std::sin(a);
		}
		dx = 1.5 * s;
		dy = size;             // = sqrt(3) * s
		half_w = s;
		half_h = size / 2.0;
	} else {
		// pointy-top: vertices at angles 30, 90, 150, 210, 270, 330
		for (int k = 0; k < 6; k++) {
			double a = M_PI / 180.0 * (60.0 * k + 30.0);
			vx_off[k] = s * std::cos(a);
			vy_off[k] = s * std::sin(a);
		}
		dx = size;             // = sqrt(3) * s
		dy = 1.5 * s;
		half_w = size / 2.0;
		half_h = s;
	}

	bool global_anchor = !std::isnan(anchor_x) && !std::isnan(anchor_y);

	double ax, ay;             // anchor (a hex is centered exactly here)
	long i_min, i_max, j_min, j_max;

	if (global_anchor) {
		ax = anchor_x;
		ay = anchor_y;
	} else {
		// anchor = lower-left first hex center (preserves legacy behaviour)
		if (flat_top) {
			ax = e.xmin + s;
			ay = e.ymin + size / 2.0;
		} else {
			ax = e.xmin + size / 2.0;
			ay = e.ymin + s;
		}
	}

	// determine the global (i, j) integer grid range covering `e`
	if (flat_top) {
		// cx = ax + i*dx
		// cy = ay + j*dy + (i mod 2) * dy/2
		i_min = (long) std::floor((e.xmin - ax - half_w) / dx);
		i_max = (long) std::ceil ((e.xmax - ax + half_w) / dx);
		j_min = (long) std::floor((e.ymin - ay - half_h - dy / 2.0) / dy);
		j_max = (long) std::ceil ((e.ymax - ay + half_h + dy / 2.0) / dy);
	} else {
		// cx = ax + i*dx + (j mod 2) * dx/2
		// cy = ay + j*dy
		i_min = (long) std::floor((e.xmin - ax - half_w - dx / 2.0) / dx);
		i_max = (long) std::ceil ((e.xmax - ax + half_w + dx / 2.0) / dx);
		j_min = (long) std::floor((e.ymin - ay - half_h) / dy);
		j_max = (long) std::ceil ((e.ymax - ay + half_h) / dy);
	}

	if (i_max < i_min || j_max < j_min) {
		out.setSRS({crs});
		return out;
	}

	out.reserve((size_t)(i_max - i_min + 1) * (size_t)(j_max - j_min + 1));

	for (long j = j_min; j <= j_max; j++) {
		long jpar = ((j % 2) + 2) % 2;             // safe mod for negatives
		double row_off_x = (!flat_top && jpar) ? dx / 2.0 : 0.0;
		double base_cy = ay + j * dy;

		for (long i = i_min; i <= i_max; i++) {
			long ipar = ((i % 2) + 2) % 2;
			double col_off_y = (flat_top && ipar) ? dy / 2.0 : 0.0;
			double cx = ax + i * dx + row_off_x;
			double cy = base_cy + col_off_y;

			// drop if hex bounding box does not intersect the extent
			if (cx + half_w < e.xmin) continue;
			if (cx - half_w > e.xmax) continue;
			if (cy + half_h < e.ymin) continue;
			if (cy - half_h > e.ymax) continue;

			std::vector<double> hx(7), hy(7);
			for (int k = 0; k < 6; k++) {
				hx[k] = cx + vx_off[k];
				hy[k] = cy + vy_off[k];
			}
			hx[6] = hx[0];
			hy[6] = hy[0];

			SpatPart p(hx, hy);
			SpatGeom g(p, polygons);
			out.addGeom(g);
		}
	}

	out.setSRS({crs});
	return out;
}


// Equal-area hexagons tessellation in a single, fixed
// global cylindrical equal-area projection (Lambert CEA on a sphere)
// "size" is the across-flats distance of a hexagon, in meters in CEA.
SpatVector SpatVector::hexagons_lonlat(SpatExtent e, double size, bool flat_top) {

	SpatVector out;
	if (size <= 0) {
		out.setError("size must be a positive number");
		return out;
	}
	if (!e.valid_notempty()) {
		out.setError("invalid extent");
		return out;
	}

	const std::string lonlat_crs = "+proj=longlat +R=6378137 +no_defs";
	const std::string cea_crs    = "+proj=cea +lon_0=0 +lat_ts=0 +R=6378137 +units=m +no_defs";

	const double R    = 6378137.0;
	const double sqrt3 = std::sqrt(3.0);
	const double circ = 2.0 * M_PI * R;     // CEA full x-extent (sphere circumference)

	// Snap "size" so that the horizontal cell spacing dx divides circ exactly.
	// Without this snap, dx*Nx != circ and cells along the antimeridian would
	// either gap or duplicate (overlap). The snap is typically <0.5%.
	double dx_orig = flat_top ? (1.5 * size / sqrt3) : size;
	double Nx_d = std::round(circ / dx_orig);
	if (Nx_d < 1.0) Nx_d = 1.0;
	double dx_snap = circ / Nx_d;
	size = flat_top ? (dx_snap * sqrt3 / 1.5) : dx_snap;

	// project the lon/lat extent rectangle to CEA. Sample many points along
	// the boundary so the projected bbox captures the curvature of parallels.
	const int nstep = 200;
	std::vector<double> bx, by;
	bx.reserve(4 * (nstep + 1));
	by.reserve(4 * (nstep + 1));
	for (int i = 0; i <= nstep; i++) {
		double t = (double)i / (double)nstep;
		bx.push_back(e.xmin + t * (e.xmax - e.xmin));
		by.push_back(e.ymin);
	}
	for (int i = 0; i <= nstep; i++) {
		double t = (double)i / (double)nstep;
		bx.push_back(e.xmax);
		by.push_back(e.ymin + t * (e.ymax - e.ymin));
	}
	for (int i = 0; i <= nstep; i++) {
		double t = (double)i / (double)nstep;
		bx.push_back(e.xmax - t * (e.xmax - e.xmin));
		by.push_back(e.ymax);
	}
	for (int i = 0; i <= nstep; i++) {
		double t = (double)i / (double)nstep;
		bx.push_back(e.xmin);
		by.push_back(e.ymax - t * (e.ymax - e.ymin));
	}

	SpatMessages tmsg = transform_coordinates(bx, by, lonlat_crs, cea_crs);
	if (tmsg.has_error) {
		out.setError(tmsg.getError());
		return out;
	}

	double cea_xmin = bx[0], cea_xmax = bx[0];
	double cea_ymin = by[0], cea_ymax = by[0];
	for (size_t i = 1; i < bx.size(); i++) {
		if (!std::isnan(bx[i])) {
			if (bx[i] < cea_xmin) cea_xmin = bx[i];
			if (bx[i] > cea_xmax) cea_xmax = bx[i];
		}
		if (!std::isnan(by[i])) {
			if (by[i] < cea_ymin) cea_ymin = by[i];
			if (by[i] > cea_ymax) cea_ymax = by[i];
		}
	}

	// Clamp the y range so that no generated hex has a vertex outside the
	// valid CEA range [-R, R]. Cells whose center lies more than half the
	// hex height from the pole would have a top/bottom vertex past CEA y=±R
	// (i.e., past the pole), and would fail to project back to lon/lat.
	const double half_h_y = flat_top ? (size / 2.0) : (size / sqrt3);
	const double clamp_ymax = R - 2.0 * half_h_y;
	const double clamp_ymin = -R + 2.0 * half_h_y;
	if (clamp_ymax <= clamp_ymin) {
		out.setError("size is too large for a global tessellation");
		return out;
	}
	if (cea_ymax > clamp_ymax) cea_ymax = clamp_ymax;
	if (cea_ymin < clamp_ymin) cea_ymin = clamp_ymin;
	if (cea_ymax <= cea_ymin) {
		out.setError("requested extent does not fit on the sphere given size");
		return out;
	}

	// If the request spans (or exceeds) the full sphere in longitude, trim
	// the x range so the bbox-cull selects exactly Nx unique cells per row.
	// Without trimming the bbox rectangle would include duplicate cells at
	// both ends (cx = -R*pi and cx = +R*pi map to the same lon = 180).
	// We pad the left side by `eps` (drops the cell at i = -Nx/2 in odd
	// rows, whose right edge sits exactly on the left bound) and trim the
	// right side by half a hex width (drops the cell at i = +Nx/2 in even
	// rows). Result: Nx unique cells per row in both row parities.
	if (cea_xmax - cea_xmin >= circ - 1.0) {
		double eps = dx_snap * 1e-6;
		double xmin0 = cea_xmin;
		cea_xmin = xmin0 + eps;
		cea_xmax = xmin0 + circ - dx_snap * 0.5 - eps;
	}

	SpatExtent cea_ext(cea_xmin, cea_xmax, cea_ymin, cea_ymax);

	// build the planar hex tessellation in CEA, anchored at the global origin
	SpatVector hex = hexagons(cea_ext, size, cea_crs, flat_top, 0.0, 0.0);
	if (hex.hasError()) {
		out.setError(hex.getError());
		return out;
	}

	// Split, in CEA, any hex that straddles the antimeridian (x = +/- R*pi).
	// We can't rely on splitting AFTER back-projection because PROJ wraps
	// vertices near x = +/- R*pi inconsistently to lon = +/- 180, producing
	// self-intersecting polygons. Instead we Sutherland-Hodgman-clip the
	// hex at x = +/- R*pi in CEA and shift the wrapping half by +/- circ
	// so both halves have x in (-R*pi, R*pi). 
	const double R_pi = M_PI * R;
	for (size_t i = 0; i < hex.geoms.size(); i++) {
		SpatGeom &g = hex.geoms[i];
		if (g.parts.empty()) continue;
		const std::vector<double> &x = g.parts[0].x;
		const std::vector<double> &y = g.parts[0].y;
		if (x.empty()) continue;
		double xmn = x[0], xmx = x[0];
		for (size_t k = 1; k < x.size(); k++) {
			if (x[k] < xmn) xmn = x[k];
			if (x[k] > xmx) xmx = x[k];
		}
		double clip_x;
		double shift;
		if (xmn < -R_pi) {
			clip_x = -R_pi;
			shift  = circ;          // shift left half (x < -R*pi) by +circ
		} else if (xmx > R_pi) {
			clip_x = R_pi;
			shift  = -circ;         // shift right half (x > R*pi) by -circ
		} else {
			continue;
		}

		// Sutherland-Hodgman: clip a closed polygon by a vertical line into
		// "inside" (x within range, sign = +1 means keep x <= clip_x with
		// shift; sign = -1 means keep x >= clip_x un-shifted).
		auto clip = [&](bool keep_le, std::vector<double> &xo, std::vector<double> &yo) {
			xo.clear(); yo.clear();
			size_t n = x.size();
			// treat as closed ring: skip a duplicate trailing vertex if present
			size_t m = (n > 1 && x.front() == x.back() && y.front() == y.back()) ? n - 1 : n;
			for (size_t k = 0; k < m; k++) {
				double x1 = x[k],         y1 = y[k];
				double x2 = x[(k + 1) % m], y2 = y[(k + 1) % m];
				bool in1 = keep_le ? (x1 <= clip_x) : (x1 >= clip_x);
				bool in2 = keep_le ? (x2 <= clip_x) : (x2 >= clip_x);
				if (in1 && in2) {
					xo.push_back(x2); yo.push_back(y2);
				} else if (in1 && !in2) {
					double t = (clip_x - x1) / (x2 - x1);
					xo.push_back(clip_x); yo.push_back(y1 + t * (y2 - y1));
				} else if (!in1 && in2) {
					double t = (clip_x - x1) / (x2 - x1);
					xo.push_back(clip_x); yo.push_back(y1 + t * (y2 - y1));
					xo.push_back(x2); yo.push_back(y2);
				}
			}
			if (xo.size() >= 3) {
				xo.push_back(xo[0]); yo.push_back(yo[0]);
			}
		};

		// Wrap-side half: x outside [-R*pi, R*pi]; will be shifted by `shift`.
		// Keep-side half: x inside; left as-is.
		bool wrap_is_le = (xmn < -R_pi);  // wrap half is the x < clip_x side
		std::vector<double> wx, wy, kx, ky;
		clip( wrap_is_le, wx, wy);   // wrap half
		clip(!wrap_is_le, kx, ky);   // keep half

		if (wx.size() < 4 || kx.size() < 4) continue;  // degenerate, leave as-is

		for (size_t k = 0; k < wx.size(); k++) wx[k] += shift;

		SpatGeom ng(g.gtype);
		ng.addPart(SpatPart(kx, ky));
		ng.addPart(SpatPart(wx, wy));
		ng.computeExtent();
		hex.geoms[i] = ng;
	}

	// densify the (planar) edges so the back-projected curves are smooth
	double densify_interval = size / 30.0;
	SpatVector hexd = hex.densify(densify_interval, true, true);
	if (hexd.hasError()) {
		out.setError(hexd.getError());
		return out;
	}

	// project to longitude/latitude
	std::vector<double> empty_aoi;
	out = hexd.project(lonlat_crs, false, "", empty_aoi, 0.0, true);
	if (out.hasError()) return out;

	for (size_t i = 0; i < out.geoms.size(); i++) {
		out.geoms[i].computeExtent();
	}

	return out;
}


// Goldberg-polyhedron tessellation of the sphere: dual of a frequency-n 
// subdivided icosahedron. The output has 10 n^2 + 2 cells 
//
// e         : longitude/latitude extent used to filter cells (by center).
// n         : subdivision frequency (positive integer).
// full_globe: if true, ignore e and return all cells.

SpatVector SpatVector::polyhedron(SpatExtent e, int n, bool full_globe) {

	SpatVector out;
	if (n < 1) {
		out.setError("n must be a positive integer");
		return out;
	}
	if (!full_globe && !e.valid_notempty()) {
		out.setError("invalid extent");
		return out;
	}

	const std::string lonlat_crs = "+proj=longlat +R=6378137 +no_defs";
	const double DEG = 180.0 / M_PI;

	// ---- icosahedron vertices and faces ------------------------------------
	double phi = (1.0 + std::sqrt(5.0)) / 2.0;
	double iv[12][3] = {
		{ 0,  1,  phi}, { 0, -1,  phi}, { 0,  1, -phi}, { 0, -1, -phi},
		{ 1,  phi,  0}, {-1,  phi,  0}, { 1, -phi,  0}, {-1, -phi,  0},
		{ phi,  0,  1}, {-phi,  0,  1}, { phi,  0, -1}, {-phi,  0, -1}
	};
	for (int i = 0; i < 12; i++) {
		double L = std::sqrt(iv[i][0]*iv[i][0] + iv[i][1]*iv[i][1] + iv[i][2]*iv[i][2]);
		iv[i][0] /= L; iv[i][1] /= L; iv[i][2] /= L;
	}
	// rotate so vertex 0 = normalize(0, 1, phi) -> (0, 0, 1) and vertex 3 ->
	// (0, 0, -1). Then the two polar pentagons are centered on the poles.
	double cs = phi / std::sqrt(1.0 + phi*phi);
	double sn = 1.0  / std::sqrt(1.0 + phi*phi);
	for (int i = 0; i < 12; i++) {
		double y = iv[i][1], z = iv[i][2];
		iv[i][1] =  cs * y + sn * z;
		iv[i][2] = -sn * y + cs * z;
	}

	// edges = pairs at the (unique) minimum pairwise distance
	double edge_d2 = 1e18;
	for (int i = 0; i < 12; i++) {
		for (int j = i+1; j < 12; j++) {
			double dx = iv[i][0]-iv[j][0], dy = iv[i][1]-iv[j][1], dz = iv[i][2]-iv[j][2];
			double d2 = dx*dx + dy*dy + dz*dz;
			if (d2 > 1e-9 && d2 < edge_d2) edge_d2 = d2;
		}
	}
	bool adj[12][12];
	for (int i = 0; i < 12; i++) for (int j = 0; j < 12; j++) adj[i][j] = false;
	for (int i = 0; i < 12; i++) {
		for (int j = i+1; j < 12; j++) {
			double dx = iv[i][0]-iv[j][0], dy = iv[i][1]-iv[j][1], dz = iv[i][2]-iv[j][2];
			double d2 = dx*dx + dy*dy + dz*dz;
			if (std::abs(d2 - edge_d2) < 1e-6) {
				adj[i][j] = true; adj[j][i] = true;
			}
		}
	}
	std::vector<std::array<int,3>> faces;
	faces.reserve(20);
	for (int i = 0; i < 10; i++) {
		for (int j = i+1; j < 11; j++) {
			if (!adj[i][j]) continue;
			for (int k = j+1; k < 12; k++) {
				if (adj[j][k] && adj[i][k]) {
					faces.push_back({i, j, k});
				}
			}
		}
	}
	// orient outward (CCW seen from outside)
	for (auto &f : faces) {
		double ax=iv[f[0]][0], ay=iv[f[0]][1], az=iv[f[0]][2];
		double bx=iv[f[1]][0], by=iv[f[1]][1], bz=iv[f[1]][2];
		double cx=iv[f[2]][0], cy=iv[f[2]][1], cz=iv[f[2]][2];
		double nx = (by-ay)*(cz-az) - (bz-az)*(cy-ay);
		double ny = (bz-az)*(cx-ax) - (bx-ax)*(cz-az);
		double nz = (bx-ax)*(cy-ay) - (by-ay)*(cx-ax);
		if (nx*ax + ny*ay + nz*az < 0) std::swap(f[1], f[2]);
	}

	// ---- subdivide and dedup vertices --------------------------------------
	const long snap = 10000000L;     // 1e7 grid for shared-edge dedup
	std::map<std::array<long,3>, int> key_map;
	std::vector<std::array<double,3>> xyz;
	xyz.reserve((size_t)10 * n * n + 12);

	int nF = (int) faces.size();
	int row_stride = n + 1;
	std::vector<std::vector<int>> face_grid(nF, std::vector<int>(row_stride * row_stride, -1));

	for (int f = 0; f < nF; f++) {
		const auto &fv = faces[f];
		double Ax=iv[fv[0]][0], Ay=iv[fv[0]][1], Az=iv[fv[0]][2];
		double Bx=iv[fv[1]][0], By=iv[fv[1]][1], Bz=iv[fv[1]][2];
		double Cx=iv[fv[2]][0], Cy=iv[fv[2]][1], Cz=iv[fv[2]][2];
		for (int i = 0; i <= n; i++) {
			for (int j = 0; j <= n - i; j++) {
				double u = (double)i / n, v = (double)j / n, w = 1.0 - u - v;
				double px = u*Ax + v*Bx + w*Cx;
				double py = u*Ay + v*By + w*Cy;
				double pz = u*Az + v*Bz + w*Cz;
				double L = std::sqrt(px*px + py*py + pz*pz);
				px /= L; py /= L; pz /= L;
				std::array<long,3> key = {
					(long) std::llround(px * snap),
					(long) std::llround(py * snap),
					(long) std::llround(pz * snap)
				};
				auto it = key_map.find(key);
				int id;
				if (it == key_map.end()) {
					id = (int) xyz.size();
					xyz.push_back({px, py, pz});
					key_map[key] = id;
				} else {
					id = it->second;
				}
				face_grid[f][i * row_stride + j] = id;
			}
		}
	}
	int nv = (int) xyz.size();

	// mark icosahedron-corner vertices (these become pentagonal cells)
	std::vector<bool> is_corner(nv, false);
	for (int f = 0; f < nF; f++) {
		is_corner[face_grid[f][0]] = true;                    // (i=0, j=0): A
		is_corner[face_grid[f][n * row_stride]] = true;       // (i=n, j=0): B
		is_corner[face_grid[f][n]] = true;                    // (i=0, j=n): C
	}

	// ---- enumerate small triangles, sphere-projected centroids -------------
	int nt_est = nF * n * n;
	std::vector<std::array<int,3>> tris;
	tris.reserve(nt_est);
	std::vector<std::array<double,3>> centers;
	centers.reserve(nt_est);
	std::vector<std::vector<int>> vert_tris(nv);

	auto add_tri = [&](int v1, int v2, int v3) {
		int t = (int) tris.size();
		tris.push_back({v1, v2, v3});
		double cx = (xyz[v1][0] + xyz[v2][0] + xyz[v3][0]) / 3.0;
		double cy = (xyz[v1][1] + xyz[v2][1] + xyz[v3][1]) / 3.0;
		double cz = (xyz[v1][2] + xyz[v2][2] + xyz[v3][2]) / 3.0;
		double L = std::sqrt(cx*cx + cy*cy + cz*cz);
		centers.push_back({cx/L, cy/L, cz/L});
		vert_tris[v1].push_back(t);
		vert_tris[v2].push_back(t);
		vert_tris[v3].push_back(t);
	};

	for (int f = 0; f < nF; f++) {
		const auto &g = face_grid[f];
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < n - i; j++) {
				add_tri(g[i*row_stride + j], g[(i+1)*row_stride + j], g[i*row_stride + (j+1)]);
			}
		}
		for (int i = 0; i < n - 1; i++) {
			for (int j = 0; j < n - 1 - i; j++) {
				add_tri(g[(i+1)*row_stride + j], g[(i+1)*row_stride + (j+1)], g[i*row_stride + (j+1)]);
			}
		}
	}

	// ---- build dual cells, densify edges, lon/lat conversion ---------------
	auto slerp = [](double ax, double ay, double az,
	                double bx, double by, double bz, double t,
	                double &ox, double &oy, double &oz) {
		double cs = ax*bx + ay*by + az*bz;
		if (cs >  1.0) cs =  1.0;
		if (cs < -1.0) cs = -1.0;
		double o = std::acos(cs);
		if (o < 1e-12) { ox = ax; oy = ay; oz = az; return; }
		double so = std::sin(o);
		double s1 = std::sin((1.0 - t) * o) / so;
		double s2 = std::sin(t * o) / so;
		ox = s1*ax + s2*bx;
		oy = s1*ay + s2*by;
		oz = s1*az + s2*bz;
		double L = std::sqrt(ox*ox + oy*oy + oz*oz);
		ox /= L; oy /= L; oz /= L;
	};

	int n_per_edge = std::max(8, (int) std::ceil(80.0 / n));

	std::vector<std::string> type_col;
	out.reserve((size_t) nv);

	for (int v = 0; v < nv; v++) {
		double cxv = xyz[v][0], cyv = xyz[v][1], czv = xyz[v][2];
		double clon = std::atan2(cyv, cxv) * DEG;
		double tmp_z = czv;
		if (tmp_z >  1.0) tmp_z =  1.0;
		if (tmp_z < -1.0) tmp_z = -1.0;
		double clat = std::asin(tmp_z) * DEG;

		if (!full_globe) {
			if (clat < e.ymin || clat > e.ymax) continue;
			if (e.xmin > -180 || e.xmax < 180) {
				bool ok = (clon >= e.xmin && clon <= e.xmax) ||
				          (clon + 360 >= e.xmin && clon + 360 <= e.xmax) ||
				          (clon - 360 >= e.xmin && clon - 360 <= e.xmax);
				if (!ok) continue;
			}
		}

		// tangent-plane basis (right-handed seen from outside the sphere)
		double Vp[3] = {cxv, cyv, czv};
		double ax[3];
		if (std::abs(Vp[2]) < 0.9) { 
			ax[0]=0; ax[1]=0; ax[2]=1; 
		} else { 
			ax[0]=1; ax[1]=0; ax[2]=0; 
		}
		double dot = ax[0]*Vp[0] + ax[1]*Vp[1] + ax[2]*Vp[2];
		double e1[3] = {ax[0] - dot*Vp[0], ax[1] - dot*Vp[1], ax[2] - dot*Vp[2]};
		double L1 = std::sqrt(e1[0]*e1[0] + e1[1]*e1[1] + e1[2]*e1[2]);
		e1[0]/=L1; e1[1]/=L1; e1[2]/=L1;
		double e2v[3] = {Vp[1]*e1[2] - Vp[2]*e1[1],
		                 Vp[2]*e1[0] - Vp[0]*e1[2],
		                 Vp[0]*e1[1] - Vp[1]*e1[0]};

		const auto &tids = vert_tris[v];
		size_t k = tids.size();
		std::vector<std::pair<double,int>> ang_id(k);
		for (size_t i = 0; i < k; i++) {
			int t = tids[i];
			double rx = centers[t][0] - Vp[0];
			double ry = centers[t][1] - Vp[1];
			double rz = centers[t][2] - Vp[2];
			double d1 = rx*e1[0]  + ry*e1[1]  + rz*e1[2];
			double d2 = rx*e2v[0] + ry*e2v[1] + rz*e2v[2];
			ang_id[i] = std::make_pair(std::atan2(d2, d1), t);
		}
		std::sort(ang_id.begin(), ang_id.end());

		// densified open ring (with closing duplicate appended)
		std::vector<double> lon, lat;
		lon.reserve(k * n_per_edge + 1);
		lat.reserve(k * n_per_edge + 1);
		for (size_t i = 0; i < k; i++) {
			int ta = ang_id[i].second;
			int tb = ang_id[(i + 1) % k].second;
			double pax = centers[ta][0], pay = centers[ta][1], paz = centers[ta][2];
			double pbx = centers[tb][0], pby = centers[tb][1], pbz = centers[tb][2];
			for (int s = 0; s < n_per_edge; s++) {
				double t = (double) s / n_per_edge;
				double px, py, pz;
				slerp(pax, pay, paz, pbx, pby, pbz, t, px, py, pz);
				double pl = std::atan2(py, px) * DEG;
				double tz = pz; if (tz > 1.0) tz = 1.0; if (tz < -1.0) tz = -1.0;
				double pa = std::asin(tz) * DEG;
				lon.push_back(pl);
				lat.push_back(pa);
			}
		}
		lon.push_back(lon[0]);
		lat.push_back(lat[0]);

		// antimeridian split
		size_t m = lon.size() - 1;
		std::vector<std::vector<double>> piece_x, piece_y;
		std::vector<double> cur_x, cur_y;
		int crossings = 0;
		for (size_t i = 0; i < m; i++) {
			cur_x.push_back(lon[i]);
			cur_y.push_back(lat[i]);
			double lon1 = lon[i], lon2 = lon[i+1];
			double lat1 = lat[i], lat2 = lat[i+1];
			if (std::abs(lon2 - lon1) > 180) {
				crossings++;
				double d, t, lat_c;
				if (lon1 > 0) {
					d = (lon2 + 360) - lon1;
					t = (180 - lon1) / d;
					lat_c = lat1 + t * (lat2 - lat1);
					cur_x.push_back( 180); cur_y.push_back(lat_c);
					piece_x.push_back(cur_x); piece_y.push_back(cur_y);
					cur_x.clear(); cur_y.clear();
					cur_x.push_back(-180); cur_y.push_back(lat_c);
				} else {
					d = (lon1 + 360) - lon2;
					t = (lon1 + 180) / d;
					lat_c = lat1 + t * (lat2 - lat1);
					cur_x.push_back(-180); cur_y.push_back(lat_c);
					piece_x.push_back(cur_x); piece_y.push_back(cur_y);
					cur_x.clear(); cur_y.clear();
					cur_x.push_back( 180); cur_y.push_back(lat_c);
				}
			}
		}
		piece_x.push_back(cur_x); piece_y.push_back(cur_y);

		SpatGeom geom(polygons);

		if (crossings == 0) {
			// single ring, already closed
			geom.addPart(SpatPart(lon, lat));
		} else if (crossings % 2 == 0) {
			// merge first and last piece (they belong to the same side of
			// the antimeridian; the walk just started in the middle of one)
			if (piece_x.size() >= 2) {
				size_t L = piece_x.size() - 1;
				piece_x[L].insert(piece_x[L].end(), piece_x[0].begin(), piece_x[0].end());
				piece_y[L].insert(piece_y[L].end(), piece_y[0].begin(), piece_y[0].end());
				std::vector<std::vector<double>> nx, ny;
				nx.push_back(piece_x[L]); ny.push_back(piece_y[L]);
				for (size_t i = 1; i < L; i++) {
					nx.push_back(piece_x[i]); ny.push_back(piece_y[i]);
				}
				piece_x = nx; piece_y = ny;
			}
			for (size_t i = 0; i < piece_x.size(); i++) {
				auto &px = piece_x[i]; auto &py = piece_y[i];
				if (px.size() < 3) continue;
				if (px.front() != px.back() || py.front() != py.back()) {
					px.push_back(px.front());
					py.push_back(py.front());
				}
				geom.addPart(SpatPart(px, py));
			}
		} else {
			// odd crossings <=> the cell encloses a pole. Rotate so the
			// merged ring runs from one antimeridian side to the other,
			// then close along the antimeridian and across the pole.
			std::vector<double> mx, my;
			for (size_t i = 1; i < piece_x.size(); i++) {
				mx.insert(mx.end(), piece_x[i].begin(), piece_x[i].end());
				my.insert(my.end(), piece_y[i].begin(), piece_y[i].end());
			}
			mx.insert(mx.end(), piece_x[0].begin(), piece_x[0].end());
			my.insert(my.end(), piece_y[0].begin(), piece_y[0].end());
			double sumlat = 0;
			for (double q : my) sumlat += q;
			double pole = (sumlat > 0) ? 90.0 : -90.0;
			double first_lon = (mx.front() < 0) ? -180.0 : 180.0;
			double last_lon  = (mx.back()  < 0) ? -180.0 : 180.0;
			mx.front() = first_lon;
			mx.back()  = last_lon;
			mx.push_back(last_lon);  my.push_back(pole);
			mx.push_back(first_lon); my.push_back(pole);
			mx.push_back(mx[0]);     my.push_back(my[0]);
			geom.addPart(SpatPart(mx, my));
		}

		if (geom.parts.empty()) continue;
		geom.computeExtent();
		out.addGeom(geom);
		type_col.push_back(is_corner[v] ? "pentagon" : "hexagon");
	}

	out.setSRS({lonlat_crs});
	if (!type_col.empty()) {
		out.add_column(type_col, "type");
	}
	return out;
}



// rectangular tessellation on a sphere.
//
// "size" is the target cell edge length in meters. Three layouts are
// supported through the align argument:
//
//   align = 0  ("fit"): latitude bands and per-band cells are stretched/shrunk
//     to fill the extent exactly. Cell size deviates a little from "size"
//
//   align = 1  ("equal"): every cell has the same constant size dlat_ideal
//     in latitude and dlon_ideal_j = dlat_ideal/cos(lat_c) in longitude per
//     band, so cell area is approximately size^2 everywhere. 
//
//   align = 2  ("cube"): like align=1 but adjusts both dlat and dlon together
//     so the cells are square in plate carree
//
// For a global longitude extent the layout is snapped to a -180..180 grid so
// no cell straddles the antimeridian.

SpatVector SpatVector::rectangles_lonlat(SpatExtent e, double size, int align) {

	SpatVector out;
	if (size <= 0) {
		out.setError("size must be a positive number");
		return out;
	}
	if (!e.valid_notempty()) {
		out.setError("invalid extent");
		return out;
	}

	const std::string lonlat_crs = "+proj=longlat +datum=WGS84 +no_defs";
	const double R = 6378137.0;
	const double DEG = 180.0 / M_PI;
	const double dlat_ideal = (size / R) * DEG;        // degrees at equator

	// extent clamped to the valid lon/lat range
	double ymin = std::max(-90.0, e.ymin);
	double ymax = std::min( 90.0, e.ymax);
	double xmin = std::max(-180.0, e.xmin);
	double xmax = std::min( 180.0, e.xmax);
	if (ymax <= ymin || xmax <= xmin) {
		out.setSRS({lonlat_crs});
		return out;
	}

	double yrange = ymax - ymin;
	double xrange = xmax - xmin;
	bool full_lat = yrange >= 180.0 - 1e-9;
	bool full_lon = xrange >= 360.0 - 1e-9;

	// Build the list of latitude bands first (lat1, lat2). The longitude
	// layout within each band is then computed in the cell-emission loop,
	// because it depends on the band center latitude.
	struct Band { double lat1, lat2; };
	std::vector<Band> bands;

	if (align == 2) {
		// "cube" mode: bands anchored at the equator, growing in height
		// toward the poles to keep dlat = dlon at each band.
		auto cube_side = [&](double lat_deg) {
			double cl = std::cos(lat_deg / DEG);
			if (cl < 1e-9) cl = 1e-9;
			return dlat_ideal / std::sqrt(cl);
		};

		auto add_band = [&](double a, double b) {
			if (b <= ymin || a >= ymax) return;        // outside extent
			bands.push_back({a, b});
		};

		// northward from equator
		double top = 0.0;
		for (int it = 0; it < 100000 && top < 90.0; it++) {
			double bottom = top;
			double s0 = cube_side(bottom);
			double lc = std::min(89.999999, bottom + 0.5 * s0);
			double s1 = cube_side(lc);
			double new_top = bottom + s1;
			if (new_top > 90.0 || new_top - bottom < 1e-9) new_top = 90.0;
			add_band(bottom, new_top);
			if (new_top >= 90.0) break;
			top = new_top;
		}
		// southward from equator
		double bot = 0.0;
		for (int it = 0; it < 100000 && bot > -90.0; it++) {
			double tp = bot;
			double s0 = cube_side(tp);
			double lc = std::max(-89.999999, tp - 0.5 * s0);
			double s1 = cube_side(lc);
			double new_bot = tp - s1;
			if (new_bot < -90.0 || tp - new_bot < 1e-9) new_bot = -90.0;
			add_band(new_bot, tp);
			if (new_bot <= -90.0) break;
			bot = new_bot;
		}

	} else if (align == 1 && !full_lat) {
		// "equal" mode: constant-height bands, centered on the extent,
		// use ceil so the stack fully covers [ymin, ymax]
		int nbands = (int) std::ceil(yrange / dlat_ideal);
		if (nbands < 1) nbands = 1;
		double y_start = ymin + 0.5 * (yrange - nbands * dlat_ideal);
		bands.reserve((size_t) nbands);
		for (int j = 0; j < nbands; j++) {
			double a = y_start + j * dlat_ideal;
			double b = a + dlat_ideal;
			if (a < -90.0) a = -90.0;
			if (b >  90.0) b =  90.0;
			if (a < b) bands.push_back({a, b});
		}

	} else {
		// "fit" mode (and equal/cube fall-throughs for full latitude range):
		// latitude bands stretched to fit the extent exactly
		int nbands = (int) std::round(yrange / dlat_ideal);
		if (nbands < 1) nbands = 1;
		double dlat = yrange / nbands;
		bands.reserve((size_t) nbands);
		for (int j = 0; j < nbands; j++) {
			double a = ymin + j * dlat;
			double b = (j == nbands - 1) ? ymax : (a + dlat);
			bands.push_back({a, b});
		}
	}

	out.reserve(bands.size() * 16);

	for (size_t j = 0; j < bands.size(); j++) {
		double lat1 = bands[j].lat1;
		double lat2 = bands[j].lat2;
		double dlat_band = lat2 - lat1;
		double lat_c = 0.5 * (lat1 + lat2);

		double cosc = std::cos(lat_c / DEG);
		if (cosc < 0) cosc = 0;

		// ideal cell width in degrees of longitude
		double dlon_ideal;
		if (align == 2) {
			// cube: dlon = dlat at this band
			dlon_ideal = dlat_band;
		} else {
			dlon_ideal = (cosc > 1e-12)
				? (size / (R * cosc)) * DEG
				: 360.0;
		}
		if (dlon_ideal > 360.0) dlon_ideal = 360.0;

		int N;
		double dlon;
		double x_start;
		if (align == 0 || full_lon) {
			// fit-style longitude (also used for any global-longitude extent)
			N = (int) std::round(xrange / dlon_ideal);
			if (N < 1) N = 1;
			dlon = xrange / N;
			x_start = xmin;
		} else {
			// equal- or cube-style longitude: constant cell width, centered,
			// ceil so the row fully covers the extent
			N = (int) std::ceil(xrange / dlon_ideal);
			if (N < 1) N = 1;
			dlon = dlon_ideal;
			x_start = xmin + 0.5 * (xrange - N * dlon);
		}

		for (int i = 0; i < N; i++) {
			double lon1 = x_start + i * dlon;
			double lon2 = lon1 + dlon;
			// snap to the extent edge for fit / full-lon modes only
			if (i == N - 1 && (align == 0 || full_lon)) {
				lon2 = xmax;
			}

			std::vector<double> rx = {lon1, lon2, lon2, lon1, lon1};
			std::vector<double> ry = {lat1, lat1, lat2, lat2, lat1};

			SpatGeom geom(SpatPart(rx, ry), polygons);
			geom.computeExtent();
			out.addGeom(geom);
		}
	}

	out.setSRS({lonlat_crs});
	return out;
}
