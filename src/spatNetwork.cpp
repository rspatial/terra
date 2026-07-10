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

#include "spatNetwork.h"
#include "distance.h"
#include "dijkstra.h"

#include <algorithm>
#include <cmath>
#include <limits>


void SpatNetwork::computeExtent() {
	if (node_x.empty()) {
		extent = SpatExtent();
		return;
	}
	double xmn = node_x[0], xmx = node_x[0];
	double ymn = node_y[0], ymx = node_y[0];
	for (size_t i = 1; i < node_x.size(); i++) {
		if (node_x[i] < xmn) xmn = node_x[i];
		if (node_x[i] > xmx) xmx = node_x[i];
		if (node_y[i] < ymn) ymn = node_y[i];
		if (node_y[i] > ymx) ymx = node_y[i];
	}
	extent.xmin = xmn; extent.xmax = xmx;
	extent.ymin = ymn; extent.ymax = ymx;
}


std::vector<size_t> SpatNetwork::node_degree() {
	std::vector<size_t> d(node_x.size(), 0);
	for (size_t i = 0; i < edge_from.size(); i++) {
		if (edge_from[i] < d.size()) d[edge_from[i]]++;
		if (edge_to[i] < d.size()) d[edge_to[i]]++;
	}
	return d;
}


std::vector<size_t> SpatNetwork::node_in_degree() {
	std::vector<size_t> d(node_x.size(), 0);
	for (size_t i = 0; i < edge_to.size(); i++) {
		if (edge_to[i] < d.size()) d[edge_to[i]]++;
	}
	return d;
}


std::vector<size_t> SpatNetwork::node_out_degree() {
	std::vector<size_t> d(node_x.size(), 0);
	for (size_t i = 0; i < edge_from.size(); i++) {
		if (edge_from[i] < d.size()) d[edge_from[i]]++;
	}
	return d;
}


// Lon-lat-aware geometric length per edge.
//   - Lon/lat CRS (srs.to_meter() == 0): geodesic distance summed along
//     the polyline, in meters (WGS84 ellipsoid via geographiclib).
//   - Otherwise: planar Euclidean, scaled to meters when the CRS unit
//     is known (`to_meter`); 1.0 (no scaling) when unknown.
void SpatNetwork::compute_edge_lengths() {
	size_t ne = edge_from.size();
	edge_length.assign(ne, 0.0);
	if (ne == 0) return;

	double m = srs.to_meter();
	bool unknown_unit = std::isnan(m);
	if (unknown_unit) m = 1.0;
	bool lonlat = (m == 0.0);

	if (lonlat) {
		for (size_t i = 0; i < ne; i++) {
			const std::vector<double> &xs = edge_x[i];
			const std::vector<double> &ys = edge_y[i];
			double L = 0.0;
			for (size_t j = 1; j < xs.size(); j++) {
				L += distance_lonlat(xs[j-1], ys[j-1], xs[j], ys[j]);
			}
			edge_length[i] = L;
		}
	} else {
		for (size_t i = 0; i < ne; i++) {
			const std::vector<double> &xs = edge_x[i];
			const std::vector<double> &ys = edge_y[i];
			double L = 0.0;
			for (size_t j = 1; j < xs.size(); j++) {
				double dx = xs[j] - xs[j-1];
				double dy = ys[j] - ys[j-1];
				L += std::sqrt(dx*dx + dy*dy);
			}
			edge_length[i] = L * m;
		}
	}
}


bool SpatNetwork::buildFromComponents(std::vector<double> nx,
                                      std::vector<double> ny,
                                      std::vector<size_t> efrom,
                                      std::vector<size_t> eto,
                                      std::vector<double> w,
                                      bool is_directed) {
	if (nx.size() != ny.size()) {
		setError("node x/y must have the same length");
		return false;
	}
	if (efrom.size() != eto.size()) {
		setError("edge from/to must have the same length");
		return false;
	}
	size_t nn = nx.size();
	size_t ne = efrom.size();
	for (size_t i = 0; i < ne; i++) {
		if (efrom[i] >= nn || eto[i] >= nn) {
			setError("edge endpoint refers to a non-existent node");
			return false;
		}
	}

	node_x = std::move(nx);
	node_y = std::move(ny);
	edge_from = std::move(efrom);
	edge_to   = std::move(eto);
	edge_x.assign(ne, std::vector<double>(2));
	edge_y.assign(ne, std::vector<double>(2));
	edge_source.assign(ne, -1);
	for (size_t i = 0; i < ne; i++) {
		edge_x[i][0] = node_x[edge_from[i]];
		edge_x[i][1] = node_x[edge_to[i]];
		edge_y[i][0] = node_y[edge_from[i]];
		edge_y[i][1] = node_y[edge_to[i]];
	}
	edge_df = SpatDataFrame();
	directed = is_directed;
	compute_edge_lengths();
	if (w.size() == ne) {
		edge_weight = std::move(w);
		weighted = true;
	} else if (!w.empty()) {
		addWarning("weight vector ignored: length must equal the number of edges");
		edge_weight.clear();
		weighted = false;
	} else {
		edge_weight.clear();
		weighted = false;
	}
	computeExtent();
	return true;
}


std::vector<double> SpatNetwork::getWeights() {
	if (weighted && edge_weight.size() == edge_from.size()) return edge_weight;
	return edge_length;
}


bool SpatNetwork::setWeights(std::vector<double> w) {
	if (w.size() != edge_from.size()) {
		addWarning("weight length must equal the number of edges");
		return false;
	}
	edge_weight = std::move(w);
	weighted = true;
	return true;
}


SpatVector SpatNetwork::as_nodes() {
	SpatVector out;
	out.srs = srs;
	out.reserve(node_x.size());
	for (size_t i = 0; i < node_x.size(); i++) {
		SpatPart p(node_x[i], node_y[i]);
		SpatGeom g(p, points);
		out.addGeom(g);
	}
	std::vector<size_t> deg = node_degree();
	std::vector<long> deg_long(deg.begin(), deg.end());
	out.df.add_column(deg_long, "degree");
	if (directed) {
		std::vector<size_t> ind = node_in_degree();
		std::vector<size_t> outd = node_out_degree();
		std::vector<long> in_long(ind.begin(),  ind.end());
		std::vector<long> ou_long(outd.begin(), outd.end());
		out.df.add_column(in_long, "in_degree");
		out.df.add_column(ou_long, "out_degree");
	}
	return out;
}


SpatVector SpatNetwork::as_edges() {
	SpatVector out;
	out.srs = srs;
	size_t n = edge_from.size();
	out.reserve(n);
	std::vector<long> from(n), to(n), src(n);
	for (size_t i = 0; i < n; i++) {
		SpatPart p(edge_x[i], edge_y[i]);
		SpatGeom g(p, lines);
		out.addGeom(g);
		from[i] = (long) edge_from[i];
		to[i]   = (long) edge_to[i];
		src[i]  = edge_source[i];           // -1 means "no source" (no shift)
	}
	out.df.add_column(from, "from_node");
	out.df.add_column(to,   "to_node");
	out.df.add_column(src,  "source_id");

	// Always include the geometric length; include weight if it differs.
	std::vector<double> L = (edge_length.size() == n) ? edge_length : std::vector<double>(n, 0.0);
	out.df.add_column(L, "length");
	if (weighted && edge_weight.size() == n) {
		out.df.add_column(edge_weight, "weight");
	}

	// Forward per-edge attributes from edge_df, if present.
	if (edge_df.nrow() == n && edge_df.ncol() > 0) {
		out.df.cbind(edge_df);
	}
	return out;
}


// Build the adjacency list once and reuse it across (potentially many)
// SSSP calls. For directed networks each input edge contributes only
// the from->to direction; for undirected ones both directions are
// added. Negative weights are clamped to 0 (Dijkstra cannot handle
// negative cycles).
static spat_dijkstra::Adjacency build_adjacency(const SpatNetwork &net,
                                                const std::vector<double> &w,
                                                bool directed) {
	size_t nn = net.node_x.size();
	size_t ne = net.edge_from.size();
	spat_dijkstra::Adjacency adj(nn);
	for (size_t i = 0; i < ne; i++) {
		size_t f = net.edge_from[i];
		size_t t = net.edge_to[i];
		double we = (i < w.size()) ? w[i] : 1.0;
		if (std::isnan(we) || we < 0.0) we = 0.0;
		if (f < nn && t < nn) {
			adj[f].push_back({t, we, i});
			if (!directed) adj[t].push_back({f, we, i});
		}
	}
	return adj;
}


// Reconstruct the polyline that follows the original edge geometries
// from `s` to `t` using the predecessor arrays. `curnode` walks
// backwards from t to s; we then play the edge list forwards,
// reversing per-edge coordinates when the edge was traversed against
// its from->to orientation.
static void reconstruct_geometry(const SpatNetwork &net,
                                 size_t s, size_t t,
                                 const std::vector<size_t> &pred_node,
                                 const std::vector<size_t> &pred_edge,
                                 std::vector<double> &X,
                                 std::vector<double> &Y) {
	const size_t NPOS = std::numeric_limits<size_t>::max();
	X.clear();
	Y.clear();
	if (s == t) return;

	std::vector<size_t> edge_path;
	size_t cur = t;
	while (cur != s) {
		size_t e = pred_edge[cur];
		size_t p = pred_node[cur];
		if (e == NPOS || p == NPOS) return;   // unreachable
		edge_path.push_back(e);
		cur = p;
	}
	std::reverse(edge_path.begin(), edge_path.end());

	size_t curnode = s;
	for (size_t step = 0; step < edge_path.size(); step++) {
		size_t e = edge_path[step];
		const std::vector<double> &ex = net.edge_x[e];
		const std::vector<double> &ey = net.edge_y[e];
		bool reverse = (net.edge_from[e] != curnode);

		std::vector<double> ox, oy;
		ox.reserve(ex.size());
		oy.reserve(ey.size());
		if (reverse) {
			for (size_t j = ex.size(); j-- > 0; ) {
				ox.push_back(ex[j]);
				oy.push_back(ey[j]);
			}
		} else {
			ox = ex;
			oy = ey;
		}

		if (step == 0) {
			X = std::move(ox);
			Y = std::move(oy);
		} else {
			// Skip first coord (duplicate of previous edge's last coord).
			for (size_t j = 1; j < ox.size(); j++) {
				X.push_back(ox[j]);
				Y.push_back(oy[j]);
			}
		}
		curnode = reverse ? net.edge_from[e] : net.edge_to[e];
	}
}


SpatVector SpatNetwork::shortest_paths(std::vector<size_t> from_nodes,
                                       std::vector<size_t> to_nodes) {
	SpatVector out;
	out.srs = srs;

	if (from_nodes.empty() || to_nodes.empty()) {
		out.setError("'from' and 'to' must each have length >= 1");
		return out;
	}
	size_t nfrom = from_nodes.size();
	size_t nto   = to_nodes.size();
	size_t npairs = std::max(nfrom, nto);
	if (nfrom != npairs && nfrom != 1) {
		out.setError("length of 'from' must equal length of 'to' or be 1");
		return out;
	}
	if (nto != npairs && nto != 1) {
		out.setError("length of 'to' must equal length of 'from' or be 1");
		return out;
	}

	size_t nn = nnodes();
	if (nn == 0) {
		out.setError("network has no nodes");
		return out;
	}
	for (size_t v : from_nodes) {
		if (v >= nn) { out.setError("'from' contains an invalid node id"); return out; }
	}
	for (size_t v : to_nodes) {
		if (v >= nn) { out.setError("'to' contains an invalid node id"); return out; }
	}

	std::vector<double> w = getWeights();
	if (w.size() != edge_from.size()) {
		// No weights at all: treat every edge as unit weight (i.e. a
		// hop count). We materialise that here so the weight vector
		// matches the edge count.
		w.assign(edge_from.size(), 1.0);
	}

	spat_dijkstra::Adjacency adj = build_adjacency(*this, w, directed);

	// Cache one SSSP per unique source -- this is the common case
	// when the user wants distances from one source to many destinations.
	std::vector<bool> have(nn, false);
	std::vector<std::vector<double>> all_dist(nn);
	std::vector<std::vector<size_t>> all_predn(nn), all_prede(nn);

	std::vector<long>   path_from(npairs), path_to(npairs);
	std::vector<double> path_dist(npairs);
	std::vector<std::vector<double>> path_x(npairs), path_y(npairs);

	for (size_t k = 0; k < npairs; k++) {
		size_t s = (nfrom == 1) ? from_nodes[0] : from_nodes[k];
		size_t t = (nto   == 1) ? to_nodes[0]   : to_nodes[k];

		if (!have[s]) {
			spat_dijkstra::sssp(adj, s, all_dist[s], all_predn[s], all_prede[s]);
			have[s] = true;
		}
		const std::vector<double> &dist = all_dist[s];
		path_from[k] = (long) s;            // 0-based; the R wrapper shifts to 1-based
		path_to[k]   = (long) t;

		if (s == t) {
			path_dist[k] = 0.0;
			// Empty geometry; users can join on `from`/`to` if they need.
		} else if (std::isinf(dist[t])) {
			path_dist[k] = NAN;
		} else {
			path_dist[k] = dist[t];
			reconstruct_geometry(*this, s, t,
			                     all_predn[s], all_prede[s],
			                     path_x[k], path_y[k]);
		}
	}

	out.reserve(npairs);
	for (size_t k = 0; k < npairs; k++) {
		if (path_x[k].empty()) {
			SpatGeom g(lines);
			out.addGeom(g);
		} else {
			SpatPart p(path_x[k], path_y[k]);
			SpatGeom g(p, lines);
			out.addGeom(g);
		}
	}
	out.df.add_column(path_from, "from");
	out.df.add_column(path_to,   "to");
	out.df.add_column(path_dist, "distance");
	return out;
}
