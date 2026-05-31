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

#include <cmath>

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


std::vector<double> SpatNetwork::edge_lengths_plane() {
	std::vector<double> out(edge_from.size(), 0.0);
	for (size_t i = 0; i < edge_from.size(); i++) {
		const std::vector<double> &xs = edge_x[i];
		const std::vector<double> &ys = edge_y[i];
		double L = 0.0;
		for (size_t j = 1; j < xs.size(); j++) {
			double dx = xs[j] - xs[j-1];
			double dy = ys[j] - ys[j-1];
			L += std::sqrt(dx*dx + dy*dy);
		}
		out[i] = L;
	}
	return out;
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
	return out;
}


SpatVector SpatNetwork::as_edges() {
	SpatVector out;
	out.srs = srs;
	size_t n = edge_from.size();
	out.reserve(n);
	std::vector<long> from(n), to(n), src(n);
	std::vector<double> len = edge_lengths_plane();
	for (size_t i = 0; i < n; i++) {
		SpatPart p(edge_x[i], edge_y[i]);
		SpatGeom g(p, lines);
		out.addGeom(g);
		from[i] = (long) edge_from[i] + 1;        // 1-based for R users
		to[i]   = (long) edge_to[i]   + 1;
		src[i]  = edge_source[i] >= 0 ? edge_source[i] + 1 : -1;
	}
	out.df.add_column(from, "from_node");
	out.df.add_column(to,   "to_node");
	out.df.add_column(src,  "source_id");
	out.df.add_column(len,  "length");

	// Forward per-edge attributes from edge_df, if present.
	if (edge_df.nrow() == n && edge_df.ncol() > 0) {
		out.df.cbind(edge_df);
	}
	return out;
}
