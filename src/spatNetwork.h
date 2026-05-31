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

#ifndef SPATNETWORK_GUARD
#define SPATNETWORK_GUARD

#include "spatVector.h"

// graph derived from a SpatVector of lines (e.g. a road network).

// Layout:
//   - Nodes: node_x[k], node_y[k] for k in 0 .. nnodes()-1
//   - Edges: for i in 0 .. nedges()-1
//       edge_from[i], edge_to[i] : indices into node_x/node_y
//       edge_x[i], edge_y[i]     : coord sequence of the edge,
//                                  including the start and end node
//       edge_source[i] : index into the source SpatVector this edge fragment came
//                        from, or -1
//
// Attribute forwarding: edge_df has one row per edge, copied from the
// row of the source SpatVector pointed to by edge_source[i].

class SpatNetwork {
public:
	std::vector<double> node_x, node_y;
	std::vector<size_t> edge_from, edge_to;
	std::vector<std::vector<double>> edge_x, edge_y;
	std::vector<long> edge_source;

	// Geometric length of each edge, in CRS units (or meters when the
	// CRS is geographic / lon-lat). Cached at construction time so that
	// `weight` can be replaced by something else (e.g. travel time)
	// without losing the underlying geometric length.
	std::vector<double> edge_length;

	// Per-edge weight used by graph algorithms (shortest path, etc.).
	// Defaults to `edge_length`. Use `setWeights()` to override or
	// `clearWeights()` to mark the network as unweighted.
	std::vector<double> edge_weight;
	bool weighted = false;

	// If true, edge orientation (from_node -> to_node) is the legal
	// travel direction; an undirected graph treats every edge as
	// bidirectional. Use this for water flow, one-way streets, etc.
	bool directed = false;

	SpatSRS srs;
	SpatExtent extent;
	SpatDataFrame edge_df;
	SpatMessages msg;

	SpatNetwork() {}
	virtual ~SpatNetwork() {}
	SpatNetwork deepCopy() { return *this; }

	size_t nnodes() { return node_x.size(); }
	size_t nedges() { return edge_from.size(); }
	bool empty() { return node_x.empty() && edge_from.empty(); }

	bool isDirected() { return directed; }
	void setDirected(bool d) { directed = d; }

	bool isWeighted() { return weighted; }
	std::vector<double> getWeights();         // edge_weight if set, else edge_length
	bool setWeights(std::vector<double> w);
	void clearWeights() { edge_weight.clear(); weighted = false; }

	// Recompute the bounding box from node positions.
	void computeExtent();

	// Per-node degree (number of incident edges, counted once per
	// occurrence so a self-loop contributes 2). For directed networks
	// this is total degree; in_degree / out_degree are separate.
	std::vector<size_t> node_degree();
	std::vector<size_t> node_in_degree();
	std::vector<size_t> node_out_degree();

	// Recompute and store per-edge geometric lengths from the current
	// edge_x / edge_y. Uses geodesic distance (in meters) when the CRS
	// is lon-lat, planar Euclidean (scaled to meters when possible)
	// otherwise.
	void compute_edge_lengths();

	// Build a network directly from primitive arrays (used by R-side
	// constructors that don't go through GEOS noding, e.g. coercion
	// from an igraph). Each edge is materialised as a straight 2-point
	// line between its two end nodes; lengths are recomputed and (if
	// `weighted`) stored as the edge weight. `edge_from`/`edge_to` are
	// 0-based indices into the node arrays. `weights` may be empty
	// (then ignored), or must have length equal to the edge count.
	bool buildFromComponents(std::vector<double> nx,
	                         std::vector<double> ny,
	                         std::vector<size_t> efrom,
	                         std::vector<size_t> eto,
	                         std::vector<double> w,
	                         bool is_directed);

	// Round-trip with GDAL's Geographic Network Model (GNM). The on-disk
	// layout is a directory (driver "GNMFile") or database (driver
	// "GNMDatabase"); the network's nodes go into a `nodes` point layer,
	// edges into an `edges` line layer, and topology into the system
	// `graph` layer that GNM manages internally. Returns false on any
	// failure (the network's error/warning state carries the details).
	// `options` accept GNM-specific create options as `"key=value"`
	// strings (e.g. `"FORMAT=ESRI Shapefile"`).
	bool write_gnm(std::string filename,
	               std::string driver_name,
	               std::vector<std::string> options);
	bool read_gnm(std::string filename);

	// Round-trips
	// as_nodes(): SpatVector of points, one per node, with a "degree" column.
	// as_edges(): SpatVector of lines, one per edge, with from_node,
	//   to_node, source_id, length, weight columns plus the forwarded
	//   per-edge attributes from edge_df.
	SpatVector as_nodes();
	SpatVector as_edges();

	// Shortest paths between pairs of node indices (0-based) using
	// Dijkstra on the current edge weights (edge_weight if set, else
	// edge_length). Pairs are taken element-wise; either side may be
	// length 1 to be recycled. Each output feature is a line that
	// follows the original (sinuous) edge geometries; an unreachable
	// pair gets an empty geometry and an NA distance. Output columns
	// are `from`, `to`, and `distance`.
	SpatVector shortest_paths(std::vector<size_t> from_nodes,
	                          std::vector<size_t> to_nodes);

	bool setSRS(std::string _srs) {
		std::string m;
		if (!srs.set(_srs, m)) {
			addWarning("cannot set SRS: " + m);
			return false;
		}
		return true;
	}
	std::string getSRS(std::string x) { return srs.get(x); }

	SpatExtent getExtent() { return extent; }

	void setError(std::string s) { msg.setError(s); }
	void addWarning(std::string s) { msg.addWarning(s); }
	bool hasError() { return msg.has_error; }
	bool hasWarning() { return msg.has_warning; }
	std::vector<std::string> getWarnings() { return msg.getWarnings(); }
	std::string getError() { return msg.getError(); }

	std::string show();
};

#endif
