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

	// Recompute the bounding box from node positions.
	void computeExtent();

	// Per-node degree (number of incident edges, counted once per
	// occurrence so a self-loop contributes 2).
	std::vector<size_t> node_degree();

	// Per-edge planar Euclidean length, computed from the stored
	// coordinate sequences.
	std::vector<double> edge_lengths_plane();

	// Round-trips
	// as_nodes(): SpatVector of points, one per node, with a "degree" column.
	// as_edges(): SpatVector of lines, one per edge, with from_node,
	//   to_node, source_id, length columns plus the forwarded
	//   per-edge attributes from edge_df.
	SpatVector as_nodes();
	SpatVector as_edges();

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
