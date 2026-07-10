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

// Generic single-source shortest-path Dijkstra on a graph given as an
// adjacency list. Used by SpatNetwork::shortest_paths() and (in
// principle) reusable by any other code that needs to operate on an
// abstract graph -- the lattice version in costDistance.cpp uses an
// implicit grid graph and is therefore separate.

#ifndef SPAT_DIJKSTRA_H
#define SPAT_DIJKSTRA_H

#include <vector>
#include <queue>
#include <limits>
#include <cstddef>
#include <utility>

namespace spat_dijkstra {

struct Edge {
	size_t to;
	double weight;
	size_t edge_id;   // index back into the original edge list
};

using Adjacency = std::vector<std::vector<Edge>>;

// Single-source shortest path. After return:
//   - dist[v]      : weight of the shortest path source -> v, or +Inf
//                    if v is unreachable.
//   - pred_node[v] : predecessor node of v on a shortest path (or
//                    NPOS if v is unreachable / is the source).
//   - pred_edge[v] : edge_id of the edge used to reach v from
//                    pred_node[v] (or NPOS).
inline void sssp(const Adjacency &adj,
                 size_t source,
                 std::vector<double> &dist,
                 std::vector<size_t> &pred_node,
                 std::vector<size_t> &pred_edge) {

	const size_t NPOS = std::numeric_limits<size_t>::max();
	const double INF  = std::numeric_limits<double>::infinity();

	size_t n = adj.size();
	dist.assign(n, INF);
	pred_node.assign(n, NPOS);
	pred_edge.assign(n, NPOS);
	if (source >= n) return;

	dist[source] = 0.0;
	using Entry = std::pair<double, size_t>;
	std::priority_queue<Entry, std::vector<Entry>, std::greater<Entry>> pq;
	pq.push({0.0, source});

	while (!pq.empty()) {
		double d = pq.top().first;
		size_t u = pq.top().second;
		pq.pop();
		if (d > dist[u]) continue;        // stale entry
		const std::vector<Edge> &out = adj[u];
		for (size_t i = 0; i < out.size(); i++) {
			const Edge &e = out[i];
			double nd = d + e.weight;
			if (nd < dist[e.to]) {
				dist[e.to]      = nd;
				pred_node[e.to] = u;
				pred_edge[e.to] = e.edge_id;
				pq.push({nd, e.to});
			}
		}
	}
}

}  // namespace spat_dijkstra

#endif  // SPAT_DIJKSTRA_H
