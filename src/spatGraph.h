// Copyright (c) 2018-2022  Robert J. Hijmans
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
		
class SpatGraph {
	public:
		virtual ~SpatGraph(){}
		std::vector<double> x;
		std::vector<double> y;
		std::vector<size_t> index;
		std::vector<size_t> edges;
		SpatDataFrame atts;
		std::string crs;
		
		SpatGraph() {};
		SpatGraph(std::vector<double> nx, std::vector<double> ny, std::vector<size_t> from, std::vector<size_t> to);

		void set_values(SpatDataFrame d);
		SpatDataFrame get_values();
		
		SpatGraph clean();
		bool writeGraph(SpatOptions opt);
		SpatGraph readGraph(SpatOptions opt);

		SpatVector shortestPath();
};

