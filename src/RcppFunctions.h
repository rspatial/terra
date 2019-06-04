// Copyright (c) 2018  Robert J. Hijmans
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

Rcpp::List getBlockSizeR(SpatRaster* r, unsigned n);
Rcpp::List getAttributes(SpatVector* v);
//bool setAttributes(SpatVector* v, Rcpp::List x, std::vector<std::string> names, std::vector<std::string> types);

Rcpp::NumericVector getGeometry(SpatVector* v);

SpatRaster rcppReclassify(SpatRaster* x, Rcpp::NumericMatrix rcl, unsigned right, bool lowest, bool othersNA, SpatOptions &opt);
//Rcpp::NumericMatrix rcppAdjacent(SpatRaster* x, std::vector<double> cells, std::string directions, bool include);
