// Copyright (c) 2018-2019  Robert J. Hijmans
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


#include <vector>
#include <string>
#include "spatMessages.h"


class SpatDataFrame {
	public:
		SpatMessages msg;
		void setError(std::string s) { msg.setError(s); }
		void addWarning(std::string s) { msg.addWarning(s); }
		bool hasError() { return msg.has_error; }
	
		std::vector<std::string> names;
		std::vector<unsigned> itype;
		std::vector<unsigned> iplace;
		std::vector< std::vector<double>> dv;
		std::vector< std::vector<long>> iv;
		std::vector< std::vector<std::string>> sv;
		std::string NAS = "____NA_+";
		
		unsigned nrow();
		unsigned ncol();
		SpatDataFrame subset_rows(std::vector<unsigned> range);
		SpatDataFrame subset_cols(std::vector<unsigned> range);
		SpatDataFrame subset_rows(unsigned i);
		SpatDataFrame subset_cols(unsigned i);
		std::vector<double> getD(unsigned i);
		std::vector<long> getI(unsigned i);
		std::vector<std::string> getS(unsigned i);
	
		void add_row();
		
		void add_column(unsigned dtype, std::string name);
		bool add_column(std::vector<double> x, std::string name);
		bool add_column(std::vector<long> x, std::string name);
		bool add_column(std::vector<std::string> x, std::string name);
		
		void insert_column(std::vector<double>);
		void insert_column(std::vector<long>);		
		void insert_column(std::vector<std::string>);

		void resize_rows(unsigned n);
		void resize_cols(unsigned n);
		void reserve(unsigned n);
		
		bool rbind(SpatDataFrame &x);
		bool cbind(SpatDataFrame &x);
		
		std::vector<std::string> get_names();
		void set_names(std::vector<std::string> nms);
};

