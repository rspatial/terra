// Copyright (c) 2018-2020  Robert J. Hijmans
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
//#include "spatMessages.h"
#include "spatBase.h"
#include "spatTime.h"

class SpatDataFrame {
	public:
		SpatDataFrame();
		virtual ~SpatDataFrame(){}

		SpatDataFrame skeleton();
	
		SpatMessages msg;
		void setError(std::string s) { msg.setError(s); }
		void addWarning(std::string s) { msg.addWarning(s); }
		bool hasError() { return msg.has_error; }
		bool hasWarning() { return msg.has_warning; }
		std::string getWarnings() { return msg.getWarnings(); }
		std::string getError() { return msg.getError(); }
	
		std::vector<std::string> names;
		std::vector<unsigned> itype; //0 double, 1 long, 2 string, 3 bool, 4 time
		std::vector<unsigned> iplace;
		std::vector< std::vector<double>> dv;
		std::vector< std::vector<long>> iv;
		std::vector< std::vector<std::string>> sv;
		std::vector< std::vector<int8_t>> bv;
		std::vector< std::vector<SpatTime_t>> tv;
		std::string NAS = "____NA_+";
		
		unsigned nrow();
		unsigned ncol();
		SpatDataFrame subset_rows(std::vector<long> range);
		SpatDataFrame subset_rows(std::vector<unsigned> range);
		SpatDataFrame subset_cols(std::vector<unsigned> range);
		SpatDataFrame subset_rows(unsigned i);
		SpatDataFrame subset_cols(unsigned i);
		std::vector<double> getD(unsigned i);
		std::vector<long> getI(unsigned i);
		std::vector<std::string> getS(unsigned i);
		std::vector<int8_t> getB(unsigned i);
		std::vector<SpatTime_t> getT(unsigned i);

		std::vector<std::string> as_string(size_t v);
		std::vector<long> as_long(size_t v);
		std::vector<double> as_double(size_t v);

		double getDvalue(unsigned i, unsigned j);
		long getIvalue(unsigned i, unsigned j);
		std::string getSvalue(unsigned i, unsigned j);
		int8_t getBvalue(unsigned i, unsigned j);
		SpatTime_t getTvalue(unsigned i, unsigned j);
	
		void add_row();
		void add_rows(size_t n);
		
		//void set_values(std::vector<double> x, std::string name);
		//void set_values(std::vector<long> x, std::string name);
		//void set_values(std::vector<std::string> x, std::string name);

		void add_column(unsigned dtype, std::string name);
		bool add_column(std::vector<double> x, std::string name);
		bool add_column(std::vector<long> x, std::string name);
		bool add_column(std::vector<int> x, std::string name);
		bool add_column(std::vector<std::string> x, std::string name);
		bool add_column(std::vector<int8_t> x, std::string name);
		bool add_column(std::vector<SpatTime_t> x, std::string name);
		bool add_column_bool(std::vector<int> x, std::string name);
		
		void insert_column(std::vector<double>, size_t i);
		void insert_column(std::vector<long>, size_t i);		
		void insert_column(std::vector<std::string>, size_t i);
		void insert_column(std::vector<int8_t>, size_t i);
		void insert_column(std::vector<SpatTime_t>, size_t i);

		bool remove_column(std::string field);
		bool remove_column(int i);		

		void resize_rows(unsigned n);
		void remove_rows(std::vector<unsigned> r);

		void resize_cols(unsigned n);
		void reserve(unsigned n);
		
		bool rbind(SpatDataFrame &x);
		bool cbind(SpatDataFrame &x);

		SpatDataFrame unique(int col);
		std::vector<int> getIndex(int col, SpatDataFrame &x);

		std::vector<std::string> get_names();
		void set_names(std::vector<std::string> nms);
		
		std::vector<std::string> get_datatypes();	
		std::string get_datatype(std::string field);
		std::string get_datatype(int field);
		int get_fieldindex(std::string field);

		bool field_exists(std::string field);
		bool write_dbf(std::string filename, bool overwrite, SpatOptions &opt);
};

