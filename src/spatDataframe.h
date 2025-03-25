// Copyright (c) 2018-2023  Robert J. Hijmans
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


#ifndef SPATDATAFRAME_GUARD
#define SPATDATAFRAME_GUARD

#include <vector>
#include <string>
//#include "spatMessages.h"
#include "spatBase.h"
#include "spatTime.h"
#include "spatFactor.h"
#include <cstdint>
#include <limits>

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
		std::vector<std::string> getWarnings() { return msg.getWarnings(); }
		std::string getError() { return msg.getError(); }
	
		std::vector<std::string> names;
		std::vector<size_t> itype; //0 double, 1 long, 2 string, 3 bool, 4 time, 5 factor
		std::vector<size_t> iplace;
		std::vector<std::vector<double>> dv;
		std::vector<std::vector<long>> iv;
		std::vector<std::vector<std::string>> sv;
		std::vector<std::vector<int8_t>> bv;
		std::vector<SpatTime_v> tv;
		std::vector<SpatFactor> fv;		
		std::string NAS = "____NA_+";
		long NAL = std::numeric_limits<long>::min();
		SpatTime_t NAT = std::numeric_limits<SpatTime_t>::min();
		
		size_t nrow();
		size_t ncol();
		SpatDataFrame subset_rows(std::vector<long> range);
		SpatDataFrame subset_rows(std::vector<size_t> range);
		SpatDataFrame subset_cols(std::vector<size_t> range);
		SpatDataFrame subset_rows(size_t i);
		SpatDataFrame subset_cols(size_t i);
		std::vector<double> getD(size_t i);
		std::vector<long> getI(size_t i);
		std::vector<std::string> getS(size_t i);
		std::vector<int8_t> getB(size_t i);
		SpatTime_v getT(size_t i);
		SpatFactor getF(size_t i);

		std::vector<std::string> as_string(size_t v);
		std::vector<long> as_long(size_t v);
		std::vector<double> as_double(size_t v);

		double getDvalue(size_t i, size_t j);
		long getIvalue(size_t i, size_t j);
		std::string getSvalue(size_t i, size_t j);
		int8_t getBvalue(size_t i, size_t j);
		SpatTime_t getTvalue(size_t i, size_t j);
		SpatFactor getFvalue(size_t i, size_t j);
	
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
		bool add_column(SpatTime_v x, std::string name);
		bool add_column(SpatFactor x, std::string name);
		bool add_column_bool(std::vector<int> x, std::string name);
		bool add_column_bool(std::vector<bool> x, std::string name);
		bool add_column_time(std::vector<SpatTime_t> x, std::string name, std::string step, std::string zone);
		void insert_column(std::vector<double>, size_t i);
		void insert_column(std::vector<long>, size_t i);		
		void insert_column(std::vector<std::string>, size_t i);
		void insert_column(std::vector<int8_t>, size_t i);
		void insert_column(SpatTime_v, size_t i);
		void insert_column(SpatFactor, size_t i);

		bool remove_column(std::string field);
		bool remove_column(int i);		

		void resize_rows(size_t n);
		void remove_rows(std::vector<size_t> r);

		void resize_cols(size_t n);
		void reserve(size_t n);
		
		bool rbind(SpatDataFrame &x);
		bool cbind(SpatDataFrame &x);

		SpatDataFrame unique_col(int col);
		std::vector<int> getIndex(int col, SpatDataFrame &x);

		std::vector<std::string> get_names();
		void set_names(std::vector<std::string> nms);
		
		std::vector<std::string> get_datatypes();	
		std::string get_datatype(std::string field);
		std::string get_datatype(int field);
		int get_fieldindex(std::string field);

		std::vector<std::string> get_timesteps();	
		std::vector<std::string> get_timezones();	

		bool field_exists(std::string field);
		bool write_dbf(std::string filename, bool overwrite, SpatOptions &opt);

		std::vector<std::vector<std::string>> to_strings();
		std::vector<std::string> one_string();
		SpatDataFrame unique();
		size_t strwidth(size_t i);

		SpatDataFrame sortby(std::string field, bool descending);
};

#endif //SPATDATAFRAME_GUARD

