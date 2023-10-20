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

#ifndef SPATTIME_GUARD
#define SPATTIME_GUARD

//#include <vector>
//#include <string>
typedef long long SpatTime_t;

class SpatTime_v {
	public:
		std::vector<SpatTime_t> x;
		std::string zone;
		std::string step;
		size_t size() { return(x.size());}
		bool empty() { return(x.empty());}
		void resize(size_t n) {x.resize(n);}
		void resize(size_t n, SpatTime_t v) {x.resize(n, v);}
		void reserve(size_t n) {x.reserve(n);}
		void push_back(SpatTime_t v) {x.push_back(v);}
};

/*
class SpatTimeP {
	long year;
	short int month;
	short int day;
	short int hour;
	short int minute;
	short int second;
};

class SpatTimeP_v {
	public:
		std::vector<SpatTimeP> x;
		std::string zone;
		std::string step;
		size_t size() { return(x.size());}
		bool empty() { return(x.empty());}
		void resize(size_t n) {x.resize(n);}
		void resize(size_t n, SpatTimeP v) {x.resize(n, v);}
		void reserve(size_t n) {x.reserve(n);}
		void push_back(SpatTimeP v) {x.push_back(v);}
};
*/

SpatTime_t get_time(long year, unsigned month, unsigned day, int hr, int min, int sec);
std::vector<int> get_date(SpatTime_t x);
std::vector<int> getymd(std::string s);
int getyear(std::string s);
SpatTime_t get_time_string(std::string s);
SpatTime_t get_time_noleap(int syear, int smonth, int sday, int shour, int smin, int ssec, double n, std::string step);
//SpatTime_t time_from_day_noleap(int syear, int smonth, int sday, double ndays);

SpatTime_t time_from_day(int syear, int smonth, int sday, double ndays);
SpatTime_t time_from_day_360(int syear, int smonth, int sday, double ndays);
SpatTime_t time_from_hour(int syear, int smonth, int sday, int shour, double nhours);
void hours_to_time(std::vector<SpatTime_t> &time, std::string origin);
SpatTime_t parse_time(std::string x);

#endif
