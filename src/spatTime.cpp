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

#include <algorithm>
#include <vector>
#include <string>
#include <cmath>
#include "string_utils.h"


typedef long long SpatTime_t;


bool isleap(const long &year) {
	return (year % 4 == 0) && ((year % 400 == 0) || (year % 100 != 0 ));
}


SpatTime_t yeartime(const long &year) {
	// seconds per year
	// 365 * 24 * 3600 = 31536000
	return isleap(year) ? 31622400 : 31536000;
}



SpatTime_t get_time(long year, unsigned month, unsigned day, int hr, int min, int sec) {

    static const unsigned mdays[2][12] = {
        {0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334},
        {0, 31, 60, 91, 121, 152, 182, 213, 244, 274, 305, 335}
    };

	if (month > 12) {
		year += month / 12;
		month = ((month-1) % 12) + 1;
	}

	// the first day does not count, we start at 1970-01-01
    // 24 * 3600 = 86400
	SpatTime_t time = -86400;

	if (year < 1970) {
		for (long y = year; y < 1970; y++){
			time -= yeartime(y);
		}
	} else {
		for (long y = 1970; y < year; y++) {
			time += yeartime(y);
		}
	}

	time += (mdays[isleap(year)][month-1] + day) * 86400;
	time += (hr * 3600) + (min * 60) + sec;
    return time;
}


SpatTime_t get_time_str(std::vector<std::string> s) {
	std::vector<long> d(6, 0);
	for (size_t i=0; i<s.size(); i++) {
		d[i] = stol(s[i]);
	}
	return get_time(d[0], d[1], d[2], d[3], d[4], d[5]);
}


std::vector<int> get_date(SpatTime_t x) {

	// seconds per month, shifted one month
    static const unsigned secdays[2][13] = {
		{0, 2678400, 5097600, 7776000, 10368000, 13046400, 15638400, 18316800, 20995200, 23587200, 26265600, 28857600, 31536000},
		{0, 2678400, 5184000, 7862400, 10454400, 13132800, 15724800, 18403200, 21081600, 23673600, 26352000, 28944000, 31622400}
	};

	// the first day does not count, we start at 1970-01-01
	int year = 1970;
	if (x < 0) {
		while (x < 0) {
			year--;
			x += yeartime(year);
		}
	} else if (x > 0) {
		while (x >= 0) {
			x -= yeartime(year);
			year++;
		}
		year--;
		x += yeartime(year);
	}
	int month;
	int leap = isleap(year);
	for (month=1; month<13; month++) {
		if (x < (secdays[leap][month])) {
			break;
		}
	}
	x -= secdays[isleap(year)][month-1];

	int day = x / 86400 + 1;
	x = x % 86400;
	int hour = x / 3600;
	x = x % 3600;
	int min = x / 60;
	int sec = x % 60;

	std::vector<int> out= {year, month, day, hour, min, sec};
    return out;
}



std::vector<std::string> splitstr(std::string s, std::string delimiter){
	std::vector<std::string> out;
	size_t pos = 0;
	while ((pos = s.find(delimiter)) != std::string::npos) {
		out.push_back(s.substr(0, pos));
		s.erase(0, pos + delimiter.length());
	}
	out.push_back(s.substr(0, pos));
	return out;
}

void replace_one_char(std::string& s, char from, char to) {
	for (size_t i = 0; i < s.size(); i++) {
		if (s[i] == from) {
			s[i] = to;
		}
	}
}


int getyear(std::string s) {
	int y;
	try {
		y = stoi(s);
	} catch(...) {
		y = 1970;
	}
	return y;
}


std::vector<int> getymd(std::string s) {

//	s = std::regex_replace(s, std::regex("T"), " ");
	lowercase(s);
	replace_one_char(s, 't', ' ');

	size_t ncolon = std::count(s.begin(), s.end(), ':');
	std::vector<std::string> x;
	std::vector<std::string> y;
	if (ncolon > 0) {
		x = splitstr(s, " ");
		s = x[0];
		if (x.size() > 1) {
			x[1].erase(std::remove(x[1].begin(), x[1].end(), 'z'), x[1].end());
			y = splitstr(x[1], ":");
		}
	}

	size_t ndash = std::count(s.begin(), s.end(), '-');
	if (ndash == 2) {
		x = splitstr(s, "-");
	}
	x.insert(x.end(), y.begin(), y.end() );
	std::vector<int> out(x.size());

	try {
		for (size_t i=0; i<out.size(); i++){
			out[i] = std::stoi(x[i]);
		}
		out.resize(6, 0);
	} catch(...) {
		out = std::vector<int>(6);
	}
	
	return out;
}



SpatTime_t get_time_string(std::string s) {

/*
	std::vector<std::string> ss;
	size_t ncolon = std::count(s.begin(), s.end(), ':');
	if (ncolon > 0) {
		ss = splitstr(s, " ");
		s = ss[0];
	}
	size_t ndash = std::count(s.begin(), s.end(), '-');
	SpatTime_t time = 0;
	if (ndash == 2) {
		ss = splitstr(s, "-");
	} else {
		return time;
	}
*/

	std::vector<int> d = getymd(s);
	return get_time(d[0], d[1], d[2], d[3], d[4], d[5]);

}

SpatTime_t time_from_hour(int syear, int smonth, int sday, int shour, double nhours) {
	SpatTime_t time = get_time(syear, smonth, sday, shour, 0, 0);
	time += nhours * 3600;
	return time;
}

void hours_to_time(std::vector<SpatTime_t> &time, std::string origin) {
	std::vector<int> ymd = getymd(origin);
	SpatTime_t otime = get_time(ymd[0], ymd[1], ymd[2], 0, 0, 0);
	for (SpatTime_t &d : time) d = otime + d * 3600;
}



SpatTime_t time_from_day(int syear, int smonth, int sday, double ndays) {
	SpatTime_t time = get_time(syear, smonth, sday, 0, 0, 0);
	time += ndays * 86400;
	return time;
}



SpatTime_t get_time_noleap(int syear, int smonth, int sday, int shour, int smin, int ssec, double n, std::string step) {

	static const int md[13] =  {0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334, 365 };

	// set start to beginning of year 
	double s = ssec + smin * 60 + shour * 3600 + (sday-1) * 24 * 3600;
	for (int i=0; i < smonth; i++) {
		s += (md[i] * 24 * 3600);
	}

	double ndays;
	if (step == "hours") {
		ndays = (n + s/3600) / 24;
	} else if (step == "minutes") {
		n += s/60;
		ndays = n / 1440;
	} else if (step == "seconds") {
		ndays = (n+s) / 86400;
	} else if (step == "days") {
		ndays = n + s/86400;
	} else {
		return 0;
	}

	int year = ndays / 365;
	double rem = ndays - year * 365;
	int month;
	for (month=1; month<13; month++) {
		if (rem < md[month]) {
			break;
		}
	}
	rem -= md[month-1];
	int day = rem;
	rem -= day;
	day++;
	rem *= 24;
	int hr = rem;
	rem -= hr;
	int mn = rem * 60;
	rem -= mn;
	int sc = rem * 60;
	
	return get_time(year+syear, month, day, hr, mn, sc);
}



SpatTime_t get_time_360(int syear, int smonth, int sday, int shour, int smin, int ssec, double n, std::string step) {

	// set start to beginning of year 
	double s = ssec + smin * 60 + shour * 3600 + (sday-1) * 24 * 3600 + (smonth-1) * 30;

	double ndays;
	if (step == "hours") {
		ndays = (n + s/3600) / 24;
	} else if (step == "minutes") {
		n += s/60;
		ndays = n / 1440;
	} else if (step == "seconds") {
		ndays = (n+s) / 86400;
	} else if (step == "days") {
		ndays = n + s/86400;
	} else {
		return 0;
	}

	int year = ndays / 360;
	double rem = ndays - year * 360;
	int month = rem / 30;
	rem -= month * 30;
	month++;
	int day = rem;
	rem -= day;
	day++;
	rem *= 24;
	int hr = rem;
	rem -= hr;
	int mn = rem * 60;
	rem -= mn;
	int sc = rem * 60;
	
	return get_time(year+syear, month, day, hr, mn, sc);
}


SpatTime_t time_from_day_360(int syear, int smonth, int sday, double ndays) {
    static const int md[13] =  {0, 30, 60, 90, 120, 150, 180, 210, 240, 270, 300, 330, 360 };
	int year = ndays / 360;
	//int doy = ndays % 360;
	int doy = ndays - (year * 360);
	int month;
	for (month=1; month<13; month++) {
		if (doy < md[month]) {
			break;
		}
	}
	month--;
	int day = doy - md[month];
	SpatTime_t time = get_time(year+syear, month+smonth, day+sday, 0, 0, 0);
	return time;
}



SpatTime_t parse_time(std::string x) {
	lrtrim(x);
	std::vector<std::string> s = strsplit(x, " ");

	std::vector<std::string> time;
	if ((!s[0].empty()) && (s[0].substr(0, 1) != "-")) {
		time = strsplit(s[0], "-");
	} else {
		time = {s[0]};		
	}

	if (time.size() == 1) {
		return stoll(time[0]);
	} else if (time.size() != 3) {
		return 0;
	}

	if (s.size() > 1) {
		std::vector<std::string> secs = strsplit(s[1], ":");
		if (secs.size() == 3) {
			time.insert(time.end(), secs.begin(), secs.end());
		}
	}

	return get_time_str(time);
}

