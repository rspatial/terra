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
		void resize(size_t n) {x.resize(n);}
		void resize(size_t n, SpatTime_t v) {x.resize(n, v);}
		void reserve(size_t n) {x.reserve(n);}
		void push_back(SpatTime_t v) {x.push_back(v);}
};


std::vector<int> get_date(SpatTime_t x);
std::vector<int> getymd(std::string s);
SpatTime_t get_time_string(std::string s);
SpatTime_t time_from_day(int syear, int smonth, int sday, double ndays);
SpatTime_t time_from_day_noleap(int syear, int smonth, int sday, double ndays);
SpatTime_t time_from_day_360(int syear, int smonth, int sday, double ndays);
SpatTime_t time_from_hour(int syear, int smonth, int sday, double nhours);
void hours_to_time(std::vector<SpatTime_t> &time, std::string origin);
SpatTime_t parse_time(std::string x);

#endif
