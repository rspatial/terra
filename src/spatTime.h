#ifndef SPATTIME_GUARD
#define SPATTIME_GUARD

//#include <vector>
//#include <string>
typedef long long SpatTime_t;
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
