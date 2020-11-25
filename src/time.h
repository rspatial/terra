#ifndef SPATTIME_GUARD
#define SPATTIME_GUARD

typedef long long SpatTime_t;
SpatTime_t get_time(long year, unsigned month, unsigned day=15, unsigned hr=0, unsigned min=0, unsigned sec=0);
std::vector<int> get_date(SpatTime_t x);
std::vector<int> getymd(std::string s);
SpatTime_t get_time_string(std::string s);
SpatTime_t time_from_day(int syear, int smonth, int sday, int ndays);
SpatTime_t time_from_day_noleap(int syear, int smonth, int sday, int ndays);
SpatTime_t time_from_day_360(int syear, int smonth, int sday, int ndays);
SpatTime_t time_from_hour(int syear, int smonth, int sday, int nhours);
void hours_to_time(std::vector<SpatTime_t> &time, std::string origin);

#endif
