#include<string>
#include<vector>

std::string double_to_string(double x);
std::vector<std::string> double_to_string(const std::vector<double> &x, std::string prep);

std::string concatenate(std::vector<std::string> v, std::string delim);
void lowercase(std::string &s);
bool is_in_set(std::string s, std::vector<std::string> ss);
std::string is_in_set_default(std::string s, std::vector<std::string> ss, std::string defvalue, bool tolower);
unsigned where_in_set(std::string s, std::vector<std::string> ss);

std::vector<std::string> strsplit(std::string s, std::string delimiter);
std::vector<double> str2dbl(std::vector<std::string> s);
std::vector<std::string> dbl2str(std::vector<double> d);
void lrtrim(std::string &s);
bool in_string(const std::string &x, std::string part);

void make_unique_names(std::vector<std::string> &s);
void make_valid_names(std::vector<std::string> &s);

void str_replace(std::string& str, const std::string& from, const std::string& to);
size_t str_replace_all(std::string& str, const std::string& from, const std::string& to);

