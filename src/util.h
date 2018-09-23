double roundn(double x, int n);
std::string concatenate(std::vector<string> v, std::string delim);
bool is_equal(double a, double b, double error_factor=1.0);
bool is_equal_range(double x, double y, double range, double tolerance);
void lowercase(std::string &s); 
bool is_in_set(string s, std::vector<string> ss);
std::string is_in_set_default(string s, std::vector<string> ss, string defvalue, bool tolower);
std::vector<std::string> strsplit(std::string s, std::string delimiter);
std::vector<double> str2dbl(std::vector<string> s);
std::vector<string> dbl2str(std::vector<double> d);
string getFileExt(const string& s);
string setFileExt(const string& s, const string& ext);
void lrtrim(std::string &s);

