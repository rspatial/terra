bool file_exists(const std::string& name);
void vector_minmax(std::vector<double> v, double &min, int &imin, double &max, int &imax);
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

template <class T> class NA {
public:
    static constexpr T value = std::is_floating_point<T>::value ? NAN : std::numeric_limits<T>::min();
};

template <typename T>
bool is_NA(const T v) {
    if (std::is_floating_point<T>::value) {
        return std::isnan(v);
    } else {
        bool b = v == NA<T>::value;
        return b;
    }
}

template <typename Iterator>
void minmax(Iterator start, Iterator end, double &min, double &max) {
    min = std::numeric_limits<double>::max();
    max = std::numeric_limits<double>::lowest();
    bool none = true;
	for (Iterator v = start; v !=end; ++v) {
		if (!std::isnan(*v)) {
			if (*v > max) {
				max = *v;
                none = false;
			}
			if (*v < min) {
				min = *v;
			}
		}
    }
    if (none) {
        min = NAN;
        max = NAN;
    }
}

