bool file_exists(const std::string& name);
void vector_minmax(std::vector<double> v, double &min, int &imin, double &max, int &imax);
double roundn(double x, int n);
std::string concatenate(std::vector<std::string> v, std::string delim);
bool is_equal(double a, double b, double error_factor=1.0);
bool is_equal_range(double x, double y, double range, double tolerance);
void lowercase(std::string &s);
bool is_in_set(std::string s, std::vector<std::string> ss);
std::string is_in_set_default(std::string s, std::vector<std::string> ss, std::string defvalue, bool tolower);
std::vector<std::string> strsplit(std::string s, std::string delimiter);
std::vector<double> str2dbl(std::vector<std::string> s);
std::vector<std::string> dbl2str(std::vector<double> d);
std::string getFileExt(const std::string& s);
std::string setFileExt(const std::string& s, const std::string& ext);
std::string basename(std::string filename);
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

template <typename T>
void setNAN(std::vector<T> &v, double naflag) {
	if (!std::isnan(naflag)) {
		T flag = naflag;
		T navalue = NA<T>::value;
		std::replace(v.begin(), v.end(), flag, navalue);
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

