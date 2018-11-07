bool is_equal(double a, double b, double error_factor=1.0);
bool is_equal_range(double x, double y, double range, double tolerance);
void vector_minmax(std::vector<double> v, double &min, int &imin, double &max, int &imax);
double roundn(double x, int n);

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

