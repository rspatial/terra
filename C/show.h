void show(SpatRaster &r);
void show(SpatVector &v);

void showValues(SpatRaster d);
void showValues(std::vector<double> v, unsigned nr, unsigned nc);
void showValues(std::vector<std::vector<double>> v);
void showValues(std::vector<std::vector<std::vector<double>>> v);

template <typename T>
void showValues(std::vector<T> v) {
    for (size_t i=0; i<v.size(); i++) {
       if (std::isnan(v[i])) {
           std::cout << "NA" << " ";
       } else {
          std::cout << v[i] << " ";
       }
    }
    std::cout << "\n";
}



