Rcpp::List getBlockSizeR(SpatRaster* r, unsigned n);
Rcpp::List getAttributes(SpatVector* v);
Rcpp::List getDataFrame(SpatDataFrame* v);

//bool setAttributes(SpatVector* v, Rcpp::List x, std::vector<std::string> names, std::vector<std::string> types);

Rcpp::NumericVector getGeometry(SpatVector* v);

SpatRaster rcppReclassify(SpatRaster* x, Rcpp::NumericMatrix rcl, unsigned right, bool lowest, bool othersNA, SpatOptions &opt);
//Rcpp::NumericMatrix rcppAdjacent(SpatRaster* x, std::vector<double> cells, std::string directions, bool include);
