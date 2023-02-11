
#ifndef VECMATHSE_GUARD
#define VECMATHSE_GUARD

double median_se_rm(const std::vector<double>& v, size_t s, size_t e);
double median_se(const std::vector<double>& v, size_t s, size_t e);
double sum_se_rm(const std::vector<double>& v, size_t s, size_t e);
double sum_se(const std::vector<double>& v, size_t s, size_t e);
double sum2_se_rm(const std::vector<double>& v, size_t s, size_t e);
double sum2_se(const std::vector<double>& v, size_t s, size_t e);
double prod_se_rm(const std::vector<double>& v, size_t s, size_t e);
double prod_se(const std::vector<double>& v, size_t s, size_t e);
double mean_se_rm(const std::vector<double>& v, size_t s, size_t e);
double mean_se(const std::vector<double>& v, size_t s, size_t e);
double sd_se_rm(const std::vector<double>& v, size_t s, size_t e);
double sd_se(const std::vector<double>& v, size_t s, size_t e);
double sdpop_se_rm(const std::vector<double>& v, size_t s, size_t e);
double sdpop_se(const std::vector<double>& v, size_t s, size_t e);
double min_se_rm(const std::vector<double>& v, size_t s, size_t e);
double min_se(const std::vector<double>& v, size_t s, size_t e);
double max_se_rm(const std::vector<double>& v, size_t s, size_t e);
double max_se(const std::vector<double>& v, size_t s, size_t e);
double first_se_rm(std::vector<double>& v, size_t s, size_t e);
double first_se(std::vector<double>& v, size_t s, size_t e);
double which_se_rm(const std::vector<double>& v, size_t s, size_t e);
double which_se(const std::vector<double>& v, size_t s, size_t e);
double whichmin_se_rm(const std::vector<double>& v, size_t s, size_t e);
double whichmin_se(const std::vector<double>& v, size_t s, size_t e);
double whichmax_se_rm(const std::vector<double>& v, size_t s, size_t e);
double whichmax_se(const std::vector<double>& v, size_t s, size_t e);
double all_se_rm(const std::vector<double>& v, size_t s, size_t e);
double all_se(const std::vector<double>& v, size_t s, size_t e);
double any_se_rm(const std::vector<double>& v, size_t s, size_t e);
double any_se(const std::vector<double>& v, size_t s, size_t e);
std::vector<double> range_se_rm(std::vector<double>& v, size_t s, size_t e);
std::vector<double> range_se(std::vector<double>& v, size_t s, size_t e);
double modal_se_rm(std::vector<double>& v, size_t s, size_t e);
double modal_se(std::vector<double>& v, size_t s, size_t e);
std::vector<bool> isna_se(const std::vector<double>& v, size_t s, size_t e);
std::vector<bool> isnotna_se(const std::vector<double>& v, size_t s, size_t e);
void cumsum_se_rm(std::vector<double>& v, size_t s, size_t e);
void cumsum_se(std::vector<double>& v, size_t s, size_t e);
void cumprod_se_rm(std::vector<double>& v, size_t s, size_t e);
void cumprod_se(std::vector<double>& v, size_t s, size_t e);
void cummax_se_rm(std::vector<double>& v, size_t s, size_t e);
void cummax_se(std::vector<double>& v, size_t s, size_t e);
void cummin_se_rm(std::vector<double>& v, size_t s, size_t e);
void cummin_se(std::vector<double>& v, size_t s, size_t e);

bool haveseFun(std::string fun);
std::function<double(std::vector<double>&, double, double)> getseFun(std::string fun, bool narm);

#endif

