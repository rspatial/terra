
#ifndef VECMATHSE_GUARD
#define VECMATHSE_GUARD

inline double median_se_rm(const std::vector<double>& v, size_t s, size_t e);
inline double median_se(const std::vector<double>& v, size_t s, size_t e);
inline double sum_se_rm(const std::vector<double>& v, size_t s, size_t e);
inline double sum_se(const std::vector<double>& v, size_t s, size_t e);
inline double sum2_se_rm(const std::vector<double>& v, size_t s, size_t e);
inline double sum2_se(const std::vector<double>& v, size_t s, size_t e);
inline double prod_se_rm(const std::vector<double>& v, size_t s, size_t e);
inline double prod_se(const std::vector<double>& v, size_t s, size_t e);
inline double mean_se_rm(const std::vector<double>& v, size_t s, size_t e);
inline double mean_se(const std::vector<double>& v, size_t s, size_t e);
inline double sd_se_rm(const std::vector<double>& v, size_t s, size_t e);
inline double sd_se(const std::vector<double>& v, size_t s, size_t e);
inline double sdpop_se_rm(const std::vector<double>& v, size_t s, size_t e);
inline double sdpop_se(const std::vector<double>& v, size_t s, size_t e);
inline double min_se_rm(const std::vector<double>& v, size_t s, size_t e);
inline double min_se(const std::vector<double>& v, size_t s, size_t e);
inline double max_se_rm(const std::vector<double>& v, size_t s, size_t e);
inline double max_se(const std::vector<double>& v, size_t s, size_t e);
inline double first_se_rm(std::vector<double>& v, size_t s, size_t e);
inline double first_se(std::vector<double>& v, size_t s, size_t e);
inline double which_se_rm(const std::vector<double>& v, size_t s, size_t e);
inline double which_se(const std::vector<double>& v, size_t s, size_t e);
inline double whichmin_se_rm(const std::vector<double>& v, size_t s, size_t e);
inline double whichmin_se(const std::vector<double>& v, size_t s, size_t e);
inline double whichmax_se_rm(const std::vector<double>& v, size_t s, size_t e);
inline double whichmax_se(const std::vector<double>& v, size_t s, size_t e);
inline double all_se_rm(const std::vector<double>& v, size_t s, size_t e);
inline double all_se(const std::vector<double>& v, size_t s, size_t e);
inline double any_se_rm(const std::vector<double>& v, size_t s, size_t e);
inline double any_se(const std::vector<double>& v, size_t s, size_t e);
inline std::vector<double> range_se_rm(std::vector<double>& v, size_t s, size_t e);
inline std::vector<double> range_se(std::vector<double>& v, size_t s, size_t e);
inline double modal_se_rm(std::vector<double>& v, size_t s, size_t e);
inline double modal_se(std::vector<double>& v, size_t s, size_t e);
inline std::vector<bool> isna_se(const std::vector<double>& v, size_t s, size_t e);
inline std::vector<bool> isnotna_se(const std::vector<double>& v, size_t s, size_t e);
inline void cumsum_se_rm(std::vector<double>& v, size_t s, size_t e);
inline void cumsum_se(std::vector<double>& v, size_t s, size_t e);
inline void cumprod_se_rm(std::vector<double>& v, size_t s, size_t e);
inline void cumprod_se(std::vector<double>& v, size_t s, size_t e);
inline void cummax_se_rm(std::vector<double>& v, size_t s, size_t e);
inline void cummax_se(std::vector<double>& v, size_t s, size_t e);
inline void cummin_se_rm(std::vector<double>& v, size_t s, size_t e);
inline void cummin_se(std::vector<double>& v, size_t s, size_t e);

bool haveseFun(std::string fun);
std::function<double(std::vector<double>&, double, double)> getseFun(std::string fun, bool narm);

#endif

