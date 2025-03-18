#include <map>
#include <vector>
#include <cmath>
#include <algorithm>


std::map<double, unsigned long long> table(std::vector<double> &v) {
	std::map<double, unsigned long long> count;
	for_each( v.begin(), v.end(), [&count]( double val ){
			if(!std::isnan(val)) count[val]++;
		}
	);
	return count;
}


std::map<double, unsigned long long int> combine_tables(std::map<double, unsigned long long int> &x, std::map<double, unsigned long long int> &y) {
	for(auto p : y) {
		x[p.first] += p.second;
	}
	return(x);
}


std::vector<double> table2vector(std::map<double, unsigned long long int> &x) {
	std::vector<std::vector<double>> out(2);
	for( auto p : x ) {
		out[0].push_back(p.first);
		out[1].push_back(p.second);
	}
	out[0].insert(out[0].end(), out[1].begin(), out[1].end());
	return out[0];
}

std::vector<std::vector<double>> table2vector2(std::map<double, unsigned long long int> &x) {
	std::vector<std::vector<double>> out(2);
	for( auto p : x ) {
		out[0].push_back(p.first);
		out[1].push_back(p.second);
	}
	return out;
}
