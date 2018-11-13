/* 
using namespace std;
using namespace Rcpp;

// xy: num[1:n, 1:2]
// x:  num[1:2, 1:n]
// y:  num[1:2, 1:n]
// v:  num[1:n, 1:4]
//     columns are: bottom-left, top-left, top-right, bottom-right

// [[Rcpp::export(name = ".doBilinear")]]
NumericVector doBilinear(NumericMatrix xy, NumericMatrix x, NumericMatrix y, NumericMatrix v) {
	size_t len = v.nrow();
	NumericVector result(len);

    for (size_t i = 0; i < len; i++) {
		double left = x(0,i);
		double right = x(1,i);
		double top = y(1,i);
		double bottom = y(0,i);
		
		double horiz = xy(i,0);
		double vert = xy(i,1);
		
		double denom = (right - left) * (top - bottom);
		
		double bottomLeftValue = v(i,0) / denom;
		double topLeftValue = v(i,1) / denom;
		double topRightValue = v(i,3) / denom;
		double bottomRightValue = v(i,2) / denom;
		result[i] = bottomLeftValue*(right-horiz)*(top-vert) + bottomRightValue*(horiz-left)*(top-vert) +
			topLeftValue*(right-horiz)*(vert-bottom) + topRightValue*(horiz-left)*(vert-bottom);
	}

	return result;
}


 */