// Copyright (c) 2018-2023  Robert J. Hijmans
//
// This file is part of the "spat" library.
//
// spat is free software: you can redistribute it and/or modify it
// under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 2 of the License, or
// (at your option) any later version.
//
// spat is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with spat. If not, see <http://www.gnu.org/licenses/>.


#include "spatRaster.h"
#include "vecmath.h"
#include "sort.h"


inline double rarea(const double &Ax, const double &Ay, const double &Bx, const double &By, const double &Cx, const double &Cy) {
   return std::abs( (Bx*Ay - Ax*By) + (Cx*By - Bx*Cy) + (Ax*Cy - Cx*Ay) ) / 2;
}


void sortvecs(std::vector<double> &X, std::vector<double> &Y, std::vector<double> &Z) {  
	std::vector<std::size_t> p = sort_order_d(X);
	permute(X, p);
	permute(Y, p);
	permute(Z, p);
  
	p = sort_order_d(Y);
	permute(X, p);
	permute(Y, p);
	permute(Z, p);
}


std::vector<std::vector<double>> SpatRaster::win_rect(std::vector<double> x, std::vector<double> y, std::vector<double> z, std::vector<double> win, SpatOptions &opt) {

	sortvecs(x, y, z);

	win[0] = std::abs(win[0]);
	win[1] = std::abs(win[1]);
    const double h = win[0] / 2; 
    const double w = win[1] / 2;

	// multiply for floating point imprecision
	const double rar = win[0] * win[1] * 1.00000001;
	double angle = std::fmod(win[2], 360.0);
	if (angle < 0) angle += 360.0;
    const bool rotated = angle != 0.0;

	double cphi=0, sphi=0, bigw=0, bigh=0, offh;
	double wcphi=0, hcphi=0, wsphi=0, hsphi=0;
	std::vector<double> ox(4);
	std::vector<double> oy(4);

	if (rotated) {
		angle = angle * M_PI / 180.0;
		cphi = cos(angle);
		sphi = sin(angle);
		wcphi= cphi * w;
		hcphi= cphi * h;
		wsphi= sphi * w;
		hsphi= sphi * h;

		ox[0] = -wcphi - hsphi;
		oy[0] = -wsphi + hcphi;
		ox[1] =  wcphi - hsphi;
		oy[1] =  wsphi + hcphi;
		ox[2] =  wcphi + hsphi;
		oy[2] =  wsphi - hcphi;
		ox[3] = -wcphi + hsphi;
		oy[3] = -wsphi - hcphi;

		bigw = (vmax(ox, false) - vmin(ox, false))/2;
		bigh = (vmax(oy, false) - vmin(oy, false))/2;
		offh = bigh * 1.00000001;
	} else {
		offh = h * 1.00000001;
	}

	const size_t nc = ncol();
	const size_t nr = nrow();

	size_t np = x.size() * 2;
	std::vector<std::vector<double>> out(2);
	out[0].reserve(np);
	out[1].reserve(np);

	size_t minpt = win[3] < 2 ? 1 : win[3];
	std::vector<double> rx(4);
	std::vector<double> ry(4);

	std::vector<int_64> cols(nc);
	std::iota(cols.begin(), cols.end(), 0);
	std::vector<double> xc = xFromCol(cols);

	if (minpt < 2) {

		if (rotated) {
			for (size_t r=0; r<nr; r++) {
				double yrow = yFromRow(r);
				double ytop = yrow + offh;
				double ybot = yrow - offh;
				double rnc = r * nc;
				ry[0] = yrow + oy[0];
				ry[1] = yrow + oy[1];
				ry[2] = yrow + oy[2];
				ry[3] = yrow + oy[3];
				np = y.size();
				for (long i=(np-1); i>=0; i--) {
					if (y[i] > ytop) {
						y.pop_back(); // above current row
					} else if (y[i] >= ybot) {
						for (long j=(nc-1); j>=0; j--) {
							double dist = x[i] - xc[j];
							if (dist > bigw) {
								break;
							} else if (std::abs(dist) <= bigw) {
								rx[0] = xc[j] + ox[0];
								rx[1] = xc[j] + ox[1];
								rx[2] = xc[j] + ox[2];
								rx[3] = xc[j] + ox[3];
								// triangles apd, dpc, cpb, bpa
								double area = rarea(rx[0], ry[0], x[i], y[i], rx[3], ry[3]);
								area += rarea(rx[3], ry[3], x[i], y[i], rx[2], ry[2]);
								area += rarea(rx[2], ry[2], x[i], y[i], rx[1], ry[1]);
								area += rarea(rx[1], ry[1], x[i], y[i], rx[0], ry[0]);
								if (area < rar) {
									out[0].push_back(rnc+j);
									out[1].push_back(z[i]);
								}
							} 
						}
					} else {
						break;
					}
				}
			}
		} else {
			for (size_t r=0; r<nr; r++) {
				double yrow = yFromRow(r);
				double ytop = yrow + offh;
				double ybot = yrow - offh;
				double rnc = r*nc;
				np = y.size();
				for (long i=(np-1); i>=0; i--) {
					if (y[i] > ytop) {
						y.pop_back(); // above current row
					} else if (y[i] >= ybot) {
						for (long j=(nc-1); j>=0; j--) {
							double dist = x[i] - xc[j];
							if (dist > w) {
								break;
							} else if (std::abs(dist) <= w) {
								out[0].push_back(rnc+j);
								out[1].push_back(z[i]);
							}
						}
					} else {
						break;
					}
				} 
			}
		}
	} else {

		if (rotated) {
			for (size_t r=0; r<nr; r++) {
				double yrow = yFromRow(r);
				double ytop = yrow + offh;
				double ybot = yrow - offh;
				double rnc = r*nc;
				ry[0] = yrow + oy[0];
				ry[1] = yrow + oy[1];
				ry[2] = yrow + oy[2];
				ry[3] = yrow + oy[3];
				np = y.size();
				std::vector<double> tmp0, tmp1;
				for (long i=(np-1); i>=0; i--) {
					if (y[i] > ytop) {
						y.pop_back(); // above current row
					} else if (y[i] >= ybot) {
						bool found = false;
						size_t minlim = 0;
						for (long j=(nc-1); j>=0; j--) {
							double dist = x[i] - xc[j];
							if (dist > bigw) {
								break;
							} else if (std::abs(dist) <= bigw) {
								rx[0] = xc[j] + ox[0];
								rx[1] = xc[j] + ox[1];
								rx[2] = xc[j] + ox[2];
								rx[3] = xc[j] + ox[3];
								// triangles apd, dpc, cpb, bpa
								double area = rarea(rx[0], ry[0], x[i], y[i], rx[3], ry[3]);
								area += rarea(rx[3], ry[3], x[i], y[i], rx[2], ry[2]);
								area += rarea(rx[2], ry[2], x[i], y[i], rx[1], ry[1]);
								area += rarea(rx[1], ry[1], x[i], y[i], rx[0], ry[0]);
								if (area < rar) {
									tmp0.push_back(rnc+j);
									tmp1.push_back(z[i]);
									found = true;
									minlim++;
								}
							} 
						}
						if (found) {
							if (minlim >= minpt) {
								out[0].insert(out[0].end(), tmp0.begin(), tmp0.end());
								out[1].insert(out[1].end(), tmp1.begin(), tmp1.end());
							}
							tmp0.resize(0);
							tmp1.resize(0);
							tmp0.reserve(10);
							tmp1.reserve(10);
						}
					} else {
						break;
					}
				}
			}
		} else {
			for (size_t r=0; r<nr; r++) {
				double yrow = yFromRow(r);
				double ytop = yrow + offh;
				double ybot = yrow - offh;
				double rnc = r*nc;
				np = y.size();
				std::vector<double> tmp0, tmp1;
				for (long i=(np-1); i>=0; i--) {
					if (y[i] > ytop) {
						y.pop_back(); // above current row
					} else if (y[i] >= ybot) {
						bool found = false;
						size_t minlim = 0;
						for (long j=(nc-1); j>=0; j--) {
							double dist = x[i] - xc[j];
							if (dist > w) {
								break;
							} else if (std::abs(dist) <= w) {
								tmp0.push_back(rnc+j);
								tmp1.push_back(z[i]);
								found = true;
								minlim++;
							}
						}
						if (found) {
							if (minlim >= minpt) {
								out[0].insert(out[0].end(), tmp0.begin(), tmp0.end());
								out[1].insert(out[1].end(), tmp1.begin(), tmp1.end());
							}
							tmp0.resize(0);
							tmp1.resize(0);
							tmp0.reserve(10);
							tmp1.reserve(10);
						}
					} else {
						break;
					}
				} 
			}
		}
	}
    return out;
}




std::vector<std::vector<double>> SpatRaster::win_circle(std::vector<double> x, std::vector<double> y, std::vector<double> z, std::vector<double> win, SpatOptions &opt) {

// the basic approach is from GDALGRID

	sortvecs(x, y, z);

    const double radius1 = win[0] * win[0]; 
    const double radius2 = win[1] * win[1];
	const double R12 = radius1 * radius2;

	double angle = std::fmod(win[2], 360.0);
	if (angle < 0) angle += 360.0;
    const bool rotated = angle != 0.0;
    angle = angle * M_PI / 180.0;

    //  coefficients for coordinate system rotation.
    const double cf1 = rotated ? cos(angle) : 0.0;
    const double cf2 = rotated ? sin(angle) : 0.0;

	size_t minpt = win[3] < 2 ? 1 : win[3];

	// for now assuming circles
    const double h = std::max(win[0], win[1]); 
//    const double w = h;

	const size_t nc = ncol();
	const size_t nr = nrow();

	size_t np = x.size() * 2;
	std::vector<std::vector<double>> out(2);
	out[0].reserve(np);
	out[1].reserve(np);

	std::vector<double> rx(4);
	std::vector<double> ry(4);

	std::vector<int_64> cols(nc);
	std::iota(cols.begin(), cols.end(), 0);
	std::vector<double> xc = xFromCol(cols);


	if (minpt < 2) {
		for (size_t r=0; r<nr; r++) {
			double yrow = yFromRow(r);
			double ytop = yrow + h;
			double ybot = yrow - h;
			double rnc = r*nc;
			np = y.size();
			for (long i=(np-1); i>=0; i--) {
				if (y[i] > ytop) {
					y.pop_back(); // above current row
				} else if (y[i] >= ybot) {
					bool found = false;
					for (long j=(nc-1); j>=0; j--) {
						double RX = x[i] - xc[j];
						double RY = y[i] - yrow;
						if (rotated) {
							RX = RX * cf1 + RY * cf2;
							RY = RY * cf1 - RX * cf2;
						}
						if ((radius2 * RX * RX + radius1 * RY * RY) <= R12) {
							out[0].push_back(rnc+j);
							out[1].push_back(z[i]);
							found = true;
						} else {
							if (found) break;
						}
					}
				} else {
					break;
				}
			}
		}
	} else {
		for (size_t r=0; r<nr; r++) {
			double yrow = yFromRow(r);
			double ytop = yrow + h;
			double ybot = yrow - h;
			double rnc = r*nc;
			np = y.size();
			std::vector<double> tmp0, tmp1;
			for (long i=(np-1); i>=0; i--) {
				if (y[i] > ytop) {
					y.pop_back(); // above current row
				} else if (y[i] >= ybot) {
					bool found = false;
					size_t minlim = 0;
					for (long j=(nc-1); j>=0; j--) {
						double RX = x[i] - xc[j];
						double RY = y[i] - yrow;
						if (rotated) {
							RX = RX * cf1 + RY * cf2;
							RY = RY * cf1 - RX * cf2;
						}
						if ((radius2 * RX * RX + radius1 * RY * RY) <= R12) {
							tmp0.push_back(rnc+j);
							tmp1.push_back(z[i]);
							found = true;
							minlim++;
						}
					}
					if (found) {
						if (minlim >= minpt) {
							out[0].insert(out[0].end(), tmp0.begin(), tmp0.end());
							out[1].insert(out[1].end(), tmp1.begin(), tmp1.end());
						}
						tmp0.resize(0);
						tmp1.resize(0);
						tmp0.reserve(10);
						tmp1.reserve(10);
					}
				} else {
					break;
				}
			}
		}
	}
    return out;
}

