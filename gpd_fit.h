/*
 * gpd_fit.h
 *
 *  Created on: 13 Aug 2010
 *      Author: auton
 */

#ifndef GPD_FIT_H_
#define GPD_FIT_H_

/*
 * gpd_fit.cpp
 *
 *  Created on: 13 Aug 2010
 *      Author: auton
 */

/*
 * gpd_fit.cpp
 *
 *  Created on: 13 Aug 2010
 *      Author: auton
 */

#include <algorithm>
#include <cmath>
#include <iostream>
#include <fstream>
#include <limits>
#include <numeric>
#include <vector>

using namespace std;

double gpd_llk(const vector<double> &data, double k, double a);
void gpd_method_of_moments(const vector<double> &data, double &k, double &a);
void gpd_ml(const vector<double> &data, double &k, double &a);
bool test_gof(const vector<double> &data, double k, double a);
double gp_cdf(double x, double k, double a);

#endif /* GPD_FIT_H_ */
