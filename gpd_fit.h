/*
 * gpd_fit.h
 *
 *  Created on: 13 Aug 2010
 *      Author: auton
 */

#ifndef GPD_FIT_H_
#define GPD_FIT_H_

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

double find_tail_approximation_p_value(const double clk_ratio, const vector<double> &clk_ratio_distribution_in);

#endif /* GPD_FIT_H_ */
