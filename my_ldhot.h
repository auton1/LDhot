/*
 * my_ldhot.h
 *
 *  Created on: Apr 16, 2010
 *      Author: auton
 */

#ifndef MY_LDHOT_H_
#define MY_LDHOT_H_

#include <ctime>
#include <iostream>
#include <limits>
#include <string>
#include <vector>

#include "output_log.h"
#include "parameters.h"
#include "tools.h"
#include "sim.h"
#include "Ran.h"
#include "gpd_fit.h"

//#include <omp.h>

using namespace std;

ofstream LOG;

const string LDHOT_VERSION="v0.3";

void output_result(ofstream &out, double lh_hotspot_pos, double rh_hotspot_pos,
		double mle_rate_hotspot,
		double lh_window_pos, double rh_window_pos,
		double mle_rate_background, double mle_rate_constant,
		unsigned int N_sims_used, double p_value, double p_value_approx);

#endif /* MY_LDHOT_H_ */
