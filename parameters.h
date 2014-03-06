/*
 * parameters.cpp
 *
 *  Created on: Nov 11, 2009
 *      Author: Adam Auton
 *      ($Revision: 148 $)
 */

// Class for reading in, checking and storing user parameters
#ifndef PARAMETERS_H_
#define PARAMETERS_H_

#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <string>
#include <vector>
#include <set>

#include "output_log.h"

using namespace std;

const string MY_LDHOT_VERSION="v0.4";

class parameters
{
public:
	string loc_filename;
	string seq_filename;
	string lk_filename;
	string res_filename;
	string output_prefix;
	double window_distance;
	double hotspot_distance;
	double start_pos;
	double end_pos;
	double pos_step;
	int N_sims;
	int seed;
	bool freq_cond;
	unsigned int freq_cond_model;
	int lk_SNP_window;

	parameters(int argc, char *argv[]);
	~parameters(){};

	void read_parameters();
	void print_help();
	void print_params();

private:
	void check_parameters();
	static void error(string err_msg, int code);

	vector<string> argv;

	string get_arg(unsigned int i);
};


#endif /* PARAMETERS_H_ */
