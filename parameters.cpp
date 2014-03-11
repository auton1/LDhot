/*
 * parameters.cpp
 *
 *  Created on: Nov 11, 2009
 *      Author: Adam Auton
 */

// Class for reading in, checking and storing user parameters

#include "parameters.h"

parameters::parameters(int argc, char *argv[])
{
	string tmp;
	for (int i=0; i<argc; i++)
	{
		tmp = argv[i];
		this->argv.push_back(tmp);
	}

	seq_filename=""; loc_filename=""; lk_filename=""; res_filename="";
	output_prefix="out";
	window_distance=50.0;	// distance either side of hotspot (kb)
	hotspot_distance=1.5;	// distance either side of centre of hotspot (kb)
	start_pos=0;
	end_pos=numeric_limits<double>::max();
	pos_step = 1;			// Step size between tests (kb)
	N_sims = 100;
	seed = 10;
	freq_cond = true;
	freq_cond_model = 0;
	lk_SNP_window = 50;
	n_threads = -1;
}

void parameters::read_parameters()
{
	print_help();
	unsigned int i=1;
	string in_str;
	seed = (long)time(NULL);
	while (i<argv.size())
	{
		in_str = argv[i];
		if (in_str ==  "--seq") { seq_filename = get_arg(i+1); i++; }
		else if (in_str ==  "--loc") { loc_filename = get_arg(i+1); i++; }
		else if (in_str ==  "--lk") { lk_filename = get_arg(i+1); i++; }
		else if (in_str ==  "--res") { res_filename = get_arg(i+1); i++; }
		else if (in_str ==  "--out") { output_prefix = get_arg(i+1); i++; }
		else if (in_str ==  "--windist") { window_distance = atof(get_arg(i+1).c_str()); i++; }
		else if (in_str ==  "--hotdist") { hotspot_distance = atof(get_arg(i+1).c_str()); i++; }
		else if (in_str ==  "--startpos") { start_pos = atof(get_arg(i+1).c_str()); i++; }
		else if (in_str ==  "--endpos") { end_pos = atof(get_arg(i+1).c_str()); i++; }
		else if (in_str ==  "--step") { pos_step = atof(get_arg(i+1).c_str()); i++; }
		else if (in_str ==  "--nsim") { N_sims = atoi(get_arg(i+1).c_str()); i++; }
		else if (in_str ==  "--seed") { seed = atoi(get_arg(i+1).c_str()); i++; }
		else if (in_str ==  "--nofreqcond") { freq_cond = false; }
		else if (in_str ==  "--freq-cond-model") { freq_cond = true; freq_cond_model = atoi(get_arg(i+1).c_str()); i++; }
		else if (in_str ==  "--lk-SNP-window") { lk_SNP_window = (unsigned)atoi(get_arg(i+1).c_str()); i++; }
		else if (in_str ==  "--nthreads") { n_threads = atoi(get_arg(i+1).c_str()); i++; }
		else
			error("Unknown option: " + string(in_str), 0);
		i++;
	}
}

string parameters::get_arg(unsigned int i)
{
	if (i>=argv.size())
		error("Requested Missing Argument",76);
	return argv[i];
}

void parameters::print_params()
{
	parameters defaults(0, 0);

	printLOG("Parameters as interpreted:\n");
	if (seq_filename != defaults.seq_filename) printLOG("\t--seq " + seq_filename + "\n");
	if (loc_filename != defaults.loc_filename) printLOG("\t--loc " + loc_filename + "\n");
	if (lk_filename != defaults.lk_filename) printLOG("\t--lk " + lk_filename + "\n");
	if (res_filename != defaults.res_filename) printLOG("\t--res " + res_filename + "\n");
	if (output_prefix != defaults.output_prefix) printLOG("\t--out " + output_prefix + "\n");
	if (start_pos != defaults.start_pos) printLOG("\t--startpos " + dbl2str_fixed(start_pos, 3) + "\n");
	if (end_pos != defaults.end_pos) printLOG("\t--endpos " + dbl2str_fixed(end_pos, 3) + "\n");
	if (pos_step != defaults.pos_step) printLOG("\t--step " + dbl2str_fixed(pos_step, 3) + "\n");
	if (window_distance != defaults.window_distance) printLOG("\t--windist " + dbl2str_fixed(window_distance, 3) + "\n");
	if (hotspot_distance != defaults.hotspot_distance) printLOG("\t--hotdist " + dbl2str_fixed(hotspot_distance, 3) + "\n");
	if (N_sims != defaults.N_sims) printLOG("\t--nsim " + int2str(N_sims) + "\n");
	if (freq_cond != defaults.freq_cond) printLOG("\t--nofreqcond\n");
	if (freq_cond_model != defaults.freq_cond_model) printLOG("\t--freq-cond-model " + int2str(freq_cond_model) + "\n");
	if (lk_SNP_window != defaults.lk_SNP_window) printLOG("\t--lk-SNP-window " + int2str(lk_SNP_window) + "\n");
	if (n_threads != defaults.n_threads) printLOG("\t--n_threads " + int2str(n_threads) + "\n");
	printLOG("\t--seed " + longint2str(seed) + "\n");
	printLOG("\n");
}

void parameters::print_help()
{
	unsigned int i;
	string in_str;

	if (argv.size() <= 1)
	{	// If there are no user parameters, display help.
		argv.push_back("--?");
		print_help();
	}

	for(i = 0; i < argv.size(); i++)
	{
		in_str = argv[i];
		if ((in_str == "-h") || (in_str == "-?") || (in_str == "-help") || (in_str == "--?") || (in_str == "--help") || (in_str == "--h"))
		{
			cout << endl << "LDhot (" << LDHOT_VERSION << ")" << endl;
			cout << "\u00A9 Adam Auton 2014" << endl << endl;

			cout << "Required Parameters: " << endl;
			cout << "--seq <filename>" << endl;
			cout << "--loc <filename>" << endl;
			cout << "--lk <filename> " << endl;
			cout << "--res <filename> " << endl;
			cout << endl;
			cout << "Important Parameters: " << endl;
			cout << "--out <prefix> " << endl;
			cout << "--nsim <int>" << endl;
			cout << endl;
			cout << "Other Parameters: " << endl;
			cout << "--startpos <double>" << endl;
			cout << "--endpos <double> " << endl;
			cout << "--step <double>" << endl;
			cout << "--windist <double>" << endl;
			cout << "--hotdist <double>" << endl;
			cout << "--seed <int>" << endl;
			cout << "--nofreqcond" << endl;
			cout << "--freq-cond-model <int>" << endl;
			cout << "--lk-SNP-window <int>" << endl;
#ifdef _OPENMP
			cout << "--n_threads <int> " << endl;
#endif
			cout << endl;
			exit(0);
		}
	}
}

void parameters::check_parameters()
{
	if (seq_filename == "") error("SEQ required.", 0);
	if (end_pos < start_pos) error("End position must be greater than Start position.", 1);
	if (hotspot_distance >= window_distance) error("Hotspot Size must be smaller Window Size.", 2);
	if (pos_step < 0.001) error("Step must be greater than 0.001.", 3);
	if (freq_cond_model > 1) error("Unknown Frequency Conditioning model number", 4);
}

void parameters::error(string err_msg, int code)
{
	printLOG("\n\nError: " + err_msg + "\n\n");
	exit(code);
}
