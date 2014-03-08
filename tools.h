/*
 * tools.h
 *
 *  Created on: Apr 16, 2010
 *      Author: auton
 */

#ifndef TOOLS_H_
#define TOOLS_H_

#include <algorithm>
#include <cassert>
#include <cmath>
#include <limits>
#include <map>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>

#include "output_log.h"

using namespace std;

typedef vector<double> vec1D;
typedef map<unsigned int, vec1D> vec2D;
typedef map<unsigned int, vec2D> vec3D;
typedef map<unsigned int, vec3D> vec4D;
typedef map<unsigned int, vec4D> vec5D;

/*
typedef vector<double> vec1D;
typedef vector<vector<double> >  vec2D;
typedef vector<vector<vector<double> > > vec3D;
typedef vector<vector<vector<vector<double> > > > vec4D;
typedef vector<vector<vector<vector<vector<double> > > > > vec5D;
*/

class data_struct
{
public:
	data_struct()
	{
		is_haplotypes = false;
		N_seq = 0;
		N_sites = 0;
	};
	~data_struct() {};

	vector<double> locs;
	vector< vector< unsigned char > > seqs;
	vector<int> freq_counts;
	vector<string> headers;
	bool is_haplotypes;
	unsigned int N_seq;
	unsigned int N_sites;

	data_struct get_subset(unsigned int start_idx, unsigned int end_idx)
	{
		assert (end_idx < locs.size());
		assert (start_idx < end_idx);

		data_struct out;
		out.is_haplotypes = is_haplotypes;
		out.N_seq = N_seq;
		out.N_sites = end_idx - start_idx;
		out.locs.resize(end_idx - start_idx);
		copy(locs.begin()+start_idx, locs.begin()+end_idx, out.locs.begin());
		out.seqs.resize(N_seq, vector< unsigned char >(end_idx - start_idx));
		for (unsigned int ui=0; ui<N_seq; ui++)
			copy(seqs[ui].begin()+start_idx, seqs[ui].begin()+end_idx, out.seqs[ui].begin());

		out.freq_counts.resize(end_idx - start_idx);
		copy(freq_counts.begin()+start_idx, freq_counts.begin()+end_idx, out.freq_counts.begin());

		//out.headers = headers;
		return out;
	}

	void calc_freq_counts()
	{
		freq_counts.resize(N_sites);
		for (unsigned int ui=0; ui<N_sites; ui++)
		{
			unsigned int count = 0;
			for (unsigned int uj=0; uj<N_seq; uj++)
			{
				count += seqs[uj][ui];
			}
			freq_counts[ui] = count;
		}
	}

	void remove_fixed_sites()
	{
		vector<double> new_locs;
		vector<int> new_freq_counts;
		vector< vector< unsigned char > > new_seqs(N_seq);
		for (unsigned int uj=0; uj<N_seq; uj++)
			new_seqs[uj].reserve(N_sites);
		new_locs.reserve(locs.size());

		calc_freq_counts();
		for (unsigned int ui=0; ui<N_sites; ui++)
		{
			if ((freq_counts[ui] > 0) && (freq_counts[ui] < (signed)N_seq))
			{
				new_locs.push_back(locs[ui]);
				new_freq_counts.push_back(freq_counts[ui]);
				for (unsigned int uj=0; uj<N_seq; uj++)
					new_seqs[uj].push_back(seqs[uj][ui]);
			}
		}

		locs = new_locs;
		seqs = new_seqs;
		N_sites = locs.size();
		freq_counts = new_freq_counts;
	}

	void remove_singleton_sites()
	{
		vector<double> new_locs;
		vector<int> new_freq_counts;
		vector< vector< unsigned char > > new_seqs(N_seq);
		for (unsigned int uj=0; uj<N_seq; uj++)
			new_seqs[uj].reserve(N_sites);
		new_locs.reserve(locs.size());

		calc_freq_counts();
		for (unsigned int ui=0; ui<N_sites; ui++)
		{
			cout << freq_counts[ui] << endl;
			if ((freq_counts[ui] > 1) && (freq_counts[ui] < ((signed)N_seq-1)))
			{
				new_locs.push_back(locs[ui]);
				new_freq_counts.push_back(freq_counts[ui]);
				for (unsigned int uj=0; uj<N_seq; uj++)
					new_seqs[uj].push_back(seqs[uj][ui]);
			}
			else
			{
				cout << "here" << endl;
			}
		}

		locs = new_locs;
		seqs = new_seqs;
		N_sites = locs.size();
		freq_counts = new_freq_counts;
	}

	void output_to_file(const string &output_filename_prefix)
	{
		string locs_file = output_filename_prefix + ".locs";
		string seqs_file = output_filename_prefix + ".seqs";

		ofstream out(locs_file.c_str());
		if (!out.good())
			error("Could not open locs file for output: " + locs_file);
		out << N_sites << "\t" << locs[N_sites-1] - locs[0] << "\tL" << endl;
		for (unsigned int ui=0; ui<N_sites; ui++)
			out << locs[ui] << endl;
		out.close();

		out.open(seqs_file.c_str());
		if (!out.good())
			error("Could not open locs file for output: " + locs_file);
		out << N_seq << "\t" << N_sites << "\t1" << endl;
		for (unsigned int ui=0; ui<N_seq; ui++)
		{
			out << ">Seq-" << ui << endl;
			for (unsigned int uj=0; uj<N_sites; uj++)
			{
				out << (int)seqs[ui][uj];
				if ((uj+1 % 4095) == 0)
					out << endl;
			}
			out << endl;
		}
		out.close();
	}
};

class lk_table_type
{
public:
	lk_table_type()
	{
		N_seq = 0; N_rho = 0; lk_rho_max = 0;
	};
	~lk_table_type() {};
	vec5D lk_table;	// N00, N01, N10, N11, rho_idx
	vector<double> lk_rho_points;
	unsigned int N_seq;
	unsigned int N_rho;
	double lk_rho_max;
};

void read_locs_file(const string &loc_filename, data_struct &data);
void read_seqs_file(const string &seq_filename, data_struct &data);
void read_lk_file(const string &lk_filename, lk_table_type &lk);
void read_res_file(const string &res_filename, vector<double> &rmap, vector<double> &rmap_pos);

void get_pair_likelihoods(const data_struct &data, const lk_table_type &lk,
							const int start_idx, const int end_idx,
							vector<vector<vector<double > > > &pair_lks,
							const int lk_SNP_window);
double calc_composite_likelihood_using_pair_likelihoods(const lk_table_type &lk, const vector<double> &rmap,
														const int start_idx, const int end_idx,
														const vector<vector<vector<double> > > &pair_lks,
														const int lk_SNP_window);

double estimate_MLE_constant_rate_across_region(const data_struct &data, const lk_table_type &lk,
													const int lk_SNP_window, double &mle_rate);
double estimate_MLE_hotspot_rate_across_region(const data_struct &data, const lk_table_type &lk,
												const int lk_SNP_window, const double lh_hotspot_pos,
												const double rh_hotspot_pos,
												double &bg_mle_rate, double &hot_mle_rate);

void find_maxima_locations_in_rates(const vector<double> &rmap_pos, const vector<double> &rmap_rate,
									vector<double> &maxima_pos_LHS, vector<double> &maxima_pos_RHS);

double get_median_bgrate_in_region(const vector<double> &rmap_pos, const vector<double> &rmap_rate,
									double lh_window_pos, double rh_window_pos, double lh_hotspot_pos, double rh_hotspot_pos);

double get_avg_rate_across_region(const vector<double> &rmap_pos, const vector<double> &rmap,
									double lh_window_pos, double rh_window_pos);

#endif /* TOOLS_H_ */
