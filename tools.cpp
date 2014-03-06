/*
 * tools.cpp
 *
 *  Created on: Apr 16, 2010
 *      Author: auton
 */

#include "tools.h"

// Return the vector of rho likelihoods for each pair in a dataset
void get_pair_likelihoods(const data_struct &data, const lk_table_type &lk,
		const int start_idx, const int end_idx, vector<vector<vector<double > > > &pair_lks,
		const int lk_SNP_window)
{
	const unsigned int N_seq = data.seqs.size();
	int N_sites = end_idx - start_idx;
	pair_lks.resize(N_sites);
	for (int i=0; i<N_sites; i++)
		pair_lks[i].resize(N_sites);

	for (int i=start_idx; i<end_idx; i++)
	{
		int j_min = max(i-lk_SNP_window, 0);
		int j_max = min(i+lk_SNP_window+1, end_idx);
		for (int j=j_min; j<j_max; j++)
		{
			if (i==j)
				continue;
			unsigned int N[2][2];
			N[0][0] = 0; N[0][1] = 0; N[1][0] = 0; N[1][1] = 0;
			for (unsigned int s=0; s<N_seq; ++s)
			{	// This is the major performance bottleneck...
				N[ data.seqs[s][i] ][ data.seqs[s][j] ]++;
			}
			unsigned int a, b, c, d;
			a=N[0][0],b=N[1][0],c=N[0][1],d=N[1][1];

			// pt and type (haploid)
			// 00: 0 = a
			// 10: 1 = b
			// 01: 2 = c
			// 11: 3 = d

			unsigned int N_types = (a>0)+(b>0)+(c>0)+(d>0);
			if (N_types == 1)
				continue;

			if (N_types == 2)
			{
				a = max(max(N[0][0],N[1][0]),max(N[0][1],N[1][1]));
				b=0; c=0; d=a;
				if ((N[0][0] > 0) && (N[0][0] != a)) d=N[0][0];
				if ((N[0][1] > 0) && (N[0][1] != a)) d=N[0][1];
				if ((N[1][0] > 0) && (N[1][0] != a)) d=N[1][0];
				if ((N[1][1] > 0) && (N[1][1] != a)) d=N[1][1];
			}

			int fl=0;
			if (2*(b+d) > data.N_seq) fl += 2; //if (2*(pt[1]+pt[3]) > nseq) fl += 2;
			if (2*(c+d) > data.N_seq) fl += 1; //if (2*(pt[2]+pt[3]) > nseq) fl += 1;
			switch(fl)
			{
			case 0 :
				break;
			case 1 :
				swap(a,c); //pswap(pt,0,2);
				swap(b,d); //pswap(pt,1,3);
				break;
			case 2 :
				swap(a,b); //pswap(pt,0,1);
				swap(c,d); //pswap(pt,2,3);
				break;
			case 3 :
				swap(a,d); //pswap(pt,0,3);
				swap(b,c); //pswap(pt,1,2);
				break;
			default :
				printf("\n\nError in sorting\n\n");
				exit(1);
			}

			if (c>b)
			{
				swap(b,c);
			}

			//pair_lks[i][j] = lk.lk_table[a][b][c][d];
			pair_lks[i][j] = ((((lk.lk_table.find(a))->second.find(b))->second.find(c))->second.find(d))->second;
		}
	}
}


// Given a set of pair likelihoods, calculate the likelihood for a given recombination map.
double calc_composite_likelihood_using_pair_likelihoods(const lk_table_type &lk,
		const vector<double> &rmap,
		const int start_idx, const int end_idx,
		const vector<vector<vector<double> > > &pair_lks,
		const int lk_SNP_window)
{
	double clk=0.0;
	const int N_rho_minus_1 = lk.N_rho-1;
	const double multiplier = N_rho_minus_1 / lk.lk_rho_max;

	for (int i=start_idx; i<end_idx; ++i)
	{
		int j_min = max(i-lk_SNP_window, 0);
		int j_max = min(i+lk_SNP_window+1, end_idx);
		int j_mid = min(max(j_min, i), j_max);
		double rmap_i = rmap[i];

		double rho_dist;
		int rho_idx, rho_idx2;
		double y1, dx, rho_idx_point;

		// Left of i
		for (int j=j_min; j<j_mid; ++j)
		{
			rho_dist = rmap_i - rmap[j];										// Map distance
			rho_idx = min((int)(rho_dist * multiplier), N_rho_minus_1);	// Lower index of map distance in table
			rho_idx2 = min(N_rho_minus_1, rho_idx+1);								// Upper index of map distance
			y1 = pair_lks[i][j][rho_idx];											// clk in table

			// Do interpolation of composite likelihood
			clk += y1;
			if (rho_idx2 > rho_idx)
			{
				rho_idx_point = lk.lk_rho_points[rho_idx];								// Rho at lower map index
				dx = lk.lk_rho_points[rho_idx2] - rho_idx_point;						// Difference in rho between indices
				clk += ((pair_lks[i][j][rho_idx2] - y1) / dx) * (rho_dist - rho_idx_point);
			}
		}

		// Right of i
		for (int j=j_mid+1; j<j_max; ++j)
		{
			rho_dist = rmap[j] - rmap_i;										// Map distance
			rho_idx = min((int)(rho_dist * multiplier), N_rho_minus_1);	// Lower index of map distance in table
			rho_idx2 = min(N_rho_minus_1, rho_idx+1);								// Upper index of map distance
			y1 = pair_lks[i][j][rho_idx];											// clk in table

			// Do interpolation of composite likelihood
			clk += y1;
			if (rho_idx2 > rho_idx)
			{
				rho_idx_point = lk.lk_rho_points[rho_idx];								// Rho at lower map index
				dx = lk.lk_rho_points[rho_idx2] - rho_idx_point;						// Difference in rho between indices
				clk += ((pair_lks[i][j][rho_idx2] - y1) / dx) * (rho_dist - rho_idx_point);
			}
		}
	}
	return clk;
}


/*
double calc_composite_likelihood_using_pair_likelihoods(const lk_table_type &lk,
		const vector<double> &rmap,
		const int start_idx, const int end_idx,
		const vector<vector<vector<double> > > &pair_lks,
		const int lk_SNP_window)
{
	double clk=0.0;
	const int window=lk_SNP_window;
	const unsigned int N_rho_minus_1 = lk.N_rho-1;
	const double multiplier = N_rho_minus_1 / lk.lk_rho_max;
	unsigned int rho_idx, rho_idx2;
	int i, j, j_min, j_max;
	double rho_dist, rmap_i, rho_idx_point, dx, y1;

	for (i=start_idx; i<end_idx; ++i)
	{
		j_min = max(i-window, 0);
		j_max = min(i+window, end_idx-1);
		rmap_i = rmap[i];

		for (j=j_min; j<=j_max; ++j)
		{
			if (i==j)
				continue;
			rho_dist = fabs(rmap_i - rmap[j]);										// Map distance
			rho_idx = min((unsigned int)(rho_dist * multiplier), N_rho_minus_1);	// Lower index of map distance in table
			rho_idx2 = min(N_rho_minus_1, rho_idx+1);								// Upper index of map distance
			rho_idx_point = lk.lk_rho_points[rho_idx];								// Rho at lower map index
			dx = lk.lk_rho_points[rho_idx2] - rho_idx_point;						// Difference in rho between indices
			y1 = pair_lks[i][j][rho_idx];											// clk in table

			// Do interpolation of composite likelihood
			clk += y1;
			if (dx > 0)
				clk += ((pair_lks[i][j][rho_idx2] - y1) / dx) * (rho_dist - rho_idx_point);
		}
	}
	return clk;
}
*/

void read_locs_file(const string &loc_filename, data_struct &data)
{
	printLOG("Reading Loci File: " + loc_filename + "\n");
	ifstream in(loc_filename.c_str());
	if (!in.good())
		error("Could not open loc file: " + loc_filename);

	double tlseq;
	char linear_or_circular;

	stringstream ss;
	string line;
	getline(in, line);
	ss.clear(); ss.str(line);
	ss >> data.N_sites >> tlseq >> linear_or_circular;
	if (linear_or_circular != 'L')
		error("Only Linear Genomes Supported.");

	data.locs.clear();
	data.locs.reserve(data.N_sites);
	double tmp=0.0, last=-numeric_limits<double>::max();

	while(!in.eof())
	{
		getline(in,line);
		if (line.size() == 0)
			continue;

		ss.clear(); ss.str(line);
		ss >> tmp;
		data.locs.push_back(tmp);

		if (tmp <= last)
			error("Loci must be strictly increasing. " + dbl2str(last, 4) + " " + dbl2str(tmp,4));
		last = tmp;
	}

	if (data.N_sites != data.locs.size())
		error("Read " +int2str(data.locs.size()) + " sites, but header says " +int2str(data.N_sites) + ".");

	printLOG("\tRead " + int2str(data.N_sites) + " sites\n");

	in.close();
}

void read_seqs_file(const string &seq_filename, data_struct &data)
{
	printLOG("Reading Sequence File:" + seq_filename + "\n");
	data.is_haplotypes = false;

	ifstream in(seq_filename.c_str());
	if (!in.good())
		error("Could not open seq file: " + seq_filename);

	unsigned int N_seq;
	unsigned int N_sites;
	int haplotype_or_genotype;

	stringstream ss;
	string line;
	getline(in, line);
	ss.clear(); ss.str(line);
	ss >> N_seq >> N_sites >> haplotype_or_genotype;

	if (haplotype_or_genotype != 1)
		error("Only Haplotypes Supported");

	data.N_seq = N_seq;
	data.N_sites = N_sites;

	if ((N_seq < 2) || (N_sites < 2))
		error("Insufficient Data");
	if ((haplotype_or_genotype != 1) && (haplotype_or_genotype != 2))
		error("Unknown Data Type (1=haplotypes, 2=genotypes). Read: " + int2str(haplotype_or_genotype));
	if (haplotype_or_genotype == 1)
		data.is_haplotypes=true;

	data.headers.resize(N_seq);
	data.seqs.resize(N_seq);
	for (unsigned int ui=0; ui<N_seq; ui++)
		data.seqs[ui].reserve(N_sites);

	data.freq_counts.resize(N_sites, 0);

	int seq_num=-1;
	while(!in.eof())
	{
		getline(in, line);
		line.erase( line.find_last_not_of(" \t\r\n") + 1);	// Trim whitespace at end of line
		if (line.size() == 0)
			continue;

		if (line[0] == '>')
		{	// Header
			seq_num++;
			if (unsigned(seq_num+1) > N_seq)
				error("Too Many Sequences.");
			data.headers[seq_num] = line.substr(1);
		}
		else
		{
			for (unsigned int ui=0; ui<line.size(); ui++)
			{
				char c = line[ui];
				switch (c)
				{
				case '0':
					data.seqs[seq_num].push_back(0);
					break;
				case '1':
					data.seqs[seq_num].push_back(1);
					data.freq_counts[data.seqs[seq_num].size()-1]++;
					break;
				case ' ':
				case '\t':
					break;
				default:
					error("Unknown base in seq " + data.headers[seq_num] + ", " + c);
					break;
				}
			}
		}
	}
	in.close();

	if (unsigned(seq_num+1) < N_seq)
		error("Too Few Sequences");

	for (unsigned int ui=0; ui<N_seq; ui++)
		if (data.seqs[ui].size() != N_sites)
			error("Incorrect number of sites in sequence " + data.headers[ui] + "\nExpected:" + int2str(data.N_sites) + " Read:" + int2str(data.seqs[ui].size()) );

	printLOG("\tRead " + int2str(N_sites) + " sites\n");
	printLOG("\tRead " + int2str(N_seq) + " haplotypes\n");
}

void read_lk_file(const string &lk_filename, lk_table_type &lk)
{
	printLOG("Reading lk file: " + lk_filename + "\n");
	ifstream in(lk_filename.c_str());
	if (!in.good())
		error("Could not open lk file: " + lk_filename);

	lk.lk_table.clear();

	stringstream ss;
	string line;

	unsigned int N_types, N_theta;
	double theta;

	line = "";
	while(line.size() == 0)
		getline(in,line);

	ss.clear(); ss.str(line);
	ss >> lk.N_seq >> N_types;

	getline(in,line); ss.clear(); ss.str(line);
	ss >> N_theta >> theta;

	getline(in,line); ss.clear(); ss.str(line);
	ss >> lk.N_rho >> lk.lk_rho_max;

	lk.lk_rho_points.resize(lk.N_rho);
	for (unsigned int i=0; i<lk.N_rho; i++)
	{
		double r = (double(i) * lk.lk_rho_max) / (lk.N_rho-1);
		lk.lk_rho_points[i] = r;
	}

	unsigned int idx, N00, N01, N10, N11;
	char c;
	double dtmp;
	unsigned int max_N00=lk.N_seq, max_N01=(lk.N_seq/2), max_N10=(lk.N_seq/2), max_N11=(lk.N_seq/2);

	printLOG("\tExpecting " + int2str(lk.N_rho) + " rho entries\n");
	printLOG("\tExpecting " + int2str(N_types) + " lk entries\n");
	printLOG("\tExpecting " + int2str(max_N00) + ", " + int2str(max_N01) + ", " + int2str(max_N10) + ", " + int2str(max_N11) + " of each pair-type\n");

	lk.lk_table.clear();

	//lk.lk_table = vec5D(max_N00+1, vec4D(max_N01+1, vec3D(max_N10+1, vec2D(max_N11+1, vec1D(lk.N_rho, 0)))));

	unsigned int rho_idx;
	unsigned int N_entries=0;
	while(!in.eof())
	{
		getline(in,line);
		if ((line.size() == 0) || (line[0] == 'T'))
			continue;	// Must be a header line.

		line.erase( line.find_last_not_of(" \t\r\n") + 1);	// Trim whitespace at end of line

		// Must be a data line
		ss.clear(); ss.str(line);
		ss >> idx >> c >> N00 >> N01 >> N10 >> N11 >> c;
		lk.lk_table[N00][N01][N10][N11].resize(lk.N_rho);
		//cout << N00 << " " << N01 << " " << N10 << " " << N11 << endl;
		rho_idx = 0;
		while (!ss.eof())
		{
			ss >> dtmp;
			lk.lk_table[N00][N01][N10][N11][rho_idx] = dtmp;
			rho_idx++;
			if (rho_idx > lk.N_rho)
				error("eh2?");
		}
		N_entries++;
	}
	printLOG("\tRead " + int2str(N_entries) + " entries\n");
}

// Use "Fewer permutations, more accurate p-values" algorithm of Knijnenburg et al. 2009
double find_tail_approximation_p_value(const double clk_ratio, const vector<double> &clk_ratio_distribution_in)
{
	double p = numeric_limits<double>::quiet_NaN();
	double k, a;

	unsigned int M=0;
	for (unsigned int ui=0; ui<clk_ratio_distribution_in.size(); ui++)
	{
		if (clk_ratio < clk_ratio_distribution_in[ui])
			M++;
	}

	if ((M < 10) && (clk_ratio_distribution_in.size() > 999))
	{	// clk_ratio is in the extreme tails of the distribution, and have enough data to use the tail approximation
		vector<double> clk_ratio_distribution = clk_ratio_distribution_in;	// Create local copy
		sort(clk_ratio_distribution.begin(), clk_ratio_distribution.end());
		vector<double> clk_distribution_tail;

		bool passed_gof=false;
		unsigned int min_idx;
		for (min_idx = (unsigned int)(double(clk_ratio_distribution.size())*0.875);	// Define the tails
			(min_idx < clk_ratio_distribution.size()-25);	// Require at least 25 points to estimate GPD.
			min_idx += 5)
		{
			double min_clk = (clk_ratio_distribution[min_idx] + clk_ratio_distribution[min_idx+1])*0.5;

			// Extract tail of the distribution
			clk_distribution_tail.resize(0);
			for (unsigned int ui=min_idx+1; ui<clk_ratio_distribution.size(); ui++)
			{
				if (clk_ratio_distribution[ui] > min_clk)
					clk_distribution_tail.push_back(clk_ratio_distribution[ui] - min_clk);
			}

			gpd_ml(clk_distribution_tail, k, a);	// Estimate GPD by Maximum Likelihood
			passed_gof = test_gof(clk_distribution_tail, k, a);	// Test Goodness of Fit
			if (passed_gof)
			{	// Estimate p-value
				p = (1.0 - gp_cdf(clk_ratio-min_clk, k, a)) * clk_distribution_tail.size() / double(clk_ratio_distribution.size());
				break;
			}
		}
	}

	return p;
}

class model
{
public:
	model()
	{
		bg_rate = 0; hot_rate = 0; clk = 0; line_coord = 0;
	};
	~model() {};
	bool operator== (const model &b)
	{
		return ((bg_rate == b.bg_rate) && (hot_rate == b.hot_rate) && (clk == b.clk));
	};
	bool operator!= (const model &b)
	{
		return (!((bg_rate == b.bg_rate) && (hot_rate == b.hot_rate) && (clk == b.clk)));
	};
	double bg_rate;
	double hot_rate;
	double clk;
	double line_coord;

	static bool compare_by_bg_rate(const model &a, const model &b)
	{
	    return a.bg_rate < b.bg_rate;
	}

	static bool compare_by_hot_rate(const model &a, const model &b)
	{
	    return a.hot_rate < b.hot_rate;
	}

	static bool compare_by_clk(const model &a, const model &b)
	{
	    return a.clk < b.clk;
	}

	static bool compare_by_line_coord(const model &a, const model &b)
	{
	    return a.line_coord < b.line_coord;
	}

	void print()
	{
		cout << "bg:" << bg_rate << " hot:" << hot_rate << " clk:" << clk << endl;
	}
};

void build_rmap(vector<double> &out, const data_struct &data,
		const double lh_hotspot_pos, const double rh_hotspot_pos,
		const double bg_rate, const double hot_rate)
{
	if ((hot_rate < 0) || (hot_rate != hot_rate))
		printLOG("Warning - negative or NaN hotspot rate.\n");
	if ((bg_rate < 0) || (bg_rate != bg_rate))
		printLOG("Warning - negative or NaN background rate.\n");
	//assert(bg_rate >= 0);
	//assert(hot_rate >= 0);
	out.resize(data.locs.size());
	out[0] = 0;
	for (unsigned int ui=1; ui<out.size(); ui++)
	{
		if (data.locs[ui] <= lh_hotspot_pos)
			out[ui] = bg_rate * (data.locs[ui] - data.locs[0]);
		else if ((data.locs[ui] > lh_hotspot_pos) && (data.locs[ui] <= rh_hotspot_pos))
		{
			out[ui] = bg_rate * (lh_hotspot_pos - data.locs[0]);
			out[ui] += hot_rate * (data.locs[ui] - lh_hotspot_pos);
		}
		else
		{
			out[ui] = bg_rate * (lh_hotspot_pos - data.locs[0]);
			out[ui] += hot_rate * (rh_hotspot_pos - lh_hotspot_pos);
			out[ui] += bg_rate * (data.locs[ui] - rh_hotspot_pos);
		}
	}
}

double find_bg_maximum(const data_struct &data,
		const lk_table_type &lk, const int lk_SNP_window,
		const double lh_hotspot_pos, const double rh_hotspot_pos,
		double &bg_mle_rate, const double hot_mle_rate, const vector<vector<vector<double> > > &pair_lks)
{
	double max_clk_difference = 0.01;	// Difference allowed when converging just background or hotspot
	double max_bg_rate_difference = 1e-6;

	vector< model > range(4);
	range[0].bg_rate = 0.0;
	range[3].bg_rate = 100.0;
	range[1].bg_rate = (range[0].bg_rate + range[3].bg_rate)*0.33;
	range[2].bg_rate = (range[0].bg_rate + range[3].bg_rate)*0.66;

	vector<double> rmap;
	for (int uj=0; uj<4; ++uj)
	{
		build_rmap(rmap, data, lh_hotspot_pos, rh_hotspot_pos, range[uj].bg_rate, hot_mle_rate);
		range[uj].clk = calc_composite_likelihood_using_pair_likelihoods(lk, rmap, 0, rmap.size(), pair_lks, lk_SNP_window);
	}

	double clk_max_estimate = -numeric_limits<double>::max();
	for (unsigned int ui=0; ui<4; ui++)
	{
		if (range[ui].clk > clk_max_estimate)
		{
			bg_mle_rate = range[ui].bg_rate;
			clk_max_estimate = range[ui].clk;
		}
	}

	// Vary background rate, keep hotspot rate fixed at mle
	unsigned int replacement_idx;
	for (unsigned int step=0; step<100; step++)
	{
		sort(range.begin(), range.end(), model::compare_by_bg_rate);
		// Step 1. Find internal point closest to maxima
		replacement_idx = 0;	// Replace left boundary
		if (range[1].clk > range[2].clk)
			replacement_idx = 3;	// Actually, point 1 is closer to the maxima, so replace right boundary

		double new_bg_rate;
		if (replacement_idx == 0)
			//new_bg_rate = 0.5*(range[1].bg_rate + range[3].bg_rate);
			if ((range[3].bg_rate-range[2].bg_rate) > (range[2].bg_rate-range[1].bg_rate))
				new_bg_rate = range[2].bg_rate + 0.38197*(range[3].bg_rate-range[2].bg_rate);
			else
				new_bg_rate = range[2].bg_rate - 0.38197*(range[2].bg_rate-range[1].bg_rate);
		else
			//new_bg_rate = 0.5*(range[0].bg_rate + range[2].bg_rate);
			if ((range[2].bg_rate-range[1].bg_rate) > (range[1].bg_rate-range[0].bg_rate))
				new_bg_rate = range[1].bg_rate + 0.38197*(range[2].bg_rate-range[1].bg_rate);
			else
				new_bg_rate = range[1].bg_rate - 0.38197*(range[1].bg_rate-range[0].bg_rate);

		if ((new_bg_rate == range[1].bg_rate) || (new_bg_rate == range[2].bg_rate) || (new_bg_rate != new_bg_rate))
			new_bg_rate = 0.5*(range[1].bg_rate + range[2].bg_rate);

		if ((new_bg_rate < 0) || (new_bg_rate != new_bg_rate))
			new_bg_rate = numeric_limits<double>::epsilon();
		build_rmap(rmap, data, lh_hotspot_pos, rh_hotspot_pos, new_bg_rate, hot_mle_rate);
		range[replacement_idx].clk = calc_composite_likelihood_using_pair_likelihoods(lk, rmap, 0, rmap.size(), pair_lks, lk_SNP_window);
		range[replacement_idx].bg_rate = new_bg_rate;

		double max_diff_bg=0.0, max_diff_clk = 0.0;
		for (unsigned int ui=0; ui<3; ui++)
			for (unsigned int uj=ui+1; uj<4; uj++)
			{
				max_diff_bg = max(max_diff_bg, fabs(range[ui].bg_rate - range[uj].bg_rate));
				max_diff_clk = max(max_diff_clk, fabs(range[ui].clk - range[uj].clk));
			}

		if 	((max_diff_clk < max_clk_difference) || (max_diff_bg < max_bg_rate_difference))
		{	// Have converged to a maximum
			break;
		}
	}

	clk_max_estimate = -numeric_limits<double>::max();
	for (unsigned int ui=0; ui<4; ui++)
	{
		if (range[ui].clk > clk_max_estimate)
		{
			bg_mle_rate = range[ui].bg_rate;
			clk_max_estimate = range[ui].clk;
		}
	}
	return clk_max_estimate;
}

double find_hot_maximum(const data_struct &data,
		const lk_table_type &lk, const int lk_SNP_window,
		const double lh_hotspot_pos, const double rh_hotspot_pos,
		const double bg_mle_rate, double &hot_mle_rate, const vector<vector<vector<double> > > &pair_lks)
{
	double max_clk_difference = 0.01;	// Difference allowed when converging just background or hotspot
	double max_hot_rate_difference = 1e-6;

	vector< model > range(4);
	// Hold background rate constant, vary hotspot rate.
	range[0].hot_rate = 0.0; range[3].hot_rate = 100.0;
	if ((hot_mle_rate > 0.0) && (hot_mle_rate < 100.0))
		range[1].hot_rate = hot_mle_rate;
	else
		range[1].hot_rate = range[3].hot_rate/3.0;
	range[2].hot_rate = (range[1].hot_rate + range[3].hot_rate)*0.5;

	vector<double> rmap;
	for (int uj=0; uj<4; ++uj)
	{
		build_rmap(rmap, data, lh_hotspot_pos, rh_hotspot_pos, bg_mle_rate, range[uj].hot_rate);
		range[uj].clk = calc_composite_likelihood_using_pair_likelihoods(lk, rmap, 0, rmap.size(), pair_lks, lk_SNP_window);
	}

	// Vary hotspot rate, keep background rate fixed
	unsigned int replacement_idx;
	for (unsigned int step=0; step<100; step++)
	{
		sort(range.begin(), range.end(), model::compare_by_hot_rate);
		replacement_idx = 0;	// Replace left boundary
		if (range[1].clk > range[2].clk)
			replacement_idx = 3;	// Actually, point 1 is closer to the maxima, so replace right boundary

		double new_hot_rate;
		if (replacement_idx == 0)
			//new_hot_rate = 0.5*(range[1].hot_rate + range[3].hot_rate);
			if ((range[3].hot_rate-range[2].hot_rate) > (range[2].hot_rate-range[1].hot_rate))
				new_hot_rate = range[2].hot_rate + 0.38197*(range[3].hot_rate-range[2].hot_rate);
			else
				new_hot_rate = range[2].hot_rate - 0.38197*(range[2].hot_rate-range[1].hot_rate);
		else
			//new_hot_rate = 0.5*(range[0].hot_rate + range[2].hot_rate);
			if ((range[2].hot_rate-range[1].hot_rate) > (range[1].hot_rate-range[0].hot_rate))
				new_hot_rate = range[1].hot_rate + 0.38197*(range[2].hot_rate-range[1].hot_rate);
			else
				new_hot_rate = range[1].hot_rate - 0.38197*(range[1].hot_rate-range[0].hot_rate);

		if ((new_hot_rate == range[1].hot_rate) || (new_hot_rate == range[2].hot_rate) || (new_hot_rate != new_hot_rate))
			new_hot_rate = 0.5*(range[1].hot_rate + range[2].hot_rate);

		if ((new_hot_rate < 0) || (new_hot_rate != new_hot_rate))
			new_hot_rate = numeric_limits<double>::epsilon();
		build_rmap(rmap, data, lh_hotspot_pos, rh_hotspot_pos, bg_mle_rate, new_hot_rate);
		range[replacement_idx].clk = calc_composite_likelihood_using_pair_likelihoods(lk, rmap, 0, rmap.size(), pair_lks, lk_SNP_window);
		range[replacement_idx].hot_rate = new_hot_rate;

		double max_diff_hot=0.0, max_diff_clk = 0.0;
		for (unsigned int ui=0; ui<3; ui++)
			for (unsigned int uj=ui+1; uj<4; uj++)
			{
				max_diff_hot = max(max_diff_hot, fabs(range[ui].hot_rate - range[uj].hot_rate));
				max_diff_clk = max(max_diff_clk, fabs(range[ui].clk - range[uj].clk));
			}

		if 	((max_diff_clk < max_clk_difference) || (max_diff_hot < max_hot_rate_difference))
		{	// Have converged to a maximum
			break;
		}
	}

	double clk_max_estimate = -numeric_limits<double>::max();
	for (unsigned int ui=0; ui<4; ui++)
	{
		if (range[ui].clk > clk_max_estimate)
		{
			hot_mle_rate = range[ui].hot_rate;
			clk_max_estimate = range[ui].clk;	// Maximum found after hot step
		}
	}
	return clk_max_estimate;
}


// Find the maximum likelihood along line (dy / dx)*bg_rate + c
double find_maximum_along_line(const data_struct &data,
		const lk_table_type &lk, const int lk_SNP_window,
		const double lh_hotspot_pos, const double rh_hotspot_pos,
		double &bg_mle_rate, double &hot_mle_rate, const vector<vector<vector<double> > > &pair_lks,
		double d_bg, double d_hot)
{
	double max_clk_difference = 0.01;	// Difference allowed when converging just background or hotspot
	double max_hot_rate_difference = 1e-6;
	double max_bg_rate_difference = 1e-6;

	vector< model > range(4);

	double m = d_hot / d_bg;
	double c = hot_mle_rate - m*bg_mle_rate;
	range[1].bg_rate = bg_mle_rate; range[1].hot_rate = hot_mle_rate;

	range[0].bg_rate = 0.0;
	range[0].hot_rate = c;

	if (range[0].hot_rate < 0.0)
	{
		range[0].hot_rate = 0.0;
		range[0].bg_rate = -c / m;
	}
	else if (range[0].hot_rate > 100.0)
	{
		range[0].hot_rate = 0.0;
		range[0].bg_rate = (100.0-c) / m;
	}

	range[3].bg_rate = 100.0;
	range[3].hot_rate = m * range[3].bg_rate + c;

	if (range[3].hot_rate < 0.0)
	{
		range[3].hot_rate = 0.0;
		range[3].bg_rate = -c / m;
	}
	else if (range[3].hot_rate > 100.0)
	{
		range[3].hot_rate = 0.0;
		range[3].bg_rate = (100.0-c) / m;
	}

	range[2].bg_rate = (range[0].bg_rate + range[3].bg_rate)*0.5;
	range[2].hot_rate = (range[0].hot_rate + range[3].hot_rate)*0.5;

	vector<double> rmap;
	for (int uj=0; uj<4; ++uj)
	{
		build_rmap(rmap, data, lh_hotspot_pos, rh_hotspot_pos, range[uj].bg_rate, range[uj].hot_rate);
		range[uj].clk = calc_composite_likelihood_using_pair_likelihoods(lk, rmap, 0, rmap.size(), pair_lks, lk_SNP_window);
	}

	// Vary along line.
	unsigned int replacement_idx;
	for (unsigned int step=0; step<100; step++)
	{
		for (int uj=0; uj<4; ++uj)
		{
			double y = range[uj].hot_rate - c;
			range[uj].line_coord = sqrt(range[uj].bg_rate*range[uj].bg_rate + y*y);
		}

		sort(range.begin(), range.end(), model::compare_by_line_coord);

		replacement_idx = 0;	// Replace left boundary
		if (range[1].clk > range[2].clk)
			replacement_idx = 3;	// Actually, point 1 is closer to the maxima, so replace right boundary

		double new_line_coord, new_bg_rate, new_hot_rate;
		if (replacement_idx == 0)
			//new_line_coord = 0.5*(range[1].line_coord + range[3].line_coord);
			if ((range[3].line_coord-range[2].line_coord) > (range[2].line_coord-range[1].line_coord))
				new_line_coord = range[2].line_coord + 0.38197*(range[3].line_coord-range[2].line_coord);
			else
				new_line_coord = range[2].line_coord - 0.38197*(range[2].line_coord-range[1].line_coord);
		else
			//new_line_coord = 0.5*(range[0].line_coord + range[2].line_coord);
			if ((range[2].line_coord-range[1].line_coord) > (range[1].line_coord-range[0].line_coord))
				new_line_coord = range[1].line_coord + 0.38197*(range[2].line_coord-range[1].line_coord);
			else
				new_line_coord = range[1].line_coord - 0.38197*(range[1].line_coord-range[0].line_coord);

		if ((new_line_coord == range[1].line_coord) || (new_line_coord == range[2].line_coord) || (new_line_coord != new_line_coord))
			new_line_coord = 0.5*(range[1].line_coord + range[2].line_coord);

		double theta = atan(d_hot / d_bg);
		new_bg_rate = new_line_coord * cos(theta);
		new_hot_rate = new_line_coord * sin(theta) + c;
		if ((new_hot_rate < 0) || (new_hot_rate != new_hot_rate))
			new_hot_rate = numeric_limits<double>::epsilon();
		if ((new_bg_rate < 0) || (new_bg_rate != new_bg_rate))
			new_bg_rate = numeric_limits<double>::epsilon();
		build_rmap(rmap, data, lh_hotspot_pos, rh_hotspot_pos, new_bg_rate, new_hot_rate);
		range[replacement_idx].clk = calc_composite_likelihood_using_pair_likelihoods(lk, rmap, 0, rmap.size(), pair_lks, lk_SNP_window);
		range[replacement_idx].bg_rate = new_bg_rate;
		range[replacement_idx].hot_rate = new_hot_rate;

		double max_diff_bg=0.0, max_diff_hot=0.0, max_diff_clk = 0.0;
		for (unsigned int ui=0; ui<3; ui++)
			for (unsigned int uj=ui+1; uj<4; uj++)
			{
				max_diff_bg = max(max_diff_bg, fabs(range[ui].bg_rate - range[uj].bg_rate));
				max_diff_hot = max(max_diff_hot, fabs(range[ui].hot_rate - range[uj].hot_rate));
				max_diff_clk = max(max_diff_clk, fabs(range[ui].clk - range[uj].clk));
			}

		if ((max_diff_clk < max_clk_difference) || ((max_diff_bg < max_bg_rate_difference) & (max_diff_hot < max_hot_rate_difference)))
		{	// Have converged to a maximum
			break;
		}
	}

	double clk_max_estimate = -numeric_limits<double>::max();
	for (unsigned int ui=0; ui<4; ui++)
	{
		if (range[ui].clk > clk_max_estimate)
		{
			hot_mle_rate = range[ui].hot_rate;
			bg_mle_rate = range[ui].bg_rate;
			clk_max_estimate = range[ui].clk;	// Maximum found after hot step
		}
	}
	return clk_max_estimate;
}



// This function is *horrible*. Basically all it does is a line search for background rate, holding the hotspot
// rate fixed. It then does a line search for hotspot rate, holding background rate fixed. It repeats until
// convergence (which can be slow). The code is a mess. Mea culpa.
double estimate_MLE_hotspot_rate_across_region(const data_struct &data,
		const lk_table_type &lk, const int lk_SNP_window,
		const double lh_hotspot_pos, const double rh_hotspot_pos,
		double &bg_mle_rate, double &hot_mle_rate)
{
	double max_clk_difference = 0.01;

	vector<vector<vector<double> > > pair_lks;
	get_pair_likelihoods(data, lk, 0, data.locs.size(), pair_lks, lk_SNP_window);
	double clk_max_estimate = -numeric_limits<double>::max();
	double prev_clk_max_estimate;
	hot_mle_rate = bg_mle_rate;
	if ((bg_mle_rate < 0) | (bg_mle_rate > 100))
	{
		bg_mle_rate = 0.0; hot_mle_rate = 0.0;
	}

	double bg_max_clk = -numeric_limits<double>::max(), hot_max_clk = -numeric_limits<double>::max();
	double line_max_clk = -numeric_limits<double>::max();
	double prev_bg_rate, prev_hot_rate;
	int outer_step;

	for (outer_step = 0; outer_step<1000; outer_step++)
	{
		prev_clk_max_estimate = clk_max_estimate;
		prev_bg_rate = bg_mle_rate; prev_hot_rate = hot_mle_rate;

		// Work out what direction to move in (if any)
		bg_max_clk = find_bg_maximum(data, lk, lk_SNP_window, lh_hotspot_pos, rh_hotspot_pos,
			bg_mle_rate, hot_mle_rate, pair_lks);

		hot_max_clk = find_hot_maximum(data, lk, lk_SNP_window, lh_hotspot_pos, rh_hotspot_pos,
			bg_mle_rate, hot_mle_rate, pair_lks);

		clk_max_estimate = max(max(bg_max_clk, hot_max_clk), clk_max_estimate);

		if (fabs(prev_clk_max_estimate - clk_max_estimate) < max_clk_difference)
			break;

		// Try moving along line in direction of previous move for quicker convergence.
		if (outer_step > 2)
		{
			double dx = prev_bg_rate-bg_mle_rate, dy = prev_hot_rate-hot_mle_rate;
			if ((dx != 0) && (dy != 0))
				for (int step=0; step<2; step++)
				{
					line_max_clk = find_maximum_along_line(data, lk, lk_SNP_window, lh_hotspot_pos, rh_hotspot_pos,
						bg_mle_rate, hot_mle_rate, pair_lks, dx, pow(-1.0, step)*dy);
				}
			clk_max_estimate = max(clk_max_estimate, line_max_clk);
		}

		if (outer_step == 999)
		{
			printLOG("Warning: Convergence failure in hotspot routine.\n");
		}
	}

	return clk_max_estimate;
}


double estimate_MLE_constant_rate_across_region(const data_struct &data,
		const lk_table_type &lk, const int lk_SNP_window, double &mle_rate)
{
	double max_clk_difference = 0.01;
	double max_bg_rate_difference = 1e-6;

	vector<vector<vector<double> > > pair_lks;
	get_pair_likelihoods(data, lk, 0, data.locs.size(), pair_lks, lk_SNP_window);
	double clk_max_estimate = -numeric_limits<double>::max();
	mle_rate = 0;

	vector< model > range(4);
	range[0].bg_rate = 0; range[3].bg_rate = 100.0;
	range[1].bg_rate = (range[0].bg_rate + range[3].bg_rate)*0.33;
	range[2].bg_rate = (range[0].bg_rate + range[3].bg_rate)*0.66;

	// Find maximum clk using line search
	// http://en.wikipedia.org/wiki/Line_search

	// Initialise
	for (unsigned int uj=0; uj<4; uj++)
	{
		vector<double> rmap(data.locs.size(), 0);
		for (unsigned int ui=1; ui<rmap.size(); ui++)
			rmap[ui] = rmap[ui-1] + range[uj].bg_rate * (data.locs[ui] - data.locs[ui-1]);

		range[uj].clk = calc_composite_likelihood_using_pair_likelihoods(lk, rmap, 0, rmap.size(), pair_lks, lk_SNP_window);
	}

	// Now perform line search
	unsigned int replacement_idx;
	int step;
	for (step=0; step<100; step++)
	{
		// Step 1. Find internal point closest to maxima
		replacement_idx = 0;	// Replace left boundary
		if (range[1].clk > range[2].clk)
			replacement_idx = 3;	// Actually, point 1 is closer to the maxima, so replace right boundary

		double new_bg_rate;
		if (replacement_idx == 0)
			//new_bg_rate = 0.5*(range[1].bg_rate + range[3].bg_rate);
			if ((range[3].bg_rate-range[2].bg_rate) > (range[2].bg_rate-range[1].bg_rate))
				new_bg_rate = range[2].bg_rate + 0.38197*(range[3].bg_rate-range[2].bg_rate);
			else
				new_bg_rate = range[2].bg_rate - 0.38197*(range[2].bg_rate-range[1].bg_rate);
		else
			//new_bg_rate = 0.5*(range[0].bg_rate + range[2].bg_rate);
			if ((range[2].bg_rate-range[1].bg_rate) > (range[1].bg_rate-range[0].bg_rate))
				new_bg_rate = range[1].bg_rate + 0.38197*(range[2].bg_rate-range[1].bg_rate);
			else
				new_bg_rate = range[1].bg_rate - 0.38197*(range[1].bg_rate-range[0].bg_rate);

		if ((new_bg_rate == range[1].bg_rate) || (new_bg_rate == range[2].bg_rate) || (new_bg_rate != new_bg_rate))
			new_bg_rate = (range[1].bg_rate + range[2].bg_rate) * 0.5;

		vector<double> rmap(data.locs.size(), 0);
		for (unsigned int ui=1; ui<rmap.size(); ui++)
			rmap[ui] = rmap[ui-1] + new_bg_rate * (data.locs[ui] - data.locs[ui-1]);

		range[replacement_idx].clk = calc_composite_likelihood_using_pair_likelihoods(lk, rmap, 0, rmap.size(), pair_lks, lk_SNP_window);
		range[replacement_idx].bg_rate = new_bg_rate;

		sort(range.begin(), range.end(), model::compare_by_bg_rate);

		double max_diff_bg=0.0, max_diff_clk = 0.0;
		for (unsigned int ui=0; ui<3; ui++)
			for (unsigned int uj=ui+1; uj<4; uj++)
			{
				max_diff_bg = max(max_diff_bg, fabs(range[ui].bg_rate - range[uj].bg_rate));
				max_diff_clk = max(max_diff_clk, fabs(range[ui].clk - range[uj].clk));
			}

		if 	((max_diff_clk < max_clk_difference) || (max_diff_bg < max_bg_rate_difference))
		{	// Have converged to a maximum
			break;
		}

		if (step == 99)
			printLOG("Warning: Convergence failure in background routine.\n");
	}

	for (unsigned int ui=0; ui<4; ui++)
	{
		if (range[ui].clk > clk_max_estimate)
		{
			mle_rate = range[ui].bg_rate;
			clk_max_estimate = range[ui].clk;
		}
	}

	return clk_max_estimate;
}



void read_res_file(const string &res_filename, vector<double> &rmap, vector<double> &rmap_pos)
{
	printLOG("Reading res: " + res_filename + "\n");
	ifstream in(res_filename.c_str());
	if (!in.good())
		error("Could not open res file: " + res_filename);

	rmap.clear();
	rmap.push_back(0.0);
	rmap_pos.push_back(-numeric_limits<double>::max());
	stringstream ss;
	string line;
	double pos=-1.0, rho_rate=-1.0;
	double last_pos=-1.0, last_rho_rate=-1.0;
	getline(in, line);	// Header

	while(!in.eof() && (pos == -1.0))
	{
		getline(in,line);
		if (line.size() == 0)
			error("No data in res file.");

		ss.clear(); ss.str(line);
		ss >> pos >> rho_rate;
	}

	rmap_pos.push_back(pos);
	rmap.push_back(0);

	last_pos = pos;
	last_rho_rate = rho_rate;

	double pos_delta, rho, rho_tot = 0.0;
	while(!in.eof())
	{
		getline(in,line);
		if (line.size() == 0)
			continue;

		ss.clear(); ss.str(line);
		ss >> pos >> rho_rate;

		if (pos < 0)
			continue;

		pos_delta = pos - last_pos;
		rho = pos_delta * last_rho_rate;

		if ((pos_delta < 0) || (rho < 0))
			error("Rmap must be increasing between: " + dbl2str(last_pos, 4) + " " + dbl2str(pos,4));

		rho_tot += rho;

		rmap.push_back(rho_tot);
		rmap_pos.push_back(pos);

		last_pos = pos;
		last_rho_rate = rho_rate;
	}
	printLOG("\tRead " + int2str(rmap.size()) + " map positions\n");
	printLOG("\tMap Covers positions " + dbl2str(rmap_pos[1], 4) + " - " + dbl2str(rmap_pos[rmap_pos.size()-1],4) + " kb\n");
	printLOG("\tMap Length = " + dbl2str(rho_tot,3) + " 4Ner\n");

	//for (unsigned int ui=0; ui<rmap.size(); ui++)
	//	cout << ui << "\t" << rmap_pos[ui] << "\t" << rmap[ui] << endl;
}

// Find peaks in the recombination rate estimates
void find_maxima_locations_in_rates(const vector<double> &rmap_pos, const vector<double> &rmap_rate,
									vector<double> &maxima_pos_LHS, vector<double> &maxima_pos_RHS)
{
	double prev_rate = -1, next_rate = -1;
	for (unsigned int ui = 1; ui < (rmap_rate.size()-1); ui++)
	{
		unsigned int uj1 = ui - 1;
		while(uj1 != 0)
		{
			prev_rate = rmap_rate[uj1];
			if (fabs(prev_rate - rmap_rate[ui]) > 1e-9)
				break;
			uj1--;
		}

		unsigned int uj2 = ui + 1;
		while(uj2 < rmap_rate.size())
		{
			next_rate = rmap_rate[uj2];
			if (fabs(next_rate - rmap_rate[ui]) > 1e-9)
				break;
			uj2++;
		}

		if ((prev_rate < rmap_rate[ui]) && (next_rate < rmap_rate[ui]))
		{
			maxima_pos_LHS.push_back(rmap_pos[ui]);
			maxima_pos_RHS.push_back(rmap_pos[ui+1]);
		}
	}
}

double get_median_bgrate_in_region(const vector<double> &rmap_pos, const vector<double> &rmap_rate,
									double lh_window_pos, double rh_window_pos, double lh_hotspot_pos, double rh_hotspot_pos)
{
	double out = 0;
	vector<double> rates;
	rates.reserve(rmap_rate.size());
	for (unsigned int ui = 0; ui < (rmap_pos.size()-1); ui++)
	{
		if ((rmap_pos[ui] >= lh_window_pos) && (rmap_pos[ui+1] <= rh_window_pos))
			if ((rmap_pos[ui+1] <= lh_hotspot_pos) || (rmap_pos[ui] >= rh_hotspot_pos))
				rates.push_back(rmap_rate[ui]);	// Rate is inside window, but outside hotspot
	}
	if (rates.size() > 0)
	{
		sort(rates.begin(), rates.end());
		if ((rates.size() % 2) == 0)
		{
			unsigned int idx = rates.size()/2;
			out = (rates[idx] + rates[idx-1])*0.5;
		}
		else
		{
			unsigned int idx = (rates.size()+1)/2;
			out = rates[idx-1];
		}
	}
	return out;
}

double get_avg_rate_across_region(const vector<double> &rmap_pos, const vector<double> &rmap,
									double lh_window_pos, double rh_window_pos)
{
	double out = 0;

	unsigned int idx1 = 0, idx2 = rmap_pos.size()-1;
	for (unsigned int ui = 0; ui < (rmap_pos.size()-1); ui++)
	{
		if (rmap_pos[ui] <= lh_window_pos)
			idx1 = ui;
		if (rmap_pos[ui] <= rh_window_pos)
			idx2 = ui;
		if (rmap_pos[ui] > rh_window_pos)
			break;
	}

	double dx1 = lh_window_pos - rmap_pos[idx1];
	double dx2 = rh_window_pos - rmap_pos[idx2];
	double rmap1 = rmap[idx1];
	if (dx1 > 0)
		rmap1 += dx1 * (rmap[idx1+1] - rmap[idx1]) / (rmap_pos[idx1+1] - rmap_pos[idx1]);

	double rmap2 = rmap[idx2];
	if (dx2 > 0)
		rmap2 += dx2 * (rmap[idx2+1] - rmap[idx2]) / (rmap_pos[idx2+1] - rmap_pos[idx2]);

	out = (rmap2 - rmap1) / (rh_window_pos - lh_window_pos);

	return out;
}


