/*
 * my_ldhot.cpp
 *
 *  Created on: Apr 16, 2010
 *      Author: Adam Auton
 */

#include "my_ldhot.h"

// Output results data to file
void output_result(ofstream &out, double lh_hotspot_pos, double rh_hotspot_pos,
		double mle_rate_hotspot,
		double lh_window_pos, double rh_window_pos,
		double mle_rate_background, double mle_rate_constant,
		unsigned int N_sims_used, double p_value, double p_value_approx)
{
	out << fixed << setprecision(3);
	out << lh_hotspot_pos << "\t" << rh_hotspot_pos;
	out << scientific;
	out << "\t" << mle_rate_hotspot;
	out << fixed << setprecision(3);
	out << "\t" << lh_window_pos << "\t" << rh_window_pos;
	out << scientific;
	out << "\t" << mle_rate_background;
	out << "\t" << mle_rate_constant;
	out << fixed << setprecision(3);
	out << "\t" << N_sims_used;
	out << scientific;
	out << "\t" << p_value;
	out << "\t" << p_value_approx;
	out << endl;
}

int main(int argc, char *argv[])
{
	time_t start,end;
	time(&start);
	// 1. Read Parameters
	parameters params(argc, argv);
	params.read_parameters();

	LOG.open((params.output_prefix + ".log").c_str());

	printLOG("\nmyLDhot - " + LDHOT_VERSION + "\n");
	printLOG("(C) Adam Auton 2014\n\n");

	params.print_params();

	data_struct data;
	Ran rangen(params.seed);	// Random number generator

	// 2. Load Locs file
	read_locs_file(params.loc_filename, data);

	// 3. Load Seq file
	read_seqs_file(params.seq_filename, data);
	if (data.is_haplotypes == false)
		error("Only Haplotype Data Supported.");

	if (data.locs.size() != data.seqs[0].size())
		error("Locs file does not match Seqs file.");

	// Remove fixed sites.
	unsigned int N_orig_sites = data.locs.size();
	data.remove_fixed_sites();
	//data.remove_singleton_sites();
	if (N_orig_sites != data.locs.size())
	{
		int diff = N_orig_sites - data.locs.size();
		printLOG("Removed " + int2str(diff) + " uninformative sites.\n");
	}

	// 4. Load rate file
	vector<double> rmap;
	vector<double> rmap_pos;
	read_res_file(params.res_filename, rmap, rmap_pos);

	vector<double> rmap_rate(rmap.size()-1, 0);
	for (unsigned int ui = 0; ui < rmap_rate.size(); ui++)
		rmap_rate[ui] = (rmap[ui+1] - rmap[ui]) / (rmap_pos[ui+1] - rmap_pos[ui]);

	// Find rate maxima locations in rate estimates
	vector<double> maxima_pos_LHS;	// LHS
	vector<double> maxima_pos_RHS;	// RHS
	find_maxima_locations_in_rates(rmap_pos, rmap_rate, maxima_pos_LHS, maxima_pos_RHS);

	// 5. Load Likelihood File
	lk_table_type lk;
	read_lk_file(params.lk_filename, lk);

	if (lk.N_seq != data.N_seq)
		error("LK file does not match Seq File.");

	if (data.N_sites < 100)
		error("Cowardly refusing to process data, as contains less than 100 SNPs.");

	ofstream out((params.output_prefix + ".hotspots.txt").c_str());
	out << "HotStart\tHotEnd\tMLE_hotspot_rate\tBgStart\tBgEnd\tMLE_bg_rate\tMLE_constant_rate\tN_used_sims\tP_ecdf\tP_tail_approx" << endl;
	// 6. Loop Over Regions
	double startpos = max(params.start_pos, data.locs[0]);
	double endpos = min(params.end_pos, data.locs[data.N_sites-1]);

	printLOG("Attempting to detect hotspots between " + dbl2str(startpos,3) + "kb and " + dbl2str(endpos,3) + "kb\n");

	const int sim_block_size = 50;

	for (double pos=startpos; pos<endpos; pos+=params.pos_step)
	{
		double lh_hotspot_pos = max(pos - params.hotspot_distance - numeric_limits<double>::epsilon(), data.locs[0]);
		double rh_hotspot_pos = min(pos + params.hotspot_distance + numeric_limits<double>::epsilon(), data.locs[data.N_sites-1]);
		// SNP index of hotspot internal sites
		int lh_hotspot_idx = distance(data.locs.begin(), lower_bound(data.locs.begin(), data.locs.end(), lh_hotspot_pos));
		int rh_hotspot_idx = distance(data.locs.begin(), upper_bound(data.locs.begin(), data.locs.end(), rh_hotspot_pos))-1;

		// Define window as +- 50 SNPs of hotspot
		int lh_window_idx = 0, rh_window_idx = data.N_sites-1;
		if (lh_hotspot_idx > params.lk_SNP_window)
			lh_window_idx = lh_hotspot_idx - params.lk_SNP_window;
		if (rh_window_idx - rh_hotspot_idx > params.lk_SNP_window)
			rh_window_idx = rh_hotspot_idx + params.lk_SNP_window;

		double lh_window_pos = data.locs[lh_window_idx] - numeric_limits<double>::epsilon();
		double rh_window_pos = data.locs[rh_window_idx] + numeric_limits<double>::epsilon();

		// If needed, expand window to meet window size requirements
		if (rh_window_pos - lh_window_pos < 2.0*params.window_distance)
		{
			double new_lh_window_pos = max(pos - params.window_distance - numeric_limits<double>::epsilon(), data.locs[0]);
			double new_rh_window_pos = min(pos + params.window_distance + numeric_limits<double>::epsilon(), data.locs[data.N_sites-1]);

			lh_window_pos = min(lh_window_pos, new_lh_window_pos);
			rh_window_pos = max(rh_window_pos, new_rh_window_pos);
		}

		lh_window_idx = distance(data.locs.begin(), lower_bound(data.locs.begin(), data.locs.end(), lh_window_pos));
		rh_window_idx = distance(data.locs.begin(), upper_bound(data.locs.begin(), data.locs.end(), rh_window_pos))-1;

		int N_SNPs_to_left = lh_hotspot_idx - lh_window_idx;
		int N_SNPs_to_right = rh_window_idx - rh_hotspot_idx;

		// Get data subset
		data_struct data_subset = data.get_subset(lh_window_idx, rh_window_idx);

		// Is this a peak in the rate file?
		bool do_test = false;
		for (unsigned int ui=0; ui<maxima_pos_LHS.size(); ui++)
		{
			if ((lh_hotspot_pos >= maxima_pos_LHS[ui]) && (lh_hotspot_pos <= maxima_pos_RHS[ui]))
				do_test = true;
			else if ((rh_hotspot_pos >= maxima_pos_LHS[ui]) && (rh_hotspot_pos <= maxima_pos_RHS[ui]))
				do_test = true;
			else if ((lh_hotspot_pos <= maxima_pos_LHS[ui]) && (rh_hotspot_pos >= maxima_pos_RHS[ui]))
				do_test = true;
			else if ((lh_hotspot_pos >= maxima_pos_LHS[ui]) && (rh_hotspot_pos <= maxima_pos_RHS[ui]))
				do_test = true;

			if ((do_test == true) || (maxima_pos_LHS[ui] > rh_hotspot_pos))
				break;
		}
		if (do_test == false)
		{
			//printLOG("\n\tNot a rate maxima.\n");
			continue;
		}

		//double res_median_rate = get_median_bgrate_in_region(rmap_pos, rmap_rate, lh_window_pos, rh_window_pos, lh_hotspot_pos, rh_hotspot_pos);
		double mean_rate = get_avg_rate_across_region(rmap_pos, rmap, lh_window_pos, rh_window_pos);

		printLOG("Hotspot: " + dbl2str_fixed(lh_hotspot_pos,3) + "-" + dbl2str_fixed(rh_hotspot_pos,3) + "kb  BG: ");
		printLOG(dbl2str_fixed(lh_window_pos,3) + "-" + dbl2str_fixed(rh_window_pos,3) + "kb ");
		printLOG("(" + int2str(data_subset.N_sites) + " sites, rate=" + dbl2str(mean_rate, 3) + ") ");

		if ((N_SNPs_to_left < params.lk_SNP_window) || (N_SNPs_to_right < params.lk_SNP_window))
		{
			printLOG("\n\tToo near boundary.\n");
			continue;
		}

		// Estimate constant background rate of real data (model 0).
		double mle_rate_constant=0.0;
		double data_clk_constant = estimate_MLE_constant_rate_across_region(data_subset, lk, params.lk_SNP_window, mle_rate_constant);
		//  Estimate constant background rate with hotspot at center (model 1).
		double mle_rate_hotspot=mle_rate_constant, mle_rate_background = mle_rate_constant;	// Take best constant rate as starting point
		double data_clk_hotspot=estimate_MLE_hotspot_rate_across_region(data_subset,
				lk, params.lk_SNP_window,
				lh_hotspot_pos, rh_hotspot_pos,
				mle_rate_background, mle_rate_hotspot);
		// Calculate likelihood ratio
		double data_clk_ratio = data_clk_hotspot - data_clk_constant;
		//printLOG("dCLK=" + dbl2str(data_clk_ratio,3) + " ");
		printLOG("Hot=" + dbl2str(mle_rate_hotspot,3) + ", bg=" + dbl2str(mle_rate_background,3) + " vs bg=" + dbl2str(mle_rate_constant,3) + " dclk=" + dbl2str(data_clk_ratio, 3) + "\n");

		unsigned int N_sims_used = 0;
		double p_value = 1.0;
		double p_value_approx = numeric_limits<double>::quiet_NaN();

		if (mle_rate_hotspot <= mle_rate_background)
		{
			printLOG("\tHotspot MLE rate <= background rate.\n");
			continue;
		}
		if (mle_rate_background > 50.0)
		{
			printLOG("\tBackground MLE rate > 50 4Ner / kb.\n");
			continue;
		}

		vector<double> clk_ratio_distribution;
		clk_ratio_distribution.reserve(params.N_sims);

		unsigned int N_lt_truth = 0;
		clock_t sim_start = clock();
		for (int sim_num=0; sim_num < params.N_sims; sim_num+=sim_block_size)
		{
			// We are going simulate data with a constant rate, but with rate drawn from exponential distribution
			// with mean mle_constant_rate (as estimated from data).

			// We're going submit simulations in batches of 50, with multithreading. In order to ensure things
			// remain deterministic, we're going to generate the random numbers/seeds we need ahead of time.
			// Update: Multithreading actually slowed things down as it appears to be bus limited :-(
			clk_ratio_distribution.resize(sim_num+sim_block_size);
			vector<int> random_seeds(sim_block_size);	// Choose random seeds for the various threads
			vector<double> sim_rate(sim_block_size);
			for (int ui=0; ui<sim_block_size; ui++)	// Do this beforehand for thread safety
			{
				random_seeds[ui] = rangen.ran_int();
				// We need to choose a null model rate to simulate with.
				// The MLE performs poorly when the background rate is very low (or zero).
				// Using the rate from LDhat seems to perform better.
				//sim_rate[ui] = rangen.rexp(res_median_rate * log(2.0)); // * log(2) for median
				//sim_rate[ui] = res_median_rate;
				//sim_rate[ui] = mean_rate;
				sim_rate[ui] = rangen.rexp(mean_rate); // * log(2) for median
				sim_rate[ui] = max(min(sim_rate[ui], 50.0), 0.0001);	// Set a minimum and maximum...
				//sim_rate[ui] = rangen.ran_double()*100;	// Uniform
			}

			for (int sub_sim=0; sub_sim<sim_block_size; sub_sim++)
			{
				// Generate null recombination map.
				vector<double> null_rmap(data_subset.N_sites, 0);
				for (unsigned int ui=1; ui < data_subset.N_sites; ui++)
					null_rmap[ui] = null_rmap[ui-1] + sim_rate[sub_sim] * (data_subset.locs[ui] - data_subset.locs[ui-1]);

				// Create a controller object for simulations under the null distribution
				sim_control con(data_subset.N_seq, data_subset.N_sites, 1.0, null_rmap, true, true, data_subset.locs, params.freq_cond, data_subset.freq_counts, params.freq_cond_model);
				// Create simulation object.
				sim simulation(con);
				simulation.set_seed(random_seeds[sub_sim]);	// Keeps things deterministic...
				simulation.run_sim();	// Generate a simulated dataset under the null recombination rate

				data_struct simulated_data;
				simulated_data.N_seq = data_subset.N_seq;
				simulated_data.N_sites = simulation.locs.size();
				simulated_data.seqs = simulation.seqs;
				simulated_data.locs = simulation.locs;

				// Estimate constant background rate of simulated data (model 0).
				double sim_mle_rate_constant=0.0;
				double sim_clk_constant = estimate_MLE_constant_rate_across_region(simulated_data, lk, params.lk_SNP_window, sim_mle_rate_constant);
				// Estimate constant background rate with hotspot of simulated data (model 1)
				double sim_mle_rate_background=0.0, sim_mle_rate_hotspot=0.0;
				double sim_clk_hotspot = estimate_MLE_hotspot_rate_across_region(simulated_data, lk, params.lk_SNP_window,
						lh_hotspot_pos, rh_hotspot_pos,	sim_mle_rate_background, sim_mle_rate_hotspot);
				// Calculate likelihood ratio
				double sim_clk_ratio = sim_clk_hotspot - sim_clk_constant;

				// Add a tiny jitter to break ties when calculating empirical p-value
				sim_clk_ratio += pow(-1.0, (double)sim_num) * 0.00001 * rangen.ran_double();

				clk_ratio_distribution[sim_num+sub_sim]=sim_clk_ratio;
			}

			N_lt_truth = 0;
			for (unsigned int ui=0; ui<clk_ratio_distribution.size(); ui++)
				if (data_clk_ratio <= clk_ratio_distribution[ui])
					N_lt_truth++;
			N_sims_used = clk_ratio_distribution.size();

			p_value = double(N_lt_truth+2.7)/(N_sims_used+5.4);	// Adjusted Wald interval (Agresti and Coull, 1998)
			double ci = 2.575829303549*sqrt(p_value*(1.0-p_value)/N_sims_used);	// 99% C.I.
			if ((p_value - ci) > 0.01)
				break;	// Don't bother carrying on - very unlikely to find anything significant.
		}
		clock_t sim_end = clock();
		double sim_time = double(sim_end - sim_start)/CLOCKS_PER_SEC;

		N_lt_truth = 0;
		for (unsigned int ui=0; ui<clk_ratio_distribution.size(); ui++)
			if (data_clk_ratio <= clk_ratio_distribution[ui])
				N_lt_truth++;
		N_sims_used = clk_ratio_distribution.size();

		p_value = double(N_lt_truth+1)/(N_sims_used+1);
		p_value_approx = find_tail_approximation_p_value(data_clk_ratio, clk_ratio_distribution);
		printLOG("\tp = " + dbl2str(p_value,4) + " (" + int2str(N_lt_truth) + "/" + int2str(N_sims_used) + ")");
		if (p_value_approx == p_value_approx)	// i.e. isn't NaN.
			printLOG(" p_tail_approx = " + dbl2str(p_value_approx,5));

		output_result(out, lh_hotspot_pos, rh_hotspot_pos, mle_rate_hotspot,
				lh_window_pos, rh_window_pos, mle_rate_background, mle_rate_constant,
				N_sims_used, p_value, p_value_approx);

		double N_sims_per_sec = N_sims_used / sim_time;
		printLOG("  " + dbl2str_fixed(N_sims_per_sec, 2) + " simulations/sec\n");
	}

	out.close();
	time(&end);
	double running_time = difftime(end,start);
	printLOG("Run Time = " + dbl2str_fixed(running_time, 2) + " seconds\n");
	printLOG("Done.\n");
	LOG.close();
	return 0;
}
