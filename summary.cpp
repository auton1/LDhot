/*
 * summarise.cpp
 *
 *  Created on: Mar 6, 2014
 *      Author: aauton
 */


#include "summary.h"


void print_help()
{
	cout << endl << "My LDhot Summarize" << endl;
	cout << "\u00A9 Adam Auton 2014" << endl << endl;

	cout << "Required Parameters: " << endl;
	cout << "--res <filename>" << endl;
	cout << "--hot <filename>" << endl;
	cout << endl;
	cout << "Other Parameters: " << endl;
	cout << "--sig <double>" << endl;
	cout << "--frac <double> " << endl;
	cout << endl;
	exit(0);
}

// Small program to combine the hotspot results, and set the widths on the basis of the
// rate estimates.
int main(int argc, char *argv[])
{
	time_t start,end;
	time(&start);

	// 1. Read Parameters
	int i=1;
	string in_str;
	string res_filename="", hot_filename="", output_prefix = "out";
	double sig = 0.001;
	double frac = 0.1;
	while (i<argc)
	{
		in_str = argv[i];
		if (in_str ==  "--res") { res_filename = string(argv[i+1]); i++; }
		else if (in_str ==  "--hot") { hot_filename = string(argv[i+1]); i++; }
		else if (in_str ==  "--out") { output_prefix = string(argv[i+1]); i++; }
		else if (in_str ==  "--sig") { sig = atof(argv[i+1]); i++; }
		else if (in_str ==  "--frac") { frac = atof(argv[i+1]); i++; }
		else if ((in_str == "-?") || (in_str == "-h")) print_help();
		else
			error("Unknown option: " + string(in_str), 0);
		i++;
	}

	if ((res_filename == "") || (hot_filename == ""))
	{
		print_help();
	}

	LOG.open((output_prefix + ".summary.log").c_str());

	printLOG("\nmyLDhot summary\n");
	printLOG("(C) Adam Auton 2014\n\n");
	printLOG("Parameters as interpreted:\n");
	printLOG("\t--res " + res_filename + "\n");
	printLOG("\t--hot " + output_prefix + "\n");
	printLOG("\t--out " + output_prefix + "\n");
	printLOG("\t--sig " + dbl2str(sig, 3) + "\n");
	printLOG("\t--frac " + dbl2str(frac, 3) + "\n");
	printLOG("\n");

	vector<double> rmap;
	vector<double> rmap_pos;
	read_res_file(res_filename, rmap, rmap_pos);

	vector<double> rmap_rate(rmap.size()-1, 0);
	for (unsigned int ui = 0; ui < rmap_rate.size(); ui++)
		rmap_rate[ui] = (rmap[ui+1] - rmap[ui]) / (rmap_pos[ui+1] - rmap_pos[ui]);

	vector<vector<double> > hotspots;	// hotStart, hotEnd, p

	printLOG("Reading " + hot_filename + "\n");
	ifstream in(hot_filename.c_str());
	if (!in.good())
		error("Could not open hot file: " + hot_filename);

	stringstream ss;
	string line;
	getline(in, line);	// Header

	double pos1, pos2, MLE_hotspot_rate, BgStart, BgEnd, MLE_bg_rate, MLE_constant_rate, N_used_sims, P_ecdf, P_tail_approx;
	double p;
	while(!in.eof())
	{
		getline(in,line);
		if (line.size() == 0)
			continue;

		ss.clear(); ss.str(line);
		ss >> pos1 >> pos2 >> MLE_hotspot_rate >> BgStart >> BgEnd >> MLE_bg_rate >> MLE_constant_rate >> N_used_sims >> P_ecdf >> P_tail_approx;

		vector<double> hotspot;
		hotspot.push_back(pos1);
		hotspot.push_back(pos2);
		p = min(P_ecdf, P_tail_approx);
		hotspot.push_back(P_ecdf);
		hotspots.push_back(hotspot);
	}

	printLOG("Read " + int2str(hotspots.size()) + " hotspot windows.\n");

	// Combine adjacent windows
	bool change=true;
	while (change == true)
	{
		change = false;
		for (unsigned int ui=0; ui<hotspots.size(); ui++)
		{
			pos1 = hotspots[ui][0];
			pos2 = hotspots[ui][1];
			p = hotspots[ui][2];
			for (unsigned int uj=0; uj<hotspots.size(); uj++)
			{
				if (ui == uj)
					continue;
				if (((pos1 >= hotspots[uj][0]) && (pos1 <= hotspots[uj][1])) ||
						((pos2 >= hotspots[uj][0]) && (pos2 <= hotspots[uj][1])))
				{
					hotspots[ui][0] = min(pos1, hotspots[uj][0]);
					hotspots[ui][1] = max(pos2, hotspots[uj][1]);
					hotspots[ui][2] = min(p, hotspots[uj][2]);
					hotspots[uj] = hotspots[ui];
					change = true;
				}
			}
		}

		sort(hotspots.begin(), hotspots.end());
		std::vector< vector<double> >::iterator it;
		it = std::unique (hotspots.begin(), hotspots.end());
		hotspots.resize( std::distance(hotspots.begin(),it) );
	}

	// Remove hotspots that aren't significant
	vector<vector<double> > new_hotspots;
	new_hotspots.reserve(hotspots.size());
	for (unsigned int ui=0; ui<hotspots.size(); ui++)
	{
		if (hotspots[ui][2] < sig)
			new_hotspots.push_back(hotspots[ui]);
	}
	hotspots = new_hotspots;

	// Find the width of each hotspot
	for (unsigned int ui=0; ui<hotspots.size(); ui++)
	{
		// Find peak rate in each hotspot
		double max_rate = -1;
		for (unsigned int uj=0; uj<(rmap_pos.size()-1); uj++)
		{
			if ((rmap_pos[uj] >= hotspots[ui][0]) && (rmap_pos[uj] < hotspots[ui][1]))
			{
				max_rate = max(max_rate, rmap_rate[uj]);
			}
		}
		// Find the positions where the rate drops to a fraction of the peak rate
		unsigned int lh_idx = 0, rh_idx = 0;
		for (unsigned int uj=0; uj<(rmap_pos.size()-1); uj++)
		{
			if ((rmap_rate[uj] <= (max_rate * frac)) && (rmap_pos[uj] <= hotspots[ui][0]))
				lh_idx = uj;
			if ((rmap_rate[uj] <= (max_rate * frac)) && (rmap_pos[uj] >= hotspots[ui][1]))
			{
				rh_idx = uj;
				break;
			}
		}
		lh_idx++;
		hotspots[ui][0] = min(rmap_pos[lh_idx], hotspots[ui][0]);
		hotspots[ui][1] = max(rmap_pos[rh_idx], hotspots[ui][1]);
	}

	// Combine adjacent windows
	change=true;
	while (change == true)
	{
		change = false;
		for (unsigned int ui=0; ui<hotspots.size(); ui++)
		{
			pos1 = hotspots[ui][0];
			pos2 = hotspots[ui][1];
			p = hotspots[ui][2];
			for (unsigned int uj=0; uj<hotspots.size(); uj++)
			{
				if (ui == uj)
					continue;
				if (((pos1 >= hotspots[uj][0]) && (pos1 <= hotspots[uj][1])) ||
						((pos2 >= hotspots[uj][0]) && (pos2 <= hotspots[uj][1])))
				{
					hotspots[ui][0] = min(pos1, hotspots[uj][0]);
					hotspots[ui][1] = max(pos2, hotspots[uj][1]);
					hotspots[ui][2] = min(p, hotspots[uj][2]);
					hotspots[uj] = hotspots[ui];
					change = true;
				}
			}
		}

		sort(hotspots.begin(), hotspots.end());
		std::vector< vector<double> >::iterator it;
		it = std::unique (hotspots.begin(), hotspots.end());
		hotspots.resize( std::distance(hotspots.begin(),it) );
	}

	printLOG("Writing " + int2str(hotspots.size()) + " hotspots.\n");

	ofstream out((output_prefix + ".hot_summary.txt").c_str());
	out << "#hotStart\thotEnd\tp" << endl;
	for (unsigned int ui=0; ui<hotspots.size(); ui++)
		out << hotspots[ui][0] << "\t" << hotspots[ui][1] << "\t" << hotspots[ui][2] << endl;
	out.close();

	time(&end);
	double running_time = difftime(end,start);
	printLOG("Run Time = " + dbl2str_fixed(running_time, 2) + " seconds\n");
	printLOG("Done.\n");
	LOG.close();
	return 0;
}
