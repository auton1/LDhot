/*
 * sim.h
 *
 *  Created on: Apr 18, 2010
 *      Author: auton
 */

#ifndef SIM_H_
#define SIM_H_

#include <cmath>
#include <cstdlib>
#include <cstdio>
#include <ctime>
#include <iostream>
#include <vector>
#include <cassert>

#include "Ran.h"

using namespace std;

struct node_tree
{
	int node_num;
	int site;
	int ndesc;
	double time;
	double time_above;
	double time_below;
	node_tree *a;		// Ancestral
	node_tree *d[2];	// Descendants
};


struct node_list
{
	bool nanc;
	int node_num;
	double rlen;
	double time;
	node_tree **asite;
};


class sim_control
{
public:
	sim_control(int Nchr, int Nsites, double theta, const vector<double> &rmap,
			bool infinite_sites, bool condition_on_seg_sites,
			const vector<double> &cond_locs, bool condition_on_freq, const vector<int> &cond_freq,
			unsigned int freq_cond_model) :
		nsamp(Nchr),
		len(Nsites),
		theta(theta),
		rmap(rmap),
		infinite_sites(infinite_sites),
		condition_on_seg_sites(condition_on_seg_sites),
		cond_locs(cond_locs),
		condition_on_freq(condition_on_freq),
		cond_freq(cond_freq),
		freq_cond_model(freq_cond_model)
	{
		w0=0.0;
		for (int i=1;i<nsamp;i++)
			w0+=1.0/double(i);
		//R = rmap[len-1];
		R = rmap[len-1] - rmap[0];
	};

	~sim_control()
	{
	};

	void set_rmap(vector<double> &rec_map)
	{
		assert((int)rec_map.size() == len);
		rmap = rec_map;
	}

	int nsamp;				// Sample Size
	int len;				// Length (= N_sites)
	double theta;			// Theta per site
	vector<double> rmap;	// Recombination Map
	double R;
	bool infinite_sites;	// Infinite sites model
	bool condition_on_seg_sites;	// Condition on segregating sites and frequency
	vector<double> cond_locs;		// Segregating sites, and frequency for conditioning
	bool condition_on_freq;
	vector<int> cond_freq;
	double w0;				// Watterson
	unsigned int freq_cond_model;	// Flag for controlling how SNP frequencies are conditioned etc.
private:
};

class sim
{
public:
	sim(sim_control &con) : con(con)
	{
		int i, site;
		node_tree *new_node;
		tree_size = 2*con.nsamp-1;
		tree_ptr = (node_tree ***) malloc((size_t) (con.len+1)*sizeof(node_tree **));
		for (site=1;site<=con.len;site++)
		{
			tree_ptr[site] = (node_tree **) malloc((size_t) (tree_size+1)*sizeof(node_tree *));
			for (i=1;i<=tree_size;i++)
			{
				new_node = (node_tree *) malloc((size_t) sizeof(node_tree));
				new_node->d[0]=NULL; new_node->d[1]=NULL;
				new_node->a=NULL;
				new_node->site=site;
				new_node->time_above=new_node->time_below=0.0;
				if (i<=con.nsamp)
				{
					new_node->time=0.0; new_node->node_num=i; new_node->ndesc=1;
				}
				tree_ptr[site][i] = new_node;
			}
		}

		segl.resize(con.len, 0);

		seqs.resize(con.nsamp, vector< unsigned char >(con.len,0));
		locs.resize(con.len,0);

		log_2 = log(2.0);

		lnfac_lookup_size = 0;
		unsigned int N=(con.nsamp*2)+1;
		lnfac_lookup.resize(N,0);
		for (unsigned int ui=1; ui<N; ui++)
			lnfac_lookup[ui] = lnfac(ui);
		lnfac_lookup_size = N;
		nrec=0, nco=0;
		rangen.set_seed(time(NULL));
	};
	~sim()
	{
		for (int site=1;site<=con.len;site++)
		{
			for (int i=1; i<=tree_size; i++)
				free(tree_ptr[site][i]);
			free(tree_ptr[site]);
		}
		free(tree_ptr);
	};

	void set_seed(unsigned long seed)
	{
		rangen.set_seed( seed );
	}

	sim_control con;

	vector<vector< unsigned char > > seqs;
	vector<double> locs;

	void run_sim();

private:
	int tree_size;
	node_tree ***tree_ptr;
	vector<int> segl;

	int nrec, nco;

	double log_2;
	unsigned int lnfac_lookup_size;
	vector<double> lnfac_lookup;

	void choose_time(double *t, int k, double rho);
	void count_rlen(node_list *nodel);
	node_list **recombine(int *k, double *rr, node_list **list, double t);
	node_list **coalesce(int *k, node_list **list, double t);
	double lnfac(int i);
	int add_mut(node_tree **tree, double fl);
	node_tree *add_mut_f_default(int fsim, node_tree **tree_site);
	node_tree *add_mut_f_alternate(int fsim, node_tree **tree_site);
	void seq_mut(node_tree *nm, int site, int base);
	void tree_summary();
	void make_tree();

	Ran rangen;
};

#endif /* SIM_H_ */
