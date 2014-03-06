/*
 * sim.cpp
 *
 *  Created on: Apr 18, 2010
 *      Author: auton
 */

#include "sim.h"

void sim::run_sim()
{
	nrec=0, nco=0;
	make_tree();

	tree_summary();
}

/**************************************************/
/*Routine to choose time for next coalescent event*/
/**************************************************/

void sim::choose_time(double *t, int k, double rho)
{
	double u=-log(rangen.ran_double());

	/*Constant population size*/
	(*t) += (double) 2*u/(k*rho+k*(k-1));
}

/***************************************/
/*Calculate potential for recombination*/
/***************************************/

void sim::count_rlen(node_list *nodel)
{
	int i;

	nodel->rlen=0.0; nodel->nanc=true;
	for (i=1; i<=con.len; i++)
		if (nodel->asite[i])
			break;
	if (i>con.len)
	{ /*If no longer any ancestral material*/
		nodel->rlen = 0; nodel->nanc=false;
		return;
	}
	nodel->rlen -= (double) con.rmap[i-1];
	for (i=con.len; i>0; i--)
		if (nodel->asite[i])
			break;
	nodel->rlen += (double) con.rmap[i-1];
}


/***********************/
/*Genereate recombinant*/
/***********************/

node_list ** sim::recombine(int *k, double *rr, node_list **list, double t)
{
	int i,j,j0, upper=con.len;
	double cump, ri;
	node_list *new_node_list;
	list = (node_list **) realloc(list, (size_t)((*k+1)*sizeof(node_list *))); /*Resize: use k+1th pos to fill in*/
	new_node_list = (node_list *) malloc((size_t) sizeof(node_list));
	new_node_list->asite = (node_tree **) malloc((size_t) (con.len+1)*sizeof(node_tree *));
	list[(*k)]=new_node_list;

	/*Choose lineage to recombine*/

	ri = (double) rangen.ran_double()*(*rr);
	cump=0;
	for (i=1; cump<ri; i++)
		cump += list[i]->rlen;
	i--;
	ri = (double) rangen.ran_double()*(list[i]->rlen);  /*Choose point along sequence to recombine*/
	for (j0=1;(list[i]->asite[j0])==NULL;j0++){}; /*Find beginning of ancestral material*/
	cump=0;
	for (j=j0;cump<ri;j++)
		cump=con.rmap[j-1]-con.rmap[j0-1]; /*Choose last point on LHS to come from parent 1*/
	j0=j-2;

	for (j=1; j<=j0; j++)
		list[(*k)]->asite[j]=NULL;	/*Copy LHS to parent*/

	/*Choose if crossing-over or gene conversion event and define upper point to copy until*/

	for (;j<=upper; j++)
	{
		list[(*k)]->asite[j]=list[i]->asite[j];
		list[i]->asite[j]=NULL;
	}
	for (;j<=con.len;j++)
		list[(*k)]->asite[j]=NULL;

	list[i]->time=list[(*k)]->time=t;

	count_rlen(list[i]);
	count_rlen(list[(*k)]);

	return list;
}


/*******************************/
/* Generate coalescent         */
/*******************************/

node_list ** sim::coalesce(int *k, node_list **list, double t)
{
	int i, j, pos, node_pos;

	i=int(1+((*k)+1)*rangen.ran_double()); /*Choose first lineage to coalesce*/
	j=i;
	while(j==i)
		j = int(1 +((*k)+1)*rangen.ran_double()); /*Choose second lineage to coalesce*/

	/*Find all sites for which coalescent event has occured and update tree accordingly*/
	for (pos=1;pos<=con.len;pos++)
	{
		if (list[i]->asite[pos] && list[j]->asite[pos])
		{/*Coalescent event at site*/
			segl[pos-1]--;	// segl[pos][1]--;
			node_pos = 2*con.nsamp-segl[pos-1]; /*Number in tree_ptr of next unused pointer for pos*/
			tree_ptr[pos][node_pos]->d[0]=list[i]->asite[pos];
			tree_ptr[pos][node_pos]->d[1]=list[j]->asite[pos];
			(tree_ptr[pos][node_pos]->d[0])->a=tree_ptr[pos][node_pos];
			(tree_ptr[pos][node_pos]->d[0])->time_above=t-(tree_ptr[pos][node_pos]->d[0])->time;
			(tree_ptr[pos][node_pos]->d[1])->a=tree_ptr[pos][node_pos];
			(tree_ptr[pos][node_pos]->d[1])->time_above=t-(tree_ptr[pos][node_pos]->d[1])->time;
			tree_ptr[pos][node_pos]->time = t;
			tree_ptr[pos][node_pos]->time_below = (double) (tree_ptr[pos][node_pos]->d[0])->time_below+\
				(tree_ptr[pos][node_pos]->d[0])->time_above+(tree_ptr[pos][node_pos]->d[1])->time_below+\
				(tree_ptr[pos][node_pos]->d[1])->time_above;
			tree_ptr[pos][node_pos]->ndesc = (tree_ptr[pos][node_pos]->d[0])->ndesc+(tree_ptr[pos][node_pos]->d[1])->ndesc;
			tree_ptr[pos][node_pos]->node_num = node_pos;
			if (segl[pos-1]>1)
				list[j]->asite[pos]=tree_ptr[pos][node_pos];
			else
				list[j]->asite[pos]=NULL;
		}
		else
			list[j]->asite[pos] = (node_tree *) (list[i]->asite[pos]?list[i]->asite[pos]:list[j]->asite[pos]); /*Just copy existing*/
	}
	count_rlen(list[j]);
	list[j]->time=t;

/*Reallocate list - check to see if can remove new node because all sites found MRCA*/
	if (!(list[j]->nanc))
	{
		free(list[i]->asite);
		free(list[i]);
		free(list[j]->asite);
		free(list[j]);
		if (j==(*k))
		{
			list[i]=list[(*k)+1];
			list[j]=list[(*k)];
		}
		else {
			list[i]=list[(*k)];
			list[j]=list[(*k)+1];
		}
		(*k)--;
	}
	else {
		free(list[i]->asite);
		free(list[i]);
		list[i]=list[(*k)+1];/*Replace i with last in (previous) list*/
	}
	list = (node_list **) realloc(list, (size_t) (*k+1)*sizeof(node_list *));
	return list;

}

double sim::lnfac(int i)
{
	if ((unsigned)i < lnfac_lookup_size)
		return lnfac_lookup[i];

	int j;
	double cp=0.0;
	for (j=2;j<=i;j++)
		cp += log((double) j);
	return cp;
}



/*********************************************/
/*Routine to place mutations on the genealogy*/
/*********************************************/
// Return the index of the branch to be mutated?
int sim::add_mut(node_tree **tree, double fl)
{
	int i;

	//for (i=1;i<=tree_size &  *fl > 0.0;i++) *fl -= tree[i]->time_above;
	for (i=1;(i<=tree_size) && (fl > 0.0);i++)
		fl -= tree[i]->time_above;
	return (--i);
}

/******************************************************************************/
/*Routine to place mutations on the genealogy conditioning on allele frequency*/
/*fsim is an integer with the MAF                                             */
/******************************************************************************/
node_tree * sim::add_mut_f_default(int fsim, node_tree **tree_site)
{
	int i, maf;
	double cump, r1, part[3], mx=0.0;

	vector<double> branch_prob(tree_size, 0);

	fsim = min(fsim, con.nsamp-fsim);

	// Heuristic for conditioning on frequency.
	//1.  Simulate an ARG for n samples.
	//2.  At the sites where you want to put mutations, list the number of descendants (j) of each internal branch.
	//3.  For each branch, weight by the branch length x the probability of seeing j in n given
	//		an earlier sample that had i in n (where i is the number you want to condition on).
	//4.  Sample the branch on which to place the mutation with probabilities given by 3.

	//Weight branches by marginal Pr(j | i)
	// Gil's original implementation
	/*
	double cons = log_2+2*lnfac(con.nsamp)-lnfac(2*con.nsamp)-lnfac(fsim-1)-lnfac(con.nsamp-fsim-1);
	for (i=1;i<tree_size;i++)
	{
		maf = min(tree_site[i]->ndesc, con.nsamp-tree_site[i]->ndesc);	// MAF of a mutation placed at this branch
		part[0] = lnfac(con.nsamp+fsim-maf-1)+lnfac(con.nsamp-fsim+maf-1);
		part[1] = lnfac(fsim+maf-1)+lnfac(2*con.nsamp-fsim-maf-1);
		part[2] = lnfac(maf)+lnfac(con.nsamp-maf);
		branch_prob[i-1] = tree_site[i]->time_above*(exp(cons+part[0]-part[2])+exp(cons+part[1]-part[2]));
		if (branch_prob[i-1]>mx)
			mx=branch_prob[i-1];
	}
	*/

	// Implementation assuming beta-binomial model. Seems to be equivalent of above.
	// Basic model assumes p ~ Beta(alpha=1, alpha=1); i ~ Bin(n, p) -> j ~ BetaBin(n, 1+i, 1+n-i).
	// but taking account of the folding that arises from using the minor allele frequency.
	double cons = log_2+lnfac(con.nsamp)+lnfac(con.nsamp+1)-lnfac(fsim)-lnfac(con.nsamp-fsim)-lnfac((2*con.nsamp)+1);

	//Weight branches by marginal Pr(j | i)
	for (i=1;i<tree_size;i++)
	{
		maf = min(tree_site[i]->ndesc, con.nsamp-tree_site[i]->ndesc);	// MAF of a mutation placed at this branch
		part[0] = lnfac(fsim+maf)+lnfac((2*con.nsamp)-fsim-maf);
		part[1] = lnfac(con.nsamp-fsim+maf)+lnfac(con.nsamp+fsim-maf);
		part[2] = lnfac(maf)+lnfac(con.nsamp-maf);
		branch_prob[i-1] = tree_site[i]->time_above*(exp(cons+part[0]-part[2])+exp(cons+part[1]-part[2]));
		if (branch_prob[i-1]>mx)
			mx=branch_prob[i-1];
	}

	//Normalise to workable numbers
	cump=0.0;
	for (i=0;i<(tree_size-1);i++)
	{
		branch_prob[i] /= mx;
		cump += branch_prob[i];
	}
	// Choose a branch
	r1 = rangen.ran_double()*cump;
	cump=0.0;
	for (i=0;(i<tree_size-1) && (cump<r1);i++)
		cump += branch_prob[i];

	return tree_site[i];
}


// Replacement of add_mut_f
// If you *can* match for frequency, then do so. Otherwise, don't try.
node_tree * sim::add_mut_f_alternate(int fsim, node_tree **tree_site)
{
	int i, maf;
	double cump, r1, mx=0.0;

	vector<double> branch_prob(tree_size, 0);

	fsim = min(fsim, con.nsamp-fsim);

	// Weight by branch length
	for (i=1;i<tree_size;i++)
	{
		maf = min(tree_site[i]->ndesc, con.nsamp-tree_site[i]->ndesc);	// MAF of a mutation placed at this branch
		if (maf == fsim)
		{
			branch_prob[i-1]= tree_site[i]->time_above;
			if (branch_prob[i-1]>mx)
				mx=branch_prob[i-1];
		}
	}

	if (mx != 0)
	{	// Could find a branch for matching allele frequency.
		// Normalise to workable numbers
		cump=0.0;
		for (i=0;i<(tree_size-1);i++)
		{
			branch_prob[i] /= mx;
			cump += branch_prob[i];
		}
		// Choose a branch
		r1 = rangen.ran_double()*cump;
		cump=0.0;
		for (i=0;(i<tree_size-1) && (cump<r1);i++)
			cump += branch_prob[i];

		return tree_site[i];
	}
	else
	{	// Could not find a branch for matching allele frequency, so select branch without freq cond.
		double tree_len = tree_site[tree_size]->time_below;
		double tt = (double) tree_len*rangen.ran_double();
		return tree_site[add_mut(tree_site, tt)];
	}
}

/**************************************************/
/*Routine to mutate sequences at tips of genealogy*/
/**************************************************/

void sim::seq_mut(node_tree *nm, int site, int base)
{
	if ((nm->d[0]==NULL) && (nm->d[1]==NULL))
	{ /*terminal*/
		seqs[nm->node_num-1][site-1] = base;
	}
	else
	{
		seq_mut(nm->d[0], site, base);
		seq_mut(nm->d[1], site, base);
	}
}



/******************************************************************/
/*Add mutations to tree and calculate summary statistics of sample*/
/******************************************************************/

void sim::tree_summary()
{
	int site, i, j;
	double tree_len=0, tt;
	node_tree *node_mut;

	for (i=0; i<con.nsamp; i++)
		for (j=0; j<con.len; j++)
			seqs[i][j]=0;

//	Add mutations to tree
	for (site=1; site<=con.len; site++)
	{
		tree_len = tree_ptr[site][tree_size]->time_below;
		if (con.condition_on_seg_sites)
		{	// Condition on seg site locations
			locs[site-1] = con.cond_locs[site-1];
			if ((con.condition_on_freq == false) || (con.cond_freq[site-1]==0))
			{	// Don't condition on frequency
				tt = (double) tree_len*rangen.ran_double();
				node_mut = tree_ptr[site][add_mut(tree_ptr[site], tt)];
			}
			else
			{	// Condition on the frequency in con.cond_freq[site-1]
				if (con.freq_cond_model == 0)
					node_mut = add_mut_f_default(con.cond_freq[site-1], tree_ptr[site]);
				else // if (con.freq_cond_model == 1)
					node_mut = add_mut_f_alternate(con.cond_freq[site-1], tree_ptr[site]);
			}
			seq_mut(node_mut, site, 1);
		}
		else
		{	// Drop mutations at random
			int nmuts=rangen.rpoiss((con.theta)*tree_len/2);
			if ((con.infinite_sites)&&(nmuts))
				nmuts=1;
			for (i=0; i<nmuts; i++)
			{
				tt = (double) ((double) tree_len*rangen.ran_double());
				j = add_mut(tree_ptr[site], tt);
				node_mut = tree_ptr[site][j];
				seq_mut(node_mut, site, i+1);
			}
		}
	}
}

void sim::make_tree()
{
	int k, i, j;
	nrec=0, nco=0;
	double rr, avR, t=0, told;
	node_list **list, *new_node_list;
	list = (node_list **) malloc((size_t)((con.nsamp+1)*sizeof(node_list *)));
	for (i=1;i<=con.nsamp;i++)
	{
		new_node_list = (node_list *) malloc((size_t) sizeof(node_list));
		new_node_list->asite = (node_tree **) malloc((size_t) (con.len+1)*sizeof(node_tree *));
		for (j=1;j<=con.len;j++)
			new_node_list->asite[j]=tree_ptr[j][i]; /*Start pointers at base of trees at each site*/
		new_node_list->nanc=true;
		new_node_list->rlen = con.R;
		new_node_list->time=0.0;
		list[i]=new_node_list;
	}
	for (i=0;i<con.len;i++)
	{
		segl[i]=con.nsamp;
	}

	k=con.nsamp;	/*Number of lineages active in list*/
	while (k>1)
	{
	/*Count lineages which can contribute to rec events*/
		rr=0;
		for(i=1; i<=k; i++)
			rr += list[i]->rlen;

		avR = (double) rr/k;
	/*Choose time for next event*/
		told=t;
		choose_time(&t, k, avR);
	/*Choose type of next event: update list and tree accordingly*/
		if (rangen.ran_double() < (double) avR/(avR+(k-1)))
		{  /*Event is a recombination*/
			k++; nrec++;
		//	cout << "rec" << endl;
			list = recombine(&k, &rr, list, t);
		}
		else
		{ /*Event is a coalescent*/
			nco++; k--;
			//cout << "coal" << endl;
			list = coalesce(&k,list,t);
		}
	}

	free(list);
}


