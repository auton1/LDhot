/*
 * gpd_fit.cpp
 *
 *  Created on: 13 Aug 2010
 *      Author: auton
 */

#include "gpd_fit.h"

// Compute the log-likelihood for the Generalized Pareto Distribution
// shape = k, scale = a
double gpd_llk(const vector<double> &data, double k, double a)
{
	double llk=numeric_limits<double>::infinity();
	double ln_a = log(a);
	unsigned int n=data.size();

	if (fabs(k) > numeric_limits<double>::epsilon())
	{
		double theta = k/a;
		double sum=0;
		for (unsigned int ui=0; ui<n; ui++)
			sum += log(1-data[ui]*theta);
		llk = -n*ln_a - (1.0-1.0/k)*sum;
	}
	else
	{
		double sum=0;
		double ai = 1.0/a;
		for (unsigned int ui=0; ui<n; ui++)
			sum += data[ui]*ai;
		llk = -n*ln_a - sum;
	}
	return -llk;
}

// Estimate shape and scale by method of moments
void gpd_method_of_moments(const vector<double> &data, double &k, double &a)
{
	double mean=0.0, var=0.0;
	unsigned int n=data.size();

	if (n<2)
	{
		cout << "Warning: Not enough data to use method of moments" << endl;
		k=0; a=mean;
		return;
	}

	double sum=accumulate(data.begin(), data.end(), 0.0);
	mean = sum / n;

	double x;
	for (unsigned int ui=0; ui<n; ui++)
	{
		x = data[ui]-mean;
		var += x*x;
	}
	var /= (n-1);

	k = 0.5 * ((mean * mean / var) - 1.0);
	a = 0.5 * mean * ((mean*mean/var) + 1.0);

	double data_max = *max_element(data.begin(), data.end());

	if ((k > 0) && (data_max >= a/k))
	{	// Method of moments failed - return exp fit
		k=0; a=mean;
	}
}

void gpd_ml(const vector<double> &data, double &k, double &a)
{
	// Get initial values
	gpd_method_of_moments(data, k, a);

	// Reparameterize into k=k, theta = k/a

	unsigned int n=data.size();
	double max_data = -numeric_limits<double>::max();
	for (unsigned int ui=0; ui<n; ui++)
		max_data = max(data[ui], max_data);

	double max_theta = (1.0/max_data)-2.0*numeric_limits<double>::epsilon();

	vector<pair<double, double> > range(3);
	range[0].first = -99999.9;
	range[1].first = min(k/a, max_theta);
	range[2].first = max_theta;

	sort(range.begin(), range.end());

	double llk = 0, theta;
	double sum;

	for (unsigned int uj=0; uj<3; uj++)
	{
		theta = range[uj].first;
		sum = 0;
		for (unsigned int ui=0; ui<n; ui++)
			sum += log(1.0-data[ui]*theta);

		llk = -n - sum - n*log(-sum / (n * theta));

		range[uj].second = llk;
	}

	llk = 0;
	double old_llk = 10.0;

	double dydx;
	for (unsigned int step=0; step<10000; step++)
	{
		old_llk = llk;

		dydx = (range[2].second - range[0].second) / (range[2].first - range[0].first);
		theta = (abs((range[2].second - range[0].second)*0.75) + dydx*range[0].first)/dydx;

		if ((theta < range[0].first) || (theta > range[2].first) || (theta != theta))
			theta = 0.5*(range[0].first + range[2].first);

		sum = 0;
		for (unsigned int ui=0; ui<n; ui++)
			sum += log(1.0-data[ui]*theta);

		llk = -n - sum - n*log(-sum / (n * theta));

		if (range[0].second < range[2].second)
		{
			range[0].first = theta;
			range[0].second = llk;
		}
		else
		{
			range[2].first = theta;
			range[2].second = llk;
		}
		sort(range.begin(), range.end());

		if (fabs(old_llk - llk) < 1e-10)
			break;
	}

	sum = 0;
	for (unsigned int ui=0; ui<n; ui++)
		sum += log(1.0-data[ui]*theta);

	k = -sum / n;
	a = k / theta;
}

double gp_cdf(double x, double k, double a)
{
	double cdf;

	if (fabs(k) < numeric_limits<double>::epsilon())
		cdf = 1.0 - exp(-x/a);
	else
	{
		if (k*x/a < 1.0)
			cdf = 1.0 - exp(log(1.0 - k*x/a)/k);
		else
			cdf = 1.0;
	}

	return cdf;
}

// Test goodness of fit at 5% level
// Return true if pass, false if fail
// This goodnees-of-fit test taken from Choulakian and Stephens 2001
bool test_gof(const vector<double> &data, double k, double a)
{
	if ((k<-0.9) || (k>=0.5))
		return false;

	vector<double> k_critvec(10);
	vector<double> W_critvec(10);
	vector<double> A_critvec(10);

	k_critvec[0] = -0.9; k_critvec[1]= -0.5; k_critvec[2] = -0.2; k_critvec[3] = -0.1; k_critvec[4] = 0;
	k_critvec[5] = 0.1; k_critvec[6] = 0.2; k_critvec[7] = 0.3; k_critvec[8] = 0.4; k_critvec[9] = 0.5;
	/*// 0.005 level
	W_critvec[0] = 0.187; W_critvec[1]= 0.204; W_critvec[2] = 0.228; W_critvec[3] = 0.240; W_critvec[4] = 0.255;
	W_critvec[5] = 0.27; W_critvec[6] = 0.291; W_critvec[7] = 0.317; W_critvec[8] = 0.349; W_critvec[9] = 0.39;

	A_critvec[0] = 1.226; A_critvec[1] = 1.336; A_critvec[2] = 1.471; A_critvec[3] = 1.532;	A_critvec[4] = 1.603;
	A_critvec[5] = 1.687; A_critvec[6] = 1.788; A_critvec[7] = 1.909; A_critvec[8] = 2.058; A_critvec[9] = 2.243;
	*/
	// 0.05 level
	W_critvec[0] = 0.115; W_critvec[1] = 0.124; W_critvec[2] = 0.137; W_critvec[3] = 0.144; W_critvec[4] = 0.153;
	W_critvec[5] = 0.16;  W_critvec[6] = 0.171; W_critvec[7] = 0.184; W_critvec[8] = 0.201; W_critvec[9] = 0.222;

	A_critvec[0] = 0.771; A_critvec[1] = 0.83;  A_critvec[2] = 0.903; A_critvec[3] = 0.935;	A_critvec[4] = 0.974;
	A_critvec[5] = 1.02;  A_critvec[6] = 1.074; A_critvec[7] = 1.14;  A_critvec[8] = 1.221; A_critvec[9] = 1.321;

	unsigned int n=data.size();
	vector<double> z(n);
	for (unsigned int ui=0; ui<n; ui++)
		z[ui] = gp_cdf(data[ui], k, a);

	double W = 0;
	double t;
	for (unsigned int ui=1; ui<=n; ui++)
	{
		t = (z[ui-1] - ((2.0*ui-1.0)/(2.0*n)));
		W += t*t + 1.0/(12.0*n);
	}

	double A=0;
	for (unsigned int ui=1; ui<=n; ui++)
		A -= (2*ui-1) * (log(z[ui-1]) + log(1.0-z[n-ui]));
	A /= n;
	A -= n;

	unsigned int idx;
	for (idx=0; idx<9; idx++)
		if ((k >= k_critvec[idx]) && (k < k_critvec[idx+1]))
			break;

	double A_crit, W_crit, dydx;
	dydx = (W_critvec[idx+1] - W_critvec[idx]) / (k_critvec[idx+1] - k_critvec[idx]);
	W_crit = W_critvec[idx] + dydx*(k-k_critvec[idx]);

	dydx = (A_critvec[idx+1] - A_critvec[idx]) / (k_critvec[idx+1] - k_critvec[idx]);
	A_crit = A_critvec[idx] + dydx*(k-k_critvec[idx]);

	//cout << k << " " << A << " " << A_crit << endl;
	//cout << k << " " << W << " " << W_crit << endl;

	bool A_pass = (A < A_crit);
	bool W_pass = (W < W_crit);

	return (A_pass && W_pass);
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
			(min_idx < clk_ratio_distribution.size()-25);	// Require at least 25 points to estimate `.
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

