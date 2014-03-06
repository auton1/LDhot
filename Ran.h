/*
 * ran.h
 *
 *  Created on: May 2, 2011
 *      Author: auton
 */

#ifndef RAN_H_
#define RAN_H_

#include <ctime>

class Ran
{
public:
	Ran() : u(1), v(4101842887655102017LL), w(1)
	{
		set_seed(time(NULL));
	}
	Ran(unsigned long seed) : u(1), v(4101842887655102017LL), w(1)
	{
		set_seed(seed);
	}
	uint64_t int64()
	{
		u = u * 2862933555777941757LL + 7046029254386353087LL;
		v ^= v >> 17; v ^= v << 31; v ^= v >> 8;
		w = 4294957665U*(w & 0xffffffff) + (w >> 32);
		uint64_t x = u ^ (u << 21); x ^= x >> 35; x ^= x << 4;
		return (x + v) ^ w;
	}
	double ran_double() { return 5.42101086242752217e-20 * int64(); }
	long ran_long() { return (long)int64(); }
	long ran_int() { return (int)int64(); }
	void set_seed(unsigned long seed)
	{
		uint64_t j = seed;
		v = 4101842887655102017LL; w = 1;
		u = j ^ v; int64();
		v = u; int64();
		w = v; int64();
	}

	int rpoiss(double x)
	{
		int i=0;
		double r1=ran_double(), cump=0.0, fn;

		fn = exp(-x);
		while (cump<r1){
			cump+=fn;
			i++;
			fn*=(double) x/i;
		}
		i--;
		return i;
	}

	double rexp(double mean)
	{	// Note rate = 1/mean
		return -(log(ran_double())*mean);
	}
private:
	uint64_t u,v,w;
};

#endif /* RAN_H_ */
