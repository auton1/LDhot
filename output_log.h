/*
 * log.h
 *
 *  Created on: Nov 11, 2009
 *      Author: Adam Auton
 *      ($Revision: 91 $)
 */

#ifndef LOG_H_
#define LOG_H_

#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>

using namespace std;

extern ofstream LOG;

void printLOG(string s);
void error(string err_msg, int error_code=0);
void error(string err_msg, double value1, double value2, int error_code=0);
void counted_warning(string err_msg);
void warning(string err_msg);
string int2str(int n);
string longint2str(long int n);
string dbl2str(double n, int prc);
string dbl2str_fixed(double n, int prc);

#endif /* LOG_H_ */
