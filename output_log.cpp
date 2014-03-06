/*
 * log.cpp
 *
 *  Created on: Nov 11, 2009
 *      Author: Adam Auton
 *      ($Revision: 66 $)
 */

#include "output_log.h"

void printLOG(string s)
{
	cout << s; cout.flush();
	LOG << s; LOG.flush();
}

void error(string err_msg, int error_code)
{
	printLOG("Error:" + err_msg + "\n");
	exit(error_code);
}


void error(string err_msg, double value1, double value2, int error_code)
{
	printLOG("Error:" + err_msg + "\n");
	stringstream ss;
	ss << "Value1=" << value1 << " Value2=" << value2 << endl;
	printLOG(ss.str());
	exit(error_code);
}


void counted_warning(string err_msg)
{
	static unsigned int warning_count = 0;
	printLOG(err_msg + "\n");
	warning_count++;
	if (warning_count > 1000)
		error("Stopping at 1000 entry-level warnings", 10);
}

void warning(string err_msg)
{
	printLOG(err_msg + "\n");
}

string int2str(int n)
{
  std::ostringstream s2( std::stringstream::out );
  s2 << n;
  return s2.str();
}

string longint2str(long int n)
{
  std::ostringstream s2( std::stringstream::out );
  s2 << n;
  return s2.str();
}

string dbl2str(double n, int prc)
{
  std::ostringstream s2;
  if ( prc > 0 )
    s2.precision(prc);
  s2 << n;
  return s2.str();
}

string dbl2str_fixed(double n, int prc)
{
  std::ostringstream s2;
  s2 << setiosflags( ios::fixed );
  if ( prc > 0 )
    s2.precision(prc);
  s2 << n;
  return s2.str();
}
