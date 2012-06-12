/*
 * sfacMeanSD.cpp -- scalar field Autocorrelation Mean and SD
 *
 * input will be a list of file names each of which has data of the form
 * sigma1,<phi(
 * sigma2,<phi(x)phi(y)>
 *   .                .
 *   .                .
 *   .                .
 * output will be ??? files 
 *
 * file1 will be called prefix-sigma-meanrp-COUNT.stsdat
 * file2 will be called prefix-sigma-sdrp-COUNT.stsdat
 *
 * where COUNT is the number of .sts2p1 files that are being processed
 * prefix is a user supplied string for identifying the spacetime
 *
 */
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <functional>
#include <iomanip>
#include <iostream>
#include <limits>
#include <list>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <utility>
#include <vector>
using namespace std;


typedef unsigned long long bignum_t;
typedef unsigned short fixnum_t;
typedef pair<fixnum_t,fixnum_t> fixnum_pair_t;

typedef multimap<fixnum_t,long double> sigma_propagator_t;
typedef multimap<fixnum_t,long double>::const_iterator sp_citer_t;
typedef multimap<fixnum_t,long double>::iterator sp_iter_t;

typedef map<fixnum_t, long double> fixnum_long_double_map_t;
typedef map<fixnum_t, long double>::const_iterator ffm_citer_t;
typedef map<fixnum_t, long double>::iterator ffm_iter_t;

fixnum_t gCOUNT = 0;
fixnum_t gSIGMA_MIN = numeric_limits<fixnum_t>::max();
fixnum_t gSIGMA_MAX = numeric_limits<fixnum_t>::min();


void parseInputFile( const char* filename, 
		     sigma_propagator_t& svsrp )
{
  ifstream fin(filename);
  if(fin.is_open())
    {
      string line;
      while(getline(fin,line))
	{
	  string::size_type begIdx = 0, 
	    endIdx = line.find_first_of(" ",begIdx);
	  fixnum_t sigma = atoi(line.substr(begIdx, endIdx).c_str());
	  if (sigma > gSIGMA_MAX)
	    gSIGMA_MAX = sigma;
	  if (sigma < gSIGMA_MIN)
	    gSIGMA_MIN = sigma;
	  long double propgat = strtod(line.substr(endIdx+1).c_str(),0);
	  svsrp.insert(make_pair(sigma,propgat));
	}
      gCOUNT++;
      fin.close();
    }
}

// file1, file2
void printOutputFiles( const char* prefix, 
		       const sigma_propagator_t& svsrp )
{
  /*
   *
   * file1 will be called prefix-sigma-meanrp-COUNT.stsdat
   * file2 will be called prefix-sigma-sdrp-COUNT.stsdat
   */
  std::ostringstream ostr1, ostr2;
  ostr1 << prefix << "-sigma-meanac-smin" << gSIGMA_MIN << "-smax" <<
    gSIGMA_MAX << "-count" << gCOUNT << ".sfacdat";
  ostr2 << prefix << "-sigma-sdac-smin" << gSIGMA_MIN << "-smax" <<
    gSIGMA_MAX << "-count" << gCOUNT << ".sfacdat";

  ofstream fout1(ostr1.str().c_str());
  ofstream fout2(ostr2.str().c_str());

  fout1 << "{";
  fout2 << "{";


  for( fixnum_t sigma = gSIGMA_MIN; sigma <= gSIGMA_MAX; ++sigma )
    {
      fout1 << "{" << sigma << ","; 
      //      fout2 << "{" << sigma << ",";


      long double mean_total = 0.0;
      for( sp_citer_t iter = svsrp.lower_bound(sigma); 
	   iter != svsrp.upper_bound(sigma); ++iter )
	mean_total += iter->second;
      long double mean = mean_total / svsrp.count(sigma);
      fout1 << std::setprecision(9) << std::fixed << mean << 
	(sigma == gSIGMA_MAX ? "}": "},");
      //rpmean[sigma] = mean; // store the mean value of the retprobs

      long double sd_total = 0.0;
      for( sp_citer_t iter = svsrp.lower_bound(sigma); 
	   iter != svsrp.upper_bound(sigma); ++iter )
	{
	  long double x = iter->second;
	  sd_total += ((x - mean)*(x - mean));
	}
      long double variance = sd_total / svsrp.count(sigma);
      long double sd = sqrtf( variance );
      fout2 << std::setprecision(9) << std::fixed << (1.0/sd) << 
	(sigma == gSIGMA_MAX ? "}": ",");
      //      rpsd[sigma] = sd; // store the sd of the retprobs

    }
  fout1 << "}\n";
  fout2 << "}\n";

  fout1.close();
  fout2.close();
}

int main(int argc, char** argv)
{
  if (argc != 3)
    {
      cerr << "usage: sfacMeanSD prefix inputfilelst" << endl;
      exit(-1);
    }

  sigma_propagator_t SVsRP;
  
  //  fixnum_long_double_map_t retprobmean;
  //fixnum_long_double_map_t retprobsd;
  
  std::ifstream fin( argv[2] );
  if( fin.is_open() )
    {
      std::string infilename;
      while( std::getline( fin, infilename ) )
	{
	  parseInputFile( infilename.c_str(), SVsRP );
	}
      fin.close();
    }    
  printOutputFiles(argv[1], SVsRP );//, retprobmean, retprobsd );
}
