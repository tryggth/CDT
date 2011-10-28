/*
 * mean_sd.cpp
 */
#include <algorithm>
#include <cassert>
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

typedef unsigned long long bignum_t;
typedef unsigned short fixnum_t;

std::vector<std::vector<long double> > gALL_RET_PROBS;
fixnum_t gCOUNT = 0;
fixnum_t gSIGMA_MIN = 0;
fixnum_t gSIGMA_MAX = 0;

void parseInputFile( const char* filename, 
		     std::vector<std::vector<long double> >& allRetProbs )
{
  std::ifstream fin(filename);
  if(fin.is_open())
    {
      std::vector<long double> currRetProb(gSIGMA_MAX - gSIGMA_MIN);
      std::string line;
      size_t ndx = 0;
      while(getline(fin,line))
	{
	  std::string::size_type begIdx = 0, 
	    endIdx = line.find_first_of("\t",begIdx);
	  fixnum_t sigma = atoi(line.substr(begIdx, endIdx).c_str());
	  float retprob = strtod(line.substr(endIdx+1).c_str(),0);
	  currRetProb[ndx++] = retprob;
	}
      std::cout << "check point 2" << std::endl;
      gCOUNT++;
      allRetProbs.push_back(currRetProb);
      std::cout << "check point 3" << std::endl;
      fin.close();
    }
}

void printOutputFiles(const char* prefix, 
		      const std::vector<std::vector<long double> >& allRetProbs,
		      const char* suffix)
{
  assert(allRetProbs.size() == gCOUNT);

  std::ostringstream ostr1, ostr2;
  ostr1 << prefix << "-sigma-meanrp-smin" << gSIGMA_MIN << "-smax" <<
    gSIGMA_MAX << "-count" << gCOUNT << suffix;
  ostr2 << prefix << "-sigma-sdrp-smin" << gSIGMA_MIN << "-smax" <<
    gSIGMA_MAX << "-count" << gCOUNT << suffix;

  std::ofstream fout1(ostr1.str().c_str());
  std::ofstream fout2(ostr2.str().c_str());

  fout1 << "{";
  fout2 << "{";

  for( fixnum_t sigma = gSIGMA_MIN; sigma <= gSIGMA_MAX; ++sigma )
    {
      fout1 << "{" << sigma << ","; 

      long double mean_total = 0.0;
      
      for(size_t ndx = 0; ndx < allRetProbs.size(); ++ndx)
	mean_total += allRetProbs[ndx][sigma - gSIGMA_MIN];

      std::cout << mean_total << std::endl;
      long double mean = mean_total / gCOUNT;

      fout1 << std::setprecision(9) << std::fixed << mean << 
	(sigma == gSIGMA_MAX ? "}": "},");

      float sd_total = 0.0;
      for(size_t ndx = 0; ndx < allRetProbs.size(); ++ndx)
	{
	  long double x = allRetProbs[ndx][sigma - gSIGMA_MIN];
	  sd_total += ((x - mean)*(x - mean));
	}
      long double variance = sd_total / gCOUNT;
      long double sd = sqrt(variance);
      fout2 << std::setprecision(9) << std::fixed << (1.0/sd) << 
	(sigma == gSIGMA_MAX ? "}": ",");
    }
  fout1 << "}\n";
  fout2 << "\n";

  fout1.close();
  fout2.close();
}

int main(int argc, char** argv)
{
  if (argc != 6)
    {
      std::cerr << 
	"usage: mean_sd prefix inputfilelst sigmin sigmax suffix" << 
	std::endl;
      exit(-1);
    }

  std::ifstream fin( argv[2] );
  gSIGMA_MIN = static_cast<fixnum_t>( atoi( argv[3] ) );
  gSIGMA_MAX = static_cast<fixnum_t>( atoi( argv[4] ) );
  gALL_RET_PROBS.resize(10);
  if( fin.is_open() )
    {
      std::string infilename;
      while( std::getline( fin, infilename ) )
	{
	  parseInputFile( infilename.c_str(), gALL_RET_PROBS );
	}
      fin.close();
    }    
  printOutputFiles(argv[1], gALL_RET_PROBS, argv[5]);
}
