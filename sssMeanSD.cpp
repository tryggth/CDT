/*
 * sssMeanSD.cpp
 *
 * input will be a list of file names each of which has data of the form
 * sigma1         retprob1
 * sigma2         retprob2
 *   .                .
 *   .                .
 *   .                .
 * sigmaN         retprobN
 *
 * output will be ??? files 
 *
 * file1 will be called prefix-sigma-meanrp-COUNT.stsdat
 * file2 will be called prefix-sigma-sdrp-COUNT.stsdat
 * file3 will be called prefix-sigma-meanrp-COUNT-FWW.stsdat
 * file4 will be called prefix-sigma-sdrp-COUNT-FWW.stsdat
 * file5 will be called prefix-sigma-Ds-COUNT.stsdat
 * file6 will be called prefix-sigma-Ds-COUNT-FWW.stsdat
 *
 * where COUNT is the number of .sss2p1 files that are being processed
 * prefix is a user supplied string for identifying the spacetime
 * FWW is the filter window width used for filtering the data
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

typedef multimap<fixnum_t,float> sigma_retprobs_t;
typedef multimap<fixnum_t,float>::const_iterator srp_citer_t;
typedef multimap<fixnum_t,float>::iterator srp_iter_t;

typedef map<fixnum_t, float> fixnum_float_map_t;
typedef map<fixnum_t, float>::const_iterator ffm_citer_t;
typedef map<fixnum_t, float>::iterator ffm_iter_t;

fixnum_t gCOUNT = 0;
fixnum_t gSIGMA_MIN = 0;
fixnum_t gSIGMA_MAX = 0;


/*
  const double CM3 = (-1.0/60.0);
  const double CM2 = (3.0/20.0);
  const double CM1 = (-3.0/4.0);
  const double CP1 = (3.0/4.0);
  const double CP2 = (-3.0/20.0);
  const double CP3 = (1.0/60.0);
*/
const double CM4 = (1.0/280.0);
const double CM3 = (-4.0/105.0);
const double CM2 = (1.0/5.0);
const double CM1 = (-4.0/5.0);
const double CP1 = (4.0/5.0);
const double CP2 = (-1.0/5.0);
const double CP3 = (4.0/105.0);
const double CP4 = (-1.0/280.0);

void parseInputFile( const char* filename, 
		     sigma_retprobs_t& svsrp, 
		     sigma_retprobs_t& sigmaVsLogRetProb )
{
  ifstream fin(filename);
  if(fin.is_open())
    {
      string line;
      while(getline(fin,line))
	{
	  string::size_type begIdx = 0, 
	    endIdx = line.find_first_of("\t",begIdx);
	  fixnum_t sigma = atoi(line.substr(begIdx, endIdx).c_str());
	  float retprob = strtod(line.substr(endIdx+1).c_str(),0);
	  float logrp = logf(retprob);
	  svsrp.insert(make_pair(sigma,retprob));
	  sigmaVsLogRetProb.insert(make_pair(sigma,retprob));
	}
      gCOUNT++;
      fin.close();
    }
}

// file1, file2, file5
void printOutputFiles( const char* prefix, 
		       const sigma_retprobs_t& svsrp,
		       const sigma_retprobs_t& sigmaVsLogRetProb,
		       fixnum_float_map_t& rpmean, 
		       fixnum_float_map_t& rpsd,
		       fixnum_float_map_t& logmeanrp, 
		       fixnum_float_map_t& logsdrp )
{
  /*
   *
   * file1 will be called prefix-sigma-meanrp-COUNT.stsdat
   * file2 will be called prefix-sigma-sdrp-COUNT.stsdat
   * file5 will be called prefix-sigma-Ds-COUNT.stsdat
   */
  std::ostringstream ostr1, ostr2, ostr5;
  ostr1 << prefix << "-sigma-meanrp" << gCOUNT << ".sssdat";
  ostr2 << prefix << "-sigma-sdrp-" << gCOUNT << ".sssdat";
  ostr5 << prefix << "-sigma-Ds-" << gCOUNT << ".sssdat";

  ofstream fout1(ostr1.str().c_str());
  ofstream fout2(ostr2.str().c_str());
  ofstream fout5(ostr5.str().c_str());

  fout1 << "{";
  fout2 << "{";


  for( fixnum_t sigma = gSIGMA_MIN; sigma <= gSIGMA_MAX; ++sigma )
    {
      fout1 << "{" << sigma << ","; 
      fout2 << "{" << sigma << ",";


      float mean_total = 0.0;
      for( srp_citer_t iter = svsrp.lower_bound(sigma); 
	   iter != svsrp.upper_bound(sigma); ++iter )
	mean_total += iter->second;
      float mean = mean_total / svsrp.count(sigma);
      fout1 << mean << (sigma == gSIGMA_MAX ? "}": "},");
      rpmean[sigma] = mean; // store the mean value of the retprobs

      float sd_total = 0.0;
      for( srp_citer_t iter = svsrp.lower_bound(sigma); 
	   iter != svsrp.upper_bound(sigma); ++iter )
	{
	  sd_total += ( (iter->second - mean)*(iter->second - mean) );
	}
      float sd = sqrtf( sd_total / svsrp.count( sigma ) );
      fout2 << sd << (sigma == gSIGMA_MAX ? "}": "},");
      rpsd[sigma] = sd; // store the sd of the retprobs

    }
  fout1 << "}\n";
  fout2 << "}\n";

  fout5 << "{";
  for( fixnum_t sigma = gSIGMA_MIN + 4; sigma <= gSIGMA_MAX - 4; ++sigma )
    {
      fout5 << "{" << sigma << ",";
      double deriv = 
	-2.0 * sigma * ( CM4*log(rpmean[sigma-4]) +
			 CM3*log(rpmean[sigma-3]) +
			 CM2*log(rpmean[sigma-2]) +
			 CM1*log(rpmean[sigma-1]) +
			 CP1*log(rpmean[sigma+1]) +
			 CP2*log(rpmean[sigma+2]) +
			 CP3*log(rpmean[sigma+3]) +
			 CP4*log(rpmean[sigma+4]) );
      fout5 << deriv << (sigma == (gSIGMA_MAX-3) ? "}" : "},");
    }
  fout5 << "}\n";

  fout1.close();
  fout2.close();
  fout5.close();
}

// filtered window width has to be an odd number
void printFilteredOutputFiles( const char* prefix,
			       fixnum_t fww, // fww = filter window width
			       const fixnum_float_map_t& rpmean, 
			       const fixnum_float_map_t& rpsd,
			       const fixnum_float_map_t& logmeanrp, 
			       const fixnum_float_map_t& logsdrp,
			       const fixnum_float_map_t& meanlogrp, 
			       const fixnum_float_map_t& sdlogrp)
{
  // gSIGMA_MAX - gSIGMA_MIN + 1 is the unfiltered data vector size
  // the filtered data window size will be the unfiltered size - (fww - 1)
  // which is equal to gSIGMA_MAX - gSIGMA_MIN - fww + 2
  float *filteredMeanRetProbs = new float[gSIGMA_MAX-gSIGMA_MIN-fww+2];
  float *filteredMeanLogRetProbs = new float[gSIGMA_MAX-gSIGMA_MIN-fww+2];
  fixnum_t *filteredSigmaVals = new fixnum_t[gSIGMA_MAX-gSIGMA_MIN-fww+2];
  fixnum_t offset = (fww - 1)/2;
  for( fixnum_t sigma = gSIGMA_MIN + offset; 
       sigma <= gSIGMA_MAX - offset; ++sigma )
    {
      filteredSigmaVals[ sigma - gSIGMA_MIN - offset ] = sigma;
      float filtered_retprob_mean_value = 0.0, filtered_meanlog_retprob = 0.0;
      for ( fixnum_t pos = sigma - offset; pos <= sigma + offset; ++pos )
	{
	  filtered_retprob_mean_value += rpmean.find(pos)->second;
	  filtered_meanlog_retprob += meanlogrp.find(pos)->second;
	}
      filtered_retprob_mean_value /= fww;
      filtered_meanlog_retprob /= fww;
      filteredMeanRetProbs[sigma-gSIGMA_MIN-offset]=filtered_retprob_mean_value;
      filteredMeanLogRetProbs[sigma-gSIGMA_MIN-offset]=filtered_meanlog_retprob;
    }

  std::ostringstream ostr3, ostr4, ostr6;
  ostr3 << prefix << "-sigma-meanrp-" << gCOUNT << "-" << fww << ".sssdat";
  ostr4 << prefix << "-sigma-sdrp-" << gCOUNT << "-" << fww << ".sssdat";
  ostr6 << prefix << "-sigma-Ds-" << gCOUNT << "-" << fww << ".sssdat";
  
  ofstream fout3(ostr3.str().c_str());
  //ofstream fout4(ostr4.str().c_str());
  ofstream fout6(ostr6.str().c_str());

  fout3 << "{";
  //fout4 << "{";

  for( fixnum_t sigma = 0; sigma < gSIGMA_MAX-gSIGMA_MIN-fww+2; ++sigma )
    {
      fout3 << "{" << filteredSigmaVals[sigma] << ",";
      //fout4 << "{" << filteredSigmaVals[sigma] << ",";
      fout3 << filteredMeanRetProbs[sigma] << 
	(sigma == (gSIGMA_MAX-gSIGMA_MIN-fww+1) ? "}" : "},");
    }
  fout3 << "}" << endl;

  fout6 << "{";
  for( fixnum_t sigma = 4; sigma < gSIGMA_MAX-gSIGMA_MIN-fww-2; ++sigma )
    {
      fout6 << "{" << filteredSigmaVals[sigma] << ",";

      double deriv = 
	-2.0 * sigma * ( CM4*log(filteredMeanRetProbs[sigma-4]) +
			 CM3*log(filteredMeanRetProbs[sigma-3]) +
			 CM2*log(filteredMeanRetProbs[sigma-2]) +
			 CM1*log(filteredMeanRetProbs[sigma-1]) +
			 CP1*log(filteredMeanRetProbs[sigma+1]) +
			 CP2*log(filteredMeanRetProbs[sigma+2]) +
			 CP3*log(filteredMeanRetProbs[sigma+3]) +
			 CP4*log(filteredMeanRetProbs[sigma+4]) );

      fout6 << deriv << (sigma == (gSIGMA_MAX-gSIGMA_MIN-fww-1) ? "}" : "},");
    }

  fout3.close();
  fout6.close();
}

int main(int argc, char** argv)
{
  if (argc != 6)
    {
      cerr << "usage: stsMeanSD prefix inputfilelst sigmin sigmax fww" << endl;
      exit(-1);
    }

  sigma_retprobs_t SVsRP, sigmaVsLogRetProb;
  
  fixnum_float_map_t retprobmean;
  fixnum_float_map_t retprobsd;
  
  fixnum_float_map_t logmeanrp; // log of mean of return probabilities
  fixnum_float_map_t logsdrp; // log of standard deviation of ret probs

  fixnum_float_map_t meanlogrp; // mean of log of return probabilities
  fixnum_float_map_t sdlogrp; // standard deviation of log of ret probs

  std::ifstream fin( argv[2] );
  if( fin.is_open() )
    {
      std::string infilename;
      while( std::getline( fin, infilename ) )
	{
	  parseInputFile( infilename.c_str(), SVsRP, sigmaVsLogRetProb );
	}
      fin.close();
    }    
  gSIGMA_MIN = static_cast<fixnum_t>( atoi( argv[3] ) );
  gSIGMA_MAX = static_cast<fixnum_t>( atoi( argv[4] ) );
  fixnum_t fww = static_cast<fixnum_t>( atoi( argv[5] ) );

  printOutputFiles(argv[1], SVsRP, sigmaVsLogRetProb, 
		   retprobmean, retprobsd,
		   logmeanrp, logsdrp );
  
  printFilteredOutputFiles( argv[1], fww, 
			    retprobmean, retprobsd,
			    logmeanrp, logsdrp, 
			    meanlogrp, sdlogrp );
}
