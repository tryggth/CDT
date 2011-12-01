/*
 * stsDiffusion3p1SI.cpp
 */
#include "stsDiffusion3p1SI.hpp"
#include <algorithm>
#include <ctime>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <list>
// boost libraries
#include <boost/generator_iterator.hpp>
#include <boost/random/linear_congruential.hpp>
#include <boost/random/uniform_int.hpp>
#include <boost/random/uniform_real.hpp>
#include <boost/random/variate_generator.hpp>
#include <boost/random/mersenne_twister.hpp>

typedef boost::mt19937 base_generator_type;

bignum_t pick_random_simplex(const std::vector<bignum_t>& sxs )
{
  static base_generator_type generator(static_cast<unsigned int>(std::time(0)));
  typedef boost::uniform_int<> distribution_type;
  typedef boost::variate_generator<base_generator_type&, distribution_type> 
    gen_type;
  gen_type die_gen(generator, distribution_type(0, sxs.size()));
  boost::generator_iterator<gen_type> die(&die_gen);
  return sxs[*die];
}

stsDiffusion3p1::stsDiffusion3p1(  const char* infilename, 
				   const char* outfilename,
				   fixnum_t sigmamin,
				   fixnum_t sigmamax,
				   fixnum_t numtrials )
  : _ifilename(infilename), _ofilename(outfilename), 
    _sigmin(sigmamin), _sigmax(sigmamax), _numtrials(numtrials)
{
  _loadSpacetimeDataFromFile();
  _computeMaxVolSandwich();
  for( fixnum_t sigma = _sigmin; sigma <= _sigmax; ++sigma )
    _sigmaVsRetProb[sigma] = 0.0;
}

stsDiffusion3p1::~stsDiffusion3p1( void )
{
}

void stsDiffusion3p1::_loadSpacetimeDataFromFile( void )
{
  std::ifstream fin( _ifilename.c_str() );
  if( fin.is_open() )
    {
      std::string line;
      getline(fin,line);
      _parseParameterLine(line);
      while (getline(fin,line))
	{
	  bignum_t id;
	  simplex4 sx;
	  _parseDataLine(line, id, sx);
	  _4sxstore[id] = sx;
	}
    }
  fin.close();
}

void stsDiffusion3p1::_parseParameterLine( std::string& line )
{
  std::string::size_type begIdx = 0, endIdx = line.find_first_of(" ",begIdx);
  _bctype = line.substr(begIdx,endIdx);

  begIdx = endIdx + 1; endIdx = line.find_first_of(" ",begIdx);
  _stopology = line.substr(begIdx,endIdx-begIdx);

  begIdx = endIdx + 1; endIdx = line.find_first_of(" ",begIdx);
  _numT = 
    static_cast<fixnum_t>(atoi(line.substr(begIdx,endIdx-begIdx).c_str()));

  begIdx = endIdx + 1; endIdx = line.find_first_of(" ",begIdx);
  _initVol = 
    static_cast<fixnum_t>(atoi(line.substr(begIdx,endIdx-begIdx).c_str()));;

  begIdx = endIdx + 1; endIdx = line.find_first_of(" ",begIdx);
  _lastUsedPoint = strtoull(line.substr(begIdx,endIdx-begIdx).c_str(),0,10);

  begIdx = endIdx + 1; endIdx = line.find_first_of(" ",begIdx);
  _lastUsedId = strtoull(line.substr(begIdx,endIdx-begIdx).c_str(),0,10);

  begIdx = endIdx + 1; endIdx = line.find_first_of(" ",begIdx);
  _N0 = strtoull(line.substr(begIdx,endIdx-begIdx).c_str(),0,10);

  begIdx = endIdx + 1; endIdx = line.find_first_of(" ",begIdx);
  _N1_SL = strtoull(line.substr(begIdx,endIdx-begIdx).c_str(),0,10);

  begIdx = endIdx + 1; endIdx = line.find_first_of(" ",begIdx);
  _N1_TL = strtoull(line.substr(begIdx,endIdx-begIdx).c_str(),0,10);

  begIdx = endIdx + 1; endIdx = line.find_first_of(" ",begIdx);
  _N2_SL = strtoull(line.substr(begIdx,endIdx-begIdx).c_str(),0,10);

  begIdx = endIdx + 1; endIdx = line.find_first_of(" ",begIdx);
  _N2_TL = strtoull(line.substr(begIdx,endIdx-begIdx).c_str(),0,10);

  begIdx = endIdx + 1; endIdx = line.find_first_of(" ",begIdx);
  _N3_TL_31 = strtoull(line.substr(begIdx,endIdx-begIdx).c_str(),0,10);

  begIdx = endIdx + 1; endIdx = line.find_first_of(" ",begIdx);
  _N3_TL_22 = strtoull(line.substr(begIdx,endIdx-begIdx).c_str(),0,10);

  begIdx = endIdx + 1; endIdx = line.find_first_of(" ",begIdx);
  _N3_SL = strtoull(line.substr(begIdx,endIdx-begIdx).c_str(),0,10);

  begIdx = endIdx + 1; endIdx = line.find_first_of(" ",begIdx);
  _N4_TL_41 = strtoull(line.substr(begIdx,endIdx-begIdx).c_str(),0,10);

  begIdx = endIdx + 1; endIdx = line.find_first_of(" ",begIdx);
  _N4_TL_32 = strtoull(line.substr(begIdx,endIdx-begIdx).c_str(),0,10);

  begIdx = endIdx + 1; endIdx = line.find_first_of(" ",begIdx);
  _eps = strtod(line.substr(begIdx,endIdx-begIdx).c_str(),0);

  begIdx = endIdx + 1; endIdx = line.find_first_of(" ",begIdx);
  _kappa0 = strtod(line.substr(begIdx,endIdx-begIdx).c_str(),0);

  begIdx = endIdx + 1; endIdx = line.find_first_of(" ",begIdx);
  _delta = strtod(line.substr(begIdx,endIdx-begIdx).c_str(),0);

  begIdx = endIdx + 1; endIdx = line.find_first_of(" ",begIdx);
  _kappa4 = strtod(line.substr(begIdx,endIdx-begIdx).c_str(),0);

}

void stsDiffusion3p1::_parseDataLine( std::string& line, 
					    bignum_t& id, 
					    simplex4& sx )
{
  std::string::size_type begIdx = 0, endIdx = line.find_first_of(" ",begIdx);
  sx._type = static_cast<fixnum_t>(atoi(line.substr(begIdx,endIdx).c_str()));
  
  begIdx = endIdx + 1; endIdx = line.find_first_of(" ",begIdx);
  sx._tmlohi.first = 
    static_cast<unsigned short>(atoi(line.substr(begIdx,endIdx-begIdx).c_str()));
  
  begIdx = endIdx + 1; endIdx = line.find_first_of(" ",begIdx);
  sx._tmlohi.second = 
    static_cast<unsigned short>(atoi(line.substr(begIdx,endIdx-begIdx).c_str()));
  
  for(size_t ndx=0; ndx<5; ++ndx)
    {
      begIdx = endIdx + 1; endIdx = line.find_first_of(" ",begIdx);
      sx._points[ndx] = (strtoull(line.substr(begIdx,endIdx-begIdx).c_str(),
				  0,10));
    }

  for(size_t ndx=0; ndx<5; ++ndx)
    {
      begIdx = endIdx + 1; endIdx = line.find_first_of(" ",begIdx);
      sx._nbrids[ndx] = (strtoull(line.substr(begIdx,endIdx-begIdx).c_str(),
				  0,10));
    }
  
  begIdx = endIdx + 1; endIdx = line.find_first_of(" ",begIdx);
  id = (strtoull(line.substr(begIdx,endIdx-begIdx).c_str(),0,10));
}

bignum_t stsDiffusion3p1::_countSimplicesInSandwich( fixnum_t tmlo, 
							   fixnum_t tmhi )
{
  bignum_t count = 0;
  fixnum_pair_t tmlohi = std::make_pair(tmlo,tmhi%_numT);
  for( S4SConstIter citer = _4sxstore.begin(); 
       citer != _4sxstore.end(); 
       ++citer )
    if(citer->second._tmlohi == tmlohi)
      count++;
  return count;
}

void stsDiffusion3p1::_computeMaxVolSandwich( void )
{
  // set up a map with (0,1), (1,2),...,(T-1,0) as the keys with all values 
  // initially equal to 0;
  typedef std::map<fixnum_pair_t,fixnum_t> foomap_t;
  typedef std::map<fixnum_pair_t,fixnum_t>::const_iterator foociter_t;
  typedef std::map<fixnum_pair_t,fixnum_t>::iterator fooiter_t;

  foomap_t lohiVmap;
  for( fixnum_t ts = 0; ts < _numT; ++ts )
    lohiVmap[std::make_pair(ts, (ts+1)%_numT)] = 0;

  // iterate through the three simplex store, incrementing the appropriate 
  // sandwich volume count
  for( S4SConstIter citer = _4sxstore.begin(); 
       citer != _4sxstore.end(); 
       ++citer )
    lohiVmap[citer->second._tmlohi]++;
  // the sandwich with the maximum volume is _maxVolSandwich;
  _maxvolume = 0;
  for( foociter_t citer = lohiVmap.begin(); citer != lohiVmap.end(); ++citer )
    if( citer->second > _maxvolume )
      {
	_maxVolSandwich = citer->first;
	_maxvolume = citer->second;
      }
}

std::vector<bignum_t> 
stsDiffusion3p1::_getSimplicesInSandwich
( fixnum_pair_t tslohi )
{
  std::vector<bignum_t> sxsInSandwich;
  for( S4SConstIter citer = _4sxstore.begin(); 
       citer != _4sxstore.end(); 
       ++citer )
    if( tslohi == citer->second._tmlohi )
      sxsInSandwich.push_back( citer->first );
  return sxsInSandwich;
}

void stsDiffusion3p1::diffusion(void)
{
  std::vector<bignum_t> sxids = _getSimplicesInSandwich(_maxVolSandwich);
  
  for( fixnum_t trial = 0; trial <= _numtrials; ++trial )
    {
      bignum_t i0 = pick_random_simplex(sxids);
      simplex4& sx = _4sxstore[i0];

      _KTstore[kt_key_t(i0,i0,0)] = 1.0;
      _KTstore[kt_key_t(i0,i0,1)] = 0.0;
      for( size_t ndx = 0; ndx < 5; ++ndx )
	_KTstore[kt_key_t(sx._nbrids[ndx], i0, 1)] = 0.2;
      for( fixnum_t sigma = _sigmin; sigma <= _sigmax; ++sigma )
	_sigmaVsRetProb[sigma] += KT(i0,i0,sigma);
    }
  for( fixnum_t sigma = _sigmin; sigma <= _sigmax; ++sigma )
    _sigmaVsRetProb[sigma] /= _numtrials;

  std::ofstream fout( _ofilename.c_str() );
  if(fout.is_open())
    {
      for(fixnum_t sigma = _sigmin; sigma <= _sigmax; ++sigma)
	fout << sigma << "\t" << std::setprecision(12) << std::fixed << 
	  _sigmaVsRetProb[sigma] << std::endl;
    }
  fout.close();
}

real_t stsDiffusion3p1::KT( bignum_t j, bignum_t i0, fixnum_t sig )
{
  if ( sig == 0 )
    {
      if( i0 == j )
	return 1.0;
      else
	return 0.0;
    }

  kt_key_t k(j,i0,sig);
  if( _KTstore.find( k ) == _KTstore.end() )
    {
      real_t kt_tot = 0.0;
      simplex4& sx = _4sxstore[j];
      for( size_t ndx = 0; ndx < 5; ++ndx )
	{
	  kt_tot += KT( sx._nbrids[ndx], i0, sig-1 );
	}
      _KTstore[k] = kt_tot/5.0;
    }
  return _KTstore[k];
}
