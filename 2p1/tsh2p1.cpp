/*
 * tsh2p1.cpp
 */
#include "tsh2p1.hpp"
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

fixnum_t pick_random_neighbor( void )
{
  static base_generator_type generator(static_cast<unsigned int>(std::time(0)));
  typedef boost::uniform_int<> distribution_type;
  typedef boost::variate_generator<base_generator_type&, distribution_type> 
    gen_type;
  gen_type die_gen(generator, distribution_type(0, 3));
  boost::generator_iterator<gen_type> die(&die_gen);
  return *die;
}

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

ThickSliceHausdorffDim::ThickSliceHausdorffDim(  const char* infilename, 
						 const char* outfilename,
						 fixnum_t threshvol )
  : _ifilename(infilename), _ofilename(outfilename), 
    _thresholdVolume( threshvol )
{
  _loadSpacetimeDataFromFile();
}

ThickSliceHausdorffDim::~ThickSliceHausdorffDim( void )
{
}

void ThickSliceHausdorffDim::_loadSpacetimeDataFromFile( void )
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
	  simplex3 sx;
	  _parseDataLine(line, id, sx);
	  _3sxstore[id] = sx;
	}
    }
  fin.close();
}

void ThickSliceHausdorffDim::_parseParameterLine( std::string& line )
{
  std::string::size_type begIdx = 0, endIdx = line.find_first_of(" ",begIdx);
  _bctype = line.substr(begIdx,endIdx);

  begIdx = endIdx + 1; endIdx = line.find_first_of(" ",begIdx);
  _stopology = line.substr(begIdx,endIdx-begIdx);

  begIdx = endIdx + 1; endIdx = line.find_first_of(" ",begIdx);
  _numT = static_cast<fixnum_t>(atoi(line.substr(begIdx,endIdx-begIdx).c_str()));

  begIdx = endIdx + 1; endIdx = line.find_first_of(" ",begIdx);
  _initVol = static_cast<fixnum_t>(atoi(line.substr(begIdx,endIdx-begIdx).c_str()));;

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
  _eps = strtod(line.substr(begIdx,endIdx-begIdx).c_str(),0);

  begIdx = endIdx + 1; endIdx = line.find_first_of(" ",begIdx);
  _k0 = strtod(line.substr(begIdx,endIdx-begIdx).c_str(),0);

  begIdx = endIdx + 1; endIdx = line.find_first_of(" ",begIdx);
  _k3 = strtod(line.substr(begIdx,endIdx-begIdx).c_str(),0);

  begIdx = endIdx + 1; endIdx = line.find_first_of(" ",begIdx);
  _alpha = strtod(line.substr(begIdx,endIdx-begIdx).c_str(),0);
}

void ThickSliceHausdorffDim::_parseDataLine( std::string& line, bignum_t& id, simplex3& sx )
{
  std::string::size_type begIdx = 0, endIdx = line.find_first_of(" ",begIdx);
  sx._type = static_cast<fixnum_t>(atoi(line.substr(begIdx,endIdx).c_str()));
  
  begIdx = endIdx + 1; endIdx = line.find_first_of(" ",begIdx);
  sx._tmlohi.first = static_cast<unsigned short>(atoi(line.substr(begIdx,endIdx-begIdx).c_str()));
  
  begIdx = endIdx + 1; endIdx = line.find_first_of(" ",begIdx);
  sx._tmlohi.second = static_cast<unsigned short>(atoi(line.substr(begIdx,endIdx-begIdx).c_str()));
  
  for(size_t ndx=0; ndx<4; ++ndx)
    {
      begIdx = endIdx + 1; endIdx = line.find_first_of(" ",begIdx);
      sx._points[ndx] = (strtoull(line.substr(begIdx,endIdx-begIdx).c_str(),0,10));
    }

  for(size_t ndx=0; ndx<4; ++ndx)
    {
      begIdx = endIdx + 1; endIdx = line.find_first_of(" ",begIdx);
      sx._nbrids[ndx] = (strtoull(line.substr(begIdx,endIdx-begIdx).c_str(),0,10));
    }
  
  begIdx = endIdx + 1; endIdx = line.find_first_of(" ",begIdx);
  id = (strtoull(line.substr(begIdx,endIdx-begIdx).c_str(),0,10));
}

void ThickSliceHausdorffDim::_computeStalkSandwiches( void )
{
  // set up a map with (0,1), (1,2),...,(T-1,0) as the keys with all values initially equal to 0;
  typedef std::map<fixnum_pair_t,fixnum_t> foomap_t;
  typedef std::map<fixnum_pair_t,fixnum_t>::const_iterator foociter_t;
  typedef std::map<fixnum_pair_t,fixnum_t>::iterator fooiter_t;

  foomap_t lohiVmap;
  for( fixnum_t ts = 0; ts < _numT; ++ts )
    lohiVmap[std::make_pair(ts, (ts+1)%_numT)] = 0;

  // iterate through the three simplex store, incrementing the appropriate 
  // sandwich volume count
  for( S3SConstIter citer = _3sxstore.begin(); citer != _3sxstore.end(); ++citer )
    lohiVmap[citer->second._tmlohi]++;

  // all the sandwiches with volume less than threshold get into the 
  // _stalkSandwiches list.
  for( foociter_t citer = lohiVmap.begin(); citer != lohiVmap.end(); ++citer )
    if( citer->second <= _thresholdVolume )
      _stalkSandwiches.push_back( citer->first );
}

void ThickSliceHausdorffDim::_getSimplicesInSandwich( fixnum_pair_t tslohi,
						      std::vector<bignum_t>& sxsInSandwich)
{
  sxsInSandwich.clear();
  for( S3SConstIter citer = _3sxstore.begin(); 
       citer != _3sxstore.end(); 
       ++citer )
    if( tslohi == citer->second._tmlohi )
      sxsInSandwich.push_back( citer->first );
}

void ThickSliceHausdorffDim::_getNeighbors(const std::set<bignum_t>& oldnbors, 
					    std::set<bignum_t>& newnbors)
{
  std::set<bignum_t>::const_iterator it;
  for( it = oldnbors.begin(); it != oldnbors.end(); ++it )
    {
      simplex3& sx = _3sxstore[*it];
      for( size_t ndx = 0; ndx < 4; ++ndx )
	{
	  simplex3& nbor = _3sxstore[sx._nbrids[ndx]];
	  if( nbor._tmlohi == sx._tmlohi )
	    newnbors.insert(sx._nbrids[ndx]);
	}
    }
}

void ThickSliceHausdorffDim::computeHausdorffDimension(void)
{
  _computeStalkSandwiches();

  for (fixnum_t ts = 0; ts < _numT; ++ts )
    {
      std::pair<fixnum_t,fixnum_t> tmlohi( ts, (ts+1)%_numT );
      if(std::find( _stalkSandwiches.begin(),
		    _stalkSandwiches.end(),
		    tmlohi ) 
	 == _stalkSandwiches.end())
	{
	  // this thick slice is not in the stalk. 
	  //get all the simplices in this thick slice
	  std::vector<bignum_t> sxsInThickSlice;
	  _getSimplicesInSandwich( tmlohi, sxsInThickSlice );
	  for ( fixnum_t ntrial = 0; ntrial < _numtrials; ++ntrial )
	    {
	      bignum_t randomSx = pick_random_simplex( sxsInThickSlice );
	      std::set<bignum_t> visited;
	      visited.insert(randomSx);
	      bignum_t n = 1;
	      std::set<bignum_t> newnbors, oldnbors;
	      oldnbors.insert(randomSx);
	      while ( visited.size() < sxsInThickSlice.size() )
		{
		  newnbors.clear();
		  _getNeighbors( oldnbors, newnbors );
		  _rVsnr[n] += newnbors.size();
		  n++;
		  // store newnbors in the set of visited simplices
		  for( std::set<bignum_t>::const_iterator it = newnbors.begin();
		       it != newnbors.end(); ++it )
		    visited.insert(*it);
		  oldnbors = newnbors;
		}
	    }
	}
    }
  std::ofstream fout( _ofilename.c_str() );
  if(fout.is_open())
    {
    }
  fout.close();
}
