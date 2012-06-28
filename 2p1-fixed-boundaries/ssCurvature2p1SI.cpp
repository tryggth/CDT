/*
 * ssCurvature2p1SI.cpp
 */
#include "ssCurvature2p1SI.hpp"
#include <ctime>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <set>
// boost libraries
#include <boost/generator_iterator.hpp>
#include <boost/random/linear_congruential.hpp>
#include <boost/random/uniform_int.hpp>
#include <boost/random/uniform_real.hpp>
#include <boost/random/variate_generator.hpp>
#include <boost/random/mersenne_twister.hpp>

typedef boost::mt19937 base_generator_type;

ssCurvature2p1::ssCurvature2p1(const char* infilename, 
			       const char* outfilename,
			       fixnum_t offset)
  : _ifilename( infilename ), _ofilename( outfilename ),
    _offset(offset)
{
  _loadDataFromFile();
}

void ssCurvature2p1::_loadDataFromFile( void )
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
	  s2simplex sx;
	  _parseDataLine(line, id, sx);
	  _s2sxstore[id] = sx;
	}
    }
  fin.close();
}

void ssCurvature2p1::_parseParameterLine( std::string& line )
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
  _N2_SL = strtoull(line.substr(begIdx,endIdx-begIdx).c_str(),0,10);

  // read it again
  begIdx = endIdx + 1; endIdx = line.find_first_of(" ",begIdx);
  _N2_SL = strtoull(line.substr(begIdx,endIdx-begIdx).c_str(),0,10);

  begIdx = endIdx + 1; endIdx = line.find_first_of(" ",begIdx);
  _k0 = strtod(line.substr(begIdx,endIdx-begIdx).c_str(),0);

  begIdx = endIdx + 1; endIdx = line.find_first_of(" ",begIdx);
  _k3 = strtod(line.substr(begIdx,endIdx-begIdx).c_str(),0);
}

void ssCurvature2p1::_parseDataLine( std::string& line, bignum_t& id, s2simplex& sx )
{
  std::string::size_type begIdx = 0, endIdx = line.find_first_of(" ",begIdx);
  id = (strtoull(line.substr(begIdx,endIdx-begIdx).c_str(),0,10));
  begIdx = endIdx + 1; endIdx = line.find_first_of(" ",begIdx);
  sx._tslice = 
    static_cast<fixnum_t>(atoi(line.substr(begIdx,endIdx-begIdx).c_str()));
  for(size_t ndx=0; ndx<3; ++ndx)
    {
      begIdx = endIdx + 1; endIdx = line.find_first_of(" ",begIdx);
      sx._points[ndx] = 
	(strtoull(line.substr(begIdx,endIdx-begIdx).c_str(),0,10));
    }
  for(size_t ndx=0; ndx<3; ++ndx)
    {
      begIdx = endIdx + 1; endIdx = line.find_first_of(" ",begIdx);
      sx._nbrids[ndx] = 
	(strtoull(line.substr(begIdx,endIdx-begIdx).c_str(),0,10));
    }
}

void ssCurvature2p1::_getS2SimplicesInSlice( fixnum_t tslice, 
					      std::vector<bignum_t>& retval )
{
  bignum_t count = 0;
  fixnum_t ts = tslice%_numT;
  retval.clear();
  for( S2SSCIter citer = _s2sxstore.begin(); 
       citer != _s2sxstore.end(); 
       ++citer )
    if(citer->second._tslice == ts)
      retval.push_back( citer->first );
}

void ssCurvature2p1::_computeMaxVolSlice( void )
{
  fixnum_t *volCounts = new fixnum_t[_numT];

  for( fixnum_t ts = 0; ts < _numT; ++ts )
    volCounts[ts] = 0;

  for (S2SSCIter citer = _s2sxstore.begin(); 
       citer != _s2sxstore.end(); 
       ++citer)
    {
      volCounts[ citer->second._tslice ]++;
    }

  bignum_t currmax = 0;
  for( fixnum_t ts = 0; ts < _numT; ++ts )
    {
      if( volCounts[ts] > currmax )
	{
	  _maxVolSlice = ts;
	  currmax = volCounts[ts];
	}
    }
  delete [] volCounts;
}

double ssCurvature2p1::curvature(void)
{
  _computeMaxVolSlice();
  std::vector<bignum_t> sxInMaxVol;
  // compute the offset from _maxVolSlice, and compute the spectral dimension
  // for this slice 
  fixnum_t ts = offset( _maxVolSlice, _offset );
  _getS2SimplicesInSlice( ts, sxInMaxVol );

  // build a set of all the points in the slice
  std::set<bignum_t> ptsInSlice;
  std::vector<bignum_t>::const_iterator citer;
  for( citer = sxInMaxVol.begin(); citer != sxInMaxVol.end(); ++citer )
    {
      s2simplex& s2sx = _s2sxstore[*citer];
      for(size_t ndx = 0; ndx < 3; ++ndx )
	ptsInSlice.insert( s2sx._points[ndx] );
    }
  
  // for each point in this set, calculate the number of triangles around this
  // point, call this N. curvature at this point is [(6-N)xpi]/3. Summing over
  // all the points in the slice, we get total curvature of the slice
  double totCurvature = 0.0;
  long N_tot = 0;
  for( std::set<bignum_t>::const_iterator siter = ptsInSlice.begin();
       siter != ptsInSlice.end(); ++siter )
    {
      long N = _trianglesAroundPoint( *siter, sxInMaxVol );
      double curvature = (6-N)*M_PI/3.0;      
      N_tot += N;
      totCurvature += curvature;
    }
  return totCurvature;
}

long ssCurvature2p1::_trianglesAroundPoint( bignum_t p, 
					    const std::vector<bignum_t>& t )
{
  long ret = 0;
  std::vector<bignum_t>::const_iterator citer;
  for( citer = t.begin(); citer != t.end(); ++citer )
    {
      s2simplex& triangle = _s2sxstore[*citer];
      for( size_t ndx = 0; ndx < 3; ++ndx )
	if( triangle._points[ndx] == p )
	  {
	    ret++;
	  }
    }
  return ret;
}

