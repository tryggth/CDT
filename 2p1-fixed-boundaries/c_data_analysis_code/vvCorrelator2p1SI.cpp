/*
 * vvCorrelator2p1SI.cpp
 */
#include "vvCorrelator2p1SI.hpp"
#include <algorithm>
#include <ctime>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <list>

vvCorrelator2p1::vvCorrelator2p1(  const char* infilename, 
				   const char* outfilename,
				   fixnum_t dmin, fixnum_t dmax )
  : _ifilename(infilename), _ofilename(outfilename), _dmin(dmin), _dmax(dmax)
{
  _loadSpatialSliceDataFromFile();
  _computeMaxVolSlice();
  _vvCorrelations.reserve( 2*_dmax + 1 );
}

vvCorrelator2p1::~vvCorrelator2p1( void )
{
}

void vvCorrelator2p1::_loadSpatialSliceDataFromFile( void )
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

void vvCorrelator2p1::_parseParameterLine( std::string& line )
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

void 
vvCorrelator2p1::_parseDataLine
( std::string& line, bignum_t& id, s2simplex& sx )
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

bignum_t 
vvCorrelator2p1::_countS2SimplicesInSlice
( fixnum_t tslice ) 
{
  bignum_t count = 0;
  fixnum_t ts = tslice%_numT;
  for( S2SSCIter citer = _s2sxstore.begin(); 
       citer != _s2sxstore.end(); 
       ++citer )
    if(citer->second._tslice == ts)
      count++;
  return count;
}

void 
vvCorrelator2p1::_computeMaxVolSlice
( void )
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

void vvCorrelator2p1::computeCorrelations( void )
{
  // compute the correlations
  for( fixnum_t delta = _dmin; delta <= _dmax; ++delta )
    {
      double tot = 0.0;
      for( fixnum_t ts = 0; ts < _numT; ++ts )
	tot += (_countS2SimplicesInSlice(ts)*
		_countS2SimplicesInSlice(offset(ts,delta)));
      tot /= (_numT*_numT);
      _vvCorrelations[delta + _dmax] = tot;
    }

  // write the correlations to file
  std::ofstream fout( _ofilename.c_str() );
  if(fout.is_open())
    {
      for(fixnum_t delta = _dmin; delta <= _dmax; ++delta)
	fout << delta << "\t" << 
	  std::setprecision(12) << std::fixed << 
	  _vvCorrelations[delta + _dmax] << std::endl;
    }
  fout.close();
}
