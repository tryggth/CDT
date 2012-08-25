/*
 * sshDiffusion2p1SI.cpp
 */
#include "sshDiffusion2p1SI.hpp"
#include <ctime>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
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

sshDiffusion2p1::sshDiffusion2p1(const char* infilename, 
				 const char* outfilename,
				 fixnum_t numtrials)
  : _ifilename( infilename ), _ofilename( outfilename ), _numtrials(numtrials) 
{
  _loadDataFromFile();
  for( fixnum_t d = 1; d <= 1000; d++ )
    _rVsN[d] = 0.0;
}

sshDiffusion2p1::~sshDiffusion2p1( void )
{
}

void sshDiffusion2p1::_loadDataFromFile( void )
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

void sshDiffusion2p1::_parseParameterLine( std::string& line )
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

void sshDiffusion2p1::_parseDataLine( std::string& line, bignum_t& id, s2simplex& sx )
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

void sshDiffusion2p1::_getS2SimplicesInSlice( fixnum_t tslice, 
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

void sshDiffusion2p1::_computeMaxVolSlice( void )
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

void sshDiffusion2p1::hausdorff(void)
{
  _computeMaxVolSlice();
  std::vector<bignum_t> sxInMaxVol;
  _getS2SimplicesInSlice( _maxVolSlice, sxInMaxVol );
  fixnum_t d = 0, d_max = 0;
  for( fixnum_t trial = 0; trial < _numtrials; ++trial )
    {
      bignum_t i0 = pick_random_simplex( sxInMaxVol );
      s2simplex& sx = _s2sxstore[i0];
      _visited.clear();
      _nbors.clear();
      for( size_t ndx = 0; ndx < 3; ++ndx )
	{
	  _visited.insert( sx._nbrids[ndx] );
	  _nbors.insert( std::make_pair( kt_key_t(i0,1), sx._nbrids[ndx] ) );
	}
      d = 2;
      while ( d <= 100 ) //_visited.size() < sxInMaxVol.size() )
	{
	  std::cout << "d = " << d << " and " << _visited.size() << " of " 
		    << sxInMaxVol.size() << " and " << _nbors.count(kt_key_t(i0,d)) << std::endl;
	  std::set<bignum_t>::const_iterator citer;
	  std::list<bignum_t> newones;
	  for( citer = _visited.begin(); citer != _visited.end(); ++citer )
	    {
	      s2simplex& sx = _s2sxstore[*citer];
	      for( size_t ndx = 0; ndx < 3; ++ndx )
		{
		  // don't want to invalidate _visited iterators by directly
		  // inserting into _visited. build a noduplicate list called
		  // newones
		  if( (_visited.find( sx._nbrids[ndx] ) == _visited.end() ) &&
		      (std::find(newones.begin(), 
				 newones.end(), 
				 sx._nbrids[ndx]) == newones.end()) )
		    newones.push_back( sx._nbrids[ndx] );
		}
	    }
	  for( std::list<bignum_t>::const_iterator it = newones.begin();
	       it != newones.end(); ++it )
	    {
	      _nbors.insert( std::make_pair( kt_key_t(i0,d), *it ) );
	      _visited.insert( *it );
	    }
	  d++;
	}
      d_max = d;
      for( d = 1; d <= d_max; ++d )
	{
	  _rVsN[d] += _nbors.count( kt_key_t(i0, d) );
	}
    }
  //  for( fixnum_t d = 1; d <= 1000; ++d )
  //  _rVsN[d] /= _numtrials;

  std::ofstream fout( _ofilename.c_str() );
  if(fout.is_open())
    {
      for(fixnum_t d = 1; d <= 1000; ++d)
	fout << d << "\t" << std::setprecision(6) << std::fixed <<
	  _rVsN[d] << std::endl;
    }
  fout.close();
}

