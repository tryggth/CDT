/*
 * sssDiffusion3p1SI.cpp
 */
#include "sssDiffusion3p1SI.hpp"
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
  
  static base_generator_type 
    generator(static_cast<unsigned int>(std::time(0)));
  typedef boost::uniform_int<> distribution_type;
  typedef 
    boost::variate_generator<base_generator_type&, distribution_type> 
    gen_type;
  gen_type die_gen(generator, distribution_type(0, sxs.size()));
  boost::generator_iterator<gen_type> die(&die_gen);
  return sxs[*die];
}

sssDiffusion3p1::sssDiffusion3p1(const char* infilename, 
				 const char* outfilename,
				 fixnum_t sigmin,
				 fixnum_t sigmax,
				 fixnum_t offset,
				 fixnum_t numtrials)
  : _ifilename( infilename ), _ofilename( outfilename ),
    _sigmin(sigmin), _sigmax(sigmax), _offset(offset), _numtrials(numtrials) 
{
  _loadDataFromFile();
  for( fixnum_t sig = sigmin; sig <= sigmax; ++sig )
    _sigVsRetProb[sig] = 0.0;
}

sssDiffusion3p1::~sssDiffusion3p1( void )
{
}

void sssDiffusion3p1::_loadDataFromFile( void )
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
	  s3simplex sx;
	  _parseDataLine(line, id, sx);
	  _s3sxstore[id] = sx;
	}
    }
  fin.close();
}

void sssDiffusion3p1::_parseParameterLine( std::string& line )
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
  _N3_SL = strtoull(line.substr(begIdx,endIdx-begIdx).c_str(),0,10);

  // read it again
  begIdx = endIdx + 1; endIdx = line.find_first_of(" ",begIdx);
  _N3_SL = strtoull(line.substr(begIdx,endIdx-begIdx).c_str(),0,10);

  begIdx = endIdx + 1; endIdx = line.find_first_of(" ",begIdx);
  _kappa0 = strtod(line.substr(begIdx,endIdx-begIdx).c_str(),0);

  begIdx = endIdx + 1; endIdx = line.find_first_of(" ",begIdx);
  _delta = strtod(line.substr(begIdx,endIdx-begIdx).c_str(),0);

  begIdx = endIdx + 1; endIdx = line.find_first_of(" ",begIdx);
  _kappa4 = strtod(line.substr(begIdx,endIdx-begIdx).c_str(),0);

}

void sssDiffusion3p1::_parseDataLine( std::string& line, 
				      bignum_t& id, 
				      s3simplex& sx )
{
  std::string::size_type begIdx = 0, endIdx = line.find_first_of(" ",begIdx);
  id = (strtoull(line.substr(begIdx,endIdx-begIdx).c_str(),0,10));
  begIdx = endIdx + 1; endIdx = line.find_first_of(" ",begIdx);
  sx._tslice = 
    static_cast<fixnum_t>(atoi(line.substr(begIdx,endIdx-begIdx).c_str()));
  for(size_t ndx=0; ndx<4; ++ndx)
    {
      begIdx = endIdx + 1; endIdx = line.find_first_of(" ",begIdx);
      sx._points[ndx] = 
	(strtoull(line.substr(begIdx,endIdx-begIdx).c_str(),0,10));
    }
  for(size_t ndx=0; ndx<4; ++ndx)
    {
      begIdx = endIdx + 1; endIdx = line.find_first_of(" ",begIdx);
      sx._nbrids[ndx] = 
	(strtoull(line.substr(begIdx,endIdx-begIdx).c_str(),0,10));
    }
}

void sssDiffusion3p1::_getS3SimplicesInSlice( fixnum_t tslice, 
					      std::vector<bignum_t>& retval )
{
  bignum_t count = 0;
  fixnum_t ts = tslice%_numT;
  retval.clear();
  for( S3SSCIter citer = _s3sxstore.begin(); 
       citer != _s3sxstore.end(); 
       ++citer )
    if(citer->second._tslice == ts)
      retval.push_back( citer->first );
}

void sssDiffusion3p1::_computeMaxVolSlice( void )
{
  fixnum_t *volCounts = new fixnum_t[_numT];

  for( fixnum_t ts = 0; ts < _numT; ++ts )
    volCounts[ts] = 0;

  for (S3SSCIter citer = _s3sxstore.begin(); 
       citer != _s3sxstore.end(); 
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

void sssDiffusion3p1::diffusion(void)
{
  _computeMaxVolSlice();
  std::vector<bignum_t> sxInMaxVol;
  // compute the offset from _maxVolSlice, and compute the spectral dimension
  // for this slice 
  fixnum_t ts = offset( _maxVolSlice, _offset );
  _getS3SimplicesInSlice( ts, sxInMaxVol );
  // do the exact diffusion for _numtrial simplices in the max volume slice
  for( size_t ndx = 0; ndx < _numtrials; ++ndx )
    {
      bignum_t i0 = sxInMaxVol[ndx];
      s3simplex& sx = _s3sxstore[i0];
      _KTstore[kt_key_t(i0,i0,0)] = 1.0;
      _KTstore[kt_key_t(i0,i0,1)] = 0.0;
      for(size_t ndx = 0; ndx < 4; ++ndx )
	_KTstore[kt_key_t(sx._nbrids[ndx], i0, 1)] = 1.0/4.0;

      for( fixnum_t sigma = _sigmin; sigma <= _sigmax; ++sigma )
	{
	  _sigVsRetProb[sigma] += KT(i0,i0,sigma);
	}
    }
  // renormailze _sigVsRetProb
  for( fixnum_t sig = _sigmin; sig <= _sigmax; ++sig )
    _sigVsRetProb[sig] /= _numtrials;

  // write the contents of sigma versus return probabilities to file
  std::ofstream fout( _ofilename.c_str() );
  if(fout.is_open())
    {
      for(fixnum_t sigma = _sigmin; sigma <= _sigmax; ++sigma)
	fout << sigma << "\t" << std::setprecision(6) << std::fixed <<
	  _sigVsRetProb[sigma] << std::endl;
    }
  fout.close();
}

real_t sssDiffusion3p1::KT( bignum_t j, bignum_t i0, fixnum_t sig )
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
      s3simplex& sx = _s3sxstore[j];
      for( size_t ndx = 0; ndx < 4; ++ndx )
	{
	  kt_tot += KT( sx._nbrids[ndx], i0, sig-1 );
	}
      _KTstore[k] = kt_tot/4.0;
    }
  return _KTstore[k];
}
