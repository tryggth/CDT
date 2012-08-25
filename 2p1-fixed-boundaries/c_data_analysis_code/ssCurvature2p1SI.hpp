#ifndef __ssCurvature__hpp__
#define __ssCurvature__hpp__

/*
 * ssCurvature2p1 --- computes the spatial slice curvature for 2p1 spacetimes. 
 * The spectral dimension is computed for the maximum volume slice +- the 
 * specified offset.
 */
#include <list>
#include <vector>
#include <map>
#include <string>
#include <boost/tuple/tuple.hpp>
#include <boost/unordered_map.hpp>
#include <boost/functional/hash.hpp>

typedef unsigned long long bignum_t;
typedef unsigned short fixnum_t;
typedef std::pair<fixnum_t,fixnum_t> fixnum_pair_t;
typedef double real_t;

// s2simplex stands for spatial 2-simplex; a generic 2-simplex could be 
// spacelike or timelike id ts p0 p1 p2 n0 n1 n2
struct s2simplex
{
  fixnum_t _tslice;
  bignum_t _points[3];
  bignum_t _nbrids[3];

  s2simplex( void ) : _tslice(-1)
  {
    for( size_t ndx = 0; ndx < 3; ++ndx )
      _points[ndx] = _nbrids[ndx] = 0;
  }

  s2simplex( fixnum_t tslice, bignum_t p0, bignum_t p1, bignum_t p2, 
	     bignum_t n0, bignum_t n1, bignum_t n2 )
    : _tslice( tslice )
  {
    _points[0] = p0; _points[1] = p1; _points[2] = p2;
    _nbrids[0] = n0; _nbrids[1] = n1; _nbrids[2] = n2;
  }

  s2simplex( const s2simplex& rhs ) : _tslice( rhs._tslice)
  {
    for( size_t ndx = 0; ndx < 3; ++ndx )
      {
	_points[ndx] = rhs._points[ndx];
	_nbrids[ndx] = rhs._nbrids[ndx];
      }
  }

  s2simplex& operator=( const s2simplex& rhs )
  {
    if( this != &rhs )
      {
	_tslice = rhs._tslice;
	for( size_t ndx = 0; ndx < 3; ++ndx )
	  {
	    _points[ndx] = rhs._points[ndx];
	    _nbrids[ndx] = rhs._nbrids[ndx];
	  }
      }
    return *this;
  }
};

typedef std::map<bignum_t, s2simplex> S2SimplexStore;
typedef std::map<bignum_t, s2simplex>::const_iterator S2SSCIter;
typedef std::map<bignum_t, s2simplex>::iterator S2SSIter;

class ssCurvature2p1
{
public:
  ssCurvature2p1( const char* infilename, 
		  const char* outfilename,
		  fixnum_t offset );
  
  double curvature(void);
    
private:
  // input and output data files
  std::string _ifilename;
  std::string _ofilename;

  // following parameters pertain to the spectral dimension calculation
  fixnum_t _offset;

  // following loaded from .s2sx2p1 file
  fixnum_t _numT;
  std::string _bctype;
  std::string _stopology;
  fixnum_t _initVol;
  bignum_t _N2_SL;
  double _k0;
  double _k3;

  inline fixnum_t offset( fixnum_t i0, fixnum_t k )
  {
    return( (i0 + (k % _numT) + _numT) % _numT );
  }

  // stores the spatial 2-simplex data
  S2SimplexStore _s2sxstore;

  // slice with the largest volume
  fixnum_t _maxVolSlice;

  // implementation methods
  void _loadDataFromFile( void );
  void _parseDataLine( std::string& line, bignum_t& id, s2simplex& sx );
  void _parseParameterLine( std::string& line );
  void _getS2SimplicesInSlice( fixnum_t tslice, std::vector<bignum_t>& sxs );
  void _computeMaxVolSlice( void );
  long _trianglesAroundPoint( bignum_t p, const std::vector<bignum_t>& t );

  // prohibit default initialization and copying
  ssCurvature2p1( void );
  ssCurvature2p1( const ssCurvature2p1& );
  ssCurvature2p1& operator=( const ssCurvature2p1& );
};

#endif
