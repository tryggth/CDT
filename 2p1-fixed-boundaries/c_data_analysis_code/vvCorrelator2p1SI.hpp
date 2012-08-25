#ifndef __vvCorrelator2p1SI__hpp__
#define __vvCorrelator2p1SI__hpp__

/*
 * Thu May  5 09:29:01 2011 
 * this is an exact simulation of the diffusion process; no random walks!
 * SI means Stalk Included
 */
#include <vector>
#include <map>
#include <string>

typedef unsigned long long bignum_t;
//typedef unsigned long fixnum_t;
typedef signed long fixnum_t;
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

class vvCorrelator2p1
{
public:
  vvCorrelator2p1( const char* infilename, 
		   const char* outfilename,
		   fixnum_t dmin, fixnum_t dmax ); 
  ~vvCorrelator2p1( void );

  inline fixnum_t maxVolume( void ) const
  {
    return _maxvolume;
  }

  void computeCorrelations(void);
private:
  // input and output data files
  std::string _ifilename;
  std::string _ofilename;

  // maximum volume sandwich
  fixnum_t _maxVolSlice;
  fixnum_t _maxvolume;
  fixnum_t _dmin;
  fixnum_t _dmax;

  // following loaded from .s2sx2p1 file
  fixnum_t _numT;
  std::string _bctype;
  std::string _stopology;
  fixnum_t _initVol;
  bignum_t _N2_SL;
  double _k0;
  double _k3;

  // the spatial 2-simplex store
  S2SimplexStore _s2sxstore;
  // the correlation values
  std::vector<real_t> _vvCorrelations;

  // implementation methods
  bignum_t _countS2SimplicesInSlice( fixnum_t tslice );
  void _loadSpatialSliceDataFromFile( void );
  void _parseDataLine( std::string& line, bignum_t& id, s2simplex& sx );
  void _parseParameterLine( std::string& line );

  void _computeMaxVolSlice( void );

  inline fixnum_t offset( fixnum_t i0, fixnum_t k )
  {
    return( (i0 + (k % _numT) + _numT) % _numT );
  }

  // prohibit copying
  vvCorrelator2p1( void );
  vvCorrelator2p1( const vvCorrelator2p1& );
  vvCorrelator2p1& operator=( const vvCorrelator2p1& );
};

#endif // __vvCorrelator2p1SI__hpp__
