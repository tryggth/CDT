#ifndef __ssHausdorf_3p1SI__hpp__
#define __ssHausdorf_3p1SI__hpp__

/*
 * ssh3p1SI.hpp --- computes the spatial slice hausdorf dimension for 3p1 
 * spacetimes. The hausdorf dimension is computed for the maximum volume slice 
 * +- the offset
 */
#include <list>
#include <vector>
#include <map>
#include <set>
#include <string>
#include <boost/tuple/tuple.hpp>
#include <boost/unordered_map.hpp>
#include <boost/functional/hash.hpp>

typedef unsigned long long bignum_t;
typedef signed short fixnum_t;
typedef std::pair<fixnum_t,fixnum_t> fixnum_pair_t;
typedef double real_t;

struct s3simplex
{
  fixnum_t _tslice;
  bignum_t _points[4];
  bignum_t _nbrids[4];

  s3simplex( void ) : _tslice(-1)
  {
    for( size_t ndx = 0; ndx < 4; ++ndx )
      _points[ndx] = _nbrids[ndx] = 0;
  }

  s3simplex( fixnum_t tslice, 
	     bignum_t p0, bignum_t p1, bignum_t p2, bignum_t p3, 
	     bignum_t n0, bignum_t n1, bignum_t n2, bignum_t n3 )
    : _tslice( tslice )
  {
    _points[0] = p0; _points[1] = p1; _points[2] = p2; _points[3] = p3;
    _nbrids[0] = n0; _nbrids[1] = n1; _nbrids[2] = n2; _nbrids[3] = n3;
  }

  s3simplex( const s3simplex& rhs ) : _tslice(rhs._tslice)
  {
    for( size_t ndx = 0; ndx < 4; ++ndx )
      {
	_points[ndx] = rhs._points[ndx];
	_nbrids[ndx] = rhs._nbrids[ndx];
      }
  }

  s3simplex& operator=( const s3simplex& rhs )
  {
    if( this != &rhs )
      {
	_tslice = rhs._tslice;
	for( size_t ndx = 0; ndx < 4; ++ndx )
	  {
	    _points[ndx] = rhs._points[ndx];
	    _nbrids[ndx] = rhs._nbrids[ndx];
	  }
      }
    return *this;
  }
};

std::ostream& operator<<(std::ostream& os, const s3simplex& sx);

typedef std::map<bignum_t, s3simplex> S3SimplexStore;
typedef std::map<bignum_t, s3simplex>::const_iterator S3SSCIter;
typedef std::map<bignum_t, s3simplex>::iterator S3SSIter;

fixnum_t pick_random_neighbor( void );

class ssHausdorf3p1
{
public:
  ssHausdorf3p1( const char* infilename, 
		 const char* outfilename,
		 fixnum_t plusminus );
  ~ssHausdorf3p1( void );
  
  void compute(void);
    
private:
  // input and output data files
  std::string _ifilename;
  std::string _ofilename;

  // following parameters pertain to the hausdorf dimension calculation
  fixnum_t _plusminus;

  // following loaded from .s3sx2p1 file
  fixnum_t _numT;
  std::string _bctype;
  std::string _stopology;
  fixnum_t _initVol;
  bignum_t _N3_SL;
  double _kappa0;
  double _delta;
  double _kappa4; 
  
  inline fixnum_t offset( fixnum_t i0, fixnum_t k )
  {
    return( (i0 + (k % _numT) + _numT) % _numT );
  }

  // stores the spatial 3-simplex data
  S3SimplexStore _s3sxstore;

  // # of neighbors at a distance of d
  std::map<fixnum_t,bignum_t> _dVsNumd;

  // slice with the largest volume
  fixnum_t _maxVolSlice;

  // implementation methods

  void _loadDataFromFile( void );
  void _parseDataLine( std::string& line, bignum_t& id, s3simplex& sx );
  void _parseParameterLine( std::string& line );
  void _getS3SimplicesInSlice( fixnum_t tslice, std::vector<bignum_t>& sxs );
  void _computeMaxVolSlice( void );

  // prohibit default initialization and copying
  ssHausdorf3p1( void );
  ssHausdorf3p1( const ssHausdorf3p1& );
  ssHausdorf3p1& operator=( const ssHausdorf3p1& );
};

#endif // __ssHausdorf_3p1SI__hpp__
