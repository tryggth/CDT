#ifndef __sss_diffusion3p1SI__hpp__
#define __sss_diffusion3p1SI__hpp__

/*
 * sss3p1 --- computes the spatial slice spectral dimension for 2p1 spacetimes. 
 * The spectral dimension is computed for the maximum volume slice +- the offset
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

  s3simplex( fixnum_t tslice, bignum_t p0, bignum_t p1, bignum_t p2, bignum_t p3, 
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

struct kt_key_t
{
  bignum_t id1; bignum_t id2; fixnum_t sig;
  kt_key_t( void ) : id1(0), id2(0), sig(0) 
  { 
  }
  kt_key_t( bignum_t i1, bignum_t i2, fixnum_t s ) : id1(i1), id2(i2), sig(s)
  {
  }
};

inline bool operator==( const kt_key_t& k1, const kt_key_t& k2 )
{
  return ( (k1.id1==k2.id1) && (k1.id2==k2.id2) && (k1.sig==k2.sig) );
}

inline std::size_t hash_value( const kt_key_t& k )
{
  std::size_t seed = 0;
  boost::hash_combine( seed, k.id1 );
  boost::hash_combine( seed, k.id2 );
  boost::hash_combine( seed, k.sig );
  return seed;
}

typedef real_t kt_val_t;

typedef 
boost::unordered_map< kt_key_t, kt_val_t, boost::hash<kt_key_t> > 
kt_store_t;

fixnum_t pick_random_neighbor( void );

class sssDiffusion3p1
{
public:
  sssDiffusion3p1( const char* infilename, 
		   const char* outfilename,
		   fixnum_t sigmin,
		   fixnum_t sigmax, 
		   fixnum_t offset,
		   fixnum_t numtrials );
  ~sssDiffusion3p1( void );
  
  void diffusion(void);
    
private:
  // input and output data files
  std::string _ifilename;
  std::string _ofilename;

  // following parameters pertain to the spectral dimension calculation
  fixnum_t _sigmin;
  fixnum_t _sigmax;
  fixnum_t _numtrials;
  fixnum_t _offset;

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
  // memoize the KT(j,k,s) values that have been computed
  kt_store_t _KTstore;

  // stores the computed sigma vs retprob values
  std::map<fixnum_t,real_t> _sigVsRetProb;

  // slice with the largest volume
  fixnum_t _maxVolSlice;

  // implementation methods
  real_t KT( bignum_t j, bignum_t i0, fixnum_t sigma );
  void _loadDataFromFile( void );
  void _parseDataLine( std::string& line, bignum_t& id, s3simplex& sx );
  void _parseParameterLine( std::string& line );
  void _getS3SimplicesInSlice( fixnum_t tslice, std::vector<bignum_t>& sxs );
  void _computeMaxVolSlice( void );

  // prohibit default initialization and copying
  sssDiffusion3p1( void );
  sssDiffusion3p1( const sssDiffusion3p1& );
  sssDiffusion3p1& operator=( const sssDiffusion3p1& );
};

#endif // __sss_diffusion3p1SI__hpp__
