#ifndef __sts_diffusion2p1SI__hpp__
#define __sts_diffusion2p1SI__hpp__

/*
 * Thu May  5 09:29:01 2011 
 * this is an exact simulation of the diffusion process; no random walks!
 * SI means Stalk Included
 */
#include <vector>
#include <map>
#include <string>
#include <boost/tuple/tuple.hpp>
#include <boost/unordered_map.hpp>
#include <boost/functional/hash.hpp>

typedef unsigned long long bignum_t;
typedef unsigned long fixnum_t;
typedef std::pair<fixnum_t,fixnum_t> fixnum_pair_t;
typedef double real_t;

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

inline std::ostream& operator<<( std::ostream& os, const kt_key_t& k )
{
  os << "(" << k.id1 << "," << k.id2 << "," << k.sig << ")";
  return os;
}

typedef real_t kt_val_t;

typedef 
boost::unordered_map< kt_key_t, kt_val_t, boost::hash<kt_key_t> > 
kt_store_t;

struct simplex3
{
  fixnum_t _type;
  fixnum_pair_t _tmlohi;
  bignum_t _points[4];
  bignum_t _nbrids[4];

  simplex3( void ) : _type(-1)
  {
    _tmlohi.first = _tmlohi.second = static_cast<fixnum_t>(-1);
    for( size_t ndx = 0; ndx < 4; ++ndx )
      _points[ndx] = _nbrids[ndx] = 0;
  }

  simplex3( fixnum_t type, fixnum_t tmlo, fixnum_t tmhi, 
	    bignum_t p0, bignum_t p1, bignum_t p2, bignum_t p3,
	    bignum_t n0, bignum_t n1, bignum_t n2, bignum_t n3 )
    : _type( type )
  {
    _tmlohi.first = tmlo; _tmlohi.second = tmhi;
    _points[0] = p0; _points[1] = p1; _points[2] = p2; _points[3] = p3;
    _nbrids[0] = n0; _nbrids[1] = n1; _nbrids[2] = n2; _nbrids[3] = n3;
  }

  simplex3( const simplex3& rhs ) : _type( rhs._type), _tmlohi( rhs._tmlohi )
  {
    for( size_t ndx = 0; ndx < 4; ++ndx )
      {
	_points[ndx] = rhs._points[ndx];
	_nbrids[ndx] = rhs._nbrids[ndx];
      }
  }

  simplex3& operator=( const simplex3& rhs )
  {
    if( this != &rhs )
      {
	_type = rhs._type;
	_tmlohi = rhs._tmlohi;
	for( size_t ndx = 0; ndx < 4; ++ndx )
	  {
	    _points[ndx] = rhs._points[ndx];
	    _nbrids[ndx] = rhs._nbrids[ndx];
	  }
      }
    return *this;
  }
};

typedef std::map<bignum_t, simplex3> Simplex3Store;
typedef std::map<bignum_t, simplex3>::const_iterator S3SConstIter;
typedef std::map<bignum_t, simplex3>::iterator S3SIter;

typedef std::vector<fixnum_pair_t> vector_of_fixnumpair_t;
typedef std::vector<fixnum_pair_t>::const_iterator lfxnp_constiter_t;
typedef std::vector<fixnum_pair_t>::iterator lfxnp_iter_t;


class stsDiffusion2p1
{
public:
  stsDiffusion2p1( const char* infilename, 
		   const char* outfilename, 
		   fixnum_t sigmin, 
		   fixnum_t sigmax,
		   fixnum_t numtrials ); 
  ~stsDiffusion2p1( void );

  inline fixnum_t maxVolume( void ) const
  {
    return _maxvolume;
  }

  void diffusion(void);


private:
  // input and output data files
  std::string _ifilename;
  std::string _ofilename;

  // following parameters pertain to the spectral dimension calculation
  fixnum_t _sigmin;
  fixnum_t _sigmax;
  fixnum_t _numtrials;
  // maximum volume sandwich
  std::pair<fixnum_t,fixnum_t> _maxVolSandwich;
  fixnum_t _maxvolume;

  // following loaded from .3sx2p1 file, none of which are used in the diffusion
  // process! we just use the _3sxstore
  fixnum_t _numT;
  std::string _bctype;
  std::string _stopology;
  fixnum_t _initVol;
  bignum_t _lastUsedPoint;
  bignum_t _lastUsedId;
  bignum_t _N0;
  bignum_t _N1_SL;
  bignum_t _N1_TL;
  bignum_t _N2_SL;
  bignum_t _N2_TL;
  bignum_t _N3_TL_31;
  bignum_t _N3_TL_22;
  double _eps; 
  double _k0;
  double _k3; 
  double _alpha;

  // the 3-simplex store
  Simplex3Store _3sxstore;
  // memoize the KT(j,k,s) values that have been computed
  kt_store_t _KTstore;
  // compute and store the return probability values
  std::map<fixnum_t, real_t> _sigmaVsRetProb;

  // implementation methods
  real_t KT( bignum_t j, bignum_t i0, fixnum_t sigma );
  bignum_t _countSimplicesInSandwich( fixnum_t tmlo, fixnum_t tmhi );
  void _loadSpacetimeDataFromFile( void );
  void _parseDataLine( std::string& line, bignum_t& id, simplex3& sx );
  void _parseParameterLine( std::string& line );

  void _computeMaxVolSandwich( void );
  std::vector<bignum_t> _getSimplicesInSandwich( fixnum_pair_t tslohi );

  // prohibit copying
  stsDiffusion2p1( void );
  stsDiffusion2p1( const stsDiffusion2p1& );
  stsDiffusion2p1& operator=( const stsDiffusion2p1& );
};

#endif // __sts_diffusion2p1SI__hpp__
