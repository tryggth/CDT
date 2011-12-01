#ifndef __ssh_diffusion2p1SI__hpp__
#define __ssh_diffusion2p1SI__hpp__

/*
 * ssh2p1 --- computes the spatial slice hausdorff dimension for 2p1 
 * spacetimes. The hausdorff dimension is computed for the maximum volume slice.
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

struct kt_key_t
{
  bignum_t id; fixnum_t sig;
  kt_key_t( void ) : id(0), sig(0) 
  { 
  }
  kt_key_t( bignum_t i, fixnum_t s ) : id(i), sig(s)
  {
  }
};

inline bool operator==( const kt_key_t& k1, const kt_key_t& k2 )
{
  return ( (k1.id==k2.id) && (k1.sig==k2.sig) );
}

inline std::size_t hash_value( const kt_key_t& k )
{
  std::size_t seed = 0;
  boost::hash_combine( seed, k.id );
  boost::hash_combine( seed, k.sig );
  return seed;
}

class sshDiffusion2p1
{
public:
  sshDiffusion2p1( const char* infilename, 
		   const char* outfilename,
		   fixnum_t numtrials );
  ~sshDiffusion2p1( void );
  
  void hausdorff(void);
    
private:
  // input and output data files
  std::string _ifilename;
  std::string _ofilename;

  fixnum_t _numtrials;

  // following loaded from .s2sx2p1 file
  fixnum_t _numT;
  std::string _bctype;
  std::string _stopology;
  fixnum_t _initVol;
  bignum_t _N2_SL;
  double _k0;
  double _k3; 

  // stores the spatial 2-simplex data
  S2SimplexStore _s2sxstore;

  std::set<bignum_t> _visited;
  boost::unordered_map< kt_key_t, bignum_t, boost::hash<kt_key_t> > _nbors;
  std::map<fixnum_t,real_t> _rVsN;

  // slice with the largest volume
  fixnum_t _maxVolSlice;
  // implementation methods
  void _loadDataFromFile( void );
  void _parseDataLine( std::string& line, bignum_t& id, s2simplex& sx );
  void _parseParameterLine( std::string& line );
  void _getS2SimplicesInSlice( fixnum_t tslice, std::vector<bignum_t>& sxs );
  void _computeMaxVolSlice( void );

  // prohibit default initialization and copying
  sshDiffusion2p1( void );
  sshDiffusion2p1( const sshDiffusion2p1& );
  sshDiffusion2p1& operator=( const sshDiffusion2p1& );
};

#endif // __ssh_diffusion2p1SI__hpp__
