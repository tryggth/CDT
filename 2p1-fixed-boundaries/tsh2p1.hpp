#ifndef __tsh2p1__hpp__
#define __tsh2p1__hpp__

#include <vector>
#include <map>
#include <set>
#include <string>

typedef unsigned long long bignum_t;
typedef unsigned short fixnum_t;
typedef std::pair<fixnum_t,fixnum_t> fixnum_pair_t;

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

fixnum_t pick_random_neighbor( void );
bignum_t pick_random_simplex( const std::vector<bignum_t>& );

class ThickSliceHausdorffDim
{
public:
  ThickSliceHausdorffDim( const char* infilename, 
			  const char* outfilename, 
			  fixnum_t thresholdvol );
  ~ThickSliceHausdorffDim( void );

  void computeHausdorffDimension(void);
    
private:
  // input and output data files
  std::string _ifilename;
  std::string _ofilename;

  // following parameters pertain to the spectral dimension calculation
  std::vector<std::pair<fixnum_t,fixnum_t> > _stalkSandwiches;
  fixnum_t _thresholdVolume;

  // following loaded from .3sx2p1 file
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
  Simplex3Store _3sxstore;

  // implementation methods
  //  bignum_t _countSimplicesInSandwich( fixnum_t tmlo, fixnum_t tmhi );
  void _loadSpacetimeDataFromFile( void );
  void _parseDataLine( std::string& line, bignum_t& id, simplex3& sx );
  void _parseParameterLine( std::string& line );

  void _computeStalkSandwiches( void );
  void _getSimplicesInSandwich( fixnum_pair_t tslohi, 
				std::vector<bignum_t>& sxs );

  void _getNeighbors( const std::set<bignum_t>& oldnbors, 
		      std::set<bignum_t>& newnbors );

  // prohibit copying
  ThickSliceHausdorffDim( void );
  ThickSliceHausdorffDim( const ThickSliceHausdorffDim& );
  ThickSliceHausdorffDim& operator=( const ThickSliceHausdorffDim& );
};

#endif // __sts2p1v2__hpp__
