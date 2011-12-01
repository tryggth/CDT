/*
 * ssh3p1SI.cpp
 */
#include "ssh3p1SI.hpp"
#include <ctime>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
// boost libraries

ssHausdorf3p1::ssHausdorf3p1(const char* infilename, 
			     const char* outfilename,
			     fixnum_t plusminus)
  : _ifilename( infilename ), _ofilename( outfilename ),
    _plusminus(plusminus)
{
  _loadDataFromFile();
  srand( time(0) );
}

ssHausdorf3p1::~ssHausdorf3p1( void )
{
}

void ssHausdorf3p1::_loadDataFromFile( void )
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

void ssHausdorf3p1::_parseParameterLine( std::string& line )
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

void ssHausdorf3p1::_parseDataLine( std::string& line, 
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

void ssHausdorf3p1::_getS3SimplicesInSlice( fixnum_t tslice, 
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

void ssHausdorf3p1::_computeMaxVolSlice( void )
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

void ssHausdorf3p1::compute(void)
{
  std::ofstream fout( _ofilename.c_str() );
  _computeMaxVolSlice();
  std::vector<bignum_t> sxInSlice;

  // go from -_offset to +_offset, in steps of 1 (so total 2*_offset+1 slices)
  for( fixnum_t pm = -1*_plusminus; pm <= _plusminus; ++pm )
    {
      fixnum_t ts = offset( _maxVolSlice, pm );
      _getS3SimplicesInSlice( ts, sxInSlice );
      double rAverage = 0.0;
      for( fixnum_t ntrial = 0; ntrial <= sxInSlice.size()/50; ++ntrial )
	{
	  unsigned long rndndx = random() % sxInSlice.size();
	  bignum_t sxid = sxInSlice[rndndx];
	  s3simplex& sx = _s3sxstore[sxid];
	  std::set<bignum_t> visited;
	  typedef std::set<bignum_t>::const_iterator visited_citer;
	  typedef std::set<bignum_t>::iterator visited_iter;
	  _dVsNumd.clear();
	  fixnum_t d = 1;
	  _dVsNumd[d] = 4;
	  for( size_t ndx = 0; ndx < 4; ++ndx )
	    {
	      visited.insert( sx._nbrids[ndx] );  // visited <- sxid.neighbors
	    }
	  while( visited.size() < sxInSlice.size() )
	    {
	      d++;
	      bignum_t num = 0;
	      std::set<bignum_t> newvis;
	      for( visited_citer c = visited.begin(); c != visited.end(); ++c )
		{
		  bignum_t currid = *c;
		  s3simplex& currsx = _s3sxstore[currid];
		  for( size_t ndx = 0; ndx < 4; ++ndx )
		    {
		      bignum_t nbrid = currsx._nbrids[ndx];
		      if( (visited.find( nbrid ) == visited.end()) &&
			  (newvis.find( nbrid ) == newvis.end()) )
			{
			  newvis.insert( nbrid ); // mark this as visited
			  num++;
			}
		    }
		}
	      // merge newvis into visited
	      for( visited_citer nv = newvis.begin(); nv != newvis.end(); ++nv )
		{
		  visited.insert( *nv );
		}
	      _dVsNumd[d] = num;
	    }
	  fixnum_t dmax = d;
	  bignum_t rTimesNr = 0;
	  for( fixnum_t r = 1; r <= dmax; ++r )
	    {
	      rTimesNr += r * _dVsNumd[r];
	    }
	  rAverage += (rTimesNr / sxInSlice.size() );
	}
      rAverage /= (sxInSlice.size()/50);
      fout << std::setprecision(6) << std::fixed << rAverage << "\t"
	   << sxInSlice.size() << std::endl;
    }
  fout.close();
}


