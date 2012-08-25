#include "vvCorrelator2p1SI.hpp"
#include <algorithm>
#include <fstream>
#include <iostream>
#include <cstdlib>

int main( int argc, char** argv )
{
  if( argc != 4 )
    {
      std::cerr << "usage: vvCorrelator2p1SI infilelist dmin dmax" << std::endl;
      exit(-1);
    }

  std::ifstream fin( argv[1] );
  if( fin.is_open() )
    {
      std::string infilename;
      while( std::getline( fin, infilename ) )
	{
	  std::string infilecopy = infilename;
	  size_t pos = infilecopy.find(".s2sx2p1");
	  std::string outfilename = infilecopy.substr(0, pos) + 
	    "-" + argv[2] + "-" + argv[3] + "-SI.vvc2p1"; 

	  vvCorrelator2p1 sd( infilename.c_str(),     // infile
			      outfilename.c_str(),    // outfile
			      atol(argv[2]),          // dmin
			      atol(argv[3]));         // dmax
	  sd.computeCorrelations();
	  std::cout << "finished processing " << outfilename << std::endl;
	}
      fin.close();
    }
}
