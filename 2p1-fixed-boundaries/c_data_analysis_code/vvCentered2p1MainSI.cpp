#include "vvCentered2p1SI.hpp"
#include <algorithm>
#include <fstream>
#include <iostream>
#include <cstdlib>

int main( int argc, char** argv )
{
  if( argc != 2 )
    {
      std::cerr << "usage: vvCentered2p1SI infilelist" << std::endl;
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
	    "-SI.cent2p1"; 

	  vvCentered2p1 sd( infilename.c_str(),     // infile
			    outfilename.c_str() );    // outfile
	  sd.center();
	  std::cout << "finished processing " << outfilename << std::endl;
	}
      fin.close();
    }
}
