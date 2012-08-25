#include "ssCurvature2p1SI.hpp"
#include <algorithm>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <cstdlib>

int main( int argc, char** argv )
{
  if( argc != 3 )
    {
      std::cerr << 
	"usage: ssCurvature2p1SI infilelist offset"
		<< std::endl;
      exit(-1);
    }

  std::ifstream fin( argv[1] );
  if( fin.is_open() )
    {
      std::string offstr = argv[4];
      size_t pos = offstr.find("-");
      if( pos != std::string::npos )
	offstr.replace(pos,1,"m");
      std::string infilename;
      while( std::getline( fin, infilename ) )
	{
	  std::string infilecopy = infilename;
	  size_t pos = infilecopy.find(".s2sx2p1");
	  std::string outfilename = infilecopy.substr(0, pos) + "-" + 
	    argv[2] + "-SI.ssR2p1"; 

	  ssCurvature2p1 sd( infilename.c_str(),     // infile
			     outfilename.c_str(),    // outfile
			      atol(argv[2]) );        // offset
	  std::cout << std::setprecision(20) << sd.curvature() << std::endl;
	}
      fin.close();
    }
}

