#include "stsDiffusion3p1SI.hpp"
#include <algorithm>
#include <fstream>
#include <iostream>
#include <cstdlib>

int main( int argc, char** argv )
{
  if( argc != 5 )
    {
      std::cerr << 
	"usage: stsDiffusion3p1SI infilelist sigma-min sigma-max numtrials"
		<< std::endl;
      exit(-1);
    }

  std::ifstream fin( argv[1] );
  if( fin.is_open() )
    {
      std::string infilename;
      while( std::getline( fin, infilename ) )
	{
	  std::string infilecopy = infilename;
	  size_t pos = infilecopy.find(".4sx3p1");
	  std::string outfilename = infilecopy.substr(0, pos) + "-" + 
	    argv[2] + "-" + argv[3] + "-" + argv[4] + "-diffusion-SI.sts3p1"; 

	  stsDiffusion3p1 sd( infilename.c_str(),     // infile
			      outfilename.c_str(),    // outfile
			      atol(argv[2]),          // sigmin
			      atol(argv[3]),          // sigmax
			      atol(argv[4]) );        // numtrials
	  sd.diffusion();
	  std::cout << "finished processing " << outfilename << std::endl;
	}
      fin.close();
    }
}

