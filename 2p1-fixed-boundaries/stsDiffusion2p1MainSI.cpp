#include "stsDiffusion2p1SI.hpp"
#include <algorithm>
#include <fstream>
#include <iostream>
#include <cstdlib>

int main( int argc, char** argv )
{
  if( argc != 5 )
    {
      std::cerr << 
	"usage: stsDiffusion2p1SI infilelist sigma-min sigma-max numtrials"
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
	  size_t pos = infilecopy.find(".3sx2p1");
	  std::string outfilename = infilecopy.substr(0, pos) + "-" + 
	    argv[2] + "-" + argv[3] + "-" + argv[4] + "-diffusion-SI.sts2p1"; 

	  stsDiffusion2p1 sd( infilename.c_str(),     // infile
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

