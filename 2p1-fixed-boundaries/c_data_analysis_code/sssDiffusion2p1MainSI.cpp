#include "sssDiffusion2p1SI.hpp"
#include <algorithm>
#include <fstream>
#include <iostream>
#include <cstdlib>

int main( int argc, char** argv )
{
  if( argc != 6 )
    {
      std::cerr << 
	"usage: sssDiffusion2p1SI infilelist sigma-min sigma-max offset numtrials"
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
	    argv[2] + "-" + argv[3] + "-" + offstr + "-" + argv[5] + 
	    "-diffusion-SI.sss2p1"; 

	  sssDiffusion2p1 sd( infilename.c_str(),     // infile
			      outfilename.c_str(),    // outfile
			      atol(argv[2]),          // sigmin
			      atol(argv[3]),          // sigmax
			      atol(argv[4]),          // offset
			      atol(argv[5]) );        // numtrials
	  sd.diffusion();
	  std::cout << "finished processing " << outfilename << std::endl;
	}
      fin.close();
    }
}

