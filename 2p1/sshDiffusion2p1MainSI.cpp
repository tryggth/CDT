#include "sshDiffusion2p1SI.hpp"
#include <algorithm>
#include <fstream>
#include <iostream>
#include <cstdlib>

int main( int argc, char** argv )
{
  if( argc != 3 )
    {
      std::cerr << 
	"usage: sshDiffusion2p1SI infilelist numtrials"
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
	  size_t pos = infilecopy.find(".s2sx2p1");
	  std::string outfilename = infilecopy.substr(0, pos) + "-" + 
	    argv[2] + "-diffusion-SI.ssh2p1"; 

	  sshDiffusion2p1 sd( infilename.c_str(),     // infile
			      outfilename.c_str(),    // outfile
			      atol(argv[2]) );        // numtrials
	  sd.hausdorff();
	  std::cout << "finished processing " << outfilename << std::endl;
	}
      fin.close();
    }
}

