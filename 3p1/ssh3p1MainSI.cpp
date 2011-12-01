#include "ssh3p1SI.hpp"
#include <algorithm>
#include <fstream>
#include <iostream>
#include <cstdlib>

int main( int argc, char** argv )
{
  if( argc != 3 )
    {
      std::cerr << "usage: ssh3p1SI infilelist plusminus" << std::endl;
      exit(-1);
    }

  std::ifstream fin( argv[1] );
  if( fin.is_open() )
    {
      std::string plusminus = argv[2];
      std::string infilename;
      while( std::getline( fin, infilename ) )
	{
	  std::string infilecopy = infilename;
	  size_t pos = infilecopy.find(".s3sx3p1");
	  std::string outfilename = infilecopy.substr(0, pos) + "-" + 
	    plusminus + "-" + "-SI.ssh3p1"; 
	  ssHausdorf3p1 hd( infilename.c_str(),     // infile
			    outfilename.c_str(),    // outfile
			    atol(argv[2]) );        // plusminus
	  hd.compute();
	  std::cout << "finished processing " << outfilename << std::endl;
	}
      fin.close();
    }
}
