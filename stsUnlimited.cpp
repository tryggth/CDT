/*
 * unlimited spectral dimension
 */
#include <cstdlib>
#include <iostream>
#include <map>
using namespace std;

typedef unsigned long bignum_t;
struct vertex
{
  bignum_t nbrids[3];

  vertex( void ) 
  {
    for( size_t ndx = 0; ndx < 3; ++ndx )
      nbrids[ndx] = 0;
  }

  vertex( bignum_t n0, bignum_t n1, bignum_t n2 )
  {
    nbrids[0] = n0; nbrids[1] = n1; nbrids[2] = n2;
  }

  vertex( const vertex& rhs ) 
  {
    for( size_t ndx = 0; ndx < 3; ++ndx )
      nbrids[ndx] = rhs.nbrids[ndx];
  }

  vertex& operator=( const vertex& rhs ) 
  {
    if( this != &rhs )
      {
	for( size_t ndx = 0; ndx < 3; ++ndx )
	  nbrids[ndx] = rhs.nbrids[ndx];
      }
    return *this;
  }
};

ostream& operator<<( ostream& os, const vertex& vtx )
{
  os << "(" << vtx.nbrids[0] << " " << vtx.nbrids[1] << " " << vtx.nbrids[2] << ")";
  return os;
}

typedef map<bignum_t, vertex> triangulation_t;
typedef map<bignum_t, vertex>::const_iterator citer_t;
typedef map<bignum_t, vertex>::iterator iter_t;

triangulation_t g_triangulation;

long double g_retprobs[56];
//unsigned long g_retprobs[56];

void initialize_retprobs( size_t i0 )
{
  for( size_t ndx = 0; ndx < 56; ++ndx )
    g_retprobs[ndx] = 0.0;
  g_retprobs[i0] = 1.0;
}

void update_retprobs( void )
{
  long double retprobs[56];
  for( size_t ndx = 0; ndx < 56; ++ndx ) 
    retprobs[ndx] = 0.0;
  for( size_t ndx = 0; ndx < 56; ++ndx )
    {
      if( g_retprobs[ndx] != 0 )
	for( size_t n2 = 0; n2 < 3; ++n2 )
	  retprobs[g_triangulation[ndx].nbrids[n2]] += (g_retprobs[ndx]/3.0);
    }
  for( size_t ndx = 0; ndx < 56; ++ndx ) 
    g_retprobs[ndx] = retprobs[ndx];
}

void show_retprobs( void )
{
  for( size_t ndx = 0; ndx < 56; ++ndx )
    cout << ndx << ":" << scientific << g_retprobs[ndx] << " ";
  cout << "\n";
}

void show_retprobsv2( bignum_t i0, size_t sigma )
{
  cout << sigma << "\t" << scientific << g_retprobs[i0] << endl;
}

void populate_triangulation( void )
{
  g_triangulation.insert( make_pair(1, vertex(2,14,44)));
  g_triangulation.insert( make_pair(2, vertex(15,3,1)));
  g_triangulation.insert( make_pair(3, vertex(2,4,46)));
  g_triangulation.insert( make_pair(4, vertex(3,5,17)));
  g_triangulation.insert( make_pair(5, vertex(6,4,48)));
  g_triangulation.insert( make_pair(6, vertex(7,19,5)));
  g_triangulation.insert( make_pair(7, vertex(8,6,50)));
  g_triangulation.insert( make_pair(8, vertex(7,9,21)));
  g_triangulation.insert( make_pair(9, vertex(8,10,52)));
  g_triangulation.insert( make_pair(10, vertex(11,23,9)));
  g_triangulation.insert( make_pair(11, vertex(12,10,54)));
  g_triangulation.insert( make_pair(12, vertex(11,13,25)));
  g_triangulation.insert( make_pair(13, vertex(12,14,56)));
  g_triangulation.insert( make_pair(14, vertex(1,27,13)));
  g_triangulation.insert( make_pair(15, vertex(2,16,28)));
  g_triangulation.insert( make_pair(16, vertex(15,17,29)));
  g_triangulation.insert( make_pair(17, vertex(18,4,16)));
  g_triangulation.insert( make_pair(18, vertex(17,19,31)));
  g_triangulation.insert( make_pair(19, vertex(20,6,18)));
  g_triangulation.insert( make_pair(20, vertex(21,33,19)));
  g_triangulation.insert( make_pair(21, vertex(8,22,20)));
  g_triangulation.insert( make_pair(22, vertex(21,35,23)));
  g_triangulation.insert( make_pair(23, vertex(10,24,22)));
  g_triangulation.insert( make_pair(24, vertex(23,37,25)));
  g_triangulation.insert( make_pair(25, vertex(26,12,24)));
  g_triangulation.insert( make_pair(26, vertex(27,25,39)));
  g_triangulation.insert( make_pair(27, vertex(14,26,28)));
  g_triangulation.insert( make_pair(28, vertex(27,41,15)));
  g_triangulation.insert( make_pair(29, vertex(42,16,30)));
  g_triangulation.insert( make_pair(30, vertex(29,31,43)));
  g_triangulation.insert( make_pair(31, vertex(18,30,32)));
  g_triangulation.insert( make_pair(32, vertex(31,33,45)));
  g_triangulation.insert( make_pair(33, vertex(20,32,34)));
  g_triangulation.insert( make_pair(34, vertex(33,35,47)));
  g_triangulation.insert( make_pair(35, vertex(22,34,36)));
  g_triangulation.insert( make_pair(36, vertex(35,37,49)));
  g_triangulation.insert( make_pair(37, vertex(24,36,38)));
  g_triangulation.insert( make_pair(38, vertex(37,39,51)));
  g_triangulation.insert( make_pair(39, vertex(26,38,40)));
  g_triangulation.insert( make_pair(40, vertex(39,41,53)));
  g_triangulation.insert( make_pair(41, vertex(28,40,42)));
  g_triangulation.insert( make_pair(42, vertex(29,41,55)));
  g_triangulation.insert( make_pair(43, vertex(30,44,56)));
  g_triangulation.insert( make_pair(44, vertex(43,1,45)));
  g_triangulation.insert( make_pair(45, vertex(44,32,46)));
  g_triangulation.insert( make_pair(46, vertex(45,3,47)));
  g_triangulation.insert( make_pair(47, vertex(46,34,48)));
  g_triangulation.insert( make_pair(48, vertex(47,5,49)));
  g_triangulation.insert( make_pair(49, vertex(48,36,50)));
  g_triangulation.insert( make_pair(50, vertex(49,7,51)));
  g_triangulation.insert( make_pair(51, vertex(50,38,52)));
  g_triangulation.insert( make_pair(52, vertex(51,9,53)));
  g_triangulation.insert( make_pair(53, vertex(52,40,54)));
  g_triangulation.insert( make_pair(54, vertex(53,11,55)));
  g_triangulation.insert( make_pair(55, vertex(42,54,56)));
  g_triangulation.insert( make_pair(56, vertex(43,55,13)));
}

void show_triangulation( void )
{
  for( citer_t ci = g_triangulation.begin(); ci != g_triangulation.end(); ++ci )
    cout << ci->first << " --> " << ci->second << "\n";
}

int main( int argc, char** argv )
{
  populate_triangulation();
  initialize_retprobs(35);

  for( size_t sigma = 0; sigma < atoi(argv[1]); ++sigma )
    {
      show_retprobsv2(35,sigma);
      //show_retprobs();
      cout << endl;
      update_retprobs();
    }
}
