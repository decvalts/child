/**************************************************************************/
/**
**  @file mathutil.cpp
**  @brief Special math routines not in math libraries. Most or all
**         from Numerical Recipes in C by Press et al.
**
**  $Id: mathutil.cpp,v 1.7 2004-06-16 13:37:26 childcvs Exp $
*/
/**************************************************************************/

#include "mathutil.h"

#include <fstream>
#include <cstdlib>
#include <ctime>

#include "mt19937ar-cok.cpp"

#include "../tInputFile/tInputFile.h"
#include "../tListInputData/tListInputData.h"

/*********************************************************\
**  ran3
**
**  Random number generator from Numerical Recipes.
**  Returns a uniform random number between 0.0 and 1.0.
**  Set idum to any negative value to initialize or
**  reinitialize the sequence.
**
**  Parameters: idum - random seed
**
\*********************************************************/

tRand::tRand(long seed)
{
  init(seed);
}

tRand::tRand( tRand const &orig )
  : inext(orig.inext), inextp(orig.inextp)
{
  for( int i=1;i<=55;i++)
    ma[i] = orig.ma[i];
}
tRand::tRand( tInputFile const &infile )
{
  initFromFile( infile );
  // read previous state if necessary
  int opt;
  if ( (opt = infile.ReadItem( opt, "OPTREADINPUT" )) == OPTREADINPUT_PREVIOUS)
    tListInputDataRand( infile, *this );
}

void tRand::initFromFile( tInputFile const &infile )
{
  int seed;
  seed = infile.ReadItem( seed, "SEED" );
  init(seed);
}

void tRand::dumpToFile( std::ofstream& outFile ){
  for(size_t i=1; i<sizeof(ma)/sizeof(ma[0]); ++i)
    outFile << ma[i] << '\n';
  outFile << inext << '\n' << inextp << '\n';
}

void tRand::readFromFile( std::ifstream& inFile ){
  for(size_t i=1; i<sizeof(ma)/sizeof(ma[0]); ++i)
    inFile >> ma[i];
  inFile >> inext;
  inFile >> inextp;
}

int tRand::numberRecords() const {
  return sizeof(ma)/sizeof(ma[0])-1+2;
}

#define MBIG 1000000000
#define MSEED 161803398
#define MZ 0
#define FAC (1.0/MBIG)

void tRand::init(long seed)
{
  int i,ii,k;

  long mj=MSEED-(seed < 0 ? -seed : seed);
  mj %= MBIG;
  ma[55]=mj;
  long mk=1;
  for (i=1;i<=54;i++) {
    ii=(21*i) % 55;
    ma[ii]=mk;
    mk=mj-mk;
    if (mk < MZ) mk += MBIG;
    mj=ma[ii];
  }
  for (k=1;k<=4;k++)
    for (i=1;i<=55;i++) {
      ma[i] -= ma[1+(i+30) % 55];
      if (ma[i] < MZ) ma[i] += MBIG;
    }
  inext=0;
  inextp=31;
}


double tRand::ran3()
{
  if (++inext == 56) inext=1;
  if (++inextp == 56) inextp=1;
  long mj=ma[inext]-ma[inextp];
  if (mj < MZ) mj += MBIG;
  ma[inext]=mj;
  return mj*FAC;
}
#undef MBIG
#undef MSEED
#undef MZ
#undef FAC

double tRand::ExpDev()
{
  double dum;
  do
    dum = ran3();
  while( dum == 0.0 );
  return -log( dum );
}

// Addition DAV 2016
// Random generator between interval, calls ran3()
double tRand::RandCustomInterval(double floatmin, double floatmax)
{
  double f = (double)ran3() / RAND_MAX;
  return floatmin + f * (floatmax - floatmin);
}

// Different implementation of generator
double tRand::RandRange2(int min, int max)
{
  std::srand(std::time(0));
  int randNum = std::rand() % (max - min + 1) + min;
  
  double ran_dbl = static_cast<double>(randNum);
  
  return ran_dbl;
}


















