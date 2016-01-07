//-*-c++-*- 

/**************************************************************************/
/**
**  @file tStorm.h
**  @brief Header for class tStorm
**
**  A tStorm object generates random storms assuming an exponential
**  distribution of rainfall intensity, storm duration, and time to the
**  next storm. It is essentially an implementation of the model of
**  P. Eagleson, 1978b, Water Resources Research. Its services include
**  reading the necessary parameters from a tInputFile, generating a new      
**  storm, and reporting its various values.
**    If you want to provide an option for NOT having storms vary
**  randomly, you can do so by setting optVariable to zero on initialization.
**    The GammaDev() function is provided for future reference; it is not
**  actually used in the current version.
**    At the user's option, the storm parameters can also be varied 
**  sinusoidally with time to simulate long-term climatic fluctuations.
**  (This is done by GenerateStorm).
**    tStorm objects could be easily modified (or inherited from) to use
**  different distributions. They can also be modified to create objects
**  for other random processes such as river flows, etc.
**
**  Created by Greg Tucker, November 1997.
**
**  Modifications:
**   - added data member "stormfile" to handle file containing history
**     of storm events
**
**  $Id: tStorm.h,v 1.31 2004-06-16 13:37:42 childcvs Exp $
*/
/**************************************************************************/

#ifndef TSTORM_H
#define TSTORM_H

#include "../tInputFile/tInputFile.h"
#include "../tTimeSeries/tTimeSeries.h"


#include "../Classes.h"
#include "../Definitions.h"
#include "../Mathutil/mathutil.h"
#include "../tArray/tArray.h"
#include "../tPtrList/tPtrList.h"
#include "../tList/tList.h"
#include "../tStreamNet/tStreamNet.h"
#include "../tLNode/tLNode.h"
#include "../MeshElements/meshElements.h"
#include "../tMesh/tMesh.h"
#include "../tInputFile/tInputFile.h"
#include "../globalFns.h"




#include <iosfwd>
#include <sstream>

class tStorm
{
public:
    tStorm( bool optVariable = true );
  tStorm( const tInputFile &, tRand *, bool no_write_mode = false );
  tStorm( const tStorm& );
    void GenerateStorm( double tm, tMesh< tLNode > *meshRef, double minp=0.0, double mind=0.0 ); //add tMesh< tLNode > &meshRef
    double getStormDuration() const;
    double interstormDur() const;
    double getRainrate() const;
    bool getOptVar() const;
	double PTLlength(double x, double y, double a, double b, double c);     // length of point (x,y) to line (ax+by+c=0) 
	bool optOroPrecip;  // Flag indicating whether orographic precipitation is used  
    void TurnOnOutput( const tInputFile& );
    void TurnOffOutput();
    inline void setRand( tRand* ptr ) {rand = ptr;}
    void setRainrate( double );

private:
    double ExpDev() const;
    double GammaDev(double) const;

    std::ofstream stormfile;// File containing history of storm events
    tRand *rand;       // Random number generator
    tTimeSeries p_ts;       // Rainfall intensity for the current storm
    tTimeSeries stdur_ts;   // Storm duration
    tTimeSeries istdur_ts;  // Time between storms
    double p;          // Actual rainfall intensity for the current storm
    double stdur;      // Actual storm duration
    double istdur;     // Actual time between storms
    double endtm;      // The end time of the run, just in case a big enough
                       // storm is never generated
    double SpeedX;     // windspeed in x direction 
	double SpeedY;     // windspeed in y direction 
	double source0;    // source water of orographic precipitation 
	double tauc;       // delay time of qc 
	double tauf;      //delay time of qs    
	double BasicP;    // the lowest value of precipitation 
	double initialqc;  // the initial value of qc near the bondary 
	double initialqs;  // the initial value of qs near the bondary 
	int avrge;        //subEdge lenght average qs and qc 
	int subEgeNum;    // sub Edge node number 
    bool optVariable;  // Flag indicating whether storms are random or not
};


inline bool tStorm::getOptVar() const {return optVariable;}


#endif
