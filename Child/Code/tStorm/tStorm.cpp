/**************************************************************************/
/**
**  @file tStorm.cpp
**  @brief Functions for class tStorm.
**
**  A tStorm object generates random storms assuming an exponential
**  distribution of rainfall intensity, storm duration, and time to the
**  next storm. It is essentially an implementation of the model of
**  P. Eagleson, 1978b, Water Resources Research. Its services include
**  reading the necessary parameters from a tInputFile, generating a new
**  storm, and reporting its various values.
**
**  $Id: tStorm.cpp,v 1.36 2004-06-16 13:37:42 childcvs Exp $
*/
/**************************************************************************/


#include <math.h>
#include <string.h>
#include "../Mathutil/mathutil.h"
#include <iostream>
#include <fstream>

#include "tStorm.h"

/**************************************************************************\
**
**  tStorm::tStorm:  Constructor for storms. The default constructor
**                   assigns a value of unity to storm depth, duration,
**                   and interstorm duration. (Note:
**                   this constructor does not allow option for sinusoidal
**                   variation in means).
**
\**************************************************************************/
tStorm::tStorm( bool optvar )
  :
  rand(0),
  p(1.0),
  stdur(1.0),
  istdur(1.0),
  endtm(1.0e9),
  optVariable(optvar)
{
   //srand( 0 );
}


/**************************************************************************\
**
**  tStorm::tStorm
**
**  Alternative constructor that reads parameters directly from a tInputFile
**  object. Reads option for variable storms (normally this is "yes"---that's
**  the point of these objects---but a user may wish to switch off variation
**  as a test), mean values for rainfall intensity, duration, and interstorm
**  period, and a random seed to initialize the random number generator.
**    Also reads an option for long-term sinusoidal variations in the mean
**  values, and if the option is selected, reads the relevant parameters.
**  Variables p0, stdur0, and istdur0 are the mean values of the means;
**  pdev, stdurdev, and istdurdev are the range of variation (e.g., if pMean
**  were to fluctuate between 5 and 10, p0 would be 7.5 and pdev 2.5).
**
**  Modifications:
**   - 3/00 initialization now includes creation of ".storm" file for storm
**     history (GT)
**
\**************************************************************************/
tStorm::tStorm( const tInputFile &infile, tRand *rand_,
        bool no_write_mode /* = false */ ) :
  rand(rand_)
{
   // Read + set parameters for storm intensity, duration, and spacing
   optVariable = infile.ReadBool( "OPTVAR" );
   optOroPrecip = infile.ReadBool( "OPTOROGRAPHICPRECIP", false ); // flag of orographic precipitation,if false run as normal without oro-rain.
   optSpatialPrecip = infile.ReadBool( "OPTSPATIALPRECIP", false ); // flag for spatial precip

   if( !no_write_mode )
     {
       infile.WarnObsoleteKeyword("PMEAN", "ST_PMEAN");
       infile.WarnObsoleteKeyword("STDUR", "ST_STDUR");
       infile.WarnObsoleteKeyword("ISTDUR", "ST_ISTDUR");
       infile.WarnObsoleteKeyword("OPTSINVAR", "ST_PMEAN");
       infile.WarnObsoleteKeyword("PERIOD", "ST_PMEAN");
       infile.WarnObsoleteKeyword("START_CYCLE_TIME", "ST_PMEAN");
       infile.WarnObsoleteKeyword("MAXPMEAN", "ST_PMEAN");
       infile.WarnObsoleteKeyword("MAXSTDURMN", "ST_STDUR");
       infile.WarnObsoleteKeyword("MAXISTDURMN", "ST_ISTDUR");
     }
   infile.ReadItem( p_ts, "ST_PMEAN");
   infile.ReadItem( stdur_ts, "ST_STDUR");
   infile.ReadItem( istdur_ts, "ST_ISTDUR");

   p = p_ts.calc(0.);
   stdur = stdur_ts.calc(0.);
   istdur = istdur_ts.calc(0.);

   endtm = infile.ReadItem( endtm, "RUNTIME" );

   if (optOroPrecip)
   {
          SpeedX = infile.ReadItem( SpeedX, "WINDSPEED_X" );
          SpeedY = infile.ReadItem( SpeedY, "WINDSPEED_Y" );
          source0 = infile.ReadItem( source0, "WATERBACKGROUND");
          tauc = infile.ReadItem( tauc, "TAU_C");
          tauf = infile.ReadItem( tauf, "TAU_F");
          avrge = infile.ReadItem( avrge, "SUBEDGE_ON");
          BasicP = infile.ReadItem( BasicP, "BASIC_P");
          initialqc = infile.ReadItem( initialqc, "INITIAL_QC");
          initialqs = infile.ReadItem( initialqs, "INITIAL_QS");
          subEgeNum = infile.ReadItem( subEgeNum, "SUBEDGE_NUM");
   }

   /// Spatial Precip - DV ///
   if (optSpatialPrecip)
   {
       // Read in the storm model code from the input file
       int cread = infile.ReadItem( cread , "SPATIAL_STORM_MODEL" );
       miStormType = IntToStormType( cread );

       // Read in the values of storm centre point, and radius. These will
       // be mean values if one of the random options is turned on
       stormcenterpoint_a = infile.ReadItem( stormcenterpoint_a, "STORMCENTER_X" );
       stormcenterpoint_b = infile.ReadItem( stormcenterpoint_b, "STORMCENTER_Y" );
       stormradius = infile.ReadItem( stormradius, "STORMRADIUS" );
   }

   double help;

   help = infile.ReadItem( help, "OPTREADINPUT" );
   if(help>0){
      help = infile.ReadItem( help, "INPUTTIME" );
      endtm += help;
   }

   // If variable storms used, create a file for writing them
   if( optVariable && !no_write_mode )
   {
      char fname[87];
#define THEEXT ".storm"
      infile.ReadItem( fname, sizeof(fname)-sizeof(THEEXT), "OUTFILENAME" );
      strcat( fname, THEEXT );
#undef THEEXT
      stormfile.open( fname );
      if( !stormfile.good() )
          std::cerr << "Warning: unable to create storm data file '"
            << fname << "'\n";
   }
 }
}

tStorm::tStorm( const tStorm& orig )
   :  stormfile(),
      rand(orig.rand),
      p_ts(orig.p_ts),
      stdur_ts(orig.stdur_ts),
      istdur_ts(orig.istdur_ts),
      p(orig.p),
      stdur(orig.stdur),
      istdur(orig.istdur),
      endtm(orig.endtm),
      optVariable(orig.optVariable)
{}
/**************************************************************************\
**
**  tStorm::TurnOnOutput, TurnOffOutput
**
**  Open output file so output will be directed to it,
**  or close output file so there won't be output.
**
**  SL, 10/2010
**
\**************************************************************************/
void tStorm::TurnOnOutput( const tInputFile& infile )
{
     // If variable storms used, create a file for writing them
  if( !stormfile.good() && optVariable )
   {
      char fname[87];
#define THEEXT ".storm"
      infile.ReadItem( fname, sizeof(fname)-sizeof(THEEXT), "OUTFILENAME" );
      strcat( fname, THEEXT );
#undef THEEXT
      stormfile.open( fname );
      if( !stormfile.good() )
          std::cerr << "Warning: unable to create storm data file '"
            << fname << "'\n";
   }
}

void tStorm::TurnOffOutput()
{
  if( stormfile.good() )
    stormfile.close();
}

/**************************************************************************\
**
**  tStorm::GenerateStorm
**
**  Generates a new storm by drawing new values of p, stdur, and istdur from
**  an exponential distribution and updating the random seed.
**    If the minp parameter is greater than zero, the function will keep
**  picking storms until it finds one with p>minp. The total elapsed time,
**  including the rejected storms and their associated interstorm periods,
**  is stored istdur.
**
**  Inputs:      minp -- minimum value of rainfall rate p (default 0)
**               mind -- minimum storm depth to produce runoff (default 0)
**               tm -- current time in simulation
**  Members updated:  p, stdur, istdur take on new random values (if optVar)
**                    pMean, stdurMean, istdurMean adjusted (if optSinVar)
**  Assumptions:  pMean > 0
**
**  Modifications:
**   - changed AND to OR in while loop, GT 5/99
**   - added to while loop an additional check to see if time has run out,
**     NG 2/00
**   - added output of time, storm intensity, & duration to ".storm" file
**     GT 3/00
**
\**************************************************************************/
void tStorm::GenerateStorm( double tm, tMesh< tLNode > *meshRef, double minp, double mind )  //////add mesh, for taking mositure of each node for orographic rainfall if it's turned on
{

   p = p_ts.calc(tm);
   stdur = stdur_ts.calc(tm);
   istdur = istdur_ts.calc(tm);

   // If option for random storms is on, pick a storm at random.
   // Keep picking and accumulating time until the storm depth or intensity
   // is greater than the minimum value needed to produce runoff.

   if( optVariable )
   {
      const double pMean = p;
      const double stdurMean = stdur;
      const double istdurMean = istdur;

      stdur = 0;
      istdur = 0;
      do
      {
         p = pMean*ExpDev();
         istdur += istdurMean*ExpDev() + stdur;
         stdur = stdurMean*ExpDev();
     if(0) { // Debug
       std::cout << "P " << p << "  ST " << stdur << "  IST " << istdur
             << "  DP " << p*stdur << " minp " << minp
             << " mind " <<mind << std::endl;
     }
      } while( (p<=minp || (p*stdur)<=mind) && (tm+istdur+stdur<endtm) );
      if( stormfile.good() )
    stormfile << istdur << " " << p << " " << stdur << std::endl;
   }

   /**************************************************************************\
   **
   ** Spatial model of rainfall distribution
   ** Based on a simple circular storm cell model
   ** The user must specify the x,y (a,b) coordinates of the eye of the storm
   ** Rainfall is currently set to be uniform throughought the cell area,
   **  and zero outside the cell.
   **
   \**************************************************************************/
   if ( optSpatialPrecip )
   {
       // Calculate which nodes should be wetted based on the spatial parameters
       // for storm centre and radius.
       tMesh< tLNode > *meshPtr;
       meshPtr = meshRef;
       tLNode *cn;
       tMesh< tLNode >::nodeListIter_t nodIter( meshPtr->getNodeList() );
       volatile double cnX;
       volatile double cnY;

       // Local variables for storm loaction and size defined here
       double center_a;
       double center_b;
       double this_radius;

       // Cases to set values of radius and x,y center should go here, with separate call out functions
       switch (miStormType)
       {
           // Different cases for generating storm points
           // Case 1: Fixed centre location
           case kStaticStormCell:
           {
                //i.e. they are all fixed
               center_a = stormcenterpoint_a;
               center_b = stormcenterpoint_b;
               this_radius = stormradius;
           }
           break;

           // Case 2: Chosen from a random distribution
           case kRandomStormCell:
           {
               do
               {
                    center_a = stormcenterpoint_a*ExpDev();
                    center_b = stormcenterpoint_b*ExpDev();
                    this_radius = stormradius*ExpDev();
           } while ()
           break;

           // Case 3: Chosen from a random weighted distribution
           case kWeightedRandomStormCell:
           {

           }
           break;
       }

       // Now set the wetted nodes based on the storm morphology values set above
       for ( cn = nodIter.FirstP(); nodIter.IsActive(); cn = nodIter.NextP() )
       {
           // Iterate over the nodes to calculate which ones get rainfall
           // and which ones don't.

           // get the node coords
           cnX = cn->getX();
           cnY = cn->getY();

           // Test if they are inside the radius of the storm
           if ( ((cnX - center_a) * (cnX - center_a) + (cnY - center_b) * (cnY - center_b) <= this_radius*this_radius) )
           {
               // We are in the storm cell, set the precip of the current node
               cn->setPreci( p );
           }
           else
           {
               // We are outside the storm cell, no rain here!
               cn->setPreci( 0.0 );
           }
       }
   }


    if( optOroPrecip )
   {
/**************************************************************************\
**
**  linear orographic precipitation model (Smith&Barstad, 2004)
**  modified by J Han
\**************************************************************************/
     #define kMaxSpokes 100
     tMesh< tLNode > *meshPtr;
       meshPtr = meshRef;
       tLNode *cn;  // pointer to the current node
       tEdge *flowedg, *flowdgeOrg, *flowedg1min, *flowedg1max, *flowedg2min, *flowedg2max;
       tMesh< tLNode >::nodeListIter_t nodIter( meshPtr->getNodeList() );
       double dhdx, dhdy, slopeOrg, angleV, angleV1max, angleV2max, angleV1min, angleV2min;
       double dx1min, dy1min, dz1min, dx1max, dy1max, dz1max, dx2min, dy2min, dz2min, dx2max, dy2max, dz2max;
       double ETS;
       double e=2.718281828;
       double molecularRation=0.622;
       double mixRation, specHum;
       double press;
       double density;
       double precipWeight;
       int ctr;
       double AngleWind, AngleWindOpp;
       double dest1qc, dest2qc, dest1qs, dest2qs, dest1L, dest2L;
       double dest1qcex, dest2qcex, dest1qsex, dest2qsex;
       const int nActiveNodes = meshPtr->getNodeList()->getActiveSize(); // # active nodes
       const int nnodes = meshPtr->getNodeList()->getSize(); // total # nodes
       tLNode **cnWind = new tLNode *[nActiveNodes];  // prt to node list ordered in wind direction
       double *plengthWind = new double [nActiveNodes]; // prt to the list of length of each point in wind direction


       int numNode;
       double minX, maxX, minY, maxY;
       double orignX, orignY;
       double k, a, b, c;  // line ax+by+c=0 which go through each point in the mesh

       /////////////////////////////////find the min, max coordinary of X and Y/////////////
       numNode=0;
        for( cn = nodIter.FirstP(); nodIter.IsActive();
        cn = nodIter.NextP() )
       {
            if (numNode == 0)
            {
                minX=cn->getX();
                maxX=cn->getX();
                minY=cn->getY();
                maxY=cn->getY();
            }
           if (numNode != 0 )
           {
               if (cn->getX()<minX)
               {
                   minX=cn->getX();
               }

               if (cn->getX()>maxX)
               {
                   maxX=cn->getX();
               }

                if (cn->getY()<minY)
               {
                   minY=cn->getY();
               }

               if (cn->getY()>maxY)
               {
                   maxY=cn->getY();
               }
           }

           numNode++;
       }

       /////////////////////////////node ordered in wind direction////////////////////////////

       if (SpeedX>=0 && SpeedY>=0)  ////////set origin point depend on wind direction
       {
           orignX=minX;
           orignY=minY;
       }

       if (SpeedX>=0 && SpeedY<0)
       {
           orignX=minX;
           orignY=maxY;
       }

       if (SpeedX<0 && SpeedY>=0)
       {
           orignX=maxX;
           orignY=minY;
       }

       if (SpeedX<0 && SpeedY<0)
       {
           orignX=maxX;
           orignY=maxY;
       }                                   /////////////////////////////////////////////////////////////



            numNode=0;
        for( cn = nodIter.FirstP(); nodIter.IsActive();  /////////////////////node order in wind direction
        cn = nodIter.NextP() )
        {
            if (SpeedX==0)
            {
                a=0;
                b=-1;
            }
            else if (SpeedY==0)
            {
                a=1;
                b=0;
            }
            else
            {
                a=-SpeedX/SpeedY;
                b=-1;
            }


            c=-a*cn->getX()-b*cn->getY();


            if (numNode == 0)
            {
                cnWind[numNode]=cn;
            }

            if (numNode != 0 )
            {
                for (int numNode2=numNode; numNode2 >= 0; numNode2--)
                {
                    if ((numNode2-1) < 0)
                    {

                        if (numNode2 < numNode)
                        {
                            for (int numNode3=numNode; numNode3 > numNode2; numNode3--)
                            {
                                cnWind[numNode3]=cnWind[numNode3-1];
                            }
                        }
                        cnWind[numNode2] = cn;
                        break;

                    }

                    else if (PTLlength(orignX, orignY, a, b, c) >= PTLlength(orignX, orignY, a, b, -a*cnWind[numNode2-1]->getX()-b*cnWind[numNode2-1]->getY()))
                    {

                        if (numNode2 < numNode)
                        {
                            for (int numNode3=numNode; numNode3 > numNode2; numNode3--)
                            {
                                cnWind[numNode3]=cnWind[numNode3-1];
                            }
                        }
                        cnWind[numNode2] = cn;
                        break;
                    }
                }
            }


            numNode++;

       }


       //////////////////////////////wind direction/////////////////////////////////////////

       if (SpeedY<0 && SpeedX/sqrt( SpeedX*SpeedX + SpeedY*SpeedY )>=0)
       {
           AngleWind=2*3.1416-acos(SpeedX/sqrt( SpeedX*SpeedX + SpeedY*SpeedY ));
       }
           if (SpeedY>=0 && SpeedX/sqrt( SpeedX*SpeedX + SpeedY*SpeedY ) >= 0)
       {
           AngleWind=acos(SpeedX/sqrt( SpeedX*SpeedX + SpeedY*SpeedY ));
       }



           if (SpeedY<0 && SpeedX/sqrt( SpeedX*SpeedX + SpeedY*SpeedY )< 0)
       {
           AngleWind=2*3.1416-(3.1416-acos(fabs(SpeedX/sqrt( SpeedX*SpeedX + SpeedY*SpeedY ))));
       }
           if (SpeedY>=0 && SpeedX/sqrt( SpeedX*SpeedX + SpeedY*SpeedY ) < 0)
       {
           AngleWind=3.1416-acos(fabs(SpeedX/sqrt( SpeedX*SpeedX + SpeedY*SpeedY )));
       }


    for( cn = nodIter.FirstP(); nodIter.IsActive();
        cn = nodIter.NextP() )
   {
       ctr=0;
       angleV1max=NULL;
       angleV2max=NULL;
       angleV1min=-3.1416*2;
       angleV2min=3.1416*2;
      flowedg = cn->getEdg();
      flowdgeOrg = cn->getEdg();

      angleV=flowedg->CalAngle();
      angleV=angleV-AngleWind;
      if (angleV<0 )
      {
          angleV1min=angleV;

          flowedg1min=flowedg;

          angleV1max=angleV;
          flowedg1max=flowedg;
      }
      if (angleV>=0)
      {
          angleV2min=angleV;
          flowedg2min=flowedg;
          angleV2max=angleV;
          flowedg2max=flowedg;
      }

      slopeOrg=flowedg->CalcSlope();
      flowedg = flowedg->getCCWEdg();

      while (flowedg!=flowdgeOrg)
      {

          angleV=flowedg->CalAngle();
          angleV=angleV-AngleWind;
          if (angleV<0 && angleV>angleV1min)  ////////////////////////calculate min angle for direction1
         {
              angleV1min=angleV;
              flowedg1min=flowedg;
              if (angleV1max==NULL)
              {
                  angleV1max=angleV1min;
                  flowedg1max=flowedg;
              }
         }

           if (angleV<0 && angleV<angleV1max) /////////////////////////// calculate max angle for direction1
         {
              angleV1max=angleV;

              flowedg1max=flowedg;
         }

          if (angleV>=0 && angleV<angleV2min) ///////////////////////// calculate min angle for direction2
          {
              angleV2min=angleV;
              flowedg2min=flowedg;
                if (angleV2max==NULL)
              {
                  angleV2max=angleV2min;
                  flowedg2max=flowedg;
              }
          }

          if (angleV>=0 && angleV>angleV2max) /////////////////////////// calculate max angle for direction2
          {
              angleV2max=angleV;
              flowedg2max=flowedg;
          }




          if( ctr>kMaxSpokes ) // Make sure to prevent std::endless loops
            {
               std::cerr << "Mesh error: node " << cn->getID()
                    << " appears to be surrounded by closed boundary nodes"
                    << std::endl;
               ReportFatalError( "Bailing out of InitFlowDirs()" );
            }

            ctr++;
            flowedg = flowedg->getCCWEdg();


      }

        if (angleV1min==-3.1416*2)
        {
            flowedg1min=flowedg2max;
            angleV1min=angleV2max;
        }
        if (angleV2min==3.1416*2)

        {
            flowedg2min=flowedg1max;
            angleV2min=angleV1max;
        }



        dhdx=(flowedg1min->CalcDz()*flowedg2min->CalcDy()-flowedg2min->CalcDz()*flowedg1min->CalcDy())/(flowedg1min->CalcDx()*flowedg2min->CalcDy()-flowedg2min->CalcDx()*flowedg1min->CalcDy());
        dhdy=(flowedg1min->CalcDz()*flowedg2min->CalcDx()-flowedg2min->CalcDz()*flowedg1min->CalcDx())/(flowedg1min->CalcDy()*flowedg2min->CalcDx()-flowedg2min->CalcDy()*flowedg1min->CalcDx());
      ETS=6.112*pow(e,17.67*(25-0.0065*cn->getZ())/(243.5+25-0.0065*cn->getZ()));
      press=101325*pow((1-0.0000225577*cn->getZ()),5.25588);
      mixRation=molecularRation*ETS/(press-ETS);
      specHum=mixRation/(1+mixRation);
      density=101325/(287.05*(25+273.15-0.0065*cn->getZ()));
      precipWeight=density*specHum*(SpeedX*dhdx+SpeedY*dhdy);  /////// unit: kg/m3
      cn->setSource(source0+precipWeight*10000/1); /// precipitation unit change (column density/time-> depth/time): kg/(m3*s) to m/(yr) need more editing

      cn->setOroqc(0);
      cn->setOroqs(0);



   }


/////////////////////////////////////////////calculate qc , qs, and precipitation //////////////////////////

        AngleWindOpp=AngleWind+3.1416; /////////get the opposite direction of wind
        if (AngleWindOpp>=3.1416*2)
        {
            AngleWindOpp=AngleWindOpp-3.1416*2;
        }
       for (int numNode=0; numNode <= nActiveNodes-1; numNode++)
        {
                   ctr=0;
       angleV1max=NULL;
       angleV2max=NULL;
       angleV1min=-3.1416*2;
       angleV2min=3.1416*2;
      flowedg = cnWind[numNode]->getEdg();
      flowdgeOrg = cnWind[numNode]->getEdg();

      angleV=flowedg->CalAngle();
      angleV=angleV-AngleWindOpp;
      if (angleV<0 )
      {
          angleV1min=angleV;

          flowedg1min=flowedg;

          angleV1max=angleV;
          flowedg1max=flowedg;
      }
      if (angleV>=0)
      {
          angleV2min=angleV;
          flowedg2min=flowedg;
          angleV2max=angleV;
          flowedg2max=flowedg;
      }

      slopeOrg=flowedg->CalcSlope();
      flowedg = flowedg->getCCWEdg();

      while (flowedg!=flowdgeOrg)
      {

          angleV=flowedg->CalAngle();
          angleV=angleV-AngleWindOpp;
          if (angleV<0 && angleV>angleV1min)  ////////////////////////calculate min angle for direction1
         {
              angleV1min=angleV;
              flowedg1min=flowedg;
              if (angleV1max==NULL)
              {
                  angleV1max=angleV1min;
                  flowedg1max=flowedg;
              }
         }

           if (angleV<0 && angleV<angleV1max) /////////////////////////// calculate max angle for direction1
         {
              angleV1max=angleV;

              flowedg1max=flowedg;
         }

          if (angleV>=0 && angleV<angleV2min) ///////////////////////// calculate min angle for direction2
          {
              angleV2min=angleV;
              flowedg2min=flowedg;
                if (angleV2max==NULL)
              {
                  angleV2max=angleV2min;
                  flowedg2max=flowedg;
              }
          }

          if (angleV>=0 && angleV>angleV2max) /////////////////////////// calculate max angle for direction2
          {
              angleV2max=angleV;
              flowedg2max=flowedg;
          }




          if( ctr>kMaxSpokes ) // Make sure to prevent std::endless loops
            {
               std::cerr << "Mesh error: node " << cnWind[numNode]->getID()
                    << " appears to be surrounded by closed boundary nodes"
                    << std::endl;
               ReportFatalError( "Bailing out of InitFlowDirs()" );
            }

            ctr++;


            flowedg = flowedg->getCCWEdg();
      }

        if (angleV1min==-3.1416*2)
        {
            flowedg1min=flowedg2max;
            angleV1min=angleV2max;
        }
        if (angleV2min==3.1416*2)

        {
            flowedg2min=flowedg1max;
            angleV2min=angleV1max;
        }

        if ( !(flowedg1min->getdest()->isNonBoundary()) && !(flowedg2min->getdest()->isNonBoundary()))
         {
              cnWind[numNode]->setOroqc(initialqc);
              cnWind[numNode]->setOroqs(initialqs);
         }


        else
        {

            if (SpeedX==0)
            {
                a=0;
                b=-1;
            }
            else if (SpeedY==0)
            {
                a=1;
                b=0;
            }
            else
            {
                a=-SpeedX/SpeedY;
                b=-1;
            }



            if (avrge==0)         //overlook the sub edge calculation of qs and qc
        {

                dest1qc=((tLNode *)flowedg1min->getdest())->getOroqc()+(((tLNode *)flowedg1min->getdest())->getSource()-((tLNode *)flowedg1min->getdest())->getOroqc()/tauc)*(PTLlength(orignX, orignY, a, b, -a*cnWind[numNode]->getX()-b*cnWind[numNode]->getY())-PTLlength(orignX, orignY, a, b, -a*flowedg1min->getdest()->getX()-b*flowedg1min->getdest()->getY()))/sqrt(SpeedX*SpeedX+SpeedY*SpeedY);

            if (dest1qc<0)
            {
                  dest1qc=0;
                  dest1qs=((tLNode *)flowedg1min->getdest())->getOroqs()+((tLNode *)flowedg1min->getdest())->getOroqc()+(((tLNode *)flowedg1min->getdest())->getSource()-((tLNode *)flowedg1min->getdest())->getOroqs()/tauf)*(PTLlength(orignX, orignY, a, b, -a*cnWind[numNode]->getX()-b*cnWind[numNode]->getY())-PTLlength(orignX, orignY, a, b, -a*flowedg1min->getdest()->getX()-b*flowedg1min->getdest()->getY()))/sqrt(SpeedX*SpeedX+SpeedY*SpeedY);
            }
            else
            {
                dest1qs=((tLNode *)flowedg1min->getdest())->getOroqs()+(((tLNode *)flowedg1min->getdest())->getOroqc()/tauc-((tLNode *)flowedg1min->getdest())->getOroqs()/tauf)*(PTLlength(orignX, orignY, a, b, -a*cnWind[numNode]->getX()-b*cnWind[numNode]->getY())-PTLlength(orignX, orignY, a, b, -a*flowedg1min->getdest()->getX()-b*flowedg1min->getdest()->getY()))/sqrt(SpeedX*SpeedX+SpeedY*SpeedY);
            }
            if (dest1qs<0)
            {
                  dest1qs=0;
            }


                dest2qc=((tLNode *)flowedg2min->getdest())->getOroqc()+(((tLNode *)flowedg2min->getdest())->getSource()-((tLNode *)flowedg2min->getdest())->getOroqc()/tauc)*(PTLlength(orignX, orignY, a, b, -a*cnWind[numNode]->getX()-b*cnWind[numNode]->getY())-PTLlength(orignX, orignY, a, b, -a*flowedg2min->getdest()->getX()-b*flowedg2min->getdest()->getY()))/sqrt(SpeedX*SpeedX+SpeedY*SpeedY);

             if (dest2qc<0)
            {
                  dest2qc=0;
                  dest2qs=((tLNode *)flowedg2min->getdest())->getOroqs()+((tLNode *)flowedg2min->getdest())->getOroqc()+(((tLNode *)flowedg2min->getdest())->getSource()-((tLNode *)flowedg2min->getdest())->getOroqs()/tauf)*(PTLlength(orignX, orignY, a, b, -a*cnWind[numNode]->getX()-b*cnWind[numNode]->getY())-PTLlength(orignX, orignY, a, b, -a*flowedg2min->getdest()->getX()-b*flowedg2min->getdest()->getY()))/sqrt(SpeedX*SpeedX+SpeedY*SpeedY);
            }
            else
            {

                dest2qs=((tLNode *)flowedg2min->getdest())->getOroqs()+(((tLNode *)flowedg2min->getdest())->getOroqc()/tauc-((tLNode *)flowedg2min->getdest())->getOroqs()/tauf)*(PTLlength(orignX, orignY, a, b, -a*cnWind[numNode]->getX()-b*cnWind[numNode]->getY())-PTLlength(orignX, orignY, a, b, -a*flowedg2min->getdest()->getX()-b*flowedg2min->getdest()->getY()))/sqrt(SpeedX*SpeedX+SpeedY*SpeedY);
            }
             if (dest2qs<0)
           {
                  dest2qs=0;
            }

        }
            if (avrge==1)                       // sub edge node calculation for qs and qc, and it might make  qs and qc more smooth
            {
                dest1qcex=((tLNode *)flowedg1min->getdest())->getOroqc();
                dest1qsex=((tLNode *)flowedg1min->getdest())->getOroqs();
                dest2qcex=((tLNode *)flowedg2min->getdest())->getOroqc();
                dest2qsex=((tLNode *)flowedg2min->getdest())->getOroqs();
                for (int loopn=0; loopn < subEgeNum; loopn++ )
                {
                    dest1qc=dest1qcex+(((tLNode *)flowedg1min->getdest())->getSource()-dest1qcex/tauc)*(PTLlength(orignX, orignY, a, b, -a*cnWind[numNode]->getX()-b*cnWind[numNode]->getY())-PTLlength(orignX, orignY, a, b, -a*flowedg1min->getdest()->getX()-b*flowedg1min->getdest()->getY()))/sqrt(SpeedX*SpeedX+SpeedY*SpeedY)/subEgeNum;

                        if (dest1qc<0)
                      {
                         dest1qc=0;
                          dest1qs=dest1qsex+dest1qcex+(((tLNode *)flowedg1min->getdest())->getSource()-dest1qsex/tauf)*(PTLlength(orignX, orignY, a, b, -a*cnWind[numNode]->getX()-b*cnWind[numNode]->getY())-PTLlength(orignX, orignY, a, b, -a*flowedg1min->getdest()->getX()-b*flowedg1min->getdest()->getY()))/sqrt(SpeedX*SpeedX+SpeedY*SpeedY)/subEgeNum;
                      }
                      else
                     {
                         dest1qs=dest1qsex+(dest1qcex/tauc-dest1qsex/tauf)*(PTLlength(orignX, orignY, a, b, -a*cnWind[numNode]->getX()-b*cnWind[numNode]->getY())-PTLlength(orignX, orignY, a, b, -a*flowedg1min->getdest()->getX()-b*flowedg1min->getdest()->getY()))/sqrt(SpeedX*SpeedX+SpeedY*SpeedY)/subEgeNum;
                     }
                     if (dest1qs<0)

                     {
                         dest1qs=0;
                     }


                        dest2qc=dest2qcex+(((tLNode *)flowedg2min->getdest())->getSource()-dest2qcex/tauc)*(PTLlength(orignX, orignY, a, b, -a*cnWind[numNode]->getX()-b*cnWind[numNode]->getY())-PTLlength(orignX, orignY, a, b, -a*flowedg2min->getdest()->getX()-b*flowedg2min->getdest()->getY()))/sqrt(SpeedX*SpeedX+SpeedY*SpeedY)/subEgeNum;

                     if (dest2qc<0)
                    {
                         dest2qc=0;
                         dest2qs=dest2qsex+dest2qcex+(((tLNode *)flowedg2min->getdest())->getSource()-dest2qsex/tauf)*(PTLlength(orignX, orignY, a, b, -a*cnWind[numNode]->getX()-b*cnWind[numNode]->getY())-PTLlength(orignX, orignY, a, b, -a*flowedg2min->getdest()->getX()-b*flowedg2min->getdest()->getY()))/sqrt(SpeedX*SpeedX+SpeedY*SpeedY)/subEgeNum;
                     }
                      else
                    {

                        dest2qs=dest2qsex+(dest2qcex/tauc-dest2qsex/tauf)*(PTLlength(orignX, orignY, a, b, -a*cnWind[numNode]->getX()-b*cnWind[numNode]->getY())-PTLlength(orignX, orignY, a, b, -a*flowedg2min->getdest()->getX()-b*flowedg2min->getdest()->getY()))/sqrt(SpeedX*SpeedX+SpeedY*SpeedY)/subEgeNum;
                    }
                      if (dest2qs<0)
                    {
                            dest2qs=0;
                    }

                    dest1qcex=dest1qc;
                    dest1qsex=dest1qs;
                    dest2qcex=dest2qc;
                    dest2qsex=dest2qs;

                }
            }




              dest1L=sqrt((cnWind[numNode]->getX()-flowedg1min->getdest()->getX())*(cnWind[numNode]->getX()-flowedg1min->getdest()->getX())+(cnWind[numNode]->getY()-flowedg1min->getdest()->getY())*(cnWind[numNode]->getY()-flowedg1min->getdest()->getY()));
              dest2L=sqrt((cnWind[numNode]->getX()-flowedg2min->getdest()->getX())*(cnWind[numNode]->getX()-flowedg2min->getdest()->getX())+(cnWind[numNode]->getY()-flowedg2min->getdest()->getY())*(cnWind[numNode]->getY()-flowedg2min->getdest()->getY()));



             if (dest1L < fabs(PTLlength(orignX, orignY, a, b, -a*cnWind[numNode]->getX()-b*cnWind[numNode]->getY())-PTLlength(orignX, orignY, a, b, -a*flowedg1min->getdest()->getX()-b*flowedg1min->getdest()->getY())))
              {
                  dest1L=fabs(PTLlength(orignX, orignY, a, b, -a*cnWind[numNode]->getX()-b*cnWind[numNode]->getY())-PTLlength(orignX, orignY, a, b, -a*flowedg1min->getdest()->getX()-b*flowedg1min->getdest()->getY()));
              }
              if (dest2L < fabs(PTLlength(orignX, orignY, a, b, -a*cnWind[numNode]->getX()-b*cnWind[numNode]->getY())-PTLlength(orignX, orignY, a, b, -a*flowedg2min->getdest()->getX()-b*flowedg2min->getdest()->getY())))
              {
                  dest2L=fabs(PTLlength(orignX, orignY, a, b, -a*cnWind[numNode]->getX()-b*cnWind[numNode]->getY())-PTLlength(orignX, orignY, a, b, -a*flowedg2min->getdest()->getX()-b*flowedg2min->getdest()->getY()));
              }


              dest1L=sqrt(dest1L*dest1L-(PTLlength(orignX, orignY, a, b, -a*cnWind[numNode]->getX()-b*cnWind[numNode]->getY())-PTLlength(orignX, orignY, a, b, -a*flowedg1min->getdest()->getX()-b*flowedg1min->getdest()->getY()))*(PTLlength(orignX, orignY, a, b, -a*cnWind[numNode]->getX()-b*cnWind[numNode]->getY())-PTLlength(orignX, orignY, a, b, -a*flowedg1min->getdest()->getX()-b*flowedg1min->getdest()->getY())));
              dest2L=sqrt(dest2L*dest2L-(PTLlength(orignX, orignY, a, b, -a*cnWind[numNode]->getX()-b*cnWind[numNode]->getY())-PTLlength(orignX, orignY, a, b, -a*flowedg2min->getdest()->getX()-b*flowedg2min->getdest()->getY()))*(PTLlength(orignX, orignY, a, b, -a*cnWind[numNode]->getX()-b*cnWind[numNode]->getY())-PTLlength(orignX, orignY, a, b, -a*flowedg2min->getdest()->getX()-b*flowedg2min->getdest()->getY())));

            if ( !(flowedg1min->getdest()->isNonBoundary()) && flowedg2min->getdest()->isNonBoundary())
             {
                  cnWind[numNode]->setOroqc(dest2qc);
                  cnWind[numNode]->setOroqs(dest2qs);
             }
             else if ( flowedg1min->getdest()->isNonBoundary() && !(flowedg2min->getdest()->isNonBoundary()))
             {
                  cnWind[numNode]->setOroqc(dest1qc);
                  cnWind[numNode]->setOroqs(dest1qs);
             }
             else
             {
              cnWind[numNode]->setOroqc(dest1L/(dest1L+dest2L)*(dest2qc-dest1qc)+dest1qc);
              cnWind[numNode]->setOroqs(dest1L/(dest1L+dest2L)*(dest2qs-dest1qs)+dest1qs);
             }

       }


             cnWind[numNode]->setPreci(cnWind[numNode]->getOroqs()/tauf);



        if (cnWind[numNode]->getOroqs()/tauf <= BasicP)
        {
            cnWind[numNode]->setPreci (BasicP);
        }


        }






      delete[] cnWind;
      delete[] plengthWind;

   }




}


/**************************************************************************\
**
**  tStorm::ExpDev:  Finds a random number with an exponential distribution
**                   (adapted from Numerical Recipes).
**
\**************************************************************************/
double tStorm::ExpDev() const
{
    double dum;

    do
        dum = rand->ran3();
    while( dum==0.0 );
    return -log(dum);
}


/**************************************************************************\
**
**  tStorm "get" routines: return various variables
**
\**************************************************************************/

double tStorm::getStormDuration() const { return stdur; }
double tStorm::interstormDur() const { return istdur; }
double tStorm::getRainrate() const { return p; }

/**************************************************************************\
**
**  tStorm "set" routines: set various variables
**
\**************************************************************************/
void tStorm::setRainrate( double pMeanNew )
{
  std::stringstream ss;
  ss << pMeanNew;
  p_ts.reconfigure( ss.str().c_str() );
}


/**************************************************************************\
**
**  GammaDev
**
**  Returns a random variable drawn from a Gamma distribution with parameter m.
**
**  (Note: not actually called; provided for future use).
**
\**************************************************************************/
double tStorm::GammaDev(double m) const
{
  double x, y,z, c,t,b,u,w,v;

  if (m<1)
    {
      c = 1/m;
      t = 0.07 + 0.75*sqrt(1-m);
      b = 1 + exp(-t)*m/t;
      bool accept = false;
      while (!accept)
        {
          u = rand->ran3();
          w = rand->ran3();
          v = b *u;
          if (v<=1)
            {
              x  = t * pow(v, c);
              accept = ((w<=((2-x)/(2+x))) || (w<=exp(-x)));
            }
          else
            {
              x = -log(c*t*(b-v));
              y = x/t;
              accept = (((w*(m + y - m*y)) <= 1) || (w<= ( pow(y, (m-1)))));
            }
        }
    }
  else
    {
      b = m-1;
      c = 3*m - 0.75;
      bool accept = false;
      while (!accept)
        {
          u = rand->ran3(); v = rand->ran3();
          w = u* ( 1-u);
          y = sqrt(c/w) * (u - 0.5);
          x = b + y;
          if ( x>= 0)
            {
              z = 64*( pow(w,3.))*v*v;
              accept = (z <= ( 1 - 2*y*y/x)) || ( log(z) <= (2*(b*log(x/b) - y)));
            }
        }
    }
  return x;
}

// length of point (x,y) to line (ax+by+c=0) hanjianwei
double tStorm::PTLlength(double x, double y, double a, double b, double c)
{
    return fabs(a*x+b*y+c)/sqrt(a*a+b*b);
}


/// Implementation of the IntToStormType function
tStorm::kStormType_t tStreamNet::IntToStormType( int c ){
  switch(c){
    case 1: return kStaticStormCell;
    case 2: return kRandomStormCell;
    case 3: return kWeightedRandomStormCell;   // Added DV

    default:
      std::cout << "You asked for spatial storm model number " << c
      << " but there is no such thing.\n"
      "Available models are:\n"
      " 1. Static Storm Cell, with specified storm radius (specified x, y coordinates of storm centre)\n"
      " 2. Randomly located storm cells, of given radius\n"
      " 3. Randomly located storm cells but specified mean radius and mean location\n";
      ReportFatalError( "Unrecognized storm model code.\n" );
  }
}
