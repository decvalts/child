/**************************************************************************\
**
**  tLNode.cpp
**
**  Functions for derived class tLNode and its member classes
**
**  Modifications:
**    - changes related to addition of Vegetation module (GT 1/2000):
**       - removed tSurface class & fn
**       - modified constructors (tLNode) to handle tau, tauc, vegCover
**    - fixed problem with layer initialization in copy constructor
**      (gt, 2/2000; see below)
** 
**  $Id: tLNode.cpp,v 1.95 2002-04-11 12:24:00 arnaud Exp $
\**************************************************************************/

#include <assert.h>
#include <math.h>
#include "../errors/errors.h"
#include "tLNode.h"
#define kBugTime 5000000

//Sets the total layer depth.  While updating depth, dgrade info is
//automatically updated to keep the same texture.  
//So if texture is going to change too - you MUST update dgrade info.
//(use setDgrade)
//which will automatically update total depth.
void tLayer::setDepth( double dep)
{
   if(dgrade.getSize()>0 && depth>0){      
      double sum=0;
      tArray<double> prop;
      prop.setSize(dgrade.getSize());
      int i=0;
      while(i<dgrade.getSize()){
         prop[i]=dgrade[i]/depth;
         dgrade[i]=prop[i]*dep;
         sum+=prop[i];
         i++;
         
      }
      if(fabs(sum-1)>0.01)
          ReportFatalError("Somewhere grain sizes got messed up");
      depth=dep;
   }
   else
       depth=dep;
}

void tLayer::setDgrade( int i, double size )
{
   assert( i<dgrade.getSize() );
   dgrade[i]=size;
   // Automatically update depth when dgrade is changed
   int j=0;
   double sum=0;
   while(j<dgrade.getSize()){
      sum += dgrade[j];
      j++;
   }
   depth=sum;
}
   
/*
tErode::tErode()                                                   //tErode
{
     //erodtype = 0;
   sedinput = zp = qs = qsp = qsin = qsinp = tau = dz = 0.0;
   nsmpts = 0;
     //cout << "  tErode()" << endl;
}

tErode::tErode( const tErode &orig )                                //tErode
{
  if( &orig != 0 )
   {
        //erodtype = orig.erodtype;
      sedinput = orig.sedinput;
      zp = orig.zp;
      qs = orig.qs;
      qsp = orig.qsp;
      qsin = orig.qsin;
      qsinp = orig.qsinp;
      tau = orig.tau;
      nsmpts = orig.nsmpts;
      dz = orig.dz;
   }
    //cout << "  tErode( orig )" << endl;
}

tErode::tErode( int numg, int nums )                                //tErode
{
   nsmpts = nums;
   sedinput = zp = qs = qsp = qsin = qsinp = tau = dz = 0.0;
     //cout << "  tErode( numg, nums )" << endl;
}

tErode::~tErode()                                                   //tErode
{
     //cout << "    ~tErode()" << endl;
}

//assignment
const tErode &tErode::operator=( const tErode &right )     //tErode
{
   if( &right != this )
   {
      sedinput = right.sedinput;             
      zp = right.zp;                         
      qs = right.qs;                         
      qsp = right.qsp;            
      qsin = right.qsin;                     
      qsinp = right.qsinp;                   
      tau = right.tau;                       
      nsmpts = right.nsmpts;                 
      dz = right.dz;                         
      smooth = right.smooth;                 
   }                                         
   return *this;                             
}
*/

tMeander::tMeander()                                              //tMeander
        : xyzd(4)
{
   meander = head = reachmember = 0;
   newx = newy = deltax = deltay = zoldright = zoldleft = bankrough = 0.0;
     //cout << "  tMeander()" << endl;
}

tMeander::tMeander( const tMeander &orig )                        //tMeander
        : xyzd( orig.xyzd )
{
   if( &orig != 0 )
   {
      meander = orig.meander;
      head = orig.head;
      reachmember = orig.reachmember;
      newx = orig.newx;
      newy = orig.newy;
      deltax = orig.deltax;
      deltay = orig.deltay;
      zoldright = orig.zoldright;
      zoldleft = orig.zoldleft;
      bankrough = orig.bankrough;
   }
     //cout << "  tMeander( orig )" << endl;
}

tMeander::tMeander( int state, double x, double y )                //tMeander
        : xyzd(4)
{
     //chanptr = 0;
   meander = state;
   newx = x;
   newy = y;
   head = reachmember = 0;
   deltax = deltay = zoldright = zoldleft = bankrough = 0.0;
     //cout << "  tMeander( state, x, y )" << endl;
}

tMeander::~tMeander()                                             //tMeander
{
     //cout << "    ~tMeander()" << endl;
}

//assignment
const tMeander &tMeander::operator=( const tMeander &right )     //tMeander
{
   if( &right != this )
   {
      meander = right.meander;
      head = right.head;
      reachmember = right.reachmember;
      newx = right.newx;
      newy = right.newy;
      deltax = right.deltax;
      deltay = right.deltay;
      zoldright = right.zoldright;
      zoldleft = right.zoldleft;
      bankrough = right.bankrough;
      xyzd = right.xyzd;
   }
   return *this;
}

tBedrock::tBedrock()                                             //tBedrock
{
   erodibility = 0;
     //cout << "  tBedrock()" << endl;
}

tBedrock::tBedrock( const tBedrock &orig )                       //tBedrock
{
   if( &orig != 0 )
   {
      erodibility = orig.erodibility;
   }
     //cout << "  tBedrock( orig )" << endl;
}

tBedrock::~tBedrock()                                            //tBedrock
{
     //cout << "    ~tBedrock()" << endl;
}

//assignment
const tBedrock &tBedrock::operator=( const tBedrock &right )     //tBedrock
{
   if( &right != this )
   {
      erodibility = right.erodibility;
   }
   return *this;
}


/***** Functions for class tRegolith **************************************/

/**************************************************************************\
**
**  Constructors:
**   - Default: initializes variables to zero.
**   - Input file: takes a tInputFile reference as an argument, and 
**              reads the default initial thickness from the file.
**   - Copy: copies all fields of tLNode and calls copy constructors for
**              embedded objects.
**
\**************************************************************************/

tRegolith::tRegolith()                                          //tRegolith
        : dgrade()
{
   thickness = 0;
     //cout << "  tRegolith()" << endl;
}

/*
tRegolith::tRegolith( tInputFile &infile )                     //tRegolith
  : dgrade( )
{
  int i;
  char add, name[20];
  double help, sum, numg;
}*/

tRegolith::tRegolith( const tRegolith &orig )                   //tRegolith
        : dgrade( orig.dgrade )
{
   if( &orig != 0 )
   {
      thickness = orig.thickness;
   }
   //cout << "  tRegolith( orig ) " << thickness << endl;
}


tRegolith::~tRegolith()                                         //tRegolith
{
     //cout << "    ~tRegolith()" << endl;
}

//assignment
const tRegolith &tRegolith::operator=( const tRegolith &right )     //tRegolith
{
   if( &right != this )
   {
      thickness = right.thickness;
      dgrade = right.dgrade;
   }
   return *this;
}

tChannel::tChannel()                                             //tChannel
{
   drarea = q = chanwidth = hydrwidth = channrough = hydrnrough =
       chandepth = hydrdepth = chanslope = hydrslope = diam = 
       mdFlowPathLength = 0;
   diam = kVeryHigh;
     //cout << "  tChannel()" << endl;
}

tChannel::tChannel( const tChannel &orig )                       //tChannel
        : /*erosion( orig.erosion ),*/ migration( orig.migration )
{
   if( &orig != 0 )
   {
      drarea = orig.drarea;
      q = orig.q;
      chanwidth = orig.chanwidth;
      chandepth = orig.chandepth;
      channrough = orig.channrough;
      chanslope = orig.chanslope;
      hydrwidth = orig.hydrwidth;
      hydrdepth = orig.hydrdepth;
      hydrnrough = orig.hydrnrough;
      hydrslope = orig.hydrslope;
      mdFlowPathLength = orig.mdFlowPathLength;
      diam = orig.diam;
   }
     //cout << "  tChannel( orig )" << endl;
}

tChannel::~tChannel()                                            //tChannel
{
     //cout << "    ~tChannel()" << endl;
}
//assignment
const tChannel &tChannel::operator=( const tChannel &right )     //tChannel
{
   if( &right != this )
   {
      //erosion = right.erosion;
       //migration = right.migration;
      drarea = right.drarea;
      q = right.q;
      chanwidth = right.chanwidth;
      chandepth = right.chandepth;
      channrough = right.channrough;
      chanslope = right.chanslope;
      hydrwidth = right.hydrwidth;
      hydrdepth = right.hydrdepth;
      hydrnrough = right.hydrnrough;
      hydrslope = right.hydrslope;
      mdFlowPathLength = right.mdFlowPathLength;
      diam = right.diam;
   }
   return *this;
}

/***** Functions for class tLNode *****************************************/

/**************************************************************************\
**  Constructors:
**   - Default: initializes several variables to zero and calls
**              constructors for embedded objects.
**   - Input file: takes a tInputFile reference as an argument, and 
**              passes it to the constructor(s) for embedded objects so
**              that they can read values for default parameters, etc.
**   - Copy: copies all fields of tLNode and calls copy constructors for
**              embedded objects.
**
**    Modifications:
**     - fixed an odd error with the layer initialization in the copy
**       constructor. The use of the embedded object initialization
**       syntax " : layerlist( orig.layerlist ) ", when compiled with
**       cxx, appeared to erroneously look for a tList constructor of
**       type tList( int ), which does not exist. Thus, the layerlist
**       was not duplicated but instead simply handed by default
**       memberwise copy, so that the new list pointed to the original 
**       one. Once the problem was discovered, the fix was simple:
**       use the assignment operator instead. This was probably the
**       source of the problem noted by nmg below. (gt, 2/2000)
**
\**************************************************************************/

// initialize static members numg and grade
int tLNode::numg = 0;
tArray<double> tLNode::grade = 1;
double tLNode::maxregdep = 1;
double tLNode::KRnew = 1.0;

tLNode::tLNode()                                                   //tLNode
        : tNode(), vegCover(), rock(), reg(), chan(), qsm(), qsinm(), 
          layerlist()
{
   //cout << "=>tLNode()" << endl;
   flood = 0;
   flowedge = 0;
   tracer = 0;
   dzdt = drdt = qs = qsin = uplift = 0.0;
   tau = tauc = 0;
}

tLNode::tLNode( tInputFile &infile )                               //tLNode
        : tNode(), vegCover(), rock(), reg(), chan(), qsm(), qsinm(), 
          layerlist()
{
   int i;
   char add[2], name[20];
   double help, extra, sum, sumbr;
   tLayer layhelp, niclay;
   tArray<double> dgradehelp;
   tArray<double> dgradebrhelp;

   tau = 0.0;
   tauc = infile.ReadItem( tauc, "TAUC" );
   
   //cout << "=>tLNode( infile )" << endl;
   flood = 0;
   flowedge = 0;
   tracer = 0;
   dzdt = drdt = qs = qsin = uplift = 0.0;
   numg = infile.ReadItem( numg, "NUMGRNSIZE" );
   grade.setSize( numg );
   maxregdep = infile.ReadItem( maxregdep, "MAXREGDEPTH" );
   KRnew = infile.ReadItem( KRnew, "KR" );
   if( KRnew<0.0 )
       ReportFatalError( "Erodibility factor KR must be positive." );
   
   i=0;
   add[0]='1';
   add[1]='\0';
   
   while ( i<numg ){
      // Reading in grain size diameter info 
      strcpy( name, "GRAINDIAM");
      strcat( name, add ); 
      help = infile.ReadItem( help, name);
      grade[i] = help;
      i++;
      add[0]++;
   }

   qsm.setSize( numg );
   qsinm.setSize( numg );      

   i = infile.ReadItem( i, "OPTREADLAYER" );

   if(i!=1){
      
      dgradehelp.setSize( numg );
      dgradebrhelp.setSize( numg );   
      sum = 0;
      sumbr = 0;
      i=0;
      add[0]='1';
      
      while ( i<numg ){
         // Reading in proportions for intital regolith and bedrock
         strcpy( name, "REGPROPORTION");
         strcat( name, add ); 
         help = infile.ReadItem( help, name);
         dgradehelp[i]=help;
         sum += help;
         strcpy( name, "BRPROPORTION");
         strcat( name, add ); 
         help = infile.ReadItem( help, name);
         dgradebrhelp[i]=help;
         sumbr += help;
         i++;
         add[0]++;
      }

      if(fabs(sum-1.0)>0.001)
          ReportFatalError("Problem with the regolith proportion of grain sizes in input file");
      if(fabs(sumbr-1.0)>0.001)
          ReportFatalError("Problem with the bedrock proportion of grain sizes in input file");
      
      help = infile.ReadItem( help, "REGINIT" );
      layhelp.setCtime( 0.0 );
      layhelp.setEtime( 0.0 );
      layhelp.setRtime( 0.0 );
      if( help > 0){
         // Make a bedrock and regolith layer, possibly two depending
         // on the depth of the regolith layer.  The code will decide
         // the total number of layers needed.  By default the regolith
         // layer(s) is/are put on top of the bedrock.
         // Bedrock layer items read in and set
         help = infile.ReadItem( help, "KB");
         if( help<0.0 )
             ReportFatalError( "Erodibility factor KB must be positive." );
         layhelp.setErody(help);
         layhelp.setSed(0);
         
         // in the regolith layer dgrade is saving
         // the proportion of grain size available from the bedrock
         
         layhelp.setDgradesize(numg);
         i=0;
         help = infile.ReadItem( help, "BEDROCKDEPTH");
         while(i<numg){
            layhelp.setDgrade(i, help*dgradebrhelp[i]);
            i++;
         }
         
         layerlist.insertAtBack( layhelp );
         
         // Regolith layer items read in and set
         layhelp.setSed(1);
         help = infile.ReadItem( help, "KR" );
         if( help<0.0 )
             ReportFatalError( "Erodibility factor KR must be positive." );
         layhelp.setErody(help);
         help = infile.ReadItem( help, "REGINIT");
         extra = 0;
         if(help > maxregdep){
            // too much regolith, create two layers the bottom layer is made here
            extra = help - maxregdep;
            //layhelp.setDepth(extra);
            //layhelp.setDgradesize(numg);
            i=0;
            while(i<numg){
               layhelp.setDgrade(i, extra*dgradehelp[i]);
               i++;
            }
            if(extra>10000) //TODO, bettter criteria for setting deep sed. flag
                layhelp.setRtime(-1);
            layerlist.insertAtFront( layhelp );
            // the top regolith layer is now made
            layhelp.setRtime(0);
            i=0;
            while(i<numg){
               layhelp.setDgrade(i, maxregdep*dgradehelp[i]);
               i++;
            }
            layerlist.insertAtFront( layhelp );
         }
         else{
            // create only one regolith layer
            i=0;
            while(i<numg){
               layhelp.setDgrade(i, help*dgradehelp[i]);
               i++;
            }
            
            layerlist.insertAtFront( layhelp );
         }
      }
      
      else{
         // no regolith, so by default everything is bedrock
         // Bedrock layer items set
         help = infile.ReadItem( help, "KB");
         if( help<0.0 )
             ReportFatalError( "Erodibility factor KB must be positive." );
         layhelp.setErody(help);
         layhelp.setSed(0);
         layhelp.setDgradesize(numg);
         i=0;
         help = infile.ReadItem( help, "BEDROCKDEPTH");
         while(i<numg){
            layhelp.setDgrade(i, help*dgradebrhelp[i]);
            i++;
         }
         layerlist.insertAtBack( layhelp );
      }
      
      //cout << layerlist.getSize() << " layers created " << endl;
   }

}

tLNode::tLNode( const tLNode &orig )                               //tLNode
        : tNode( orig ),
          vegCover( orig.vegCover ), rock( orig.rock ), 
          reg( orig.reg ), chan( orig.chan ), qsm( orig.qsm),
          qsinm( orig.qsinm )
    //Xlayerlist( orig.layerlist )
{

//Be aware that the copy constructor should be called using the 
//following syntax
//tLNode *newnode = new tLNode( *oldtLNode );
//NOT 
//tLNode newnode(*oldtLNode)
//The latter seems to call a default constructor which does not
//properly copy the layerlist
   flowedge = orig.flowedge;
   flood = orig.flood;
   tracer = orig.tracer;
   dzdt = orig.dzdt;
   drdt = orig.drdt;
   qs = orig.qs;
   qsin = orig.qsin;
   uplift = orig.uplift;
   tau = orig.tau;
   tauc = orig.tauc;
   layerlist = orig.layerlist;
   //cout << "=>tLNode( orig )" << endl;
}

tLNode::~tLNode()                                                  //tLNode
{
     //cout << "    ~tLNode()" << endl;
}

//assignment
const tLNode &tLNode::operator=( const tLNode &right )                  //tNode
{
   if( &right != this )
   {
      tNode::operator=( right );
      vegCover = right.vegCover;
      rock = right.rock;
      //surf = right.surf;
      reg = right.reg;
      chan = right.chan;
      flowedge = right.flowedge;
      flood = right.flood;
      tracer = right.tracer;
      dzdt = right.dzdt;
      drdt = right.drdt;
      qs = right.qs;
      qsin = right.qsin;
      tau = right.tau;
      tauc = right.tauc;
      uplift = right.uplift;
      qsm = right.qsm;
      qsinm = right.qsinm;
      layerlist = right.layerlist;
   }
   return *this;
}

//"get" and "set" functions; most simply return or set a data value, respectively:

const tBedrock &tLNode::getRock() const {return rock;}
//Xconst tSurface &tLNode::getSurf() const {return surf;}
const tRegolith &tLNode::getReg() const {return reg;}
const tChannel &tLNode::getChan() const {return chan;}

void tLNode::setRock( const tBedrock & val ) {rock = val;}
//Xvoid tLNode::setSurf( const tSurface & val ) {surf = val;}
void tLNode::setReg( const tRegolith & val ) {reg = val;}
void tLNode::setChan( const tChannel & val ) {chan = val;}

int tLNode::getFloodStatus() {   return flood; }

void tLNode::setFloodStatus( int status )
{
   flood = status;
}

tEdge * tLNode::getFlowEdg() 
{
   return flowedge;
}

void tLNode::setFlowEdg( tEdge * newflowedge )
{
   assert( newflowedge > 0 );  // Fails when passed an invalid edge
   flowedge = newflowedge;
}

void tLNode::setDrArea( double val ) {chan.drarea = val;}
void tLNode::AddDrArea( double val ) {chan.drarea += val;}
void tLNode::AddDischarge( double val ) {chan.q += val;}

tLNode * tLNode::getDownstrmNbr()
{
   //assert( flowedge!=0 );
   if( flowedge == 0 ) return 0;
   return (tLNode *)flowedge->getDestinationPtrNC();     
}

int tLNode::Meanders() const {return chan.migration.meander;}
void tLNode::setMeanderStatus( int val )
{chan.migration.meander = ( val == 0 || val == 1 ) ? val : 0;}


double tLNode::getHydrWidth() const {return chan.hydrwidth;}
double tLNode::getChanWidth() const {return chan.chanwidth;}
double tLNode::getHydrDepth() const {return chan.hydrdepth;}
double tLNode::getChanDepth() const {return chan.chandepth;}
double tLNode::getHydrRough() const {return chan.hydrnrough;}
double tLNode::getChanRough() const {return chan.channrough;}
double tLNode::getHydrSlope() const {return chan.hydrslope;}
double tLNode::getChanSlope() const {return chan.chanslope;}
//NG changed getDiam 02/1999
double tLNode::getDiam() const {
   int i;
   double di = 0;
   for(i=0; i<numg; i++){
      di+=grade[i]*getLayerDgrade(0,i)/getLayerDepth(0);
   }
   return di;
}
double tLNode::getBankRough() const {return chan.migration.bankrough;}

//TODO: suggest doing away with the zero test for performance enhancement
void tLNode::setHydrWidth( double val )  {chan.hydrwidth = ( val > 0 ) ? val : 0;}
void tLNode::setChanWidth( double val )  {chan.chanwidth = ( val > 0 ) ? val : 0;}
void tLNode::setHydrDepth( double val )  {chan.hydrdepth = ( val > 0 ) ? val : 0;}
void tLNode::setChanDepth( double val )  {chan.chandepth = ( val > 0 ) ? val : 0;}
void tLNode::setHydrRough( double val )  {chan.hydrnrough = ( val > 0 ) ? val : 0;}
void tLNode::setChanRough( double val )  {chan.channrough = ( val > 0 ) ? val : 0;}
void tLNode::setHydrSlope( double val )  {chan.hydrslope = ( val > 0 ) ? val : 0;}
void tLNode::setChanSlope( double val )  {chan.chanslope = ( val > 0 ) ? val : 0;}
void tLNode::setBankRough( double val )
{chan.migration.bankrough = ( val > 0 ) ? val : 0;}

double tLNode::getDrArea() const {return chan.drarea;}

tArray< double >
tLNode::getZOld() const
{
   tArray< double > rl(2);
   rl[0] = chan.migration.zoldright;
   rl[1] = chan.migration.zoldleft;
   return rl;
}

void tLNode::setZOld( double right, double left )
{
   chan.migration.zoldright = right;
   chan.migration.zoldleft = left;
}

tArray< double >                                                   //tNode
tLNode::getNew2DCoords() const
{
   tArray< double > xy(2);
   if( Meanders() )
   {
      xy[0] = chan.migration.newx;
      xy[1] = chan.migration.newy;
   }
   else
   {
      xy[0] = x;
      xy[1] = y;
   }
   return xy;
}

tArray< double >                                                   //tNode
tLNode::getNew3DCoords() const
{
   tArray< double > xyz(3);
   if( Meanders() )
   {
      xyz[0] = chan.migration.newx;
      xyz[1] = chan.migration.newy;
   }
   else
   {
      xyz[0] = x;
      xyz[1] = y;
   }
   xyz[2] = z;
   return xyz;
}

void tLNode::setNew2DCoords( double val1, double val2 )        //tNode
{
   chan.migration.newx = val1;
   chan.migration.newy = val2;
}

tArray< double > tLNode::
getLatDisplace() const
{
   tArray< double > xy(2);
   xy[0] = chan.migration.deltax;
   xy[1] = chan.migration.deltay;
   return xy;
}

void tLNode::setLatDisplace( double dx, double dy )
{
   chan.migration.deltax = dx;
   chan.migration.deltay = dy;
}

void tLNode::addLatDisplace( double dx, double dy )
{
   chan.migration.deltax += dx;
   chan.migration.deltay += dy;
}



void tLNode::setDischarge( double val ) {chan.q = ( val > 0 ) ? val : 0;}

void tLNode::RevertToOldCoords()
{
   chan.migration.newx = x;
   chan.migration.newy = y;
}

void tLNode::UpdateCoords()
{
   x = chan.migration.newx;
   y = chan.migration.newy;
}

// nb: if channel is integrated into node, change this
double tLNode::getQ()
{
   return chan.q;
}

/************************************************************************\
**  GetSlope: Computes and returns the slope of the node's flowedg, or
**  zero if the slope is less than zero.
**
**  The name of this function makes NG cringe.
**  TODO Change to CALCSLOPE!!!!
**
**  Assumptions: edge lengths up to date and nonzero, flowedg's up to
**    date.
**  Modified, 2/19/98: to return a reach slope over some distance
**   independent of the discretization for meandering nodes.
**   Nominally, we set that distance to 10 channel widths; we limit this
**   reach slope calculation to meandering nodes because they are the only
**   ones for which hydraulic geometry is now calculated; if we want to
**   use a reach slope everywhere, then we need to also calculate the
**   hydraulic geometry everywhere.
**   AND make sure that the node 'crn' is not at a lower elevation than
**   the outlet; in the latter case, it returns slope=0; if it runs into
**   an infinite loop either searching for the 10-widths-downstream node
**   or the outlet node, it returns a negative number (-1)
\************************************************************************/
#define kLargeNumber 1000000
double tLNode::getSlope()
{
   int ctr;
   double rlen, curlen, slp, delz, downz;
   tLNode *dn, *on;

   assert( flowedge != 0 );
   assert( flowedge->getLength()>0 ); // failure means lengths not init'd

   if( Meanders() )
   {
      rlen = 10.0 * chan.chanwidth;
      curlen = 0.0;
      dn = this;
      ctr = 0;
      //built a couple of 'firewalls' in these loops; first, quit loop if
      //we have a flowedge with zero length; second, if we've somehow gotten
      //an infinite loop, return a negative number as a flag to the calling
      //routine that it needs to update the network (tStreamNet::UpdateNet())
      //and reaches (tStreamMeander::MakeReaches). We've run across this loop
      //bug in a call from tStreamMeander::InterpChannel, which is called by
      //tStreamMeander::MakeReaches.
      while( curlen < rlen && dn->getBoundaryFlag() == kNonBoundary &&
             dn->flowedge->getLength() > 0 )
      {
         ctr++;
         if( ctr > kLargeNumber )
         {
            return (-1);
            //ReportFatalError("infinite loop in tLNode::GetSlope(), 1st loop");
         }
         curlen += dn->flowedge->getLength();
         dn = dn->getDownstrmNbr();
         assert( dn != 0 );
      }
      if(curlen <= 0){
         cout<<"going to die in getSlope(), curlen is "<<curlen<<endl;
         TellAll();
         assert(curlen>0.0);
      }
      
      assert( curlen > 0 );
      downz = dn->z;
      //if( timetrack >= kBugTime ) cout << "GetSlope 1; " << flush;
      delz = z - downz;
      //if( timetrack >= kBugTime ) cout << "GS 2; " << flush;
      slp = delz / curlen;
      //if( timetrack >= kBugTime ) cout << "GS 3; " << flush;
      on = dn;
      ctr = 0;
      while( on->getBoundaryFlag() == kNonBoundary &&
             on->flowedge->getLength() > 0 )
      {
         ctr++;
         if( ctr > kLargeNumber )
         {
            return (-1);
            //ReportFatalError("infinite loop in tLNode::GetSlope(), 2nd loop");
         }
         on = on->getDownstrmNbr();
      }
      if( z - on->z < 0.0 ) slp = 0.0;
   }
   else{
      slp = (z - getDownstrmNbr()->z ) / flowedge->getLength();
   }
   
   //if( timetrack >= kBugTime ) cout << "GS 4; " << endl << flush;
   if( slp>=0.0 ) return slp;
   else return 0.0;
}
#undef kLargeNumber

double tLNode::getDSlopeDt()
{
   double rlen, curlen, slp;
   tLNode *dn;
   assert( flowedge != 0 );
   assert( flowedge->getLength()>0 ); // failure means lengths not init'd
   if( Meanders() )
   {
      rlen = 10.0 * chan.chanwidth;
      curlen = 0.0;
      dn = this;
      while( curlen < rlen && dn->getBoundaryFlag() == kNonBoundary )
      {
         curlen += dn->flowedge->getLength();
         dn = dn->getDownstrmNbr();
      }
      assert( curlen > 0 );
      //slp = (dzdt - dn->dzdt + uplift - dn->uplift ) / curlen;
   }
   
   else
   {
      curlen = flowedge->getLength();
      assert( curlen > 0.0 );
      dn = getDownstrmNbr();
      //slp = (dzdt - dn->dzdt + uplift - dn->uplift ) / curlen;
   }
   slp = ( dzdt - dn->dzdt + uplift - dn->uplift ) / curlen;
   if( slp>=0.0 ) return slp;
   else return 0.0;
}




//void tLNode::insertFrontSpokeList( tEdge *eptr )             //tNode
//{spokeList.insertAtFront( eptr );}

//void tLNode::insertBackSpokeList( tEdge *eptr )              //tNode
//{spokeList.insertAtBack( eptr );}

/*****************************************************************************\
**
**
**      DistNew: replaces disnew; the distance of the node's (newx, newy)
**              from the line 
**              formed by points  p0->(newx, newy) and p1->(newx, newy)
**      Data members updated: 
**      Called by: 
**      Calls:  
**        
**
\*****************************************************************************/
double tLNode::DistNew(tLNode * p0,tLNode * p1 )
{
  double a,b,c,res, xp, yp, x0, y0, x1, y1;

  x0 = ( p0->Meanders() ) ? p0->chan.migration.newx : p0->x;
  y0 = ( p0->Meanders() ) ? p0->chan.migration.newy : p0->y;
  x1 = ( p1->Meanders() ) ? p1->chan.migration.newx : p1->x;
  y1 = ( p1->Meanders() ) ? p1->chan.migration.newy : p0->y;
  xp = ( Meanders() ) ? chan.migration.newx : x;
  yp = ( Meanders() ) ? chan.migration.newy : y;
  
  a = y1 - y0; //(p1->chan.migration.newy) - (p0->chan.migration.newy);  
  b = x0 - x1; //-((p1->chan.migration.newx) - (p0->chan.migration.newx));
  c = -( a * x0 + b * y0 );
    //-((a*(p0->chan.migration.newx)) + (b*(p0->chan.migration.newy)));
  res = (a * xp + b * yp + c) / sqrt( a * a + b * b );
    //res = (a*chan.migration.newx + b*chan.migration.newy + c) / sqrt(a*a+b*b);
  if( res < 0 ) res = -res;
  return res;
}



/**************************************************************************\
**
**  TellAll
**
**  ...is a debugging routine that prints out a bunch of info about a node.
**
\**************************************************************************/
#ifndef NDEBUG
void tLNode::TellAll()
{
   tLNode * nbr;
   int i;
   
   cout << " NODE " << id << ":\n";
   cout << "  x=" << x << " y=" << y << " z=" << z;
   if( edg ) {
      cout << "  points to edg #" << edg->getID() << endl;
      cout << "  dr area: " << getDrArea() << "  disch: " << getQ()
           << "  boundary: " << boundary << "  flood: " << flood
           << "\n  varea: " << varea << endl;
      
      if( flowedge ) {
         nbr = (tLNode *)flowedge->getDestinationPtrNC();
         cout << "  Flows along edg " << flowedge->getID() << " to node "
              << nbr->getID() << " at (" << nbr->getX() << ","
              << nbr->getY() << "," << nbr->getZ() << ")\n    with vedglen "
              << flowedge->getVEdgLen() << endl;
         flowedge->TellCoords();
         //cout<<"  ccwedge of flowedge is "<<flowedge->getCCWEdg()->getID();
         //cout<<" originates at "<<flowedge->getCCWEdg()->getOriginPtrNC()->getID()<<endl;
         cout << "  qs: " << qs << "  qsin: " << qsin << "  slp: "
              << flowedge->getSlope() << "  reg: " << reg.thickness << endl;
         for(i=0; i<numg; i++)
             cout<<"  qsi "<<i<<" "<<qsm[i];
         cout<<endl;
         //for(i=0; i<numg; i++)
         //cout<<"  qsini "<<i<<" "<<qsinm[i];
         //cout<<endl;         
         cout<<" creation time top "<<getLayerCtime(0);
         cout<<"numlayers is "<<getNumLayer()<<endl;
         int j;
         for( j=0; j<getNumLayer(); j++ )
             for(i=0; i<numg; i++)
                 cout<<"  dgrade "<<i<<" "<<getLayerDgrade(j,i);
         cout << "  dzdt: " << dzdt << "  drdt: " << drdt;
         cout<<" meanders "<< Meanders()<<endl;
         cout << "  width: " << getHydrWidth() << " tau: " << getTau();
	 cout << " taucrit: " << getTauCrit() << endl;
      }
      else cout << "  Flowedg is undefined\n";
      
   }
   else cout << "  edg is undefined!\n";
   
   cout << "layerlist addresses are:\n";
   layerlist.DebugTellPtrs();
}
#endif


/**************************************************************************\
**
**  EroDep
**
**  Erodes or Deposits material of depth dz.
**
**  Note: as of now, just changes elev. Should also change alluvial
**  thickness and set to zero if negative (bedrock erosion).
**  Alluvial thickness is outdated - just use other EroDep.
**
**  1/15/98 gt
\**************************************************************************/
void tLNode::EroDep( double dz )
{
   z += dz;
   //cout << "  eroding " << id << " by " << dz << endl << flush;
   
   //sed += dz;
   //if( sed<0 ) sed=0;
   reg.thickness += dz;
   if( reg.thickness < 0 ) reg.thickness = 0.0;
}

void tLNode::InsertLayerBack( tLayer lyr )
{
   layerlist.insertAtBack( lyr );
}


void tLNode::setAlluvThickness( double val )
{
   //reg.thickness = ( val >= 0.0 ) ? val : 0.0;
   // gt changed for performance speedup (if stmt shouldn't ever be needed)
   assert( val>=0 );
   reg.thickness = val;
}

double tLNode::getAlluvThickness() const {return reg.thickness;}

tArray< double >
tLNode::getAlluvThicknessm( ) const 
{
   return reg.dgrade;
}

/*Xvoid tLNode::setVegErody( double val )
{surf.vegerody = ( val >= 0.0 ) ? val : 0.0;}

double tLNode::getVegErody() const {return surf.vegerody;}
*/

void tLNode::setBedErody( double val )
{rock.erodibility = ( val >= 0.0 ) ? val : 0.0;}

double tLNode::getBedErody() const {return rock.erodibility;}

void tLNode::setReachMember( int val )
{chan.migration.reachmember = ( val == 0 || val == 1 ) ? val : 0;}

int tLNode::getReachMember() const {return chan.migration.reachmember;}

void tLNode::setQs( double val ) {qs = val;}

void tLNode::setQs( int i, double val ) 
{
   if(i>=numg)
       ReportFatalError( "Trying to index sediment sizes that don't exist ");
   qsm[i] = val;
   qs += val;
}

double tLNode::getQs() const {return qs;}

double tLNode::getQs( int i)
{
   if(i>=numg)
       ReportFatalError( "Trying to index sediment sizes that don't exist ");
   return qsm[i];
}

tArray< double >
tLNode::getQsm( ) const
{
   return qsm;
}

void tLNode::setQsin( double val ) {qsin = val;}

void tLNode::setQsin( int i, double val ) 
{
   if(i>=numg)
       ReportFatalError( "Trying to index sediment sizes that don't exist ");
   if(i>=qsinm.getSize()){
      cout<<"trying to set index "<<i<<" but size of array is "<<qsinm.getSize()<<endl<<flush;
      TellAll();
   }
   qsinm[i]=val;
   double tot=0;
   for(int j=0; j<numg; j++)
       tot+=qsinm[j];   
   qsin=tot;
}

void tLNode::addQsin( double val ) 
{
   qsin += val;
}

void tLNode::addQsin( int i, double val )
{
   if(i>=numg)
       ReportFatalError( "Trying to index sediment sizes that don't exist ");
   qsinm[i] += val;
   qsin += val;
   
}

void tLNode::addQsin( tArray< double > val )
{
   int i;
   for(i=0; i<val.getSize(); i++){
      qsin +=  val[i];
      qsinm[i] += val[i];
   }
   
}

void tLNode::addQs( double val ) 
{
   qs += val;
}

void tLNode::addQs( int i, double val )
{
   if(i>=numg)
       ReportFatalError( "Trying to index sediment sizes that don't exist ");
   qsm[i] += val;
   qs += val;
   
}

void tLNode::addQs( tArray< double > val )
{
   int i;
   
   for(i=0; i<val.getSize(); i++){
      qsm[i] += val[i];
      qs += val[i];
   }

   
}


double tLNode::getQsin() const {return qsin;}

double tLNode::getQsin( int i )
{
   if(i>=numg)
       ReportFatalError( "Trying to index sediment sizes that don't exist ");
   return qsinm[i];
}

tArray< double >
tLNode::getQsinm( ) const
{
   return qsinm;
}

void tLNode::setGrade( int i, double size )
{
   if(i>=numg)
       ReportFatalError("Trying to set a grain size for an index which is too large");
   grade[i] = size;
}

double tLNode::getMaxregdep() const
{
   return maxregdep;
}

int tLNode::getNumg() const 
{
   return numg;
}

void tLNode::setNumg( int size )
{
   numg = size;
}

double tLNode::getGrade( int i )
{
   return grade[i];
}

tArray< double >
tLNode::getGrade( ) const
{
   return grade;
}  

double tLNode::getLayerCtime( int i ) const
{
   tLayer hlp;
   hlp = layerlist.getIthData(i);
   return hlp.getCtime();
}

void tLNode::setLayerCtime( int i, double tt)
 {
    tListIter<tLayer> ly ( layerlist );
    tLayer  * hlp;
    hlp=ly.FirstP();

    int n=0;

    while(n<i){
       n++;
       hlp=ly.NextP();
    }

    hlp->setCtime( tt );
}

double tLNode::getLayerRtime( int i ) const
{
   tLayer hlp;
   hlp = layerlist.getIthData(i);
   return hlp.getRtime();
}

void tLNode::setLayerRtime( int i, double tt)
 {
    tListIter<tLayer> ly ( layerlist );
    tLayer  * hlp;
    hlp=ly.FirstP();
    
    int n=0;

    while(n<i){
       n++;
       hlp=ly.NextP();
    }

    hlp->setRtime( tt );
}

double tLNode::getLayerEtime( int i ) const
{
   tLayer hlp;
   hlp = layerlist.getIthData(i);
   return hlp.getEtime();
}

void tLNode::setLayerEtime( int i, double tt)
 {
    tListIter<tLayer> ly ( layerlist );
    tLayer  * hlp;
    hlp=ly.FirstP();
    
    int n=0;

    while(n<i){
       n++;
       hlp=ly.NextP();
    }

    hlp->setEtime( tt );
}

void tLNode::addLayerEtime( int i, double tt)
 {
    tListIter<tLayer> ly ( layerlist );
    tLayer  * hlp;
    hlp=ly.FirstP();
    
    int n=0;

    while(n<i){
       n++;
       hlp=ly.NextP();
    }

    hlp->addEtime( tt );
 }

double tLNode::getLayerDepth( int i ) const
{
   tLayer hlp;
   if( layerlist.isEmpty() )
     {
       cout << "** WARNING lyr list empty\n";
   cout << " NODE " << id << ":\n";
   cout << "  x=" << x << " y=" << y << " z=" << z;
       
     }
   hlp = layerlist.getIthData(i);
   return hlp.getDepth();
}

void tLNode::setLayerDepth( int i, double dep)
 {
    tListIter<tLayer> ly ( layerlist );
    tLayer  * hlp;
    hlp=ly.FirstP();
    
    int n=0;

    while(n<i){
       n++;
       hlp=ly.NextP();
    }

    hlp->setDepth( dep );
}

double tLNode::getLayerErody( int i ) const
{
   tLayer hlp;
   hlp = layerlist.getIthData(i);
   return hlp.getErody();
}

void tLNode::setLayerErody( int i, double ero)
 {
    tListIter<tLayer> ly ( layerlist );
    tLayer  * hlp;
    hlp=ly.FirstP();
    
    int n=0;

    while(n<i){
       n++;
       hlp=ly.NextP();
    }

    hlp->setErody( ero );
}


int tLNode::getLayerSed( int i ) const
{
   tLayer hlp;
   hlp = layerlist.getIthData(i);
   return hlp.getSed();
}

void tLNode::setLayerSed( int i, int s)
 {
    tListIter<tLayer> ly ( layerlist );
    tLayer  * hlp;
    hlp=ly.FirstP();
    
    int n=0;

    while(n<i){
       n++;
       hlp=ly.NextP();
    }

    hlp->setSed( s );
}

double tLNode::getLayerDgrade( int i, int num ) const
{
   //cout << "newgetlayerdgrade\n";
   //tLayer * lyr;
   //lyr = layerlist.getIthDataPtr(i);
   //return lyr->getDgrade(num);
   tLayer hlp;
   hlp = layerlist.getIthData(i);
   return hlp.getDgrade(num);
}

void tLNode::setLayerDgrade( int i, int g, double val)
 {
    tListIter<tLayer> ly ( layerlist );
    tLayer  * hlp;
    hlp=ly.FirstP();
    
    int n=0;

    while(n<i){
       n++;
       hlp=ly.NextP();
    }

    hlp->setDgrade(g, val );
}

/************************************************************************\
  tLNode::LayerInterpolation

  Called by tGrid when new tLNode is added to the grid, or a node is
  moved.
  tri = The old triangle in which the new node will be contained.
        The attributes of its vertices are used for the interpolation.
  tx = The new x location at which the point will be located.
  ty = The new y location at which the point will be located.
  time = current time

  You may wonder why x and y need to be passed.  This is just in case a
  node is moving, in which case you need to reference the newx and newy
  value.  If the node is just being added, only need to reference x and y.
  To make things generic, I just decided to pass in the x and y values.

  This function first creates the new list of layers and then sets the
  layerlist of the node equal to the list which was just created.
  This is done in case a node is moving and its old layers are needed
  for the interpolation.

  LayerInterpolation now checks to make sure that you are interpolating
  only between non-boundary nodes.  
  4/1999 Makes sure that top layer depth = maxregdep

  Designed by GT, Coded by NG, last update 4/1999 
\*************************************************************************/

#define kAncient -999
void tLNode::LayerInterpolation( tTriangle * tri, double tx, double ty, double time )
{
   assert(tri!=0);
   
   //cout<<endl<<"tLNode::LayerInterpolation....";
//   cout<<" current x = "<<x<<" current y = "<<y;
//   cout<<" newx= "<<tx<<" newy= "<<ty<<endl<<flush;
   
   tLNode *lnds[3];
   int i;

   int numnodes=0;
   for(i=0; i<=2; i++){
      if(!tri->pPtr(i)->getBoundaryFlag()){
         lnds[numnodes] = (tLNode *) tri->pPtr(i);
         numnodes++;
      }
   }
   //cout<<"numnodes = "<<numnodes<<" newx= "<<tx<<" newy= "<<ty<<endl;

   tList< tLayer > helplist; //Make the layer list first.  When
   //the list is made, then set the nodes layerlist equal to helplist.

   if( numnodes==3 ){
      tArray<int> layindex(3);//What layer are you at for each node
      tArray<double> dep(3);//The depth of the layer with the matching age
      tArray<double> age(3);//The age of the current layer at each node
      tArray<int> sed(3);//The sediment flag of the current layer at each node
      
      double dist[3]; //distance b/w new point and the three triangle points
      dist[0]=DistanceBW2Points(tx, ty, lnds[0]->getX(), lnds[0]->getY());
      dist[1]=DistanceBW2Points(tx, ty, lnds[1]->getX(), lnds[1]->getY());
      dist[2]=DistanceBW2Points(tx, ty, lnds[2]->getX(), lnds[2]->getY());

      double CA=-1;
      
      for(i=0; i<=2; i++){
         if(lnds[i]->getLayerRtime(0)>CA){
            CA=lnds[i]->getLayerRtime(0);
         }
         sed[i]=(lnds[i]->getLayerSed(0)<1 ? 0 : 1);
         if(sed[i]>0)
             age[i]=lnds[i]->getLayerRtime(0);
         else
             age[i]=kAncient;
         layindex[i]=0;//Initialize layindex
      }
      
      //cout<<"Current age is "<<CA<<endl;
      //cout<<"Sed Type is "<<SD<<endl;
      //CA now contains the youngest surface layer time of the three nodes.
      //Remember that LayerRtime is the most recent time visited, which
      //implies larger time = younger layer
      double newtex; //texture value of new layer
      double sum;  //helper used to calculate the new texture & erody value
      double newerody; //erodability of new layer
      double newdep; //interpolated depth value
      double newetime; //interpolated exposure time
      tLayer layhelp; //set values in this layer then insert in back of layerlist
      layhelp.setDgradesize(numg);
      layhelp.setCtime(time);

      do
      {
         newtex=0;
         sum=0;
         newerody=0;
         newetime=0;
         
         for(i=0; i<=2; i++){
            if(age[i]<=CA+10&&age[i]>=CA-10&&sed[i]>0){
               //layer is in the window of acceptable ages
               //set node texture helper to texture of that layer
               //NOTE - This only works for two sizes right now.
               //set node depth helper to depth of that layer
               newtex+=lnds[i]->getLayerDgrade(layindex[i],0)/dist[i];
               newerody+=lnds[i]->getLayerErody(layindex[i])*lnds[i]->getLayerDepth(layindex[i])/dist[i];
               newetime+=lnds[i]->getLayerEtime(layindex[i])*lnds[i]->getLayerDepth(layindex[i])/dist[i];
               sum+=lnds[i]->getLayerDepth(layindex[i])/dist[i];
               dep[i]=lnds[i]->getLayerDepth(layindex[i]);
               layhelp.setSed(sed[i]);
               if( layindex[i] != lnds[i]->getNumLayer()-1 ){
                  //More layers below, continue to search
                  layindex[i]+=1;
                  sed[i]=(lnds[i]->getLayerSed(layindex[i])<1 ? 0 : 1);
                  if(sed[i]>0 && lnds[i]->getLayerRtime(layindex[i])>=0)
                      age[i]=lnds[i]->getLayerRtime(layindex[i]);
                  else
                      age[i]=kAncient;
               }
            }
            else{
               dep[i]=0;
               //cout<<"not in correct range, dep set to 0"<<endl;
            }
         }
         
         
         //Now find the values of the texture and layer depth in the new
         //location, given what we found out above.  Set everything, and
         //then insert the new layer.
         
         newdep=PlaneFit(tx, ty, lnds[0]->get2DCoords(), lnds[1]->get2DCoords(),
                         lnds[2]->get2DCoords(), dep);
         if(newdep>grade[0]){
            newtex=newtex/sum;   
            layhelp.setDepth(newdep);
            layhelp.setDgrade(0,newtex*newdep);
            layhelp.setErody(newerody/sum);
            layhelp.setEtime(newetime/sum);
            layhelp.setRtime(CA);
            if(numg>1)
                layhelp.setDgrade(1,(1-newtex)*newdep);
            helplist.insertAtBack( layhelp );
         }
         CA=age[0];
         for(i=1; i<=2; i++){
            if(age[i]>CA){
               CA=age[i];
            }
         }
         //cout<<endl<<"after iter, current age is set to "<<CA<<endl;
         
      }while(CA>kAncient);
      
      //If there is a deep sediment layer, interpolate it
      if(lnds[0]->getLayerRtime(layindex[0])<0){//here assuming if one
         //node has a deep sed layer, they all do
         //TODO fix this!
         newtex=0;
         newerody=0;
         sum=0;
         for(i=0; i<=2; i++){
            //NOTE - This only works for two sizes right now.
            newtex+=lnds[i]->getLayerDgrade(layindex[i],0)/dist[i];
            newerody+=lnds[i]->getLayerErody(layindex[i])*lnds[i]->getLayerDepth(layindex[i])/dist[i];
            sum+=lnds[i]->getLayerDepth(layindex[i])/dist[i];
            dep[i]=lnds[i]->getLayerDepth(layindex[i]);
            layindex[i]++;
         }
         newtex=newtex/sum;
         newdep=PlaneFit(tx, ty, lnds[0]->get2DCoords(),lnds[1]->get2DCoords(),
                         lnds[2]->get2DCoords(), dep);
         layhelp.setDepth(newdep);
         layhelp.setDgrade(0,newtex*newdep);
         if(numg>1)
             layhelp.setDgrade(1,(1-newtex)*newdep);
         layhelp.setSed(1); 
         layhelp.setErody(newerody/sum);
         layhelp.setRtime(-1);
         layhelp.setEtime(0);
         helplist.insertAtBack( layhelp );
      }

      newtex=0;
      newerody=0;
      newetime=0;
      sum=0;
      //Now, you need to interpolate bedrock layer here
      for(i=0; i<=2; i++){
         //NOTE - This only works for two sizes right now.
         //debugging routine
//           if(lnds[i]->getNumLayer()<=layindex[i]){
//              lnds[i]->TellAll();
//              for(int j=0; j<lnds[i]->getNumLayer(); j++){
//                 cout << "layer " << j+1<<endl ;
//                 cout << lnds[i]->getLayerCtime(j);
//                 cout << " " << lnds[i]->getLayerRtime(j);
//                 cout << " "<<lnds[i]->getLayerEtime(j)<<endl;
//                 cout << lnds[i]->getLayerDepth(j);
//                 cout << " " << lnds[i]->getLayerErody(j);
//                 cout << " " << lnds[i]->getLayerSed(j) << endl;
//                 cout << lnds[i]->getLayerDgrade(j,0);
//              }
//           } 
         newtex+=lnds[i]->getLayerDgrade(layindex[i],0)/dist[i];
         newerody+=lnds[i]->getLayerErody(layindex[i])*lnds[i]->getLayerDepth(layindex[i])/dist[i];
         newetime+=lnds[i]->getLayerEtime(layindex[i])*lnds[i]->getLayerDepth(layindex[i])/dist[i];
         sum+=lnds[i]->getLayerDepth(layindex[i])/dist[i];
         dep[i]=lnds[i]->getLayerDepth(layindex[i]);
      }
      newtex=newtex/sum;
      newdep=PlaneFit(tx, ty, lnds[0]->get2DCoords(),lnds[1]->get2DCoords(),
                      lnds[2]->get2DCoords(), dep);
      layhelp.setDepth(newdep);
      layhelp.setDgrade(0,newtex*newdep);
      if(numg>1)
          layhelp.setDgrade(1,(1-newtex)*newdep);
      layhelp.setSed(0);
      layhelp.setErody(newerody/sum);
      layhelp.setRtime(0);
      layhelp.setEtime(newetime/sum);

      helplist.insertAtBack( layhelp );

      //This is the new part that nmg added where the top layers are
      //truncated according to the new elevation and the elevation
      //around the node
      double theoryz;
      tArray<double> elevs(3);
      elevs[0]=lnds[0]->getZ();
      elevs[1]=lnds[1]->getZ();
      elevs[2]=lnds[2]->getZ();
      
      theoryz=PlaneFit(tx, ty, lnds[0]->get2DCoords(), lnds[1]->get2DCoords(),
                         lnds[2]->get2DCoords(), elevs );
      double diff=theoryz-getZ();

      /*cout<<"new z "<<getZ()<<" theoretical z "<<theoryz<<" diff "<<diff<<endl;
      cout<<"before removal, layer info"<<endl;
      cout<<"numlayers = "<<helplist.getSize()<<endl;
      cout<<"depth of first layer = "<<helplist.FirstP()->getDepth()<<endl;*/
      
      while(diff>0){
         if(helplist.FirstP()->getDepth()<diff){
            diff-=helplist.FirstP()->getDepth();
            helplist.removeFromFront(layhelp);
         }
         else{
            double dpth = helplist.FirstP()->getDepth();
            tLayer * firstlay = helplist.FirstP();
            firstlay->setDepth(dpth-diff);
            diff=0;
         }
      }
      
   }
   else if( numnodes ==2 ){
      tArray<int> layindex(2);//What layer are you at for each node
      tArray<double> age(2);//The age of the current layer at each node
      tArray<int> sed(2);//The sediment flag of the current layer at each node
      
      double dist[2]; //distance b/w new point and the three triangle points
      dist[0]=DistanceBW2Points(tx, ty, lnds[0]->getX(), lnds[0]->getY());
      dist[1]=DistanceBW2Points(tx, ty, lnds[1]->getX(), lnds[1]->getY());
      double distsum=(1/dist[0])+(1/dist[1]);
      
      double CA=-1;
      
      for(i=0; i<=1; i++){
         if(lnds[i]->getLayerRtime(0)>CA){
            CA=lnds[i]->getLayerRtime(0);
         }
         sed[i]=(lnds[i]->getLayerSed(0)<1 ? 0 : 1);
         if(sed[i]>0)
             age[i]=lnds[i]->getLayerRtime(0);
         else
             age[i]=kAncient;
         layindex[i]=0;//Initialize layindex
      }
      
      //CA now contains the youngest surface layer time of the three nodes.
      //Remember that LayerRtime is the most recent time visited, which
      //implies larger time = younger layer
      double newtex; //texture value of new layer
      double sum;  //helper used to calculate the new texture & erody value
      double newerody; //erodability of new layer
      double newetime; //exposure time of new layer
      double newdep; //interpolated depth value
      tLayer layhelp; //set values in this layer then insert in back of layerlist
      layhelp.setDgradesize(numg);
      layhelp.setCtime(time);

      do
      {
         newtex=0;
         sum=0;
         newerody=0;
         newetime=0;
         newdep=0;
         
         for(i=0; i<=1; i++){
            if(age[i]<=CA+10&&age[i]>=CA-10&&sed[i]>0){
               //layer is in the window of acceptable ages
               //set node texture helper to texture of that layer
               //NOTE - This only works for two sizes right now.
               newtex+=lnds[i]->getLayerDgrade(layindex[i],0)/dist[i];
               newerody+=lnds[i]->getLayerErody(layindex[i])*lnds[i]->getLayerDepth(layindex[i])/dist[i];
               newetime+=lnds[i]->getLayerEtime(layindex[i])*lnds[i]->getLayerDepth(layindex[i])/dist[i];
               sum+=lnds[i]->getLayerDepth(layindex[i])/dist[i];
               newdep+=lnds[i]->getLayerDepth(layindex[i])/dist[i]/distsum;
               
               layhelp.setSed(sed[i]);
               if( layindex[i] != lnds[i]->getNumLayer()-1 ){
                  //More layers below, continue to search
                  layindex[i]+=1;
                  sed[i]=(lnds[i]->getLayerSed(layindex[i])<1 ? 0 : 1);
                  if(sed[i]>0 && lnds[i]->getLayerRtime(layindex[i])>=0)
                      age[i]=lnds[i]->getLayerRtime(layindex[i]);
                  else
                      age[i]=kAncient;
               }
            }
         }         
         
         //Now find the values of the texture and layer depth in the new
         //location, given what we found out above.  Set everything, and
         //then insert the new layer.
            
         
         if(newdep>grade[0]){
            newtex=newtex/sum;   
            layhelp.setDepth(newdep);
            layhelp.setDgrade(0,newtex*newdep);
            layhelp.setErody(newerody/sum);
            layhelp.setEtime(newetime/sum);
            layhelp.setRtime(CA);
            if(numg>1)
                layhelp.setDgrade(1,(1-newtex)*newdep);
            helplist.insertAtBack( layhelp );
         }
         CA=age[0];
         if(age[1]>CA){
            CA=age[1];
         }
         
      }while(CA>kAncient);

      //Now, you need to interpolate deep sediment layer (if it exists)
      if(lnds[0]->getLayerRtime(layindex[0])<0){//Just assume that if one
         //node has the deep sediment layer that they all do.
         //TODO fix this!
         newtex=0;
         newerody=0;
         sum=0;
         newdep=0;
         for(i=0; i<=1; i++){
            //NOTE - This only works for two sizes right now.
            newtex+=lnds[i]->getLayerDgrade(layindex[i],0)/dist[i];
            newerody+=lnds[i]->getLayerErody(layindex[i])*lnds[i]->getLayerDepth(layindex[i])/dist[i];
            sum+=lnds[i]->getLayerDepth(layindex[i])/dist[i];
            newdep+=lnds[i]->getLayerDepth(layindex[i])/dist[i]/distsum;
            layindex[i]++;
         }
         newtex=newtex/sum;
         layhelp.setDepth(newdep);
         layhelp.setDgrade(0,newtex*newdep);
         if(numg>1)
             layhelp.setDgrade(1,(1-newtex)*newdep);
         layhelp.setSed(1); 
         layhelp.setErody(newerody/sum);
         layhelp.setRtime(-1);
         layhelp.setEtime(0);
         helplist.insertAtBack( layhelp );
      }

      //Finally, interpolate bedrock
      newtex=0;
      newerody=0;
      newetime=0;
      sum=0;
      newdep=0;
      for(i=0; i<=1; i++){
         //NOTE - This only works for two sizes right now.
         newtex+=lnds[i]->getLayerDgrade(layindex[i],0)/dist[i];
         newerody+=lnds[i]->getLayerErody(layindex[i])*lnds[i]->getLayerDepth(layindex[i])/dist[i];
         newetime+=lnds[i]->getLayerEtime(layindex[i])*lnds[i]->getLayerDepth(layindex[i])/dist[i];
         sum+=lnds[i]->getLayerDepth(layindex[i])/dist[i];
         newdep+=lnds[i]->getLayerDepth(layindex[i])/dist[i]/distsum;
      }
      newtex=newtex/sum;
      layhelp.setDepth(newdep);
      layhelp.setDgrade(0,newtex*newdep);
      if(numg>1)
          layhelp.setDgrade(1,(1-newtex)*newdep);
      layhelp.setSed(1); 
      layhelp.setErody(newerody/sum);
      layhelp.setRtime(0);
      layhelp.setEtime(newetime/sum);
      helplist.insertAtBack( layhelp );
   

      //This is the new part that nmg added where the top layers are
      //truncated according to the new elevation and the elevation
      //around the node
      double theoryz;
      theoryz=LineFit(0, lnds[0]->getZ(), dist[0]+dist[1],
                      lnds[1]->getZ(), dist[0]);
      double diff=theoryz-getZ();
      while(diff>0){
         if(helplist.FirstP()->getDepth()<diff){
            diff-=helplist.FirstP()->getDepth();
            helplist.removeFromFront(layhelp);
         }
         else{
            double dpth = helplist.FirstP()->getDepth();
            tLayer * firstlay = helplist.FirstP();
            firstlay->setDepth(dpth-diff);
            diff=0;
         }
      }
      
   }
   else{
      //only one case left - only one non-boundary node 
      //check done w/assertion at beginning
      //just set layerlist = to layerlist of non-boundary node
      tLayer layhelp; //set values in this layer then insert in back of layerlist
      layhelp.setDgradesize(numg);
      layhelp.setCtime(time);

      for(i=0; i<lnds[0]->getNumLayer(); i++){
         layhelp.setDepth(lnds[0]->getLayerDepth(i));
         layhelp.setDgrade(0,lnds[0]->getLayerDgrade(i,0));
         if(numg>1)
          layhelp.setDgrade(1,lnds[0]->getLayerDgrade(i,1));
         layhelp.setSed(lnds[0]->getLayerSed(i));
         layhelp.setErody(lnds[0]->getLayerErody(i));
         layhelp.setEtime(lnds[0]->getLayerEtime(i));
         layhelp.setRtime(lnds[0]->getLayerRtime(i));
         helplist.insertAtBack( layhelp );
      }

      //This is the new part that nmg added where the top layers are
      //truncated according to the new elevation and the elevation
      //around the node
      double diff=lnds[0]->getZ()-getZ();
      while(diff>0){
         if(helplist.FirstP()->getDepth()<diff){
            diff-=helplist.FirstP()->getDepth();
            helplist.removeFromFront(layhelp);
         }
         else{
            double dpth = helplist.FirstP()->getDepth();
            tLayer * firstlay = helplist.FirstP();
            firstlay->setDepth(dpth-diff);
            diff=0;
         }
      }  
   }

   if(maxregdep-helplist.FirstP()->getDepth()>0.001 && helplist.getIthData(1).getSed()>0){
      //top layer is too small, want to maintain maxregdep for erosion reasons
      //Because NG doesn't really know how to manipulate lists,
      //this is probably done in a cockeyed way.  First remove
      //top layer from list, update it, then add it back to front of list.
      tLayer firstlay;
      helplist.removeFromFront(firstlay);
      double diff=maxregdep-firstlay.getDepth();      
      double newerody=firstlay.getErody()*firstlay.getDepth();
      double newetime=firstlay.getEtime()*firstlay.getDepth();
      double newrtime=firstlay.getRtime()*firstlay.getDepth();
      do{
         tLayer * nextlay=helplist.FirstP();
         if(nextlay->getDepth()<diff){
            //add entire contents of layer below to top layer
            for(i=0; i<numg; i++)
                firstlay.addDgrade(i,nextlay->getDgrade(i));
            newerody+=nextlay->getErody()*nextlay->getDepth();
            newetime+=nextlay->getEtime()*nextlay->getDepth();
            newrtime+=nextlay->getRtime()*nextlay->getDepth();
            diff-=nextlay->getDepth();
            helplist.removeFromFront(*nextlay);
         }
         else{
            //don't remove entire layer below, only take some material
            double prevdep = nextlay->getDepth();
            for(i=0; i<numg; i++){
               double moving = diff*nextlay->getDgrade(i)/prevdep;
               firstlay.addDgrade(i,moving);
               nextlay->addDgrade(i,-1*moving);
            }
            newerody+=nextlay->getErody()*diff;
            newetime+=nextlay->getEtime()*diff;
            //don't want to include negative exposure time flag
            if(nextlay->getRtime()>0)
                newrtime+=nextlay->getRtime()*diff;
            diff=0;
         }
      }while(diff>0.001);
      firstlay.setErody(newerody/firstlay.getDepth());
      firstlay.setEtime(newetime/firstlay.getDepth());
      firstlay.setRtime(newrtime/firstlay.getDepth());
      helplist.insertAtFront(firstlay);
   }
   else if(helplist.FirstP()->getRtime()<0){
      //top layer is the deep sediment layer.  This messes up things 
      //for layering in general.  Form a new top layer with material from
      //the bottom layer.
      tLayer layhelp;
      layhelp.setSed(1);
      layhelp.setEtime(0);
      layhelp.setCtime(time);
      layhelp.setRtime(0);
      layhelp.setDgrade(0,maxregdep*helplist.FirstP()->getDgrade(0)/helplist.FirstP()->getDepth());
      if(numg>1)
          layhelp.setDgrade(1,maxregdep*helplist.FirstP()->getDgrade(1)/helplist.FirstP()->getDepth());
      layhelp.setErody(helplist.FirstP()->getErody());
      helplist.FirstP()->setDepth( helplist.FirstP()->getDepth()-maxregdep);
      helplist.insertAtFront( layhelp );
      
   }
   
   layerlist=helplist;
//     if(getNumLayer()<=2)
//         cout<<"layerinterp made 2 layers at "<<x<<", "<<y<<endl;

   //Below  if for debugging purposes
   /*if(getLayerEtime(0)<0 || (tx>505.0 && tx<506.0 && ty>331.0 && ty<332.0)  ){
    */
   /*i=0; 
   cout<<"x "<<tx<<" y "<<ty; 
   while(i<layerlist.getSize()){ 
      cout << "layer " << i+1 <<endl; 
      cout << getLayerCtime(i); 
      cout << " " << getLayerRtime(i); 
      cout << " " << getLayerEtime(i)<<endl; 
      cout << getLayerDepth(i); 
      cout << " " << getLayerErody(i); 
      cout << " " << getLayerSed(i) << endl; 
      cout << getLayerDgrade(i,0) ; 
      if( numg>1 ) cout << " " << getLayerDgrade(i,1);
      i++;
      cout<<endl;
   }
   
   int j;
   for(j=0; j<numnodes; j++){
      cout<<endl;
      cout<<"node "<<j<<" x "<<lnds[j]->getX()<<" y "<<lnds[j]->getY();
      cout<<" boundary "<<lnds[j]->getBoundaryFlag()<<endl;
      tLNode * nicn = lnds[j];
      i=0;
      while(i<nicn->getNumLayer()){
         cout << "layer " << i+1<<endl ;
         cout << nicn->getLayerCtime(i);
         cout << " " << nicn->getLayerRtime(i);
         cout << " "<<nicn->getLayerEtime(i)<<endl;
         cout << nicn->getLayerDepth(i);
         cout << " " << nicn->getLayerErody(i);
         cout << " " << nicn->getLayerSed(i) << endl;
         cout << nicn->getLayerDgrade(i,0);
         if( numg>1 ) cout << " " << nicn->getLayerDgrade(i,1) << endl;
         i++;
      }
      }*/ 
   
}
#undef kAncient


/*************************************************************\
  tLNode::WarnSpokeLeaving(tEdge * edglvingptr)

  This function is called when an edge is being removed from the edge list.
  If flowedge is pointing to the edge
  which will be removed, this flowedge must be updated.

  edglvingptr is, as it sais, a pointer to the edge which will be
  removed.

  Called from tGrid::ExtricateEdge

  9/98 NG and GT
\*******************************************************************/
void tLNode::WarnSpokeLeaving(tEdge * edglvingptr)
{
   //cout<<"tLNode::WarnSpokeLeaving..... node #"<<id<<endl;

   //Make sure that edg pointer in tNode won't be affected
   tNode::WarnSpokeLeaving( edglvingptr );

   if( edglvingptr == flowedge ){
      do{
         flowedge = flowedge->getCCWEdg();
      }while( (flowedge->getBoundaryFlag()==kClosedBoundary) && (flowedge != edglvingptr) );
      assert( flowedge->getOriginPtr() == this );
      
      //After looping around edges, if flow is along a non-flow edge,
      //make this a closedboundary node.  This is potentially a bad
      //thing to do since a node will be marked as a boundary node 
      //but not put at the end of the node list.  For now we ignore this
      //TODO fix this - probably in tGrid::ExtricateEdge
//        if(flowedge->getBoundaryFlag()==kClosedBoundary){
//           boundary = kClosedBoundary;
//           cout<<"node "<<getID()<<" x "<<getX()<<" y "<<getY()<<" set to boundary in WarnSpokeLeaving"<<endl<<flush;         
//        }
   }
}

/**********************************************************************\
 **
 ** tLNode::InitializeNode()
 **
 ** A virtual function.  Used when nodes are creating during the evolution
 ** process, not at the start of a run.  
 ** - Assigns one of the spokes as the flow edge.
 ** - sets size of Qsi and Qsini in case not properly set.
 ** 
 **
 ** Called from tGrid::AddNode
 **
 ** Created 1/1999 NG
 ***********************************************************************/
void tLNode::InitializeNode()
{
   int debugcount=0;

   // If we're not a boundary node, make sure we have a valid flowedge
   if( boundary==kNonBoundary )
   {
      flowedge=edg;
      do {
         flowedge = flowedge->getCCWEdg();
         debugcount++;
         assert( debugcount<10000 );
      } 
      while( (flowedge->getBoundaryFlag()==kClosedBoundary) && flowedge!=edg );
   }
   
   // Set size of sediment influx and outflux arrays to number of grain sizes
   if( qsm.getSize()!=numg ) {
      qsm.setSize( numg );
      qsinm.setSize( numg );
   }

   // cout<<"tLNode::InitializeNode node "<<id<<" flow edge "<<flowedge->getID()<<endl;
}


/***********************************************************************\
  tArray<double> tLNode::EroDep

  This function erodes and deposits material and updates
  the layering at the same time.  When eroding, layers are updated
  to maintain the maximum layer depth as long as there is sediment
  to refill them.  Bedrock layers can never receive material.
  
  Takes : - i is layer to erode/dep from/into
          - valgrd contains the depth of each grain size to erode
          (negative value) or deposit (positive value)
          - tt is the current time - for use if depositing

  Returns : an array containing the actual depth of each grain
            size that was actually ero'd/dep'd.  Only an issue
            if eroding because you might be limited in what you
            can erode (stuff might not be there, only erode from
            one layer at a time).

  Calls : - addtoLayer (both versions)-  for updating the layers below,
          and (if eroding) finding out what to add to the layer you
          are eroding from.
          - removeLayer - sometimes wierd things happens and removing
          everything doesn't really remove everything.  (Why????)
          Don't want really small layers, so remove them.
          - makeNewLayerBelow - have a maximum layer depth.  If you
          are depositing into a layer and it will surpass maxlayerdepth,
          need to move stuff into layer below.  If there is no room
          below, just make a new layer.

  Created 6/98 ng

\***********************************************************************/


tArray<double> tLNode::EroDep( int i, tArray<double> valgrd, double tt)
{
   int g;
   double amt, val, olddep;
   tArray<double> update;
   update.setSize(numg);
   tArray<double> hupdate;
   hupdate.setSize(numg);
   
   //NIC these are for testing
   //Xbefore=getLayerDepth(i);
   tArray<double> helper=valgrd;
   // int numlay=getNumLayer();
   //Xint h;
   //int nh = 0;
   
   //if(getLayerDgrade(i,0)/getLayerDepth(i)>0.99)
   //{
      //if(x<560.418 && x>560.416){
   // nh=1;
   // cout<<"layer "<<i<<" size 0 "<< getLayerDgrade(i,0)<<" size 1 "<<getLayerDgrade(i,1)<<endl;
   // if(getNumLayer()==2){
   //    cout<<"ERODEP x "<<x<<" y "<<y<<" numlayers "<<getNumLayer();
   //    cout<<" to erode b4 "<<valgrd[0]<< " " <<valgrd[1]<<endl;
   // }
   //}
   
   g=0;
   val=0;
   double max, min;
   max = -10000.;
   min = 10000.;
   while(g<numg){
      if(-1*valgrd[g]>getLayerDgrade(i,g))
          valgrd[g]=-1*getLayerDgrade(i,g);
      // Checking to see that there is enough stuff
      if(valgrd[g]>max) max = valgrd[g];
      if(valgrd[g]<min) min = valgrd[g];
      val+=valgrd[g];
      g++;
   }

//   if(nh){
//     if(x<560.418 && x>560.416){
//      cout<<"amt to ero/dep "<<val<<endl;
//      cout<<"individual amts "<<valgrd[0]<< " " <<valgrd[1]<<endl;
//         for(g=0; g<getNumLayer(); g++){
//            cout<<getLayerCtime(g)<<" "<<getLayerRtime(g)<<" "<<getLayerEtime(g)<<" "<<getLayerSed(g)<<" "<<getLayerDepth(g)<<" "<<getLayerErody<<endl;
//         }
//   }
          

   // val is now set to the total amount of erosion or deposition
   z += val;
   // total elevation has now been changed

   double sume, sumd, olde; //used for calculation of exposure time
   
   // Branch according to whether there is:
   // (1) erosion of all size-classes
   // (2) deposition of all size-classes
   // (3) erosion of some, deposition of others
   if(max <= 0 && min < -0.0000000001)
   {
      // CASE OF EROSION IN ALL SIZE CLASSES
      if(getLayerSed(i) != 0 && getLayerSed(i) == getLayerSed(i+1) && getLayerDepth(i)+val<=maxregdep)
      {
         // Updating will also be done if entering this statement
         // No updating of Bedrock layers.
         // Only update if layer below has the same material
         sume=(getLayerDepth(i)+val)*getLayerEtime(i); //for averaging
         //of the exposure time.
         olde=getLayerEtime(i+1);
         while(val<-0.000000001 && getLayerSed(i) == getLayerSed(i+1))
         {
            // keep eroding until you either get all the material you
            // need to refill the top layer, or you run out of material
            hupdate = addtoLayer(i+1, val);//remove stuff from lower layer
            g=0;
            sumd=0;
            while(g<numg)
            {
               sumd-=hupdate[g];
               val-=hupdate[g];//hupdate stores texture of material that will
               //refil the top layer
               update[g] += hupdate[g];//need update cause might need
               //to get material from more than one layer
               g++;
            }
            sume+=sumd*olde;
         }
         g=0;
         while(g<numg){
            addtoLayer(i, g, valgrd[g], -1); // Erosion 
            addtoLayer(i,g,-1*update[g],-1);//Updating with material from below
            g++;
         }
         setLayerEtime(i, sume/getLayerDepth(i));
      }
      else
      {
         // No updating, just eroding, don't need to change exposure time
         //Do this if you have only bedrock below, or if layer you are eroding
         //from is >maxregdepth-val (could be if lots of deposition)
         g=0;
         while(g<numg){
            addtoLayer(i, g, valgrd[g], -1); // Erosion done on this line
            g++;
         }
      }
      if(getLayerDepth(i)<1e-7)
          removeLayer(i);
   }
   else if(min >= 0.0 && max > 0.0000000001)
   {
      // CASE OF DEPOSITION IN ALL SIZE CLASSES
      // method seems to make good sense for surface layers
      // but may not be as appropriate for lower layers.
      // You may want to either make a provision for lower layers
      // or else write a new algorithm for them (will you ever deposit
      // into lower layers?  I don't know)
      // Also, no test done to make sure that you are depositing the right
      // material into the layer, would need to pass the flag for this.
      // For now assume that the test will be done in another place.
      if(getLayerSed(i)>0 && val < maxregdep){
         // top layer is sediment, so no issues
         // depositing less than an entire layer of stuff
         olde=getLayerEtime(i);
         if(getLayerDepth(i)+val>maxregdep){
            // Need to move stuff out of top layer to make room for deposited material
            if(getLayerSed(i) == getLayerSed(i+1) && getLayerDepth(i+1)+val<maxregdep)
            {
               // The layer below is of the appropriate material and has space
               amt = getLayerDepth(i)+val-maxregdep;//how much to move out
               setLayerEtime(i+1, (1/(getLayerDepth(i+1)+amt))*
                             (olde*amt + getLayerDepth(i+1)*getLayerEtime(i+1)));//lower layer's etime is now properly set.
               setLayerEtime(i, (1/maxregdep)*(olde*(getLayerDepth(i)-amt)));
               //now etime is set in the layer you are depositing into
               //note deposited material has an etime of 0
               olddep = getLayerDepth(i);
               g=0;
               while(g<numg){
                  addtoLayer(i+1,g,amt*getLayerDgrade(i,g)/olddep, -1);
                  // putting material from top layer to layer below
                  // nic, at this point you have decided not to change
                  // the recent time on the lower layer when you move
                  // stuff down into it.  Might think about this.
                  addtoLayer(i,g,-1*amt*getLayerDgrade(i,g)/olddep+valgrd[g], tt);
                  // changing top layer composition by depositing and removing
                  g++;
               }
            }
            else
            {
               // Need to create new layer
               amt = getLayerDepth(i)+val-maxregdep;
               setLayerEtime(i, (1/maxregdep)*(olde*(getLayerDepth(i)-amt)));
               //now etime is set in the layer you are depositing into
               //note deposited material has an etime of 0
               olddep = getLayerDepth(i);
               g=0;
               while(g<numg){
                  update[g]=amt*getLayerDgrade(i,g)/olddep;
                  // material which will be moved from top layer
                  addtoLayer(i,g,-1*amt*getLayerDgrade(i,g)/olddep+valgrd[g], tt);
                  // changing top layer composition by depositing
                  g++;
               }
               makeNewLayerBelow(i, getLayerSed(i), getLayerErody(i), update, tt);
               //When new layer is created then you change the recent time.
               // put material into a new layer which is made below
               setLayerEtime(i+1, olde); //now layer below has the proper etime
            }
         }
         
         else
         {
            // No need to move stuff out of top layer, just deposit
            setLayerEtime(i, (1/(getLayerDepth(i)+val))*olde*getLayerDepth(i));
            //set top layers etime
            g=0;
            while(g<numg){
               addtoLayer(i,g,valgrd[g],tt);
               g++;
            }
         }
      }
      else{
         // depositing sediment on top of br - create a new surface layer baby.
         // setting all new sediment layers to have erodibility of KRnew,
         // value read in at begining
         // Also use this if the amount of material deposited is 
         // greater than maxregdep
         makeNewLayerBelow(-1,1,KRnew,valgrd,tt);
      }
   }
   else if(max>0.0000000001 && min<-0.0000000001)
   {
      // CASE OFSIMULTANEOUS EROSION AND DEPOSITION:
      // Erosion of one or more sizes, deposition of other(s).
      // First, test whether the net effect is erosion or deposition.
      //Need if for the zero erosion case, in which nothing is done.
      if(val < 0)
      {
         // Case of net erosion
         if(getLayerSed(i) != 0)
         {
            //layer is sediment
            if(getLayerSed(i) == getLayerSed(i+1))
            {
               sume=(getLayerDepth(i)+val)*getLayerEtime(i); //for averaging
               //of the exposure time.
               olde=getLayerEtime(i+1);
               //There is material below to update with
               while(val<-0.000000001 && getLayerSed(i) == getLayerSed(i+1))
               {
                  // keep getting material from below  until you
                  // either get all the material you
                  // need to refill the top layer, or you run out of material
                  hupdate = addtoLayer(i+1, val);//remove stuff from
                  //lower layer, hupdate stores texture of material that will
                  //refil the top layer
                  g=0;
                  sumd=0;
                  while(g<numg)
                  {
                     sumd-=hupdate[g];
                     val-=hupdate[g];
                     update[g] += hupdate[g];
                     g++;
                  }
                  sume+=sumd*olde;
               }
               g=0;
               while(g<numg){
                  addtoLayer(i, g, valgrd[g], tt); // Erosion and deposition
                  addtoLayer(i,g,-1*update[g], tt);//Updating with material from below
                  //Set layer recent time because some deposition was done
                  g++;
               }
               setLayerEtime(i, sume/getLayerDepth(i));
            }
            else
            {
               //No updating will be done, but can put the
               //deposited material into the surface layer
               g=0;
               sumd=0;
               while(g<numg){
                  if(valgrd[g]>0)
                      sumd+=valgrd[g];
                  addtoLayer(i, g, valgrd[g], tt); // Erosion/Deposition
                  //Set layer recent time because some deposition was done
                  g++;
               }
               setLayerEtime(i, getLayerEtime(i)*(getLayerDepth(i)-sumd)/
                             getLayerDepth(i));
            }
         }
         else
         {
            //Layer is bedrock
            //First remove material from bedrock, then create a new layer
            //for the deposited material.
            for(g=0; g<numg; g++){
               update[g]=valgrd[g];//update stores the composition of new layer
               if(valgrd[g]<0){
                  addtoLayer(i, g, valgrd[g], -1);
                  update[g]=0.0;
               }
            }
            makeNewLayerBelow(i-1, 1, KRnew, update, tt);
            //New layer made with deposited material
         }
         if(getLayerDepth(i)<1e-7)
             removeLayer(i);
      }
      else
      {
         // Case of net deposition
         // First, we loop over all size-fractions. If erosion in a given
         // fraction is occurring, remove the material from the top layer;
         // if deposition, add it to val, which at the end of the loop
         // will record the total depth to be deposited.
         val=0; //val will now contain total amt to deposit
         for(g=0; g<numg; g++) {
            if(valgrd[g]<0) {
               // Here, we erode material in size fraction g
               update[g]=0;//If new layer needs to be made on top of BR
               //update will be composition.
               //Not necessarily used unles on bedrock, but I threw it
               //here anyway since you one should rarely enter this
               //neck of the woods.
               addtoLayer(i, g, valgrd[g], -1);
            }
            else {
               // Size fraction g is going to be deposited, so add it to the
               // total amount deposited (val)
               update[g]=valgrd[g];
               val+=valgrd[g];
            }
         }
         if(getLayerSed(i)>0 && val<maxregdep) {
            // Case in which top layer is "sediment" and depth to be
            // deposited is less than nominal layer thickness.
            // top layer is sediment, so no issues
            if(getLayerDepth(i)+val>maxregdep){
               // Need to move stuff out of top layer to make room for deposited mat
               if(getLayerSed(i) == getLayerSed(i+1) && getLayerDepth(i+1)+val<maxregdep)
               {
                  // The layer below is of the appropriate material and has space
                  amt = getLayerDepth(i)+val-maxregdep;//how much to move out
                  olddep = getLayerDepth(i);
                  olde=getLayerEtime(i);
                  setLayerEtime(i, (1/maxregdep)*olde*(getLayerDepth(i)-amt));
                  //exposure time of layer which is getting ero'd/dep'd is set
                  setLayerEtime(i+1, (1/(amt+getLayerDepth(i+1)))*
                                (amt*olde+getLayerDepth(i+1)*getLayerEtime(i+1)));//now exposure time of lower layer is set
                  g=0;
                  while(g<numg){
                     addtoLayer(i+1,g,amt*getLayerDgrade(i,g)/olddep, -1);
                     // putting material from top layer to layer below
                     // nic, at this point you have decided not to change
                     // the recent time on the lower layer when you move
                     // stuff down into it.  Might think about this.
                     addtoLayer(i,g,-1*amt*getLayerDgrade(i,g)/olddep+valgrd[g], tt);
                     // changing top layer composition by depositing and removing
                     g++;
                  }
               }
               else
               {
                  // Need to create new layer
                  amt = getLayerDepth(i)+val-maxregdep;
                  olddep = getLayerDepth(i);
                  olde=getLayerEtime(i);
                  setLayerEtime(i, (1/maxregdep)*olde*(getLayerDepth(i)-amt));
                  g=0;
                  while(g<numg){
                     update[g]=amt*getLayerDgrade(i,g)/olddep;
                     // material which will be moved from top layer
                     addtoLayer(i,g,-1*amt*getLayerDgrade(i,g)/olddep+valgrd[g], tt);
                     // changing top layer composition by depositing
                     g++;
                  }
                  makeNewLayerBelow(i, getLayerSed(i), getLayerErody(i), update, tt);
                  setLayerEtime(i+1, olde);
                  //When new layer is created then you change the time.
                  // put material into a new layer which is made below
               }
            }
            else
            {
               // No need to move stuff out of top layer, just deposit
               setLayerEtime(i, (1/(getLayerDepth(i)+val))*getLayerDepth(i)*
                             getLayerEtime(i));
               //now exposure time of top layer is properly set
               g=0;
               while(g<numg){
                  if(valgrd[g]>0)
                      addtoLayer(i,g,valgrd[g],tt);
                  g++;
               }
            }
         }
         else
         {
            //Layer is bedrock, so make a new layer on top to deposit into
            //or, depositing more than maxregdep
            makeNewLayerBelow(i-1, 1, KRnew, update, tt);
         }
      }
   }

   //Check to make sure that top layer is not too deep
   if(getLayerDepth(i)>1.1*maxregdep && getLayerSed(i) !=0 ){
      //Make a top layer that is maxregdep deep so that further erosion
      //is not screwed up
      hupdate = addtoLayer(i, -1*maxregdep);
      for(g=0; g<numg; g++)
          hupdate[g]=-1* hupdate[g];
      makeNewLayerBelow(-1,1,getLayerErody(i),hupdate, tt);
      setLayerRtime(i,0);
   }
   
   //   if(getLayerDepth(0)>1.1*maxregdep){
   //   cout<<"TOO MUCH SEDIMENT IN TOP LAYER"<<endl;
   //   TellAll();
   //}

//   if(nh)
//   {
//      cout<<"layer "<<i<<" size 0 "<< getLayerDgrade(i,0)<<" size 1 "<<getLayerDgrade(i,1)<<endl;
//      cout<<"at end, valgrd[0] "<<valgrd[0]<<" [1] "<<valgrd[1]<<endl;
//   }
   
          
   
   return valgrd;
}



/**************************************************************
 ** tLNode::addtoLayer(int i, int g, double val, double tt)
 **
 ** This function is called from the EroDep which updates layering.
 ** Used when material is added/removed to/from a layer, grain size 
 ** by grain size. i.e. texture may change
 ** i = layer to add to
 ** g = grain size class
 ** val = amount of that size to deposit
 ** tt = time of deposition/erosion
 **
 ** created by NG
 **
 **    Modifications:
 **      - 3/30/00: if erosion causes a layer to have zero 
 **        thickness, the layer is removed. (GT)
 *****************************************************************/
void tLNode::addtoLayer(int i, int g, double val, double tt)
{
   // For adding material to a layer, grain size by grain size
    tListIter<tLayer> ly ( layerlist );
    tLayer  * hlp;
    hlp=ly.FirstP();
    
    int n=0;
    while(n<i){
       n++;
       hlp=ly.NextP();
    }

    if(tt>0)
       hlp->setRtime( tt );
    hlp->addDgrade(g,val);

    //Although check really should be here, I think there may
    //be an issue because this function is called from a loop
    //and if the layer is removed before the loop is through -> big trouble
    // GT re-added this line, 3/00
    if(hlp->getDepth() <= 0)
        removeLayer(i);
}

/******************************************************************
 **  tLNode::addtoLayer(int i, double val)
 **  This addtolayer is only used for removing material from a layer.  
 **  This erosion will not change the texture of a layer, only depth.
 **  As always, the depth of each grain size class is also updated.
 **  i = layer to deplete
 **  val = amount to deplete by
 **  Returns an array containing the texture of material removed.
 **
 **  created NG
 **  Since only for erosion, nic modified this so that the time
 **  is not passed, since time will not be reset for erosion.
 ******************************************************************/
tArray<double> tLNode::addtoLayer(int i, double val)
{
   if(val>0)
       ReportFatalError("Using wrong function for depositing into a layer");

   tListIter<tLayer> ly ( layerlist );
   tLayer  * hlp;
   hlp=ly.FirstP();
   tArray<double> ret;
   ret.setSize(numg);
   double amt;
   
   int n=0;
   while(n<i){
      n++;
      hlp=ly.NextP();
   }

   if(hlp->getDepth()+val>1e-7)
   {
      // have enough material in this layer to fufill all erosion needs
      amt=hlp->getDepth();
      n=0;
      while(n<numg)
      {
         ret[n]=hlp->getDgrade(n)*val/amt;
         hlp->addDgrade(n,hlp->getDgrade(n)*val/amt);
         n++;
      }
      return ret;
   }
   else
   {
      // need to remove entire layer
      n=0;
      while(n<numg)
      {
         ret[n]=-1*hlp->getDgrade(n);
         n++;  
      }
      removeLayer(i);
      return ret;
   }
   
}

/*******************************************************
 ** tLNode::removeLayer(int i)
 **
 ** Removes layer i.
 **
 ** called by EroDep which updates layers.
 **
 **    Modifications:
 **      - 3/30/00: apparently the existing code was
 **        actually removing layer i+1. Bug fixed (GT).
 *******************************************************/
void tLNode::removeLayer(int i)
{
   //tLayer  * hlp;
   //tLayer niclay;
    //hlp=ly.FirstP();

    // GT added code in bug fix attempt 3/00:
    tLayer lay;

    if( i==0 )
        layerlist.removeFromFront( lay );
    else if( (i+1)==layerlist.getSize() )
        layerlist.removeFromBack( lay );
    else 
    {
       tListIter<tLayer> ly ( layerlist );
       int n;
       for( n=1; n<i; n++ )
           ly.Next();
       n = layerlist.removeNext( lay, ly.NodePtr() );
       assert( n>0 );
    }
    
    // end of GT added code

/*    int n=0;

    while(n<i){
       n++;
       hlp=ly.NextP();
    }

    if(i+1==layerlist.getSize()){
	layerlist.removeFromBack(*hlp );
	}
    else{	
       n=layerlist.removeNext((*hlp), layerlist.getListNode(hlp) );
       }
    if(n==0){
//       n=0;
//         while(n<layerlist.getSize()){
//            cout << "layer " << n+1 << " node ID "<< getID()<< endl;
//            niclay = layerlist.getIthData(n);
//            cout << "layer creation time is " << getLayerCtime(n) << endl;
//            cout << "layer recent time is " << getLayerRtime(n) << endl;
//            cout << "layer depth is " << getLayerDepth(n) << endl;
//            cout << "layer erodibility is " << getLayerErody(n) << endl;
//            cout << "is layer sediment? " << getLayerSed(n) << endl;
//            cout << "dgrade 1 is " << getLayerDgrade(n,0) << endl;
//            cout << "dgrade 2 is " << getLayerDgrade(n,1) << endl;
//            n++;  
//         }
       
       ReportFatalError("couldn't remove next layer");
       }*/
    
}

/*****************************************************************
 ** tLNode::makeNewLayerBelow(int i, int sd, double erd, tArray<double> sz, double tt)
 ** Makes a new layer below layer i.
 ** if i<0 then make new top layer.
 ** sd = sediment flag
 ** erd = erodibility
 ** sz = array of depths of each grain size class
 ** tt = time
 ** Creation and recent time set to current time, exposure time updated
 ** in erodep.
 ********************************************************************/
void tLNode::makeNewLayerBelow(int i, int sd, double erd, tArray<double> sz, double tt)
{
   tLayer hlp, niclay;
   int n;

   //cout << "making a new layer below " << i << " - layers before are:" << endl;

   hlp.setCtime(tt);
   hlp.setRtime(tt);
   hlp.setEtime(0);
   n=0;
   hlp.setDgradesize(numg);
   while(n<numg){
      hlp.setDgrade(n, sz[n]);
      n++;
   }
   hlp.setErody(erd);
   hlp.setSed(sd);

   if(i>=0){
      
      tListIter<tLayer> ly ( layerlist );
      tLayer  * pt;
      pt=ly.FirstP();
      
      n=0;
      
      while(n<i){
         n++;
         pt=ly.NextP();
      }  
      
      layerlist.insertAtNext(hlp, layerlist.getListNode(pt));
   }
   else{
      layerlist.insertAtFront(hlp);
   }
}

int tLNode::getNumLayer() const
{
   return layerlist.getSize();
}

//Returns depth of all the layers - pretty useless for the current model
double tLNode::getTotalLayerDepth() const
{
   double sz = layerlist.getSize();
   int i=0;
   double elev = 0;
   while(i<sz)
   {
      elev += getLayerDepth(i);
      i++;
   }
   return elev;
}


void tLNode::setDzDt( double val ) {dzdt = val;}

double tLNode::getDzDt() {return dzdt;}

void tLNode::setDrDt( double val ) {drdt = val;}

void tLNode::addDrDt( double val ) 
{
   drdt += val;
}

double tLNode::getDrDt() {return drdt;}

void tLNode::setXYZD( tArray< double > arr )
{
   chan.migration.xyzd = ( arr.getSize() == 4 ) ? arr : tArray< double >(4);
     //cout << "setXYZD: " << chan.migration.xyzd[0] 
     //   << " " << chan.migration.xyzd[1] << " " << chan.migration.xyzd[2]
     //   << " " << chan.migration.xyzd[3] << endl;
}

tArray< double >
tLNode::getXYZD() const {return chan.migration.xyzd;}


double tLNode::DistFromOldXY() const
{
   double xo, yo;
   tArray< double > oldpos( chan.migration.xyzd );
   
   xo = oldpos[0];
   yo = oldpos[1];
   
   return sqrt( (x-xo) * (x-xo) + (y-yo) * (y-yo) );
}

// Tests whether bedrock is exposed at a node
int tLNode::OnBedrock()
{
   // For multi-size model, criterion might be active layer thickness less
   // than a nominal thickness; here, it's just an arbitrary alluvial depth
   return ( reg.thickness<0.1 );
}

void tLNode::setUplift( double val ) {uplift = val;}

double tLNode::getUplift() const {return uplift;}
#undef kBugTime


/******************************************************************
 **  tLNode::CopyLayerList
 **
 **  Assigns the layer list in fromNode to this node. Information
 **  about the pre-existing layer list is lost.
 **
 **  Created in order to assign stratigraphy from a footwall point
 **  to the surface node when the fault plane is exhumed.
 **
 **  Created:  12/99 GT
 **
 ********************************************************************/
void tLNode::CopyLayerList( tLNode * fromNode )
{
   assert( fromNode!=0 );
   layerlist = fromNode->layerlist;
}
