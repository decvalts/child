//-*-c++-*- 

/**************************************************************************/
/**
**  @file tStormGrid.cpp
**  @brief Functions for class tStormGrid
**
**
**  Created by Declan Valters, October 2015.
**
**
**  $Id: tStormGrid.cpp ,v 1.31 2004-06-16 13:37:42 childcvs Exp $
*/
/**************************************************************************/

#include "../tInputFile/tInputFile.h"
#include "../tTimeSeries/tTimeSeries.h"
#include "../tStorm/tStorm.h"

#include <iosfwd>
#include <sstream>

/**************************************************************************\
 **								 tSTORMGRID
 **  @fn tStormGrid( tInputFile &infile, tMesh<tLNode *mp)
 **  @brief Main constructor for tStormGrid
 **
 **  @param infile Input file from which parameters are read
 **  @param mp Pointer to the mesh
 **
 **  Takes the tInputFile as an argument and reads from it the various
 **  necessary parameters and spans the grid. The necessary input values are:
 **
 **  xllc = x(m) lower left corner
 **  yllc = y(m) lower left corner
 **
 **  griddx   = distance(m) between the tStormGrid nodes
 **  grwidth  = width (m) of the model domain
 **  grlength = length(m) of the model domain
 **
 **
\**************************************************************************/
tStormGrid::tStormGrid( tInputFile const &infile, tMesh<tLNode> *mp_) 
  : mp(mp_),
    StormNodeMatrix(0),
    StormConnect(0),
    imax(0),
    jmax(0)
{
  // Keep this pointer to the list of nodes for node list access
  assert( mp!=0 );
  
  int i, j, k;    // x and y position indices
  double x, y;    // x and y coordinates
  int step_x, step_y;  
  
  //iterator for tLNodes
  tMesh<tLNode>::nodeListIter_t ni( mp->getNodeList() ) ;  
  
  i = 0;
  x = y = 0;
  
  // Read in values related to dimensions and resolution of the mesh
  // and desired output format of the rainfall grid sections
  xcorner  = infile.ReadItem(xcorner,"XCORNER");
  ycorner  = infile.ReadItem(ycorner,"YCORNER");
  griddx   = infile.ReadItem(griddx,"GRIDDX");
  grid_xwidth  = infile.ReadItem(grid_xwidth,"GR_XWIDTH"); // need to add these
  grid_ywidth = infile.ReadItem(grid_ywidth,"GR_YWIDTH");
  
  int endtime  = infile.ReadItem(endtime,"RUNTIME");
  nWrite   = infile.ReadItem(nWrite,"OPINTRVL");
  
  // Calculate the max dims of the rectangle that encompasses the mesh
  imax = int(grid_xwidth/griddx);
  jmax = int(grid_ywidth/griddx);
  
  std::cout<<"     "<<std::endl;
  std::cout <<"StormGrid: number of nodes in x-direction = " << imax <<'\n';
  std::cout <<"StormGrid: number of nodes in y-direction = " << jmax <<'\n';
  
  // Call the constructor for the matrix of StormNodes 
  StormNodeMatrix = new tMatrix<tStormNode>(imax, jamx);
  StormConnect = new tMatrix<tTriangle*>(imax, jmax);
  
  // Fill one StormNode with the initial storm properties
  const tStormNode a_StormNode(infile)

  //Construct the grid by assigning coordinates and the initial layerlist 
  // to the nodes present in the tStormNode Matrix.
  for (i=0; i<max; i++)
  {
    for(j=0; j<jmax; j++)
    {
      x = xcorner + i * griddx;
      y = ycorner + j * griddx;
      
      (*StormNodeMatrix)(i,j).setX(x);
      (*StormNodeMatrix)(i,j).setY(y);
      (*StormNodeMatrix)(i,j).setI(i);
      (*StormNodeMatrix)(i,j).setJ(j);
      (*StormNodeMatrix)(i,j).setZ(0.0);
    }
  }
  
  // Code to handle making cross sections would go here:
  
  // Code
  
  // Build StormConnect
  if (1) // For DEBUG
  {
    std::cout
      <<"   \n"
      <<"Building the StormConnect table of triangles in constructor...."
      <<std::endl;
  }
  updateConnect();
  if (1) 
  {
    std::cout
      <<"Building the StormConnect table of triangles in constructor finished"
      <<"\n    "<<std::endl;
  }
  
  // Initialise the heights of StormGrid by interpolating between the triangles
  // of the tMesh
  if (1)
  { 
    std::cout<<" Initializing tStormGrid elevations by interpolation, for the first Time "<<std::endl;
  }
  InterpolateElevations();

  setSectionBase(); 
  // DEBUG FUNCTION, all Stormnodes have to know their initial, Stormigraphy basis    
  if (1) 
  {
    std::cout
      <<" Finished Initializing tStormGrid elevations by interpolation, for the first Time "
      <<"\n    "<<std::endl;
  }
}

/// Construct a tStormGrid from another tStormGrid object
tStormGrid::tStormGrid( const tStormGrid& orig )
  : xcorner(orig.xcorner),
    ycorner(orig.ycorner),
    griddx(orig.griddx),
    grid_xwidth(orig.grid_xwidth),
    grid_ywidth(orig.grid_ywidth),
    mp(orig.mp),              	 // ptr to triangular mesh
    imax(orig.imax),
    jmax(orig.jmax),
    //optSurferFiles(orig.optSurferFiles),
    //nWrite(orig.nWrite),
    //section(orig.section),          // array with 10 section locations, 5x, and 5y
    //surface(orig.surface),	 // timeslice specific surface area
    //subsurface(orig.subsurface),    // timeslice specific subsurface cummulative height in the entire floodplain
    //subsurface_mbelt(orig.subsurface_mbelt),    // timeslice specific subsurface cummulative height in meander belt
    outputTime(orig.outputTime)    // array with all the exact timings of t    
{
  const int nrM = orig.StormNodeMatrix->getNumRows();
  const int ncM = orig.StormNodeMatrix->getNumCols();
  StormNodeMatrix = new tMatrix<tStormNode>( nrM, ncM );
  
  const int nrC = orig.StormNodeMatrix->getNumRows();
  const int ncC = orig.StormNodeMatrix->getNumCols();
  StormConnect = new tMatrix<tTriangle*>( nrC, ncC );
}

/// Operator overloader
tStormGrid& tStormGrid::operator=(const tStormGrid& right )
{
  if ( &right != this )
  {
    xcorner = right.xcorner;
    ycorner = right.ycorner;
    griddx = right.griddx;
    grid_xlength = right.grid_xlength;
    grid_ylength = right.grid_ylength;
    mp = right.mp;
    imax = right.imax;
    jmax = right.jmax;
    
    outputTime = right.outputTime;
    
    const int nrM = right.StormNodeMatrix->getNumRows();
    const int ncM = right.StormNodeMatirx->getNumCols();
    
    StormNodeMatrix = new tMatrix<tStormNode>( nrM, nrC );
    
    for( int i=0; i<nrM; ++i)
    {
      for( int j=0; j<ncM; ++j)
      {
        (*StormNodeMatrix)(i,j) = (*right.StormNodeMatrix)(i,j);
       }
    }    
    // Connectivity to to tMesh
    const int nrC = right.StormNodeMatrix->getNumRows();
    const int nrM = right.StormNodeMatrix->getNumCols();
    
    StormConnect = new tMatrix<tTriangle*>( nrC, ncC );
    
    for( int i=0; i<nrC; ++i)
    {
      for( int j=0; i<ncC; ++j )
      {
        (*StormConnect)(i, j) = (*right.StormConnect)(i, j);
      }
    }
  }
  return *this;
}
    
/**************************************************************************\
 **
 **  @fn tStormGrid
 **  @brief destructor
 **
\**************************************************************************/        
tStormGrid::~tStormGrid()
{
  mp = 0;
  delete StormNodeMatrix;
  delete StormConnect;
}

// Set section base went here (removed)

/**************************************************************************\
 **
 **  @class
 **  @brief Find a rectangular box in the Storm grid contained within
 **  a given triangle.
 **
\**************************************************************************/
class TriBox
{
  TriBox();
public:
  int imin;
  int imax;
  int jmin;
  int jmax;
  
  TriBox(tTriangle const *, double, double, double);
  
  bool constainsNone() const 
  {
    return (imin > imax) || (jmin > jmax);
  }
};

TriBox::TriBox(tTriangle const *ct,
                double xcorner,
                double ycorner,
                double griddx)
{
  tNode const * const node1Ptr = ct->pPtr(0);
  
  double maxx = node1Ptr->getX();
  double maxy = node1Ptr->getY();
  double minx = maxx;
  double miny = maxy;
  
#define COMPUTE_MINMAX(NODEID) \
    do { \
      tNode const * const nodePtr = ct->pPtr(NODEID); \
      const double xx = nodePtr->getX(); \
      const double yy = nodePtr->getY(); \
      maxx = max( maxx, xx ); \
      maxy = max( maxy, yy ); \
      minx = min( minx, xx ); \
      miny = min( miny, yy ); \
    } while(0)
    
  COMPUTE_MINMAX(1);
  COMPUTE_MINMAX(2);
  
#undef COMPUTE_MINMAX

  imin = int( ceil((minx-xcorner)/griddx));
  imax = int(floor((maxx-xcorner)/griddx));
  jmin = int( ceil((miny-ycorner)/griddx));
  jmax = int(floor((maxy-ycorner)/griddx));
}

/**************************************************************************\
 **
 **  tStormGrid::updateConnect
 **  @brief update connectivity table StormConnect
 **
\**************************************************************************/
void tStormGrid::updateConnect()
{
  // nullify table
  for (int i=0; i<imax; ++i)
  {
    for(int j=0; j<jmax; ++j)
    {
      (*StormConnect)(i,j) = NULL;
    }
  }
  
  // Loop over the TIN triangles by iterating over a linked list
  // 'triIter' containing pointers to triangles in the mesh
  tTriangle *ct;
  
  tMesh< tLNode >::triListIter_t triIter( mp->getTriList() );
  for( ct = triIter.FirstP(); !( triIter.AtEnd() ); ct = triIter.NextP() )
  {
    // Find the box containing the current triangle
    const TriBox thisBox( ct, xcorner, ycorner, griddx );
    // but what if triangle does not contain any StormNode
    if (thisBox.containsNone())
      continue;
    
    // Set StormConnect for the StormNode withinthe current triangle
    {
      // clip the bounds within the actual bounds of the StormGrid
      const int bimin = max(0, thisBox.imin);
      const int bimax = min(getImax() - 1, thisBox.imax);
      const int bjmin = max(0, thisBox.imax); 
      const int bjmax = min(getJmax() - 1, thisBox.jmax);
      
      // find which StormNode is contained within the current triangle
      for (int i=bimin; i<=bimax; ++i)
      {
        for(int j=bjmin; j<=bjmax; ++j)
        {
          if ((*StormConnect)(i,j) != NULL)
            continue;
            
          // Assign i,j a pointer to the Triangle  
          if (ct->containsPoint( (*StormNodeMatrix)(i,j).getX(),
                      (*StormNodeMatrix)(i,j).getY() )
              )
          {
            (*StormConnect)(i, j) = ct;
          }
        }
      }
    }
  }
}
                
                
//--------------------CLASS STORMNODE----------------------------------------    
    
/**************************************************************************\
 **
 **  fn tStormNode( tInputFile &infile, tMesh<tLNode *mp );
 **  @brief Main constructor for tStormNode
 **
\**************************************************************************/    

// initialise static members of the class
//int StormNode::numg = 0;
//tArray<double> tStormNode::grade = 1;
//double tStormNode::maxregdepth = 1.0;
//double tStormNode::

// What information should we store in a single rain node?
// rain on or off
// rate
// duration(?)
// can these values be pulled from tStorm?

double tStormNode::rainrate = 0.0;
bool tStormNode::rain_flag = true;


// Default tStormNode constructor
tStormNode::tStormNode() :
  layerlist(),
  closest_node(0),
  x(0.0), y(0.0), z(0.0), sectionZ(0.0), newZ(0.0),
  i(0), j(0)
{
  // all done!
}

tStormNode::tStormNode( int ) :
  layerlist(),
  ClosestNode(0),
  x(0.0), y(0.0), z(0.0), sectionZ(0.0), newZ(0.0),
  i(0), j(0)
{
  
}

// StormNode constructor for indices and coordinates
tStormNode::tStormNode( double x_, double y_, const StormNode &orig) :
  layerlist(),
  ClosestNode(0),
  x(x_), y(y_), z(0.0), sectionZ(0.0), newZ(0.0),
  i(0), j(0)
{
  layerlist = orig.layerlist;
}

// tStormNode constructor for layerlist initialisation by inputfile
tStormNode::tStormNode( tInputFile const &infile )
{
  // To do
}

// So, when do I use this one...?
// Constructs a tStormNode from another tStormNode object
tStormNode::tStormNode( const tStormNode &orig )
  :
  layerlist(),
  ClosestNode(orig.ClosestNode),
  x(orig.x),
  y(orig.y),
  z(orig.z),
  sectionZ(orig.sectionZ),
  newz(orig.newz),
  i(orig.i),
  j(orig.j)
{
  layerlist = orig.layerlist;
}

// operator
tStormNode &tStormNode::operator=( const tStormNode &right )
{
  if( &right != this )
    {
      layerlist = right.layerlist;
      ClosestNode = right.ClosestNode;
      x = right.x;
      y = right.y;
      z = right.z;
      sectionZ = right.sectionZ;
      newz = right.newz;
      i = right.i;
      j = right.j;
    }
  return *this;
}

// Destructor
tStormNode::~tStormNode()
{
  if (0) //DEBUG
    std::cout << "    ~STORMNODE()" << std::endl;
}

//"set" and "get" functions for the coordinates:
void tStormNode::setX( double val ) {x = val;}
void tStormNode::setY( double val ) {y = val;}
void tStormNode::setZ( double val ) {z = val;}
void tStormNode::setNewZ( double val ) { newz = val;}
void tStormNode::setSectionBase( double val ) { sectionZ = val; }

void tStormNode::setI( int val ) {i = val;}
void tStormNode::setJ( int val ) {j = val;}

