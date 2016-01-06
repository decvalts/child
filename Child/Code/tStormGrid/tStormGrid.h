//-*-c++-*- 

/**************************************************************************/
/**
**  @file tStormGrid.h
**  @brief Header for class tStormGrid
**
**
**  Created by Declan Valters, October 2015.
* 
*   Class tStormGrid creates a second mesh type, an equidistant rectangular grid,
*   used to store gridded rainfall data. When tStreamNet runs in 'spatial rain-
*   fall mode' (name tbc), it pulls in the unique values of rainfall at each 
*   node in the TIN. Thus, giving the user a method for using spatially variable
*   precipitation or climate data over the landscape model domain.
**
**
**  $Id: tStormGrid.h,v 1.31 2004-06-16 13:37:42 childcvs Exp $
*/
/**************************************************************************/


#ifndef TSTORMGRID_H
#define TSTORMGRID_H

#include "../tInputFile/tInputFile.h"
#include "../tMatrix/tMatrix.h"
#include "../tList/tList.h"
#include "../tLNode/tLNode.h"

#include <iosfwd>
#include <sstream>

/**************************************************************************/
/**
 **  @class tStormGrid
 **
 **  Class tStormGrid contains the values and functions for a equidistant
 **  rectangular grid, which carries the storm areal distribution.
 **
 */
/**************************************************************************/

class tStormGrid
{
public:
  /// A function that reads in a rainfall grid from an external file
  tStormGrid( tInputFile const &infile);
  
  /// A function that creates a rainfall grid (storms) from specified parameters
  /// @details The user would provide the storm areal characteristics in the
  ///   input file, e.g. average diameter, morphology, or specifiy x, y, coords?
  /// @param Takes the input file and a pointer to the mesh
  tStormGrid( tInputFile const &infile, <tMesh<tLnode> *mp);
  
  /// @brief Creates a StormGrid from another StormGrid
  tStormGrid( const tStormGrid& );
  
  /// @brief The Operator function
  tStormGrid& operator=(const tStormGrid&);
  
  /// @brief The Destructor function
  ~tStormGrid();
  
  /// @brief Take a pointer to the mesh of nodes and assign pointer
  void setMesh( tMesh<tLNode>* ptr ){ mp = ptr;}

private:
  
  int xcorner;  // x-corner (lower left) for the StormGrid
  int ycorner;  // y-corner (lower left) for the StormGrid
  double griddx; // Inter-node distance of the StormGrid
  double grid_xwidth; // x-width of the StormGrid (metres)
  double grid_ywidth; // y-width of the StormGrid (metres)
  
  tMesh<tLNode> *mp;    // pointer to the triangular mesh
  
  tMatrix<tStormNode> *StormNodeMatrix; // pointer to the matrix of StormNodes
  tMatrix<tTriangle*> *StormConnect;    // pointer to matrix of triangles
  
  int imax, jmax,       // maximum bounding box coords of the triangle
  
  
};


/**************************************************************************/
/**
 **  @class tStormNode
 **
 **  Class tStormNode represents the functions
 **  and data members associated with an individual rect. grid node.
 **
 */
/**************************************************************************/
class tStormNode
{
public:

private:
};
