/**************************************************************************\
**
**  tEolian.cpp: functions for class tEolian
**
**  (Created 1/99 by GT)
**
**  $Id: tEolian.cpp,v 1.2 1999-02-05 14:58:25 gtucker Exp $
\**************************************************************************/

#include "tEolian.h"

/**************************************************************************\
**
**  Constructor:
**
**  Takes a tInputFile as an argument and reads from it two parameters:
**  the deposition rate, and the # of grain sizes. The depositDepth
**  array is initialized to have the same size as the # of grain sizes,
**  with all values initialized to zero (done automatically by tArray).
**  tEolian::DepositLoess uses the first element depositDepth array,
**  which is then passed to tLNode::EroDep.
**
\**************************************************************************/
tEolian::tEolian( tInputFile &infile )
{
   int numg=1;
   
   loessDepRate = infile.ReadItem( loessDepRate, "LOESS_DEP_RATE" );
   numg = infile.ReadItem( numg, "NUMGRNSIZE" );
   depositDepth.setSize( numg );
   
}


/**************************************************************************\
**
**  tEolian::DepositLoess
**
**  Implements spatially uniform loess deposition, using EroDep to
**  perform the deposition. Note that if multiple sizes are used, the
**  deposits consist solely of the 1st grain size fraction, which is
**  assumed to be the finest.
**
**    Parameters:
**      gp -- pointer to the mesh to access the nodes
**      delt -- duration of current time interval
**      ctime -- current time
**    Notes:
**      - could be done as a polynomial surface
**
\**************************************************************************/
void tEolian::DepositLoess( tGrid<tLNode> *gp, double delt, double ctime )
{
   tGridListIter<tLNode> ni( gp->getNodeList() ); // iterator for nodes
   tLNode *cn;
   
   depositDepth[0] = loessDepRate*delt;
   for( cn=ni.FirstP(); ni.IsActive(); cn=ni.NextP() )
       cn->EroDep( 0, depositDepth, ctime );
}


