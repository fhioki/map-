/***************************************************************************
 *   Copyright (C) 2014 mdm                                                *
 *   marco.masciola@gmail.com                                              *
 *                                                                         *
 *   MAP++ is free software; you can redistribute it and/or modify it      *
 *   under the terms of the GNU General Public License as published by     *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.           *
 ***************************************************************************/

/**
 * @file 
 */

#include "map.h"
#include "maperror.h"
#include "freedata.h"


MAP_EXTERNCALL void MAP_InitInput_Delete(InitializationData* init_data)
{
  MAPFREE(init_data); 
};


void MAP_OtherState_Delete(ModelData* model_data)
{
  MAPFREE(model_data);
};


/** @addtogroup FortranCall */
/*@{*/
MAP_EXTERNCALL int free_init_data (InitializationData* init, char* map_msg, MAP_ERROR_CODE* ierr) 
{
  int sizeOfString = 0;
  int i = 0;
  
  MAPFREE(init->summaryFileName);

  for(i=0 ; i<init->sizeOfFullNodeString ; i++) {
    MAPFREE(init->expandedNodeInputString[i]);
  };
  MAPFREE(init->expandedNodeInputString);

  sizeOfString = init->sizeOfFullElementString;  
  for(i=0 ; i<sizeOfString ; i++) {
    MAPFREE(init->expandedElementInputString[i]);
  };
  MAPFREE(init->expandedElementInputString);  
  
  for(i=0 ; i<init->librarySize ; i++) {
    MAPFREE(init->libraryInputString[i]);
  };
  MAPFREE(init->libraryInputString);
  
  for(i=0 ; i<init->nodeSize ; i++) {
    MAPFREE(init->nodeInputString[i]);
  };
  MAPFREE(init->nodeInputString);
  
  for(i=0 ; i<init->elementSize ; i++) {
    MAPFREE(init->elementInputString[i]);
  };
  MAPFREE(init->elementInputString);  
  
  for(i=0 ; i<init->solverOptionsSize ; i++) {
    MAPFREE(init->solverOptionsString[i]);
  };
  MAPFREE(init->solverOptionsString);
  return MAP_SAFE;
};
/*@}*/


MAP_ERROR_CODE free_element(list_t *restrict element) 
{
  Element* iterElem = NULL;
  list_iterator_start(element); /* starting an iteration "session" */
  while (list_iterator_hasnext(element)) { /* tell whether more values available */
    iterElem = (Element*)list_iterator_next(element);
    MAPFREE( iterElem->psi.name ); 
    MAPFREE( iterElem->psi.units );
    MAPFREE( iterElem->alpha.name );
    MAPFREE( iterElem->alpha.units );
    MAPFREE( iterElem->alphaAtAnchor.name );
    MAPFREE( iterElem->alphaAtAnchor.units );
    MAPFREE( iterElem->l.name ); 
    MAPFREE( iterElem->l.units );
    MAPFREE( iterElem->lb.name ); 
    MAPFREE( iterElem->lb.units );
    MAPFREE( iterElem->Lu.name ); 
    MAPFREE( iterElem->Lu.units );
    MAPFREE( iterElem->h.name ); 
    MAPFREE( iterElem->h.units );
    MAPFREE( iterElem->H.name ); 
    MAPFREE( iterElem->H.units );
    MAPFREE( iterElem->V.name ); 
    MAPFREE( iterElem->V.units );
    MAPFREE( iterElem->HAtAnchor.name ); 
    MAPFREE( iterElem->HAtAnchor.units );
    MAPFREE( iterElem->VAtAnchor.name ); 
    MAPFREE( iterElem->VAtAnchor.units );
    MAPFREE( iterElem->forceAtFairlead.fx.name ); 
    MAPFREE( iterElem->forceAtFairlead.fx.units );
    MAPFREE( iterElem->forceAtFairlead.fy.name );
    MAPFREE( iterElem->forceAtFairlead.fy.units );
    MAPFREE( iterElem->forceAtFairlead.fz.name ); 
    MAPFREE( iterElem->forceAtFairlead.fz.units );
    MAPFREE( iterElem->forceAtAnchor.fx.name );
    MAPFREE( iterElem->forceAtAnchor.fx.units );
    MAPFREE( iterElem->forceAtAnchor.fy.name );
    MAPFREE( iterElem->forceAtAnchor.fy.units );
    MAPFREE( iterElem->forceAtAnchor.fz.name ); 
    MAPFREE( iterElem->forceAtAnchor.fz.units );
    MAPFREE( iterElem->T.name );
    MAPFREE( iterElem->T.units );
    MAPFREE( iterElem->TAtAnchor.name );
    MAPFREE( iterElem->TAtAnchor.units );
  
    /* don't let any pointers dangle */
    iterElem->lineProperty = NULL;      
    iterElem->label = NULL;
    iterElem->lineTension = NULL;
    iterElem->anchor = NULL; 
    iterElem->fairlead = NULL;
  };
  list_iterator_stop(element); /* ending the iteration "session" */  
};


MAP_ERROR_CODE free_node(list_t *restrict node)
{
  Node* iterNode = NULL;
  list_iterator_start(node);            /* starting an iteration "session" */
  while (list_iterator_hasnext(node)) { /* tell whether more values available */ 
    iterNode = (Node*)list_iterator_next(node);
    MAPFREE(iterNode->MApplied.name); 
    MAPFREE(iterNode->MApplied.units);
    MAPFREE(iterNode->BApplied.name);
    MAPFREE(iterNode->BApplied.units);
    
    MAPFREE(iterNode->externalForce.fx.name);
    MAPFREE(iterNode->externalForce.fx.units);
    MAPFREE(iterNode->externalForce.fy.name);
    MAPFREE(iterNode->externalForce.fy.units);
    MAPFREE(iterNode->externalForce.fz.name);
    MAPFREE(iterNode->externalForce.fz.units);
    
    MAPFREE(iterNode->positionPtr.x.name); 
    MAPFREE(iterNode->positionPtr.x.units);
    MAPFREE(iterNode->positionPtr.y.name); 
    MAPFREE(iterNode->positionPtr.y.units);
    MAPFREE(iterNode->positionPtr.z.name); 
    MAPFREE(iterNode->positionPtr.z.units);
    
    MAPFREE(iterNode->sumForcePtr.fx.name); 
    MAPFREE(iterNode->sumForcePtr.fx.units);
    MAPFREE(iterNode->sumForcePtr.fy.name); 
    MAPFREE(iterNode->sumForcePtr.fy.units);
    MAPFREE(iterNode->sumForcePtr.fz.name); 
    MAPFREE(iterNode->sumForcePtr.fz.units);
  };
  list_iterator_stop(node);             /* ending the iteration "session" */  

};


MAP_ERROR_CODE free_vessel(Vessel* floater) 
{
    /* Now delete the vessel information */
  MAPFREE(floater->xi);
  MAPFREE(floater->yi);
  MAPFREE(floater->zi);

  MAPFREE(floater->displacement.x.name);
  MAPFREE(floater->displacement.x.units);
  MAPFREE(floater->displacement.y.name);
  MAPFREE(floater->displacement.y.units);
  MAPFREE(floater->displacement.z.name);
  MAPFREE(floater->displacement.z.units);
          
  MAPFREE(floater->refOrigin.x.name);
  MAPFREE(floater->refOrigin.x.units);
  MAPFREE(floater->refOrigin.y.name);
  MAPFREE(floater->refOrigin.y.units);
  MAPFREE(floater->refOrigin.z.name);
  MAPFREE(floater->refOrigin.z.units);
          
  MAPFREE(floater->lineSumForce.fx.name);
  MAPFREE(floater->lineSumForce.fx.units);
  MAPFREE(floater->lineSumForce.fy.name);
  MAPFREE(floater->lineSumForce.fy.units);
  MAPFREE(floater->lineSumForce.fz.name);
  MAPFREE(floater->lineSumForce.fz.units);
          
  MAPFREE(floater->orientation.phi.name);
  MAPFREE(floater->orientation.phi.units);
  MAPFREE(floater->orientation.the.name);
  MAPFREE(floater->orientation.the.units);
  MAPFREE(floater->orientation.psi.name);
  MAPFREE(floater->orientation.psi.units);
  return MAP_SAFE;
}
