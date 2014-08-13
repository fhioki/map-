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
MAP_EXTERNCALL int free_init_data (InitializationData* init_data, char* map_msg, MAP_ERROR_CODE* ierr) 
{
  int ret = 0;

  ret = bdestroy(init_data->summaryFileName);
  // MAPFREE(init->summaryFileName);
  ret = bstrListDestroy(init_data->expandedNodeInputString);
  ret = bstrListDestroy(init_data->expandedElementInputString);
  ret = bstrListDestroy(init_data->libraryInputString);  
  ret = bstrListDestroy(init_data->nodeInputString);
  ret = bstrListDestroy(init_data->elementInputString);
  ret = bstrListDestroy(init_data->solverOptionsString);
  return MAP_SAFE;
};
/*@}*/


MAP_ERROR_CODE free_outlist(ModelData* data, char* map_msg, MAP_ERROR_CODE* ierr)
{
  MAPFREE(data->yList);
  return MAP_SAFE;
};


// MAP_ERROR_CODE free_element(list_t *restrict element) 
// {
//   Element* iterElem = NULL;
//   list_iterator_start(element); /* starting an iteration "session" */
//   while (list_iterator_hasnext(element)) { /* tell whether more values available */
//     iterElem = (Element*)list_iterator_next(element);
//     MAPFREE( iterElem->psi.name ); 
//     MAPFREE( iterElem->psi.units );
//     MAPFREE( iterElem->alpha.name );
//     MAPFREE( iterElem->alpha.units );
//     MAPFREE( iterElem->alphaAtAnchor.name );
//     MAPFREE( iterElem->alphaAtAnchor.units );
//     MAPFREE( iterElem->l.name ); 
//     MAPFREE( iterElem->l.units );
//     MAPFREE( iterElem->lb.name ); 
//     MAPFREE( iterElem->lb.units );
//     MAPFREE( iterElem->Lu.name ); 
//     MAPFREE( iterElem->Lu.units );
//     MAPFREE( iterElem->h.name ); 
//     MAPFREE( iterElem->h.units );
//     MAPFREE( iterElem->H.name ); 
//     MAPFREE( iterElem->H.units );
//     MAPFREE( iterElem->V.name ); 
//     MAPFREE( iterElem->V.units );
//     MAPFREE( iterElem->HAtAnchor.name ); 
//     MAPFREE( iterElem->HAtAnchor.units );
//     MAPFREE( iterElem->VAtAnchor.name ); 
//     MAPFREE( iterElem->VAtAnchor.units );
//     MAPFREE( iterElem->forceAtFairlead.fx.name ); 
//     MAPFREE( iterElem->forceAtFairlead.fx.units );
//     MAPFREE( iterElem->forceAtFairlead.fy.name );
//     MAPFREE( iterElem->forceAtFairlead.fy.units );
//     MAPFREE( iterElem->forceAtFairlead.fz.name ); 
//     MAPFREE( iterElem->forceAtFairlead.fz.units );
//     MAPFREE( iterElem->forceAtAnchor.fx.name );
//     MAPFREE( iterElem->forceAtAnchor.fx.units );
//     MAPFREE( iterElem->forceAtAnchor.fy.name );
//     MAPFREE( iterElem->forceAtAnchor.fy.units );
//     MAPFREE( iterElem->forceAtAnchor.fz.name ); 
//     MAPFREE( iterElem->forceAtAnchor.fz.units );
//     MAPFREE( iterElem->T.name );
//     MAPFREE( iterElem->T.units );
//     MAPFREE( iterElem->TAtAnchor.name );
//     MAPFREE( iterElem->TAtAnchor.units );
//   
//     /* don't let any pointers dangle */
//     iterElem->lineProperty = NULL;      
//     iterElem->label = NULL;
//     iterElem->lineTension = NULL;
//     iterElem->anchor = NULL; 
//     iterElem->fairlead = NULL;
//   };
//   list_iterator_stop(element); /* ending the iteration "session" */  
// };


MAP_ERROR_CODE free_node(list_t *restrict node)
{
  Node* iterNode = NULL;
  MAP_ERROR_CODE success = MAP_SAFE;
  list_iterator_start(node);            /* starting an iteration "session" */
  while (list_iterator_hasnext(node)) { /* tell whether more values available */ 
    iterNode = (Node*)list_iterator_next(node);

    success = bdestroy(iterNode->MApplied.name); 
    success = bdestroy(iterNode->MApplied.units);
    success = bdestroy(iterNode->BApplied.name);
    success = bdestroy(iterNode->BApplied.units);

    success = bdestroy(iterNode->externalForce.fx.name);
    success = bdestroy(iterNode->externalForce.fx.units);
    success = bdestroy(iterNode->externalForce.fy.name);
    success = bdestroy(iterNode->externalForce.fy.units);
    success = bdestroy(iterNode->externalForce.fz.name);
    success = bdestroy(iterNode->externalForce.fz.units);

    success = bdestroy(iterNode->positionPtr.x.name); 
    success = bdestroy(iterNode->positionPtr.x.units);
    success = bdestroy(iterNode->positionPtr.y.name); 
    success = bdestroy(iterNode->positionPtr.y.units);
    success = bdestroy(iterNode->positionPtr.z.name); 
    success = bdestroy(iterNode->positionPtr.z.units);

    success = bdestroy(iterNode->sumForcePtr.fx.name); 
    success = bdestroy(iterNode->sumForcePtr.fx.units);
    success = bdestroy(iterNode->sumForcePtr.fy.name); 
    success = bdestroy(iterNode->sumForcePtr.fy.units);
    success = bdestroy(iterNode->sumForcePtr.fz.name); 
    success = bdestroy(iterNode->sumForcePtr.fz.units);
  };
  list_iterator_stop(node);             /* ending the iteration "session" */  
  return MAP_SAFE;
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


MAP_ERROR_CODE map_free_types(MAP_InputType_t* u_type, MAP_ParameterType_t* p_type, MAP_ContinuousStateType_t* x_type, MAP_ConstraintStateType_t* z_type, MAP_OtherStateType_t* other_type, MAP_OutputType_t* y_type)
{
  /* inputs */
  MAPFREE(u_type->x);
  MAPFREE(u_type->y);
  MAPFREE(u_type->z);

  /* parameters are skipped for now; they are set in fortran since depth, gravity and sea density are set by glue code */

  /* continuous state */

  /* constraint state */  
  MAPFREE(z_type->H);     
  MAPFREE(z_type->V);     
  MAPFREE(z_type->x);     
  MAPFREE(z_type->y);     
  MAPFREE(z_type->z);     

  /* other state */
  MAPFREE(other_type->H); 
  MAPFREE(other_type->V); 
  MAPFREE(other_type->Ha);
  MAPFREE(other_type->Va);
  MAPFREE(other_type->x); 
  MAPFREE(other_type->y); 
  MAPFREE(other_type->z); 
  MAPFREE(other_type->xa);
  MAPFREE(other_type->ya);
  MAPFREE(other_type->za);
  MAPFREE(other_type->Fx_connect); 
  MAPFREE(other_type->Fy_connect); 
  MAPFREE(other_type->Fz_connect); 
  MAPFREE(other_type->Fx_anchor); 
  MAPFREE(other_type->Fy_anchor); 
  MAPFREE(other_type->Fz_anchor); 

  /* outputs */
  MAPFREE(y_type->Fx);    
  MAPFREE(y_type->Fy);    
  MAPFREE(y_type->Fz);    
  MAPFREE(y_type->wrtOutput);
  
  return MAP_SAFE;
};
