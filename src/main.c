/*
 * Copyright (c) 2014 mdm <marco.masciola@gmail.com>
 *
 * This file is part of MAP++.
 *
 * MAP++ is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * MAP++ is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with MAP++. If not, see <http://www.gnu.org/licenses/>.
 */

#include "map.h"
#include "maperror.h"
#include "protos.h"
#include "pyprotos.h"

extern const char MAP_ERROR_STRING[][1024];


/**
 * Frees data internal state data allcoated in the mapcall_msqs_init( ) function
 *
 * @acceses: none
 * @calledby: mapcall_msqs_init( )
 * @see: repeat_nodes( ), repeat_elements( )
 */
MAP_EXTERNCALL int free_init_data ( InitializationData* init, char *map_msg, MAP_ERROR_CODE *ierr ) 
{
  int sizeOfString = 0;
  int i = 0;
  
  MAPFREE(init->summaryFileName);

  for( i=0 ; i<init->sizeOfFullNodeString ; i++ ) {
    MAPFREE( init->fullNodeInputString[i] );
  };
  MAPFREE( init->fullNodeInputString );

  sizeOfString = init->sizeOfFullElementString;  
  for( i=0 ; i<sizeOfString ; i++ ) {
    MAPFREE( init->fullElementInputString[i] );
  };
  MAPFREE( init->fullElementInputString );  
  
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


/**
 *
 */
MAP_EXTERNCALL void mapcall_msqs_init ( MAP_InitInputType_t* initType, 
                                        MAP_InputType_t* uType,
                                        MAP_ParameterType_t* pType,
                                        MAP_ContinuousStateType_t* xType,
                                        void* none,
                                        MAP_ConstraintStateType_t* zType,
                                        MAP_OtherStateType_t* otherType,
                                        MAP_OutputType_t* yType,
                                        MAP_InitOutputType_t* ioType,
                                        MAP_ERROR_CODE *ierr,
                                        char *map_msg ) {    
  InitializationData* init = initType->object;   
  ModelData* data = otherType->object;
  MAP_ERROR_CODE success = MAP_SAFE;

  data->MAP_SOLVE_TYPE = -999;
  data->yList=NULL;  
  map_reset_universal_error(map_msg, ierr);
  do { 
    success = initialize_solver_data(&data->outerSolveData, &data->solverData, map_msg, ierr);    
    success = map_get_version(ioType);
    success = set_vessel_to_null(&data->vessel, map_msg, ierr);    
    /* 
       I've decided to call this later on in this function to initialize the line connection
       position relative to the body frame (this is simply uType, which isn't initialize at 
       this point).     
       // @rm initialize_vessel(&data->vessel);     
    */

    /* 
       initialize types; set doubles to -999.9, int=0, pointers=NULL 
       @todo: add other variables as neccessary. This needs to be fixed each time the registry
              is run an new varaibles are introdued 
    */  
    success = initialize_fortran_types(uType,pType,xType,zType,otherType,yType,ioType); 
    success = initialize_init_input(init); CHECKERRQ(MAP_FATAL_48);
    
    /* create a cable library link list */
    list_init(&data->cableLibrary);  
    list_attributes_copy(&data->cableLibrary, cable_library_meter, 1);
      
    /* create a node link list */
    list_init(&data->node);  
    list_attributes_copy(&data->node, node_meter, 1);
    
    /* create an element link list */
    list_init(&data->element);  
    list_attributes_copy(&data->element, cable_element_meter, 1);    
    
    data->sizeOfCableLibrary = init->librarySize;
    data->sizeOfNodes = init->nodeSize;
    data->sizeOfElements = init->elementSize;

    success = allocate_internal_states(data,map_msg,ierr); CHECKERRQ(MAP_FATAL_47);
    list_init(&data->yList->outList);  
    list_init(&data->yList->outListPtr);  

    success = set_model_options_list(data, init, map_msg, ierr); CHECKERRQ(MAP_FATAL_33);    
    success = set_cable_library_list(data, init, map_msg, ierr); CHECKERRQ(MAP_FATAL_16);

    success = repeat_nodes(data, init, map_msg, ierr);
    success = repeat_elements(data, init, map_msg, ierr);

    success = set_node_list(pType, uType, zType, otherType, yType, data, init->fullNodeInputString, map_msg, ierr); CHECKERRQ(MAP_FATAL_16);    
    success = set_element_list(zType, data, init->fullElementInputString, map_msg, ierr); CHECKERRQ(MAP_FATAL_16);    

    /* now create an output list to print to and output file. */
    list_attributes_copy(&data->yList->outList, vartype_meter, 1);  
    list_attributes_copy(&data->yList->outListPtr, vartype_ptr_meter, 1);  
    success = set_output_list(data, ioType, map_msg, ierr); 
    success = initialize_vessel(&data->vessel, uType); CHECKERRQ(MAP_FATAL_69);

    /*
      @todo: possible include hook for LM model here 
             eventually place MAP_SOVE = LM;
    */
    if (zType->x_Len!=0) {
      success = allocate_outer_solve_data(&data->outerSolveData, zType->x_Len, map_msg, ierr); CHECKERRQ(MAP_FATAL_72);
      data->MAP_SOLVE_TYPE = PARTITIONED;      
    } else {
      data->MAP_SOLVE_TYPE = MONOLITHIC;
    };

    printf("MAP environment properties (set externally)...\n");
    printf("    Gravity constant          [m/s^2]  : %1.2f\n", pType->g ); 
    printf("    Sea density               [kg/m^3] : %1.2f\n", pType->rhoSea );
    printf("    Water depth               [m]      : %1.2f\n", pType->depth );
    printf("    Vessel reference position [m]      : %1.2f , %1.2f , %1.2f\n", data->vessel.refOrigin.x.value, data->vessel.refOrigin.y.value, data->vessel.refOrigin.z.value); 

    success = initialize_cable_library_variables(data, pType, map_msg, ierr); CHECKERRQ(MAP_FATAL_41);
    success = set_line_variables_pre_solve(data, map_msg, ierr);
    success = reset_node_force_to_zero(data, map_msg, ierr);
    success = set_element_initial_guess(data, map_msg, ierr);
    success = first_solve(data, uType, zType, otherType, yType, map_msg, ierr); CHECKERRQ(MAP_FATAL_39);
    success = set_line_variables_post_solve(data, map_msg, ierr);    
    success = write_summary_file(init, pType, data, map_msg, ierr); CHECKERRQ(MAP_FATAL_37);           
    success = get_iteration_output_stream(yType, otherType, map_msg, ierr); // @todo CHECKERRQ()
    free_init_data(init, map_msg, ierr); /* @todo: does this need to be removed and called independently elsewhere after mapcall_msqs_init? */    
  } while (0);  
};


/**
 *
 */
MAP_EXTERNCALL void mapcall_msqs_update_states(double t,
                                               int interval,
                                               MAP_InputType_t* uType,
                                               MAP_ParameterType_t* pType,
                                               MAP_ContinuousStateType_t* xType,
                                               void* none,
                                               MAP_ConstraintStateType_t* zType,
                                               MAP_OtherStateType_t* otherType,
                                               MAP_ERROR_CODE* ierr,
                                               char* map_msg ) {
  ModelData* data = otherType->object;
  MAP_ERROR_CODE success = MAP_SAFE;

  map_reset_universal_error(map_msg, ierr);
  do {
    if (data->MAP_SOLVE_TYPE==MONOLITHIC) {
      success = line_solve_sequence(data, 0.0, map_msg, ierr);
    } else {
      success = node_solve_sequence(data, uType, zType, otherType, map_msg, ierr); // @todo CHECKERRQ()
    };    
  } while (0);
};    


/**
 *       
 */
MAP_EXTERNCALL void mapcall_msqs_calc_output(double t,
                                             MAP_InputType_t* uType,
                                             MAP_ParameterType_t* pType,
                                             MAP_ContinuousStateType_t* xType,
                                             void* none,
                                             MAP_ConstraintStateType_t* zType,
                                             MAP_OtherStateType_t* otherType,
                                             MAP_OutputType_t* yType,
                                             MAP_ERROR_CODE* ierr,
                                             char* map_msg ) {
  ModelData* data = otherType->object;
  MAP_ERROR_CODE success = MAP_SAFE;
  map_reset_universal_error(map_msg, ierr);
  success = get_iteration_output_stream(yType, otherType, map_msg, ierr); // @todo: CHECKERRQ();
};


/**
 *
 */
MAP_EXTERNCALL void mapcall_msqs_end(MAP_InputType_t* uType,
                                     MAP_ParameterType_t* pType,
                                     MAP_ContinuousStateType_t* xType,
                                     void* none,
                                     MAP_ConstraintStateType_t* zType,
                                     MAP_OtherStateType_t* otherType,
                                     MAP_OutputType_t* yType,                                                                           
                                     MAP_ERROR_CODE* ierr,
                                     char* map_msg ) {  
  Element* iterElem = NULL;
  Node* iterNode = NULL;
  CableLibrary* iterCableLibrary = NULL;
  ModelData* data = otherType->object;
  MAP_ERROR_CODE success = MAP_SAFE;

  do {
    success = free_outer_solve_data(&data->outerSolveData, zType->x_Len, map_msg, ierr); CHECKERRQ(MAP_FATAL_73);
    success = free_fortran_types(uType,pType,xType,zType,otherType,yType); 
    list_destroy(&data->yList->outList);    /* destroy output lists for writting information to output file */
    list_destroy(&data->yList->outListPtr); /* destroy output lists for writting information to output file */
    success = free_internal_states(data,map_msg,ierr); CHECKERRQ(MAP_FATAL_47);//@rm, should be replaced with a MAPFREE(data->yList)   
  } while (0);
  
  list_iterator_start(&data->element);          /* starting an iteration "session" */
  while (list_iterator_hasnext(&data->element)) { /* tell whether more values available */
    iterElem = (Element*)list_iterator_next(&data->element);
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
  list_iterator_stop(&data->element);            /* ending the iteration "session" */  
   
  list_iterator_start( &data->node);            /* starting an iteration "session" */
  while (list_iterator_hasnext(&data->node)) { /* tell whether more values available */ 
    iterNode = (Node*)list_iterator_next( &data->node );
    MAPFREE( iterNode->MApplied.name ); 
    MAPFREE( iterNode->MApplied.units );
    MAPFREE( iterNode->BApplied.name );
    MAPFREE( iterNode->BApplied.units );
    
    MAPFREE( iterNode->externalForce.fx.name );
    MAPFREE( iterNode->externalForce.fx.units );
    MAPFREE( iterNode->externalForce.fy.name );
    MAPFREE( iterNode->externalForce.fy.units );
    MAPFREE( iterNode->externalForce.fz.name );
    MAPFREE( iterNode->externalForce.fz.units );
    
    MAPFREE( iterNode->positionPtr.x.name ); 
    MAPFREE( iterNode->positionPtr.x.units );
    MAPFREE( iterNode->positionPtr.y.name ); 
    MAPFREE( iterNode->positionPtr.y.units );
    MAPFREE( iterNode->positionPtr.z.name ); 
    MAPFREE( iterNode->positionPtr.z.units );
    
    MAPFREE( iterNode->sumForcePtr.fx.name ); 
    MAPFREE( iterNode->sumForcePtr.fx.units );
    MAPFREE( iterNode->sumForcePtr.fy.name ); 
    MAPFREE( iterNode->sumForcePtr.fy.units );
    MAPFREE( iterNode->sumForcePtr.fz.name ); 
    MAPFREE( iterNode->sumForcePtr.fz.units );
  };
  list_iterator_stop(&data->node);             /* ending the iteration "session" */  
  
  list_iterator_start( &data->cableLibrary );          /* starting an iteration "session" */
  while ( list_iterator_hasnext(&data->cableLibrary) ) { /* tell whether more values available */
    iterCableLibrary = (CableLibrary*)list_iterator_next( &data->cableLibrary );
    MAPFREE( iterCableLibrary->label );
  };
  list_iterator_stop(&data->cableLibrary);             /* ending the iteration "session" */  
   
  list_destroy( &data->element );
  list_destroy( &data->node );
  list_destroy( &data->cableLibrary );
  MAPFREE( data->modelOptions.repeatAngles );
    
  /* Now delete the vessel information */
  MAPFREE(data->vessel.xi);
  MAPFREE(data->vessel.yi);
  MAPFREE(data->vessel.zi);

  MAPFREE( data->vessel.displacement.x.name );
  MAPFREE( data->vessel.displacement.x.units );
  MAPFREE( data->vessel.displacement.y.name );
  MAPFREE( data->vessel.displacement.y.units );
  MAPFREE( data->vessel.displacement.z.name );
  MAPFREE( data->vessel.displacement.z.units );
  
  MAPFREE( data->vessel.refOrigin.x.name );
  MAPFREE( data->vessel.refOrigin.x.units );
  MAPFREE( data->vessel.refOrigin.y.name );
  MAPFREE( data->vessel.refOrigin.y.units );
  MAPFREE( data->vessel.refOrigin.z.name );
  MAPFREE( data->vessel.refOrigin.z.units );
  
  MAPFREE( data->vessel.lineSumForce.fx.name );
  MAPFREE( data->vessel.lineSumForce.fx.units );
  MAPFREE( data->vessel.lineSumForce.fy.name );
  MAPFREE( data->vessel.lineSumForce.fy.units );
  MAPFREE( data->vessel.lineSumForce.fz.name );
  MAPFREE( data->vessel.lineSumForce.fz.units );
  
  MAPFREE( data->vessel.orientation.phi.name );
  MAPFREE( data->vessel.orientation.phi.units );
  MAPFREE( data->vessel.orientation.the.name );
  MAPFREE( data->vessel.orientation.the.units );
  MAPFREE( data->vessel.orientation.psi.name );
  MAPFREE( data->vessel.orientation.psi.units );

  // @rm
  // delete_input_data( u );     // MAPFREE( u );
  // delete_parameter_data( p ); // MAPFREE( p );
  // delete_continous_data( x ); // MAPFREE( x );
  // delete_constraint_data( z );// MAPFREE( z );
  // delete_output_data( y );    // MAPFREE( y );
  // delete_model_data( data );  // MAPFREE( data );
};
