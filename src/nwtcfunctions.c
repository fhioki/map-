/****************************************************************
 *   Copyright (C) 2014 mdm                                     *
 *   marco[dot]masciola[at]gmail                                *
 *                                                              *
 * Licensed to the Apache Software Foundation (ASF) under one   *
 * or more contributor license agreements.  See the NOTICE file *
 * distributed with this work for additional information        *
 * regarding copyright ownership.  The ASF licenses this file   *
 * to you under the Apache License, Version 2.0 (the            *
 * "License"); you may not use this file except in compliance   *
 * with the License.  You may obtain a copy of the License at   *
 *                                                              *
 *   http://www.apache.org/licenses/LICENSE-2.0                 *
 *                                                              *
 * Unless required by applicable law or agreed to in writing,   *
 * software distributed under the License is distributed on an  *
 * "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY       *
 * KIND, either express or implied.  See the License for the    *
 * specific language governing permissions and limitations      *      
 * under the License.                                           *  
 ****************************************************************/


/**
 * @file 
 * NWTC required functions necessary to hook MAP into FAST. These functions provide the binding to 
 * solve constraints for subsequent time-steps.
 */


#include "initialization.h"
#include "freedata.h"
#include "mapinit.h"
#include "nwtcfunctions.h"
#include "lineroutines.h"
#include "outputstream.h"


MAP_EXTERNCALL void map_init(MAP_InitInputType_t* init_type, 
                             MAP_InputType_t* u_type,
                             MAP_ParameterType_t* p_type,
                             MAP_ContinuousStateType_t* x_type,
                             void* none,
                             MAP_ConstraintStateType_t* z_type,
                             MAP_OtherStateType_t* other_type,
                             MAP_OutputType_t* y_type,
                             MAP_InitOutputType_t* io_type,
                             MAP_ERROR_CODE *ierr,
                             char *map_msg) {    
  InitializationData* init_data = init_type->object;   
  ModelData* model_data = other_type->object;
  MAP_ERROR_CODE success = MAP_SAFE;
  bstring user_msg = NULL;

  size_t len = strlen(map_msg);
  
  map_reset_universal_error(map_msg, ierr);
  do { 
    /*  initialize types; set doubles to -999.9, int=0, pointers=NULL 
     *  @todo: add other variables as neccessary. This needs to be fixed each time the registry
     *         is run an new varaibles are introduced 
     */  
    success = initialize_fortran_types(u_type, p_type, x_type, z_type, other_type, y_type, io_type);     
    success = map_get_version(io_type);
    
    /* create a cable library link lists for:
     *  - nodes
     *  - lines
     *  - cable library (properties)
     * The following are simclist routines 
     */
    list_init(&model_data->cableLibrary); 
    list_init(&model_data->node); 
    list_init(&model_data->element);  
    list_attributes_copy(&model_data->cableLibrary, cable_library_meter, 1); 
    list_attributes_copy(&model_data->node, node_meter, 1); 
    list_attributes_copy(&model_data->element, cable_element_meter, 1);    

    success = allocate_outlist(model_data, map_msg, ierr); CHECKERRQ(MAP_FATAL_47);
    list_init(&model_data->yList->out_list); /* simclist routine */
    list_init(&model_data->yList->out_list_ptr); /* simclist routine */

    /* The follow routines expand the input file contents based on the number of repeat
     * angles. If not repeat angles are declared, then expandedNodeInputString=nodeInputString
     * and expandedElementInputString=elementInputString. This is just a convenient way
     * to duplicate lines if a symetric mooring is employed
     */
    success = set_model_options_list(model_data, init_data, map_msg, ierr); CHECKERRQ(MAP_FATAL_33);
    success = set_cable_library_list(model_data, init_data, map_msg, ierr); CHECKERRQ(MAP_FATAL_16);

    success = repeat_nodes(model_data, init_data, map_msg, ierr);
    success = repeat_elements(model_data, init_data, map_msg, ierr);
     
    success = set_node_list(p_type, u_type, z_type, other_type, y_type, model_data, init_data->expandedNodeInputString, map_msg, ierr); CHECKERRQ(MAP_FATAL_16);    
    success = set_element_list(z_type, model_data, init_data->expandedElementInputString, map_msg, ierr); CHECKERRQ(MAP_FATAL_16);    
    
    /* now create an output list to print to and output file. */
    list_attributes_copy(&model_data->yList->out_list, vartype_meter, 1);  
    list_attributes_copy(&model_data->yList->out_list_ptr, vartype_ptr_meter, 1);  
    success = set_output_list(model_data, io_type, map_msg, ierr); 
    success = set_vessel(&model_data->vessel, u_type, map_msg, ierr); CHECKERRQ(MAP_FATAL_69);
    
    /* @todo: possible include hook for LM model here eventually place MAP_SOVE = LM */
    if (z_type->x_Len!=0) { /* this means there are no connect nodes. This does NOT mean z_type->H_len==0 */
      success = allocate_outer_solve_data(&model_data->outer_loop, z_type->x_Len, map_msg, ierr); CHECKERRQ(MAP_FATAL_72);
      model_data->MAP_SOLVE_TYPE = PARTITIONED;      
    } else {
      model_data->MAP_SOLVE_TYPE = MONOLITHIC;
    };

    /* if DEBUG is raised in CCFLAGS, then MAP version number is printed to screen */    
#   ifdef DEBUG
    print_machine_name_to_screen( );
#   endif 
     
    printf("MAP environment properties (set externally)...\n");
    printf("    Gravity constant          [m/s^2]  : %1.2f\n", p_type->g ); 
    printf("    Sea density               [kg/m^3] : %1.2f\n", p_type->rhoSea );
    printf("    Water depth               [m]      : %1.2f\n", p_type->depth );
    printf("    Vessel reference position [m]      : %1.2f , %1.2f , %1.2f\n", model_data->vessel.refOrigin.x.value, model_data->vessel.refOrigin.y.value, model_data->vessel.refOrigin.z.value); 
   
    success = initialize_cable_library_variables(model_data, p_type, map_msg, ierr); CHECKERRQ(MAP_FATAL_41);
    success = set_line_variables_pre_solve(model_data, map_msg, ierr); CHECKERRQ(MAP_FATAL_86);// @rm, not needed. This is called in line_solve_sequence
    success = reset_node_force_to_zero(model_data, map_msg, ierr); // @rm, not needed. This is called in line_solve_sequence
    success = set_element_initial_guess(model_data, map_msg, ierr);
    success = first_solve(model_data, u_type, z_type, other_type, y_type, map_msg, ierr); CHECKERRQ(MAP_FATAL_39);
    success = set_line_variables_post_solve(model_data, map_msg, ierr);    // @rm, not needed. This is called in line_solve_sequence
    success = write_summary_file(init_data, p_type, model_data, map_msg, ierr); CHECKERRQ(MAP_FATAL_37);           
    // success = get_iteration_output_stream(y_type, other_type, map_msg, ierr); // @todo CHECKERRQ()
  } while (0);  
  free_init_data(init_data, map_msg, ierr); 
  MAP_InitInput_Delete(init_data);
};


MAP_EXTERNCALL void map_update_states(double t,
                                      int interval,
                                      MAP_InputType_t* u_type,
                                      MAP_ParameterType_t* p_type,
                                      MAP_ContinuousStateType_t* x_type,
                                      void* none,
                                      MAP_ConstraintStateType_t* z_type,
                                      MAP_OtherStateType_t* other_type,
                                      MAP_ERROR_CODE* ierr,
                                      char* map_msg ) {
   ModelData* model_data = other_type->object;
   MAP_ERROR_CODE success = MAP_SAFE;
   
   map_reset_universal_error(map_msg, ierr);
   do {
     if (model_data->MAP_SOLVE_TYPE==MONOLITHIC) {
       success = line_solve_sequence(model_data, 0.0, map_msg, ierr);
     } else {
       success = node_solve_sequence(model_data, u_type, z_type, other_type, map_msg, ierr); // @todo CHECKERRQ()
     };    
   } while (0);
};    


MAP_EXTERNCALL void map_calc_output(double t,
                                    MAP_InputType_t* u_type,
                                    MAP_ParameterType_t* p_type,
                                    MAP_ContinuousStateType_t* x_type,
                                    void* none,
                                    MAP_ConstraintStateType_t* z_type,
                                    MAP_OtherStateType_t* other_type,
                                    MAP_OutputType_t* y_type,
                                    MAP_ERROR_CODE* ierr,
                                    char* map_msg ) {
   ModelData* model_data = other_type->object;
   MAP_ERROR_CODE success = MAP_SAFE;
   map_reset_universal_error(map_msg, ierr);
   success = get_iteration_output_stream(y_type, other_type, map_msg, ierr); // @todo: CHECKERRQ();
};


MAP_EXTERNCALL void map_end(MAP_InputType_t* u_type,
                            MAP_ParameterType_t* p_type,
                            MAP_ContinuousStateType_t* x_type,
                            void* none,
                            MAP_ConstraintStateType_t* z_type,
                            MAP_OtherStateType_t* other_type,
                            MAP_OutputType_t* y_type,                                                                           
                            MAP_ERROR_CODE* ierr,
                            char* map_msg ) {  
  Node* iterNode = NULL;
  CableLibrary* iterCableLibrary = NULL;
  ModelData* model_data = other_type->object;
  MAP_ERROR_CODE success = MAP_SAFE;
  bstring user_msg = NULL;

   map_reset_universal_error(map_msg, ierr);  
   do {
     success = free_outer_solve_data(&model_data->outer_loop, z_type->x_Len, map_msg, ierr); CHECKERRQ(MAP_FATAL_73);
     success = map_free_types(u_type, p_type, x_type, z_type, other_type, y_type); 
     list_destroy(&model_data->yList->out_list);    /* destroy output lists for writting information to output file */
     list_destroy(&model_data->yList->out_list_ptr); /* destroy output lists for writting information to output file */
     success = free_outlist(model_data,map_msg,ierr); CHECKERRQ(MAP_FATAL_47);//@rm, should be replaced with a MAPFREE(data->yList)   
     success = free_element(&model_data->element);
     success = free_node(&model_data->node);
     success = free_vessel(&model_data->vessel);
   
   
     list_iterator_start(&model_data->cableLibrary);          /* starting an iteration "session" */
     while ( list_iterator_hasnext(&model_data->cableLibrary)) { /* tell whether more values available */
       iterCableLibrary = (CableLibrary*)list_iterator_next(&model_data->cableLibrary);
       bdestroy(iterCableLibrary->label);
     };
     list_iterator_stop(&model_data->cableLibrary);             /* ending the iteration "session" */  
     
     list_destroy(&model_data->element);
     list_destroy(&model_data->node);
     list_destroy(&model_data->cableLibrary);
     MAPFREE(model_data->modelOptions.repeat_angle);
     MAP_OtherState_Delete(model_data);
   } while (0);
};
