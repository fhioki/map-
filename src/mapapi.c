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


#include "map.h"
#include "maperror.h"
#include "MAP_Types.h"
#include "mapapi.h"
#include "lineroutines.h"
#include "freedata.h"
#include "mapinit.h"
#include "outputstream.h"
#include "numeric.h"


extern const char MAP_ERROR_STRING[][1024];


/**
 * @file 
 * NWTC required functions necessary to hook MAP into FAST. These functions provide the binding to 
 * solve constraints for subsequent time-steps.
 */




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

  model_data->HEAD_U_TYPE = u_type;  
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
    list_init(&model_data->library); 
    list_init(&model_data->node); 
    list_init(&model_data->line);  
    list_init(&model_data->u_update_list);  
    list_attributes_copy(&model_data->library, cable_library_meter, 1); 
    list_attributes_copy(&model_data->node, node_meter, 1); 
    list_attributes_copy(&model_data->line, cable_line_meter, 1);    
    list_attributes_copy(&model_data->u_update_list, u_list_meter, 1);    

    success = allocate_outlist(model_data, map_msg, ierr); CHECKERRQ(MAP_FATAL_47);
    list_init(&model_data->y_list->out_list); /* simclist routine */
    list_init(&model_data->y_list->out_list_ptr); /* simclist routine */

    /* The follow routines expand the input file contents based on the number of repeat
     * angles. If not repeat angles are declared, then expanded_node_input_string=node_input_string
     * and expanded_line_input_string=line_input_string. This is just a convenient way
     * to duplicate lines if a symetric mooring is employed
     */
    success = set_model_options_list(model_data, init_data, map_msg, ierr); CHECKERRQ(MAP_FATAL_33);
    success = set_cable_library_list(model_data, init_data, map_msg, ierr); CHECKERRQ(MAP_FATAL_16);

    success = repeat_nodes(model_data, init_data, map_msg, ierr);
    success = repeat_lines(model_data, init_data, map_msg, ierr);
     
    success = set_node_list(p_type, u_type, z_type, other_type, y_type, model_data, init_data->expanded_node_input_string, map_msg, ierr); CHECKERRQ(MAP_FATAL_16);    
    success = set_line_list(z_type, model_data, init_data->expanded_line_input_string, map_msg, ierr); CHECKERRQ(MAP_FATAL_16);    
    
    /* now create an output list to print to and output file. */
    list_attributes_copy(&model_data->y_list->out_list, vartype_meter, 1);  
    list_attributes_copy(&model_data->y_list->out_list_ptr, vartype_ptr_meter, 1);  
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
    printf("    Sea density               [kg/m^3] : %1.2f\n", p_type->rho_sea );
    printf("    Water depth               [m]      : %1.2f\n", p_type->depth );
    printf("    Vessel reference position [m]      : %1.2f , %1.2f , %1.2f\n", model_data->vessel.ref_origin.x.value, model_data->vessel.ref_origin.y.value, model_data->vessel.ref_origin.z.value); 
   
    success = initialize_cable_library_variables(model_data, p_type, map_msg, ierr); CHECKERRQ(MAP_FATAL_41);
    success = set_line_variables_pre_solve(model_data, map_msg, ierr); CHECKERRQ(MAP_FATAL_86);// @rm, not needed. This is called in line_solve_sequence
    success = reset_node_force_to_zero(model_data, map_msg, ierr); // @rm, not needed. This is called in line_solve_sequence
    success = set_line_initial_guess(model_data, map_msg, ierr);
    success = first_solve(model_data, p_type, u_type, z_type, other_type, y_type, map_msg, ierr); CHECKERRQ(MAP_FATAL_39);
    success = set_line_variables_post_solve(model_data, map_msg, ierr);    // @rm, not needed. This is called in line_solve_sequence
    success = write_summary_file(init_data, p_type, model_data, map_msg, ierr); CHECKERRQ(MAP_FATAL_37);           
    success = get_iteration_output_stream(y_type, other_type, map_msg, ierr); // @todo CHECKERRQ()
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
   ReferencePoint* point_iter = NULL;
   Node* node_iter = NULL;
   int i = 0;
   int j = 0;
   
   map_reset_universal_error(map_msg, ierr);
   do {
     /* If the reference to u_type changes, then we have to update the location MAP internal states are pointing 
      * to. This is accomplished in the following code. The issue here is when this is called in Fortran:
      *
      *    CALL MAP_CopyInput(u(1), u_interp, MESH_NEWCOPY, ErrStat, ErrMsg)      
      *
      * u_interp is passed into into the argument for map_update_states(); however, the internal states are not
      * pointing to data in u_interp. We address this below. Note that the initial reference for point_iter is set
      * in set_node_list(...)
      */
     if (u_type!=model_data->HEAD_U_TYPE) { /* this is intended to be triggered when couled to FAST */
       list_iterator_start(&model_data->u_update_list);  
       while (list_iterator_hasnext(&model_data->u_update_list)) { 
         point_iter = (ReferencePoint*)list_iterator_next(&model_data->u_update_list);               
         point_iter->x->value = &(u_type->x[i]);
         point_iter->y->value = &(u_type->y[i]);
         point_iter->z->value = &(u_type->z[i]);                 
         i++;
       };
       list_iterator_stop(&model_data->u_update_list);

       // /* This proves the node position is updated with the Fortran interpolated input */
       // list_iterator_start(&model_data->node);  
       // while (list_iterator_hasnext(&model_data->node)) { 
       //   node_iter = (Node*)list_iterator_next(&model_data->node);               
       //   printf("After update>  %1.2f  %1.2f  %1.2f\n", *node_iter->position_ptr.x.value, *node_iter->position_ptr.y.value, *node_iter->position_ptr.z.value);
       // };
       // list_iterator_stop(&model_data->node);
       
       if (i!=u_type->x_Len) { /* raise error if the input array are exceeded */
         set_universal_error_with_message(map_msg, ierr, MAP_FATAL_89, "u_type range: <%d>. Updated array range: <%d>", u_type->x_Len, i);
         break;
       };
     };
     
     if (model_data->MAP_SOLVE_TYPE==MONOLITHIC) {
       success = line_solve_sequence(model_data, p_type, 0.0, map_msg, ierr);
     } else {
       success = node_solve_sequence(model_data, p_type, u_type, z_type, other_type, map_msg, ierr); // @todo CHECKERRQ()
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

  map_reset_universal_error(map_msg, ierr);  
  do {
    success = free_outer_solve_data(&model_data->outer_loop, z_type->x_Len, map_msg, ierr); CHECKERRQ(MAP_FATAL_73);
    success = map_free_types(u_type, p_type, x_type, z_type, other_type, y_type); 
    list_destroy(&model_data->y_list->out_list);    /* destroy output lists for writting information to output file */
    list_destroy(&model_data->y_list->out_list_ptr); /* destroy output lists for writting information to output file */
    success = free_outlist(model_data,map_msg,ierr); CHECKERRQ(MAP_FATAL_47);//@rm, should be replaced with a MAPFREE(data->y_list)   
    success = free_line(&model_data->line);
    success = free_node(&model_data->node);
    success = free_vessel(&model_data->vessel);
    success = free_cable_library(&model_data->library);
    success = free_update_list(&model_data->u_update_list);
     
    list_destroy(&model_data->line);
    list_destroy(&model_data->node);
    list_destroy(&model_data->library);
    list_destroy(&model_data->u_update_list);  
    MAPFREE(model_data->model_options.repeat_angle);
    MAP_OtherState_Delete(model_data);
  } while (0);
};


MAP_EXTERNCALL void map_offset_vessel(MAP_OtherStateType_t* other_type, MAP_InputType_t* u_type, double x, double y, double z, double phi, double the, double psi, char* map_msg, MAP_ERROR_CODE* ierr)
{
  ModelData* data = other_type->object;
  Vessel* vessel = &data->vessel;
  int i = 0;
  const int u_size = u_type->x_Len;
  double R[3][3];
  double cphi = 0.0;
  double sphi = 0.0;
  double cthe = 0.0;
  double sthe = 0.0;
  double cpsi = 0.0;
  double spsi = 0.0;
  double rx = 0.0;
  double ry = 0.0;
  double rz = 0.0;

  map_reset_universal_error(map_msg, ierr);
  
  /* define angles */
  cphi = cos(phi*DEG2RAD);
  sphi = sin(phi*DEG2RAD);
  cthe = cos(the*DEG2RAD);
  sthe = sin(the*DEG2RAD);
  cpsi = cos(psi*DEG2RAD);
  spsi = sin(psi*DEG2RAD);

  /* define transformation matrix */
  R[0][0] = cpsi*cthe;    R[0][1] = cpsi*sthe*sphi - spsi*cphi;   R[0][2] = cpsi*sthe*cphi + spsi*sphi;
  R[1][0] = spsi*cthe;    R[1][1] = sphi*sthe*sphi + cpsi*cphi;   R[1][2] = spsi*sthe*cphi - cpsi*sphi;
  R[2][0] = -sthe;        R[2][1] = cthe*sphi;                    R[2][2] = cthe*cphi;

  for (i=0 ; i<u_size ; i++) { 
    /* @todo: need to include the reference position for non-zero reference origins , i.e. r = (xi-ref) 
     *        xi, yi, and zi are the original node position. We are adding the new displacement to it.
     */
    rx = vessel->xi[i];
    ry = vessel->yi[i];
    rz = vessel->zi[i];

    /* matrix-vector product */
    u_type->x[i] = x + rx*R[0][0] + ry*R[0][1] + rz*R[0][2];
    u_type->y[i] = y + rx*R[1][0] + ry*R[1][1] + rz*R[1][2];
    u_type->z[i] = z + rx*R[2][0] + ry*R[2][1] + rz*R[2][2];
  };
};



MAP_EXTERNCALL double** map_linearize_matrix(MAP_InputType_t* u_type, MAP_ParameterType_t* p_type, MAP_OtherStateType_t* other_type, MAP_OutputType_t* y_type, MAP_ConstraintStateType_t* z_type, double epsilon, MAP_ERROR_CODE* ierr, char* map_msg)
{
  MapReal* x_original = NULL;
  MapReal* y_original = NULL;
  MapReal* z_original = NULL;
  ModelData* data = other_type->object;
  MAP_ERROR_CODE success = MAP_SAFE;
  const int n = u_type->x_Len;
  const int SIX = 6;
  int i = 0;
  int j = 0;
  int k = 0;
  Fd force;
  double** K;
 
  map_reset_universal_error(map_msg, ierr);

  K = (double**)malloc(SIX*sizeof(double*));
  for (i=0 ; i<SIX ; i++) {
    K[i] = (double*)malloc(SIX*sizeof(double));
    
    /* initialize K(6x6) allocated above to zero, row-by-row */
    K[i][0] = 0.0;
    K[i][1] = 0.0;
    K[i][2] = 0.0;
    K[i][3] = 0.0;
    K[i][4] = 0.0;
    K[i][5] = 0.0;
  };

   force.fx = (double*)malloc(n*sizeof(double));
   force.fy = (double*)malloc(n*sizeof(double));
   force.fz = (double*)malloc(n*sizeof(double));
   force.mx = (double*)malloc(n*sizeof(double));
   force.my = (double*)malloc(n*sizeof(double));
   force.mz = (double*)malloc(n*sizeof(double));  
   x_original = (double*)malloc(n*sizeof(double));
   y_original = (double*)malloc(n*sizeof(double));
   z_original = (double*)malloc(n*sizeof(double));
   
   /* initialize stuff allocated above to zero */
   for (i=0 ; i<n ; i++) {
    force.fx[i] = 0.0;
    force.fy[i] = 0.0;
    force.fz[i] = 0.0;
    force.mx[i] = 0.0;
    force.my[i] = 0.0;
    force.mz[i] = 0.0;
    x_original[i] = 0.0;
    y_original[i] = 0.0;
    z_original[i] = 0.0;
   };
    
   do {    
     /* first get the original values for the displacements */
     for (k=0 ; k<n ; k++) {
       x_original[k] = u_type->x[k];
       y_original[k] = u_type->y[k];
       z_original[k] = u_type->z[k];      
     };
    
     for (i=0 ; i<SIX ; i++) { /* down, force direction changes */
       success = reset_force_to_zero(force.fx, force.fy, force.fz, force.mx, force.my, force.mz, n);
       if (i==0) {        
         success = fd_x_sequence(other_type, p_type, u_type, y_type, z_type, &force, epsilon, n, x_original, map_msg, ierr); CHECKERRQ(MAP_FATAL_62);
         success = calculate_stiffness(K[0], &force, epsilon, n); CHECKERRQ(MAP_FATAL_62);
       } else if (i==1) {
         success = fd_y_sequence(other_type, p_type, u_type, y_type, z_type, &force, epsilon, n, y_original, map_msg, ierr); CHECKERRQ(MAP_FATAL_63);
         success = calculate_stiffness(K[1], &force, epsilon, n); CHECKERRQ(MAP_FATAL_63);
       } else if (i==2) {
         success = fd_z_sequence(other_type, p_type, u_type, y_type, z_type, &force, epsilon, n, z_original, map_msg, ierr); CHECKERRQ(MAP_FATAL_64);
         success = calculate_stiffness(K[2], &force, epsilon, n); CHECKERRQ(MAP_FATAL_64);
       } else if (i==3) {
         success = fd_phi_sequence(other_type, p_type, u_type, y_type, z_type, &force, epsilon, n, x_original, y_original, z_original, map_msg, ierr); CHECKERRQ(MAP_FATAL_65);
         success = calculate_stiffness(K[3], &force, epsilon, n); CHECKERRQ(MAP_FATAL_65);
       } else if (i==4) {
         success = fd_the_sequence(other_type, p_type, u_type, y_type, z_type, &force, epsilon, n, x_original, y_original, z_original, map_msg, ierr); CHECKERRQ(MAP_FATAL_66);
         success = calculate_stiffness(K[4], &force, epsilon, n); CHECKERRQ(MAP_FATAL_66);
       } else if (i==5) {
         success = fd_psi_sequence(other_type, p_type, u_type, y_type, z_type, &force, epsilon, n, x_original, y_original, z_original, map_msg, ierr); CHECKERRQ(MAP_FATAL_67);
         success = calculate_stiffness(K[5], &force, epsilon, n); CHECKERRQ(MAP_FATAL_67);
       };
     };
   } while (0);  
 
   success = reset_force_to_zero(force.fx, force.fy, force.fz, force.mx, force.my, force.mz, n);
   success = restore_original_displacement(u_type->x, x_original, n);
   success = restore_original_displacement(u_type->y, y_original, n);
   success = restore_original_displacement(u_type->z, z_original, n);
   success = line_solve_sequence(data, p_type, 0.0, map_msg, ierr); 
   
  MAPFREE(force.fx);
  MAPFREE(force.fy);
  MAPFREE(force.fz);
  MAPFREE(force.mx);
  MAPFREE(force.my);
  MAPFREE(force.mz);
  MAPFREE(x_original);
  MAPFREE(y_original);
  MAPFREE(z_original);  
  return K;
};


MAP_EXTERNCALL void map_free_linearize_matrix(double** array)
{
  int i = 0;
  for(i=0 ; i<6 ; i++) {
    MAPFREE(array[i]);
  };
  MAPFREE(array);
};


MAP_EXTERNCALL double* map_plot_x_array(MAP_OtherStateType_t* other_type, int i, int num_points, char *map_msg, MAP_ERROR_CODE *ierr)
{
  ModelData* model_data = other_type->object;
  Line* line = NULL;
  MapReal H = 0.0;
  MapReal V = 0.0;
  MapReal EA = 0.0;
  MapReal Lu = 0.0;
  MapReal Lb = 0.0;
  MapReal w = 0.0;
  MapReal S = 0.0;
  MapReal dS = 0.0;
  MapReal fairlead_x = 0.0;
  MapReal anchor_x = 0.0;
  MapReal cb = 0.0;
  MapReal lambda = 0.0;
  double* array_x = NULL;
  int ret = 0;
  int s = 0;

  map_reset_universal_error(map_msg, ierr);
  line = (Line*)list_get_at(&model_data->line, i);  
  
  if (line==NULL) {    
    set_universal_error_with_message(map_msg, ierr, MAP_FATAL_42, "Line out of range: <%d>.", i);
  } else if (line->options.linear_spring) {
    fairlead_x = *(line->fairlead->position_ptr.x.value);
    anchor_x = *(line->anchor->position_ptr.x.value);
    array_x = (double*)malloc(num_points*sizeof(double));
    dS = (fairlead_x-anchor_x)/(MapReal)(num_points-1);
    for (s=0 ; s<num_points ; s++) {
      array_x[s] = fairlead_x - S;
      S += dS;
    };    
  } else {
    fairlead_x = *(line->fairlead->position_ptr.x.value);
    anchor_x = *(line->anchor->position_ptr.x.value);
    array_x = (double*)malloc(num_points*sizeof(double));
    H = *(line->H.value);
    V = *(line->V.value);  
    EA = line->line_property->EA;
    Lu = line->Lu.value;   
    w = line->line_property->omega;
    cb = line->line_property->cb;
    dS = Lu/(MapReal)(num_points-1) ;
    
    /* If the cable is not resting on the seabed, we use the classic catenary equation
     * for a hanging chain to plot the mooring line profile. Otherwise if it is, we 
     * the modified version as done in the FAST wind turbine program. 
     *     
     * @ref : J. Jonkman, November 2007. "Dynamic Modeling and Loads Analysis of an 
     *        Offshore Floating Wind Turbine." NREL Technical Report NREL/TP-500-41958.
     */        
    if (line->options.omit_contact==true || w<0.0 || (V-w*Lu)>0.0) { /* true when no portion of the line rests on the seabed */
      for (s=0 ; s<num_points ; s++) {
        array_x[s] = fairlead_x - ( (H/w)*ARCSINH( V/H ) - (H/w)*ARCSINH( (V-S*w)/H ) + (H*S)/(EA) )*cos(line->psi.value);
        S += dS;
      };
    } else {
      Lb = Lu - (V/w);      
      lambda = (Lb-H/(cb*w))>0 ? (Lb-H/(cb*w)) : 0; 
      for (s=0 ; s<num_points ; s++) {
        if (0<=S && S<=(Lb-H/(cb*w))) { /* for 0 <= s <= Lb - H/(Cb*w) */
          array_x[s] = S*cos(line->psi.value) + anchor_x; 
        } else if ((Lb-H/(cb/w))<S && S<=Lb) { /* for Lb - H/(Cb*w) < s <= Lb */
          array_x[s] = (S + ((cb*w)/(2*EA)) * (S*S - 2*(Lb-H/(cb*w))*S + (Lb- H/(cb*w))*lambda))*cos(line->psi.value) + anchor_x; 
        } else { /* for Lb < s <= L */
          array_x[s] = (Lb + (H/w)*ARCSINH((w*(S-Lb))/H))*cos(line->psi.value)    
            - ( ((H*S)/(EA)) + ((cb*w)/(2*EA))*(-Lb*Lb + (Lb-H/(cb*w))*lambda))*cos(line->psi.value) + anchor_x; 
        };
        S += dS;
      };
    };
  };
  return array_x;
};


MAP_EXTERNCALL double* map_plot_y_array(MAP_OtherStateType_t* other_type, int i, int num_points, char* map_msg, MAP_ERROR_CODE* ierr)
{
  ModelData* data = other_type->object;
  Line* line = NULL;
  MapReal H = 0.0;
  MapReal V = 0.0;
  MapReal EA = 0.0;
  MapReal Lu = 0.0;
  MapReal Lb = 0.0;
  MapReal w = 0.0;
  MapReal S = 0.0;
  MapReal dS = 0.0;
  MapReal fairlead_y = 0.0;
  MapReal anchor_y = 0.0;
  MapReal cb = 0.0;
  MapReal lambda = 0.0;
  double* array_y = NULL;
  int s = 0;
  int ret = 0;

  map_reset_universal_error(map_msg, ierr);  
  line = (Line*)list_get_at(&data->line, i);
  
  if (line==NULL) {
    set_universal_error_with_message(map_msg, ierr, MAP_FATAL_42, "Line out of range: <%d>.", i);
  } else if (line->options.linear_spring) {
    fairlead_y = *(line->fairlead->position_ptr.y.value);
    anchor_y = *(line->anchor->position_ptr.y.value);
    array_y = (double*)malloc(num_points*sizeof(double));
    dS = (fairlead_y-anchor_y)/(MapReal)(num_points-1);
    for (s=0 ; s<num_points ; s++) {
      array_y[s] = fairlead_y - S;
      S += dS;
    };    
  } else {    
    fairlead_y = *(line->fairlead->position_ptr.y.value);
    anchor_y = *(line->anchor->position_ptr.y.value);
    array_y = (double*)malloc(num_points*sizeof(double));
    H = *(line->H.value);
    V = *(line->V.value);  
    EA = line->line_property->EA;
    Lu = line->Lu.value;
    w = line->line_property->omega;
    cb = line->line_property->cb;
    dS = Lu/(MapReal)(num_points-1) ;
    
    /* If the cable is not resting on the seabed, we use the classic catenary equation
     * for a hanging chain to plot the mooring line profile. Otherwise if it is, we 
     * the modified version as done in the FAST wind turbine program. 
     *     
     * @ref : J. Jonkman, November 2007. "Dynamic Modeling and Loads Analysis of an 
     *        Offshore Floating Wind Turbine." NREL Technical Report NREL/TP-500-41958.
     */        
    if (line->options.omit_contact==true || w<0.0 || (V-w*Lu)>0.0) { /* true when no portion of the line rests on the seabed */
      for (s=0 ; s<num_points ; s++) {
        array_y[s] = fairlead_y - ( (H/w)*ARCSINH( V/H ) - (H/w)*ARCSINH( (V-S*w)/H ) + (H*S)/(EA) )*sin(line->psi.value);
        S += dS; 
      };
    } else {
      Lb = Lu - (V/w);      
      lambda = (Lb-H/(cb*w))>0 ? (Lb-H/(cb*w)) : 0; 
      for (s=0 ; s<num_points ; s++) {
        if (0<=S && S<=(Lb-H/(cb*w))) { /* for 0 <= s <= Lb - H/(Cb*w) */
          array_y[s] = S*sin(line->psi.value) + anchor_y; 
        } else if ((Lb-H/(cb/w))<S && S<=Lb) { /* for Lb - H/(Cb*w) < s <= Lb */
          array_y[s] = (S + ((cb*w)/(2*EA))*(S*S - 2*(Lb-H/(cb*w))*S + (Lb- H/(cb*w))*lambda))*sin(line->psi.value) + anchor_y; 
        } else { /* for Lb < s <= L */
          array_y[s] = (Lb + (H/w)*ARCSINH((w*(S-Lb))/H))*sin(line->psi.value)    
            - (((H*S)/(EA)) + ((cb*w)/(2*EA))*(-Lb*Lb + (Lb-H/(cb*w))*lambda))*sin(line->psi.value) + anchor_y; 
        };
        S += dS;
      };
    };
  };
  return array_y;
};


MAP_EXTERNCALL double* map_plot_z_array(MAP_OtherStateType_t* other_type, int i, int num_points, char* map_msg, MAP_ERROR_CODE* ierr)
{
  ModelData* data = other_type->object;
  Line* line = NULL;
  MapReal H = 0.0;
  MapReal V = 0.0;
  MapReal Va = 0.0;
  MapReal EA = 0.0;
  MapReal Lu = 0.0;
  MapReal Lb = 0.0;
  MapReal w = 0.0;
  MapReal S = 0.0;
  MapReal dS = 0.0;
  MapReal fairlead_z = 0.0;
  MapReal anchor_z = 0.0;
  double* array_z = NULL;
  int ret = 0;
  int s = 0;

  map_reset_universal_error(map_msg, ierr);  
  line = (Line*)list_get_at(&data->line, i);
  
  if (line==NULL){
    set_universal_error_with_message(map_msg, ierr, MAP_FATAL_42, "Line out of range: <%d>.", i);
  } else if (line->options.linear_spring) {
    fairlead_z = *(line->fairlead->position_ptr.z.value);
    anchor_z = *(line->anchor->position_ptr.z.value);
    array_z = (double*)malloc(num_points*sizeof(double));
    dS = fabs(fairlead_z-anchor_z)/(MapReal)(num_points-1);
    for (s=0 ; s<num_points ; s++) {
      array_z[s] = fairlead_z - S;
      S += dS;
    };    
  } else {
    fairlead_z = *(line->fairlead->position_ptr.z.value);
    anchor_z = *(line->anchor->position_ptr.z.value);    
    array_z = (double*)malloc(num_points*sizeof(double));
    H = *(line->H.value);
    V = *(line->V.value);  
    EA = line->line_property->EA;
    Lu = line->Lu.value;
    w = line->line_property->omega;
    dS = Lu/(MapReal)(num_points-1) ;

    /* If the cable is not resting on the seabed, we use the classic catenary equation
     * for a hanging chain to plot the mooring line profile. Otherwise if it is, we 
     * the modified version as done in the FAST wind turbine program. 
     *     
     * @ref : J. Jonkman, November 2007. "Dynamic Modeling and Loads Analysis of an 
     *        Offshore Floating Wind Turbine." NREL Technical Report NREL/TP-500-41958.
     */        
    if (line->options.omit_contact==true || w<0.0 || (V-w*Lu)>0.0) { /* true when no portion of the line rests on the seabed */
      for (s=0 ; s<num_points ; s++) {
        array_z[s] =  fairlead_z - ((H/w)*(sqrt(1+pow(V/H,2)) - sqrt(1+pow((V-w*S)/H,2))) + (1/EA)*(V*S+w*S*S/2)); /* Z position of line in global coordinates */
        S += dS;
      };
    } else {
      Lb = Lu - (V/w);      
      for (s=0 ; s<num_points ; s++) {        
        if (0<=S && S<=Lb) {           
          array_z[s] = anchor_z; 
        } else {        
          /* @todo: verify this equation before someone uses this program to  design something that matters */           
          array_z[s] = ((H/w)*(sqrt(1 + pow((w*(S-Lb)/H),2)) - 1) + ((w*pow((S-Lb),2))/(2*EA))) + anchor_z; 
        };
        S += dS;
      };      
    };
  };
  return array_z;
};


MAP_EXTERNCALL void map_plot_array_free(MapReal* array) 
{
  MAPFREE(array);
}


MAP_EXTERNCALL MapReal map_residual_function_length(MAP_OtherStateType_t* other_type, int i, char* map_msg, MAP_ERROR_CODE* ierr)
{ 
  ModelData* model_data = other_type->object;
  Line* line = NULL;
  MapReal Fh = 0.0;
  MapReal Fv = 0.0;
  MapReal EA = 0.0;
  MapReal Lu = 0.0;
  MapReal length = 0.0;
  MapReal omega = 0.0;
  MapReal cb = 0.0;
  bool contact_flag = false;

  map_reset_universal_error(map_msg, ierr);  
  line = (Line*)list_get_at(&model_data->line, i);

  if (line==NULL) {
    set_universal_error_with_message(map_msg, ierr, MAP_FATAL_42, "Line out of range: <%d>.", i);
    return -999.9;
  };

  Fh = *(line->H.value);
  Fv = *(line->V.value);  
  EA = line->line_property->EA;
  Lu = line->Lu.value;
  length = line->l.value;
  omega = line->line_property->omega;
  contact_flag = line->options.omit_contact;
  cb = line->line_property->cb;

  if (contact_flag==true || omega<0.0 || (Fv-omega*Lu)>0.0) { /* true when no portion of the line rests on the seabed */
    return residual_function_length_no_contact(Fv, Fh, omega, Lu, EA, length);
  } else { /* true when a portion of the line rests on the seabed and the anchor tension is nonzero */
    return residual_function_length_contact(Fv, Fh, omega, Lu, EA, length, cb);
  };
};


MAP_EXTERNCALL MapReal map_residual_function_height(MAP_OtherStateType_t* other_type, int i, char* map_msg, MAP_ERROR_CODE* ierr)
{
  ModelData* model_data = other_type->object;
  Line* line = NULL;
  MapReal Fh = 0.0;
  MapReal Fv = 0.0;
  MapReal EA = 0.0;
  MapReal Lu = 0.0;
  MapReal height = 0.0;
  MapReal omega = 0.0;
  MapReal cb = 0.0;
  bool contact_flag = false;


  map_reset_universal_error(map_msg, ierr);  
  line = (Line*)list_get_at(&model_data->line, i);

  if (line==NULL) {    
    set_universal_error_with_message(map_msg, ierr, MAP_FATAL_42, "Line out of range: <%d>.", i);
    return -999.9;
  };

  Fh = *(line->H.value);
  Fv = *(line->V.value);  
  EA = line->line_property->EA;
  Lu = line->Lu.value;
  height = line->h.value;
  omega = line->line_property->omega;
  contact_flag = line->options.omit_contact;
  cb = line->line_property->cb;

  if (contact_flag==true || omega<0.0 || (Fv-omega*Lu)>0.0) { /* true when no portion of the line rests on the seabed */
    return residual_function_height_no_contact(Fv, Fh, omega, Lu, EA, height);
  } else { /* true when a portion of the line rests on the seabed and the anchor tension is nonzero */
    return residual_function_height_contact(Fv, Fh, omega, Lu, EA, height, cb);
  };
};


MAP_EXTERNCALL MapReal map_jacobian_dxdh(MAP_OtherStateType_t* other_type, int i, char* map_msg, MAP_ERROR_CODE* ierr)
{
  ModelData* model_data = other_type->object;
  Line* line = NULL;
  MapReal Fh = 0.0;
  MapReal Fv = 0.0;
  MapReal EA = 0.0;
  MapReal Lu = 0.0;
  MapReal omega = 0.0;
  MapReal cb = 0.0;
  bool contact_flag = false;

  map_reset_universal_error(map_msg, ierr);  
  line = (Line*)list_get_at(&model_data->line, i);

  if (line==NULL) {    
    set_universal_error_with_message(map_msg, ierr, MAP_FATAL_42, "Line out of range: <%d>.", i);
    return -999.9;
  };

  Fh = *(line->H.value);
  Fv = *(line->V.value);  
  EA = line->line_property->EA;
  Lu = line->Lu.value;
  omega = line->line_property->omega;
  contact_flag = line->options.omit_contact;
  cb = line->line_property->cb;

  if (contact_flag==true || omega<0.0 || (Fv-omega*Lu)>0.0) { /* true when no portion of the line rests on the seabed */
    return jacobian_dxdh_no_contact(Fv, Fh, omega, Lu, EA);
  } else { /* true when a portion of the line rests on the seabed and the anchor tension is nonzero */
    return jacobian_dxdh_contact(Fv, Fh, omega, Lu, EA, cb);
  };
};


MAP_EXTERNCALL MapReal map_jacobian_dxdv(MAP_OtherStateType_t* other_type, int i, char* map_msg, MAP_ERROR_CODE* ierr)
{
  ModelData* model_data = other_type->object;
  Line* line = NULL;
  MapReal Fh = 0.0;
  MapReal Fv = 0.0;
  MapReal EA = 0.0;
  MapReal Lu = 0.0;
  MapReal omega = 0.0;
  MapReal cb = 0.0;
  bool contact_flag = false;

  map_reset_universal_error(map_msg, ierr);  
  line = (Line*)list_get_at(&model_data->line, i);
  
  if (line==NULL) {    
    set_universal_error_with_message(map_msg, ierr, MAP_FATAL_42, "Line out of range: <%d>.", i);
    return -999.9;
  };

  Fh = *(line->H.value);
  Fv = *(line->V.value);  
  EA = line->line_property->EA;
  Lu = line->Lu.value;
  omega = line->line_property->omega;
  contact_flag = line->options.omit_contact;
  cb = line->line_property->cb;
  
  if (contact_flag==true || omega<0.0 || (Fv-omega*Lu)>0.0) { /* true when no portion of the line rests on the seabed */
    return jacobian_dxdv_no_contact(Fv, Fh, omega, Lu, EA);
  } else { /* true when a portion of the line rests on the seabed and the anchor tension is nonzero */
    return jacobian_dxdv_contact(Fv, Fh, omega, Lu, EA, cb);
  };
};


MAP_EXTERNCALL MapReal map_jacobian_dzdh(MAP_OtherStateType_t* other_type, int i, char* map_msg, MAP_ERROR_CODE* ierr)
{
  ModelData* model_data = other_type->object;
  Line* line = NULL;
  MapReal Fh = 0.0;
  MapReal Fv = 0.0;
  MapReal EA = 0.0;
  MapReal Lu = 0.0;
  MapReal omega = 0.0;
  MapReal cb = 0.0;
  bool contact_flag = false;

  map_reset_universal_error(map_msg, ierr);  
  line = (Line*)list_get_at(&model_data->line, i);

  if (line==NULL) {    
    set_universal_error_with_message(map_msg, ierr, MAP_FATAL_42, "Line out of range: <%d>.", i);
    return -999.9;
  };

  Fh = *(line->H.value);
  Fv = *(line->V.value);  
  EA = line->line_property->EA;
  Lu = line->Lu.value;
  omega = line->line_property->omega;
  contact_flag = line->options.omit_contact;
  cb = line->line_property->cb;
  
  if (contact_flag==true || omega<0.0 || (Fv-omega*Lu)>0.0) { /* true when no portion of the line rests on the seabed */
    return jacobian_dzdh_no_contact(Fv, Fh, omega, Lu, EA);
  } else { /* true when a portion of the line rests on the seabed and the anchor tension is nonzero */
    return jacobian_dzdh_contact(Fv, Fh, omega, Lu, EA, cb);
  };
};


MAP_EXTERNCALL MapReal map_jacobian_dzdv(MAP_OtherStateType_t* other_type, int i, char* map_msg, MAP_ERROR_CODE* ierr)
{
  ModelData* model_data = other_type->object;
  Line* line = NULL;
  MapReal Fh = 0.0;
  MapReal Fv = 0.0;
  MapReal EA = 0.0;
  MapReal Lu = 0.0;
  MapReal omega  = 0.0;
  MapReal cb = 0.0;
  bool contact_flag = false;

  map_reset_universal_error(map_msg, ierr);  
  line = (Line*)list_get_at(&model_data->line, i);

  if (line==NULL) {    
    set_universal_error_with_message(map_msg, ierr, MAP_FATAL_42, "Line out of range: <%d>.", i);
    return -999.9;
  };

  Fh = *(line->H.value);
  Fv = *(line->V.value);  
  EA = line->line_property->EA;
  Lu = line->Lu.value;
  omega  = line->line_property->omega;
  contact_flag = line->options.omit_contact;
  cb = line->line_property->cb;

  if (contact_flag==true || omega<0.0 || (Fv-omega*Lu)>0.0) { /* true when no portion of the line rests on the seabed */
    return jacobian_dzdv_no_contact(Fv, Fh, omega, Lu, EA);
  } else { /* true when a portion of the line rests on the seabed and the anchor tension is nonzero */
    return jacobian_dzdv_contact(Fv, Fh, omega, Lu, EA, cb);
  };
};


MAP_EXTERNCALL void map_get_fairlead_force_2d(double* H, double* V, MAP_OtherStateType_t* other_type, int index, char* map_msg, MAP_ERROR_CODE* ierr)
{
  Line* iter_line = NULL;
  ModelData* model_data = other_type->object;

  map_reset_universal_error(map_msg, ierr);  

  if (index<=list_size(&model_data->line)-1) {
    iter_line = (Line*)list_get_at(&model_data->line, index);
    *H = *(iter_line->H.value);
    *V = *(iter_line->V.value);
  } else {
    /* throw error: line out of range */
    set_universal_error_with_message(map_msg, ierr, MAP_FATAL_42, "Line out of range: %d.", index);
  };
}


MAP_EXTERNCALL void map_get_fairlead_force_3d(double* fx, double* fy, double* fz, MAP_OtherStateType_t* other_type, int index, char* map_msg, MAP_ERROR_CODE* ierr)
{
  Line* iter_line = NULL;
  ModelData* model_data = other_type->object;
  double psi = 0.0;

  if (index<=list_size(&model_data->line)-1) {
    iter_line = (Line*)list_get_at(&model_data->line, index);
    psi = iter_line->psi.value;
    *fx = *(iter_line->H.value)*cos(psi);
    *fy = *(iter_line->H.value)*sin(psi);
    *fz = *(iter_line->V.value);
  } else {
    /* throw error: line out of range */
    set_universal_error_with_message(map_msg, ierr, MAP_FATAL_42, "Line out of range: %d.", index);
  };
}


MAP_EXTERNCALL int map_size_lines(MAP_OtherStateType_t* other_type, MAP_ERROR_CODE* ierr, char* map_msg)
{
  ModelData* model_data = other_type->object;
  return list_size(&model_data->line);
}


MAP_EXTERNCALL void map_set_summary_file_name(MAP_InitInputType_t* init_type, char *map_msg, MAP_ERROR_CODE *ierr) 
{  
  InitializationData* init_data = init_type->object;   
  init_data->summary_file_name = bformat("%s", init_type->summary_file_name);
};


MAP_EXTERNCALL void map_get_header_string(int* n, char** str_array, MAP_OtherStateType_t* other_type)
{ 
  int count = 0;    
  ModelData* model_data = other_type->object;
  VarTypePtr* vartype_ptr = NULL;
  VarType* vartype = NULL;

  list_iterator_start(&model_data->y_list->out_list_ptr);
  while (list_iterator_hasnext(&model_data->y_list->out_list_ptr)) { 
    vartype_ptr = (VarTypePtr*)list_iterator_next(&model_data->y_list->out_list_ptr);
    strcpy(str_array[count],vartype_ptr->name->data);
    count++;
  };
  list_iterator_stop(&model_data->y_list->out_list_ptr);     

  list_iterator_start(&model_data->y_list->out_list);
  while (list_iterator_hasnext(&model_data->y_list->out_list)) { 
    vartype = (VarType*)list_iterator_next(&model_data->y_list->out_list);
    strcpy(str_array[count],vartype->name->data);
    count++;
  };
  list_iterator_stop(&model_data->y_list->out_list);     
  /* @todo this should raise and error when count != n */
};


MAP_EXTERNCALL void map_get_unit_string(int* n, char** str_array, MAP_OtherStateType_t* other_type)
{ 
  int count = 0;    
  ModelData* model_data = other_type->object;
  VarTypePtr* vartype_ptr = NULL;
  VarType* vartype = NULL;

  list_iterator_start(&model_data->y_list->out_list_ptr);
  while (list_iterator_hasnext(&model_data->y_list->out_list_ptr)) { 
    vartype_ptr = (VarTypePtr*)list_iterator_next(&model_data->y_list->out_list_ptr );
    strcpy(str_array[count],vartype_ptr->units->data);
    count++;
  };
  list_iterator_stop(&model_data->y_list->out_list_ptr);     

  list_iterator_start(&model_data->y_list->out_list);
  while (list_iterator_hasnext(&model_data->y_list->out_list)) { 
    vartype = (VarType*)list_iterator_next(&model_data->y_list->out_list );
    strcpy(str_array[count],vartype->units->data);
    count++;
  };
  list_iterator_stop(&model_data->y_list->out_list);     
};



/**
 * Fortran binding routine
 * SUBROUTINE MAP_set_initinput_to_null(interf,msg,err) bind(C,name='set_init_to_null') 
 */
MAP_EXTERNCALL void set_init_to_null(MAP_InitInputType_t* init_type, char* map_msg, MAP_ERROR_CODE* ierr)
{
  InitializationData* init = init_type->object; 

  init->summary_file_name = NULL;  
};


/** @addtogroup FortranCall */
/* @{ */
MAP_EXTERNCALL void map_set_sea_depth(MAP_ParameterType_t* p_type, const MapReal depth)
{
  p_type->depth = depth;
};


MAP_EXTERNCALL void map_set_sea_density(MAP_ParameterType_t* p_type, const MapReal rho)
{
  p_type->rho_sea = rho;
};


MAP_EXTERNCALL void map_set_gravity(MAP_ParameterType_t* p_type, double gravity)
{
  p_type->g = gravity;
};
/* @} */


MAP_EXTERNCALL void map_add_cable_library_input_text(MAP_InitInputType_t* init_type)
{
  InitializationData* init_data = init_type->object; 
  const int n = init_data->library_input_string->qty;
  int ret = 0;
  
  ret = bstrListAlloc(init_data->library_input_string, n+1);
  init_data->library_input_string->entry[n] = bfromcstr(init_type->library_input_str);
  init_data->library_input_string->qty++;
};


MAP_EXTERNCALL void map_add_node_input_text(MAP_InitInputType_t* init_type)
{
  InitializationData* init_data = init_type->object; 
  const int n = init_data->node_input_string->qty;
  int ret = 0;
  
  ret = bstrListAlloc(init_data->node_input_string, n+1);
  init_data->node_input_string->entry[n] = bfromcstr(init_type->node_input_str);
  init_data->node_input_string->qty++;
};


MAP_EXTERNCALL void map_add_line_input_text(MAP_InitInputType_t* init_type)
{
  InitializationData* init_data = init_type->object; 
  const int n = init_data->line_input_string->qty;
  int ret = 0;
  
  ret = bstrListAlloc(init_data->line_input_string, n+1);
  init_data->line_input_string->entry[n] = bfromcstr(init_type->line_input_str);
  init_data->line_input_string->qty++;
};


MAP_EXTERNCALL void map_add_options_input_text(MAP_InitInputType_t* init_type)
{
  InitializationData* init_data = init_type->object; 
  const int n = init_data->solver_options_string->qty;
  int ret = 0;
  
  ret = bstrListAlloc(init_data->solver_options_string, n+1);
  init_data->solver_options_string->entry[n] = bfromcstr(init_type->option_input_str);
  init_data->solver_options_string->qty++;
};


MAP_EXTERNCALL InitializationData* MAP_InitInput_Create(char* map_msg, MAP_ERROR_CODE* ierr)
{
  InitializationData* new_data = NULL;

  map_reset_universal_error(map_msg, ierr);
  new_data = malloc(sizeof(InitializationData));
  if (new_data == NULL) {
    set_universal_error(map_msg, ierr, MAP_FATAL_4);    
    return new_data;
  } else {
    initialize_init_data_to_null(new_data);    
    return new_data;
  };
};


MAP_EXTERNCALL MAP_InitInputType_t* map_create_init_type(char* map_msg, MAP_ERROR_CODE* ierr)
{
  MAP_InitInputType_t* new_data = NULL;  

  map_reset_universal_error(map_msg, ierr);
  new_data = malloc(sizeof(MAP_InitInputType_t));
  if (new_data==NULL) {
    set_universal_error(map_msg, ierr, MAP_FATAL_4);    
  } else {    
    initialize_init_type_to_null(new_data); /* set F2C types to null */
    new_data->object = NULL;
    new_data->object = (InitializationData*)(uintptr_t)MAP_InitInput_Create(map_msg, ierr);
  };
  return new_data;
};


MAP_EXTERNCALL ModelData* MAP_OtherState_Create(char* map_msg, MAP_ERROR_CODE* ierr)
{
  ModelData* new_data = NULL;

  map_reset_universal_error(map_msg, ierr);
  new_data = malloc(sizeof(ModelData));
  if (new_data==NULL) {
    set_universal_error(map_msg, ierr, MAP_FATAL_43);    
    return new_data;
  } else {
    initialize_model_data_to_null(new_data);    
    return new_data;
  };
};

MAP_EXTERNCALL MAP_OtherStateType_t* map_create_other_type(char* map_msg, MAP_ERROR_CODE* ierr)
{
  MAP_OtherStateType_t* new_data = NULL;  

  map_reset_universal_error(map_msg, ierr);
  new_data = malloc(sizeof(MAP_OtherStateType_t)); 
  if (new_data==NULL) {
    set_universal_error(map_msg, ierr, MAP_FATAL_43);    
  } else {
    new_data->object = NULL;
    new_data->object = (ModelData*)(uintptr_t)MAP_OtherState_Create(map_msg, ierr);
  };
  return new_data;    
};


MAP_EXTERNCALL MAP_InitOutputType_t* map_create_initout_type(char* map_msg, MAP_ERROR_CODE* ierr)
{
  MAP_InitOutputType_t* new_type = NULL;  

  map_reset_universal_error(map_msg, ierr);
  new_type = malloc(sizeof(MAP_InitOutputType_t));
  if (new_type==NULL) {
    set_universal_error(map_msg, ierr, MAP_FATAL_11);    
  } else {
    new_type->object = NULL;
  };
  return new_type;    
};


MAP_EXTERNCALL MAP_InputType_t* map_create_input_type(char* map_msg, MAP_ERROR_CODE* ierr)
{
  MAP_InputType_t* new_type = NULL;  

  map_reset_universal_error(map_msg, ierr);
  new_type = malloc(sizeof(MAP_InputType_t)); 
  if (new_type == NULL) {
    set_universal_error(map_msg, ierr, MAP_FATAL_5);        
  } else { 
    new_type->object = NULL;
  };
  return new_type;    
};


MAP_EXTERNCALL MAP_ParameterType_t* map_create_parameter_type(char* map_msg, MAP_ERROR_CODE* ierr)
{
  MAP_ParameterType_t* new_type = NULL;  

  map_reset_universal_error(map_msg, ierr);
  new_type = malloc(sizeof(MAP_ParameterType_t));
  if (new_type==NULL) {
    set_universal_error(map_msg, ierr, MAP_FATAL_6);    
  } else {
    new_type->object = NULL;
  };
  return new_type;
};


MAP_EXTERNCALL MAP_ConstraintStateType_t* map_create_constraint_type(char* map_msg, MAP_ERROR_CODE* ierr)
{
  MAP_ConstraintStateType_t* new_type = NULL;  

  map_reset_universal_error(map_msg, ierr);
  new_type = malloc(sizeof(MAP_ConstraintStateType_t)); 
  if (new_type==NULL) {
    set_universal_error(map_msg, ierr, MAP_FATAL_8);    
  } else {
    new_type->object = NULL;
  };
  return new_type;      
};


MAP_EXTERNCALL MAP_OutputType_t* map_create_output_type(char* map_msg, MAP_ERROR_CODE* ierr)
{
  MAP_OutputType_t* new_type = NULL;  

  map_reset_universal_error(map_msg, ierr);
  new_type = malloc(sizeof(MAP_OutputType_t));
  if (new_type==NULL) {
    set_universal_error(map_msg, ierr, MAP_FATAL_10);    
  } else {
    new_type->object = NULL;
  };
  return new_type;    
};


MAP_EXTERNCALL MAP_ContinuousStateType_t* map_create_continuous_type(char* map_msg, MAP_ERROR_CODE* ierr)
{
  MAP_ContinuousStateType_t* new_data = NULL;  

  map_reset_universal_error(map_msg, ierr);
  new_data = malloc(sizeof(MAP_ContinuousStateType_t)); 
  if (new_data==NULL) {
    set_universal_error(map_msg, ierr, MAP_FATAL_7);    
  } else {
    new_data->object = NULL;    
  };
  return new_data;    
};


MAP_EXTERNCALL void MAP_InitInput_Delete(InitializationData* init_data)
{
  MAPFREE(init_data); 
};


MAP_EXTERNCALL int free_init_data(InitializationData* init_data, char* map_msg, MAP_ERROR_CODE* ierr) 
{
  int ret = 0;

  ret = bdestroy(init_data->summary_file_name);
  ret = bstrListDestroy(init_data->expanded_node_input_string);
  ret = bstrListDestroy(init_data->expanded_line_input_string);
  ret = bstrListDestroy(init_data->library_input_string);  
  ret = bstrListDestroy(init_data->node_input_string);
  ret = bstrListDestroy(init_data->line_input_string);
  ret = bstrListDestroy(init_data->solver_options_string);


  return MAP_SAFE;
};
