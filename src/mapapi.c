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
  R[1][0] = sphi*cthe;    R[1][1] = sphi*sthe*sphi + cpsi*cphi;   R[1][2] = spsi*sthe*cphi - cpsi*sphi;
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



MAP_EXTERNCALL double** map_linearize_matrix(MAP_InputType_t* u_type, MAP_OtherStateType_t* other_type, MAP_OutputType_t* y_type, MAP_ConstraintStateType_t* z_type, double epsilon, MAP_ERROR_CODE* ierr, char* map_msg)
{
  checkpoint();
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
         success = fd_x_sequence(other_type, u_type, y_type, z_type, &force, epsilon, n, x_original, map_msg, ierr); CHECKERRQ(MAP_FATAL_62);
         success = calculate_stiffness_2(K[0], &force, epsilon, n); CHECKERRQ(MAP_FATAL_62);
       } else if (i==1) {
         success = fd_y_sequence(other_type, u_type, y_type, z_type, &force, epsilon, n, y_original, map_msg, ierr); CHECKERRQ(MAP_FATAL_63);
         success = calculate_stiffness_2(K[1], &force, epsilon, n); CHECKERRQ(MAP_FATAL_63);
       } else if (i==2) {
         success = fd_z_sequence(other_type, u_type, y_type, z_type, &force, epsilon, n, z_original, map_msg, ierr); CHECKERRQ(MAP_FATAL_64);
         success = calculate_stiffness_2(K[2], &force, epsilon, n); CHECKERRQ(MAP_FATAL_64);
       } else if (i==3) {
         success = fd_phi_sequence(other_type, u_type, y_type, z_type, &force, epsilon, n, x_original, y_original, z_original, map_msg, ierr); //CHECKERRQ(MAP_FATAL_65);
         success = calculate_stiffness_2(K[3], &force, epsilon, n); CHECKERRQ(MAP_FATAL_65);
       } else if (i==4) {
         success = fd_the_sequence(other_type, u_type, y_type, z_type, &force, epsilon, n, x_original, y_original, z_original, map_msg, ierr); //CHECKERRQ(MAP_FATAL_66);
         success = calculate_stiffness_2(K[4], &force, epsilon, n); CHECKERRQ(MAP_FATAL_66);
       } else if (i==5) {
         success = fd_psi_sequence(other_type, u_type, y_type, z_type, &force, epsilon, n, x_original, y_original, z_original, map_msg, ierr); //CHECKERRQ(MAP_FATAL_67);
         success = calculate_stiffness_2(K[5], &force, epsilon, n); CHECKERRQ(MAP_FATAL_67);
       };
     };
   } while (0);  
 
   success = reset_force_to_zero(force.fx, force.fy, force.fz, force.mx, force.my, force.mz, n);
   success = restore_original_displacement(u_type->x, x_original, n);
   success = restore_original_displacement(u_type->y, y_original, n);
   success = restore_original_displacement(u_type->z, z_original, n);
   success = line_solve_sequence(data, 0.0, map_msg, ierr); 
   
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


MAP_EXTERNCALL double* map_plot_x_array(MAP_OtherStateType_t* otherType, int i, int numPlotPoints, char *map_msg, MAP_ERROR_CODE *ierr)
{
  ModelData* data = otherType->object;
  const Element* element = (Element*)list_get_at(&data->element, i);  
  MapReal H = 0.0;
  MapReal V = 0.0;
  MapReal EA = 0.0;
  MapReal Lu = 0.0;
  MapReal Lb = 0.0;
  MapReal w = 0.0;
  MapReal S = 0.0;
  MapReal dS = 0.0;
  MapReal fairleadX=0.0;
  MapReal anchorX=0.0;
  MapReal cb = 0.0;
  MapReal lambda = 0.0;
  double* arrayX = NULL;
  char buffer[64] = "";
  int ret = 0;
  int cx = 0;
  int s = 0;
  bstring user_msg = NULL;

  *ierr = MAP_SAFE;
  map_msg[0] = 0;
  // element = (Element*)list_get_at(&data->element, i);

  if (element==NULL) {    
    user_msg = bformat("Element out of range: <%d>.", i);
    *ierr = map_set_universal_error(user_msg, map_msg, *ierr, MAP_FATAL_42);        
    ret = bdestroy(user_msg); 
    user_msg = NULL;
  } else {
    fairleadX = *(element->fairlead->positionPtr.x.value);
    anchorX = *(element->anchor->positionPtr.x.value);
    arrayX = (double*)malloc(numPlotPoints*sizeof(double));
    H = *(element->H.value);
    V = *(element->V.value);  
    EA = element->lineProperty->ea;
    Lu = element->Lu.value;   
    w = element->lineProperty->omega;
    cb = element->lineProperty->cb;
    dS = Lu/(MapReal)(numPlotPoints-1) ;
    
    /* If the cable is not resting on the seabed, we use the classic catenary equation
     * for a hanging chain to plot the mooring line profile. Otherwise if it is, we 
     * the modified version as done in the FAST wind turbine program. 
     *     
     * @ref : J. Jonkman, November 2007. "Dynamic Modeling and Loads Analysis of an 
     *        Offshore Floating Wind Turbine." NREL Technical Report NREL/TP-500-41958.
     */        
    if (element->options.omitContact==true || w<0.0 || (V-w*Lu)>0.0) { /* true when no portion of the line rests on the seabed */
      for (s=0 ; s<numPlotPoints ; s++) {
        arrayX[s] = fairleadX - ( (H/w)*ARCSINH( V/H ) - (H/w)*ARCSINH( (V-S*w)/H ) + (H*S)/(EA) )*cos(element->psi.value);
        S += dS;
      };
    } else {
      Lb = Lu - (V/w);      
      lambda = (Lb-H/(cb*w))>0 ? (Lb-H/(cb*w)) : 0; 
      for (s=0 ; s<numPlotPoints ; s++) {
        if (0<=S && S<=(Lb-H/(cb*w))) { /* for 0 <= s <= Lb - H/(Cb*w) */
          arrayX[s] = S*cos(element->psi.value) + anchorX; 
        } else if ((Lb-H/(cb/w))<S && S<=Lb) { /* for Lb - H/(Cb*w) < s <= Lb */
          arrayX[s] = (S + ((cb*w)/(2*EA)) * (S*S - 2*(Lb-H/(cb*w))*S + (Lb- H/(cb*w))*lambda))*cos(element->psi.value) + anchorX; 
        } else { /* for Lb < s <= L */
          arrayX[s] = (Lb + (H/w)*ARCSINH((w*(S-Lb))/H))*cos(element->psi.value)    
            - ( ((H*S)/(EA)) + ((cb*w)/(2*EA))*(-Lb*Lb + (Lb-H/(cb*w))*lambda))*cos(element->psi.value) + anchorX; 
        };
        S += dS;
      };
    };
  };
  return arrayX;
};


MAP_EXTERNCALL double* map_plot_y_array(MAP_OtherStateType_t* otherType, int i, int numPlotPoints, char* map_msg, MAP_ERROR_CODE* ierr)
{
  ModelData* data = otherType->object;
  Element* element = NULL;
  MapReal H = 0.0;
  MapReal V = 0.0;
  MapReal EA = 0.0;
  MapReal Lu = 0.0;
  MapReal Lb = 0.0;
  MapReal w = 0.0;
  MapReal S = 0.0;
  MapReal dS = 0.0;
  MapReal fairleadY=0.0;
  MapReal anchorY=0.0;
  MapReal cb = 0.0;
  MapReal lambda = 0.0;
  double* arrayY = NULL;
  char buffer[64] = "";
  int cx = 0;
  int s = 0;
  int ret = 0;
  bstring user_msg = NULL;
  
  *ierr = MAP_SAFE;
  map_msg[0] = 0;
  element = (Element*)list_get_at(&data->element, i);
  
  if (element==NULL) {
    user_msg = bformat("Element out of range: <%d>.", i);
    *ierr = map_set_universal_error(user_msg, map_msg, *ierr, MAP_FATAL_42);        
    ret = bdestroy(user_msg); 
    user_msg = NULL;
    //cx = map_snprintf(buffer, 64, "Element out of range: %d.", i); assert(cx>=0);
    //*ierr = map_set_universal_error(buffer, map_msg, ierr, MAP_FATAL_42);
  } else {    
    fairleadY = *(element->fairlead->positionPtr.y.value);
    anchorY = *(element->anchor->positionPtr.y.value);
    arrayY = (double*)malloc(numPlotPoints*sizeof(double));
    H = *(element->H.value);
    V = *(element->V.value);  
    EA = element->lineProperty->ea;
    Lu = element->Lu.value;
    w = element->lineProperty->omega;
    cb = element->lineProperty->cb;
    dS = Lu/(MapReal)(numPlotPoints-1) ;
    
    /* If the cable is not resting on the seabed, we use the classic catenary equation
     * for a hanging chain to plot the mooring line profile. Otherwise if it is, we 
     * the modified version as done in the FAST wind turbine program. 
     *     
     * @ref : J. Jonkman, November 2007. "Dynamic Modeling and Loads Analysis of an 
     *        Offshore Floating Wind Turbine." NREL Technical Report NREL/TP-500-41958.
     */        
    if (element->options.omitContact==true || w<0.0 || (V-w*Lu)>0.0) { /* true when no portion of the line rests on the seabed */
      for (s=0 ; s<numPlotPoints ; s++) {
        arrayY[s] = fairleadY - ( (H/w)*ARCSINH( V/H ) - (H/w)*ARCSINH( (V-S*w)/H ) + (H*S)/(EA) )*sin(element->psi.value);
        S += dS; 
      };
    } else {
      Lb = Lu - (V/w);      
      lambda = (Lb-H/(cb*w))>0 ? (Lb-H/(cb*w)) : 0; 
      for (s=0 ; s<numPlotPoints ; s++) {
        if (0<=S && S<=(Lb-H/(cb*w))) { /* for 0 <= s <= Lb - H/(Cb*w) */
          arrayY[s] = S*sin(element->psi.value) + anchorY; 
        } else if ((Lb-H/(cb/w))<S && S<=Lb) { /* for Lb - H/(Cb*w) < s <= Lb */
          arrayY[s] = (S + ((cb*w)/(2*EA))*(S*S - 2*(Lb-H/(cb*w))*S + (Lb- H/(cb*w))*lambda))*sin(element->psi.value) + anchorY; 
        } else { /* for Lb < s <= L */
          arrayY[s] = (Lb + (H/w)*ARCSINH((w*(S-Lb))/H))*sin(element->psi.value)    
            - (((H*S)/(EA)) + ((cb*w)/(2*EA))*(-Lb*Lb + (Lb-H/(cb*w))*lambda))*sin(element->psi.value) + anchorY; 
        };
        S += dS;
      };
    };
  };
  return arrayY;
};


MAP_EXTERNCALL double* map_plot_z_array(MAP_OtherStateType_t* otherType, int i, int numPlotPoints, char* map_msg, MAP_ERROR_CODE* ierr)
{
  ModelData* data = otherType->object;
  Element* element = NULL;
  MapReal H = 0.0;
  MapReal V = 0.0;
  MapReal Va = 0.0;
  MapReal EA = 0.0;
  MapReal Lu = 0.0;
  MapReal Lb = 0.0;
  MapReal w = 0.0;
  MapReal S = 0.0;
  MapReal dS = 0.0;
  MapReal fairleadZ = 0.0;
  MapReal anchorZ = 0.0;
  double* arrayZ = NULL;
  char buffer[64] = "";
  int ret = 0;
  int cx = 0;
  int s = 0;
  bstring user_msg = NULL;

  *ierr = MAP_SAFE;
  map_msg[0] = 0;
  element = (Element*)list_get_at(&data->element, i);
  
  if (element==NULL){
    user_msg = bformat("Element out of range: <%d>.", i);
    *ierr = map_set_universal_error(user_msg, map_msg, *ierr, MAP_FATAL_42);        
    ret = bdestroy(user_msg); 
    user_msg = NULL;
  } else {
    fairleadZ = *(element->fairlead->positionPtr.z.value);
    anchorZ = *(element->anchor->positionPtr.z.value);    
    arrayZ = (double*)malloc(numPlotPoints*sizeof(double));
    H = *(element->H.value);
    V = *(element->V.value);  
    EA = element->lineProperty->ea;
    Lu = element->Lu.value;
    w = element->lineProperty->omega;
    dS = Lu/(MapReal)(numPlotPoints-1) ;

    /* If the cable is not resting on the seabed, we use the classic catenary equation
     * for a hanging chain to plot the mooring line profile. Otherwise if it is, we 
     * the modified version as done in the FAST wind turbine program. 
     *     
     * @ref : J. Jonkman, November 2007. "Dynamic Modeling and Loads Analysis of an 
     *        Offshore Floating Wind Turbine." NREL Technical Report NREL/TP-500-41958.
     */        
    if (element->options.omitContact==true || w<0.0 || (V-w*Lu)>0.0) { /* true when no portion of the line rests on the seabed */
      for (s=0 ; s<numPlotPoints ; s++) {
        arrayZ[s] =  fairleadZ - ((H/w)*(sqrt(1+pow(V/H,2)) - sqrt(1+pow((V-w*S)/H,2))) + (1/EA)*(V*S+w*S*S/2)); /* Z position of element in global coordinates */
        S += dS;
      };
    } else {
      Lb = Lu - (V/w);      
      for (s=0 ; s<numPlotPoints ; s++) {        
        if (0<=S && S<=Lb) {           
          arrayZ[s] = anchorZ; 
        } else {        
          /* @todo: verify this equation before someone uses this program to  design something that matters */           
          arrayZ[s] = ((H/w)*(sqrt(1 + pow((w*(S-Lb)/H),2)) - 1) + ((w*pow((S-Lb),2))/(2*EA))) + anchorZ; 

        };
        S += dS;
      };      
    };
  };
  return arrayZ;
};


MAP_EXTERNCALL void map_plot_array_free(MapReal* array) 
{
  MAPFREE(array);
}


MAP_EXTERNCALL MapReal map_residual_function_length(MAP_OtherStateType_t* other_type, int i, char* map_msg, MAP_ERROR_CODE* ierr)
{ 
  ModelData* model_data = other_type->object;
  Element* element = NULL;
  MapReal Fh = 0.0;
  MapReal Fv = 0.0;
  MapReal EA = 0.0;
  MapReal Lu = 0.0;
  MapReal length = 0.0;
  MapReal omega = 0.0;
  MapReal cb = 0.0;
  bool contact_flag = false;
  bstring user_msg = NULL;

  map_reset_universal_error(map_msg, ierr);  
  element = (Element*)list_get_at(&model_data->element, i);

  if (element==NULL) {
    user_msg = bformat("Element out of range: %d.", i);
    *ierr = map_set_universal_error(user_msg, map_msg, *ierr, MAP_FATAL_42);
    bdestroy(user_msg);
    return -999.9;
  };

  Fh = *(element->H.value);
  Fv = *(element->V.value);  
  EA = element->lineProperty->ea;
  Lu = element->Lu.value;
  length = element->l.value;
  omega = element->lineProperty->omega;
  contact_flag = element->options.omitContact;
  cb = element->lineProperty->cb;

  if (contact_flag==true || omega<0.0 || (Fv-omega*Lu)>0.0) { /* true when no portion of the line rests on the seabed */
    return residual_function_length_no_contact(Fv, Fh, omega, Lu, EA, length);
  } else { /* true when a portion of the line rests on the seabed and the anchor tension is nonzero */
    return residual_function_length_contact(Fv, Fh, omega, Lu, EA, length, cb);
  };
};


MAP_EXTERNCALL MapReal map_residual_function_height(MAP_OtherStateType_t* other_type, int i, char* map_msg, MAP_ERROR_CODE* ierr)
{
  ModelData* model_data = other_type->object;
  Element* element = NULL;
  MapReal Fh = 0.0;
  MapReal Fv = 0.0;
  MapReal EA = 0.0;
  MapReal Lu = 0.0;
  MapReal height = 0.0;
  MapReal omega = 0.0;
  MapReal cb = 0.0;
  bool contact_flag = false;
  bstring user_msg = NULL;

  map_reset_universal_error(map_msg, ierr);  
  element = (Element*)list_get_at(&model_data->element, i);

  if (element==NULL) {    
    user_msg = bformat("Element out of range: %d.", i);
    *ierr = map_set_universal_error(user_msg, map_msg, *ierr, MAP_FATAL_42);
    bdestroy(user_msg);
    return -999.9;
  };

  Fh = *(element->H.value);
  Fv = *(element->V.value);  
  EA = element->lineProperty->ea;
  Lu = element->Lu.value;
  height = element->h.value;
  omega = element->lineProperty->omega;
  contact_flag = element->options.omitContact;
  cb = element->lineProperty->cb;

  if (contact_flag==true || omega<0.0 || (Fv-omega*Lu)>0.0) { /* true when no portion of the line rests on the seabed */
    return residual_function_height_no_contact(Fv, Fh, omega, Lu, EA, height);
  } else { /* true when a portion of the line rests on the seabed and the anchor tension is nonzero */
    return residual_function_height_contact(Fv, Fh, omega, Lu, EA, height, cb);
  };
};


MAP_EXTERNCALL MapReal map_jacobian_dxdh(MAP_OtherStateType_t* other_type, int i, char* map_msg, MAP_ERROR_CODE* ierr)
{
  ModelData* model_data = other_type->object;
  Element* element = NULL;
  MapReal Fh = 0.0;
  MapReal Fv = 0.0;
  MapReal EA = 0.0;
  MapReal Lu = 0.0;
  MapReal omega = 0.0;
  MapReal cb = 0.0;
  bool contact_flag = false;
  bstring user_msg = NULL;

  map_reset_universal_error(map_msg, ierr);  
  element = (Element*)list_get_at(&model_data->element, i);

  if (element==NULL) {    
    user_msg = bformat("Element out of range: %d.", i);
    *ierr = map_set_universal_error(user_msg, map_msg, *ierr, MAP_FATAL_42);
    bdestroy(user_msg);
    return -999.9;
  };

  Fh = *(element->H.value);
  Fv = *(element->V.value);  
  EA = element->lineProperty->ea;
  Lu = element->Lu.value;
  omega = element->lineProperty->omega;
  contact_flag = element->options.omitContact;
  cb = element->lineProperty->cb;

  if (contact_flag==true || omega<0.0 || (Fv-omega*Lu)>0.0) { /* true when no portion of the line rests on the seabed */
    return jacobian_dxdh_no_contact(Fv, Fh, omega, Lu, EA);
  } else { /* true when a portion of the line rests on the seabed and the anchor tension is nonzero */
    return jacobian_dxdh_contact(Fv, Fh, omega, Lu, EA, cb);
  };
};


MAP_EXTERNCALL MapReal map_jacobian_dxdv(MAP_OtherStateType_t* other_type, int i, char* map_msg, MAP_ERROR_CODE* ierr)
{
  ModelData* model_data = other_type->object;
  Element* element = NULL;
  MapReal Fh = 0.0;
  MapReal Fv = 0.0;
  MapReal EA = 0.0;
  MapReal Lu = 0.0;
  MapReal omega = 0.0;
  MapReal cb = 0.0;
  bool contact_flag = false;
  bstring user_msg = NULL;

  map_reset_universal_error(map_msg, ierr);  
  element = (Element*)list_get_at(&model_data->element, i);
  
  if (element==NULL) {    
    user_msg = bformat("Element out of range: %d.", i);
    *ierr = map_set_universal_error(user_msg, map_msg, *ierr, MAP_FATAL_42);
    bdestroy(user_msg);
    return -999.9;
  };

  Fh = *(element->H.value);
  Fv = *(element->V.value);  
  EA = element->lineProperty->ea;
  Lu = element->Lu.value;
  omega = element->lineProperty->omega;
  contact_flag = element->options.omitContact;
  cb = element->lineProperty->cb;
  
  if (contact_flag==true || omega<0.0 || (Fv-omega*Lu)>0.0) { /* true when no portion of the line rests on the seabed */
    return jacobian_dxdv_no_contact(Fv, Fh, omega, Lu, EA);
  } else { /* true when a portion of the line rests on the seabed and the anchor tension is nonzero */
    return jacobian_dxdv_contact(Fv, Fh, omega, Lu, EA, cb);
  };
};


MAP_EXTERNCALL MapReal map_jacobian_dzdh(MAP_OtherStateType_t* other_type, int i, char* map_msg, MAP_ERROR_CODE* ierr)
{
  ModelData* model_data = other_type->object;
  Element* element = NULL;
  MapReal Fh = 0.0;
  MapReal Fv = 0.0;
  MapReal EA = 0.0;
  MapReal Lu = 0.0;
  MapReal omega = 0.0;
  MapReal cb = 0.0;
  bool contact_flag = false;
  bstring user_msg = NULL;

  map_reset_universal_error(map_msg, ierr);  
  element = (Element*)list_get_at(&model_data->element, i);

  if (element==NULL) {    
    user_msg = bformat("Element out of range: %d.", i);
    *ierr = map_set_universal_error(user_msg, map_msg, *ierr, MAP_FATAL_42);
    bdestroy(user_msg);
    return -999.9;
  };

  Fh = *(element->H.value);
  Fv = *(element->V.value);  
  EA = element->lineProperty->ea;
  Lu = element->Lu.value;
  omega = element->lineProperty->omega;
  contact_flag = element->options.omitContact;
  cb = element->lineProperty->cb;
  
  if (contact_flag==true || omega<0.0 || (Fv-omega*Lu)>0.0) { /* true when no portion of the line rests on the seabed */
    return jacobian_dzdh_no_contact(Fv, Fh, omega, Lu, EA);
  } else { /* true when a portion of the line rests on the seabed and the anchor tension is nonzero */
    return jacobian_dzdh_contact(Fv, Fh, omega, Lu, EA, cb);
  };
};


MAP_EXTERNCALL MapReal map_jacobian_dzdv(MAP_OtherStateType_t* other_type, int i, char* map_msg, MAP_ERROR_CODE* ierr)
{
  ModelData* model_data = other_type->object;
  Element* element = NULL;
  MapReal Fh = 0.0;
  MapReal Fv = 0.0;
  MapReal EA = 0.0;
  MapReal Lu = 0.0;
  MapReal omega  = 0.0;
  MapReal cb = 0.0;
  bool contact_flag = false;
  bstring user_msg = NULL;

  map_reset_universal_error(map_msg, ierr);  
  element = (Element*)list_get_at(&model_data->element, i);

  if (element==NULL) {    
    user_msg = bformat("Element out of range: %d.", i);
    *ierr = map_set_universal_error(user_msg, map_msg, *ierr, MAP_FATAL_42);
    bdestroy(user_msg);
    return -999.9;
  };

  Fh = *(element->H.value);
  Fv = *(element->V.value);  
  EA = element->lineProperty->ea;
  Lu = element->Lu.value;
  omega  = element->lineProperty->omega;
  contact_flag = element->options.omitContact;
  cb = element->lineProperty->cb;

  if (contact_flag==true || omega<0.0 || (Fv-omega*Lu)>0.0) { /* true when no portion of the line rests on the seabed */
    return jacobian_dzdv_no_contact(Fv, Fh, omega, Lu, EA);
  } else { /* true when a portion of the line rests on the seabed and the anchor tension is nonzero */
    return jacobian_dzdv_contact(Fv, Fh, omega, Lu, EA, cb);
  };
};


// MAP_EXTERNCALL void map_get_fairlead_force_2d(double* H, double* V, MAP_OtherStateType_t* other_type, int index, char* map_msg, MAP_ERROR_CODE* ierr)
// {
//   Element* iter_element = NULL;
//   ModelData* model_data = other_type->object;
//   if (index<=list_size(&model_data->element)-1) {
//     iter_element = (Element*)list_get_at(&model_data->element, index);
//     *H = *(iter_element->H.value);
//     *V = *(iter_element->V.value);
//   } else {
//     /* @todo: throw error: element out of range */
//   };
// }


// MAP_EXTERNCALL void map_get_fairlead_force_3d(double* fx, double* fy, double* fz, MAP_OtherStateType_t* other_type, int index, char* map_msg, MAP_ERROR_CODE* ierr)
// {
//   Element* iter_element = NULL;
//   ModelData* model_data = other_type->object;
//   double psi = 0.0;
//   if (index<=list_size(&model_data->element)-1) {
//     iter_element = (Element*)list_get_at(&model_data->element, index);
//     psi = iter_element->psi.value;
//     *fx = *(iter_element->H.value)*cos(psi);
//     *fy = *(iter_element->H.value)*sin(psi);
//     *fz = *(iter_element->V.value);
//   } else {
//     /* @todo: throw error: element out of range */
//   };
// }


MAP_EXTERNCALL int map_size_elements(MAP_OtherStateType_t* other_type, MAP_ERROR_CODE* ierr, char* map_msg)
{
  ModelData* model_data = other_type->object;
  return list_size(&model_data->element);
}


MAP_EXTERNCALL void map_set_summary_file_name(MAP_InitInputType_t* init_type, char *map_msg, MAP_ERROR_CODE *ierr) 
{  
  InitializationData* init_data = init_type->object;   
  init_data->summaryFileName = bformat("%s", init_type->summaryFileName);
};


// MAP_EXTERNCALL void map_get_header_string(int* n, char** str_array, MAP_OtherStateType_t* other_type)
// { 
//   int count = 0;    
//   ModelData* model_data = other_type->object;
//   VarTypePtr* vartype_ptr = NULL;
//   VarType* vartype = NULL;
// 
//   list_iterator_start(&model_data->yList->out_list_ptr);
//   while (list_iterator_hasnext(&model_data->yList->out_list_ptr)) { 
//     vartype_ptr = (VarTypePtr*)list_iterator_next(&model_data->yList->out_list_ptr);
//     strcpy(str_array[count],vartype_ptr->name);
//     count++;
//   };
//   list_iterator_stop(&model_data->yList->out_list_ptr);     
// 
//   list_iterator_start(&model_data->yList->out_list);
//   while (list_iterator_hasnext(&model_data->yList->out_list)) { 
//     vartype = (VarType*)list_iterator_next(&model_data->yList->out_list);
//     strcpy(str_array[count],vartype->name);
//     count++;
//   };
//   list_iterator_stop(&model_data->yList->out_list);     
//   /* @todo this should raise and error when count != n */
// };


// MAP_EXTERNCALL void map_get_unit_string(int* n, char** str_array, MAP_OtherStateType_t* other_type)
// { 
//   int count = 0;    
//   ModelData* model_data = other_type->object;
//   VarTypePtr* vartype_ptr = NULL;
//   VarType* vartype = NULL;
// 
//   list_iterator_start(&model_data->yList->out_list_ptr);
//   while (list_iterator_hasnext(&model_data->yList->out_list_ptr)) { 
//     vartype_ptr = (VarTypePtr*)list_iterator_next(&model_data->yList->out_list_ptr );
//     strcpy(str_array[count],vartype_ptr->units);
//     count++;
//   };
//   list_iterator_stop(&model_data->yList->out_list_ptr);     
// 
//   list_iterator_start(&model_data->yList->out_list);
//   while (list_iterator_hasnext(&model_data->yList->out_list)) { 
//     vartype = (VarType*)list_iterator_next(&model_data->yList->out_list );
//     strcpy(str_array[count],vartype->units);
//     count++;
//   };
//   list_iterator_stop(&model_data->yList->out_list);     
// };



/**
 * Fortran binding routine
 * SUBROUTINE MAP_set_initinput_to_null(interf,msg,err) bind(C,name='set_init_to_null') 
 */
MAP_EXTERNCALL void set_init_to_null(MAP_InitInputType_t* init_type, char* map_msg, MAP_ERROR_CODE* ierr)
{
  InitializationData* init = init_type->object; 

  init->summaryFileName = NULL;  
};


/** @addtogroup FortranCall */
/* @{ */
MAP_EXTERNCALL void map_set_sea_depth(MAP_ParameterType_t* p_type, const MapReal depth)
{
  p_type->depth = depth;
};


MAP_EXTERNCALL void map_set_sea_density(MAP_ParameterType_t* p_type, const MapReal rho)
{
  p_type->rhoSea = rho;
};


MAP_EXTERNCALL void map_set_gravity(MAP_ParameterType_t* p_type, double gravity)
{
  p_type->g = gravity;
};
/* @} */


MAP_EXTERNCALL void map_add_cable_library_input_text(MAP_InitInputType_t* init_type)
{
  InitializationData* init_data = init_type->object; 
  const int n = init_data->libraryInputString->qty;
  int ret = 0;
  
  ret = bstrListAlloc(init_data->libraryInputString, n+1);
  init_data->libraryInputString->entry[n] = bfromcstr(init_type->libraryInputLine);
  init_data->libraryInputString->qty++;
};


MAP_EXTERNCALL void map_add_node_input_text(MAP_InitInputType_t* init_type)
{
  InitializationData* init_data = init_type->object; 
  const int n = init_data->nodeInputString->qty;
  int ret = 0;
  
  ret = bstrListAlloc(init_data->nodeInputString, n+1);
  init_data->nodeInputString->entry[n] = bfromcstr(init_type->nodeInputLine);
  init_data->nodeInputString->qty++;
};


MAP_EXTERNCALL void map_add_element_input_text(MAP_InitInputType_t* init_type)
{
  InitializationData* init_data = init_type->object; 
  const int n = init_data->elementInputString->qty;
  int ret = 0;
  
  ret = bstrListAlloc(init_data->elementInputString, n+1);
  init_data->elementInputString->entry[n] = bfromcstr(init_type->elementInputLine);
  init_data->elementInputString->qty++;
};


MAP_EXTERNCALL void map_add_options_input_text(MAP_InitInputType_t* init_type)
{
  InitializationData* init_data = init_type->object; 
  const int n = init_data->solverOptionsString->qty;
  int ret = 0;
  
  ret = bstrListAlloc(init_data->solverOptionsString, n+1);
  init_data->solverOptionsString->entry[n] = bfromcstr(init_type->optionInputLine);
  init_data->solverOptionsString->qty++;
};


MAP_EXTERNCALL InitializationData* MAP_InitInput_Create(char* map_msg, MAP_ERROR_CODE* ierr)
{
  InitializationData* new_data = NULL;

  map_reset_universal_error(map_msg, ierr);
  new_data = malloc(sizeof(InitializationData));
  if (new_data == NULL) {
    *ierr = map_set_universal_error(NULL, map_msg, *ierr, MAP_FATAL_4);    
    return new_data;
  } else {
    initialize_init_data_to_null(new_data);
    *ierr = MAP_SAFE;
    return new_data;
  };
};


MAP_EXTERNCALL MAP_InitInputType_t* map_create_init_type(char* map_msg, MAP_ERROR_CODE* ierr)
{
  MAP_InitInputType_t* new_data = NULL;  

  *ierr = MAP_SAFE;
  new_data = malloc(sizeof(MAP_InitInputType_t));
  if (new_data==NULL) {
    *ierr = map_set_universal_error(NULL, map_msg, *ierr, MAP_FATAL_4);    
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
    *ierr = map_set_universal_error(NULL, map_msg, *ierr, MAP_FATAL_43);    
    return new_data;
  } else {
    initialize_model_data_to_null(new_data);
    *ierr = MAP_SAFE;
    return new_data;
  };
};

MAP_EXTERNCALL MAP_OtherStateType_t* map_create_other_type(char* map_msg, MAP_ERROR_CODE* ierr)
{
  MAP_OtherStateType_t* new_data = NULL;  

  *ierr = MAP_SAFE;
  new_data = malloc(sizeof(MAP_OtherStateType_t)); 
  if (new_data==NULL) {
    *ierr = map_set_universal_error(NULL, map_msg, *ierr, MAP_FATAL_43);    
  } else {
    new_data->object = NULL;
    new_data->object = (ModelData*)(uintptr_t)MAP_OtherState_Create(map_msg, ierr);
  };
  return new_data;    
};


MAP_EXTERNCALL MAP_InitOutputType_t* map_create_initout_type(char* map_msg, MAP_ERROR_CODE* ierr)
{
  MAP_InitOutputType_t* new_type = NULL;  

  *ierr = MAP_SAFE;
  new_type = malloc(sizeof(MAP_InitOutputType_t));
  if (new_type==NULL) {
    *ierr = map_set_universal_error(NULL, map_msg, *ierr, MAP_FATAL_11);    
  } else {
    new_type->object = NULL;
  };
  return new_type;    
};


MAP_EXTERNCALL MAP_InputType_t* map_create_input_type(char* map_msg, MAP_ERROR_CODE* ierr)
{
  MAP_InputType_t* new_type = NULL;  

  *ierr = MAP_SAFE;
  new_type = malloc(sizeof(MAP_InputType_t)); 
  if (new_type == NULL) {
    *ierr = map_set_universal_error(NULL, map_msg, *ierr, MAP_FATAL_5);        
  } else { 
    new_type->object = NULL;
  };
  return new_type;    
};


MAP_EXTERNCALL MAP_ParameterType_t* map_create_parameter_type(char* map_msg, MAP_ERROR_CODE* ierr)
{
  MAP_ParameterType_t* new_type = NULL;  

  *ierr = MAP_SAFE;
  new_type = malloc(sizeof(MAP_ParameterType_t));
  if (new_type==NULL) {
    *ierr = map_set_universal_error(NULL, map_msg, *ierr, MAP_FATAL_6);    
  } else {
    new_type->object = NULL;
  };
  return new_type;
};


MAP_EXTERNCALL MAP_ConstraintStateType_t* map_create_constraint_type(char* map_msg, MAP_ERROR_CODE* ierr)
{
  MAP_ConstraintStateType_t* new_type = NULL;  

  *ierr = MAP_SAFE;
  new_type = malloc(sizeof(MAP_ConstraintStateType_t)); 
  if (new_type==NULL) {
    *ierr = map_set_universal_error(NULL, map_msg, *ierr, MAP_FATAL_8);    
  } else {
    new_type->object = NULL;
  };
  return new_type;      
};


MAP_EXTERNCALL MAP_OutputType_t* map_create_output_type(char* map_msg, MAP_ERROR_CODE* ierr)
{
  MAP_OutputType_t* new_type = NULL;  

  *ierr = MAP_SAFE;
  new_type = malloc(sizeof(MAP_OutputType_t));
  if (new_type==NULL) {
    *ierr = map_set_universal_error(NULL, map_msg, *ierr, MAP_FATAL_10);    
  } else {
    new_type->object = NULL;
  };
  return new_type;    
};


MAP_EXTERNCALL MAP_ContinuousStateType_t* map_create_continuous_type(char* map_msg, MAP_ERROR_CODE* ierr)
{
  MAP_ContinuousStateType_t* new_data = NULL;  

  *ierr = MAP_SAFE;
  new_data = malloc(sizeof(MAP_ContinuousStateType_t)); 
  if (new_data==NULL) {
    *ierr = map_set_universal_error(NULL, map_msg, *ierr, MAP_FATAL_7);    
  } else {
    new_data->object = NULL;    
  };
  return new_data;    
};
