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


#include "lineroutines.h"
#include "mapinit.h"

extern const char MAP_ERROR_STRING[][1024];


MAP_ERROR_CODE reset_node_force_to_zero(ModelData* model_data, char* map_msg, MAP_ERROR_CODE* ierr)
{
  Node* node_iter = NULL;

  list_iterator_start(&model_data->node);            /* starting an iteration "session" */
  while (list_iterator_hasnext(&model_data->node)) { /* tell whether more values available */
    node_iter = (Node*)list_iterator_next(&model_data->node);
    *(node_iter->sumForcePtr.fx.value) = 0.0;
    *(node_iter->sumForcePtr.fy.value) = 0.0;
    *(node_iter->sumForcePtr.fz.value) = 0.0;    
  };
  list_iterator_stop(&model_data->node); /* ending the iteration "session" */    
  return MAP_SAFE;
};


MAP_ERROR_CODE calculate_node_sum_force(ModelData* model_data)
{
  Element* element_iter = NULL;
  Node* node_iter = NULL;
  double psi = 0.0;
  double fx = 0.0;
  double fy = 0.0;
  double fz = 0.0;
  double fx_a = 0.0;
  double fy_a = 0.0;
  double fz_a = 0.0;

  list_iterator_start(&model_data->element);         /* starting an iteration "session" */
  while (list_iterator_hasnext(&model_data->element)) { /* tell whether more values available */
    element_iter = (Element*)list_iterator_next(&model_data->element);
    psi = element_iter->psi.value;
    fx = *(element_iter->H.value)*cos(psi);
    fx_a = -(element_iter->HAtAnchor.value)*cos(psi);
  
    fy = *(element_iter->H.value)*sin(psi);
    fy_a = -(element_iter->HAtAnchor.value)*sin(psi);
    
    fz = *(element_iter->V.value);
    fz_a = -(element_iter->VAtAnchor.value);
  
    add_to_sum_fx(element_iter->fairlead, fx);
    add_to_sum_fx(element_iter->anchor, fx_a);
  
    add_to_sum_fy(element_iter->fairlead, fy);
    add_to_sum_fy(element_iter->anchor, fy_a);
  
    add_to_sum_fz(element_iter->fairlead, fz);
    add_to_sum_fz(element_iter->anchor, fz_a);
  }; 
  list_iterator_stop(&model_data->element); /* ending the iteration "session" */    

  list_iterator_start(&model_data->node);            /* starting an iteration "session" */
  while (list_iterator_hasnext(&model_data->node)) { /* tell whether more values available */
    node_iter = (Node*)list_iterator_next(&model_data->node);
    if (node_iter->type==CONNECT) {
      add_to_sum_fx(node_iter, -(node_iter->externalForce.fx.value));
      add_to_sum_fy(node_iter, -(node_iter->externalForce.fy.value));
      add_to_sum_fz(node_iter, -(node_iter->externalForce.fz.value));
     };
  }; 
  list_iterator_stop(&model_data->node); /* ending the iteration "session" */    
  return MAP_SAFE;
}


void add_to_sum_fx(Node* node, const MapReal fx)
{
  *(node->sumForcePtr.fx.value) += fx;
};


void add_to_sum_fy(Node* node, const MapReal fy)
{
  *(node->sumForcePtr.fy.value) += fy;
};


void add_to_sum_fz(Node* node, const MapReal fz)
{
  *(node->sumForcePtr.fz.value) += fz;
};


MAP_ERROR_CODE set_line_variables_pre_solve(ModelData* model_data, char* map_msg, MAP_ERROR_CODE* ierr)
{
  MAP_ERROR_CODE success = MAP_SAFE;
  Element* element_iter = NULL;
  bstring user_msg = NULL;
  int i = 0;

  list_iterator_start(&model_data->element);            /* starting an iteration "session" */
  while (list_iterator_hasnext(&model_data->element)) { /* tell whether more values available */ 
    element_iter = (Element*)list_iterator_next(&model_data->element);    
    
    /* horizontal cable excursion*/
    element_iter->l.value = set_horizontal_excursion(element_iter);

    /* vertical cable excursion */
    element_iter->h.value = set_vertical_excursion(element_iter);

    /* angle of azimuth */
    success = set_psi(element_iter, map_msg, ierr); 
    if (success!=MAP_SAFE) {
      user_msg = bformat("Element number %d", i);
      *ierr = map_set_universal_error(user_msg, map_msg, *ierr, MAP_WARNING_6);        
      success = bdestroy(user_msg); 
      user_msg = NULL;
    };
    i++;
  };
  list_iterator_stop(&model_data->element); /* ending the iteration "session" */    
  return MAP_SAFE;
};


MAP_ERROR_CODE set_line_variables_post_solve(ModelData* model_data, char* map_msg, MAP_ERROR_CODE* ierr)
{
  Element* element_iter = NULL;
  double H = 0.0;
  double V = 0.0;
  double Ha = 0.0;
  double Va = 0.0;
  double w = 0.0;
  double Lu = 0.0;
  double Lb = 0.0;
  double cb = 0.0;
  bool contact_flag = false;
   
  list_iterator_start(&model_data->element); /* starting an iteration "session" */
  while (list_iterator_hasnext(&model_data->element)) { /* tell whether more values available */ 
    element_iter = (Element*)list_iterator_next(&model_data->element);
    
    /* altitude angle at fairlead */
    H = *(element_iter->H.value);
    V = *(element_iter->V.value);
    element_iter->alpha.value = atan2(V, H);
    
    /* tension at fairlead */
    element_iter->T.value = sqrt(H*H + V*V);
     
    /* line in contact with seabed options */
    w = element_iter->lineProperty->omega;
    cb = element_iter->lineProperty->cb;
    contact_flag = element_iter->options.omitContact;
    Lu = element_iter->Lu.value;
    if (contact_flag==true || w<0.0 || (V-w*Lu)>0.0) { 
      element_iter->lb.value = 0.0;
      Ha = H;
      Va = V - w*Lu;
    } else { /* line is touching the seabed */
      Lb = Lu - (V/w);
      element_iter->lb.value = Lb;
      ((H-cb*w*Lb)>0.0) ? (Ha = (H-cb*w*Lb)) : (Ha = 0.0);
      Va = 0.0;
    };
    
    /* tension at anchor */
    element_iter->HAtAnchor.value = Ha;
    element_iter->VAtAnchor.value = Va;
    element_iter->TAtAnchor.value = sqrt(Ha*Ha + Va*Va);
    element_iter->alphaAtAnchor.value = atan2(Va, Ha);    
  };
  list_iterator_stop(&model_data->element); /* ending the iteration "session" */        
  return MAP_SAFE;
};


MAP_ERROR_CODE set_psi(Element* element, char* map_msg, MAP_ERROR_CODE* ierr)
{
  double overlapValue = 0.0;
  MapReal fairleadX = *(element->fairlead->positionPtr.x.value);
  MapReal fairleadY = *(element->fairlead->positionPtr.y.value);
  MapReal fairleadZ = *(element->fairlead->positionPtr.z.value);
  MapReal anchorX = *(element->anchor->positionPtr.x.value);
  MapReal anchorY = *(element->anchor->positionPtr.y.value);
  MapReal anchorZ = *(element->anchor->positionPtr.z.value);

  overlapValue = sqrt(pow((fairleadX-anchorX),2) + pow((fairleadY-anchorY),2) +  pow((fairleadZ-anchorZ),2));
  
  /* make sure the demoninator is not zero (or close to it) */
  if (overlapValue<=1e-2) {    
    // @todo: what happens when this returns prematurely?
    return MAP_WARNING;
  };
  
  /* 
     find the angle psi element simply finds the angle psi between the local and global reference frames simply by 
     evaluating trig relationships
  */
  element->psi.value = atan2((fairleadY-anchorY), (fairleadX-anchorX));
  return MAP_SAFE;
};


MapReal set_horizontal_excursion(Element* element)
{
  const MapReal fairleadX = *(element->fairlead->positionPtr.x.value);
  const MapReal fairleadY = *(element->fairlead->positionPtr.y.value);
  const MapReal anchorX = *(element->anchor->positionPtr.x.value);
  const MapReal anchorY = *(element->anchor->positionPtr.y.value);
  return sqrt(pow((fairleadX-anchorX),2) + pow((fairleadY-anchorY),2));
};


MapReal set_vertical_excursion(Element* element)
{
  const MapReal fairleadZ = *(element->fairlead->positionPtr.z.value);
  const MapReal anchorZ = *(element->anchor->positionPtr.z.value);
  return fabs(fairleadZ - anchorZ);
};


MAP_ERROR_CODE set_element_initial_guess(ModelData* model_data, char* map_msg, MAP_ERROR_CODE* ierr)
{
  Element* element_iter = NULL;
  MapReal lambda = 0.0;
  MapReal length = 0.0;
  MapReal height = 0.0;
  MapReal w = 0.0;
  MapReal Lu = 0.0;

  list_iterator_start(&model_data->element);            /* starting an iteration "session" */
  while (list_iterator_hasnext(&model_data->element)) { /* tell whether more values available */ 
    element_iter = (Element*)list_iterator_next(&model_data->element);    
    w = element_iter->lineProperty->omega;
    length = element_iter->l.value;
    height = element_iter->h.value;
    Lu = element_iter->Lu.value;
    
    /* @todo: this is legacy from MAP v > 1.0 for the case of a perfectly vertical cable        
     * if ( norm <= TOL ) {
     *   psi = 1.0;//checkpoint();
     * } else {
     *   if ( (fairlead->Y.value - anchor->Y.value) >= 0) {
     *     // this simply finds the angle psi between the local and global reference frames simply by evaluation trig relationships
     *     psi = acos( (fairlead->X.value - anchor->X.value)/norm );
     *   } else {
     *     psi = -acos( (fairlead->X.value - anchor->X.value)/norm );
     *   }
     * }
     */

    /* note: the element horizontal (l) and vertical (h) excursion are previously initialized in set_line_variables_pre_solve(...) */ 
    if (length<=1e-2) {
      lambda = 1000000;
    } else if (sqrt(length*length+height*height)>=Lu) {
      lambda = 0.2;
    } else {
      lambda = sqrt(3*((Lu*Lu - height*height)/(length*length) - 1));
    };

    /* note: omega is previously initialized in function initialize_cable_library_variables(.. ) in mapinit.c */    
    *(element_iter->H.value) = fabs(w*length/(2*lambda));
    *(element_iter->V.value) = (w/2)*((height/tanh(lambda)) + Lu);
  };
  list_iterator_stop(&model_data->element); /* ending the iteration "session" */    
  return MAP_SAFE;
};



MAP_ERROR_CODE line_solve_sequence(ModelData* model_data, double t, char* map_msg, MAP_ERROR_CODE* ierr) 
{
  MAP_ERROR_CODE success = MAP_SAFE;

  do { 
    success = set_line_variables_pre_solve(model_data, map_msg, ierr);
    success = reset_node_force_to_zero(model_data, map_msg, ierr);    
    success = solve_line(model_data, t, map_msg, ierr); 
    // /* prematurely terminating the the line solve routine to return back to the
    //    caller function. Do a back tracking on the solution and attempt recovery */
    // if (success==MAP_FATAL) {
    //   return MAP_FATAL_59;
    // };   
    // CHECKERRQ(MAP_FATAL_39); /* @todo: this should be called elsewhere to notifiy 
    //                             users ofpremature termination. This won't be caught
    //                             when solve_line fails. */    
    success = set_line_variables_post_solve(model_data, map_msg, ierr);
    success = calculate_node_sum_force(model_data);
  } while (0);

  MAP_RETURN;
};    


MAP_ERROR_CODE node_solve_sequence(ModelData* model_data, MAP_InputType_t* u_type, MAP_ConstraintStateType_t* z_type, MAP_OtherStateType_t* other_type, char* map_msg, MAP_ERROR_CODE* ierr)
{
  OuterSolveAttributes* ns = &model_data->outer_loop;
  MAP_ERROR_CODE success = MAP_SAFE;
  Element* element_iter = NULL;
  const int THREE = 3;
  const int z_size = z_type->z_Len; //N
  const int m = THREE*(other_type->Fz_connect_Len); /* rows */
  const int n = THREE*(z_type->z_Len);              /* columns */
  double error = 0.0;
  int SIZE = THREE*z_size;
  int col = 0;
  int row = 0;
  int i = 0;
  int j = 0;
  int lineCounter = 0;
  bstring user_msg = NULL;

  ns->iterationCount = 1;
  do {
    error = 0.0;
    success = line_solve_sequence(model_data, 0.0, map_msg, ierr); CHECKERRQ(MAP_FATAL_79);
    switch (ns->fd) {
    case BACKWARD_DIFFERENCE :
      success = backward_difference_jacobian(other_type, z_type, model_data, map_msg, ierr); CHECKERRQ(MAP_FATAL_75);
      break;
    case CENTRAL_DIFFERENCE :
      success = central_difference_jacobian(other_type, z_type, model_data, map_msg, ierr); CHECKERRQ(MAP_FATAL_76);
      break;
    case FORWARD_DIFFERENCE :
      success = forward_difference_jacobian(other_type, z_type, model_data, map_msg, ierr); CHECKERRQ(MAP_FATAL_77);
      break;
    };
    success = line_solve_sequence(model_data, 0.0, map_msg, ierr);
    success = lu(ns, SIZE, map_msg, ierr); CHECKERRQ(MAP_FATAL_74);
    success = lu_back_substitution(ns, SIZE, map_msg, ierr); CHECKERRQ(MAP_FATAL_74);
    
    /* Note that: ns->x = J^(-1) * F
     *  [x,y,z]_i+1 =  [x,y,z]_i - J^(-1) * F        
     */   
    for (i=0 ; i<z_size ; i++) { 
      z_type->x[i] -= ns->x[THREE*i];
      z_type->y[i] -= ns->x[THREE*i+1];
      z_type->z[i] -= ns->x[THREE*i+2];
      error += (pow(other_type->Fx_connect[i],2)+ pow(other_type->Fy_connect[i],2) + pow(other_type->Fz_connect[i],2));
    };
    ns->iterationCount++;
    if (ns->iterationCount>ns->maxIts) {
      *ierr = map_set_universal_error(user_msg, map_msg, *ierr, MAP_FATAL_80);
      break;
    };
    
    // printf("Error: %f, tol=%f\n",sqrt(error),ns->tol);
    
    /* @todo: end when iterations is exceeded. need some way to indicate that simulation did not suuficiently 
     * meet termination criteria
     */
  } while (sqrt(error)>ns->tol);

  if (success==MAP_SAFE) {
    return MAP_SAFE;
  } else if (success==MAP_ERROR) {
    return MAP_ERROR;
  } else {
    return MAP_FATAL;
  };
};


MAP_ERROR_CODE solve_line(ModelData* model_data, double time, char* map_msg, MAP_ERROR_CODE* ierr)
{
  MAP_ERROR_CODE success = MAP_SAFE;
  Element* element_iter = NULL;
  bstring user_msg = NULL;
  int n = 1; /* line counter */

  list_iterator_start(&model_data->element);            /* starting an iteration "session" */
  while (list_iterator_hasnext(&model_data->element)) { /* tell whether more values available */ 
    element_iter = (Element*)list_iterator_next(&model_data->element);

    /* check for fatal errors 
     * @todo: this should handle the exception of a perfectly vertical cable. Added to the priority list
     */    
    if (element_iter->l.value<=0.0) {
      user_msg = bformat("Line segment %d, l = %d [m].", n, element_iter->l.value);
      *ierr = map_set_universal_error(user_msg, map_msg, *ierr, MAP_FATAL_54);        
      success = bdestroy(user_msg); 
      user_msg = NULL;
      break; 
    } else if (element_iter->h.value<=-MACHINE_EPSILON) {
      user_msg = bformat("Line segment %d, h = %d [m].", n, element_iter->h.value);
      *ierr = map_set_universal_error(user_msg, map_msg, *ierr, MAP_FATAL_55);        
      success = bdestroy(user_msg); 
      user_msg = NULL;
      break; 
    } else if (element_iter->lineProperty->omega>0.0) {
      success = check_maximum_line_length(element_iter, element_iter->options.omitContact, map_msg, ierr);
      if (success) {        
        user_msg = bformat("Line segment %d.", n);
        *ierr = map_set_universal_error(user_msg, map_msg, *ierr, MAP_FATAL_59);        
        success = bdestroy(user_msg); 
        user_msg = NULL;
        break;
      };
    };    
    success = call_minpack_lmder(element_iter, &model_data->inner_loop, &model_data->modelOptions, n, time, map_msg, ierr); CHECKERRQ(MAP_FATAL_79);
    n++;
  };
  list_iterator_stop(&model_data->element); /* ending the iteration "session" */    

  if (*ierr==MAP_SAFE) {
    return MAP_SAFE;
  } else if (*ierr==MAP_ERROR) {
    return MAP_ERROR;
  } else {
    return MAP_FATAL;
  };
};


MAP_ERROR_CODE check_maximum_line_length(Element* element, const bool contact_flag, char *map_msg, MAP_ERROR_CODE *ierr)
{
  MAP_ERROR_CODE success = MAP_SAFE;
  MapReal LMax = 0.0;
  MapReal l = element->l.value;
  MapReal h = element->h.value;
  MapReal EA = element->lineProperty->ea;
  MapReal W = element->lineProperty->omega;
  MapReal Lu = element->Lu.value;
  bstring user_msg = NULL;

  LMax = l - EA/W + sqrt(pow((EA/W),2) + 2.0*h*EA/W);
  if (Lu>=LMax && contact_flag==false) {
    user_msg = bformat("LMax = %f [m].", LMax);
    // printf("%s\n",user_msg->data);
    *ierr = map_set_universal_error(user_msg, map_msg, *ierr, MAP_FATAL_59);        
    success = bdestroy(user_msg); 
    user_msg = NULL;
    return MAP_FATAL;
  };
  return MAP_SAFE;
};


// MAP_ERROR_CODE increment_dof_by_delta(double* inputType, const MapReal delta, const int size)
// {
//   int i = 0;
//   for (i=0 ; i<size ; i++) {
//     inputType[i] += delta;
//   };
//   return MAP_SAFE;
// };


// MAP_ERROR_CODE restore_original_displacement(double* inputType, const double* initialValue, const int size)
// {
//   int i = 0;
//   for (i=0 ; i<size ; i++) {
//     inputType[i] = initialValue[i];
//   };
//   return MAP_SAFE;
// };


// MAP_ERROR_CODE reset_force_to_zero(double* fx, double* fy, double* fz, double* mx, double* my, double* mz, const int size)
// {
//   int i = 0;
//   for (i=0 ; i<size ; i++) {
//     fx[i] = 0.0;
//     fy[i] = 0.0;
//     fz[i] = 0.0;
//     mx[i] = 0.0;
//     my[i] = 0.0;
//     mz[i] = 0.0;
//   };
//   return MAP_SAFE;
// };



// MAP_ERROR_CODE set_force_minus(const double* inputType, double* force, const int size)
// {
//   int i = 0;
//   for (i=0 ; i<size ; i++) {
//     force[i] -= -inputType[i];
//   };
//   return MAP_SAFE;
// };


// MAP_ERROR_CODE set_force_plus(const double* inputType, double* force, const int size)
// {
//   int i = 0;
//   for (i=0 ; i<size ; i++) {
//     force[i] += -inputType[i];
//   };
//   return MAP_SAFE;
// };



// MAP_ERROR_CODE set_moment_minus(const MAP_OutputType_t* outputType, const Vessel* vessel, double* mx, double* my, double* mz, const int size)
// {
//   int i = 0;
//   double rx = 0.0;
//   double ry = 0.0;
//   double rz = 0.0;
// 
//   for (i=0 ; i<size ; i++) {    
//     /* 
//        @todo: this is not consistent with set_moment_minus_2... 
//               ensure moments are in the global frame
//     */
//     rx = (vessel->xi[i] - vessel->refOrigin.x.value);
//     ry = (vessel->yi[i] - vessel->refOrigin.y.value);
//     rz = (vessel->zi[i] - vessel->refOrigin.z.value);    
//     mx[i] -= ((-outputType->Fz[i]*ry) - (-outputType->Fy[i]*rz));
//     my[i] -= ((-outputType->Fx[i]*rz) - (-outputType->Fz[i]*rx));
//     mz[i] -= ((-outputType->Fy[i]*rx) - (-outputType->Fx[i]*ry));
//   };
//   return MAP_SAFE;
// }



// MAP_ERROR_CODE set_moment_plus(const MAP_OutputType_t* outputType, const Vessel* vessel, double* mx, double* my, double* mz, const int size)
// {
//   int i = 0;
//   double rx = 0.0;
//   double ry = 0.0;
//   double rz = 0.0;
// 
//   for (i=0 ; i<size ; i++) {    
//     rx = (vessel->xi[i] - vessel->refOrigin.x.value);
//     ry = (vessel->yi[i] - vessel->refOrigin.y.value);
//     rz = (vessel->zi[i] - vessel->refOrigin.z.value);    
//     mx[i] += ((-outputType->Fz[i]*ry) - (-outputType->Fy[i]*rz));
//     my[i] += ((-outputType->Fx[i]*rz) - (-outputType->Fz[i]*rx));
//     mz[i] += ((-outputType->Fy[i]*rx) - (-outputType->Fx[i]*ry));
//   };
//   return MAP_SAFE;
// };



// MAP_ERROR_CODE increment_phi_dof_by_delta(MAP_InputType_t* uType, const Vessel* vessel, const double delta, const int size)
// {
//   int i = 0;
//   double R[3][3];
//   double phi = vessel->orientation.phi.value;
//   double rx = 0.0;
//   double ry = 0.0;
//   double rz = 0.0;
//   R[0][0] = 1.0;   R[0][1] = 0.0;              R[0][2] = 0.0;
//   R[1][0] = 0.0;   R[1][1] = cos(phi+delta);   R[1][2] = -sin(phi+delta);  
//   R[2][0] = 0.0;   R[2][1] = sin(phi+delta);   R[2][2] = cos(phi+delta);
// 
//   for (i=0 ; i<size ; i++) {
//     /* @todo: should also include the reference origin location? */
//     rx = vessel->xi[i];
//     ry = vessel->yi[i];
//     rz = vessel->zi[i];
//     // uType->y[i] = ry - rz*delta;
//     // uType->z[i] = rz + ry*delta;
//     // printf("%f  %f  %f\n",uType->x[i],uType->y[i],uType->z[i]);
//     uType->x[i] = R[0][0]*rx;
//     uType->y[i] = R[1][1]*ry + R[1][2]*rz;
//     uType->z[i] = R[2][1]*ry + R[2][2]*rz;
//   };
//   return MAP_SAFE;
// };


// MAP_ERROR_CODE increment_the_dof_by_delta(MAP_InputType_t* uType, const Vessel* vessel, const double delta, const int size)
// {
//   int i = 0;
//   double R[3][3];
//   double the = vessel->orientation.the.value;
//   double rx = 0.0;
//   double ry = 0.0;
//   double rz = 0.0;
//   R[0][0] = cos(the+delta);   R[0][1] = 0.0;   R[0][2] = sin(the+delta);
//   R[1][0] = 0.0;              R[1][1] = 1.0;   R[1][2] = 0.0;
//   R[2][0] = -sin(the+delta);  R[2][1] = 0.0;   R[2][2] = cos(the+delta);
//   
//   for (i=0 ; i<size ; i++) {
//     /* @todo: should also include the reference origin location? */
//     rx = vessel->xi[i];
//     ry = vessel->yi[i];
//     rz = vessel->zi[i];
//     // uType->x[i] = rx + rz*delta;
//     // uType->z[i] = rz - rx*delta;    
//     uType->x[i] = R[0][0]*rx + R[0][2]*rz;
//     uType->y[i] = R[1][1]*ry;
//     uType->z[i] = R[2][0]*rx + R[2][2]*rz;
//   };
//   return MAP_SAFE;
// };


// MAP_ERROR_CODE increment_psi_dof_by_delta(MAP_InputType_t* uType, const Vessel* vessel, const double delta, const int size)
// {
//   int i = 0;
//   double R[3][3];
//   double psi = vessel->orientation.psi.value;
//   double rx = 0.0;
//   double ry = 0.0;
//   double rz = 0.0;
//   R[0][0] = cos(psi+delta);   R[0][1] = -sin(psi+delta);   R[0][2] = 0.0;
//   R[1][0] = sin(psi+delta);   R[1][1] = cos(psi+delta);    R[1][2] = 0.0;
//   R[2][0] = 0.0;              R[2][1] = 0.0;               R[2][2] = 1.0;
// 
//   for (i=0 ; i<size ; i++) {
//     /* @todo: should also include the reference origin location? */
//     rx = vessel->xi[i];
//     ry = vessel->yi[i];
//     rz = vessel->zi[i];
//     // uType->x[i] = rx - ry*delta;
//     // uType->y[i] = ry + rx*delta;    
//     uType->x[i] = R[0][0]*rx + R[0][1]*ry;
//     uType->y[i] = R[1][0]*rx + R[1][1]*ry;
//     uType->z[i] = R[2][2]*rz;
//   };
//   return MAP_SAFE;
// };
// 


// MAP_ERROR_CODE fd_x_sequence(MAP_OtherStateType_t* otherType, MAP_InputType_t* uType, MAP_OutputType_t* yType, MAP_ConstraintStateType_t* zType, Fd* force, const double epsilon, const int size, const double* originalPos, char* map_msg, MAP_ERROR_CODE* ierr)
// {
//   MAP_ERROR_CODE success = MAP_SAFE;
//   ModelData* data = otherType->object;
//   Vessel* vessel = &data->vessel;
// 
//   do {
//     /* minus epsilon sequence */
//     success = increment_dof_by_delta(uType->x, -epsilon, size); CHECKERRQ(MAP_FATAL_61);
//     if (data->MAP_SOLVE_TYPE==MONOLITHIC) {
//       success = line_solve_sequence(data, 0.0, map_msg, ierr);
//     } else {
//       success = node_solve_sequence(data, uType, zType, otherType, map_msg, ierr); // @todo CHECKERRQ()
//     };    
//     success = set_force_plus(yType->Fx, force->fx, size); CHECKERRQ(MAP_FATAL_61);
//     success = set_force_plus(yType->Fy, force->fy, size); CHECKERRQ(MAP_FATAL_61);
//     success = set_force_plus(yType->Fz, force->fz, size); CHECKERRQ(MAP_FATAL_61);
//     success = set_moment_plus(yType, vessel, force->mx, force->my, force->mz, size); CHECKERRQ(MAP_FATAL_61);
//     success = restore_original_displacement(uType->x, originalPos, size); CHECKERRQ(MAP_FATAL_61);
//   
//     /* plus epsilon sequence */
//     success = increment_dof_by_delta(uType->x, epsilon, size); CHECKERRQ(MAP_FATAL_61);
//     if (data->MAP_SOLVE_TYPE==MONOLITHIC) {
//       success = line_solve_sequence(data, 0.0, map_msg, ierr);
//     } else {
//       success = node_solve_sequence(data, uType, zType, otherType, map_msg, ierr); // @todo CHECKERRQ()
//     };    
//     success = set_force_minus(yType->Fx, force->fx, size); CHECKERRQ(MAP_FATAL_61);
//     success = set_force_minus(yType->Fy, force->fy, size); CHECKERRQ(MAP_FATAL_61);
//     success = set_force_minus(yType->Fz, force->fz, size); CHECKERRQ(MAP_FATAL_61);
//     success = set_moment_minus(yType, vessel, force->mx, force->my, force->mz, size); CHECKERRQ(MAP_FATAL_61);
//     success = restore_original_displacement(uType->x, originalPos, size); CHECKERRQ(MAP_FATAL_61);
//   } while (0);
//   MAP_RETURN;
// };


// MAP_ERROR_CODE fd_y_sequence(MAP_OtherStateType_t* otherType, MAP_InputType_t* uType, MAP_OutputType_t* yType, MAP_ConstraintStateType_t* zType, Fd* force, const double epsilon, const int size, const double* originalPos, char* map_msg, MAP_ERROR_CODE* ierr)
// {
//   MAP_ERROR_CODE success = MAP_SAFE;
//   ModelData* data = otherType->object;
//   Vessel* vessel = &data->vessel;
// 
//   do {
//     /* minus epsilon sequence */
//     success = increment_dof_by_delta(uType->y, -epsilon, size); CHECKERRQ(MAP_FATAL_61);
//     if (data->MAP_SOLVE_TYPE==MONOLITHIC) {
//       success = line_solve_sequence(data, 0.0, map_msg, ierr);
//     } else {
//       success = node_solve_sequence(data, uType, zType, otherType, map_msg, ierr); // @todo CHECKERRQ()
//     };    
//     success = set_force_plus(yType->Fx, force->fx, size); CHECKERRQ(MAP_FATAL_61);
//     success = set_force_plus(yType->Fy, force->fy, size); CHECKERRQ(MAP_FATAL_61);
//     success = set_force_plus(yType->Fz, force->fz, size); CHECKERRQ(MAP_FATAL_61);
//     success = set_moment_plus(yType, vessel, force->mx, force->my, force->mz, size); CHECKERRQ(MAP_FATAL_61);
//     success = restore_original_displacement(uType->y, originalPos, size); CHECKERRQ(MAP_FATAL_61);
//         
//     /* plus epsilon sequence */
//     success = increment_dof_by_delta(uType->y, epsilon, size); CHECKERRQ(MAP_FATAL_61);
//     if (data->MAP_SOLVE_TYPE==MONOLITHIC) {
//       success = line_solve_sequence(data, 0.0, map_msg, ierr);
//     } else {
//       success = node_solve_sequence(data, uType, zType, otherType, map_msg, ierr); // @todo CHECKERRQ()
//     };    
//     success = set_force_minus(yType->Fx, force->fx, size); CHECKERRQ(MAP_FATAL_61);
//     success = set_force_minus(yType->Fy, force->fy, size); CHECKERRQ(MAP_FATAL_61);
//     success = set_force_minus(yType->Fz, force->fz, size); CHECKERRQ(MAP_FATAL_61);
//     success = set_moment_minus(yType, vessel, force->mx, force->my, force->mz, size); CHECKERRQ(MAP_FATAL_61);
//     success = restore_original_displacement(uType->y, originalPos, size); CHECKERRQ(MAP_FATAL_61);
//   } while (0);
//   
//   MAP_RETURN;
// };


// MAP_ERROR_CODE fd_z_sequence(MAP_OtherStateType_t* otherType, MAP_InputType_t* uType, MAP_OutputType_t* yType, MAP_ConstraintStateType_t* zType, Fd* force, const double epsilon, const int size, const double* originalPos, char* map_msg, MAP_ERROR_CODE* ierr)
// {
//   MAP_ERROR_CODE success = MAP_SAFE;
//   ModelData* data = otherType->object;
//   Vessel* vessel = &data->vessel;
// 
//   do {
//     /* minus epsilon sequence */
//     success = increment_dof_by_delta(uType->z, -epsilon, size); CHECKERRQ(MAP_FATAL_61);
//     if (data->MAP_SOLVE_TYPE==MONOLITHIC) {
//       success = line_solve_sequence(data, 0.0, map_msg, ierr);
//     } else {
//       success = node_solve_sequence(data, uType, zType, otherType, map_msg, ierr); // @todo CHECKERRQ()
//     };    
//     success = set_force_plus(yType->Fx, force->fx, size); CHECKERRQ(MAP_FATAL_61);
//     success = set_force_plus(yType->Fy, force->fy, size); CHECKERRQ(MAP_FATAL_61);
//     success = set_force_plus(yType->Fz, force->fz, size); CHECKERRQ(MAP_FATAL_61);
//     success = set_moment_plus(yType, vessel, force->mx, force->my, force->mz, size); CHECKERRQ(MAP_FATAL_61);
//     success = restore_original_displacement(uType->z, originalPos, size); CHECKERRQ(MAP_FATAL_61);
//         
//     /* plus epsilon sequence */
//     success = increment_dof_by_delta(uType->z, epsilon, size); CHECKERRQ(MAP_FATAL_61);
//     if (data->MAP_SOLVE_TYPE==MONOLITHIC) {
//       success = line_solve_sequence(data, 0.0, map_msg, ierr);
//     } else {
//       success = node_solve_sequence(data, uType, zType, otherType, map_msg, ierr); // @todo CHECKERRQ()
//     };    
//     success = set_force_minus(yType->Fx, force->fx, size); CHECKERRQ(MAP_FATAL_61);
//     success = set_force_minus(yType->Fy, force->fy, size); CHECKERRQ(MAP_FATAL_61);
//     success = set_force_minus(yType->Fz, force->fz, size); CHECKERRQ(MAP_FATAL_61);
//     success = set_moment_minus(yType, vessel, force->mx, force->my, force->mz, size); CHECKERRQ(MAP_FATAL_61);
//     success = restore_original_displacement(uType->z, originalPos, size); CHECKERRQ(MAP_FATAL_61);
//   } while (0);
// 
//   MAP_RETURN;
// };


// MAP_ERROR_CODE fd_phi_sequence(MAP_OtherStateType_t* otherType, MAP_InputType_t* uType, MAP_OutputType_t* yType, MAP_ConstraintStateType_t* zType, Fd* force, const double epsilon, const int size, const double* originalX, const double* originalY, const double* originalZ, char* map_msg, MAP_ERROR_CODE* ierr)
// {
//   MAP_ERROR_CODE success = MAP_SAFE;
//   ModelData* data = otherType->object;
//   Vessel* vessel = &data->vessel;
// 
//   do {
//     /* minus epsilon sequence */
//     success = increment_phi_dof_by_delta(uType, vessel, -epsilon, size); CHECKERRQ(MAP_FATAL_61);        
//     if (data->MAP_SOLVE_TYPE==MONOLITHIC) {
//        success = line_solve_sequence(data, 0.0, map_msg, ierr);
//     } else {
//       success = node_solve_sequence(data, uType, zType, otherType, map_msg, ierr); // @todo CHECKERRQ()
//     };    
//     success = set_force_plus(yType->Fx, force->fx, size); CHECKERRQ(MAP_FATAL_61);
//     success = set_force_plus(yType->Fy, force->fy, size); CHECKERRQ(MAP_FATAL_61);
//     success = set_force_plus(yType->Fz, force->fz, size); CHECKERRQ(MAP_FATAL_61);
//     success = set_moment_plus_2(uType, yType, vessel, force->mx, force->my, force->mz, size); CHECKERRQ(MAP_FATAL_61);
//     success = restore_original_displacement(uType->x, originalX, size); CHECKERRQ(MAP_FATAL_61);        
//     success = restore_original_displacement(uType->y, originalY, size); CHECKERRQ(MAP_FATAL_61);        
//     success = restore_original_displacement(uType->z, originalZ, size); CHECKERRQ(MAP_FATAL_61);        
//     
//     /* plus epsilon sequence */
//     success = increment_phi_dof_by_delta(uType, vessel, epsilon, size); CHECKERRQ(MAP_FATAL_61);        
//     if (data->MAP_SOLVE_TYPE==MONOLITHIC) {
//       success = line_solve_sequence(data, 0.0, map_msg, ierr);
//     } else {
//       success = node_solve_sequence(data, uType, zType, otherType, map_msg, ierr); // @todo CHECKERRQ()
//     };    
//     success = set_force_minus(yType->Fx, force->fx, size); CHECKERRQ(MAP_FATAL_61);
//     success = set_force_minus(yType->Fy, force->fy, size); CHECKERRQ(MAP_FATAL_61);
//     success = set_force_minus(yType->Fz, force->fz, size); CHECKERRQ(MAP_FATAL_61);
//     success = set_moment_minus_2(uType, yType, vessel, force->mx, force->my, force->mz, size); CHECKERRQ(MAP_FATAL_61);
//     success = restore_original_displacement(uType->x, originalX, size); CHECKERRQ(MAP_FATAL_61);        
//     success = restore_original_displacement(uType->y, originalY, size); CHECKERRQ(MAP_FATAL_61);        
//     success = restore_original_displacement(uType->z, originalZ, size); CHECKERRQ(MAP_FATAL_61);                
//   } while (0);
//   
//   MAP_RETURN;
// };


// MAP_ERROR_CODE fd_the_sequence(MAP_OtherStateType_t* otherType, MAP_InputType_t* uType, MAP_OutputType_t* yType, MAP_ConstraintStateType_t* zType, Fd* force, const double epsilon, const int size, const double* originalX, const double* originalY, const double* originalZ, char* map_msg, MAP_ERROR_CODE* ierr)
// {
//   MAP_ERROR_CODE success = MAP_SAFE;
//   ModelData* data = otherType->object;
//   Vessel* vessel = &data->vessel;
// 
//   do {
//     /* minus epsilon sequence */
//     success = increment_the_dof_by_delta(uType, vessel, -epsilon, size); CHECKERRQ(MAP_FATAL_61);
//     if (data->MAP_SOLVE_TYPE==MONOLITHIC) {
//       success = line_solve_sequence(data, 0.0, map_msg, ierr);
//     } else {
//       success = node_solve_sequence(data, uType, zType, otherType, map_msg, ierr); // @todo CHECKERRQ()
//     };    
//     success = set_force_plus(yType->Fx, force->fx, size); CHECKERRQ(MAP_FATAL_61);
//     success = set_force_plus(yType->Fy, force->fy, size); CHECKERRQ(MAP_FATAL_61);
//     success = set_force_plus(yType->Fz, force->fz, size); CHECKERRQ(MAP_FATAL_61);
//     success = set_moment_plus_2(uType, yType, vessel, force->mx, force->my, force->mz, size); CHECKERRQ(MAP_FATAL_61);
//     success = restore_original_displacement(uType->x, originalX, size); CHECKERRQ(MAP_FATAL_61);        
//     success = restore_original_displacement(uType->y, originalY, size); CHECKERRQ(MAP_FATAL_61);        
//     success = restore_original_displacement(uType->z, originalZ, size); CHECKERRQ(MAP_FATAL_61);        
//     
//     /* plut epsilon sequence */
//     success = increment_the_dof_by_delta(uType, vessel, epsilon, size); CHECKERRQ(MAP_FATAL_61);        
//     if (data->MAP_SOLVE_TYPE==MONOLITHIC) {
//       success = line_solve_sequence(data, 0.0, map_msg, ierr);
//     } else {
//       success = node_solve_sequence(data, uType, zType, otherType, map_msg, ierr); // @todo CHECKERRQ()
//     };    
//     success = set_force_minus(yType->Fx, force->fx, size); CHECKERRQ(MAP_FATAL_61);
//     success = set_force_minus(yType->Fy, force->fy, size); CHECKERRQ(MAP_FATAL_61);
//     success = set_force_minus(yType->Fz, force->fz, size); CHECKERRQ(MAP_FATAL_61);
//     success = set_moment_minus_2(uType, yType, vessel, force->mx, force->my, force->mz, size); CHECKERRQ(MAP_FATAL_61);
//     success = restore_original_displacement(uType->x, originalX, size); CHECKERRQ(MAP_FATAL_61);        
//     success = restore_original_displacement(uType->y, originalY, size); CHECKERRQ(MAP_FATAL_61);        
//     success = restore_original_displacement(uType->z, originalZ, size); CHECKERRQ(MAP_FATAL_61);        
//   } while (0);
// 
//   MAP_RETURN;
// };


// MAP_ERROR_CODE fd_psi_sequence(MAP_OtherStateType_t* otherType, MAP_InputType_t* uType, MAP_OutputType_t* yType, MAP_ConstraintStateType_t* zType, Fd* force, const double epsilon, const int size, const double* originalX, const double* originalY, const double* originalZ, char* map_msg, MAP_ERROR_CODE* ierr)
// {
//   MAP_ERROR_CODE success = MAP_SAFE;
//   ModelData* data = otherType->object;
//   Vessel* vessel = &data->vessel;
//   do {
//     /* minus epsilon sequence */
//     success = increment_psi_dof_by_delta(uType, vessel, -epsilon, size); CHECKERRQ(MAP_FATAL_61);
//     if (data->MAP_SOLVE_TYPE==MONOLITHIC) {
//       success = line_solve_sequence(data, 0.0, map_msg, ierr);
//     } else {
//       success = node_solve_sequence(data, uType, zType, otherType, map_msg, ierr); // @todo CHECKERRQ()
//     };    
//     success = set_force_plus(yType->Fx, force->fx, size); CHECKERRQ(MAP_FATAL_61);
//     success = set_force_plus(yType->Fy, force->fy, size); CHECKERRQ(MAP_FATAL_61);
//     success = set_force_plus(yType->Fz, force->fz, size); CHECKERRQ(MAP_FATAL_61);
//     success = set_moment_plus_2(uType, yType, vessel, force->mx, force->my, force->mz, size); CHECKERRQ(MAP_FATAL_61);
//     success = restore_original_displacement(uType->x, originalX, size); CHECKERRQ(MAP_FATAL_61);
//     success = restore_original_displacement(uType->y, originalY, size); CHECKERRQ(MAP_FATAL_61);
//     success = restore_original_displacement(uType->z, originalZ, size); CHECKERRQ(MAP_FATAL_61);
//         
//     /* plut epsilon sequence */
//     success = increment_psi_dof_by_delta(uType, vessel, epsilon, size); CHECKERRQ(MAP_FATAL_61);
//     if (data->MAP_SOLVE_TYPE==MONOLITHIC) {
//       success = line_solve_sequence(data, 0.0, map_msg, ierr);
//     } else {
//       success = node_solve_sequence(data, uType, zType, otherType, map_msg, ierr); // @todo CHECKERRQ()
//     };    
//     success = set_force_minus(yType->Fx, force->fx, size); CHECKERRQ(MAP_FATAL_61);
//     success = set_force_minus(yType->Fy, force->fy, size); CHECKERRQ(MAP_FATAL_61);
//     success = set_force_minus(yType->Fz, force->fz, size); CHECKERRQ(MAP_FATAL_61);
//     success = set_moment_minus_2(uType, yType, vessel, force->mx, force->my, force->mz, size); CHECKERRQ(MAP_FATAL_61);
//     success = restore_original_displacement(uType->x, originalX, size); CHECKERRQ(MAP_FATAL_61);
//     success = restore_original_displacement(uType->y, originalY, size); CHECKERRQ(MAP_FATAL_61);
//     success = restore_original_displacement(uType->z, originalZ, size); CHECKERRQ(MAP_FATAL_61);
//   } while (0);
// 
//   MAP_RETURN;
// };



// MAP_ERROR_CODE calculate_stiffness_2(double* K, Fd* force, const double delta, const int size)
// {  
//   int i = 0;
//   for (i=0 ; i<size ; i++) {
//     K[0] += (force->fx[i]/(2*delta));
//     K[1] += (force->fy[i]/(2*delta));
//     K[2] += (force->fz[i]/(2*delta));
//     K[3] += (force->mx[i]/(2*delta));
//     K[4] += (force->my[i]/(2*delta));
//     K[5] += (force->mz[i]/(2*delta));
//   };
//   return MAP_SAFE;
// };


// MAP_ERROR_CODE set_moment_plus_2(const MAP_InputType_t* uType, const MAP_OutputType_t* outputType, const Vessel* vessel, double* mx, double* my, double* mz, const int size)
// {
//   int i = 0;
//   double rx = 0.0;
//   double ry = 0.0;
//   double rz = 0.0;
//   double m1 = 0.0;
//   double m2 = 0.0;
//   double m3 = 0.0;
//   double phi = vessel->orientation.phi.value;
// 
//   double delta = 1e-4;
//   for (i=0 ; i<size ; i++) {    
//     /* 
//        @todo forces are in global frame. Make sure everything is consistent. Should also add
//        reference origin. This need s to be converted to the global frame to be consistent with
//        uType. 
//     */
//     rx = uType->x[i];// (vessel->xi[i] - vessel->refOrigin.x.value);
//     ry = uType->y[i];// (vessel->yi[i] - vessel->refOrigin.y.value);
//     rz = uType->z[i];// (vessel->zi[i] - vessel->refOrigin.z.value);    
//     m1 = ((-outputType->Fz[i]*ry) - (-outputType->Fy[i]*rz));
//     m2 = ((-outputType->Fx[i]*rz) - (-outputType->Fz[i]*rx));
//     m3 = ((-outputType->Fy[i]*rx) - (-outputType->Fx[i]*ry));
//     mx[i] += m1;
//     my[i] += m2*cos(phi+delta) - m3*sin(phi+delta);
//     mz[i] += m2*sin(phi+delta) + m3*cos(phi+delta);    
//     // mx[i] += ((-outputType->Fz[i]*ry) - (-outputType->Fy[i]*rz));
//     // my[i] += ((-outputType->Fx[i]*rz) - (-outputType->Fz[i]*rx));
//     // mz[i] += ((-outputType->Fy[i]*rx) - (-outputType->Fx[i]*ry));
//   };
//   return MAP_SAFE;
// };


// MAP_ERROR_CODE set_moment_minus_2(const MAP_InputType_t* uType, const MAP_OutputType_t* outputType, const Vessel* vessel, double* mx, double* my, double* mz, const int size)
// {
//   int i = 0;
//   double rx = 0.0;
//   double ry = 0.0;
//   double rz = 0.0;
//   double m1 = 0.0;
//   double m2 = 0.0;
//   double m3 = 0.0;
//   double phi = vessel->orientation.phi.value;
// 
//   double delta = 1e-4;
//   for (i=0 ; i<size ; i++) {    
//     rx = uType->x[i];// (vessel->xi[i] - vessel->refOrigin.x.value);
//     ry = uType->y[i];// (vessel->yi[i] - vessel->refOrigin.y.value);
//     rz = uType->z[i];// (vessel->zi[i] - vessel->refOrigin.z.value);    
//     m1 = ((-outputType->Fz[i]*ry) - (-outputType->Fy[i]*rz));
//     m2 = ((-outputType->Fx[i]*rz) - (-outputType->Fz[i]*rx));
//     m3 = ((-outputType->Fy[i]*rx) - (-outputType->Fx[i]*ry));
//     mx[i] -= m1;
//     my[i] -= m2*cos(phi+delta) - m3*sin(phi+delta);
//     mz[i] -= m2*sin(phi+delta) + m3*cos(phi+delta);    
//     // mx[i] -= ((-outputType->Fz[i]*ry) - (-outputType->Fy[i]*rz));
//     // my[i] -= ((-outputType->Fx[i]*rz) - (-outputType->Fz[i]*rx));
//     // mz[i] -= ((-outputType->Fy[i]*rx) - (-outputType->Fx[i]*ry));
//   };
//   return MAP_SAFE;
// };


// /**
//  * success = calculate_stiffness(stiffness, fx, epsilon, N);
//  */
// MAP_ERROR_CODE calculate_stiffness(double* K, double* force, const double delta, const int size)
// {  
//   *K = 0.0;
//   int i = 0;
//   for (i=0 ; i<size ; i++) {
//     *K += (force[i]/(2*delta));
//   };
//   return MAP_SAFE;
// };
