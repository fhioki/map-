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


#include "lineroutines.h"


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
  // Node* node_iter = NULL;
  // Element* element_iter = NULL;
  // MapReal psi = 0.0;
  // MapReal fx = 0.0;
  // MapReal fy = 0.0;
  // MapReal fz = 0.0;
  // MapReal fx_a = 0.0;
  // MapReal fy_a = 0.0;
  // MapReal fz_a = 0.0;
  // 
  // list_iterator_start(&data->element);         /* starting an iteration "session" */
  // while (list_iterator_hasnext(&data->element)) { /* tell whether more values available */
  //   element_iter = (Element*)list_iterator_next(&data->element);
  //   psi = element_iter->psi.value;
  //   fx = *(element_iter->H.value)*cos(psi);
  //   fx_a = -(element_iter->HAtAnchor.value)*cos(psi);
  // 
  //   fy = *(element_iter->H.value)*sin(psi);
  //   fy_a = -(element_iter->HAtAnchor.value)*sin(psi);
  //   
  //   fz = *(element_iter->V.value);
  //   fz_a = -(element_iter->VAtAnchor.value);
  // 
  //   add_to_sum_fx(element_iter->fairlead, fx);
  //   add_to_sum_fx(element_iter->anchor, fx_a);
  // 
  //   add_to_sum_fy(element_iter->fairlead, fy);
  //   add_to_sum_fy(element_iter->anchor, fy_a);
  // 
  //   add_to_sum_fz(element_iter->fairlead, fz);
  //   add_to_sum_fz(element_iter->anchor, fz_a);
  // }; 
  // list_iterator_stop(&model_data->element); /* ending the iteration "session" */    
  // 
  // list_iterator_start(&model_data->node);            /* starting an iteration "session" */
  // while (list_iterator_hasnext(&model_data->node)) { /* tell whether more values available */
  //   node_iter = (Node*)list_iterator_next(&model_data->node);
  //   if (node_iter->type==CONNECT) {
  //     add_to_sum_fx(node_iter, -(node_iter->externalForce.fx.value));
  //     add_to_sum_fy(node_iter, -(node_iter->externalForce.fy.value));
  //     add_to_sum_fz(node_iter, -(node_iter->externalForce.fz.value));
  //   };
  // }; 
  // list_iterator_stop(&model_data->node); /* ending the iteration "session" */    
  return MAP_SAFE;
}


/**
 * increment sum force value by (f) if node is fairlead; (-f) is node is anchor. 
 */
// void add_to_sum_fx(Node* node, const MapReal fx)
// {
//   *(node->sumForcePtr.fx.value) += fx;
// };
// 
// 
// /**
//  * increment sum force value by (f) if node is fairlead; (-f) is node is anchor. 
//  */
// void add_to_sum_fy(Node* node, const MapReal fy)
// {
//   *(node->sumForcePtr.fy.value) += fy;
// };
// 
// 
// /**
//  * increment sum force value by (f) if node is fairlead; (-f) is node is anchor. 
//  */
// void add_to_sum_fz(Node* node, const MapReal fz)
// {
//   *(node->sumForcePtr.fz.value) += fz;
// };
// 
// 
// 
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
// 
// 
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
    checkpoint();
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
    
    printf("Error: %f, tol=%f\n",sqrt(error),ns->tol);
    
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
      user_msg = bformat("Line segment %s, l = %d [m].", n, element_iter->l.value);
      *ierr = map_set_universal_error(user_msg, map_msg, *ierr, MAP_FATAL_54);        
      success = bdestroy(user_msg); 
      user_msg = NULL;
      break; 
    } else if (element_iter->h.value<=-MACHINE_EPSILON) {
      user_msg = bformat("Line segment %s, h = %d [m].", n, element_iter->h.value);
      *ierr = map_set_universal_error(user_msg, map_msg, *ierr, MAP_FATAL_55);        
      success = bdestroy(user_msg); 
      user_msg = NULL;
      break; 
    } else if (element_iter->lineProperty->omega>0.0) {
      success = check_maximum_line_length(element_iter, element_iter->options.omitContact, map_msg, ierr);
      if (success) {        
        user_msg = bformat("Line segment %s.", n);
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
    *ierr = map_set_universal_error(user_msg, map_msg, *ierr, MAP_FATAL_59);        
    success = bdestroy(user_msg); 
    user_msg = NULL;
    return MAP_FATAL;
  };
  return MAP_SAFE;
};
