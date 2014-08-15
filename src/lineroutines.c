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
#include "lineroutines.h"


// MAP_ERROR_CODE reset_node_force_to_zero(ModelData* data, char* map_msg, MAP_ERROR_CODE* ierr)
// {
//   Node* iterNode = NULL;
// 
//   list_iterator_start(&data->node);            /* starting an iteration "session" */
//   while (list_iterator_hasnext(&data->node)) { /* tell whether more values available */
//     iterNode = (Node*)list_iterator_next(&data->node);
//     *(iterNode->sumForcePtr.fx.value) = 0.0;
//     *(iterNode->sumForcePtr.fy.value) = 0.0;
//     *(iterNode->sumForcePtr.fz.value) = 0.0;    
//   };
//   list_iterator_stop(&data->node); /* ending the iteration "session" */    
//   return MAP_SAFE;
// };
// 
// 
// /**
//  * call immediately after set_line_variables_post_solve(); this added H and V
//  */
// MAP_ERROR_CODE calculate_node_sum_force(ModelData* data)
// {
//   Node* iterNode = NULL;
//   Element* iterElement = NULL;
//   MapReal psi = 0.0;
//   MapReal fx = 0.0;
//   MapReal fy = 0.0;
//   MapReal fz = 0.0;
//   MapReal fx_a = 0.0;
//   MapReal fy_a = 0.0;
//   MapReal fz_a = 0.0;
//   
//   list_iterator_start(&data->element);         /* starting an iteration "session" */
//   while (list_iterator_hasnext(&data->element)) { /* tell whether more values available */
//     iterElement = (Element*)list_iterator_next(&data->element);
//     psi = iterElement->psi.value;
//     fx = *(iterElement->H.value)*cos(psi);
//     fx_a = -(iterElement->HAtAnchor.value)*cos(psi);
// 
//     fy = *(iterElement->H.value)*sin(psi);
//     fy_a = -(iterElement->HAtAnchor.value)*sin(psi);
//     
//     fz = *(iterElement->V.value);
//     fz_a = -(iterElement->VAtAnchor.value);
// 
//     add_to_sum_fx(iterElement->fairlead, fx);
//     add_to_sum_fx(iterElement->anchor, fx_a);
// 
//     add_to_sum_fy(iterElement->fairlead, fy);
//     add_to_sum_fy(iterElement->anchor, fy_a);
// 
//     add_to_sum_fz(iterElement->fairlead, fz);
//     add_to_sum_fz(iterElement->anchor, fz_a);
//   }; 
//   list_iterator_stop(&data->element); /* ending the iteration "session" */    
// 
//   list_iterator_start(&data->node);            /* starting an iteration "session" */
//   while (list_iterator_hasnext(&data->node)) { /* tell whether more values available */
//     iterNode = (Node*)list_iterator_next(&data->node);
//     if (iterNode->type==CONNECT) {
//       add_to_sum_fx(iterNode, -(iterNode->externalForce.fx.value));
//       add_to_sum_fy(iterNode, -(iterNode->externalForce.fy.value));
//       add_to_sum_fz(iterNode, -(iterNode->externalForce.fz.value));
//     };
//   }; 
//   list_iterator_stop(&data->node); /* ending the iteration "session" */    
//   return MAP_SAFE;
// }
// 
// 
// /**
//  * increment sum force value by (f) if node is fairlead; (-f) is node is anchor. 
//  */
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


// /**
//  * set: alpha, alpha at anchor, H, V and T at anchor, T at fairlead
//  */
// MAP_ERROR_CODE set_line_variables_post_solve(ModelData* data, char* map_msg, MAP_ERROR_CODE* ierr)
// {
//   Element* iterElement = NULL;
//   MapReal H = 0.0;
//   MapReal V = 0.0;
//   MapReal Ha = 0.0;
//   MapReal Va = 0.0;
//   MapReal w = 0.0;
//   MapReal Lu = 0.0;
//   MapReal Lb = 0.0;
//   MapReal cb = 0.0;
//   bool contactFlag = false;
// 
//   list_iterator_start(&data->element); /* starting an iteration "session" */
//   while (list_iterator_hasnext(&data->element)) { /* tell whether more values available */ 
//     iterElement = (Element*)list_iterator_next(&data->element);
//     
//     /* altitude angle at fairlead */
//     H = *(iterElement->H.value);
//     V = *(iterElement->V.value);
//     iterElement->alpha.value = atan2(V, H);
//     
//     /* tension at fairlead */
//     iterElement->T.value = sqrt(H*H + V*V);
//     
//     /* line in contact with seabed options */
//     w = iterElement->lineProperty->omega;
//     cb = iterElement->lineProperty->cb;
//     contactFlag = iterElement->options.omitContact;
//     Lu = iterElement->Lu.value;
//     if (contactFlag==true || w<0.0 || (V-w*Lu)>0.0) { 
//       iterElement->lb.value = 0.0;
//       Ha = H;
//       Va = V - w*Lu;
//     } else { /* line is touching the seabed */
//       Lb = Lu - (V/w);
//       iterElement->lb.value = Lb;
//       ((H-cb*w*Lb)>0.0) ? (Ha = (H-cb*w*Lb)) : (Ha = 0.0);
//       Va = 0.0;
//     };
//     
//     /* tension at anchor */
//     iterElement->HAtAnchor.value = Ha;
//     iterElement->VAtAnchor.value = Va;
//     iterElement->TAtAnchor.value = sqrt(Ha*Ha + Va*Va);
//     iterElement->alphaAtAnchor.value = atan2(Va, Ha);    
//   };
//   list_iterator_stop(&data->element); /* ending the iteration "session" */        
//   return MAP_SAFE;
// };


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


// /**
//  *
//  */
// MAP_ERROR_CODE set_element_initial_guess(ModelData* data, char* map_msg, MAP_ERROR_CODE* ierr)
// {
//   Element* iterElement = NULL;
//   MapReal lambda = 0.0;
//   MapReal length = 0.0;
//   MapReal height = 0.0;
//   MapReal w = 0.0;
//   MapReal Lu = 0.0;
// 
//   list_iterator_start(&data->element);            /* starting an iteration "session" */
//   while (list_iterator_hasnext(&data->element)) { /* tell whether more values available */ 
//     iterElement = (Element*)list_iterator_next(&data->element);    
//     w = iterElement->lineProperty->omega;
//     length = iterElement->l.value;
//     height = iterElement->h.value;
//     Lu = iterElement->Lu.value;
// 
//     /* 
//        @todo: this is legacy from MAP v > 1.0 for the case of a perfectly vertical cable        
//        if ( norm <= TOL ) {
//          psi = 1.0;//checkpoint();
//        } else {
//          if ( (fairlead->Y.value - anchor->Y.value) >= 0) {
//            // this simply finds the angle psi between the local and global reference frames simply by evaluation trig relationships
//            psi = acos( (fairlead->X.value - anchor->X.value)/norm );
//          } else {
//            psi = -acos( (fairlead->X.value - anchor->X.value)/norm );
//          }
//        }
//     */
// 
//     /* note: the element horizontal (l) and vertical (h) excursion are previously initialized in set_line_variables_pre_solve(...) */ 
//     if (length<=1e-2) {
//       lambda = 1000000;
//     } else if (sqrt(length*length+height*height)>=Lu) {
//       lambda = 0.2;
//     } else {
//       lambda = sqrt(3*((Lu*Lu - height*height)/(length*length) - 1));
//     };
// 
//     /* note: omega is previously initialized in function initialize_cable_library_variables(.. ) in mapinit.c */    
//     *(iterElement->H.value) = fabs(w*length/(2*lambda));
//     *(iterElement->V.value) = (w/2)*((height/tanh(lambda)) + Lu);
//   };
//   list_iterator_stop(&data->element); /* ending the iteration "session" */    
//   return MAP_SAFE;
// };
// 
// 
// 
// /**
//  *
//  */
// MAP_ERROR_CODE line_solve_sequence(ModelData* data, double t, char* map_msg, MAP_ERROR_CODE* ierr) 
// {
//   MAP_ERROR_CODE success = MAP_SAFE;
// 
//   do { 
//     success = set_line_variables_pre_solve(data, map_msg, ierr);
//     success = reset_node_force_to_zero(data, map_msg, ierr);    
//     success = solve_line(data, t, map_msg, ierr); 
//     // /* prematurely terminating the the line solve routine to return back to the
//     //    caller function. Do a back tracking on the solution and attempt recovery */
//     // if (success==MAP_FATAL) {
//     //   return MAP_FATAL_59;
//     // };   
//     // CHECKERRQ(MAP_FATAL_39); /* @todo: this should be called elsewhere to notifiy 
//     //                             users ofpremature termination. This won't be caught
//     //                             when solve_line fails. */    
//     success = set_line_variables_post_solve(data, map_msg, ierr);
//     success = calculate_node_sum_force(data);
//   } while (0);
// 
//   MAP_RETURN;
// };    
