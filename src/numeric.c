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


#include "numeric.h"
#include "lineroutines.h"


int inner_function_evals(void* line_ptr, int m, int n, const __cminpack_real__* x, __cminpack_real__* fvec, __cminpack_real__* fjac, int ldfjac, int iflag) 
{
  Line* line = (Line*)line_ptr;
  double Fh = x[0]; /* Fh is not const because it is artificially set to epsilon when the line is perfectly vertical */
  const double Fv = x[1];  
  const double EA = line->line_property->EA;
  const double Lu = line->Lu.value;
  const double height = line->h.value;
  const double length = line->l.value;
  const double omega = line->line_property->omega;
  const double cb = line->line_property->cb;
  const bool contactFlag = line->options.omit_contact;
  
  if (iflag==0) {
    return 0;
  };
 
  if (Fh<MAP_HORIZONTAL_TOL) { /* perfectly vertical case */
    Fh = MAP_HORIZONTAL_TOL;
  };

  if (iflag!=2) {
    if (contactFlag==true || omega<0.0 || (Fv-omega*Lu)>0.0) { /* true when no portion of the line rests on the seabed */
      fvec[0] = residual_function_length_no_contact(Fv, Fh, omega, Lu, EA, length);
      fvec[1] = residual_function_height_no_contact(Fv, Fh, omega, Lu, EA, height); 
    } else { /* true when a portion of the line rests on the seabed and the anchor tension is nonzero */
      fvec[0] = residual_function_length_contact(Fv, Fh, omega, Lu, EA, length, cb);
      fvec[1] = residual_function_height_contact(Fv, Fh, omega, Lu, EA, height, cb); 
    };
  } else {
    if (contactFlag==true || omega<0.0 || (Fv-omega*Lu)>0.0) { /* true when no portion of the line rests on the seabed */
      fjac[0] = jacobian_dxdh_no_contact(Fv, Fh, omega, Lu, EA);
      fjac[1] = jacobian_dxdv_no_contact(Fv, Fh, omega, Lu, EA);
      fjac[2] = jacobian_dzdh_no_contact(Fv, Fh, omega, Lu, EA);
      fjac[3] = jacobian_dzdv_no_contact(Fv, Fh, omega, Lu, EA);
    } else { /* true when a portion of the line rests on the seabed and the anchor tension is nonzero */
      fjac[0] = jacobian_dxdh_contact(Fv, Fh, omega, Lu, EA, cb);
      fjac[1] = jacobian_dxdv_contact(Fv, Fh, omega, Lu, EA, cb);
      fjac[2] = jacobian_dzdh_contact(Fv, Fh, omega, Lu, EA, cb);
      fjac[3] = jacobian_dzdv_contact(Fv, Fh, omega, Lu, EA, cb);
    };
  };
  return 0;
};


MapReal residual_function_length_no_contact(const MapReal V, const MapReal H, const MapReal w, const MapReal Lu, const MapReal EA, const MapReal l)
{  
  return (H/w)*asinh(V/H) - (H/w)*asinh( (V-w*Lu)/H ) + ((H*Lu)/(EA)) - l;
};


MapReal residual_function_height_no_contact(const MapReal V, const MapReal H, const MapReal w, const MapReal Lu, const MapReal EA, const MapReal h)
{
  return (H/w)* sqrt(1 + pow((V/H), 2)) - (H/w)*sqrt(1 + pow(((V-w*Lu)/H), 2)) + 1/(EA)*(V*Lu - (w*Lu*Lu)/2) - h;
};


MapReal jacobian_dxdh_no_contact(const MapReal V, const MapReal H, const MapReal w, const MapReal Lu, const MapReal EA)
{
  return (asinh(V/H) - asinh((V-w*Lu)/H))/w - ((V/H + pow(V/H, 2)/sqrt(1.0 + pow(V/H, 2)))/(V/H + sqrt(1.0 + pow(V/H, 2))) 
                                               - ((V-w*Lu)/H + pow((V-w*Lu)/H, 2)/sqrt(1.0 + pow((V-w*Lu)/H, 2)))
                                               /((V-w*Lu)/H + sqrt(1.0 + pow((V-w*Lu)/H, 2))))/w + (Lu/(EA));      
};


MapReal jacobian_dxdv_no_contact(const MapReal V, const MapReal H, const MapReal w, const MapReal Lu, const MapReal EA)
{
  return ((1.0 + V/H /sqrt(1.0 + pow(V/H, 2)))/(V/H + sqrt(1.0 + pow(V/H, 2))) 
          - (1.0 + (V-w*Lu)/H /sqrt(1.0 + pow( (V-w*Lu)/H , 2)))
          /((V-w*Lu)/H + sqrt(1.0 + pow((V-w*Lu)/H, 2))))/w;
};


MapReal jacobian_dzdh_no_contact(const MapReal V, const MapReal H, const MapReal w, const MapReal Lu, const MapReal EA)
{
  return ( sqrt( 1.0 + pow( V/H , 2) ) - sqrt( 1.0 + pow( (V-w*Lu)/H , 2) ) )/w  
    - ( pow( V/H , 2 )/sqrt( 1.0 + pow( V/H , 2) ) - pow( (V-w*Lu)/H , 2)/sqrt( 1.0 + pow( (V-w*Lu)/H , 2) ) )/w;    
};


MapReal jacobian_dzdv_no_contact(const MapReal V, const MapReal H, const MapReal w, const MapReal Lu, const MapReal EA)
{
  return ( V/H/sqrt( 1.0 + pow( V/H , 2) ) - (V-w*Lu)/H /sqrt( 1.0 + pow( (V-w*Lu)/H , 2) ) )/w + (Lu/(EA));
};


MapReal residual_function_length_contact(const MapReal V, const MapReal H, const MapReal w, const MapReal Lu, const MapReal EA, const MapReal l, const MapReal cb)
{  
  /* Note that Lb = Lu - V/w */
  if (-cb*(V-w*Lu)<H) { /* true when a portion of the line rests on the seabed and the anchor tension is nonzero */
    return log((V/H) + sqrt(1.0 + pow(V/H,2)))*(H/w) - 0.5*(cb/EA)*w*(Lu-V/w)*(Lu-V/w) + (Lu/EA)*H + (Lu-V/w) - l;
  } else { /* 0.0<H<=-CB*(V-w*Lu), A  portion of the line must rest on the seabed and the anchor tension is zero */
    return log((V/H) + sqrt(1.0 + pow(V/H,2)))*(H/w) - 0.5*(cb/EA)*w*((Lu-V/w)*(Lu-V/w) - ((Lu-V/w) - (H/w)/cb)*((Lu-V/w) - (H/w)/cb)) + (Lu/EA)*H + (Lu-V/w) - l;
  };
};


MapReal residual_function_height_contact(const MapReal V, const MapReal H, const MapReal w, const MapReal Lu, const MapReal EA, const MapReal h, const MapReal cb)
{
  // @todo remove this conditional statement because the equations are the same regardless of the outcome
  if (-cb*(V-w*Lu)<H) { /* true when a portion of the line rests on the seabed and the anchor tension is nonzero */
    return (sqrt(1.0 + pow(V/H,2)) - 1.0)*(H/w) + 0.5*V*(V/(w*EA)) - h;
  } else { /* 0.0<H<=-CB*(V-w*Lu), A  portion of the line must rest on the seabed and the anchor tension is zero */
    return (sqrt(1.0 + pow(V/H,2)) - 1.0)*(H/w) + 0.5*V*(V/(w*EA)) - h;
  };
};


MapReal jacobian_dxdh_contact(const MapReal V, const MapReal H, const MapReal w, const MapReal Lu, const MapReal EA, const MapReal cb)
{
  if (-cb*(V-w*Lu)<H) { /* true when a portion of the line rests on the seabed and the anchor tension is nonzero */
    return log((V/H) + sqrt(1.0 + pow(V/H,2)))/w - (((V/H) + (V/H)*(V/H)/sqrt(1.0 + pow(V/H,2)))/((V/H) + sqrt(1.0 + pow(V/H,2))))/w + (Lu/EA);
  } else { /* 0.0<H<=-CB*(V-w*Lu), A  portion of the line must rest on the seabed and the anchor tension is zero */
    return log((V/H) + sqrt(1.0 + pow(V/H,2)))/w - (((V/H) + (V/H)*(V/H)/sqrt(1.0 + pow(V/H,2)))/((V/H) + sqrt(1.0 + pow(V/H,2))))/w + (Lu/EA) - ((Lu-V/w) - (H/w)/cb)/EA;
  };
};


MapReal jacobian_dxdv_contact(const MapReal V, const MapReal H, const MapReal w, const MapReal Lu, const MapReal EA, const MapReal cb)
{
  if (-cb*(V-w*Lu)<H) { /* true when a portion of the line rests on the seabed and the anchor tension is nonzero */
    return ((1.0 + (V/H)/sqrt(1.0 + pow(V/H,2)))/((V/H) + sqrt(1.0 + pow(V/H,2))))/w + (cb/EA)*(Lu-V/w) - 1.0/w;
  } else { /* 0.0<H<=-CB*(V-w*Lu), A  portion of the line must rest on the seabed and the anchor tension is zero */
    return ((1.0 + (V/H)/sqrt(1.0 + pow(V/H,2)))/((V/H) + sqrt(1.0 + pow(V/H,2))))/w + (H/(w*EA)) - 1.0/w;
  };
};


MapReal jacobian_dzdh_contact(const MapReal V, const MapReal H, const MapReal w, const MapReal Lu, const MapReal EA, const MapReal cb)
{
  // @todo remove this conditional statement because the equations are the same regardless of the outcome
  if (-cb*(V-w*Lu)<H) { /* true when a portion of the line rests on the seabed and the anchor tension is nonzero */
    return (sqrt(1.0 + pow(V/H,2)) - 1.0 - (V/H)*(V/H)/sqrt(1.0 + pow(V/H,2)))/w;
  } else { /* 0.0<H<=-CB*(V-w*Lu), A  portion of the line must rest on the seabed and the anchor tension is zero */
    return (sqrt(1.0 + pow(V/H,2)) - 1.0 - (V/H)*(V/H)/sqrt(1.0 + pow(V/H,2)))/w;
  };
};



MapReal jacobian_dzdv_contact(const MapReal V, const MapReal H, const MapReal w, const MapReal Lu, const MapReal EA, const MapReal cb)
{
  // @todo remove this conditional statement because the equations are the same regardless of the outcome
  if (-cb*(V-w*Lu)<H) { /* true when a portion of the line rests on the seabed and the anchor tension is nonzero */
    return ((V/H)/sqrt(1.0 + pow(V/H,2)))/w + (V/(w*EA));
  } else { /* 0.0<H<=-CB*(V-w*Lu), A  portion of the line must rest on the seabed and the anchor tension is zero */
    return ((V/H)/sqrt(1.0 + pow(V/H,2)))/w + (V/(w*EA));
  };  
};




double get_maximum_line_length(Line* line)
{
  const double l = line->l.value;
  const double h = line->h.value;
  const double EA = line->line_property->EA;
  const double w = line->line_property->omega;
  const double Lu = line->Lu.value;

  return (l - EA/w + sqrt(pow((EA/w),2) + 2.0*h*EA/w));
};


MAP_ERROR_CODE lu_back_substitution(OuterSolveAttributes* ns, const int n, char* map_msg, MAP_ERROR_CODE* ierr)
{
  int i = 0;
  int j = 0;

  /* Forward solve Ly = b */
  for (i=0 ; i<n ; i++) {
    ns->y[i] = ns->b[i];
    for (j=0 ; j<i ; j++) {
      ns->y[i] -= (ns->l[i][j])*(ns->y[j]);
    };
    if (fabs(ns->l[i][i])<MACHINE_EPSILON) {
      return MAP_FATAL;
    };
    ns->y[i] /= ns->l[i][i];    
  };

  /* Backward solve Ux = y */
  for (i=n-1 ; i>=0 ; i--) {
    ns->x[i] = ns->y[i];
    for (j=i+1 ; j<n ; j++) {
      ns->x[i] -= (ns->u[i][j])*(ns->x[j]);
    };    
    if (fabs(ns->u[i][i])<MACHINE_EPSILON) {
      return MAP_FATAL;
    };
    ns->x[i] /= ns->u[i][i];
  };
  return MAP_SAFE;
};



MAP_ERROR_CODE lu(OuterSolveAttributes* ns, const int n, char* map_msg, MAP_ERROR_CODE* ierr)
{
  int i = 0;
  int j = 0;
  int k = 0;

  for (i=0 ; i<n ; i++) {
    ns->l[i][i] = 1.0;
    for (j=i+1 ; j<n ; j++) {
      if (fabs(ns->jac[i][i])<MACHINE_EPSILON) {
        return MAP_FATAL;
      };
      ns->l[j][i] = (ns->jac[j][i])/(ns->jac[i][i]);
      /* ns->jac[j][j] = ns->l[j][j] */
      for (k=i+1 ; k<n ; k++) {
        ns->jac[j][k] = ns->jac[j][k] - (ns->l[j][i])*(ns->jac[i][k]);
      };
    };
    
    for (k=i ; k<n ; k++) {
      ns->u[i][k] = ns->jac[i][k] ;
    };
  };

  return MAP_SAFE;
};



MAP_ERROR_CODE forward_difference_jacobian(MAP_OtherStateType_t* other_type, MAP_ParameterType_t* p_type, MAP_ConstraintStateType_t* z_type, ModelData* model_data, char* map_msg, MAP_ERROR_CODE* ierr)
{
  OuterSolveAttributes* ns = &model_data->outer_loop;
  MAP_ERROR_CODE success = MAP_SAFE;
  double original_displacement = 0.0;
  const int THREE = 3;
  const int z_size = z_type->z_Len; 
  const int m = THREE*(other_type->Fz_connect_Len); // rows
  const int n = THREE*(z_type->z_Len);              // columns
  int i = 0;
  int j = 0;
  
  for (i=0 ; i<z_size ; i++) {
    ns->b[THREE*i] = other_type->Fx_connect[i];
    ns->b[THREE*i+1] = other_type->Fy_connect[i];
    ns->b[THREE*i+2] = other_type->Fz_connect[i];      
  }

  /* First store the connect node force before applying ns->epsilon displacement */
  for (i=0 ; i<n ; i++) {            
    for (j=0 ; j<z_size ; j++) {            
      ns->jac[THREE*j][i] = -other_type->Fx_connect[j];
      ns->jac[THREE*j+1][i] = -other_type->Fy_connect[j];
      ns->jac[THREE*j+2][i] = -other_type->Fz_connect[j];
    };
  };
    
  for (j=0 ; j<z_size ; j++) {   // column               
    for (i=0 ; i<z_size ; i++) { // rows           
      original_displacement = z_type->x[j];
      z_type->x[j] += ns->epsilon;
      success = line_solve_sequence(model_data, p_type, 0.0, map_msg, ierr);
      if (success) {
        set_universal_error_with_message(map_msg, ierr, MAP_FATAL_78, "Forward difference, x[%d]+delta, row %d, col %d.", j+1, THREE*i, THREE*j);
        return MAP_FATAL;
      };
      ns->jac[THREE*i][THREE*j] += other_type->Fx_connect[i];
      ns->jac[THREE*i][THREE*j] /= ns->epsilon;
      ns->jac[THREE*i+1][THREE*j] += other_type->Fy_connect[i];
      ns->jac[THREE*i+1][THREE*j] /= ns->epsilon;
      ns->jac[THREE*i+2][THREE*j] += other_type->Fz_connect[i];
      ns->jac[THREE*i+2][THREE*j] /= ns->epsilon;
      z_type->x[j] = original_displacement;
        
      original_displacement = z_type->y[j];
      z_type->y[j] += ns->epsilon;
      success = line_solve_sequence(model_data, p_type, 0.0, map_msg, ierr);
      if (success) {
        set_universal_error_with_message(map_msg, ierr, MAP_FATAL_78, "Forward difference, x[%d]+delta, row %d, col %d.", j+1, THREE*i, THREE*j);
        return MAP_FATAL;
      };
      ns->jac[THREE*i][THREE*j+1] += other_type->Fx_connect[i];
      ns->jac[THREE*i][THREE*j+1] /= ns->epsilon;
      ns->jac[THREE*i+1][THREE*j+1] += other_type->Fy_connect[i];
      ns->jac[THREE*i+1][THREE*j+1] /= ns->epsilon;
      ns->jac[THREE*i+2][THREE*j+1] += other_type->Fz_connect[i];
      ns->jac[THREE*i+2][THREE*j+1] /= ns->epsilon;
      z_type->y[j] = original_displacement;
    
      original_displacement = z_type->z[j];
      z_type->z[j] += ns->epsilon;
      success = line_solve_sequence(model_data, p_type, 0.0, map_msg, ierr);
      if (success) {
        set_universal_error_with_message(map_msg, ierr, MAP_FATAL_78, "Forward difference, x[%d]+delta, row %d, col %d.", j+1, THREE*i, THREE*j);
        return MAP_FATAL;
      };
      ns->jac[THREE*i][THREE*j+2] += other_type->Fx_connect[i];
      ns->jac[THREE*i][THREE*j+2] /= ns->epsilon;
      ns->jac[THREE*i+1][THREE*j+2] += other_type->Fy_connect[i];
      ns->jac[THREE*i+1][THREE*j+2] /= ns->epsilon;
      ns->jac[THREE*i+2][THREE*j+2] += other_type->Fz_connect[i];
      ns->jac[THREE*i+2][THREE*j+2] /= ns->epsilon;
      z_type->z[j] = original_displacement;
    };
  };

  /* read flag to set scaling parameter */
  if (ns->pg) {
    for (i=0 ; i<THREE*z_size ; i++) { 
      ns->jac[i][i] += (ns->ds/pow(ns->iteration_count,1.5)+ns->d);
    };
  };

  return MAP_SAFE;
};


MAP_ERROR_CODE backward_difference_jacobian(MAP_OtherStateType_t* other_type, MAP_ParameterType_t* p_type, MAP_ConstraintStateType_t* z_type, ModelData* model_data, char* map_msg, MAP_ERROR_CODE* ierr)
{
  OuterSolveAttributes* ns = &model_data->outer_loop;
  MAP_ERROR_CODE success = MAP_SAFE;
  double original_displacement = 0.0;
  const int THREE = 3;
  const int z_size = z_type->z_Len; // N
  const int m = THREE*(other_type->Fz_connect_Len); // rows
  const int n = THREE*(z_type->z_Len);              // columns
  int i = 0;
  int j = 0;
  
  for (i=0 ; i<z_size ; i++) {
    ns->b[THREE*i] = other_type->Fx_connect[i];
    ns->b[THREE*i+1] = other_type->Fy_connect[i];
    ns->b[THREE*i+2] = other_type->Fz_connect[i];      
  }

  /* First store the connect node force before applying ns->epsilon displacement */
  for (i=0 ; i<n ; i++) {            
    for (j=0 ; j<z_size ; j++) {            
      ns->jac[THREE*j][i] = other_type->Fx_connect[j];
      ns->jac[THREE*j+1][i] = other_type->Fy_connect[j];
      ns->jac[THREE*j+2][i] = other_type->Fz_connect[j];
    };
  };

  for (j=0 ; j<z_size ; j++) {   // column               
    for (i=0 ; i<z_size ; i++) { // rows           
      original_displacement = z_type->x[j];
      z_type->x[j] -= ns->epsilon;
      success = line_solve_sequence(model_data, p_type, 0.0, map_msg, ierr);
      if (success) {
        set_universal_error_with_message(map_msg, ierr, MAP_FATAL_78, "Backward difference, x[%d]+delta, row %d, col %d.", j+1, THREE*i, THREE*j);
        return MAP_FATAL;
      };
      ns->jac[THREE*i][THREE*j] -= other_type->Fx_connect[i];
      ns->jac[THREE*i][THREE*j] /= ns->epsilon;
      ns->jac[THREE*i+1][THREE*j] -= other_type->Fy_connect[i];
      ns->jac[THREE*i+1][THREE*j] /= ns->epsilon;
      ns->jac[THREE*i+2][THREE*j] -= other_type->Fz_connect[i];
      ns->jac[THREE*i+2][THREE*j] /= ns->epsilon;
      z_type->x[j] = original_displacement;
        
      original_displacement = z_type->y[j];
      z_type->y[j] -= ns->epsilon;
      success = line_solve_sequence(model_data, p_type, 0.0, map_msg, ierr);
      if (success) {
        set_universal_error_with_message(map_msg, ierr, MAP_FATAL_78, "Backward difference, x[%d]+delta, row %d, col %d.", j+1, THREE*i, THREE*j);
        return MAP_FATAL;
      };
      ns->jac[THREE*i][THREE*j+1] -= other_type->Fx_connect[i];
      ns->jac[THREE*i][THREE*j+1] /= ns->epsilon;
      ns->jac[THREE*i+1][THREE*j+1] -= other_type->Fy_connect[i];
      ns->jac[THREE*i+1][THREE*j+1] /= ns->epsilon;
      ns->jac[THREE*i+2][THREE*j+1] -= other_type->Fz_connect[i];
      ns->jac[THREE*i+2][THREE*j+1] /= ns->epsilon;
      z_type->y[j] = original_displacement;
    
      original_displacement = z_type->z[j];
      z_type->z[j] -= ns->epsilon;
      success = line_solve_sequence(model_data, p_type, 0.0, map_msg, ierr);
      if (success) {
        set_universal_error_with_message(map_msg, ierr, MAP_FATAL_78, "Backward difference, x[%d]+delta, row %d, col %d.", j+1, THREE*i, THREE*j);
        return MAP_FATAL;
      };
      ns->jac[THREE*i][THREE*j+2] -= other_type->Fx_connect[i];
      ns->jac[THREE*i][THREE*j+2] /= ns->epsilon;
      ns->jac[THREE*i+1][THREE*j+2] -= other_type->Fy_connect[i];
      ns->jac[THREE*i+1][THREE*j+2] /= ns->epsilon;
      ns->jac[THREE*i+2][THREE*j+2] -= other_type->Fz_connect[i];
      ns->jac[THREE*i+2][THREE*j+2] /= ns->epsilon;
      z_type->z[j] = original_displacement;
    };
  };

  /* read flag to set scaling parameter */
  if (ns->pg) {
    for (i=0 ; i<THREE*z_size ; i++) { 
      ns->jac[i][i] += (ns->ds/pow(ns->iteration_count,1.5)+ns->d);
    };
  };

  return MAP_SAFE;
};


MAP_ERROR_CODE central_difference_jacobian(MAP_OtherStateType_t* other_type, MAP_ParameterType_t* p_type, MAP_ConstraintStateType_t* z_type, ModelData* model_data, char* map_msg, MAP_ERROR_CODE* ierr)
{
  OuterSolveAttributes* ns = &model_data->outer_loop;
  MAP_ERROR_CODE success = MAP_SAFE;
  double original_displacement = 0.0;
  const int THREE = 3;
  const int z_size = z_type->z_Len; //M
  const int m = THREE*(other_type->Fz_connect_Len); // rows
  const int n = THREE*(z_type->z_Len);              // columns
  int i = 0;
  int j = 0;
  
  for (i=0 ; i<z_size ; i++) {
    ns->b[THREE*i] = other_type->Fx_connect[i];
    ns->b[THREE*i+1] = other_type->Fy_connect[i];
    ns->b[THREE*i+2] = other_type->Fz_connect[i];      
  }

  /* First store the connect node force before applying ns->epsilon displacement */
  for (j=0 ; j<z_size ; j++) {   // column               
    for (i=0 ; i<z_size ; i++) { // rows           
      original_displacement = z_type->x[j];
      z_type->x[j] += ns->epsilon;
      success = line_solve_sequence(model_data, p_type, 0.0, map_msg, ierr);
      if (success) {
        set_universal_error_with_message(map_msg, ierr, MAP_FATAL_78, "Central difference, x[%d]+delta, row %d, col %d.", j+1, THREE*i, THREE*j);
        return MAP_FATAL;
      };
      ns->jac[THREE*i][THREE*j] = other_type->Fx_connect[i];      
      ns->jac[THREE*i+1][THREE*j] = other_type->Fy_connect[i];
      ns->jac[THREE*i+2][THREE*j] = other_type->Fz_connect[i];
      z_type->x[j] = original_displacement;
        
      original_displacement = z_type->y[j];
      z_type->y[j] += ns->epsilon;
      if (success) {
        set_universal_error_with_message(map_msg, ierr, MAP_FATAL_78, "Central difference, x[%d]+delta, row %d, col %d.", j+1, THREE*i, THREE*j);
        return MAP_FATAL;
      };
      success = line_solve_sequence(model_data, p_type, 0.0, map_msg, ierr);
      ns->jac[THREE*i][THREE*j+1] = other_type->Fx_connect[i];
      ns->jac[THREE*i+1][THREE*j+1] = other_type->Fy_connect[i];
      ns->jac[THREE*i+2][THREE*j+1] = other_type->Fz_connect[i];
      z_type->y[j] = original_displacement;
    
      original_displacement = z_type->z[j];
      z_type->z[j] += ns->epsilon;
      if (success) {
        set_universal_error_with_message(map_msg, ierr, MAP_FATAL_78, "Central difference, x[%d]+delta, row %d, col %d.", j+1, THREE*i, THREE*j);
        return MAP_FATAL;
      };
      success = line_solve_sequence(model_data, p_type, 0.0, map_msg, ierr);
      ns->jac[THREE*i][THREE*j+2] = other_type->Fx_connect[i];
      ns->jac[THREE*i+1][THREE*j+2] = other_type->Fy_connect[i];
      ns->jac[THREE*i+2][THREE*j+2] = other_type->Fz_connect[i];
      z_type->z[j] = original_displacement;
    };
  };
    
  for (j=0 ; j<z_size ; j++) {   // column               
    for (i=0 ; i<z_size ; i++) { // rows           
      original_displacement = z_type->x[j];
      z_type->x[j] -= ns->epsilon;
      if (success) {
        set_universal_error_with_message(map_msg, ierr, MAP_FATAL_78, "Central difference, x[%d]+delta, row %d, col %d.", j+1, THREE*i, THREE*j);
        return MAP_FATAL;
      };
      success = line_solve_sequence(model_data, p_type, 0.0, map_msg, ierr);
      ns->jac[THREE*i][THREE*j] -= other_type->Fx_connect[i];
      ns->jac[THREE*i][THREE*j] /= (2*ns->epsilon);
      ns->jac[THREE*i+1][THREE*j] -= other_type->Fy_connect[i];
      ns->jac[THREE*i+1][THREE*j] /= (2*ns->epsilon);
      ns->jac[THREE*i+2][THREE*j] -= other_type->Fz_connect[i];
      ns->jac[THREE*i+2][THREE*j] /= (2*ns->epsilon);
      z_type->x[j] = original_displacement;
        
      original_displacement = z_type->y[j];
      z_type->y[j] -= ns->epsilon;
      success = line_solve_sequence(model_data, p_type, 0.0, map_msg, ierr);
      if (success) {
        set_universal_error_with_message(map_msg, ierr, MAP_FATAL_78, "Central difference, x[%d]+delta, row %d, col %d.", j+1, THREE*i, THREE*j);
        return MAP_FATAL;
      };
      ns->jac[THREE*i][THREE*j+1] -= other_type->Fx_connect[i];
      ns->jac[THREE*i][THREE*j+1] /= (2*ns->epsilon);
      ns->jac[THREE*i+1][THREE*j+1] -= other_type->Fy_connect[i];
      ns->jac[THREE*i+1][THREE*j+1] /= (2*ns->epsilon);
      ns->jac[THREE*i+2][THREE*j+1] -= other_type->Fz_connect[i];
      ns->jac[THREE*i+2][THREE*j+1] /= (2*ns->epsilon);
      z_type->y[j] = original_displacement;
    
      original_displacement = z_type->z[j];
      z_type->z[j] -= ns->epsilon;
      success = line_solve_sequence(model_data, p_type, 0.0, map_msg, ierr);
      if (success) {
        set_universal_error_with_message(map_msg, ierr, MAP_FATAL_78, "Central difference, x[%d]+delta, row %d, col %d.", j+1, THREE*i, THREE*j);
        return MAP_FATAL;
      };
      ns->jac[THREE*i][THREE*j+2] -= other_type->Fx_connect[i];
      ns->jac[THREE*i][THREE*j+2] /= (2*ns->epsilon);
      ns->jac[THREE*i+1][THREE*j+2] -= other_type->Fy_connect[i];
      ns->jac[THREE*i+1][THREE*j+2] /= (2*ns->epsilon);
      ns->jac[THREE*i+2][THREE*j+2] -= other_type->Fz_connect[i];
      ns->jac[THREE*i+2][THREE*j+2] /= (2*ns->epsilon);
      z_type->z[j] = original_displacement;
    };
  };

  /* read flag to set scaling parameter */
  if (ns->pg) {
    ns->coef = pow(ns->iteration_count,1.5);
    for (i=0 ; i<THREE*z_size ; i++) { 
      ns->jac[i][i] += (ns->ds/ns->coef + ns->d);
    };
  };

  return MAP_SAFE;
};


MAP_ERROR_CODE call_minpack_lmder(Line* line, InnerSolveAttributes* inner_opt, ModelOptions* opt, const int line_num, const double time, char* map_msg, MAP_ERROR_CODE* ierr)
{
  MAP_ERROR_CODE success = MAP_SAFE;

  /* initial guess vector is set in set_line_initial_guess(..); otherwise, the previous solution is used as the initial guess */
  inner_opt->x[0] = *(line->H.value);
  inner_opt->x[1] = *(line->V.value);

  line->evals = 0;
  line->njac_evals = 0;

  inner_opt->info = __cminpack_func__(lmder)(inner_function_evals, 
                                      line, 
                                      inner_opt->m, 
                                      inner_opt->n, 
                                      inner_opt->x, 
                                      inner_opt->fvec, 
                                      inner_opt->fjac, 
                                      inner_opt->ldfjac, 
                                      inner_opt->f_tol, 
                                      inner_opt->x_tol, 
                                      inner_opt->g_tol, 
                                      inner_opt->max_its, 
                                      inner_opt->diag,
                                      inner_opt->mode, 
                                      inner_opt->factor, 
                                      inner_opt->nprint, 
                                      &line->evals, 
                                      &line->njac_evals, 
                                      inner_opt->ipvt, 
                                      inner_opt->qtf, 
                                      inner_opt->wa1 ,
                                      inner_opt->wa2 ,
                                      inner_opt->wa3 , 
                                      inner_opt->wa4);
  
  line->residual_norm = (MapReal)__minpack_func__(enorm)(&inner_opt->m, inner_opt->fvec);

  if (line->residual_norm>1e-3) {
    success = MAP_WARNING;
    set_universal_error_with_message(map_msg, ierr, MAP_WARNING_11, "Line segment %d.", line_num);
  };

  if (line->options.diagnostics_flag && (double)line->diagnostic_type>time ) { 
    printf("\n      %4.3f [sec]  Line %d\n",time, line_num+1);
    printf("      ----------------------------------------------------\n");
    printf("      Residual l2 norm at solution:  %15.7g\n", line->residual_norm);
    printf("      Function evaluations:         %10i\n", line->evals);
    printf("      Jacobian evaluations:         %10i\n", line->njac_evals);
    printf("      Exit parameter                %10i\n\n", inner_opt->info);
  };
  
  *(line->H.value) = inner_opt->x[0];
  *(line->V.value) = inner_opt->x[1];
  line->converge_reason = inner_opt->info;
  
  switch (inner_opt->info) {
  case 0 :
    success = MAP_FATAL;
    set_universal_error_with_message(map_msg, ierr, MAP_FATAL_39, "Line segment %d.", line_num);
    break;
  case 1 :
    success = MAP_SAFE;
    break;
  case 2 :
    success = MAP_SAFE;
    break;
  case 3 :
    success = MAP_SAFE;
    break;
  case 4 :
    success = MAP_SAFE;
    break;
  case 5 :
    success = MAP_FATAL;
    set_universal_error_with_message(map_msg, ierr, MAP_FATAL_40, "Line segment %d.", line_num);
    break;
  case 6 :
    success = MAP_FATAL;
    set_universal_error_with_message(map_msg, ierr, MAP_ERROR_11, "Line segment %d.", line_num);
    break;
  case 7 :
    success = MAP_FATAL;
    set_universal_error_with_message(map_msg, ierr, MAP_ERROR_13, "Line segment %d.", line_num);
    break;
  case 8 :
    success = MAP_FATAL;
    set_universal_error_with_message(map_msg, ierr, MAP_ERROR_12, "Line segment %d.", line_num);
    break;
  default :
    success = MAP_SAFE;
    break;
  };
  return MAP_SAFE;
};
