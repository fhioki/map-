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
#include "numeric.h"

extern const char MAP_ERROR_STRING[][1024];


MAP_ERROR_CODE reset_node_force_to_zero(Domain* domain, char* map_msg, MAP_ERROR_CODE* ierr)
{
  Node* node_iter = NULL;

  list_iterator_start(&domain->node);            /* starting an iteration "session" */
  while (list_iterator_hasnext(&domain->node)) { /* tell whether more values available */
    node_iter = (Node*)list_iterator_next(&domain->node);
    *(node_iter->sum_force_ptr.fx.value) = 0.0;
    *(node_iter->sum_force_ptr.fy.value) = 0.0;
    *(node_iter->sum_force_ptr.fz.value) = 0.0;    
  };
  list_iterator_stop(&domain->node); /* ending the iteration "session" */    
  return MAP_SAFE;
};


MAP_ERROR_CODE calculate_node_sum_force(Domain* domain, MAP_ParameterType_t* p_type)
{
  Line* line_iter = NULL;
  Node* node_iter = NULL;
  const double g = p_type->g;
  const double rho = p_type->rho_sea;
  double psi = 0.0;
  double sum_fx = 0.0;
  double sum_fy = 0.0;
  double sum_fz = 0.0;
  double fx = 0.0;
  double fy = 0.0;
  double fz = 0.0;
  double fx_a = 0.0;
  double fy_a = 0.0;
  double fz_a = 0.0;

  list_iterator_start(&domain->line);         /* starting an iteration "session" */
  while (list_iterator_hasnext(&domain->line)) { /* tell whether more values available */
    line_iter = (Line*)list_iterator_next(&domain->line);
    psi = line_iter->psi;
    fx = *(line_iter->H.value)*cos(psi);
    fx_a = -(line_iter->H_at_anchor)*cos(psi);
  
    fy = *(line_iter->H.value)*sin(psi);
    fy_a = -(line_iter->H_at_anchor)*sin(psi);
    
    fz = *(line_iter->V.value);
    fz_a = -(line_iter->V_at_anchor);
  
    add_to_sum_fx(line_iter->fairlead, fx);
    add_to_sum_fx(line_iter->anchor, fx_a);
  
    add_to_sum_fy(line_iter->fairlead, fy);
    add_to_sum_fy(line_iter->anchor, fy_a);
  
    add_to_sum_fz(line_iter->fairlead, fz);
    add_to_sum_fz(line_iter->anchor, fz_a);
  }; 
  list_iterator_stop(&domain->line); /* ending the iteration "session" */    

  /* This is where we include the externally applied forces on the node. Note that
   *     \sum F_x= \left \{ \mathbf{f}_\textup{lines} \right \}_x-\left \{ \mathbf{f}_\textup{ext} \right \}_x \\
   *     \sum F_y= \left \{ \mathbf{f}_\textup{lines} \right \}_y-\left \{ \mathbf{f}_\textup{ext} \right \}_y \ \
   *     \sum F_z= \left \{ \mathbf{f}_\textup{lines} \right \}_z-\left \{ \mathbf{f}_\textup{ext} \right \}_z - gM_{\textup{app}} + \rho gB_{\textup{app}}
   * The \mathbf{f}_{\textup{lines}} portion is summed in the line iterator above. 
   */
  list_iterator_start(&domain->node);  /* starting an iteration "session" */
  while (list_iterator_hasnext(&domain->node)) { /* tell whether more values available */
    node_iter = (Node*)list_iterator_next(&domain->node);
    if (node_iter->type==CONNECT) {
      sum_fx = -(node_iter->external_force.fx.value);
      sum_fy = -(node_iter->external_force.fy.value);
      sum_fz = -(node_iter->external_force.fz.value - node_iter->M_applied.value*g + node_iter->B_applied.value*rho*g);

      add_to_sum_fx(node_iter, sum_fx);
      add_to_sum_fy(node_iter, sum_fy);
      add_to_sum_fz(node_iter, sum_fz);
     };
  }; 
  list_iterator_stop(&domain->node); /* ending the iteration "session" */    
  return MAP_SAFE;
}


void add_to_sum_fx(Node* node, const double fx)
{
  *(node->sum_force_ptr.fx.value) += fx;
};


void add_to_sum_fy(Node* node, const double fy)
{
  *(node->sum_force_ptr.fy.value) += fy;
};


void add_to_sum_fz(Node* node, const double fz)
{
  *(node->sum_force_ptr.fz.value) += fz;
};


MAP_ERROR_CODE set_line_variables_pre_solve(Domain* domain, char* map_msg, MAP_ERROR_CODE* ierr)
{
  MAP_ERROR_CODE success = MAP_SAFE;
  Line* line_iter = NULL;
  int i = 0;

  list_iterator_start(&domain->line);            /* starting an iteration "session" */
  while (list_iterator_hasnext(&domain->line)) { /* tell whether more values available */ 
    line_iter = (Line*)list_iterator_next(&domain->line);    
        
    /* no fairlead or anchor was not set. End program gracefully. 
     * there is likely an error in the input file 
     */
    if (!line_iter->fairlead || !line_iter->anchor) {
      *ierr = MAP_FATAL;      
      break;
    };

    /* horizontal cable excursion */
    line_iter->l = set_horizontal_excursion(line_iter);
    
    /* vertical cable excursion */
    line_iter->h = set_vertical_excursion(line_iter);

    /* angle of azimuth */
    success = set_psi(line_iter, map_msg, ierr); 
    if (success!=MAP_SAFE) {
      set_universal_error_with_message(map_msg, ierr, MAP_WARNING_6, "Line number %d", i);
    };
    i++;
  };
  list_iterator_stop(&domain->line); /* ending the iteration "session" */    
  MAP_RETURN;//return MAP_SAFE;
};


MAP_ERROR_CODE set_line_variables_post_solve(Domain* domain, char* map_msg, MAP_ERROR_CODE* ierr)
{
  Line* line_iter = NULL;
  double H = 0.0;
  double V = 0.0;
  double Ha = 0.0;
  double Va = 0.0;
  double w = 0.0;
  double Lu = 0.0;
  double Lb = 0.0;
  double cb = 0.0;
  bool contact_flag = false;
   
  list_iterator_start(&domain->line); /* starting an iteration "session" */
  while (list_iterator_hasnext(&domain->line)) { /* tell whether more values available */ 
    line_iter = (Line*)list_iterator_next(&domain->line);
    
    /* altitude angle at fairlead */
    H = *(line_iter->H.value);
    V = *(line_iter->V.value);
    line_iter->alpha = atan2(V, H);
    
    /* tension at fairlead */
    line_iter->T = sqrt(H*H + V*V);

    if (line_iter->options.linear_spring) {
      Ha = H;
      Va = V;
    } else {
      /* line in contact with seabed options */
      w = line_iter->line_property->omega;
      cb = line_iter->line_property->cb;
      contact_flag = line_iter->options.omit_contact;
      Lu = line_iter->Lu.value;
      if (contact_flag==true || w<0.0 || (V-w*Lu)>0.0) { 
        line_iter->Lb = 0.0;
        Ha = H;
        Va = V - w*Lu;
      } else { /* line is touching the seabed */
        Lb = Lu - (V/w);
        line_iter->Lb = Lb;
        ((H-cb*w*Lb)>0.0) ? (Ha = (H-cb*w*Lb)) : (Ha = 0.0);
        Va = 0.0;
      };
    };
    
    /* tension at anchor */
    line_iter->H_at_anchor = Ha;
    line_iter->V_at_anchor = Va;
    line_iter->T_at_anchor = sqrt(Ha*Ha + Va*Va);
    line_iter->alpha_at_anchor = atan2(Va, Ha);    

    /* set global x,y,z force at fairlead and achor */
    line_iter->fx_fairlead = H*cos(line_iter->psi);
    line_iter->fy_fairlead = H*sin(line_iter->psi);
    line_iter->fz_fairlead = V;

    line_iter->fx_anchor = Ha*cos(line_iter->psi);
    line_iter->fy_anchor = Ha*cos(line_iter->psi);
    line_iter->fz_anchor = Va;
  };
  list_iterator_stop(&domain->line); /* ending the iteration "session" */        
  return MAP_SAFE;
};


MAP_ERROR_CODE set_psi(Line* line, char* map_msg, MAP_ERROR_CODE* ierr)
{
  double overlap_value = 0.0;
  const double x_fair = *(line->fairlead->position_ptr.x.value);
  const double y_fair = *(line->fairlead->position_ptr.y.value);
  const double z_fair = *(line->fairlead->position_ptr.z.value);
  const double x_anch = *(line->anchor->position_ptr.x.value);
  const double y_anch = *(line->anchor->position_ptr.y.value);
  const double z_anch = *(line->anchor->position_ptr.z.value);

  overlap_value = sqrt(pow((x_fair-x_anch),2) + pow((y_fair-y_anch),2) +  pow((z_fair-z_anch),2));
  
  /* make sure the demoninator is not zero (or close to it) */
  if (overlap_value<=1e-2) {    
    // @todo: what happens when this returns prematurely?
    return MAP_WARNING;
  };
  
  /* find the angle psi line simply finds the angle psi between the local and global reference frames simply by 
   * evaluating trig relationships
   */
  line->psi = atan2((y_fair-y_anch), (x_fair-x_anch));
  return MAP_SAFE;
};


double set_horizontal_excursion(Line* line)
{
  const double x_fair = *(line->fairlead->position_ptr.x.value);
  const double y_fair = *(line->fairlead->position_ptr.y.value);
  const double x_anch = *(line->anchor->position_ptr.x.value);
  const double y_anch = *(line->anchor->position_ptr.y.value);  

  return sqrt(pow((x_fair-x_anch),2) + pow((y_fair-y_anch),2));
};


double set_vertical_excursion(Line* line)
{
  const double z_fair = *(line->fairlead->position_ptr.z.value);
  const double z_anch = *(line->anchor->position_ptr.z.value);

  return fabs(z_fair - z_anch);
};


MAP_ERROR_CODE set_line_initial_guess(Domain* domain, char* map_msg, MAP_ERROR_CODE* ierr)
{
  Line* line_iter = NULL;
  double lambda = 0.0;
  double length = 0.0;
  double height = 0.0;
  double w = 0.0;
  double Lu = 0.0;

  list_iterator_start(&domain->line);            /* starting an iteration "session" */
  while (list_iterator_hasnext(&domain->line)) { /* tell whether more values available */ 
    line_iter = (Line*)list_iterator_next(&domain->line);    
    w = line_iter->line_property->omega;
    length = line_iter->l;
    height = line_iter->h;
    Lu = line_iter->Lu.value;
    
    /* set initial guess to user defined value (horizontal fairlead force) 
       in the MAP input file 
    */
    if (line_iter->fairlead->external_force.fx.user_initial_guess && line_iter->fairlead->external_force.fy.user_initial_guess) {      
      *(line_iter->H.value) = sqrt(pow(line_iter->fairlead->external_force.fx.value,2) +       
                                   pow(line_iter->fairlead->external_force.fy.value,2));
    } else {    
      /* note: the line horizontal (l) and vertical (h) excursion are previously initialized in set_line_variables_pre_solve(...) */ 
      if (length<=1e-2) {
        lambda = 1000000;
      } else if (sqrt(length*length+height*height)>=Lu) {
        lambda = 0.2;
      } else {
        lambda = sqrt(3*((Lu*Lu - height*height)/(length*length) - 1));
      };

      /* note: omega is previously initialized in function initialize_cable_library_variables(.. ) in mapinit.c */    
      *(line_iter->H.value) = fabs(w*length/(2*lambda));
    };

    /* set initial guess to user defined value (vertical fairlead force) 
       in the MAP input file 
    */
    if (line_iter->fairlead->external_force.fz.user_initial_guess) {
      *(line_iter->V.value) = line_iter->fairlead->external_force.fz.value;      
    } else {
      *(line_iter->V.value) = (w/2)*((height/tanh(lambda)) + Lu);
    };
    
  };
  list_iterator_stop(&domain->line); /* ending the iteration "session" */    
  return MAP_SAFE;
};



MAP_ERROR_CODE line_solve_sequence(Domain* domain, MAP_ParameterType_t* p_type, double t, char* map_msg, MAP_ERROR_CODE* ierr) 
{
  MAP_ERROR_CODE success = MAP_SAFE;
  
  do { 
    success = set_line_variables_pre_solve(domain, map_msg, ierr);
    success = reset_node_force_to_zero(domain, map_msg, ierr);    
    success = solve_line(domain, t, map_msg, ierr); CHECKERRQ(MAP_FATAL_88);
    // /* prematurely terminating the the line solve routine to return back to the
    //    caller function. Do a back tracking on the solution and attempt recovery */
    // if (success==MAP_FATAL) {
    //   return MAP_FATAL_59;
    // };   
    // CHECKERRQ(MAP_FATAL_39); /* @todo: this should be called elsewhere to notifiy 
    //                             users ofpremature termination. This won't be caught
    //                             when solve_line fails. */    
    success = set_line_variables_post_solve(domain, map_msg, ierr);
    success = calculate_node_sum_force(domain, p_type);
  } while (0);

  if (success==MAP_SAFE) {
    return MAP_SAFE;
  } else if (success==MAP_ERROR) {
    return MAP_ERROR;
  } else {
    return MAP_FATAL;
  };
};    


MAP_ERROR_CODE node_solve_sequence(Domain* domain, MAP_ParameterType_t* p_type, MAP_InputType_t* u_type, MAP_ConstraintStateType_t* z_type, MAP_OtherStateType_t* other_type, char* map_msg, MAP_ERROR_CODE* ierr)
{
  OuterSolveAttributes* ns = &domain->outer_loop;
  MAP_ERROR_CODE success = MAP_SAFE;
  Line* line_iter = NULL;
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

  ns->iteration_count = 1;
  do {
    // if (domain->outer_loop.krylov_accelerator) {
    //   succuss = checkpoint();
    // } else {

    error = 0.0;
    success = line_solve_sequence(domain, p_type, 0.0, map_msg, ierr); CHECKERRQ(MAP_FATAL_79);
    switch (ns->fd) {
    case BACKWARD_DIFFERENCE :
      success = backward_difference_jacobian(other_type, p_type, z_type, domain, map_msg, ierr); CHECKERRQ(MAP_FATAL_75);
      break;
    case CENTRAL_DIFFERENCE :
      success = central_difference_jacobian(other_type, p_type, z_type, domain, map_msg, ierr); CHECKERRQ(MAP_FATAL_76);
      break;
    case FORWARD_DIFFERENCE :
      success = forward_difference_jacobian(other_type, p_type, z_type, domain, map_msg, ierr); CHECKERRQ(MAP_FATAL_77);
      break;
    };

    success = line_solve_sequence(domain, p_type, 0.0, map_msg, ierr); CHECKERRQ(MAP_FATAL_78);
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
    ns->iteration_count++;
    if (ns->iteration_count>ns->max_its) {
      set_universal_error(map_msg, ierr, MAP_FATAL_80);
      break;
    };
    
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


MAP_ERROR_CODE solve_line(Domain* domain, double time, char* map_msg, MAP_ERROR_CODE* ierr)
{
  MAP_ERROR_CODE success = MAP_SAFE;
  Line* line_iter = NULL;
  int n = 1; /* line counter */

  list_iterator_start(&domain->line);            /* starting an iteration "session" */
  while (list_iterator_hasnext(&domain->line)) { /* tell whether more values available */ 
    line_iter = (Line*)list_iterator_next(&domain->line);

    if (line_iter->l<MAP_HORIZONTAL_TOL && line_iter->l>=0.0) { /* perfectly vertical */
      /* this should be triggered for  perfectly vertical cable */
      line_iter->l = MAP_HORIZONTAL_TOL;
    } else if (line_iter->l<0.0) {
      set_universal_error_with_message(map_msg, ierr, MAP_FATAL_54, "Line segment %d, l = %d [m].", n, line_iter->l);
      break; 
    } else if (line_iter->h<=-MACHINE_EPSILON) {
      set_universal_error_with_message(map_msg, ierr, MAP_FATAL_55, "Line segment %d, h = %d [m].", n, line_iter->h);
      break; 
    } else if (line_iter->line_property->omega>0.0) {
      success = check_maximum_line_length(line_iter, line_iter->options.omit_contact, map_msg, ierr);
      if (success) {        
        set_universal_error_with_message(map_msg, ierr, MAP_FATAL_59, "Line segment %d.", n);
        break;
      };
    };    
    if (line_iter->options.linear_spring) {
      success = solve_linear_spring_cable(line_iter, map_msg, ierr); CHECKERRQ(MAP_FATAL_87);
    } else {
      success = call_minpack_lmder(line_iter, &domain->inner_loop, &domain->model_options, n, time, map_msg, ierr); CHECKERRQ(MAP_FATAL_79);
    };

    /* 
       @todo: this should be moved to outside the loop 
    */
    // /* check if L^2 norm is small. If not, MAP converged prematurely */
    // if (line_iter->residual_norm>1e-3) {
    //   set_universal_error_with_message(map_msg, ierr, MAP_FATAL_90, "Line segment %d.", n);
    //   break;      
    // };
    n++;
  };
  list_iterator_stop(&domain->line); /* ending the iteration "session" */    

  if (*ierr==MAP_SAFE) {
    return MAP_SAFE;
  } else if (*ierr==MAP_WARNING) {
    return MAP_SAFE;
  } else if (*ierr==MAP_ERROR) {
    return MAP_ERROR;
  } else {
    return MAP_FATAL;
  };
};


MAP_ERROR_CODE solve_linear_spring_cable(Line* line, char* map_msg, MAP_ERROR_CODE* ierr)
{
  MAP_ERROR_CODE success = MAP_SAFE;
  double d_norm = 0.0;
  Vector d;  /* r1-r2 */
  Vector r1; /* anchor position */
  Vector r2; /* fairlead position */
  Vector force;
  Vector d_unit;
  const double EA = line->line_property->EA;
  const double Lu = line->Lu.value;  

  /* line anchor position */
  r1.x = *line->anchor->position_ptr.x.value; 
  r1.y = *line->anchor->position_ptr.y.value; 
  r1.z = *line->anchor->position_ptr.z.value; 

  /* line fairlead position */
  r2.x = *line->fairlead->position_ptr.x.value; 
  r2.y = *line->fairlead->position_ptr.y.value; 
  r2.z = *line->fairlead->position_ptr.z.value; 

  /* difference vector */
  d.x =  r1.x - r2.x;
  d.y =  r1.y - r2.y;
  d.z =  r1.z - r2.z;
  
  /* calculate norm of \mathbf{d} */
  d_norm = sqrt(d.x*d.x + d.y*d.y + d.z*d.z);

  if (d_norm<MACHINE_EPSILON) {
    return MAP_FATAL;
  }

  /* \hat{u}_{\mathbf{d}} = \frac{\mathbf{d}} {\left \| \mathbf{d} \right \|} */
  d_unit.x = d.x/d_norm; 
  d_unit.y = d.y/d_norm; 
  d_unit.z = d.z/d_norm; 

  if (d_norm<=Lu) { /* a cable cannot support a compression; force will be negative */
    force.x = 0.0;
    force.y = 0.0;
    force.z = 0.0;
  } else {
    force.x = EA/Lu*(d_norm - Lu)*d_unit.x;
    force.y = EA/Lu*(d_norm - Lu)*d_unit.y;
    force.z = EA/Lu*(d_norm - Lu)*d_unit.z;
  };

  *(line->H.value) = sqrt(force.x*force.x + force.y*force.y);
  *(line->V.value) = fabs(force.z);
  return MAP_SAFE;
};


MAP_ERROR_CODE check_maximum_line_length(Line* line, const bool contact_flag, char *map_msg, MAP_ERROR_CODE *ierr)
{
  MAP_ERROR_CODE success = MAP_SAFE;
  double LMax = 0.0;
  const double l = line->l;
  const double h = line->h;
  const double EA = line->line_property->EA;
  const double W = line->line_property->omega;
  const double Lu = line->Lu.value;

  LMax = l - EA/W + sqrt(pow((EA/W),2) + 2.0*h*EA/W);
  if (Lu>=LMax && contact_flag==false) {
    set_universal_error_with_message(map_msg, ierr, MAP_FATAL_59, "LMax = %f [m].", LMax);
    return MAP_FATAL;
  };
  return MAP_SAFE;
};


MAP_ERROR_CODE increment_dof_by_delta(double* u_type, const double delta, const int size)
{
  int i = 0;
  for (i=0 ; i<size ; i++) {
    u_type[i] += delta;
  };
  return MAP_SAFE;
};


MAP_ERROR_CODE restore_original_displacement(double* u_ype, const double* initial_value, const int size)
{
  int i = 0;
  for (i=0 ; i<size ; i++) {
    u_ype[i] = initial_value[i];
  };
  return MAP_SAFE;
};


MAP_ERROR_CODE reset_force_to_zero(double* fx, double* fy, double* fz, double* mx, double* my, double* mz, const int size)
{
  int i = 0;
  for (i=0 ; i<size ; i++) {
    fx[i] = 0.0;
    fy[i] = 0.0;
    fz[i] = 0.0;
    mx[i] = 0.0;
    my[i] = 0.0;
    mz[i] = 0.0;
  };
  return MAP_SAFE;
};



MAP_ERROR_CODE set_force_minus(const double* u_type, double* force, const int size)
{
  int i = 0;
  for (i=0 ; i<size ; i++) {
    force[i] -= -u_type[i];
  };
  return MAP_SAFE;
};


MAP_ERROR_CODE set_force_plus(const double* input_type, double* force, const int size)
{
  int i = 0;
  for (i=0 ; i<size ; i++) {
    force[i] += -input_type[i];
  };
  return MAP_SAFE;
};



MAP_ERROR_CODE set_moment_minus(const MAP_OutputType_t* y_type, const Vessel* vessel, double* mx, double* my, double* mz, const int size)
{
  int i = 0;
  double rx = 0.0;
  double ry = 0.0;
  double rz = 0.0;

  for (i=0 ; i<size ; i++) {    
    /* @todo: this is not consistent with set_moment_minus_2... 
     *        ensure moments are in the global frame
     */
    rx = (vessel->xi[i] - vessel->ref_origin.x.value);
    ry = (vessel->yi[i] - vessel->ref_origin.y.value);
    rz = (vessel->zi[i] - vessel->ref_origin.z.value);    

    /* cross  product: 
     * \mathbf{m}=\mathbf{r} \times \mathbf{F} 
     */
    mx[i] -= ((-y_type->Fz[i]*ry) - (-y_type->Fy[i]*rz));
    my[i] -= ((-y_type->Fx[i]*rz) - (-y_type->Fz[i]*rx));
    mz[i] -= ((-y_type->Fy[i]*rx) - (-y_type->Fx[i]*ry));
  };
  return MAP_SAFE;
}



MAP_ERROR_CODE set_moment_plus(const MAP_OutputType_t* output_type, const Vessel* vessel, double* mx, double* my, double* mz, const int size)
{
  int i = 0;
  double rx = 0.0;
  double ry = 0.0;
  double rz = 0.0;

  for (i=0 ; i<size ; i++) {    
    rx = (vessel->xi[i] - vessel->ref_origin.x.value);
    ry = (vessel->yi[i] - vessel->ref_origin.y.value);
    rz = (vessel->zi[i] - vessel->ref_origin.z.value);    

    /* cross  product: 
     * \mathbf{m}=\mathbf{r} \times \mathbf{F} 
     */
    mx[i] += ((-output_type->Fz[i]*ry) - (-output_type->Fy[i]*rz));
    my[i] += ((-output_type->Fx[i]*rz) - (-output_type->Fz[i]*rx));
    mz[i] += ((-output_type->Fy[i]*rx) - (-output_type->Fx[i]*ry));
  };
  return MAP_SAFE;
};



MAP_ERROR_CODE increment_phi_dof_by_delta(MAP_InputType_t* u_type, const Vessel* vessel, const double delta, const int size)
{
  int i = 0;
  double R[3][3];
  double phi = vessel->orientation.phi.value;
  double rx = 0.0;
  double ry = 0.0;
  double rz = 0.0;
  R[0][0] = 1.0;   R[0][1] = 0.0;              R[0][2] = 0.0;
  R[1][0] = 0.0;   R[1][1] = cos(phi+delta);   R[1][2] = -sin(phi+delta);  
  R[2][0] = 0.0;   R[2][1] = sin(phi+delta);   R[2][2] = cos(phi+delta);

  for (i=0 ; i<size ; i++) {
    /* @todo: should also include the reference origin location? */
    rx = vessel->xi[i];
    ry = vessel->yi[i];
    rz = vessel->zi[i];
    // uType->y[i] = ry - rz*delta;
    // uType->z[i] = rz + ry*delta;
    // printf("%f  %f  %f\n",uType->x[i],uType->y[i],uType->z[i]);
    u_type->x[i] = R[0][0]*rx;
    u_type->y[i] = R[1][1]*ry + R[1][2]*rz;
    u_type->z[i] = R[2][1]*ry + R[2][2]*rz;
  };
  return MAP_SAFE;
};


MAP_ERROR_CODE increment_the_dof_by_delta(MAP_InputType_t* u_type, const Vessel* vessel, const double delta, const int size)
{
  int i = 0;
  double R[3][3];
  double the = vessel->orientation.the.value;
  double rx = 0.0;
  double ry = 0.0;
  double rz = 0.0;

  R[0][0] = cos(the+delta);   R[0][1] = 0.0;   R[0][2] = sin(the+delta);
  R[1][0] = 0.0;              R[1][1] = 1.0;   R[1][2] = 0.0;
  R[2][0] = -sin(the+delta);  R[2][1] = 0.0;   R[2][2] = cos(the+delta);
  
  for (i=0 ; i<size ; i++) {
    /* @todo: should also include the reference origin location? */
    rx = vessel->xi[i];
    ry = vessel->yi[i];
    rz = vessel->zi[i];
    // u_type->x[i] = rx + rz*delta;
    // u_type->z[i] = rz - rx*delta;    
    u_type->x[i] = R[0][0]*rx + R[0][2]*rz;
    u_type->y[i] = R[1][1]*ry;
    u_type->z[i] = R[2][0]*rx + R[2][2]*rz;
  };
  return MAP_SAFE;
};


MAP_ERROR_CODE increment_psi_dof_by_delta(MAP_InputType_t* u_type, const Vessel* vessel, const double delta, const int size)
{
  int i = 0;
  double R[3][3];
  double psi = vessel->orientation.psi.value;
  double rx = 0.0;
  double ry = 0.0;
  double rz = 0.0;
  R[0][0] = cos(psi+delta);   R[0][1] = -sin(psi+delta);   R[0][2] = 0.0;
  R[1][0] = sin(psi+delta);   R[1][1] = cos(psi+delta);    R[1][2] = 0.0;
  R[2][0] = 0.0;              R[2][1] = 0.0;               R[2][2] = 1.0;

  for (i=0 ; i<size ; i++) {
    /* @todo: should also include the reference origin location? */
    rx = vessel->xi[i];
    ry = vessel->yi[i];
    rz = vessel->zi[i];
    // u_type->x[i] = rx - ry*delta;
    // u_type->y[i] = ry + rx*delta;    
    u_type->x[i] = R[0][0]*rx + R[0][1]*ry;
    u_type->y[i] = R[1][0]*rx + R[1][1]*ry;
    u_type->z[i] = R[2][2]*rz;
  };
  return MAP_SAFE;
};



MAP_ERROR_CODE fd_x_sequence(MAP_OtherStateType_t* other_type, MAP_ParameterType_t* p_type, MAP_InputType_t* u_type, MAP_OutputType_t* y_type, MAP_ConstraintStateType_t* z_type, Fd* force, const double epsilon, const int size, const double* original_pos, char* map_msg, MAP_ERROR_CODE* ierr)
{
  MAP_ERROR_CODE success = MAP_SAFE;
  Domain* domain = other_type->object;
  Vessel* vessel = &domain->vessel;

  do {
    /* minus epsilon sequence */
    success = increment_dof_by_delta(u_type->x, -epsilon, size); CHECKERRQ(MAP_FATAL_61);
    if (domain->MAP_SOLVE_TYPE==MONOLITHIC) {
      success = line_solve_sequence(domain, p_type, 0.0, map_msg, ierr);
    } else {
      success = node_solve_sequence(domain, p_type, u_type, z_type, other_type, map_msg, ierr); // @todo CHECKERRQ()
    };    
    success = set_force_plus(y_type->Fx, force->fx, size); CHECKERRQ(MAP_FATAL_61);
    success = set_force_plus(y_type->Fy, force->fy, size); CHECKERRQ(MAP_FATAL_61);
    success = set_force_plus(y_type->Fz, force->fz, size); CHECKERRQ(MAP_FATAL_61);
    success = set_moment_plus(y_type, vessel, force->mx, force->my, force->mz, size); CHECKERRQ(MAP_FATAL_61);
    success = restore_original_displacement(u_type->x, original_pos, size); CHECKERRQ(MAP_FATAL_61);
  
    /* plus epsilon sequence */
    success = increment_dof_by_delta(u_type->x, epsilon, size); CHECKERRQ(MAP_FATAL_61);
    if (domain->MAP_SOLVE_TYPE==MONOLITHIC) {
      success = line_solve_sequence(domain, p_type, 0.0, map_msg, ierr);
    } else {
      success = node_solve_sequence(domain, p_type, u_type, z_type, other_type, map_msg, ierr); // @todo CHECKERRQ()
    };    
    success = set_force_minus(y_type->Fx, force->fx, size); CHECKERRQ(MAP_FATAL_61);
    success = set_force_minus(y_type->Fy, force->fy, size); CHECKERRQ(MAP_FATAL_61);
    success = set_force_minus(y_type->Fz, force->fz, size); CHECKERRQ(MAP_FATAL_61);
    success = set_moment_minus(y_type, vessel, force->mx, force->my, force->mz, size); CHECKERRQ(MAP_FATAL_61);
    success = restore_original_displacement(u_type->x, original_pos, size); CHECKERRQ(MAP_FATAL_61);
  } while (0);
  MAP_RETURN;
};


MAP_ERROR_CODE fd_y_sequence(MAP_OtherStateType_t* other_type, MAP_ParameterType_t* p_type, MAP_InputType_t* u_type, MAP_OutputType_t* y_type, MAP_ConstraintStateType_t* z_type, Fd* force, const double epsilon, const int size, const double* original_pos, char* map_msg, MAP_ERROR_CODE* ierr)
{
  MAP_ERROR_CODE success = MAP_SAFE;
  Domain* domain = other_type->object;
  Vessel* vessel = &domain->vessel;

  do {
    /* minus epsilon sequence */
    success = increment_dof_by_delta(u_type->y, -epsilon, size); CHECKERRQ(MAP_FATAL_61);
    if (domain->MAP_SOLVE_TYPE==MONOLITHIC) {
      success = line_solve_sequence(domain, p_type, 0.0, map_msg, ierr);
    } else {
      success = node_solve_sequence(domain, p_type, u_type, z_type, other_type, map_msg, ierr); // @todo CHECKERRQ()
    };    
    success = set_force_plus(y_type->Fx, force->fx, size); CHECKERRQ(MAP_FATAL_61);
    success = set_force_plus(y_type->Fy, force->fy, size); CHECKERRQ(MAP_FATAL_61);
    success = set_force_plus(y_type->Fz, force->fz, size); CHECKERRQ(MAP_FATAL_61);
    success = set_moment_plus(y_type, vessel, force->mx, force->my, force->mz, size); CHECKERRQ(MAP_FATAL_61);
    success = restore_original_displacement(u_type->y, original_pos, size); CHECKERRQ(MAP_FATAL_61);
        
    /* plus epsilon sequence */
    success = increment_dof_by_delta(u_type->y, epsilon, size); CHECKERRQ(MAP_FATAL_61);
    if (domain->MAP_SOLVE_TYPE==MONOLITHIC) {
      success = line_solve_sequence(domain, p_type, 0.0, map_msg, ierr);
    } else {
      success = node_solve_sequence(domain, p_type, u_type, z_type, other_type, map_msg, ierr); // @todo CHECKERRQ()
    };    
    success = set_force_minus(y_type->Fx, force->fx, size); CHECKERRQ(MAP_FATAL_61);
    success = set_force_minus(y_type->Fy, force->fy, size); CHECKERRQ(MAP_FATAL_61);
    success = set_force_minus(y_type->Fz, force->fz, size); CHECKERRQ(MAP_FATAL_61);
    success = set_moment_minus(y_type, vessel, force->mx, force->my, force->mz, size); CHECKERRQ(MAP_FATAL_61);
    success = restore_original_displacement(u_type->y, original_pos, size); CHECKERRQ(MAP_FATAL_61);
  } while (0);
  
  MAP_RETURN;
};


MAP_ERROR_CODE fd_z_sequence(MAP_OtherStateType_t* other_type, MAP_ParameterType_t* p_type, MAP_InputType_t* u_type, MAP_OutputType_t* y_type, MAP_ConstraintStateType_t* z_type, Fd* force, const double epsilon, const int size, const double* original_pos, char* map_msg, MAP_ERROR_CODE* ierr)
{
  MAP_ERROR_CODE success = MAP_SAFE;
  Domain* domain = other_type->object;
  Vessel* vessel = &domain->vessel;

  do {
    /* minus epsilon sequence */
    success = increment_dof_by_delta(u_type->z, -epsilon, size); CHECKERRQ(MAP_FATAL_61);
    if (domain->MAP_SOLVE_TYPE==MONOLITHIC) {
      success = line_solve_sequence(domain, p_type, 0.0, map_msg, ierr);
    } else {
      success = node_solve_sequence(domain, p_type, u_type, z_type, other_type, map_msg, ierr); // @todo CHECKERRQ()
    };    
    success = set_force_plus(y_type->Fx, force->fx, size); CHECKERRQ(MAP_FATAL_61);
    success = set_force_plus(y_type->Fy, force->fy, size); CHECKERRQ(MAP_FATAL_61);
    success = set_force_plus(y_type->Fz, force->fz, size); CHECKERRQ(MAP_FATAL_61);
    success = set_moment_plus(y_type, vessel, force->mx, force->my, force->mz, size); CHECKERRQ(MAP_FATAL_61);
    success = restore_original_displacement(u_type->z, original_pos, size); CHECKERRQ(MAP_FATAL_61);
        
    /* plus epsilon sequence */
    success = increment_dof_by_delta(u_type->z, epsilon, size); CHECKERRQ(MAP_FATAL_61);
    if (domain->MAP_SOLVE_TYPE==MONOLITHIC) {
      success = line_solve_sequence(domain, p_type, 0.0, map_msg, ierr);
    } else {
      success = node_solve_sequence(domain, p_type, u_type, z_type, other_type, map_msg, ierr); // @todo CHECKERRQ()
    };    
    success = set_force_minus(y_type->Fx, force->fx, size); CHECKERRQ(MAP_FATAL_61);
    success = set_force_minus(y_type->Fy, force->fy, size); CHECKERRQ(MAP_FATAL_61);
    success = set_force_minus(y_type->Fz, force->fz, size); CHECKERRQ(MAP_FATAL_61);
    success = set_moment_minus(y_type, vessel, force->mx, force->my, force->mz, size); CHECKERRQ(MAP_FATAL_61);
    success = restore_original_displacement(u_type->z, original_pos, size); CHECKERRQ(MAP_FATAL_61);
  } while (0);

  MAP_RETURN;
};


MAP_ERROR_CODE fd_phi_sequence(MAP_OtherStateType_t* other_type, MAP_ParameterType_t* p_type, MAP_InputType_t* u_type, MAP_OutputType_t* y_type, MAP_ConstraintStateType_t* z_type, Fd* force, const double epsilon, const int size, const double* original_x, const double* original_y, const double* original_z, char* map_msg, MAP_ERROR_CODE* ierr)
{
  MAP_ERROR_CODE success = MAP_SAFE;
  Domain* domain = other_type->object;
  Vessel* vessel = &domain->vessel;

  do {
    /* minus epsilon sequence */
    success = increment_phi_dof_by_delta(u_type, vessel, -epsilon, size); CHECKERRQ(MAP_FATAL_61);        
    if (domain->MAP_SOLVE_TYPE==MONOLITHIC) {
      success = line_solve_sequence(domain, p_type, 0.0, map_msg, ierr);
    } else {
      success = node_solve_sequence(domain, p_type, u_type, z_type, other_type, map_msg, ierr); // @todo CHECKERRQ()
    };    
    success = set_force_plus(y_type->Fx, force->fx, size); CHECKERRQ(MAP_FATAL_61);
    success = set_force_plus(y_type->Fy, force->fy, size); CHECKERRQ(MAP_FATAL_61);
    success = set_force_plus(y_type->Fz, force->fz, size); CHECKERRQ(MAP_FATAL_61);
    success = set_moment_plus_2(u_type, y_type, vessel, force->mx, force->my, force->mz, size); CHECKERRQ(MAP_FATAL_61);
    success = restore_original_displacement(u_type->x, original_x, size); CHECKERRQ(MAP_FATAL_61);        
    success = restore_original_displacement(u_type->y, original_y, size); CHECKERRQ(MAP_FATAL_61);        
    success = restore_original_displacement(u_type->z, original_z, size); CHECKERRQ(MAP_FATAL_61);        
    
    /* plus epsilon sequence */
    success = increment_phi_dof_by_delta(u_type, vessel, epsilon, size); CHECKERRQ(MAP_FATAL_61);        
    if (domain->MAP_SOLVE_TYPE==MONOLITHIC) {
      success = line_solve_sequence(domain, p_type, 0.0, map_msg, ierr);
    } else {
      success = node_solve_sequence(domain, p_type, u_type, z_type, other_type, map_msg, ierr); // @todo CHECKERRQ()
    };    
    success = set_force_minus(y_type->Fx, force->fx, size); CHECKERRQ(MAP_FATAL_61);
    success = set_force_minus(y_type->Fy, force->fy, size); CHECKERRQ(MAP_FATAL_61);
    success = set_force_minus(y_type->Fz, force->fz, size); CHECKERRQ(MAP_FATAL_61);
    success = set_moment_minus_2(u_type, y_type, vessel, force->mx, force->my, force->mz, size); CHECKERRQ(MAP_FATAL_61);
    success = restore_original_displacement(u_type->x, original_x, size); CHECKERRQ(MAP_FATAL_61);        
    success = restore_original_displacement(u_type->y, original_y, size); CHECKERRQ(MAP_FATAL_61);        
    success = restore_original_displacement(u_type->z, original_y, size); CHECKERRQ(MAP_FATAL_61);                
  } while (0);
  
  MAP_RETURN;
};


MAP_ERROR_CODE fd_the_sequence(MAP_OtherStateType_t* other_type, MAP_ParameterType_t* p_type, MAP_InputType_t* u_type, MAP_OutputType_t* y_type, MAP_ConstraintStateType_t* z_type, Fd* force, const double epsilon, const int size, const double* original_x, const double* original_y, const double* original_z, char* map_msg, MAP_ERROR_CODE* ierr)
{
  MAP_ERROR_CODE success = MAP_SAFE;
  Domain* domain = other_type->object;
  Vessel* vessel = &domain->vessel;

  do {
    /* minus epsilon sequence */
    success = increment_the_dof_by_delta(u_type, vessel, -epsilon, size); CHECKERRQ(MAP_FATAL_61);
    if (domain->MAP_SOLVE_TYPE==MONOLITHIC) {
      success = line_solve_sequence(domain, p_type, 0.0, map_msg, ierr);
    } else {
      success = node_solve_sequence(domain, p_type, u_type, z_type, other_type, map_msg, ierr); // @todo CHECKERRQ()
    };    
    success = set_force_plus(y_type->Fx, force->fx, size); CHECKERRQ(MAP_FATAL_61);
    success = set_force_plus(y_type->Fy, force->fy, size); CHECKERRQ(MAP_FATAL_61);
    success = set_force_plus(y_type->Fz, force->fz, size); CHECKERRQ(MAP_FATAL_61);
    success = set_moment_plus_2(u_type, y_type, vessel, force->mx, force->my, force->mz, size); CHECKERRQ(MAP_FATAL_61);
    success = restore_original_displacement(u_type->x, original_x, size); CHECKERRQ(MAP_FATAL_61);        
    success = restore_original_displacement(u_type->y, original_y, size); CHECKERRQ(MAP_FATAL_61);        
    success = restore_original_displacement(u_type->z, original_z, size); CHECKERRQ(MAP_FATAL_61);        
    
    /* plut epsilon sequence */
    success = increment_the_dof_by_delta(u_type, vessel, epsilon, size); CHECKERRQ(MAP_FATAL_61);        
    if (domain->MAP_SOLVE_TYPE==MONOLITHIC) {
      success = line_solve_sequence(domain, p_type, 0.0, map_msg, ierr);
    } else {
      success = node_solve_sequence(domain, p_type, u_type, z_type, other_type, map_msg, ierr); // @todo CHECKERRQ()
    };    
    success = set_force_minus(y_type->Fx, force->fx, size); CHECKERRQ(MAP_FATAL_61);
    success = set_force_minus(y_type->Fy, force->fy, size); CHECKERRQ(MAP_FATAL_61);
    success = set_force_minus(y_type->Fz, force->fz, size); CHECKERRQ(MAP_FATAL_61);
    success = set_moment_minus_2(u_type, y_type, vessel, force->mx, force->my, force->mz, size); CHECKERRQ(MAP_FATAL_61);
    success = restore_original_displacement(u_type->x, original_x, size); CHECKERRQ(MAP_FATAL_61);        
    success = restore_original_displacement(u_type->y, original_y, size); CHECKERRQ(MAP_FATAL_61);        
    success = restore_original_displacement(u_type->z, original_z, size); CHECKERRQ(MAP_FATAL_61);        
  } while (0);

  MAP_RETURN;
};


MAP_ERROR_CODE fd_psi_sequence(MAP_OtherStateType_t* other_type, MAP_ParameterType_t* p_type, MAP_InputType_t* u_type, MAP_OutputType_t* y_type, MAP_ConstraintStateType_t* z_type, Fd* force, const double epsilon, const int size, const double* original_x, const double* original_y, const double* original_z, char* map_msg, MAP_ERROR_CODE* ierr)
{
  MAP_ERROR_CODE success = MAP_SAFE;
  Domain* domain = other_type->object;
  Vessel* vessel = &domain->vessel;
  do {
    /* minus epsilon sequence */
    success = increment_psi_dof_by_delta(u_type, vessel, -epsilon, size); CHECKERRQ(MAP_FATAL_61);
    if (domain->MAP_SOLVE_TYPE==MONOLITHIC) {
      success = line_solve_sequence(domain, p_type, 0.0, map_msg, ierr);
    } else {
      success = node_solve_sequence(domain, p_type, u_type, z_type, other_type, map_msg, ierr); // @todo CHECKERRQ()
    };    
    success = set_force_plus(y_type->Fx, force->fx, size); CHECKERRQ(MAP_FATAL_61);
    success = set_force_plus(y_type->Fy, force->fy, size); CHECKERRQ(MAP_FATAL_61);
    success = set_force_plus(y_type->Fz, force->fz, size); CHECKERRQ(MAP_FATAL_61);
    success = set_moment_plus_2(u_type, y_type, vessel, force->mx, force->my, force->mz, size); CHECKERRQ(MAP_FATAL_61);
    success = restore_original_displacement(u_type->x, original_x, size); CHECKERRQ(MAP_FATAL_61);
    success = restore_original_displacement(u_type->y, original_y, size); CHECKERRQ(MAP_FATAL_61);
    success = restore_original_displacement(u_type->z, original_z, size); CHECKERRQ(MAP_FATAL_61);
        
    /* plut epsilon sequence */
    success = increment_psi_dof_by_delta(u_type, vessel, epsilon, size); CHECKERRQ(MAP_FATAL_61);
    if (domain->MAP_SOLVE_TYPE==MONOLITHIC) {
      success = line_solve_sequence(domain, p_type, 0.0, map_msg, ierr);
    } else {
      success = node_solve_sequence(domain, p_type, u_type, z_type, other_type, map_msg, ierr); // @todo CHECKERRQ()
    };    
    success = set_force_minus(y_type->Fx, force->fx, size); CHECKERRQ(MAP_FATAL_61);
    success = set_force_minus(y_type->Fy, force->fy, size); CHECKERRQ(MAP_FATAL_61);
    success = set_force_minus(y_type->Fz, force->fz, size); CHECKERRQ(MAP_FATAL_61);
    success = set_moment_minus_2(u_type, y_type, vessel, force->mx, force->my, force->mz, size); CHECKERRQ(MAP_FATAL_61);
    success = restore_original_displacement(u_type->x, original_x, size); CHECKERRQ(MAP_FATAL_61);
    success = restore_original_displacement(u_type->y, original_y, size); CHECKERRQ(MAP_FATAL_61);
    success = restore_original_displacement(u_type->z, original_z, size); CHECKERRQ(MAP_FATAL_61);
  } while (0);

  MAP_RETURN;
};



MAP_ERROR_CODE calculate_stiffness(double* K, Fd* force, const double delta, const int size)
{  
  int i = 0;
  for (i=0 ; i<size ; i++) {
    K[0] += (force->fx[i]/(2*delta));
    K[1] += (force->fy[i]/(2*delta));
    K[2] += (force->fz[i]/(2*delta));
    K[3] += (force->mx[i]/(2*delta));
    K[4] += (force->my[i]/(2*delta));
    K[5] += (force->mz[i]/(2*delta));
  };
  return MAP_SAFE;
};


MAP_ERROR_CODE set_moment_plus_2(const MAP_InputType_t* u_type, const MAP_OutputType_t* y_type, const Vessel* vessel, double* mx, double* my, double* mz, const int size)
{
  int i = 0;
  double rx = 0.0;
  double ry = 0.0;
  double rz = 0.0;
  double m1 = 0.0;
  double m2 = 0.0;
  double m3 = 0.0;
  double phi = vessel->orientation.phi.value;

  double delta = 1e-4;
  for (i=0 ; i<size ; i++) {    
    /* @todo: forces are in global frame. Make sure everything is consistent. Should also add
     *        reference origin. This need s to be converted to the global frame to be consistent with
     *       uType. 
     */
    rx = u_type->x[i];// (vessel->xi[i] - vessel->ref_origin.x.value);
    ry = u_type->y[i];// (vessel->yi[i] - vessel->ref_origin.y.value);
    rz = u_type->z[i];// (vessel->zi[i] - vessel->ref_origin.z.value);    
    m1 = ((-y_type->Fz[i]*ry) - (-y_type->Fy[i]*rz));
    m2 = ((-y_type->Fx[i]*rz) - (-y_type->Fz[i]*rx));
    m3 = ((-y_type->Fy[i]*rx) - (-y_type->Fx[i]*ry));
    mx[i] += m1;
    my[i] += m2*cos(phi+delta) - m3*sin(phi+delta);
    mz[i] += m2*sin(phi+delta) + m3*cos(phi+delta);    
    // mx[i] += ((-output_type->Fz[i]*ry) - (-output_type->Fy[i]*rz));
    // my[i] += ((-output_type->Fx[i]*rz) - (-output_type->Fz[i]*rx));
    // mz[i] += ((-output_type->Fy[i]*rx) - (-output_type->Fx[i]*ry));
  };
  return MAP_SAFE;
};


MAP_ERROR_CODE set_moment_minus_2(const MAP_InputType_t* u_type, const MAP_OutputType_t* y_type, const Vessel* vessel, double* mx, double* my, double* mz, const int size)
{
  int i = 0;
  double rx = 0.0;
  double ry = 0.0;
  double rz = 0.0;
  double m1 = 0.0;
  double m2 = 0.0;
  double m3 = 0.0;
  double phi = vessel->orientation.phi.value;

  double delta = 1e-4;
  for (i=0 ; i<size ; i++) {    
    rx = u_type->x[i];// (vessel->xi[i] - vessel->ref_origin.x.value);
    ry = u_type->y[i];// (vessel->yi[i] - vessel->ref_origin.y.value);
    rz = u_type->z[i];// (vessel->zi[i] - vessel->ref_origin.z.value);    
    m1 = ((-y_type->Fz[i]*ry) - (-y_type->Fy[i]*rz));
    m2 = ((-y_type->Fx[i]*rz) - (-y_type->Fz[i]*rx));
    m3 = ((-y_type->Fy[i]*rx) - (-y_type->Fx[i]*ry));
    mx[i] -= m1;
    my[i] -= m2*cos(phi+delta) - m3*sin(phi+delta);
    mz[i] -= m2*sin(phi+delta) + m3*cos(phi+delta);    
    // mx[i] -= ((-y_type->Fz[i]*ry) - (-y_type->Fy[i]*rz));
    // my[i] -= ((-y_type->Fx[i]*rz) - (-y_type->Fz[i]*rx));
    // mz[i] -= ((-y_type->Fy[i]*rx) - (-y_type->Fx[i]*ry));
  };
  return MAP_SAFE;
};
