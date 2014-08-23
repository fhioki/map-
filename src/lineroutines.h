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


#ifndef _LINEROUTINES_H
#define _LINEROUTINES_H


#include "map.h"


/**
 * success = incremenet_x_dof_by_delta(uType, -epsilon);
 */
MAP_ERROR_CODE increment_dof_by_delta(double* u_type, const MapReal delta, const int size);


/**
 * success = restore_original_displacement(uType->x, xOriginal, N);
 */
MAP_ERROR_CODE restore_original_displacement(double* u_ype, const double* initial_value, const int size);


/**
 * success = reset_force_to_zero(yType, N);
 */
MAP_ERROR_CODE reset_force_to_zero(double* fx, double* fy, double* fz, double* mx, double* my, double* mz, const int size);


/**
 * success = set_force_plus(yType->Fx, fx, N);
 */
MAP_ERROR_CODE set_force_minus(const double* inputType, double* force, const int size);


/**
 * success = set_force_plus(yType->Fx, fx, N);
 */
MAP_ERROR_CODE set_force_plus(const double* inputType, double* force, const int size);


MAP_ERROR_CODE calculate_stiffness(double* K, double* force, const double delta, const int size);


/**
 * success = set_moment_minus(yType, vessel, mx, my, mz, N)
 */
MAP_ERROR_CODE set_moment_minus(const MAP_OutputType_t* outputType, const Vessel* vessel, double* mx, double* my, double* mz, const int size);


MAP_ERROR_CODE set_moment_plus(const MAP_OutputType_t* outputType, const Vessel* vessel, double* mx, double* my, double* mz, const int size);


MAP_ERROR_CODE increment_phi_dof_by_delta(MAP_InputType_t* u_type, const Vessel* vessel, const double delta, const int size);


MAP_ERROR_CODE increment_the_dof_by_delta(MAP_InputType_t* u_type, const Vessel* vessel, const double delta, const int size);


MAP_ERROR_CODE increment_psi_dof_by_delta(MAP_InputType_t* u_type, const Vessel* vessel, const double delta, const int size);


MAP_ERROR_CODE fd_x_sequence(MAP_OtherStateType_t* other_type, MAP_InputType_t* u_type, MAP_OutputType_t* y_type, MAP_ConstraintStateType_t* z_type, Fd* force, const double epsilon, const int size, const double* original_pos, char* map_msg, MAP_ERROR_CODE* ierr);
MAP_ERROR_CODE fd_y_sequence(MAP_OtherStateType_t* other_type, MAP_InputType_t* u_type, MAP_OutputType_t* y_type, MAP_ConstraintStateType_t* z_type, Fd* force, const double epsilon, const int size, const double* original_pos, char* map_msg, MAP_ERROR_CODE* ierr);
MAP_ERROR_CODE fd_z_sequence(MAP_OtherStateType_t* other_type, MAP_InputType_t* u_type, MAP_OutputType_t* y_type, MAP_ConstraintStateType_t* z_type, Fd* force, const double epsilon, const int size, const double* original_pos, char* map_msg, MAP_ERROR_CODE* ierr);
MAP_ERROR_CODE fd_phi_sequence(MAP_OtherStateType_t* other_type, MAP_InputType_t* u_type, MAP_OutputType_t* y_type, MAP_ConstraintStateType_t* z_type, Fd* force, const double epsilon, const int size, const double* original_x, const double* original_y, const double* original_z, char* map_msg, MAP_ERROR_CODE* ierr);
MAP_ERROR_CODE fd_the_sequence(MAP_OtherStateType_t* other_type, MAP_InputType_t* u_type, MAP_OutputType_t* y_type, MAP_ConstraintStateType_t* z_type, Fd* force, const double epsilon, const int size, const double* original_x, const double* original_y, const double* original_z, char* map_msg, MAP_ERROR_CODE* ierr);
MAP_ERROR_CODE fd_psi_sequence(MAP_OtherStateType_t* other_type, MAP_InputType_t* u_type, MAP_OutputType_t* y_type, MAP_ConstraintStateType_t* z_type, Fd* force, const double epsilon, const int size, const double* original_x, const double* original_y, const double* original_z, char* map_msg, MAP_ERROR_CODE* ierr);
MAP_ERROR_CODE calculate_stiffness_2(double* K, Fd* force, const double delta, const int size);
MAP_ERROR_CODE set_moment_plus_2(const MAP_InputType_t* u_type, const MAP_OutputType_t* y_type, const Vessel* vessel, double* mx, double* my, double* mz, const int size);
MAP_ERROR_CODE set_moment_minus_2(const MAP_InputType_t* u_type, const MAP_OutputType_t* y_type, const Vessel* vessel, double* mx, double* my, double* mz, const int size);

/**
 * sets cable excursions (l and h) and reference frame psi rotation
 */
MAP_ERROR_CODE set_line_variables_pre_solve(ModelData* model_data, char* map_msg, MAP_ERROR_CODE* ierr);


/**
 * set: alpha, alpha at anchor, H, V and T at anchor, T at fairlead
 */
MAP_ERROR_CODE set_line_variables_post_solve(ModelData* model_data, char* map_msg, MAP_ERROR_CODE* ierr);

MapReal set_vertical_excursion(Element* element);
MapReal set_horizontal_excursion(Element* element);


/**
 *  'psi' is the angle between the element x-axis (local frame) and X-axis (global frame). This essentially produces this rotation matrix:
 *  
 *  \mathbf{R}(\psi) = \begin{bmatrix}
 *                     \cos\psi & -\sin\psi & 0 \\ 
 *                     \sin\psi &  \cos\psi & 0 \\ 
 *                            0 &         0 & 1
 *                     \end{bmatrix}
 *  
 *       1) first find psi - the angle of rotation between the element frame and the global reference frame
 *       2) r_j = fairlead displacement
 *       3) r_i = anchor displacement
 *       4) cos^{-1}( dot( (r_j-r_i) , (u_i) ) / ( norm(r_j-r_i) ) )
 */
MAP_ERROR_CODE set_psi(Element* element, char* map_msg, MAP_ERROR_CODE* ierr);

MAP_ERROR_CODE reset_node_force_to_zero(ModelData* model_data, char* map_msg, MAP_ERROR_CODE* ierr);


/**
 * Initialized omega (weight per unit length) and cross-section area of a cable. The formula is 
 *
 *   A=\pi*\frac{radius^{2}}{4}
 *   \omega=g*(\mu-A*\rho_{seawater})
 *
 * @acceses: none
 * @calledby: mapcall_msqs_init( )
 */                                                   
MAP_ERROR_CODE set_element_initial_guess(ModelData* model_data, char* map_msg, MAP_ERROR_CODE* ierr);


/**
 * MAP_InputType_t* uType,
 * MAP_ConstraintStateType_t* zType,
 * MAP_OtherStateType_t* otherType,
 * MAP_OutputType_t* yType,
 */
MAP_ERROR_CODE node_solve_sequence(ModelData* model_data, MAP_InputType_t* u_type, MAP_ConstraintStateType_t* z_type, MAP_OtherStateType_t* other_type, char* map_msg, MAP_ERROR_CODE* ierr);
MAP_ERROR_CODE line_solve_sequence(ModelData* model_data, double t, char* map_msg, MAP_ERROR_CODE* ierr);
MAP_ERROR_CODE solve_line(ModelData* model_data, double time, char* map_msg, MAP_ERROR_CODE* ierr);


/**
 * check to see if the line is double backing
 ELSEIF ( W  >  0.0_DbKi )  THEN   ! .TRUE. when the line will sink in fluid

         LMax      = XF - EA/W + SQRT( (EA/W)*(EA/W) + 2.0_DbKi*ZF*EA/W )  ! Compute the maximum stretched length of the line with seabed interaction beyond which the line would have to double-back on itself; here the line forms an "L" between the anchor and fairlead (i.e. it is horizontal along the seabed from the anchor, then vertical to the fairlead)

         IF ( ( L  >=  LMax   ) .AND. ( CB >= 0.0_DbKi ) )  &  ! .TRUE. if the line is as long or longer than its maximum possible value with seabed interaction
            CALL ProgAbort ( ' Unstretched mooring line length too large. '// &
                         ' Routine Catenary() cannot solve quasi-static mooring line solution.' )


      ENDIF
 */
MAP_ERROR_CODE check_maximum_line_length(Element* element, const bool contact_flag, char* map_msg, MAP_ERROR_CODE* ierr);

/**
 * call immediately after set_line_variables_post_solve(); this added H and V
 */
MAP_ERROR_CODE calculate_node_sum_force(ModelData* model_data);



/**
 * increment sum force value by (f) if node is fairlead; (-f) is node is anchor. 
 */
void add_to_sum_fx(Node* node, const MapReal fx);


/**
 * increment sum force value by (f) if node is fairlead; (-f) is node is anchor. 
 */
void add_to_sum_fy(Node* node, const MapReal fy);


/**
 * increment sum force value by (f) if node is fairlead; (-f) is node is anchor. 
 */
void add_to_sum_fz(Node* node, const MapReal fz);


#endif /* _LINEROUTINES_H */