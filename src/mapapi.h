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


#ifndef _MAPAPI_H
#define _MAPAPI_H


MAP_EXTERNCALL double* map_plot_x_array(MAP_OtherStateType_t* otherFortType, int i, int numPlotPoints, char* map_msg, MAP_ERROR_CODE* ierr);
MAP_EXTERNCALL double* map_plot_y_array(MAP_OtherStateType_t* otherFortType, int i, int numPlotPoints, char* map_msg, MAP_ERROR_CODE* ierr);
MAP_EXTERNCALL double* map_plot_z_array(MAP_OtherStateType_t* otherFortType, int i, int numPlotPoints, char* map_msg, MAP_ERROR_CODE* ierr);
MAP_EXTERNCALL void map_plot_array_free(MapReal * array) ;
MAP_EXTERNCALL MapReal map_residual_function_length(MAP_OtherStateType_t* otherFortType, int i, char* map_msg, MAP_ERROR_CODE* ierr);
MAP_EXTERNCALL MapReal map_residual_function_height(MAP_OtherStateType_t* otherFortType, int i, char* map_msg, MAP_ERROR_CODE* ierr);
MAP_EXTERNCALL MapReal map_jacobian_dxdh(MAP_OtherStateType_t* otherFortType, int i, char* map_msg, MAP_ERROR_CODE* ierr);
MAP_EXTERNCALL MapReal map_jacobian_dxdv(MAP_OtherStateType_t* otherFortType, int i, char* map_msg, MAP_ERROR_CODE* ierr);
MAP_EXTERNCALL MapReal map_jacobian_dzdh(MAP_OtherStateType_t* otherFortType, int i, char* map_msg, MAP_ERROR_CODE* ierr);
MAP_EXTERNCALL MapReal map_jacobian_dzdv(MAP_OtherStateType_t* otherFortType, int i, char* map_msg, MAP_ERROR_CODE* ierr);


/**
 * @brief   Returns vertical and horizontal fairlead force along element plane
 * @details 
 * @param   H, reference to horizontal fairlead force magnitude
 * @param   V, reference to vertical fairlead force magnitude
 * @param   other_type, other state type fortran derived
 * @param   index, element number we are requesting the data for
 * @param   map_msg, error message
 * @param   ierr, error code
 * @see     
 */
MAP_EXTERNCALL void map_get_fairlead_force_2d(double* H, double* V, MAP_OtherStateType_t* other_type, int index, char* map_msg, MAP_ERROR_CODE* ierr);


/**
 * @brief   Returns X, Y, and Z fairlead force in global reference frame
 * @details
 * @param   fx, reference to horizontal X fairlead force in global frame
 * @param   fy, reference to horizontal Y fairlead force in global frame
 * @param   fz, reference to vertical Z fairlead force in global frame
 * @param   other_type, other state type fortran derived
 * @param   index, element number we are requesting the data for
 * @param   map_msg, error message
 * @param   ierr, error code
 * @see     
 */
MAP_EXTERNCALL void map_get_fairlead_force_3d(double* fx, double* fy, double* fz, MAP_OtherStateType_t* other_type, int index, char* map_msg, MAP_ERROR_CODE* ierr);


/**
 * call this in python: offset_vessel().argtypes = [MapData_Type, MapInput_Type, c_char_p, POINTER(c_int)]        
 * angles are radians
 * 
 *     lib.map_offset_vessel.argtypes = [MapData_Type, MapInput_Type, c_double, c_double, c_double, c_double, c_double, c_double, c_char_p, POINTER(c_int)]        
 *
 * Angles are in degrees. This routine converts angles from deg to rad
 * Per Fossen (Fossen, Thor I. Guidance and control of ocean vehicles. Vol. 199. No. 4. New York: Wiley, 1994), this 
 * matrix converts vector from the body frame to the global reference frame:
 *
 * R = | cos(ψ)*cos(θ)    cos(ψ)*sin(θ)*sin(φ) − sin(ψ)*cos(φ)       cos(ψ)*sin(θ)*cos(φ) + sin(ψ)*sin(φ)  |
 *     | sin(ψ)*cos(θ)    sin(ψ)*sin(θ)*sin(φ) + cos(ψ)*cos(φ)       sin(ψ)*sin(θ)*cos(φ) − cos(ψ)*sin(φ)  |
 *     |   −sin(θ)                   cos(θ)*sin(φ)                                cos(θ)*cos(φ)            |
 *
 * We need to invoke this procedure to move the vessel nodes with body rotations factored:
 *
 * u_type = x + [R]*r       
 *   ▲     ▲       ▲
 *   |     |       |
 * global  |       |
 *      global     |
 *               local
 */
MAP_EXTERNCALL void map_offset_vessel(MAP_OtherStateType_t* other_type, MAP_InputType_t* u_type, double x, double y, double z, double phi, double the, double psi, char* map_msg, MAP_ERROR_CODE* ierr);


/**
 * lib.linearize_matrix.argtypes = [MapInput_Type, MapData_Type, MapOutnput_Type, c_double, c_char_p, POINTER(c_int)]        
 */
MAP_EXTERNCALL double** map_linearize_matrix(MAP_InputType_t* uType, MAP_OtherStateType_t* otherType, MAP_OutputType_t* yType, MAP_ConstraintStateType_t* zType, double epsilon, char* map_msg, MAP_ERROR_CODE* ierr);


/**
 * lib.py_free_linearize_matrix.argtypes = [POINTER(POINTER(c_double))]
 */
MAP_EXTERNCALL void map_free_linearize_matrix(double** array);

#endif /* _MAPAPI_H */
