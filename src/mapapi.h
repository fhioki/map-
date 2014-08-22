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
MAP_EXTERNCALL int map_size_elements(MAP_OtherStateType_t* other_type, MAP_ERROR_CODE* ierr, char* map_msg);

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
MAP_EXTERNCALL double** map_linearize_matrix(MAP_InputType_t* u_type, MAP_OtherStateType_t* other_type, MAP_OutputType_t* y_type, MAP_ConstraintStateType_t* z_type, double epsilon, MAP_ERROR_CODE* ierr, char* map_msg);


/**
 * lib.py_free_linearize_matrix.argtypes = [POINTER(POINTER(c_double))]
 */
MAP_EXTERNCALL void map_free_linearize_matrix(double** array);


/**
 * @brief     Initializes the MAP model and allocates memory. Inconsistencies with the input file is reported here. 
 *            This should only be called once. 
 * @param     init_type initialization type, F2C FAST-native derived type
 * @param     u_type input type, F2C FAST-native derived type. 
 * @param     p_type parameter type, F2C FAST-native derived type
 * @param     x_type continuous-state type, F2C FAST-native derived type
 * @param     none void type. It exists because MAP does not use discrete types (it would fill this place otherwise)
 * @param     z_type constraint-state type, F2C FAST-native derived tpe
 * @param     other_type other-state type, F2C FAST-native derived type. This is core structure that houses all internal MAP states that cannot be mapped to Fortran
 * @param     y_type output type, F2C FAST-native derived type
 * @param     ioType init-output types, F2C FAST-native derived type
 * @param     ierr error code
 * @param     map_msg error string
 */
MAP_EXTERNCALL void map_init(MAP_InitInputType_t* init_type, 
                             MAP_InputType_t* u_type,
                             MAP_ParameterType_t* p_type,
                             MAP_ContinuousStateType_t* x_type,
                             void* none,
                             MAP_ConstraintStateType_t* z_type,
                             MAP_OtherStateType_t* other_type,
                             MAP_OutputType_t* y_type,
                             MAP_InitOutputType_t* ioType,
                             MAP_ERROR_CODE* ierr,
                             char* map_msg);


/**
 * @brief     Solves the statics problem for the MSQS system and should be called at each time step, or after updating the vessel displacement. 
 *            Can be called multiple times, but must be called between {@link map_init()} and {@link map_end()}            
 * @param     t current (global) time
 * @param     interval coupling interval
 * @param     u_type input type, F2C FAST-native derived type. 
 * @param     p_type parameter type, F2C FAST-native derived type
 * @param     x_type continuous-state type, F2C FAST-native derived type
 * @param     none void type. It exists because MAP does not use discrete types (it would fill this place otherwise)
 * @param     z_type constraint-state type, F2C FAST-native derived tpe
 * @param     other_type other-state type, F2C FAST-native derived type. This is core structure that houses all internal MAP states that cannot be mapped to Fortran
 * @param     ierr error code
 * @param     map_msg error string
 * @see       {@link map_calc_output()}
 */
MAP_EXTERNCALL void map_update_states(double t,
                                      int interval,
                                      MAP_InputType_t* u_type,
                                      MAP_ParameterType_t* p_type,
                                      MAP_ContinuousStateType_t* x_type,
                                      void* none,
                                      MAP_ConstraintStateType_t* z_type,
                                      MAP_OtherStateType_t* other_type,
                                      MAP_ERROR_CODE* ierr,
                                      char* map_msg);


/**
 * @brief     Retrieves the values after the statics problem is solved. This function should be called immediately after map_update_states. 
 *            Can be called multiple times, but must be called between {@link map_init()} and {@link map_end()}
 * @param     t current (global) time
 * @param     u_type input type, F2C FAST-native derived type. 
 * @param     p_type parameter type, F2C FAST-native derived type
 * @param     x_type continuous-state type, F2C FAST-native derived type
 * @param     none void type. It exists because MAP does not use discrete types (it would fill this place otherwise)
 * @param     z_type constraint-state type, F2C FAST-native derived tpe
 * @param     other_type other-state type, F2C FAST-native derived type. This is core structure that houses all internal MAP states that cannot be mapped to Fortran
 * @param     y_type output type, F2C FAST-native derived type
 * @param     ierr error code
 * @param     map_msg error string
 * @see       {@link map_update_states()}
 */
MAP_EXTERNCALL void map_calc_output(double t,
                                    MAP_InputType_t* u_type,
                                    MAP_ParameterType_t* p_type,
                                    MAP_ContinuousStateType_t* x_type,
                                    void* none,
                                    MAP_ConstraintStateType_t* z_type,
                                    MAP_OtherStateType_t* other_type,
                                    MAP_OutputType_t* y_type,
                                    MAP_ERROR_CODE* ierr,
                                    char* map_msg);


/**
 * @brief     Deallocates all memory. Must be called after {@link map_init()}. This is called once. 
 * @param     u_type input type, F2C FAST-native derived type. 
 * @param     p_type parameter type, F2C FAST-native derived type
 * @param     x_type continuous-state type, F2C FAST-native derived type
 * @param     none void type. It exists because MAP does not use discrete types (it would fill this place otherwise)
 * @param     z_type constraint-state type, F2C FAST-native derived tpe
 * @param     other_type other-state type, F2C FAST-native derived type. This is core structure that houses all internal MAP states that cannot be mapped to Fortran
 * @param     y_type output type, F2C FAST-native derived type
 * @param     ierr error code
 * @param     map_msg error string
 * @see       {@link map_update_states()}
 */
MAP_EXTERNCALL void map_end(MAP_InputType_t* u_type,
                            MAP_ParameterType_t* p_type,
                            MAP_ContinuousStateType_t* x_type,
                            void* none,
                            MAP_ConstraintStateType_t* z_type,
                            MAP_OtherStateType_t* other_type,
                            MAP_OutputType_t* y_type,                                                                           
                            MAP_ERROR_CODE* ierr,
                            char* map_msg);

#endif /* _MAPAPI_H */
