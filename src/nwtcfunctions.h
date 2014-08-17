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


#ifndef _NWTC_FUNCTIONS_H
#define _NWTC_FUNCTIONS_H


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

#endif /* _NWTC_FUNCTIONS_H */
