/***************************************************************************
 *   Copyright (C) 2014 mdm                                                *
 *   marco.masciola@gmail.com                                              *
 *                                                                         *
 *   MAP++ is free software; you can redistribute it and/or modify it      *
 *   under the terms of the GNU General Public License as published by     *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.           *
 ***************************************************************************/


#ifndef _EXTERN_PROTOS_H
#define _EXTERN_PROTOS_H


/**
 * @brief     Initializes the MAP model and allocates memory. Inconsistencies with the input file is reported here. 
 *            This should only be called once. 
 * @param     initType initialization type, F2C FAST-native derived type
 * @param     uType input type, F2C FAST-native derived type. 
 * @param     pType parameter type, F2C FAST-native derived type
 * @param     xType continuous-state type, F2C FAST-native derived type
 * @param     none void type. It exists because MAP does not use discrete types (it would fill this place otherwise)
 * @param     zType constraint-state type, F2C FAST-native derived tpe
 * @param     otherType other-state type, F2C FAST-native derived type. This is core structure that houses all internal MAP states that cannot be mapped to Fortran
 * @param     yType output type, F2C FAST-native derived type
 * @param     ioType init-output types, F2C FAST-native derived type
 * @param     ierr error code
 * @param     map_msg error string
 */
MAP_EXTERNCALL void map_init(MAP_InitInputType_t* initType, 
                             MAP_InputType_t* uType,
                             MAP_ParameterType_t* pType,
                             MAP_ContinuousStateType_t* xType,
                             void* none,
                             MAP_ConstraintStateType_t* zType,
                             MAP_OtherStateType_t* otherType,
                             MAP_OutputType_t* yType,
                             MAP_InitOutputType_t* ioType,
                             MAP_ERROR_CODE* ierr,
                             char* map_msg );

/**
 * @brief     Solves the statics problem for the MSQS system and should be called at each time step, or after updating the vessel displacement. 
 *            Can be called multiple times, but must be called between {@link map_init()} and {@link map_end()}            
 * @param     t current (global) time
 * @param     interval coupling interval
 * @param     uType input type, F2C FAST-native derived type. 
 * @param     pType parameter type, F2C FAST-native derived type
 * @param     xType continuous-state type, F2C FAST-native derived type
 * @param     none void type. It exists because MAP does not use discrete types (it would fill this place otherwise)
 * @param     zType constraint-state type, F2C FAST-native derived tpe
 * @param     otherType other-state type, F2C FAST-native derived type. This is core structure that houses all internal MAP states that cannot be mapped to Fortran
 * @param     ierr error code
 * @param     map_msg error string
 * @see       {@link map_calc_output()}
 */
MAP_EXTERNCALL void map_update_states(double t,
                                      int interval,
                                      MAP_InputType_t* uType,
                                      MAP_ParameterType_t* pType,
                                      MAP_ContinuousStateType_t* xType,
                                      void* none,
                                      MAP_ConstraintStateType_t* zType,
                                      MAP_OtherStateType_t* otherType,
                                      MAP_ERROR_CODE* ierr,
                                      char* map_msg );

/**
 * @brief     Retrieves the values after the statics problem is solved. This function should be called immediately after map_update_states. 
 *            Can be called multiple times, but must be called between {@link map_init()} and {@link map_end()}
 * @param     t current (global) time
 * @param     uType input type, F2C FAST-native derived type. 
 * @param     pType parameter type, F2C FAST-native derived type
 * @param     xType continuous-state type, F2C FAST-native derived type
 * @param     none void type. It exists because MAP does not use discrete types (it would fill this place otherwise)
 * @param     zType constraint-state type, F2C FAST-native derived tpe
 * @param     otherType other-state type, F2C FAST-native derived type. This is core structure that houses all internal MAP states that cannot be mapped to Fortran
 * @param     yType output type, F2C FAST-native derived type
 * @param     ierr error code
 * @param     map_msg error string
 * @see       {@link map_update_states()}
 */
MAP_EXTERNCALL void map_calc_output(double t,
                                    MAP_InputType_t* uType,
                                    MAP_ParameterType_t* pType,
                                    MAP_ContinuousStateType_t* xType,
                                    void* none,
                                    MAP_ConstraintStateType_t* zType,
                                    MAP_OtherStateType_t* otherType,
                                    MAP_OutputType_t* yType,
                                    MAP_ERROR_CODE* ierr,
                                    char* map_msg);

/**
 * @brief     Deallocates all memory. Must be called after {@link map_init()}. This is called once. 
 * @param     uType input type, F2C FAST-native derived type. 
 * @param     pType parameter type, F2C FAST-native derived type
 * @param     xType continuous-state type, F2C FAST-native derived type
 * @param     none void type. It exists because MAP does not use discrete types (it would fill this place otherwise)
 * @param     zType constraint-state type, F2C FAST-native derived tpe
 * @param     otherType other-state type, F2C FAST-native derived type. This is core structure that houses all internal MAP states that cannot be mapped to Fortran
 * @param     yType output type, F2C FAST-native derived type
 * @param     ierr error code
 * @param     map_msg error string
 * @see       {@link map_update_states()}
 */
MAP_EXTERNCALL void map_end(MAP_InputType_t* uType,
                            MAP_ParameterType_t* pType,
                            MAP_ContinuousStateType_t* xType,
                            void* none,
                            MAP_ConstraintStateType_t* zType,
                            MAP_OtherStateType_t* otherType,
                            MAP_OutputType_t* yType,                                                                           
                            MAP_ERROR_CODE* ierr,
                            char* map_msg);

#endif /* _EXTERN_PROTOS_H */
