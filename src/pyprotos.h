/**
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


#ifndef _PYPROTOS_H
#define _PYPROTOS_H

#include "MAP_Types.h"

/**
 *
 */
MAP_EXTERNCALL MAP_InitInputType_t* py_create_init_data( char* msg, MAP_ERROR_CODE* status );
MAP_EXTERNCALL MAP_InitOutputType_t* py_create_initout_data( char* msg, MAP_ERROR_CODE* status );
MAP_EXTERNCALL MAP_InputType_t* py_create_input_data( char* msg, MAP_ERROR_CODE* status );
MAP_EXTERNCALL MAP_ParameterType_t* py_create_parameter_data( char* msg, MAP_ERROR_CODE* status );
MAP_EXTERNCALL MAP_ConstraintStateType_t* py_create_constraint_data( char* msg, MAP_ERROR_CODE* status );
MAP_EXTERNCALL MAP_OtherStateType_t* py_create_model_data( char* msg, MAP_ERROR_CODE* status );
MAP_EXTERNCALL MAP_OutputType_t* py_create_output_data( char* msg, MAP_ERROR_CODE* status );
MAP_EXTERNCALL MAP_ContinuousStateType_t* py_create_continuous_data( char* msg, MAP_ERROR_CODE* status );


/**
 */
MAP_EXTERNCALL int set_gravity( ModelData* data, MapReal g );
MAP_EXTERNCALL int set_sea_depth( ModelData* data, MapReal d );
MAP_EXTERNCALL int set_seaDensity( ModelData* data, MapReal rho );

/**
 * Frees data internal state data allcoated in the mapcall_msqs_init( ) function
 *
 * @acceses: none
 * @calledby: mapcall_msqs_init( )
 * @see: repeat_nodes( ), repeat_elements( )
 */
MAP_EXTERNCALL int free_init_data( InitializationData* init, char *map_msg, MAP_ERROR_CODE *ierr );


/**
 *
 */
MAP_EXTERNCALL MapReal get_residual_function_length( MAP_OtherStateType_t* otherFortType, int i, char *map_msg, MAP_ERROR_CODE *ierr );
MAP_EXTERNCALL MapReal get_residual_function_height( MAP_OtherStateType_t* otherFortType, int i, char *map_msg, MAP_ERROR_CODE *ierr );
MAP_EXTERNCALL MapReal pyget_jacobian_dxdh( MAP_OtherStateType_t* otherFortType, int i, char *map_msg, MAP_ERROR_CODE *ierr );
MAP_EXTERNCALL MapReal pyget_jacobian_dxdv( MAP_OtherStateType_t* otherFortType, int i, char *map_msg, MAP_ERROR_CODE *ierr );
MAP_EXTERNCALL MapReal pyget_jacobian_dzdh( MAP_OtherStateType_t* otherFortType, int i, char *map_msg, MAP_ERROR_CODE *ierr );
MAP_EXTERNCALL MapReal pyget_jacobian_dzdv( MAP_OtherStateType_t* otherFortType, int i, char *map_msg, MAP_ERROR_CODE *ierr );


/**
 *
 */
MAP_EXTERNCALL double* pyget_plot_x( MAP_OtherStateType_t* otherFortType, int i, int numPlotPoints, char *map_msg, MAP_ERROR_CODE *ierr );
MAP_EXTERNCALL double* pyget_plot_y( MAP_OtherStateType_t* otherFortType, int i, int numPlotPoints, char *map_msg, MAP_ERROR_CODE *ierr );
MAP_EXTERNCALL double* pyget_plot_z( MAP_OtherStateType_t* otherFortType, int i, int numPlotPoints, char *map_msg, MAP_ERROR_CODE *ierr );
MAP_EXTERNCALL void pyget_plot_array_free( MapReal *array ) ;


#endif /* _PYPROTOS_H */
