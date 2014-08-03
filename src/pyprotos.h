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
