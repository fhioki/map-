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


#ifndef _NUMERIC_H
#define _NUMERIC_H

#include "map.h"
#include "maperror.h"
#include "cminpack/minpack.h"
#include "cminpack/cminpack.h"
// #include "MAP_Types.h"


MAP_ERROR_CODE call_minpack_lmder(Element* element, InnerSolveAttributes* inner_opt, ModelOptions* opt, const int line_num, const double time, char* map_msg, MAP_ERROR_CODE* ierr);
double get_maximum_line_length(Element* element);

int inner_function_evals(void* element_ptr, int m, int n, const __cminpack_real__* x, __cminpack_real__* fvec, __cminpack_real__* fjac, int ldfjac, int iflag);

MapReal residual_function_length_no_contact(const MapReal V, const MapReal H, const MapReal w, const MapReal Lu, const MapReal EA, const MapReal l);
MapReal residual_function_height_no_contact(const MapReal V, const MapReal H, const MapReal w, const MapReal Lu, const MapReal EA, const MapReal h);
MapReal residual_function_length_contact(const MapReal V, const MapReal H, const MapReal w, const MapReal Lu, const MapReal EA, const MapReal l, const MapReal cb);
MapReal residual_function_height_contact(const MapReal V, const MapReal H, const MapReal w, const MapReal Lu, const MapReal EA, const MapReal h, const MapReal cb);
MapReal jacobian_dxdh_no_contact(const MapReal V, const MapReal H, const MapReal w, const MapReal Lu, const MapReal EA);
MapReal jacobian_dxdv_no_contact(const MapReal V, const MapReal H, const MapReal w, const MapReal Lu, const MapReal EA);
MapReal jacobian_dzdh_no_contact(const MapReal V, const MapReal H, const MapReal w, const MapReal Lu, const MapReal EA);
MapReal jacobian_dzdv_no_contact(const MapReal V, const MapReal H, const MapReal w, const MapReal Lu, const MapReal EA);
MapReal jacobian_dxdh_contact(const MapReal V, const MapReal H, const MapReal w, const MapReal Lu, const MapReal EA, const MapReal cb);
MapReal jacobian_dxdv_contact(const MapReal V, const MapReal H, const MapReal w, const MapReal Lu, const MapReal EA, const MapReal cb);
MapReal jacobian_dzdh_contact(const MapReal V, const MapReal H, const MapReal w, const MapReal Lu, const MapReal EA, const MapReal cb);
MapReal jacobian_dzdv_contact(const MapReal V, const MapReal H, const MapReal w, const MapReal Lu, const MapReal EA, const MapReal cb);

MAP_ERROR_CODE forward_difference_jacobian(MAP_OtherStateType_t* other_type, MAP_ConstraintStateType_t* z_type, ModelData* model_data, char* map_msg, MAP_ERROR_CODE* ierr);
MAP_ERROR_CODE backward_difference_jacobian(MAP_OtherStateType_t* other_type, MAP_ConstraintStateType_t* z_type, ModelData* model_data, char* map_msg, MAP_ERROR_CODE* ierr);
MAP_ERROR_CODE central_difference_jacobian(MAP_OtherStateType_t* other_type, MAP_ConstraintStateType_t* z_type, ModelData* model_data, char* map_msg, MAP_ERROR_CODE* ierr);
MAP_ERROR_CODE lu(OuterSolveAttributes* ns, const int n, char* map_msg, MAP_ERROR_CODE* ierr);


/**
 * Ax = b -> LUx = b. Then y is defined to be Ux
 */
MAP_ERROR_CODE lu_back_substitution(OuterSolveAttributes* ns, const int n, char* map_msg, MAP_ERROR_CODE* ierr);


#endif // _NUMERIC_H
