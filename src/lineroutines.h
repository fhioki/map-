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


#ifndef _LINEROUTINES_H
#define _LINEROUTINES_H

#include "cminpack/cminpack.h"
#include "MAP_Types.h"


/**
 * sets cable excursions (l and h) and reference frame psi rotation
 */
MAP_ERROR_CODE set_line_variables_pre_solve(ModelData* model_data, char* map_msg, MAP_ERROR_CODE* ierr);


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

#endif /* _LINEROUTINES_H */
