/***************************************************************************
 *   Copyright (C) 2014 mdm                                                *
 *   marco[dot]masciola at gmail                                           *
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


#ifndef _INITIALIZATION_H
#define _INITIALIZATION_H

#include "MAP_Types.h"


/**
 * @brief Set the water depth. Should be called before {@link map_init()}
 * @param p_type paramter type, native C struct {@link MAP_ParameterType_t}
 * @param depth water depth [m]
 *
 * Example Fortran usage:
 * @code
 * INTERFACE                                                                
 *    SUBROUTINE mapextern_set_depth(interf, fc_val) BIND(C,name='set_sea_depth')  
 *      IMPORT                            
 *      IMPLICIT NONE                     
 *      TYPE(MAP_ParameterType_C) interf
 *      REAL(C_DOUBLE), VALUE :: fc_val  
 *    END SUBROUTINE mapextern_set_depth    
 * END INTERFACE                      
 *
 *   ...
 *
 * ! access the function using this subroutine call: 
 * CALL mapextern_set_depth(p%C_obj, depth)
 * @endcode
 */
MAP_EXTERNCALL void set_sea_depth(MAP_ParameterType_t* p_type, const MapReal depth);


/**
 * @brief Set the water density. Should be called before {@link map_init()}
 * @param p_type paramter type, native C struct {@link MAP_ParameterType_t}
 * @param rho water density [kg/m^3]
 *
 * Example Fortran usage:
 * @code
 * INTERFACE                                                                         
 *    SUBROUTINE mapextern_set_density(interf, fc_val) BIND(C,name='set_sea_density') 
 *      IMPORT                                            
 *      IMPLICIT NONE                                     
 *      TYPE(MAP_ParameterType_C) interf                
 *      REAL(C_DOUBLE), VALUE :: fc_val                      
 *    END SUBROUTINE mapextern_set_density                      
 * END INTERFACE                                          
 *
 *   ...
 *
 * ! access the function using this subroutine call: 
 * CALL mapextern_set_density(p%C_obj, rho)
 * @endcode
 */
MAP_EXTERNCALL void set_sea_density(MAP_ParameterType_t* p_type, const MapReal rho);


/**
 * @brief Set the gravitational constant. Should be called before {@link map_init()}
 * @param p_type paramter type, native C struct {@link MAP_ParameterType_t}
 * @param grtavity gravitational acceleration [m/s^2]
 *
 * Example Fortran usage:
 * @code
 * INTERFACE                                                               
 *    SUBROUTINE MAP_set_gravity(interf, fc_val) BIND(C,name='set_gravity')
 *      IMPORT                                                         
 *      IMPLICIT NONE                                                  
 *      TYPE(MAP_ParameterType_C) interf                               
 *      REAL(C_DOUBLE), VALUE :: fc_val                                
 *    END SUBROUTINE MAP_set_gravity                                   
 * END INTERFACE                                                       
 *
 *   ...
 *
 * ! access the function using this subroutine call: 
 * CALL mapextern_set_gravity(p%C_obj, g)
 * @endcode
 */
MAP_EXTERNCALL void set_gravity(MAP_ParameterType_t* p_type, const MapReal gravity);


/**
 * @brief Initializes the Fortran/C iteroperability types
 * @details This is called in the py_create_init_data routine following successful allocation 
 *          of memory. This should not be called directly by the user code. 
 * @param init_type, Fortran/C interoperable type {@link MAP_InitInputType_t}
 * @see py_create_init_data()
 */
void initialize_init_type_to_null(MAP_InitInputType_t* init_type);


/**
 * @brief Initializes MAP internal initialization data
 * @details This is called in the py_create_init_data routine following successful allocation 
 *          of memory. This should not be called directly by the user code. The internal states
 *          are nullified and set to zero. 
 * @param init_data, internal MAP initialization data {@link InitializationData}
 * @see py_create_init_data()
 */
void initialize_init_data_to_null(InitializationData* init_data);


#endif /* _INITIALIZATION_H */
