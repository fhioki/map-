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

#include "map.h"
#include "initialization.h"


/**
 * @file 
 */


/** @addtogroup FortranCall */
/* @{ */
MAP_EXTERNCALL void set_sea_depth(MAP_ParameterType_t* p_type, const MapReal depth)
{
  p_type->depth = depth;
};


MAP_EXTERNCALL void set_sea_density(MAP_ParameterType_t* p_type, const MapReal rho)
{
  p_type->rhoSea = rho;
};


MAP_EXTERNCALL void set_gravity(MAP_ParameterType_t* p_type, const MapReal gravity)
{
  p_type->g = gravity;
};
/* @} */
