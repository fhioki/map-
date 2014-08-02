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
#include "maperror.h"
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


MAP_EXTERNCALL void set_gravity(MAP_ParameterType_t* p_type, double gravity)
{
  p_type->g = gravity;
};
/* @} */


void initialize_init_type_to_null(MAP_InitInputType_t* init_type)
{
  /* initialize the native Fortran/C types */
  init_type->gravity = -999.9;
  init_type->seaDensity = -999.9;
  init_type->depth = -999.9;
  init_type->fileName[0] = 0;
  init_type->summaryFileName[0] = 0;
  init_type->libraryInputLine[0] = 0;
  init_type->nodeInputLine[0] = 0;
  init_type->elementInputLine[0] = 0;
  init_type->optionInputLine[0] = 0;
};


void initialize_init_data_to_null(InitializationData* init_data)
{
  /* initialize the MAP initialization internal data strcture */
  init_data->libraryInputString = NULL;
  init_data->nodeInputString = NULL;
  init_data->elementInputString = NULL;
  init_data->solverOptionsString = NULL;
  init_data->expandedNodeInputString = NULL; 
  init_data->expandedElementInputString = NULL;
  init_data->summaryFileName = NULL;  
  init_data->sizeOfFullNodeString = 0; 
  init_data->sizeOfFullElementString = 0; 
  init_data->librarySize = 0;
  init_data->nodeSize = 0;
  init_data->elementSize = 0;
  init_data->solverOptionsSize = 0;
};
