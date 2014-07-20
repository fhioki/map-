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

#include "map.h"
#include "maperror.h"
#include "protos.h"
#include "pyprotos.h"
#include "MAP_Types.h"

extern const char MAP_ERROR_STRING[][1024];


/* 
   @rm this is now handled by the fcall routines from python
  
   MAP_EXTERNCALL int set_gravity( ModelData* data, MapReal g )
  {  
    data->gravity = g;
    return MAP_SAFE;
  };
  
  
  MAP_EXTERNCALL int set_sea_depth( ModelData* data, MapReal d )
  {
    data->depth = d;
    return MAP_SAFE;
  };
  
  
  MAP_EXTERNCALL int set_sea_density( ModelData* data, MapReal rho )
  {
    data->seaDensity = rho;
    return MAP_SAFE;
  };
*/



