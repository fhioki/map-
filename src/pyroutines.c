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


/**
 * MAP_EXTERNCALL InitializationData* MAP_InitInput_Create( char* map_msg, MAP_ERROR_CODE* ierr )
 * obj = Map.lib.MAP_InitInput_Create( self.status, pointer(self.ierr) )
 */
MAP_EXTERNCALL MAP_InitInputType_t* py_create_init_data(char* map_msg, MAP_ERROR_CODE* ierr)
{
  MAP_InitInputType_t* new_data = NULL;  
  *ierr = MAP_SAFE;
  new_data = (MAP_InitInputType_t*)malloc(sizeof(MAP_InitInputType_t));
  if (new_data==NULL) {
    *ierr = map_set_universal_error("", map_msg, ierr, MAP_FATAL_4);    
    return new_data;
  } else {    
    new_data->object = NULL;
    new_data->object = (InitializationData*)malloc(sizeof(InitializationData));
    if (new_data->object==NULL) {
      *ierr = map_set_universal_error("", map_msg, ierr, MAP_FATAL_4);          
    };
    return new_data;
  };
};


/**
 * MAP_EXTERNCALL InitializationOutputData* MAP_InitOutput_Create( char *map_msg, MAP_ERROR_CODE *ierr )
 * obj = Map.lib.MAP_InitOutput_Create( self.status, pointer(self.ierr) )
 */
MAP_EXTERNCALL MAP_InitOutputType_t* py_create_initout_data(char* map_msg, MAP_ERROR_CODE* ierr)
{
  MAP_InitOutputType_t* new_data = NULL;  
  *ierr = MAP_SAFE;
  new_data = (MAP_InitOutputType_t*)malloc(sizeof(MAP_InitOutputType_t));
  if (new_data==NULL) {
    *ierr = map_set_universal_error("", map_msg, ierr, MAP_FATAL_11);    
  }; 
  new_data->object = NULL;
  return new_data;    
};


/**
 * MAP_EXTERNCALL InputData* MAP_Input_Create( char* map_msg, MAP_ERROR_CODE *ierr) 
 * obj = Map.lib.MAP_Input_Create( self.status, pointer(self.ierr) )
 */
MAP_EXTERNCALL MAP_InputType_t* py_create_input_data( char* map_msg, MAP_ERROR_CODE* ierr )
{
  MAP_InputType_t* new_data=NULL;  
  *ierr = MAP_SAFE;
  new_data = (MAP_InputType_t*)malloc( sizeof(MAP_InputType_t) ); // @todo: needs to be freed
  if ( new_data == NULL ) 
  {
    *ierr = map_set_universal_error( "", map_msg, ierr, MAP_FATAL_5 );    
    return new_data;
  }; 
  new_data->object=NULL;
  return new_data;    
};


/**
 * MAP_EXTERNCALL ParameterData* MAP_Param_Create( char* map_msg, MAP_ERROR_CODE *ierr )
 * obj = Map.lib.MAP_Param_Create( self.status, pointer(self.ierr) )
 */
MAP_EXTERNCALL MAP_ParameterType_t* py_create_parameter_data(char* map_msg, MAP_ERROR_CODE* ierr)
{
  MAP_ParameterType_t* new_data = NULL;  
  *ierr = MAP_SAFE;
  new_data = (MAP_ParameterType_t*)malloc(sizeof(MAP_ParameterType_t));
  if (new_data==NULL) {
    *ierr = map_set_universal_error("", map_msg, ierr, MAP_FATAL_6);    
    return new_data;
  }; 
  new_data->object = NULL;
  return new_data;
};


/**
 * MAP_EXTERNCALL ConstraintData* MAP_ConstrState_Create( char* map_msg, MAP_ERROR_CODE *ierr )
 * obj = Map.lib.MAP_ConstrState_Create( self.status, pointer(self.ierr) )
 */
MAP_EXTERNCALL MAP_ConstraintStateType_t* py_create_constraint_data(char* map_msg, MAP_ERROR_CODE* ierr)
{
  MAP_ConstraintStateType_t* new_data = NULL;  
  *ierr = MAP_SAFE;
  new_data = (MAP_ConstraintStateType_t*)malloc(sizeof(MAP_ConstraintStateType_t)); 
  if (new_data==NULL) {
    *ierr = map_set_universal_error("", map_msg, ierr, MAP_FATAL_8);    
    return new_data;    
  }; 
  new_data->object = NULL;
  return new_data;      
};


/**
 * MAP_EXTERNCALL ModelData *MAP_OtherState_Create( char *map_msg, MAP_ERROR_CODE *ierr )
 * obj = Map.lib.MAP_OtherState_Create( self.status, pointer(self.ierr) )
 */
MAP_EXTERNCALL MAP_OtherStateType_t* py_create_model_data(char* map_msg, MAP_ERROR_CODE* ierr)
{
  MAP_OtherStateType_t* new_data = NULL;  
  *ierr = MAP_SAFE;
  new_data = (MAP_OtherStateType_t*)malloc(sizeof(MAP_OtherStateType_t)); 
  if (new_data==NULL) {
    *ierr = map_set_universal_error("", map_msg, ierr, MAP_FATAL_43);    
    return new_data;
  } 
  new_data->object = NULL;
  new_data->object = (ModelData*)malloc(sizeof(ModelData)); 
  if (new_data->object==NULL) {
    *ierr = map_set_universal_error("", map_msg, ierr, MAP_FATAL_43);            
  };
  return new_data;    
};


/**
 * MAP_EXTERNCALL OutputData *MAP_Output_Create( char *map_msg, MAP_ERROR_CODE *ierr )
 * obj = Map.lib.MAP_Output_Create( self.status, pointer(self.ierr) )
 */
MAP_EXTERNCALL MAP_OutputType_t* py_create_output_data(char* map_msg, MAP_ERROR_CODE* ierr)
{
  MAP_OutputType_t* new_data = NULL;  
  *ierr = MAP_SAFE;
  new_data = (MAP_OutputType_t*)malloc(sizeof(MAP_OutputType_t));
  if (new_data==NULL) {
    *ierr = map_set_universal_error("", map_msg, ierr, MAP_FATAL_10);    
    return new_data;
  };
  new_data->object = NULL;
  return new_data;    
};


/**
 * MAP_EXTERNCALL ContinuousData* MAP_ContState_Create( char* map_msg, MAP_ERROR_CODE *ierr )
 * obj = Map.lib.MAP_ContState_Create( self.status, pointer(self.ierr) )
 */
MAP_EXTERNCALL MAP_ContinuousStateType_t* py_create_continuous_data(char* map_msg, MAP_ERROR_CODE* ierr)
{
  MAP_ContinuousStateType_t* new_data = NULL;  
  *ierr = MAP_SAFE;
  new_data = (MAP_ContinuousStateType_t*)malloc(sizeof(MAP_ContinuousStateType_t)); 
  if (new_data==NULL) {
    *ierr = map_set_universal_error("", map_msg, ierr, MAP_FATAL_7);    
    return new_data;
  } 
  new_data->object = NULL;    
  return new_data;    
};


/**
 * call this in python: offset_vessel().argtypes = [MapData_Type, MapInput_Type, c_char_p, POINTER(c_int)]        
 * angles are radians
 * 
 *     lib.py_offset_vessel.argtypes = [MapData_Type, MapInput_Type, c_double, c_double, c_double, c_double, c_double, c_double, c_char_p, POINTER(c_int)]        
 *
 * Angles are in degrees. This routine converts angles from deg to rad
 */
MAP_EXTERNCALL void py_offset_vessel(MAP_OtherStateType_t* otherType, MAP_InputType_t* uType, double x, double y, double z, double phi, double the, double psi, char* map_msg, MAP_ERROR_CODE* ierr)
{
  ModelData* data = otherType->object;
  Vessel* vessel = &data->vessel;
  int i = 0;
  int N = uType->x_Len;
  double R[3][3];
  double cphi = 0.0;
  double sphi = 0.0;
  double cthe = 0.0;
  double sthe = 0.0;
  double cpsi = 0.0;
  double spsi = 0.0;
  double rx = 0.0;
  double ry = 0.0;
  double rz = 0.0;

  map_reset_universal_error(map_msg, ierr);
  
  /* define angles */
  cphi = cos(phi*DEG2RAD);
  sphi = sin(phi*DEG2RAD);
  cthe = cos(the*DEG2RAD);
  sthe = sin(the*DEG2RAD);
  cpsi = cos(psi*DEG2RAD);
  spsi = sin(psi*DEG2RAD);

  /*
    Per Fossen (Fossen, Thor I. Guidance and control of ocean vehicles. Vol. 199. No. 4. New York: Wiley, 1994), this 
    matrix converts vector from the body frame to the global reference frame:

    R = | cos(ψ)*cos(θ)    cos(ψ)*sin(θ)*sin(φ) − sin(ψ)*cos(φ)       cos(ψ)*sin(θ)*cos(φ) + sin(ψ)*sin(φ)  |
        | sin(ψ)*cos(θ)    sin(ψ)*sin(θ)*sin(φ) + cos(ψ)*cos(φ)       sin(ψ)*sin(θ)*cos(φ) − cos(ψ)*sin(φ)  |
        |   −sin(θ)                   cos(θ)*sin(φ)                                cos(θ)*cos(φ)            |

    We need to invoke this procedure to move the vessel nodes with body rotations factored:

    uType = x + [R]*r       
      ▲     ▲       ▲
      |     |       |
    global  |       |
         global     |
                  local
  */
  R[0][0] = cpsi*cthe;    R[0][1] = cpsi*sthe*sphi - spsi*cphi;   R[0][2] = cpsi*sthe*cphi + spsi*sphi;
  R[1][0] = sphi*cthe;    R[1][1] = sphi*sthe*sphi + cpsi*cphi;   R[1][2] = spsi*sthe*cphi - cpsi*sphi;
  R[2][0] = -sthe;        R[2][1] = cthe*sphi;                    R[2][2] = cthe*cphi;

  for (i=0 ; i<N ; i++) { 
    /* 
       @todo: need to include the reference position for non-zero reference origins , i.e. r = (xi-ref) 
       xi, yi, and zi are the original node position. We are adding the new displacement to it.
    */
    rx = vessel->xi[i];
    ry = vessel->yi[i];
    rz = vessel->zi[i];

    /* matrix-vector product */
    uType->x[i] = x + rx*R[0][0] + ry*R[0][1] + rz*R[0][2];
    uType->y[i] = y + rx*R[1][0] + ry*R[1][1] + rz*R[1][2];
    uType->z[i] = z + rx*R[2][0] + ry*R[2][1] + rz*R[2][2];
  };
};



/**
 *
 */
MAP_ERROR_CODE fd_x_sequence(MAP_OtherStateType_t* otherType, MAP_InputType_t* uType, MAP_OutputType_t* yType, MAP_ConstraintStateType_t* zType, Fd* force, const double epsilon, const int size, const double* originalPos, char* map_msg, MAP_ERROR_CODE* ierr)
{
  MAP_ERROR_CODE success = MAP_SAFE;
  ModelData* data = otherType->object;
  Vessel* vessel = &data->vessel;

  do {
    /* minus epsilon sequence */
    success = increment_dof_by_delta(uType->x, -epsilon, size); CHECKERRQ(MAP_FATAL_61);
    if (data->MAP_SOLVE_TYPE==MONOLITHIC) {
      success = line_solve_sequence(data, 0.0, map_msg, ierr);
    } else {
      success = node_solve_sequence(data, uType, zType, otherType, map_msg, ierr); // @todo CHECKERRQ()
    };    
    success = set_force_plus(yType->Fx, force->fx, size); CHECKERRQ(MAP_FATAL_61);
    success = set_force_plus(yType->Fy, force->fy, size); CHECKERRQ(MAP_FATAL_61);
    success = set_force_plus(yType->Fz, force->fz, size); CHECKERRQ(MAP_FATAL_61);
    success = set_moment_plus(yType, vessel, force->mx, force->my, force->mz, size); CHECKERRQ(MAP_FATAL_61);
    success = restore_original_displacement(uType->x, originalPos, size); CHECKERRQ(MAP_FATAL_61);
  
    /* plus epsilon sequence */
    success = increment_dof_by_delta(uType->x, epsilon, size); CHECKERRQ(MAP_FATAL_61);
    if (data->MAP_SOLVE_TYPE==MONOLITHIC) {
      success = line_solve_sequence(data, 0.0, map_msg, ierr);
    } else {
      success = node_solve_sequence(data, uType, zType, otherType, map_msg, ierr); // @todo CHECKERRQ()
    };    
    success = set_force_minus(yType->Fx, force->fx, size); CHECKERRQ(MAP_FATAL_61);
    success = set_force_minus(yType->Fy, force->fy, size); CHECKERRQ(MAP_FATAL_61);
    success = set_force_minus(yType->Fz, force->fz, size); CHECKERRQ(MAP_FATAL_61);
    success = set_moment_minus(yType, vessel, force->mx, force->my, force->mz, size); CHECKERRQ(MAP_FATAL_61);
    success = restore_original_displacement(uType->x, originalPos, size); CHECKERRQ(MAP_FATAL_61);
  } while (0);
  MAP_RETURN;
};


/**
 *
 */
MAP_ERROR_CODE fd_y_sequence(MAP_OtherStateType_t* otherType, MAP_InputType_t* uType, MAP_OutputType_t* yType, MAP_ConstraintStateType_t* zType, Fd* force, const double epsilon, const int size, const double* originalPos, char* map_msg, MAP_ERROR_CODE* ierr)
{
  MAP_ERROR_CODE success = MAP_SAFE;
  ModelData* data = otherType->object;
  Vessel* vessel = &data->vessel;

  do {
    /* minus epsilon sequence */
    success = increment_dof_by_delta(uType->y, -epsilon, size); CHECKERRQ(MAP_FATAL_61);
    if (data->MAP_SOLVE_TYPE==MONOLITHIC) {
      success = line_solve_sequence(data, 0.0, map_msg, ierr);
    } else {
      success = node_solve_sequence(data, uType, zType, otherType, map_msg, ierr); // @todo CHECKERRQ()
    };    
    success = set_force_plus(yType->Fx, force->fx, size); CHECKERRQ(MAP_FATAL_61);
    success = set_force_plus(yType->Fy, force->fy, size); CHECKERRQ(MAP_FATAL_61);
    success = set_force_plus(yType->Fz, force->fz, size); CHECKERRQ(MAP_FATAL_61);
    success = set_moment_plus(yType, vessel, force->mx, force->my, force->mz, size); CHECKERRQ(MAP_FATAL_61);
    success = restore_original_displacement(uType->y, originalPos, size); CHECKERRQ(MAP_FATAL_61);
        
    /* plus epsilon sequence */
    success = increment_dof_by_delta(uType->y, epsilon, size); CHECKERRQ(MAP_FATAL_61);
    if (data->MAP_SOLVE_TYPE==MONOLITHIC) {
      success = line_solve_sequence(data, 0.0, map_msg, ierr);
    } else {
      success = node_solve_sequence(data, uType, zType, otherType, map_msg, ierr); // @todo CHECKERRQ()
    };    
    success = set_force_minus(yType->Fx, force->fx, size); CHECKERRQ(MAP_FATAL_61);
    success = set_force_minus(yType->Fy, force->fy, size); CHECKERRQ(MAP_FATAL_61);
    success = set_force_minus(yType->Fz, force->fz, size); CHECKERRQ(MAP_FATAL_61);
    success = set_moment_minus(yType, vessel, force->mx, force->my, force->mz, size); CHECKERRQ(MAP_FATAL_61);
    success = restore_original_displacement(uType->y, originalPos, size); CHECKERRQ(MAP_FATAL_61);
  } while (0);
  
  MAP_RETURN;
};


/**
 *
 */
MAP_ERROR_CODE fd_z_sequence(MAP_OtherStateType_t* otherType, MAP_InputType_t* uType, MAP_OutputType_t* yType, MAP_ConstraintStateType_t* zType, Fd* force, const double epsilon, const int size, const double* originalPos, char* map_msg, MAP_ERROR_CODE* ierr)
{
  MAP_ERROR_CODE success = MAP_SAFE;
  ModelData* data = otherType->object;
  Vessel* vessel = &data->vessel;

  do {
    /* minus epsilon sequence */
    success = increment_dof_by_delta(uType->z, -epsilon, size); CHECKERRQ(MAP_FATAL_61);
    if (data->MAP_SOLVE_TYPE==MONOLITHIC) {
      success = line_solve_sequence(data, 0.0, map_msg, ierr);
    } else {
      success = node_solve_sequence(data, uType, zType, otherType, map_msg, ierr); // @todo CHECKERRQ()
    };    
    success = set_force_plus(yType->Fx, force->fx, size); CHECKERRQ(MAP_FATAL_61);
    success = set_force_plus(yType->Fy, force->fy, size); CHECKERRQ(MAP_FATAL_61);
    success = set_force_plus(yType->Fz, force->fz, size); CHECKERRQ(MAP_FATAL_61);
    success = set_moment_plus(yType, vessel, force->mx, force->my, force->mz, size); CHECKERRQ(MAP_FATAL_61);
    success = restore_original_displacement(uType->z, originalPos, size); CHECKERRQ(MAP_FATAL_61);
        
    /* plus epsilon sequence */
    success = increment_dof_by_delta(uType->z, epsilon, size); CHECKERRQ(MAP_FATAL_61);
    if (data->MAP_SOLVE_TYPE==MONOLITHIC) {
      success = line_solve_sequence(data, 0.0, map_msg, ierr);
    } else {
      success = node_solve_sequence(data, uType, zType, otherType, map_msg, ierr); // @todo CHECKERRQ()
    };    
    success = set_force_minus(yType->Fx, force->fx, size); CHECKERRQ(MAP_FATAL_61);
    success = set_force_minus(yType->Fy, force->fy, size); CHECKERRQ(MAP_FATAL_61);
    success = set_force_minus(yType->Fz, force->fz, size); CHECKERRQ(MAP_FATAL_61);
    success = set_moment_minus(yType, vessel, force->mx, force->my, force->mz, size); CHECKERRQ(MAP_FATAL_61);
    success = restore_original_displacement(uType->z, originalPos, size); CHECKERRQ(MAP_FATAL_61);
  } while (0);

  MAP_RETURN;
};


/**
 *
 */
MAP_ERROR_CODE fd_phi_sequence(MAP_OtherStateType_t* otherType, MAP_InputType_t* uType, MAP_OutputType_t* yType, MAP_ConstraintStateType_t* zType, Fd* force, const double epsilon, const int size, const double* originalX, const double* originalY, const double* originalZ, char* map_msg, MAP_ERROR_CODE* ierr)
{
  MAP_ERROR_CODE success = MAP_SAFE;
  ModelData* data = otherType->object;
  Vessel* vessel = &data->vessel;

  do {
    /* minus epsilon sequence */
    success = increment_phi_dof_by_delta(uType, vessel, -epsilon, size); CHECKERRQ(MAP_FATAL_61);        
    if (data->MAP_SOLVE_TYPE==MONOLITHIC) {
       success = line_solve_sequence(data, 0.0, map_msg, ierr);
    } else {
      success = node_solve_sequence(data, uType, zType, otherType, map_msg, ierr); // @todo CHECKERRQ()
    };    
    success = set_force_plus(yType->Fx, force->fx, size); CHECKERRQ(MAP_FATAL_61);
    success = set_force_plus(yType->Fy, force->fy, size); CHECKERRQ(MAP_FATAL_61);
    success = set_force_plus(yType->Fz, force->fz, size); CHECKERRQ(MAP_FATAL_61);
    success = set_moment_plus_2(uType, yType, vessel, force->mx, force->my, force->mz, size); CHECKERRQ(MAP_FATAL_61);
    success = restore_original_displacement(uType->x, originalX, size); CHECKERRQ(MAP_FATAL_61);        
    success = restore_original_displacement(uType->y, originalY, size); CHECKERRQ(MAP_FATAL_61);        
    success = restore_original_displacement(uType->z, originalZ, size); CHECKERRQ(MAP_FATAL_61);        
    
    /* plus epsilon sequence */
    success = increment_phi_dof_by_delta(uType, vessel, epsilon, size); CHECKERRQ(MAP_FATAL_61);        
    if (data->MAP_SOLVE_TYPE==MONOLITHIC) {
      success = line_solve_sequence(data, 0.0, map_msg, ierr);
    } else {
      success = node_solve_sequence(data, uType, zType, otherType, map_msg, ierr); // @todo CHECKERRQ()
    };    
    success = set_force_minus(yType->Fx, force->fx, size); CHECKERRQ(MAP_FATAL_61);
    success = set_force_minus(yType->Fy, force->fy, size); CHECKERRQ(MAP_FATAL_61);
    success = set_force_minus(yType->Fz, force->fz, size); CHECKERRQ(MAP_FATAL_61);
    success = set_moment_minus_2(uType, yType, vessel, force->mx, force->my, force->mz, size); CHECKERRQ(MAP_FATAL_61);
    success = restore_original_displacement(uType->x, originalX, size); CHECKERRQ(MAP_FATAL_61);        
    success = restore_original_displacement(uType->y, originalY, size); CHECKERRQ(MAP_FATAL_61);        
    success = restore_original_displacement(uType->z, originalZ, size); CHECKERRQ(MAP_FATAL_61);                
  } while (0);
  
  MAP_RETURN;
};


/**
 *
 */
MAP_ERROR_CODE fd_the_sequence(MAP_OtherStateType_t* otherType, MAP_InputType_t* uType, MAP_OutputType_t* yType, MAP_ConstraintStateType_t* zType, Fd* force, const double epsilon, const int size, const double* originalX, const double* originalY, const double* originalZ, char* map_msg, MAP_ERROR_CODE* ierr)
{
  MAP_ERROR_CODE success = MAP_SAFE;
  ModelData* data = otherType->object;
  Vessel* vessel = &data->vessel;

  do {
    /* minus epsilon sequence */
    success = increment_the_dof_by_delta(uType, vessel, -epsilon, size); CHECKERRQ(MAP_FATAL_61);
    if (data->MAP_SOLVE_TYPE==MONOLITHIC) {
      success = line_solve_sequence(data, 0.0, map_msg, ierr);
    } else {
      success = node_solve_sequence(data, uType, zType, otherType, map_msg, ierr); // @todo CHECKERRQ()
    };    
    success = set_force_plus(yType->Fx, force->fx, size); CHECKERRQ(MAP_FATAL_61);
    success = set_force_plus(yType->Fy, force->fy, size); CHECKERRQ(MAP_FATAL_61);
    success = set_force_plus(yType->Fz, force->fz, size); CHECKERRQ(MAP_FATAL_61);
    success = set_moment_plus_2(uType, yType, vessel, force->mx, force->my, force->mz, size); CHECKERRQ(MAP_FATAL_61);
    success = restore_original_displacement(uType->x, originalX, size); CHECKERRQ(MAP_FATAL_61);        
    success = restore_original_displacement(uType->y, originalY, size); CHECKERRQ(MAP_FATAL_61);        
    success = restore_original_displacement(uType->z, originalZ, size); CHECKERRQ(MAP_FATAL_61);        
    
    /* plut epsilon sequence */
    success = increment_the_dof_by_delta(uType, vessel, epsilon, size); CHECKERRQ(MAP_FATAL_61);        
    if (data->MAP_SOLVE_TYPE==MONOLITHIC) {
      success = line_solve_sequence(data, 0.0, map_msg, ierr);
    } else {
      success = node_solve_sequence(data, uType, zType, otherType, map_msg, ierr); // @todo CHECKERRQ()
    };    
    success = set_force_minus(yType->Fx, force->fx, size); CHECKERRQ(MAP_FATAL_61);
    success = set_force_minus(yType->Fy, force->fy, size); CHECKERRQ(MAP_FATAL_61);
    success = set_force_minus(yType->Fz, force->fz, size); CHECKERRQ(MAP_FATAL_61);
    success = set_moment_minus_2(uType, yType, vessel, force->mx, force->my, force->mz, size); CHECKERRQ(MAP_FATAL_61);
    success = restore_original_displacement(uType->x, originalX, size); CHECKERRQ(MAP_FATAL_61);        
    success = restore_original_displacement(uType->y, originalY, size); CHECKERRQ(MAP_FATAL_61);        
    success = restore_original_displacement(uType->z, originalZ, size); CHECKERRQ(MAP_FATAL_61);        
  } while (0);

  MAP_RETURN;
};


/**
 *
 */
MAP_ERROR_CODE fd_psi_sequence(MAP_OtherStateType_t* otherType, MAP_InputType_t* uType, MAP_OutputType_t* yType, MAP_ConstraintStateType_t* zType, Fd* force, const double epsilon, const int size, const double* originalX, const double* originalY, const double* originalZ, char* map_msg, MAP_ERROR_CODE* ierr)
{
  MAP_ERROR_CODE success = MAP_SAFE;
  ModelData* data = otherType->object;
  Vessel* vessel = &data->vessel;
  do {
    /* minus epsilon sequence */
    success = increment_psi_dof_by_delta(uType, vessel, -epsilon, size); CHECKERRQ(MAP_FATAL_61);
    if (data->MAP_SOLVE_TYPE==MONOLITHIC) {
      success = line_solve_sequence(data, 0.0, map_msg, ierr);
    } else {
      success = node_solve_sequence(data, uType, zType, otherType, map_msg, ierr); // @todo CHECKERRQ()
    };    
    success = set_force_plus(yType->Fx, force->fx, size); CHECKERRQ(MAP_FATAL_61);
    success = set_force_plus(yType->Fy, force->fy, size); CHECKERRQ(MAP_FATAL_61);
    success = set_force_plus(yType->Fz, force->fz, size); CHECKERRQ(MAP_FATAL_61);
    success = set_moment_plus_2(uType, yType, vessel, force->mx, force->my, force->mz, size); CHECKERRQ(MAP_FATAL_61);
    success = restore_original_displacement(uType->x, originalX, size); CHECKERRQ(MAP_FATAL_61);
    success = restore_original_displacement(uType->y, originalY, size); CHECKERRQ(MAP_FATAL_61);
    success = restore_original_displacement(uType->z, originalZ, size); CHECKERRQ(MAP_FATAL_61);
        
    /* plut epsilon sequence */
    success = increment_psi_dof_by_delta(uType, vessel, epsilon, size); CHECKERRQ(MAP_FATAL_61);
    if (data->MAP_SOLVE_TYPE==MONOLITHIC) {
      success = line_solve_sequence(data, 0.0, map_msg, ierr);
    } else {
      success = node_solve_sequence(data, uType, zType, otherType, map_msg, ierr); // @todo CHECKERRQ()
    };    
    success = set_force_minus(yType->Fx, force->fx, size); CHECKERRQ(MAP_FATAL_61);
    success = set_force_minus(yType->Fy, force->fy, size); CHECKERRQ(MAP_FATAL_61);
    success = set_force_minus(yType->Fz, force->fz, size); CHECKERRQ(MAP_FATAL_61);
    success = set_moment_minus_2(uType, yType, vessel, force->mx, force->my, force->mz, size); CHECKERRQ(MAP_FATAL_61);
    success = restore_original_displacement(uType->x, originalX, size); CHECKERRQ(MAP_FATAL_61);
    success = restore_original_displacement(uType->y, originalY, size); CHECKERRQ(MAP_FATAL_61);
    success = restore_original_displacement(uType->z, originalZ, size); CHECKERRQ(MAP_FATAL_61);
  } while (0);

  MAP_RETURN;
};


/**
 * success = calculate_stiffness(stiffness, fx, epsilon, N);
 */
MAP_ERROR_CODE calculate_stiffness_2(double* K, Fd* force, const double delta, const int size)
{  
  int i = 0;
  for (i=0 ; i<size ; i++) {
    K[0] += (force->fx[i]/(2*delta));
    K[1] += (force->fy[i]/(2*delta));
    K[2] += (force->fz[i]/(2*delta));
    K[3] += (force->mx[i]/(2*delta));
    K[4] += (force->my[i]/(2*delta));
    K[5] += (force->mz[i]/(2*delta));
  };
  return MAP_SAFE;
};


/**
 * lib.linearize_matrix.argtypes = [MapInput_Type, MapData_Type, MapOutnput_Type, c_double, c_char_p, POINTER(c_int)]        
 */
MAP_EXTERNCALL double** py_linearize_matrix(MAP_InputType_t* uType, MAP_OtherStateType_t* otherType, MAP_OutputType_t* yType, MAP_ConstraintStateType_t* zType, double epsilon, char* map_msg, MAP_ERROR_CODE* ierr)
{
  MapReal* xOriginal = NULL;
  MapReal* yOriginal = NULL;
  MapReal* zOriginal = NULL;
  ModelData* data = otherType->object;
  MAP_ERROR_CODE success = MAP_SAFE;
  const int N = uType->x_Len;
  const int SIX = 6;
  int i = 0;
  int j = 0;
  int k = 0;
  Fd force;
  double** K;

  map_reset_universal_error(map_msg, ierr);
  K = (double**)malloc(SIX*sizeof(double*));
  for (i=0 ; i<SIX ; i++) {
    K[i] = (double*)malloc(SIX*sizeof(double));

    /* initialize K(6x6) allocated above to zero, row-by-row */
    K[i][0] = 0.0;
    K[i][1] = 0.0;
    K[i][2] = 0.0;
    K[i][3] = 0.0;
    K[i][4] = 0.0;
    K[i][5] = 0.0;
  };
  
  force.fx = (double*)malloc(N*sizeof(double));
  force.fy = (double*)malloc(N*sizeof(double));
  force.fz = (double*)malloc(N*sizeof(double));
  force.mx = (double*)malloc(N*sizeof(double));
  force.my = (double*)malloc(N*sizeof(double));
  force.mz = (double*)malloc(N*sizeof(double));  
  xOriginal = (double*)malloc(N*sizeof(double));
  yOriginal = (double*)malloc(N*sizeof(double));
  zOriginal = (double*)malloc(N*sizeof(double));
  
  /* initialize stuff allocated above to zero */
  for (i=0 ; i<N ; i++) {
    force.fx[i] = 0.0;
    force.fy[i] = 0.0;
    force.fz[i] = 0.0;
    force.mx[i] = 0.0;
    force.my[i] = 0.0;
    force.mz[i] = 0.0;
    xOriginal[i] = 0.0;
    yOriginal[i] = 0.0;
    zOriginal[i] = 0.0;
  };
   
  do {    
    /* first get the original values for the displacements */
    for (k=0 ; k<N ; k++) {
      xOriginal[k] = uType->x[k];
      yOriginal[k] = uType->y[k];
      zOriginal[k] = uType->z[k];      
    };
   
    for (i=0 ; i<SIX ; i++) { /* down, force direction changes */
      success = reset_force_to_zero(force.fx, force.fy, force.fz, force.mx, force.my, force.mz, N);
      if (i==0) {        
        success = fd_x_sequence(otherType, uType, yType, zType, &force, epsilon, N, xOriginal, map_msg, ierr); CHECKERRQ(MAP_FATAL_62);
        success = calculate_stiffness_2(K[0], &force, epsilon, N); CHECKERRQ(MAP_FATAL_62);
        // printf("This is the stiffness: %12.2f, %12.2f, %12.2f, %12.2f, %12.2f, %12.2f\n", K[0][0], K[0][1], K[0][2], K[0][3], K[0][4], K[0][5]);
      } else if (i==1) {
        success = fd_y_sequence(otherType, uType, yType, zType, &force, epsilon, N, yOriginal, map_msg, ierr); CHECKERRQ(MAP_FATAL_63);
        success = calculate_stiffness_2(K[1], &force, epsilon, N); CHECKERRQ(MAP_FATAL_63);
        // printf("This is the stiffness: %12.2f, %12.2f, %12.2f, %12.2f, %12.2f, %12.2f\n", K[1][0], K[1][1], K[1][2], K[1][3], K[1][4], K[1][5]);
      } else if (i==2) {
        success = fd_z_sequence(otherType, uType, yType, zType, &force, epsilon, N, zOriginal, map_msg, ierr); CHECKERRQ(MAP_FATAL_64);
        success = calculate_stiffness_2(K[2], &force, epsilon, N); CHECKERRQ(MAP_FATAL_64);
        // printf("This is the stiffness: %12.2f, %12.2f, %12.2f, %12.2f, %12.2f, %12.2f\n", K[2][0], K[2][1], K[2][2], K[2][3], K[2][4], K[2][5]);
      } else if (i==3) {
        success = fd_phi_sequence(otherType, uType, yType, zType, &force, epsilon, N, xOriginal, yOriginal, zOriginal, map_msg, ierr); //CHECKERRQ(MAP_FATAL_65);
        success = calculate_stiffness_2(K[3], &force, epsilon, N); CHECKERRQ(MAP_FATAL_65);
        // printf("This is the stiffness: %12.2f, %12.2f, %12.2f, %12.2f, %12.2f, %12.2f\n", K[3][0], K[3][1], K[3][2], K[3][3], K[3][4], K[3][5]);
      } else if (i==4) {
        success = fd_the_sequence(otherType, uType, yType, zType, &force, epsilon, N, xOriginal, yOriginal, zOriginal, map_msg, ierr); //CHECKERRQ(MAP_FATAL_66);
        success = calculate_stiffness_2(K[4], &force, epsilon, N); CHECKERRQ(MAP_FATAL_66);
        // printf("This is the stiffness: %12.2f, %12.2f, %12.2f, %12.2f, %12.2f, %12.2f\n", K[4][0], K[4][1], K[4][2], K[4][3], K[4][4], K[4][5]);
      } else if (i==5) {
        success = fd_psi_sequence(otherType, uType, yType, zType, &force, epsilon, N, xOriginal, yOriginal, zOriginal, map_msg, ierr); //CHECKERRQ(MAP_FATAL_67);
        success = calculate_stiffness_2(K[5], &force, epsilon, N); CHECKERRQ(MAP_FATAL_67);
        // printf("This is the stiffness: %12.2f, %12.2f, %12.2f, %12.2f, %12.2f, %12.2f\n", K[5][0], K[5][1], K[5][2], K[5][3], K[5][4], K[5][5]);
      };
    };
  } while (0);  

  success = reset_force_to_zero(force.fx, force.fy, force.fz, force.mx, force.my, force.mz, N);
  success = restore_original_displacement(uType->x, xOriginal, N);
  success = restore_original_displacement(uType->y, yOriginal, N);
  success = restore_original_displacement(uType->z, zOriginal, N);
  success = line_solve_sequence(data, 0.0, map_msg, ierr); 
  
  MAPFREE(force.fx);
  MAPFREE(force.fy);
  MAPFREE(force.fz);
  MAPFREE(force.mx);
  MAPFREE(force.my);
  MAPFREE(force.mz);
  MAPFREE(xOriginal);
  MAPFREE(yOriginal);
  MAPFREE(zOriginal);  
  return K;
};



/**
 * lib.py_free_linearize_matrix.argtypes = [POINTER(POINTER(c_double))]
 */
MAP_EXTERNCALL void py_free_linearize_matrix(double** array)
{
  int i = 0;
  for(i=0 ; i<6 ; i++) {
    MAPFREE(array[i]);
  };
  MAPFREE(array);
};


/**
 * success = reset_force_to_zero(yType, N);
 */
MAP_ERROR_CODE reset_force_to_zero(double* fx, double* fy, double* fz, double* mx, double* my, double* mz, const int size)
{
  int i = 0;
  for (i=0 ; i<size ; i++) {
    fx[i] = 0.0;
    fy[i] = 0.0;
    fz[i] = 0.0;
    mx[i] = 0.0;
    my[i] = 0.0;
    mz[i] = 0.0;
  };
  return MAP_SAFE;
};


/**
 * success = restore_original_displacement(uType->x, xOriginal, N);
 */
MAP_ERROR_CODE restore_original_displacement(double* inputType, const double* initialValue, const int size)
{
  int i = 0;
  for (i=0 ; i<size ; i++) {
    inputType[i] = initialValue[i];
  };
  return MAP_SAFE;
};


/**
 * success = incremenet_x_dof_by_delta(uType, -epsilon);
 */
MAP_ERROR_CODE increment_dof_by_delta(double* inputType, const MapReal delta, const int size)
{
  int i = 0;
  for (i=0 ; i<size ; i++) {
    inputType[i] += delta;
  };
  return MAP_SAFE;
};


/**
 * success = set_force_plus(yType->Fx, fx, N);
 */
MAP_ERROR_CODE set_force_plus(const double* inputType, double* force, const int size)
{
  int i = 0;
  for (i=0 ; i<size ; i++) {
    force[i] += -inputType[i];
  };
  return MAP_SAFE;
};


/**
 * success = set_force_plus(yType->Fx, fx, N);
 */
MAP_ERROR_CODE set_force_minus(const double* inputType, double* force, const int size)
{
  int i = 0;
  for (i=0 ; i<size ; i++) {
    force[i] -= -inputType[i];
  };
  return MAP_SAFE;
};


// /**
//  * success = calculate_stiffness(stiffness, fx, epsilon, N);
//  */
// MAP_ERROR_CODE calculate_stiffness(double* K, double* force, const double delta, const int size)
// {  
//   *K = 0.0;
//   int i = 0;
//   for (i=0 ; i<size ; i++) {
//     *K += (force[i]/(2*delta));
//   };
//   return MAP_SAFE;
// };


/**
 * success = set_moment_minus(yType, vessel, mx, my, mz, N)
 */
MAP_ERROR_CODE set_moment_minus(const MAP_OutputType_t* outputType, const Vessel* vessel, double* mx, double* my, double* mz, const int size)
{
  int i = 0;
  double rx = 0.0;
  double ry = 0.0;
  double rz = 0.0;

  for (i=0 ; i<size ; i++) {    
    /* 
       @todo: this is not consistent with set_moment_minus_2... 
              ensure moments are in the global frame
    */
    rx = (vessel->xi[i] - vessel->refOrigin.x.value);
    ry = (vessel->yi[i] - vessel->refOrigin.y.value);
    rz = (vessel->zi[i] - vessel->refOrigin.z.value);    
    mx[i] -= ((-outputType->Fz[i]*ry) - (-outputType->Fy[i]*rz));
    my[i] -= ((-outputType->Fx[i]*rz) - (-outputType->Fz[i]*rx));
    mz[i] -= ((-outputType->Fy[i]*rx) - (-outputType->Fx[i]*ry));
  };
  return MAP_SAFE;
}


/**
 * success = set_moment_minus(yType, vessel, mx, my, mz, N)
 */
MAP_ERROR_CODE set_moment_plus(const MAP_OutputType_t* outputType, const Vessel* vessel, double* mx, double* my, double* mz, const int size)
{
  int i = 0;
  double rx = 0.0;
  double ry = 0.0;
  double rz = 0.0;

  for (i=0 ; i<size ; i++) {    
    rx = (vessel->xi[i] - vessel->refOrigin.x.value);
    ry = (vessel->yi[i] - vessel->refOrigin.y.value);
    rz = (vessel->zi[i] - vessel->refOrigin.z.value);    
    mx[i] += ((-outputType->Fz[i]*ry) - (-outputType->Fy[i]*rz));
    my[i] += ((-outputType->Fx[i]*rz) - (-outputType->Fz[i]*rx));
    mz[i] += ((-outputType->Fy[i]*rx) - (-outputType->Fx[i]*ry));
  };
  return MAP_SAFE;
};


/**
 * success = increment_phi_dof_by_delta(uType, vessel, -epsilon, N);
 */
MAP_ERROR_CODE increment_phi_dof_by_delta(MAP_InputType_t* uType, const Vessel* vessel, const double delta, const int size)
{
  int i = 0;
  double R[3][3];
  double phi = vessel->orientation.phi.value;
  double rx = 0.0;
  double ry = 0.0;
  double rz = 0.0;
  R[0][0] = 1.0;   R[0][1] = 0.0;              R[0][2] = 0.0;
  R[1][0] = 0.0;   R[1][1] = cos(phi+delta);   R[1][2] = -sin(phi+delta);  
  R[2][0] = 0.0;   R[2][1] = sin(phi+delta);   R[2][2] = cos(phi+delta);

  for (i=0 ; i<size ; i++) {
    /* @todo: should also include the reference origin location? */
    rx = vessel->xi[i];
    ry = vessel->yi[i];
    rz = vessel->zi[i];
    // uType->y[i] = ry - rz*delta;
    // uType->z[i] = rz + ry*delta;
    // printf("%f  %f  %f\n",uType->x[i],uType->y[i],uType->z[i]);
    uType->x[i] = R[0][0]*rx;
    uType->y[i] = R[1][1]*ry + R[1][2]*rz;
    uType->z[i] = R[2][1]*ry + R[2][2]*rz;
  };
  return MAP_SAFE;
};


/**
 * success = increment_the_dof_by_delta(uType, vessel, -epsilon, N);
 */
MAP_ERROR_CODE increment_the_dof_by_delta(MAP_InputType_t* uType, const Vessel* vessel, const double delta, const int size)
{
  int i = 0;
  double R[3][3];
  double the = vessel->orientation.the.value;
  double rx = 0.0;
  double ry = 0.0;
  double rz = 0.0;
  R[0][0] = cos(the+delta);   R[0][1] = 0.0;   R[0][2] = sin(the+delta);
  R[1][0] = 0.0;              R[1][1] = 1.0;   R[1][2] = 0.0;
  R[2][0] = -sin(the+delta);  R[2][1] = 0.0;   R[2][2] = cos(the+delta);
  
  for (i=0 ; i<size ; i++) {
    /* @todo: should also include the reference origin location? */
    rx = vessel->xi[i];
    ry = vessel->yi[i];
    rz = vessel->zi[i];
    // uType->x[i] = rx + rz*delta;
    // uType->z[i] = rz - rx*delta;    
    uType->x[i] = R[0][0]*rx + R[0][2]*rz;
    uType->y[i] = R[1][1]*ry;
    uType->z[i] = R[2][0]*rx + R[2][2]*rz;
  };
  return MAP_SAFE;
};


/**
 * success = increment_the_dof_by_delta(uType, vessel, -epsilon, N);
 */
MAP_ERROR_CODE increment_psi_dof_by_delta(MAP_InputType_t* uType, const Vessel* vessel, const double delta, const int size)
{
  int i = 0;
  double R[3][3];
  double psi = vessel->orientation.psi.value;
  double rx = 0.0;
  double ry = 0.0;
  double rz = 0.0;
  R[0][0] = cos(psi+delta);   R[0][1] = -sin(psi+delta);   R[0][2] = 0.0;
  R[1][0] = sin(psi+delta);   R[1][1] = cos(psi+delta);    R[1][2] = 0.0;
  R[2][0] = 0.0;              R[2][1] = 0.0;               R[2][2] = 1.0;

  for (i=0 ; i<size ; i++) {
    /* @todo: should also include the reference origin location? */
    rx = vessel->xi[i];
    ry = vessel->yi[i];
    rz = vessel->zi[i];
    // uType->x[i] = rx - ry*delta;
    // uType->y[i] = ry + rx*delta;    
    uType->x[i] = R[0][0]*rx + R[0][1]*ry;
    uType->y[i] = R[1][0]*rx + R[1][1]*ry;
    uType->z[i] = R[2][2]*rz;
  };
  return MAP_SAFE;
};



MAP_ERROR_CODE set_moment_plus_2(const MAP_InputType_t* uType, const MAP_OutputType_t* outputType, const Vessel* vessel, double* mx, double* my, double* mz, const int size)
{
  int i = 0;
  double rx = 0.0;
  double ry = 0.0;
  double rz = 0.0;
  double m1 = 0.0;
  double m2 = 0.0;
  double m3 = 0.0;
  double phi = vessel->orientation.phi.value;

  double delta = 1e-4;
  for (i=0 ; i<size ; i++) {    
    /* 
       @todo forces are in global frame. Make sure everything is consistent. Should also add
       reference origin. This need s to be converted to the global frame to be consistent with
       uType. 
    */
    rx = uType->x[i];// (vessel->xi[i] - vessel->refOrigin.x.value);
    ry = uType->y[i];// (vessel->yi[i] - vessel->refOrigin.y.value);
    rz = uType->z[i];// (vessel->zi[i] - vessel->refOrigin.z.value);    
    m1 = ((-outputType->Fz[i]*ry) - (-outputType->Fy[i]*rz));
    m2 = ((-outputType->Fx[i]*rz) - (-outputType->Fz[i]*rx));
    m3 = ((-outputType->Fy[i]*rx) - (-outputType->Fx[i]*ry));
    mx[i] += m1;
    my[i] += m2*cos(phi+delta) - m3*sin(phi+delta);
    mz[i] += m2*sin(phi+delta) + m3*cos(phi+delta);    
    // mx[i] += ((-outputType->Fz[i]*ry) - (-outputType->Fy[i]*rz));
    // my[i] += ((-outputType->Fx[i]*rz) - (-outputType->Fz[i]*rx));
    // mz[i] += ((-outputType->Fy[i]*rx) - (-outputType->Fx[i]*ry));
  };
  return MAP_SAFE;
};


MAP_ERROR_CODE set_moment_minus_2(const MAP_InputType_t* uType, const MAP_OutputType_t* outputType, const Vessel* vessel, double* mx, double* my, double* mz, const int size)
{
  int i = 0;
  double rx = 0.0;
  double ry = 0.0;
  double rz = 0.0;
  double m1 = 0.0;
  double m2 = 0.0;
  double m3 = 0.0;
  double phi = vessel->orientation.phi.value;

  double delta = 1e-4;
  for (i=0 ; i<size ; i++) {    
    rx = uType->x[i];// (vessel->xi[i] - vessel->refOrigin.x.value);
    ry = uType->y[i];// (vessel->yi[i] - vessel->refOrigin.y.value);
    rz = uType->z[i];// (vessel->zi[i] - vessel->refOrigin.z.value);    
    m1 = ((-outputType->Fz[i]*ry) - (-outputType->Fy[i]*rz));
    m2 = ((-outputType->Fx[i]*rz) - (-outputType->Fz[i]*rx));
    m3 = ((-outputType->Fy[i]*rx) - (-outputType->Fx[i]*ry));
    mx[i] -= m1;
    my[i] -= m2*cos(phi+delta) - m3*sin(phi+delta);
    mz[i] -= m2*sin(phi+delta) + m3*cos(phi+delta);    
    // mx[i] -= ((-outputType->Fz[i]*ry) - (-outputType->Fy[i]*rz));
    // my[i] -= ((-outputType->Fx[i]*rz) - (-outputType->Fz[i]*rx));
    // mz[i] -= ((-outputType->Fy[i]*rx) - (-outputType->Fx[i]*ry));
  };
  return MAP_SAFE;
};
