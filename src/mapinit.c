/*
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
//#include "cminpack/minpack.h"
//#include "pyprotos.h"
#include "bstring/bstrlib.h"
#include "mapinit.h"
#include "lineroutines.h"


MAP_ERROR_CODE set_vessel(Vessel* floater, const MAP_InputType_t* u_type, char* map_msg, MAP_ERROR_CODE* ierr)
{
  MAP_ERROR_CODE success = MAP_SAFE;
  int i = 0;
  int n = u_type->x_Len;
  bstring user_msg = NULL;

  do {
    /* vessel displacement */
    success = set_vartype_float("[m]", "Vessel_X", -999, &floater->displacement.x, 0.0); CHECKERRQ(MAP_FATAL_68);
    success = set_vartype_float("[m]", "Vessel_Y", -999, &floater->displacement.y, 0.0); CHECKERRQ(MAP_FATAL_68);
    success = set_vartype_float("[m]", "Vessel_Z", -999, &floater->displacement.z, 0.0); CHECKERRQ(MAP_FATAL_68);
     
    /* vessel reference origin. When ==[0.0, 0.0, 0.0], then the reference origin is aligned with the SWL */
    success = set_vartype_float("[m]", "Vessel_Xref", -999, &floater->refOrigin.x, 0.0); CHECKERRQ(MAP_FATAL_68);
    success = set_vartype_float("[m]", "Vessel_Yref", -999, &floater->refOrigin.y, 0.0); CHECKERRQ(MAP_FATAL_68);
    success = set_vartype_float("[m]", "Vessel_Zref", -999, &floater->refOrigin.z, 0.0); CHECKERRQ(MAP_FATAL_68);
    
    /* sum force of all fairleads connecte to the vessel */
    success = set_vartype_float("[N]", "Vessel_fx", -999, &floater->lineSumForce.fx, 0.0); CHECKERRQ(MAP_FATAL_68);
    success = set_vartype_float("[N]", "Vessel_fy", -999, &floater->lineSumForce.fy, 0.0); CHECKERRQ(MAP_FATAL_68);
    success = set_vartype_float("[N]", "Vessel_fz", -999, &floater->lineSumForce.fz, 0.0); CHECKERRQ(MAP_FATAL_68);
    
    /* orientation of the vessel. This is used as input from the user */
    success = set_vartype_float("[deg]", "Vessel_phi", -999, &floater->orientation.phi, 0.0); CHECKERRQ(MAP_FATAL_68);
    success = set_vartype_float("[deg]", "Vessel_the", -999, &floater->orientation.the, 0.0); CHECKERRQ(MAP_FATAL_68);
    success = set_vartype_float("[deg]", "Vessel_psi", -999, &floater->orientation.psi, 0.0); CHECKERRQ(MAP_FATAL_68);
  } while(0);

  floater->xi = (double*)malloc(n*sizeof(double));  
  floater->yi = (double*)malloc(n*sizeof(double));  
  floater->zi = (double*)malloc(n*sizeof(double));  

  if (floater->xi==NULL || floater->yi==NULL || floater->zi==NULL) {
    return MAP_FATAL;
  };
  
  for (i=0 ; i<n ; i++) {
    floater->xi[i] = u_type->x[i];
    floater->yi[i] = u_type->y[i];
    floater->zi[i] = u_type->z[i];
  };
    
  return MAP_SAFE;
};


MAP_ERROR_CODE first_solve(ModelData* model_data, MAP_InputType_t* u_type, MAP_ConstraintStateType_t* z_type, MAP_OtherStateType_t* other_type, MAP_OutputType_t* y_type, char* map_msg, MAP_ERROR_CODE* ierr)
{
  MAP_ERROR_CODE success = MAP_SAFE;

  if (model_data->MAP_SOLVE_TYPE==MONOLITHIC) {
    success = line_solve_sequence(model_data, 0.0, map_msg, ierr); /* @todo CHECKERRQ() */
  } else {
    success = node_solve_sequence(model_data, u_type, z_type, other_type, map_msg, ierr); /* @todo CHECKERRQ() */
  };

  /* @todo replace with MAP_RETURN? */
  if (success==MAP_SAFE) {
    return MAP_SAFE;
  } else if (success==MAP_ERROR) {
    return MAP_ERROR;
  } else {
    return MAP_FATAL;
  };
};


MAP_ERROR_CODE allocate_outer_solve_data(OuterSolveAttributes* ns, const int size, char* map_msg, MAP_ERROR_CODE* ierr)
{
  int ret = 0;
  bstring user_msg = NULL;
  const int THREE = 3;  
  const int SIZE = THREE*size;
  int i = 0;

  ns->jac = (double**)malloc(SIZE*sizeof(double*));
  ns->l = (double**)malloc(SIZE*sizeof(double*));  
  ns->u = (double**)malloc(SIZE*sizeof(double*));  
  ns->x = (double*)malloc(SIZE*sizeof(double));
  ns->b = (double*)malloc(SIZE*sizeof(double));
  ns->y = (double*)malloc(SIZE*sizeof(double*));  
  
  if (ns->jac==NULL) {
    user_msg = bformat("allocate_outer_solve_data(...).");
    *ierr = map_set_universal_error(user_msg, map_msg, *ierr, MAP_FATAL_8);        
    ret = bdestroy(user_msg); 
    user_msg = NULL;
    return MAP_FATAL;
  };

  if (ns->x==NULL) {
    user_msg = bformat("allocate_outer_solve_data(...).");
    *ierr = map_set_universal_error(user_msg, map_msg, *ierr, MAP_FATAL_8);        
    ret = bdestroy(user_msg); 
    user_msg = NULL;
    return MAP_FATAL;
  };

  if (ns->b==NULL) {
    user_msg = bformat("allocate_outer_solve_data(...).");
    *ierr = map_set_universal_error(user_msg, map_msg, *ierr, MAP_FATAL_8);        
    ret = bdestroy(user_msg); 
    user_msg = NULL;
    return MAP_FATAL;
  };

  if (ns->l==NULL) {
    user_msg = bformat("allocate_outer_solve_data(...).");
    *ierr = map_set_universal_error(user_msg, map_msg, *ierr, MAP_FATAL_8);        
    ret = bdestroy(user_msg); 
    user_msg = NULL;
    return MAP_FATAL;
  };

  if (ns->u==NULL) {
    user_msg = bformat("allocate_outer_solve_data(...).");
    *ierr = map_set_universal_error(user_msg, map_msg, *ierr, MAP_FATAL_8);        
    ret = bdestroy(user_msg); 
    user_msg = NULL;
    return MAP_FATAL;
  };

  if (ns->y==NULL) {
    user_msg = bformat("allocate_outer_solve_data(...).");
    *ierr = map_set_universal_error(user_msg, map_msg, *ierr, MAP_FATAL_8);        
    ret = bdestroy(user_msg); 
    user_msg = NULL;
    return MAP_FATAL;
  };

  for(i=0 ; i<SIZE ; i++) {
    ns->jac[i] = (double*)malloc(SIZE*sizeof(double));    
    ns->l[i] = (double*)malloc(SIZE*sizeof(double));    
    ns->u[i] = (double*)malloc(SIZE*sizeof(double));    
  };

  return MAP_SAFE;
};


// /**
//  *
//  */
// MAP_ERROR_CODE free_outer_solve_data(OuterSolveAttributes* ns, const int size, char* map_msg, MAP_ERROR_CODE* ierr)
// {
//   const int SIZE = 3*size;
//   int i = 0;
// 
//   if (ns->jac) { /* is it allocated? */
//     for(i=0 ; i<SIZE ; i++) {
//       MAPFREE(ns->jac[i]);
//     };
//   };
//  
//   if (ns->l) { /* is it allocated? */
//     for(i=0 ; i<SIZE ; i++) {
//      MAPFREE(ns->l[i]);
//     };
//   };
// 
//   if (ns->u) { /* is it allocated? */
//     for(i=0 ; i<SIZE ; i++) {
//       MAPFREE(ns->u[i]);
//    };  
//   };
// 
//   MAPFREE(ns->jac);
//   MAPFREE(ns->l);
//   MAPFREE(ns->u);
//   MAPFREE(ns->b);
//   MAPFREE(ns->x);  
//   MAPFREE(ns->y);
//   return MAP_SAFE;
// };


MAP_ERROR_CODE check_help_flag(bstring list)
{
  MAP_ERROR_CODE success = 0;

  success = biseqcstrcaseless(list,"HELP"); /* string compare */
  if (success) { 
    print_help_to_screen();
  }; 
  return MAP_SAFE;
};


MAP_ERROR_CODE check_inner_f_tol_flag(struct bstrList* list, double* ftol)
{
  int success = 0;
  int n = 0;
  const char* word = NULL;

  success = biseqcstrcaseless(list->entry[0],"INNER_FTOL"); /* string compare */
  if (success==BSTR_ERR) {
    return MAP_FATAL;
  } else if (success) { 
    while (n<list->qty-1) { /* iterating through all strings */      
      if (list->entry[n+1]->slen) { /* if the string length is not 0 */
        word = list->entry[n+1]->data;
        if (is_numeric(word)) { 
          *ftol = (double)atof(word);
          return MAP_SAFE;
        } else {
          return MAP_FATAL;
        };
      };
      n++;
    };
  }; 
  word = NULL;
  return MAP_SAFE;
};


MAP_ERROR_CODE check_inner_g_tol_flag(struct bstrList* list, double* gtol)
{
  int success = 0;
  int n = 0;
  const char* word = NULL;

  success = biseqcstrcaseless(list->entry[0],"INNER_GTOL"); /* string compare */
  if (success==BSTR_ERR) {
    return MAP_FATAL;
  } else if (success) { 
    while (n<list->qty-1) { /* iterating through all strings */      
      if (list->entry[n+1]->slen) { /* if the string length is not 0 */
        word = list->entry[n+1]->data;
        if (is_numeric(word)) { 
          *gtol = (double)atof(word);
          return MAP_SAFE;
        } else {
          return MAP_FATAL;
        };
      };
      n++;
    };
  }; 
  word = NULL;
  return MAP_SAFE;
};


MAP_ERROR_CODE check_inner_x_tol_flag(struct bstrList* list, double* xtol)
{
  int success = 0;
  int n = 0;
  const char* word = NULL;

  success = biseqcstrcaseless(list->entry[0],"INNER_XTOL"); /* string compare */
  if (success==BSTR_ERR) {
    return MAP_FATAL;
  } else if (success) { 
    while (n<list->qty-1) { /* iterating through all strings */      
      if (list->entry[n+1]->slen) { /* if the string length is not 0 */
        word = list->entry[n+1]->data;
        if (is_numeric(word)) { 
          *xtol = (double)atof(word);
          return MAP_SAFE;
        } else {
          return MAP_FATAL;
        };
      };
      n++;
    };
  }; 
  word = NULL;
  return MAP_SAFE;
};



MAP_ERROR_CODE check_inner_max_its_flag(struct bstrList* list, int* max_its)
{
  int success = 0;
  int n = 0;
  const char* word = NULL;

  success = biseqcstrcaseless(list->entry[0],"INNER_MAX_ITS"); /* string compare */
  if (success==BSTR_ERR) {
    return MAP_FATAL;
  } else if (success) { 
    while (n<list->qty-1) { /* iterating through all strings */      
      if (list->entry[n+1]->slen) { /* if the string length is not 0 */
        word = list->entry[n+1]->data;
        if (is_numeric(word)) { 
          *max_its = (int)atof(word);
          return MAP_SAFE;
        } else {
          return MAP_FATAL;
        };
      };
      n++;
    };
  }; 
  word = NULL;
  return MAP_SAFE;
};


MAP_ERROR_CODE check_outer_max_its_flag(struct bstrList* list, int* max_its)
{
  int success = 0;
  int n = 0;
  const char* word = NULL;

  success = biseqcstrcaseless(list->entry[0],"OUTER_MAX_ITS"); /* string compare */
  if (success==BSTR_ERR) {
    return MAP_FATAL;
  } else if (success) { 
    while (n<list->qty-1) { /* iterating through all strings */      
      if (list->entry[n+1]->slen) { /* if the string length is not 0 */
        word = list->entry[n+1]->data;
        if (is_numeric(word)) { 
          *max_its = (int)atof(word);
          return MAP_SAFE;
        } else {
          return MAP_FATAL;
        };
      };
      n++;
    };
  }; 
  word = NULL;
  return MAP_SAFE;
};


MAP_ERROR_CODE check_outer_tol_flag(struct bstrList* list, double* outer_tol)
{
  int success = 0;
  int n = 0;
  const char* word = NULL;

  success = biseqcstrcaseless(list->entry[0],"OUTER_TOL"); /* string compare */
  if (success==BSTR_ERR) {
    return MAP_FATAL;
  } else if (success) { 
    while (n<list->qty-1) { /* iterating through all strings */      
      if (list->entry[n+1]->slen) { /* if the string length is not 0 */
        word = list->entry[n+1]->data;
        if (is_numeric(word)) { 
          *outer_tol = (double)atof(word);
          return MAP_SAFE;
        } else {
          return MAP_FATAL;
        };
      };
      n++;
    };
  }; 
  word = NULL;
  return MAP_SAFE;
};


MAP_ERROR_CODE check_outer_epsilon_flag(struct bstrList* list, double* epsilon)
{
  int success = 0;
  int n = 0;
  const char* word = NULL;

  success = biseqcstrcaseless(list->entry[0],"OUTER_EPSILON"); /* string compare */
  if (success==BSTR_ERR) {
    return MAP_FATAL;
  } else if (success) { 
    while (n<list->qty-1) { /* iterating through all strings */      
      if (list->entry[n+1]->slen) { /* if the string length is not 0 */
        word = list->entry[n+1]->data;
        if (is_numeric(word)) { 
          *epsilon = (double)atof(word);
          return MAP_SAFE;
        } else {
          return MAP_FATAL;
        };
      };
      n++;
    };
  }; 
  word = NULL;
  return MAP_SAFE;
};


MAP_ERROR_CODE check_integration_dt_flag(struct bstrList* list, double* dt)
{
  int success = 0;
  int n = 0;
  const char* word = NULL;
  success = biseqcstrcaseless(list->entry[0],"INTEGRATION_DT"); /* string compare */
  if (success==BSTR_ERR) {
    return MAP_FATAL;
  } else if (success) { 
    while (n<list->qty-1) { /* iterating through all strings */      
      if (list->entry[n+1]->slen) { /* if the string length is not 0 */
        word = list->entry[n+1]->data;
        if (is_numeric(word)) { 
          *dt = (double)atof(word);
          return MAP_WARNING;
        } else {
          return MAP_FATAL;
        };
      };
      n++;
    };
  }; 
  word = NULL;
  return MAP_SAFE;
};


MAP_ERROR_CODE check_kb_default_flag(struct bstrList* list, double* kb)
{
  int success = 0;
  int n = 0;
  const char* word = NULL;

  success = biseqcstrcaseless(list->entry[0],"KB_DEFAULT"); /* string compare */
  if (success==BSTR_ERR) {
    return MAP_FATAL;
  } else if (success) { 
    while (n<list->qty-1) { /* iterating through all strings */      
      if (list->entry[n+1]->slen) { /* if the string length is not 0 */
        word = list->entry[n+1]->data;
        if (is_numeric(word)) { 
          *kb = (double)atof(word);
          return MAP_WARNING;
        } else {
          return MAP_FATAL;
        };
      };
      n++;
    };
  }; 
  word = NULL;
  return MAP_SAFE;
};


MAP_ERROR_CODE check_cb_default_flag(struct bstrList* list, double* cb)
{
  int success = 0;
  int n = 0;
  const char* word = NULL;

  success = biseqcstrcaseless(list->entry[0],"CB_DEFAULT"); /* string compare */
  if (success==BSTR_ERR) {
    return MAP_FATAL;
  } else if (success) { 
    while (n<list->qty-1) { /* iterating through all strings */      
      if (list->entry[n+1]->slen) { /* if the string length is not 0 */
        word = list->entry[n+1]->data;
        if (is_numeric(word)) { 
          *cb = (double)atof(word);
          return MAP_WARNING;
        } else {
          return MAP_FATAL;
        };
      };
      n++;
    };
  }; 
  word = NULL;
  return MAP_SAFE;
};


MAP_ERROR_CODE check_outer_bd_flag(struct bstrList* list, FdType* bd)
{
  int success = 0;

  success = biseqcstrcaseless(list->entry[0],"OUTER_BD"); /* string compare */
  if (success) {
    *bd = BACKWARD_DIFFERENCE;
  };
  return MAP_SAFE;
};


MAP_ERROR_CODE check_outer_cd_flag(struct bstrList* list, FdType* cd)
{
  int success = 0;

  success = biseqcstrcaseless(list->entry[0],"OUTER_CD"); /* string compare */
  if (success) {
    *cd = CENTRAL_DIFFERENCE;
  };
  return MAP_SAFE;
};


MAP_ERROR_CODE check_outer_fd_flag(struct bstrList* list, FdType* fd)
{
  int success = 0;

  success = biseqcstrcaseless(list->entry[0],"OUTER_FD"); /* string compare */
  if (success) {
    *fd = FORWARD_DIFFERENCE;
  };
  return MAP_SAFE;
};


MAP_ERROR_CODE check_wave_kinematics_flag(struct bstrList* list, bool* wave)
{
  int success = 0;

  success = biseqcstrcaseless(list->entry[0],"WAVE_KINEMATICS"); /* string compare */
  if (success) {
    *wave = false;
    return MAP_FATAL;
  };
  return MAP_SAFE;
};


MAP_ERROR_CODE check_pg_cooked_flag(struct bstrList* list, OuterSolveAttributes* solver)
{
  int n = 0;
  int success = 0;
  int next = 0; 
  const char* word = NULL;

  success = biseqcstrcaseless(list->entry[0],"PG_COOKED"); /* string compare */
  if (success==BSTR_ERR) {
    return MAP_FATAL;
  } else if (success) {
    while (n<list->qty-1) { /* iterating through all strings */      
      if (list->entry[n+1]->slen) { /* if the string length is not 0 */
        word = list->entry[n+1]->data;
        if (is_numeric(word)) {         
          if (!next) {
            solver->d = (double)atof(word);
            next++;
          } else {
            solver->ds = (double)atof(word);
            solver->pg = true;
            return MAP_SAFE;
          };
        } else {
          return MAP_FATAL;
        };
      };
      n++;
    };
  }; 
  word = NULL;
  if (!next) {
    return MAP_SAFE;
  } else {
    return MAP_WARNING;
  };
};


MAP_ERROR_CODE check_repeat_flag(struct bstrList* list, ModelOptions* options)
{
  double* more_angles = NULL;
  char* current = NULL;
  int success = 0;
  int n = 0; /* word interator in the list */
  int i = 0; /* synonym for repeat_angle_size */

  success = biseqcstrcaseless(list->entry[0],"REPEAT"); /* string compare */
  if (success==BSTR_ERR) {
    return MAP_FATAL;
  } else if (success) {
    while (n<list->qty-1) { /* iterating through all strings */      
      if (list->entry[n+1]->slen) { /* if the string length is not 0 */
        current = list->entry[n+1]->data;
        i = options->repeat_angle_size;
        more_angles = realloc(options->repeat_angle, (i+1)*sizeof(double));
        if (more_angles) {
          options->repeat_angle = more_angles;
          if (is_numeric(current)) { /* add repeat angle if word is numeric */
            options->repeat_angle[i] = atof(current);
            options->repeat_angle_size++;
          } else { 
            MAPFREE(more_angles);
            return MAP_FATAL;
          };
        } else {
          MAPFREE(more_angles);
          return MAP_FATAL;
        };        
      };
      n++;
    };
  };
  return MAP_SAFE;
};


MAP_ERROR_CODE check_ref_position_flag(struct bstrList* list, Point* ref_position)
{
  int n = 0;
  int success = 0;
  int next = 0; 
  const char* word = NULL;

  success = biseqcstrcaseless(list->entry[0],"REF_POSITION"); /* string compare */
  if (success==BSTR_ERR) {
    return MAP_FATAL;
  } else if (success) {
    while (n<list->qty-1) { /* iterating through all strings */      
      if (list->entry[n+1]->slen) { /* if the string length is not 0 */
        word = list->entry[n+1]->data;
        if (is_numeric(word)) {         
          if (!next) {
            ref_position->x.value = (double)atof(word);
            next++;
          } else if (next==1) {
            ref_position->y.value = (double)atof(word);
            next++;
          } else {
            ref_position->z.value = (double)atof(word);
            return MAP_SAFE;
          };
        };
      };
      n++;
    };
    return MAP_WARNING;
  };
  word = NULL;
  return MAP_SAFE;
};


MAP_ERROR_CODE check_uncaught_flag(struct bstrList* list)
{
  int success = 0;  

  if (biseqcstrcaseless(list->entry[0],"")) {
    return MAP_SAFE;
  } else if (biseqcstrcaseless(list->entry[0],"HELP")) {
    return MAP_SAFE;
  } else if (biseqcstrcaseless(list->entry[0],"INNER_FTOL")) { 
    return MAP_SAFE;
  } else if (biseqcstrcaseless(list->entry[0],"INNER_GTOL")) {
    return MAP_SAFE;
  } else if (biseqcstrcaseless(list->entry[0],"INNER_XTOL")) {
    return MAP_SAFE;
  } else if (biseqcstrcaseless(list->entry[0],"INNER_MAX_ITS")) { 
    return MAP_SAFE;
  } else if (biseqcstrcaseless(list->entry[0],"OUTER_MAX_ITS")) {
    return MAP_SAFE;
  } else if (biseqcstrcaseless(list->entry[0],"OUTER_TOL")) {
    return MAP_SAFE;
  } else if (biseqcstrcaseless(list->entry[0],"OUTER_EPSILON")) {
    return MAP_SAFE;
  } else if (biseqcstrcaseless(list->entry[0],"INTEGRATION_DT")) {
    return MAP_SAFE;
  } else if (biseqcstrcaseless(list->entry[0],"KB_DEFAULT")) {
    return MAP_SAFE;
  } else if (biseqcstrcaseless(list->entry[0],"CB_DEFAULT")) {
    return MAP_SAFE;
  } else if (biseqcstrcaseless(list->entry[0],"OUTER_CD")) {
    return MAP_SAFE;
  } else if (biseqcstrcaseless(list->entry[0],"OUTER_BD")) {
    return MAP_SAFE;
  } else if (biseqcstrcaseless(list->entry[0],"OUTER_FD")) {
    return MAP_SAFE;
  } else if (biseqcstrcaseless(list->entry[0],"WAVE_KINEMATICS")) {
    return MAP_SAFE;
  } else if (biseqcstrcaseless(list->entry[0],"PG_COOKED")) {
    return MAP_SAFE;
  } else if (biseqcstrcaseless(list->entry[0],"REPEAT")) {
    return MAP_SAFE;
  } else if (biseqcstrcaseless(list->entry[0],"REF_POSITION")) {
    return MAP_SAFE;
  }; 
  return MAP_WARNING;
};


MAP_ERROR_CODE set_library_diameter(bstring word, CableLibrary* library_ptr)
{
  if (is_numeric(word->data)) { 
    library_ptr->diam = (double)atof(word->data);
  } else { 
    return MAP_FATAL;
  };  
  return MAP_SAFE;
};


MAP_ERROR_CODE set_library_mass_density(bstring word, CableLibrary* library_ptr)
{  
  if (is_numeric(word->data)) { 
    library_ptr->massDensityInAir = (double)atof(word->data);
  } else { 
    return MAP_FATAL;
  };  
  return MAP_SAFE;
};


MAP_ERROR_CODE set_library_ea(bstring word, CableLibrary* library_ptr)
{
  if (is_numeric(word->data)) { 
    library_ptr->ea = (double)atof(word->data);
  } else { 
    return MAP_FATAL;
  };  
  return MAP_SAFE;
};


MAP_ERROR_CODE set_library_cb(bstring word, CableLibrary* library_ptr)
{
  if (is_numeric(word->data)) { 
    library_ptr->cb = (double)atof(word->data);
  } else { 
    return MAP_FATAL;
  };  
  return MAP_SAFE;
};


MAP_ERROR_CODE set_library_internal_damping(bstring word, CableLibrary* library_ptr)
{
  if (is_numeric(word->data)) { 
    library_ptr->cIntDamp = (double)atof(word->data);
  } else { 
    return MAP_FATAL;
  };  
  return MAP_SAFE;
};


MAP_ERROR_CODE set_library_added_mass_coefficient(bstring word, CableLibrary* library_ptr)
{
  if (is_numeric(word->data)) { 
    library_ptr->cAdded = (double)atof(word->data);
  } else { 
    return MAP_FATAL;
  };  
  return MAP_SAFE;
};


MAP_ERROR_CODE set_library_cross_flow_drag_coefficient(bstring word, CableLibrary* library_ptr)
{
  if (is_numeric(word->data)) { 
    library_ptr->cDragNormal = (double)atof(word->data);
  } else { 
    return MAP_FATAL;
  };  
  return MAP_SAFE;
};


MAP_ERROR_CODE set_library_tangent_drag_coefficient(bstring word, CableLibrary* library_ptr)
{
  if (is_numeric(word->data)) { 
    library_ptr->cDragTangent = (double)atof(word->data);
  } else { 
    return MAP_FATAL;
  };  
  return MAP_SAFE;
};



MAP_ERROR_CODE set_model_options_list(ModelData* model_data, InitializationData* init_data, char* map_msg, MAP_ERROR_CODE* ierr)
{
  MAP_ERROR_CODE success = MAP_SAFE;
  bstring user_msg = NULL;
  int i = 0;
  const int n_lines = (init_data->solverOptionsString->qty)-1;
  struct bstrList* parsed = NULL;
  struct tagbstring tokens; 

  cstr2tbstr(tokens," \t\n\r"); /* token for splitting line into indivdual words is a tab and space */   
  for (i=0 ; i<=n_lines ; i++) { 
    parsed = bsplits(init_data->solverOptionsString->entry[i], &tokens);
    do {
      success = check_help_flag(parsed->entry[0]); CHECKERRQ(MAP_FATAL_85);
      success = check_inner_f_tol_flag(parsed, &model_data->inner_loop.f_tol); CHECKERRK(MAP_ERROR_2);
      success = check_outer_max_its_flag(parsed, &model_data->outer_loop.maxIts); CHECKERRK(MAP_ERROR_3);
      success = check_inner_max_its_flag(parsed, &model_data->inner_loop.max_its); CHECKERRK(MAP_ERROR_4);
      success = check_inner_g_tol_flag(parsed, &model_data->inner_loop.g_tol); CHECKERRK(MAP_ERROR_9);
      success = check_inner_x_tol_flag(parsed, &model_data->inner_loop.x_tol); CHECKERRK(MAP_ERROR_10);
      success = check_outer_tol_flag(parsed, &model_data->outer_loop.tol); CHECKERRK(MAP_ERROR_3);
      success = check_outer_epsilon_flag(parsed, &model_data->outer_loop.epsilon); CHECKERRK(MAP_ERROR_3);
      success = check_integration_dt_flag(parsed, &model_data->modelOptions.integrationDt); CHECKERRK(MAP_ERROR_15); 
      success = check_kb_default_flag(parsed, &model_data->modelOptions.kbLm); CHECKERRK(MAP_ERROR_16); 
      success = check_cb_default_flag(parsed, &model_data->modelOptions.cbLm); CHECKERRK(MAP_ERROR_17); 
      success = check_outer_bd_flag(parsed, &model_data->outer_loop.fd);
      success = check_outer_cd_flag(parsed, &model_data->outer_loop.fd);
      success = check_outer_fd_flag(parsed, &model_data->outer_loop.fd);      
      success = check_wave_kinematics_flag(parsed, &model_data->modelOptions.waveKinematics); CHECKERRK(MAP_WARNING_10);
      success = check_pg_cooked_flag(parsed, &model_data->outer_loop); CHECKERRK(MAP_WARNING_8);
      success = check_repeat_flag(parsed, &model_data->modelOptions); CHECKERRQ(MAP_FATAL_34);
      success = check_ref_position_flag(parsed, &model_data->vessel.refOrigin); CHECKERRQ(MAP_FATAL_36);
      success = check_uncaught_flag(parsed);       
      if (success) {
        user_msg = bformat("word: <%s>", parsed->entry[0]->data);
        *ierr = map_set_universal_error(user_msg, map_msg, *ierr, MAP_WARNING_1);        
        success = bdestroy(user_msg); /* clear user_msg to not inavertendly picked up elsewhere */
        user_msg = NULL;
      };
    } while (0);   
    success = bstrListDestroy(parsed);
  };
  MAP_RETURN;
};


MAP_ERROR_CODE reset_cable_library(CableLibrary* library_ptr)
{
  MAP_ERROR_CODE success = MAP_SAFE;

  library_ptr->diam = 0.0;
  library_ptr->massDensityInAir = 0.0;
  library_ptr->ea = 0.0;          
  library_ptr->omega = 0.0;       
  library_ptr->a = 0.0;           
  library_ptr->cb = 0.0;          
  library_ptr->cIntDamp = 0.0;    
  library_ptr->cAdded = 0.0;      
  library_ptr->cDragNormal = 0.0; 
  library_ptr->cDragTangent = 0.0; 
  // if (library_ptr->label) {
  //   success = bdestroy(library_ptr->label);
  // };
  return MAP_SAFE;
};


MAP_ERROR_CODE set_cable_library_list(ModelData* model_data, InitializationData* init_data, char* map_msg, MAP_ERROR_CODE* ierr)
{
  MAP_ERROR_CODE success = MAP_SAFE;
  bstring user_msg = NULL;
  int i = 0;
  int n = 0;
  int next = 0; 
  int ret = 0;
  const int n_lines = (init_data->libraryInputString->qty)-1;
  struct bstrList* parsed = NULL;
  struct tagbstring tokens; 
  CableLibrary new_cable_library;
  CableLibrary* library_iter = NULL;
  
  cstr2tbstr(tokens," \t\n\r"); /* token for splitting line into indivdual words is a tab and space */   
  success = reset_cable_library(&new_cable_library);

  for (i=0 ; i<=n_lines ; i++) { 
    list_append(&model_data->cableLibrary, &new_cable_library);
    library_iter = (CableLibrary*)list_get_at(&model_data->cableLibrary, i);

    parsed = bsplits(init_data->libraryInputString->entry[i], &tokens);
    n = 0;
    next = 0;
    do {  
      while (n<parsed->qty-1) { /* iterating through all strings */              
        if (parsed->entry[n]->slen) { /* if the string length is not 0 */
          if (next==0) {
            library_iter->label = bstrcpy(parsed->entry[n]);                         
            next++;
          } else if (next==1) {
             success = set_library_diameter(parsed->entry[n], library_iter); CHECKERRQ(MAP_FATAL_12);
             next++;            
          } else if (next==2) {
             success = set_library_mass_density(parsed->entry[n], library_iter); CHECKERRQ(MAP_FATAL_13);
             next++;
          } else if (next==3) {
            success = set_library_ea(parsed->entry[n], library_iter); CHECKERRQ(MAP_FATAL_14);
            next++;
          } else if (next==4) {
            success = set_library_cb(parsed->entry[n], library_iter); CHECKERRQ(MAP_FATAL_15);
            next++;
          } else if (next==5) {
            success = set_library_internal_damping(parsed->entry[n], library_iter); CHECKERRQ(MAP_FATAL_81);
            next++;
          } else if (next==6) {
            success = set_library_added_mass_coefficient(parsed->entry[n], library_iter); CHECKERRQ(MAP_FATAL_82);
            next++;
          } else if (next==7) {
            success = set_library_cross_flow_drag_coefficient(parsed->entry[n], library_iter); CHECKERRQ(MAP_FATAL_83);
            next++;
          } else if (next==8) {
            success = set_library_tangent_drag_coefficient(parsed->entry[n], library_iter); CHECKERRQ(MAP_FATAL_84);
            next++;
          };
        };
        n++;
      };
    } while (0);   
    // list_append(&model_data->cableLibrary, &new_cable_library);
    // success = reset_cable_library(&new_cable_library);
    success = bstrListDestroy(parsed);
  };
  model_data->sizeOfCableLibrary = list_size(&model_data->cableLibrary); /* SimCList routine */
  MAP_RETURN;
};


MAP_ERROR_CODE initialize_cable_library_variables(ModelData* model_data, MAP_ParameterType_t* p_type, char* map_msg, MAP_ERROR_CODE* ierr)
{
  MAP_ERROR_CODE success = MAP_SAFE;
  double radius = 0.0;
  double area = 0.0;
  double muCable = 0.0;
  double rhoFluid = 0.0; 
  const double g = p_type->g;
  const double PI = 3.14159264;
  CableLibrary* library_iter = NULL;
  bstring user_msg = NULL;

  list_iterator_start(&model_data->cableLibrary); /* starting an iteration "session" */
  while ( list_iterator_hasnext(&model_data->cableLibrary)) { /* tell whether more values available */ 
    library_iter = (CableLibrary*)list_iterator_next(&model_data->cableLibrary);
    radius = library_iter->diam/2;
    area = PI*pow(radius,2);
    muCable = library_iter->massDensityInAir;
    rhoFluid = p_type->rhoSea;
    library_iter->omega = g*(muCable-area*rhoFluid);

    library_iter->a = area;
    if (fabs(library_iter->omega)<=1) {
      user_msg = bformat("omega = %f <= 1.0", library_iter->omega);
      *ierr = map_set_universal_error(user_msg, map_msg, *ierr, MAP_WARNING_5);
      success = bdestroy(user_msg); /* clear user_msg to not inavertendly picked up elsewhere */
      user_msg = NULL;
    };
  };
  list_iterator_stop(&model_data->cableLibrary); /* ending the iteration "session" */    
  
  if (fabs(library_iter->omega)<=1e-3) {
    return MAP_FATAL;
  }
  return MAP_SAFE;
};


MAP_ERROR_CODE expand_node_number(const int n_line, bstring line)
{
  bstring current_entry = NULL;
  int ret = 0;

  current_entry = bformat("%d   ", n_line);              
  ret = bconcat(line, current_entry);
  ret = bdestroy(current_entry);               
  return MAP_SAFE;
};

MAP_ERROR_CODE expand_node_type(const char* word, bstring line)
{
  bstring current_entry = NULL;
  int ret = 0;

  current_entry = bformat("%s   ", word);              
  ret = bconcat(line, current_entry);
  ret = bdestroy(current_entry);
  return MAP_SAFE;
};


MAP_ERROR_CODE expand_node_position_x(double* x, const char* word)
{
  if (is_numeric(word)) { /* if number is numeric */
    *x = (double)atof(word);                
  } else {
    return MAP_FATAL;
  };
  return MAP_SAFE;
};


MAP_ERROR_CODE expand_node_position_y(double* y, const char* word)
{
  if (is_numeric(word)) { /* if number is numeric */
    *y = (double)atof(word);                
  } else {
    return MAP_FATAL;
  };
  return MAP_SAFE;
};


MAP_ERROR_CODE expand_node_position_z(Vector* position, const double angle, const double x, const double y, const char* word, bstring line)
{
  bstring current_entry = NULL;
  int ret = 0;

  position->x =  x*cos(angle) + y*sin(angle);
  position->y = -x*sin(angle) + y*cos(angle);                    
  if (is_numeric(word)) { /* if number is numeric */
    position->z = (double)atof(word);                
  } else {
    return MAP_FATAL;
  };
  current_entry = bformat("%1.4f   %1.4f   %1.4f   ",position->x, position->y, position->z);
  ret = bconcat(line, current_entry);
  ret = bdestroy(current_entry);
  return MAP_SAFE;
};


MAP_ERROR_CODE expand_node_mass(const char* word, bstring line)
{
  bstring current_entry = NULL;
  int ret = 0;

  current_entry = bformat("%s   ", word);              
  ret = bconcat(line, current_entry);
  ret = bdestroy(current_entry);
  return MAP_SAFE;
};


MAP_ERROR_CODE expand_node_buoyancy(const char* word, bstring line)
{
  bstring current_entry = NULL;
  int ret = 0;

  current_entry = bformat("%s   ", word);              
  ret = bconcat(line, current_entry);
  ret = bdestroy(current_entry);
  return MAP_SAFE;
};


MAP_ERROR_CODE expand_node_force_x(double* fx, const char* word)
{
  if (is_numeric(word)) { /* if number is numeric */
    *fx = (double)atof(word);                
  } else if (word[0]=='#') { /* if the nuymber is iterated */
    if (is_numeric(remove_first_character(word))) { 
      *fx = (double)atof(remove_first_character(word));
    } else { /* in this case, it is presumed the force is just '#' */
      *fx = 0.0;                  
    };
  } else {
    return MAP_FATAL;
  };
};


MAP_ERROR_CODE expand_node_force_y(double* fy, const char* word)
{
  if (is_numeric(word)) { /* if number is numeric */
    *fy = (double)atof(word);                
  } else if (word[0]=='#') { /* if the nuymber is iterated */
    if (is_numeric(remove_first_character(word))) { 
      *fy = (double)atof(remove_first_character(word));
    } else { /* in this case, it is presumed the force is just '#' */
      *fy = 0.0;                  
    };
  } else {
    return MAP_FATAL;
  };
};


MAP_ERROR_CODE expand_node_force_z(Vector* force, const double angle, const double fx, const double fy, const char* word, bstring line)
{
  bstring current_entry = NULL;
  int ret = 0;

  force->x =  fx*cos(angle) + fy*sin(angle);
  force->y = -fx*sin(angle) + fy*cos(angle);                    
  if (is_numeric(word)) { /* if number is numeric */
    force->z = (double)atof(word);                
  } else if (word[0]=='#') { /* if the nuymber is iterated */
    if (is_numeric(remove_first_character(word))) { 
      force->z = (double)atof(remove_first_character(word));
    } else { /* in this case, it is presumed the force is just '#' */
      force->z = 0.0;                  
    };
  } else {
    return MAP_FATAL;
  };
  if (word[0]=='#') { 
    force->z = (MapReal)atof(remove_first_character(word));
    current_entry = bformat("#%1.4f   #%1.4f   #%1.4f\n",force->x, force->y, force->z);              
  } else {
    force->z = (MapReal)atof(word);
    current_entry = bformat("%1.4f   %1.4f   %1.4f\n",force->x, force->y, force->z);              
  };
  ret = bconcat(line, current_entry);
  ret = bdestroy(current_entry);               
  return MAP_SAFE;
};


MAP_ERROR_CODE repeat_nodes(ModelData* model_data, InitializationData* init_data, char* map_msg, MAP_ERROR_CODE* ierr)
{
  MAP_ERROR_CODE success = MAP_SAFE;
  int i = 0;
  int j = 0;
  int next = 0; 
  int i_parsed = 0;
  int n_line = 0;
  const int num_repeat = model_data->modelOptions.repeat_angle_size; 
  const int num_node = init_data->nodeInputString->qty;
  const int n = (num_node)*(num_repeat+1);
  const char* word = NULL;  
  double x_position = 0.0;
  double y_position = 0.0;
  double x_force = 0.0;
  double y_force = 0.0;
  double current_angle = 0.0;
  bstring line = bformat("");            
  struct bstrList* parsed = NULL;
  struct tagbstring tokens; 
  Vector force;     
  Vector position;
  Node new_node; 

  cstr2tbstr(tokens," \t\n\r"); /* token for splitting line into indivdual words is a tab and space */   

  force.x = 0.0;
  force.y = 0.0;
  force.z = 0.0;
  position.x = 0.0;
  position.y = 0.0;
  position.z = 0.0;
   
  /* allocate space needed to expand the number of node lines */
  success = bstrListAlloc(init_data->expandedNodeInputString, n+1); 
  init_data->expandedNodeInputString->qty = 0;
  
  for(i=0 ; i<num_node ; i++) {     
    init_data->expandedNodeInputString->entry[i] = bfromcstr(init_data->nodeInputString->entry[i]->data);// bstrcpy(init_data->nodeInputString->entry[i]);
    init_data->expandedNodeInputString->qty++;
  };

  for(i=0 ; i<num_repeat ; i++) { /* this is skipped if not repeat angles are declared */
    for(j=0 ; j<num_node ; j++) { 
      success = reset_node(&new_node);  
      n_line = (i+1)*num_node + j;
      current_angle = model_data->modelOptions.repeat_angle[i]*(DEG2RAD);
      parsed = bsplits(init_data->nodeInputString->entry[j], &tokens);
      next = 0;
      i_parsed = 0;
      do {  
        while (i_parsed<parsed->qty-1) { /* iterating through all strings */              
          if (parsed->entry[i_parsed]->slen) { /* if the string length is not 0 */
            word = parsed->entry[i_parsed]->data;
            if (next==0) {
              success = expand_node_number(n_line+1, line);/* @todo: checkerrq */
              next++;
            } else if (next==1) {
              success = expand_node_type(word, line);/* @todo: checkerrq */
              next++;
            } else if (next==2) {
              success = expand_node_position_x(&x_position, word);/* @todo: checkerrq */
              next++;
            } else if (next==3) {
              success = expand_node_position_y(&y_position, word);/* @todo: checkerrq */
              next++;
            } else if (next==4) {
              success = expand_node_position_z(&position, current_angle, x_position, y_position, word, line);/* @todo: checkerrq */
              next++;
            } else if (next==5) { /* node mass */
              success = expand_node_mass(word, line);/* @todo: checkerrq */
              next++;
            } else if (next==6) { /* node buoyancy */
              success = expand_node_buoyancy(word, line);/* @todo: checkerrq */
              next++;
            } else if (next==7) {
              success = expand_node_force_x(&x_force, word);/* @todo: checkerrq */
              next++;
            } else if (next==8) {
              success = expand_node_force_y(&y_force, word);/* @todo: checkerrq */
              next++;
            } else if (next==9) {
              success = expand_node_force_z(&force, current_angle, x_force, y_force, word, line);/* @todo: checkerrq */
              next++;
            };
          };
          i_parsed++;
        };
        init_data->expandedNodeInputString->qty++;
        init_data->expandedNodeInputString->entry[n_line] = bstrcpy(line);
        success = bassigncstr(line, "");
      } while (0);   
      success = bstrListDestroy(parsed);
    };  
  };
  success = bdestroy(line);               

  MAP_RETURN;
};


MAP_ERROR_CODE expand_element_number(const int n_line, bstring line)
{
  bstring current_entry = NULL;
  int ret = 0;
              
  current_entry = bformat("%d   ", n_line);              
  ret = bconcat(line, current_entry);
  ret = bdestroy(current_entry);               
  return MAP_SAFE;
};


MAP_ERROR_CODE expand_element_property_name(const char* word, bstring line)
{
  int ret = 0;
  bstring current_entry = NULL;

  current_entry = bformat("%s   ", word);              
  ret = bconcat(line, current_entry);
  ret = bdestroy(current_entry);
  return MAP_SAFE;
};


MAP_ERROR_CODE expand_element_length(const char* word, bstring line)
{
  int ret = 0;
  bstring current_entry = NULL;
  
  current_entry = bformat("%s   ", word);              
  ret = bconcat(line, current_entry);
  ret = bdestroy(current_entry);
  return MAP_SAFE;
};


MAP_ERROR_CODE expand_element_anchor_number(const char* word, const int index, const int n, bstring line)
{
  int ret = 0;
  bstring current_entry = NULL;
  int anchor_num = 0;

  if (is_numeric(word)) {
    anchor_num = (int)atoi(word);
  } else {
    return MAP_SAFE;
  };  
  current_entry = bformat("%d   ", (index+1)*n+anchor_num);
  ret = bconcat(line, current_entry);
  ret = bdestroy(current_entry);
  return MAP_SAFE;
};


MAP_ERROR_CODE expand_element_fairlead_number(const char* word, const int index, const int n, bstring line)
{
  int ret = 0;
  bstring current_entry = NULL;
  int fairlead_num = 0;

  if (is_numeric(word)) {
    fairlead_num = (int)atoi(word);
  } else {
    return MAP_SAFE;
  };  
  current_entry = bformat("%d   ", (index+1)*n+fairlead_num);
  ret = bconcat(line, current_entry);
  ret = bdestroy(current_entry);
  return MAP_SAFE;
};

MAP_ERROR_CODE expand_element_flag(const char* word, bstring line)
{
  int ret = 0;
  bstring current_entry = NULL;
  
  current_entry = bformat("%s   ", word);              
  ret = bconcat(line, current_entry);
  ret = bdestroy(current_entry);
  return MAP_SAFE;
};


MAP_ERROR_CODE repeat_elements(ModelData* model_data, InitializationData* init_data, char* map_msg, MAP_ERROR_CODE* ierr)
{
  MAP_ERROR_CODE success = MAP_SAFE;
  int i = 0;
  int j = 0;
  int next = 0; 
  int i_parsed = 0;
  int n_line = 0;
  const int num_repeat = model_data->modelOptions.repeat_angle_size; 
  const int num_element = init_data->elementInputString->qty;
  const int num_node = init_data->nodeInputString->qty;
  const int n = (num_element)*(num_repeat+1);
  const char* word = NULL;  
  double current_angle = 0.0;
  bstring line = bformat("");            
  struct bstrList* parsed = NULL;
  struct tagbstring tokens; 
  Element new_element;

  cstr2tbstr(tokens," \t\n\r"); /* token for splitting line into indivdual words is a tab and space */   

  /* allocate space needed to expand the number of node lines */
  success = bstrListAlloc(init_data->expandedElementInputString, n+1); 
  init_data->expandedElementInputString->qty = 0;
  
  for(i=0 ; i<num_element ; i++) {     
    init_data->expandedElementInputString->entry[i] = bfromcstr(init_data->elementInputString->entry[i]->data);
    init_data->expandedElementInputString->qty++;
  };
  

  for(i=0 ; i<num_repeat ; i++) { /* this is skipped if not repeat angles are declared */
    for(j=0 ; j<num_element ; j++) { 
      success = reset_element(&new_element);  
      n_line = (i+1)*num_element + j;
      current_angle = model_data->modelOptions.repeat_angle[i]*(DEG2RAD);
      parsed = bsplits(init_data->elementInputString->entry[j], &tokens);
      next = 0;
      i_parsed = 0;
      do {  
        while (i_parsed<parsed->qty-1) { /* iterating through all strings */              
          if (parsed->entry[i_parsed]->slen) { /* if the string length is not 0 */
            word = parsed->entry[i_parsed]->data;
            if (next==0) {
              success = expand_element_number(n_line+1, line);
              next++;
            } else if (next==1) {
              success = expand_element_property_name(word, line);
              next++;
            } else if (next==2) {
              success = expand_element_length(word, line);
              next++;
            } else if (next==3) {
              success = expand_element_anchor_number(word, i, num_node, line);
              next++;
            } else if (next==4) {
              success = expand_element_fairlead_number(word, i, num_node, line);
              next++;
            } else {
              success = expand_element_flag(word, line);
              next++;
            };
          };
          i_parsed++;
        };
        init_data->expandedElementInputString->qty++;
        init_data->expandedElementInputString->entry[n_line] = bstrcpy(line);
        success = bassigncstr(line, "");
      } while (0);      
      success = bstrListDestroy(parsed);
    };
  };
  success = bdestroy(line);               
  
  // for(i=0 ; i<init_data->expandedElementInputString->qty ; i++) {
  //   printf("%s\n",init_data->expandedElementInputString->entry[i]->data);
  // };

  MAP_RETURN;
};



MAP_ERROR_CODE allocate_types_for_nodes(MAP_InputType_t* u_type, MAP_ConstraintStateType_t* z_type, MAP_OtherStateType_t* other_type, MAP_OutputType_t* y_type, ModelData* model_data, struct bstrList* node_input_string, char* map_msg, MAP_ERROR_CODE* ierr)
{
  int i = 0;
  int i_parsed = 0;
  int next = 0;
  int fix_num = 0;
  int vessel_num = 0;
  int connect_num = 0;
  MAP_ERROR_CODE success = MAP_SAFE;
  bstring user_msg = NULL;
  const int num_nodes = node_input_string->qty;
  struct bstrList* parsed = NULL;
  struct tagbstring tokens; 
  
  cstr2tbstr(tokens," \t\n\r"); /* token for splitting line into indivdual words is a tab and space */   

   /* First need to run loop to find the number of inputs, outputs, constraints.
    * Basically we are just counting the number of VESSEL, FIX, and CONNECT nodes     
    */
  for(i=0 ; i<num_nodes ; i++) {             
    i_parsed = 0;
    next = 0;
    parsed = bsplits(node_input_string->entry[i], &tokens);
    do {  
      while (i_parsed<parsed->qty-1) { /* iterating through all strings */              
        if (parsed->entry[i_parsed]->slen) { /* if the string length is not 0 */
          if (next==1) {
            if (biseqcstrcaseless(parsed->entry[i_parsed],"FIX")) {
              fix_num++;
              break; /* break the while-loop because the agenda is reached */
            } else if (biseqcstrcaseless(parsed->entry[i_parsed],"CONNECT")) {
              connect_num++;
              break; /* break the while-loop because the agenda is reached */
            } else if (biseqcstrcaseless(parsed->entry[i_parsed],"VESSEL")) {
              vessel_num++;
              break; /* break the while-loop because the agenda is reached */
            } else {
              user_msg = bformat("Value: <%s>", parsed->entry[i_parsed]->data);
              *ierr = map_set_universal_error(user_msg, map_msg, *ierr, MAP_FATAL_25);        
              success = bdestroy(user_msg); /* clear user_msg to not inavertendly picked up elsewhere */
              user_msg = NULL;
            };
          };          
          next++;
        };
        i_parsed++;
      };
    } while (0);
    success = bstrListDestroy(parsed);
  }; 

  other_type->x_Len = fix_num;
  other_type->y_Len = fix_num;
  other_type->z_Len = fix_num;
  other_type->x = (double*)malloc(other_type->x_Len*sizeof(double));
  other_type->y = (double*)malloc(other_type->y_Len*sizeof(double));
  other_type->z = (double*)malloc(other_type->z_Len*sizeof(double));

  /* If the node is VESSEL, then the applied force is an output state. Otherwise, 
   * it has to be an other state because it can't be associated with any other type. 
   * This is what is done below.
   */
  other_type->Fx_connect_Len = connect_num;
  other_type->Fy_connect_Len = connect_num;
  other_type->Fz_connect_Len = connect_num;
  other_type->Fx_connect = (double*)malloc(other_type->Fx_connect_Len*sizeof(double));
  other_type->Fy_connect = (double*)malloc(other_type->Fy_connect_Len*sizeof(double));
  other_type->Fz_connect = (double*)malloc(other_type->Fz_connect_Len*sizeof(double));

  other_type->Fx_anchor_Len = fix_num;
  other_type->Fy_anchor_Len = fix_num;
  other_type->Fz_anchor_Len = fix_num;
  other_type->Fx_anchor = (double*)malloc(other_type->Fx_anchor_Len*sizeof(double));
  other_type->Fy_anchor = (double*)malloc(other_type->Fy_anchor_Len*sizeof(double));
  other_type->Fz_anchor = (double*)malloc(other_type->Fz_anchor_Len*sizeof(double));
  
  z_type->x_Len = connect_num;          
  z_type->y_Len = connect_num;          
  z_type->z_Len = connect_num;          
  z_type->x = (double*)malloc(z_type->x_Len*sizeof(double));
  z_type->y = (double*)malloc(z_type->y_Len*sizeof(double));
  z_type->z = (double*)malloc(z_type->z_Len*sizeof(double));
  
  u_type->x_Len = vessel_num;
  u_type->y_Len = vessel_num;
  u_type->z_Len = vessel_num;
  u_type->x = (double*)malloc(u_type->x_Len*sizeof(double));
  u_type->y = (double*)malloc(u_type->y_Len*sizeof(double));
  u_type->z = (double*)malloc(u_type->z_Len*sizeof(double));

  y_type->Fx_Len = vessel_num;
  y_type->Fy_Len = vessel_num;
  y_type->Fz_Len = vessel_num;
  y_type->Fx = (double*)malloc(y_type->Fx_Len*sizeof(double));
  y_type->Fy = (double*)malloc(y_type->Fy_Len*sizeof(double));
  y_type->Fz = (double*)malloc(y_type->Fz_Len*sizeof(double));

  return MAP_SAFE;
};


MAP_ERROR_CODE compare_length(int a, int b)
{
  if (a!=b) {
    return MAP_FATAL;
  }; 
  return MAP_SAFE;
};


MAP_ERROR_CODE set_node_list(const MAP_ParameterType_t* p_type,  MAP_InputType_t* u_type, MAP_ConstraintStateType_t* z_type, MAP_OtherStateType_t* other_type, MAP_OutputType_t* y_type, ModelData* model_data, struct bstrList* node_input_string, char* map_msg, MAP_ERROR_CODE* ierr)
{
  MAP_ERROR_CODE success = MAP_SAFE;
  int i = 0;
  int i_parsed = 0;
  int next = 0;
  const int num_nodes = node_input_string->qty;
  int fix_num = 0;
  int vessel_num = 0;
  int connect_num = 0;
  Node new_node;
  Node* node_iter = NULL;
  struct bstrList* parsed = NULL;
  struct tagbstring tokens; 
  bstring user_msg = NULL;
  bstring alias = NULL;
  bstring value_string = NULL;
  const double depth = p_type->depth;

  cstr2tbstr(tokens," \t\n\r"); /* token for splitting line into indivdual words is a tab and space */   

  success = allocate_types_for_nodes(u_type, z_type, other_type, y_type, model_data, node_input_string, map_msg, ierr);
  success = reset_node(&new_node); /* create an empty node */
   
  for(i=0 ; i<num_nodes ; i++) {         
    list_append(&model_data->node, &new_node); /* append node to list */
    node_iter = (Node*)list_get_at(&model_data->node, i);

    i_parsed = 0;
    next = 0;
    parsed = bsplits(node_input_string->entry[i], &tokens);
    do {  
      while (i_parsed<parsed->qty-1) { /* iterating through all strings */              
        if (parsed->entry[i_parsed]->slen) { /* if the string length is not 0 */
          if (next==0) {            
            next++;
          } else if (next==1) {
            if (biseqcstrcaseless(parsed->entry[i_parsed],"FIX")) {
              node_iter->type = FIX;
              fix_num++;                   /* VarTypePtr              FAST derived  array index */
              success = associate_vartype_ptr(&node_iter->positionPtr.x, other_type->x, fix_num);
              success = associate_vartype_ptr(&node_iter->positionPtr.y, other_type->y, fix_num);
              success = associate_vartype_ptr(&node_iter->positionPtr.z, other_type->z, fix_num);
              success = associate_vartype_ptr(&node_iter->sumForcePtr.fx, other_type->Fx_anchor, fix_num);
              success = associate_vartype_ptr(&node_iter->sumForcePtr.fy, other_type->Fy_anchor, fix_num);
              success = associate_vartype_ptr(&node_iter->sumForcePtr.fz, other_type->Fz_anchor, fix_num);
            } else if (biseqcstrcaseless(parsed->entry[i_parsed],"CONNECT")) {
              node_iter->type = CONNECT;
              connect_num++;
              success = associate_vartype_ptr(&node_iter->positionPtr.x, z_type->x, connect_num);
              success = associate_vartype_ptr(&node_iter->positionPtr.y, z_type->y, connect_num);
              success = associate_vartype_ptr(&node_iter->positionPtr.z, z_type->z, connect_num);
              success = associate_vartype_ptr(&node_iter->sumForcePtr.fx, other_type->Fx_connect, connect_num);
              success = associate_vartype_ptr(&node_iter->sumForcePtr.fy, other_type->Fy_connect, connect_num);
              success = associate_vartype_ptr(&node_iter->sumForcePtr.fz, other_type->Fz_connect, connect_num);
            } else if (biseqcstrcaseless(parsed->entry[i_parsed],"VESSEL")) {
              node_iter->type = VESSEL;
              vessel_num++;
              success = associate_vartype_ptr(&node_iter->positionPtr.x, u_type->x, vessel_num);
              success = associate_vartype_ptr(&node_iter->positionPtr.y, u_type->y, vessel_num);
              success = associate_vartype_ptr(&node_iter->positionPtr.z, u_type->z, vessel_num);
              success = associate_vartype_ptr(&node_iter->sumForcePtr.fx, y_type->Fx, vessel_num);
              success = associate_vartype_ptr(&node_iter->sumForcePtr.fy, y_type->Fy, vessel_num);
              success = associate_vartype_ptr(&node_iter->sumForcePtr.fz, y_type->Fz, vessel_num);
            } else {
              user_msg = bformat("Value: <%s>", parsed->entry[i_parsed]->data);
              *ierr = map_set_universal_error(user_msg, map_msg, *ierr, MAP_FATAL_25);        
              success = bdestroy(user_msg); /* clear user_msg to not inavertendly picked up elsewhere */
              user_msg = NULL;
            };
            next++;
          } else if (next==2) { /* set initial X node position values */
            alias = bformat("X[%d]", i+1);                          
            success = set_vartype_ptr("[m]", alias, i, &node_iter->positionPtr.x, parsed->entry[i_parsed]); CHECKERRQ(MAP_FATAL_17);
            bdestroy(alias);
            next++;
          } else if (next==3) { /* set initial Y node position values */
            alias = bformat("Y[%d]", i+1);                          
            success = set_vartype_ptr("[m]", alias, i, &node_iter->positionPtr.y, parsed->entry[i_parsed]); CHECKERRQ(MAP_FATAL_18);
            bdestroy(alias);
            next++;
          } else if (next==4) { /* set initial Z node position values */
            alias = bformat("Z[%d]", i+1);                          
            if (biseqcstrcaseless(parsed->entry[i_parsed],"DEPTH")) {          
              if (node_iter->type!=FIX) { /* can only use 'DEPTH' flag in input file for FIX (anchor) nodes */
                user_msg = bformat("Value: <%s>", parsed->entry[i_parsed]->data);
                *ierr = map_set_universal_error(user_msg, map_msg, *ierr, MAP_FATAL_71); 
                success = bdestroy(user_msg); /* clear user_msg to not inavertendly picked up elsewhere */
                user_msg = NULL;
              } else {
                value_string = bformat("%f", depth);                          
                success = set_vartype_ptr("[m]", alias, i, &node_iter->positionPtr.z, value_string); CHECKERRQ(MAP_FATAL_19);
              };
            } else { /* all other nodes not using the 'DEPTH' flag */
              success = set_vartype_ptr("[m]", alias, i, &node_iter->positionPtr.z, parsed->entry[i_parsed]); CHECKERRQ(MAP_FATAL_19);
            };        
            bdestroy(alias);
            next++;
          } else if (next==5) { /* set the node mass */            
            alias = bformat("M[%d]", i+1);                          
            success = set_vartype("[kg]", alias, i, &node_iter->MApplied, parsed->entry[i_parsed]); CHECKERRQ(MAP_FATAL_20);
            bdestroy(alias);
            next++;  
          } else if (next==6) { /* set the node buoyancy */
            alias = bformat("B[%d]", i+1);                          
            success = set_vartype("[m^3]", alias, i, &node_iter->BApplied, parsed->entry[i_parsed]); CHECKERRQ(MAP_FATAL_21);
            bdestroy(alias);
            next++; 
          } else if (next==7) { /* set applied X external force (or user guess) of the node */                    
            alias = bformat("FX[%d]", i+1);                          
            success = set_vartype("[N]", alias, i, &node_iter->externalForce.fx, parsed->entry[i_parsed]); CHECKERRQ(MAP_FATAL_22);
            bdestroy(alias);
            next++;
          } else if (next==8) { /* set applied Y external force (or user guess) of the node */            
            alias = bformat("FY[%d]", i+1);                          
            success = set_vartype("[N]", alias, i, &node_iter->externalForce.fy, parsed->entry[i_parsed]); CHECKERRQ(MAP_FATAL_23);
            bdestroy(alias);
            next++;
          } else if (next==9) { /* set applied Z external force (or user guess) of the node */
            alias = bformat("FZ[%d]", i+1);                          
            success = set_vartype("[N]", alias, i, &node_iter->externalForce.fz, parsed->entry[i_parsed]); CHECKERRQ(MAP_FATAL_24);
            bdestroy(alias);
            next++;
          } else {            
            next++;
          };
        };
        i_parsed++;
      };
      // printf("Node is fixed: %d\n", new_node.type);
      // printf("Name: %s\n", new_node.positionPtr.z.name->data);
      // printf("Name: %f\n", *new_node.positionPtr.z.value);
      // printf("units: %f\n\n", new_node.externalForce.fz.value);
        
      // init_data->expandedNodeInputString->qty++;
      // init_data->expandedNodeInputString->entry[n_line] = bstrcpy(line);
      //success = bassigncstr(line, "");
    } while (0);   
    success = bstrListDestroy(parsed);
    /* @todo: need to make sure next==9; otherwise not enough inputs and an error should
     *        be thrown
     */
    // list_append(&model_data->node, &new_node);    
  };  

  model_data->sizeOfNodes = list_size(&model_data->node);

  /* check to make sure the number of allocated array spaces for fortran derived types matches 
   * what was actually set in the node initialization front end.
   */     
  do {
    success = compare_length(other_type->Fx_connect_Len, connect_num); CHECKERRQ(MAP_FATAL_49);
    success = compare_length(other_type->Fy_connect_Len, connect_num); CHECKERRQ(MAP_FATAL_49);
    success = compare_length(other_type->Fz_connect_Len, connect_num); CHECKERRQ(MAP_FATAL_49);

    success = compare_length(other_type->Fx_anchor_Len, fix_num); CHECKERRQ(MAP_FATAL_49); // @todo: change error code
    success = compare_length(other_type->Fy_anchor_Len, fix_num); CHECKERRQ(MAP_FATAL_49); // @todo: change error code
    success = compare_length(other_type->Fz_anchor_Len, fix_num); CHECKERRQ(MAP_FATAL_49); // @todo: change error code

    success = compare_length(other_type->x_Len, fix_num); CHECKERRQ(MAP_FATAL_49);
    success = compare_length(other_type->y_Len, fix_num); CHECKERRQ(MAP_FATAL_49);
    success = compare_length(other_type->z_Len, fix_num); CHECKERRQ(MAP_FATAL_49);
              
    success = compare_length(u_type->x_Len, vessel_num); CHECKERRQ(MAP_FATAL_50);
    success = compare_length(u_type->y_Len, vessel_num); CHECKERRQ(MAP_FATAL_50);
    success = compare_length(u_type->z_Len, vessel_num); CHECKERRQ(MAP_FATAL_50);    
              
    success = compare_length(y_type->Fx_Len, vessel_num); CHECKERRQ(MAP_FATAL_51);
    success = compare_length(y_type->Fy_Len, vessel_num); CHECKERRQ(MAP_FATAL_51);
    success = compare_length(y_type->Fz_Len, vessel_num); CHECKERRQ(MAP_FATAL_51);    
              
    success = compare_length(z_type->x_Len, connect_num); CHECKERRQ(MAP_FATAL_52);
    success = compare_length(z_type->y_Len, connect_num); CHECKERRQ(MAP_FATAL_52);
    success = compare_length(z_type->z_Len, connect_num); CHECKERRQ(MAP_FATAL_52);    
  } while (0);  

  MAP_RETURN;
};


MAP_ERROR_CODE set_vartype_float(const char* unit, const char* alias, const int num, VarType* type, const double value)
{
  type->name = bfromcstr(alias);
  type->units = bfromcstr(unit); 
  type->referenceCounter = 0;
  type->id = num;
  type->value = value;

  return MAP_SAFE;
};


MAP_ERROR_CODE set_vartype(const char* unit, bstring alias, const int num, VarType* type, bstring property)
{
  type->name = bstrcpy(alias);
  type->units = bfromcstr(unit); 
  type->referenceCounter = 0;
  type->id = num;
  
  if (!property) { /* this option should only be called for setting element vartypes */
    type->value = -999.9;
  } else {
    if (property->data[0]=='#') { /* this variable is an iterated parameter */      
      type->isFixed = false;
      if (property->slen==1) { /* implies that property->data = "#" */
        type->value = -999.9;
      } else if (is_numeric(remove_first_character(property->data))) { 
        type->value = (double)atof(remove_first_character(property->data));
      } else {
        return MAP_FATAL;
      };
    } else { /* this variable is constant */    
      type->isFixed = true;
      if (is_numeric(property->data)) { 
        type->value = (double)atof(property->data);
      } else {
        return MAP_FATAL;
      };
    };
  };
  return MAP_SAFE;
};


MAP_ERROR_CODE set_vartype_ptr(const char* unit, bstring alias, const int num, VarTypePtr* type, bstring property)
{
  type->name = bstrcpy(alias);
  type->units = bfromcstr(unit); 
  type->referenceCounter = 0;
  type->id = num;  

  if (property->data[0]=='#') { /* this variable is an iterated parameter */      
    type->isFixed = false;
    if (property->slen==1) { /* implies that property->data = "#" */
      *type->value = -999.9;
    } else if (is_numeric(remove_first_character(property->data))) { 
      *type->value = (double)atof(remove_first_character(property->data));
    } else {
      return MAP_FATAL;
    };
  } else { /* this variable is constant */    
    type->isFixed = true;
    if (is_numeric(property->data)) { 
      *type->value = (double)atof(property->data);
    } else {
      return MAP_FATAL;
    };
  };
  return MAP_SAFE;
};


MAP_ERROR_CODE set_element_option_flags(bstring word, Element* element_ptr, char* map_msg, MAP_ERROR_CODE* ierr)
{
  MAP_ERROR_CODE success = MAP_SAFE;
  bstring user_msg = NULL;

  if (biseqcstrcaseless(word,"PLOT")) {    
    element_ptr->options.plotFlag = true;
  } else if (biseqcstrcaseless(word, "GX_POS")) {
    element_ptr->options.gxPosFlag = true;
  } else if (biseqcstrcaseless(word, "GY_POS")) {
    element_ptr->options.gyPosFlag = true;
  } else if (biseqcstrcaseless(word, "GZ_POS")) {
    element_ptr->options.gzPosFlag = true;
  } else if (biseqcstrcaseless(word, "GX_A_POS")) {
    element_ptr->options.gxAnchorPosFlag = true;
  } else if (biseqcstrcaseless(word, "GY_A_POS")) {
    element_ptr->options.gyAnchorPosFlag = true;
  } else if (biseqcstrcaseless(word, "GZ_A_POS")) {
    element_ptr->options.gzAnchorPosFlag = true;
  } else if (biseqcstrcaseless(word, "GX_FORCE")) {
    element_ptr->options.gxForceFlag= true;
  } else if (biseqcstrcaseless(word, "GY_FORCE")) {
    element_ptr->options.gyForceFlag = true;
  } else if (biseqcstrcaseless(word, "GZ_FORCE")) {
    element_ptr->options.gzForceFlag = true;
  } else if (biseqcstrcaseless(word, "H_FAIR")) {
    element_ptr->options.HFlag = true;
  } else if (biseqcstrcaseless(word, "H_ANCH")) {
    element_ptr->options.HAnchorFlag = true;
  } else if (biseqcstrcaseless(word, "V_FAIR")) {
    element_ptr->options.VFlag = true;
  } else if (biseqcstrcaseless(word, "V_ANCH")) {
    element_ptr->options.VAnchorFlag = true;
  } else if (biseqcstrcaseless(word, "TENSION_FAIR")) {
    element_ptr->options.fairleadTensionFlag = true;
  } else if (biseqcstrcaseless(word, "TENSION_ANCH")) {
    element_ptr->options.anchorTensionFlag = true;
  } else if (biseqcstrcaseless(word, "X_EXCURSION")) {
    element_ptr->options.horizontalExcursionFlag = true;
  } else if (biseqcstrcaseless(word, "Z_EXCURSION")) {
    element_ptr->options.verticalExcursionFlag = true;
  } else if (biseqcstrcaseless(word, "AZIMUTH")) {
    element_ptr->options.azimuthFlag = true;
  } else if (biseqcstrcaseless(word, "ALTITUDE")) {
    element_ptr->options.altitudeFlag = true;
  } else if (biseqcstrcaseless(word, "ALTITUDE_A")) {
    element_ptr->options.altitudeAnchorFlag = true;
  } else if (biseqcstrcaseless(word, "LINE_TENSION")) {
    element_ptr->options.lineTensionFlag = true;
  } else if (biseqcstrcaseless(word, "OMIT_CONTACT")) {
    element_ptr->options.omitContact = true;
  } else if (biseqcstrcaseless(word, "SEG_SIZE")) {
    if (is_numeric(word->data)) {
      element_ptr->segmentSize = (MapReal)atof(word->data);
    } else { /* should not cancel the simulation; simply ignore it */      
      user_msg = bformat("Option <%s>", word->data);
      *ierr = map_set_universal_error(user_msg, map_msg, *ierr, MAP_ERROR_18);
      success = bdestroy(user_msg); /* clear user_msg to not inavertendly picked up elsewhere */
      user_msg = NULL;
      // *ierr = map_set_universal_error(buffer, map_msg, ierr, MAP_ERROR_18); 
    };          
  } else if (biseqcstrcaseless(word, "LAY_LENGTH")) {
    element_ptr->options.layLengthFlag = true;
  // } else if (biseqcstrcaseless("DAMAGE_TIME", word)) {
  //   word = strtok(NULL, " ,\n\t\r\0");
  //   success = is_numeric(word);
  //   if (success==MAP_SAFE) {
  //     elem->options.damageTimeFlag = true;
  //     elem->damageTime = (MapReal)atof(word);
  //   } else { /* should not cancel the simulation; simply ignore it */     
  // *ierr=map_set_universal_error(buffer, map_msg, ierr, MAP_ERROR_1);
  //   };          
  // } else if (biseqcstrcaseless("DIAGNOSTIC", word)) {
  //   word = strtok(NULL, " ,\n\t\r\0");
  //   success = is_numeric(word);
  //   if (success==MAP_SAFE) {
  //     elem->options.diagnosticsFlag = true;
  //     elem->diagnosticType = (int)atoi(word);
  //   } else {
  //     /* should not cancel the simulation; simply ignore it */
  //     *ierr = map_set_universal_error(buffer, map_msg, ierr, MAP_ERROR_14);
  //     elem->options.diagnosticsFlag = true;
  //     elem->diagnosticType = 0;
  //   };          
  } else {
    /* should not cancel the simulation; simply ignore it */
    user_msg = bformat("Option <%s>", word->data);
    *ierr = map_set_universal_error(user_msg, map_msg, *ierr, MAP_WARNING_3);
    success = bdestroy(user_msg); /* clear user_msg to not inavertendly picked up elsewhere */
    user_msg = NULL;
  };
  return MAP_SAFE;
};


MAP_ERROR_CODE set_element_list(MAP_ConstraintStateType_t* z_type, ModelData* model_data, struct bstrList* element_input_string, char* map_msg, MAP_ERROR_CODE* ierr)
{
  MAP_ERROR_CODE success = MAP_SAFE;
  int i = 0;
  int i_parsed = 0;
  int next = 0;
  const int num_elements = element_input_string->qty;
  Element new_element;
  Element* element_iter = NULL;
  struct bstrList* parsed = NULL;
  struct tagbstring tokens; 
  bstring user_msg = NULL;
  bstring alias = NULL;

  cstr2tbstr(tokens," \t\n\r"); /* token for splitting line into indivdual words is a tab and space */   
  success = reset_element(&new_element);

  model_data->sizeOfElements = num_elements;

  z_type->H_Len = model_data->sizeOfElements;          
  z_type->V_Len = model_data->sizeOfElements;          
  z_type->H = (double*)malloc(z_type->H_Len*sizeof(double));
  z_type->V = (double*)malloc(z_type->V_Len*sizeof(double));

  if (z_type->H==NULL || z_type->V==NULL) {
    user_msg = bformat("Failed allocation of a z_type");
    *ierr = map_set_universal_error(user_msg, map_msg, *ierr, MAP_FATAL_53);
    success = bdestroy(user_msg); /* clear user_msg to not inavertendly picked up elsewhere */
    user_msg = NULL;
    return MAP_FATAL;
  };
  
  for(i=0 ; i<num_elements ; i++) {         
    list_append(&model_data->element, &new_element);
    element_iter = (Element*)list_get_at(&model_data->element, i);
    success = set_element_vartype(element_iter); /* @todo: check error */

    i_parsed = 0;
    next = 0;
    parsed = bsplits(element_input_string->entry[i], &tokens);
    do {  
      while (i_parsed<parsed->qty-1) { /* iterating through all strings */              
        if (parsed->entry[i_parsed]->slen) { /* if the string length is not 0 */
          if (next==0) { /* use this first option as an opportunity to set the run-time flags to false */             
            success = associate_vartype_ptr(&element_iter->H, z_type->H, i+1);
            success = associate_vartype_ptr(&element_iter->V, z_type->V, i+1);

            element_iter->H.isFixed = false;
            alias = bformat("H[%d]", i+1);
            success = set_vartype_ptr("[N]", alias, i, &element_iter->H, parsed->entry[i_parsed]); CHECKERRQ(MAP_FATAL_32);            
            success = bdestroy(alias);

            element_iter->V.isFixed = false;
            alias = bformat("V[%d]", i+1);
            success = set_vartype_ptr("[N]", alias, i, &element_iter->V, parsed->entry[i_parsed]); CHECKERRQ(MAP_FATAL_32);                        
            success = bdestroy(alias);
             
            next++;
          } else if (next==1) {
            success = associate_element_with_cable_property(element_iter, model_data, parsed->entry[i_parsed]->data, map_msg, ierr); CHECKERRQ(MAP_FATAL_32);           
            next++;
          } else if (next==2) { 
            alias = bformat("Lu[%d]", i+1);
            success = set_vartype("[m]", alias, i, &element_iter->Lu, parsed->entry[i_parsed]); CHECKERRQ(MAP_FATAL_26);
            success = bdestroy(alias);
            next++;
          } else if (next==3) { 
            success = associate_element_with_anchor_node(element_iter, model_data, i+1, parsed->entry[i_parsed]->data,  map_msg, ierr); CHECKERRQ(MAP_FATAL_32);        
            next++;
          } else if (next==4) { 
            success = associate_element_with_fairlead_node(element_iter, model_data, i+1, parsed->entry[i_parsed]->data,  map_msg, ierr); CHECKERRQ(MAP_FATAL_32);        
            next++;
          } else { /* set the node mass */            
            success = set_element_option_flags(parsed->entry[i_parsed], element_iter, map_msg, ierr);
            // next++;  
          };
        };
        i_parsed++;
      };
      // printf("Node is fixed: %d\n", new_node.type);
      // printf("Name: %s\n", new_node.positionPtr.z.name->data);
      // printf("Name: %f\n", *new_node.positionPtr.z.value);
      // printf("units: %f\n\n", new_node.externalForce.fz.value);
        
      // init_data->expandedNodeInputString->qty++;
      // init_data->expandedNodeInputString->entry[n_line] = bstrcpy(line);
      //success = bassigncstr(line, "");
    } while (0);   
    success = bstrListDestroy(parsed);
  };  





//   int i = 0;
//   int cx = 0;
//   int cnt = 0;
//   int success = 0;
//   int sizeOfString = 0;
//   char buffer[64] = "";
//   char* line = NULL;
//   char* word = NULL;
//   Element newElement;
//   Element* elementIter = NULL;
//   int lineCounter = 0;
// 
//   zType->H_Len = data->sizeOfElements;          
//   zType->V_Len = data->sizeOfElements;          
//   zType->H = (double*)malloc(zType->H_Len*sizeof(double));
//   zType->V = (double*)malloc(zType->V_Len*sizeof(double));
// 
//   if (zType->H==NULL) {
//     *ierr = map_set_universal_error(buffer, map_msg, ierr, MAP_FATAL_53);
//     return MAP_FATAL;
//   };
// 
//   if (zType->V==NULL) {
//     *ierr = map_set_universal_error(buffer, map_msg, ierr, MAP_FATAL_53);
//     return MAP_FATAL;
//   };
// 
//   for(i=0 ; i<=data->sizeOfElements-1 ; i++) {
//     sizeOfString = strlen(elementInputString[i]);
//     line = (char*)realloc(line, (sizeOfString+1)*sizeof(char));    
//     strcpy(line, elementInputString[i]);    
//     word = strtok(line, " ,\n\t\r\0");
//     while (word!=NULL) {
//       do {
//         if (cnt==0) {
//           /* use this first option as an opportunity to set the run-time flags to false */
//           initialize_element(&newElement, i+1);
//           success = associate_vartype_ptr(&newElement.H, zType->H, i+1);
//           success = associate_vartype_ptr(&newElement.V, zType->V, i+1);          
//           cx = map_snprintf(buffer, 64, "V[%d]", i+1);assert(cx>=0);                  
//           success = set_vartype_ptr("[N]", buffer, i, &newElement.V, word); CHECKERRQ(MAP_FATAL_32);
//           cx = map_snprintf(buffer, 64, "H[%d]", i+1); assert(cx>=0);                            
//           success = set_vartype_ptr("[N]", buffer, i, &newElement.H, word); CHECKERRQ(MAP_FATAL_32);
//         } else if (cnt==1) {
//           success = associate_element_with_cable_property(&newElement, data, word, map_msg, ierr); CHECKERRQ(MAP_FATAL_32);
//         } else if (cnt==2) { 
//           cx = map_snprintf(buffer, 64, "Lu[%d]", i+1); assert(cx>=0);                
//           success = set_vartype("[m]", buffer, i, &newElement.Lu, word); CHECKERRQ(MAP_FATAL_26);
//         } else if (cnt==3) {
//           success = associate_element_with_anchor_node(data, &newElement, i+1, word, map_msg, ierr); CHECKERRQ(MAP_FATAL_32);        
//         } else if (cnt==4) {
//           success = associate_element_with_fairlead_node(data, &newElement, i+1, word, map_msg, ierr); CHECKERRQ(MAP_FATAL_32); 
//         } else { 
//           success = map_set_option_flags (word, &newElement, map_msg, ierr);
//         };
//       } while (0);
//       cnt++;
//       word = strtok(NULL, " ,\n\t\r\0");
//     };  
//     list_append(&data->element, &newElement);
//     cnt=0;
//   };
// 
//   /* check EA, Lu errors. These values cannot be negative. */
//   list_iterator_start(&data->element);            /* starting an iteration "session" */
//   while (list_iterator_hasnext(&data->element)) { /* tell whether more values available */ 
//     elementIter = (Element*)list_iterator_next(&data->element);
//     if(elementIter->Lu.value<=0.0) {
//       success = MAP_FATAL;
//       cx = map_snprintf(buffer, 64, "Line segment %d, Lu = %f [m].", lineCounter, elementIter->Lu.value); assert(cx>=0);
//       *ierr = map_set_universal_error(buffer, map_msg, ierr, MAP_FATAL_56);      
//     } else if(elementIter->lineProperty->ea<=0.0) {
//       success = MAP_FATAL;
//       cx = map_snprintf(buffer, 64, "Line segment %d, EA = %f [N].", lineCounter, elementIter->lineProperty->ea); assert(cx>=0);
//       *ierr = map_set_universal_error(buffer, map_msg, ierr, MAP_FATAL_57);      
//     };
//     lineCounter++;
//   };
//   list_iterator_stop(&data->element); /* ending the iteration "session" */    
// 
//   MAPFREE(line);
// 
//   if (*ierr==MAP_FATAL) {
//     return MAP_FATAL;
//   };
//   return MAP_SAFE;
  return MAP_SAFE;
};



MAP_ERROR_CODE set_output_list(ModelData* model_data, MAP_InitOutputType_t* io_type, char* map_msg, MAP_ERROR_CODE* ierr)
{
  Element* element_iter = NULL;
  OutputList* yList = model_data->yList;
  
  list_iterator_start(&model_data->element); /* starting an iteration "session" */
  while (list_iterator_hasnext(&model_data->element)) { /* tell whether more values available */ 
    element_iter = (Element*)list_iterator_next(&model_data->element);    
    
    if (element_iter->options.gxAnchorPosFlag) {
      list_append(&yList->out_list_ptr, &element_iter->anchor->positionPtr.x);      
      io_type->writeOutputHdr_Len++;
      io_type->writeOutputUnt_Len++;
    };

    if (element_iter->options.gyAnchorPosFlag) {
      list_append(&yList->out_list_ptr, &element_iter->anchor->positionPtr.y);
      io_type->writeOutputHdr_Len++;
      io_type->writeOutputUnt_Len++;
    };

    if (element_iter->options.gzAnchorPosFlag) {
      list_append(&yList->out_list_ptr, &element_iter->anchor->positionPtr.z);
      io_type->writeOutputHdr_Len++;
      io_type->writeOutputUnt_Len++;
    };

    if (element_iter->options.gxPosFlag) {
      list_append(&yList->out_list_ptr, &element_iter->fairlead->positionPtr.x);
      io_type->writeOutputHdr_Len++;
      io_type->writeOutputUnt_Len++;
    };

    if (element_iter->options.gyPosFlag) {
      list_append(&yList->out_list_ptr, &element_iter->fairlead->positionPtr.y);
      io_type->writeOutputHdr_Len++;
      io_type->writeOutputUnt_Len++;
    };

    if (element_iter->options.gzPosFlag) {
      list_append(&yList->out_list_ptr, &element_iter->fairlead->positionPtr.z);
      io_type->writeOutputHdr_Len++;
      io_type->writeOutputUnt_Len++;
    };

    if (element_iter->options.HFlag) {
      list_append(&yList->out_list_ptr, &element_iter->H);
      io_type->writeOutputHdr_Len++;
      io_type->writeOutputUnt_Len++;
    };

    if (element_iter->options.VFlag) {
      list_append(&yList->out_list_ptr, &element_iter->V);
      io_type->writeOutputHdr_Len++;
      io_type->writeOutputUnt_Len++;
    };
    
    if (element_iter->options.HAnchorFlag) {
      list_append(&yList->out_list, &element_iter->HAtAnchor);
      io_type->writeOutputHdr_Len++;
      io_type->writeOutputUnt_Len++;
    };

    if (element_iter->options.VAnchorFlag) {
      list_append(&yList->out_list, &element_iter->VAtAnchor);
      io_type->writeOutputHdr_Len++;
      io_type->writeOutputUnt_Len++;
    };

    if (element_iter->options.anchorTensionFlag) {
      list_append(&yList->out_list, &element_iter->TAtAnchor);
      io_type->writeOutputHdr_Len++;
      io_type->writeOutputUnt_Len++;
    };

    if (element_iter->options.altitudeAnchorFlag) {
      list_append(&yList->out_list, &element_iter->alphaAtAnchor);
      io_type->writeOutputHdr_Len++;
      io_type->writeOutputUnt_Len++;
    };

    if (element_iter->options.gxForceFlag) {
      list_append(&yList->out_list, &element_iter->forceAtFairlead.fx); /* @todo: this is not correct. Should point to fairlead->sumForce.fx */
      io_type->writeOutputHdr_Len++;
      io_type->writeOutputUnt_Len++;
    };

    if (element_iter->options.gyForceFlag) {
      list_append(&yList->out_list, &element_iter->forceAtFairlead.fy); /* @todo: this is not correct. Should point to fairlead->sumForce.fy */
      io_type->writeOutputHdr_Len++;
      io_type->writeOutputUnt_Len++;
    };

    if (element_iter->options.gzForceFlag) {
      list_append(&yList->out_list, &element_iter->forceAtFairlead.fz); /* @todo: this is not correct. Should point to fairlead->sumForce.fz */
      io_type->writeOutputHdr_Len++;
      io_type->writeOutputUnt_Len++;
    };

    if (element_iter->options.VFlag) {
      list_append(&yList->out_list, &element_iter->forceAtFairlead.fz); /* @todo: this is not correct. Doubled up with above */
      io_type->writeOutputHdr_Len++;
      io_type->writeOutputUnt_Len++;
    };

    if (element_iter->options.fairleadTensionFlag) {
      list_append(&yList->out_list, &element_iter->T);
      io_type->writeOutputHdr_Len++;
      io_type->writeOutputUnt_Len++;
    };

    if (element_iter->options.horizontalExcursionFlag) {
      list_append(&yList->out_list, &element_iter->l);
      io_type->writeOutputHdr_Len++;
      io_type->writeOutputUnt_Len++;
    };

    if (element_iter->options.verticalExcursionFlag) {
      list_append(&yList->out_list, &element_iter->h);
      io_type->writeOutputHdr_Len++;
      io_type->writeOutputUnt_Len++;
    };

    if (element_iter->options.layLengthFlag) {
      list_append(&yList->out_list, &element_iter->lb);
      io_type->writeOutputHdr_Len++;
      io_type->writeOutputUnt_Len++;
    };

    if (element_iter->options.azimuthFlag) {
      list_append(&yList->out_list, &element_iter->psi);
      io_type->writeOutputHdr_Len++;
      io_type->writeOutputUnt_Len++;
    };

    if (element_iter->options.altitudeFlag) {
      list_append(&yList->out_list, &element_iter->alpha);
      io_type->writeOutputHdr_Len++;
      io_type->writeOutputUnt_Len++;
    };
  };
  list_iterator_stop(&model_data->element); /* ending the iteration session */  

  return MAP_SAFE;
};


/**
 * this function returns the size of elements 
 */
size_t vartype_meter(const void* el) 
{
  /* every element has the constant size of a rectangle structure */
  return sizeof(VarType);
};


/**
 * this function returns the size of elements 
 */
size_t vartype_ptr_meter(const void* el) 
{
  /* every element has the constant size of a rectangle structure */
  return sizeof(VarTypePtr);
};


MAP_ERROR_CODE reset_element(Element* element_ptr)
{
  element_ptr->event.LuRestore = 0.0;
  element_ptr->event.dLu = 0.0;
  element_ptr->event.LuMax = 0.0;
  element_ptr->event.payinFlag = false;

  /* run-time flags */
  element_ptr->options.plotFlag = false;
  element_ptr->options.gxPosFlag = false;
  element_ptr->options.gyPosFlag = false;
  element_ptr->options.gzPosFlag = false;
  element_ptr->options.gxAnchorPosFlag = false;
  element_ptr->options.gyAnchorPosFlag = false;
  element_ptr->options.gzAnchorPosFlag = false;
  element_ptr->options.gxForceFlag = false;
  element_ptr->options.gyForceFlag = false;
  element_ptr->options.gzForceFlag = false;
  element_ptr->options.HFlag = false;
  element_ptr->options.HAnchorFlag = false;
  element_ptr->options.VFlag = false;
  element_ptr->options.VAnchorFlag = false;
  element_ptr->options.fairleadTensionFlag = false;
  element_ptr->options.anchorTensionFlag = false;
  element_ptr->options.horizontalExcursionFlag = false;
  element_ptr->options.verticalExcursionFlag = false;
  element_ptr->options.azimuthFlag = false;
  element_ptr->options.altitudeFlag = false;
  element_ptr->options.altitudeAnchorFlag = false;
  element_ptr->options.lineTensionFlag = false;
  element_ptr->options.omitContact = false;
  element_ptr->options.layLengthFlag = false;
  element_ptr->options.damageTimeFlag = false;
  element_ptr->options.diagnosticsFlag = false;
  
  element_ptr->lineProperty = NULL;      
  element_ptr->label = NULL;
  element_ptr->lineTension  = NULL;
  element_ptr->anchor = NULL;             /* Anchor node */
  element_ptr->fairlead = NULL;           /* Fairlead node */
  
  element_ptr->psi.name = NULL;
  element_ptr->psi.units = NULL;
  element_ptr->alpha.name = NULL;
  element_ptr->alpha.units = NULL;
  element_ptr->alphaAtAnchor.name = NULL;
  element_ptr->alphaAtAnchor.units = NULL;
  element_ptr->l.name = NULL;
  element_ptr->l.units = NULL;
  element_ptr->lb.name = NULL;
  element_ptr->lb.units = NULL;
  element_ptr->h.name = NULL;
  element_ptr->h.units = NULL;
  element_ptr->H.name = NULL;
  element_ptr->H.units = NULL;
  element_ptr->V.name = NULL;
  element_ptr->V.units = NULL;
  element_ptr->HAtAnchor.name = NULL;
  element_ptr->HAtAnchor.units = NULL;
  element_ptr->VAtAnchor.name = NULL;
  element_ptr->VAtAnchor.units = NULL;
  element_ptr->forceAtFairlead.fx.name = NULL;
  element_ptr->forceAtFairlead.fx.units = NULL;
  element_ptr->forceAtFairlead.fy.name = NULL;
  element_ptr->forceAtFairlead.fy.units = NULL;
  element_ptr->forceAtFairlead.fz.name = NULL;
  element_ptr->forceAtFairlead.fz.units = NULL;
  element_ptr->forceAtAnchor.fx.name = NULL;
  element_ptr->forceAtAnchor.fx.units = NULL;
  element_ptr->forceAtAnchor.fy.name = NULL;
  element_ptr->forceAtAnchor.fy.units = NULL;
  element_ptr->forceAtAnchor.fz.name = NULL;
  element_ptr->forceAtAnchor.fz.units = NULL;
  element_ptr->T.name = NULL;
  element_ptr->T.units = NULL;
  element_ptr->TAtAnchor.name = NULL;
  element_ptr->TAtAnchor.units = NULL;

  element_ptr->damageTime = -999.9;
  element_ptr->diagnosticType = -9999;
  element_ptr->segmentSize = 10;

};


MAP_ERROR_CODE set_element_vartype(Element* element_ptr)
{
  MAP_ERROR_CODE success = MAP_SAFE;
  bstring alias = NULL;

  alias = bformat("psi");                          
  success = set_vartype("[deg]", alias, 0, &element_ptr->psi, NULL); /* @todo: check error */
  bdestroy(alias);

  alias = bformat("alpha");                          
  success = set_vartype("[deg]", alias, 0, &element_ptr->alpha, NULL); /* @todo: check error */
  bdestroy(alias);

  alias = bformat("alpha_a");                          
  success = set_vartype("[deg]", alias, 0, &element_ptr->alphaAtAnchor, NULL); /* @todo: check error */
  bdestroy(alias);

  alias = bformat("l");                          
  success = set_vartype("[m]", alias, 0, &element_ptr->l, NULL); /* @todo: check error */
  bdestroy(alias);

  alias = bformat("lb");                          
  success = set_vartype("[m]", alias, 0, &element_ptr->lb, NULL); /* @todo: check error */
  bdestroy(alias);

  alias = bformat("h");                          
  success = set_vartype("[m]", alias, 0, &element_ptr->h, NULL); /* @todo: check error */
  bdestroy(alias);

  alias = bformat("H_a");                          
  success = set_vartype("[N]", alias, 0, &element_ptr->HAtAnchor, NULL); /* @todo: check error */
  bdestroy(alias);

  alias = bformat("fx");                          
  success = set_vartype("[N]", alias, 0, &element_ptr->forceAtFairlead.fx, NULL); /* @todo: check error */
  bdestroy(alias);

  alias = bformat("fy");                          
  success = set_vartype("[N]", alias, 0, &element_ptr->forceAtFairlead.fy, NULL); /* @todo: check error */
  bdestroy(alias);

  alias = bformat("fz");                          
  success = set_vartype("[N]", alias, 0, &element_ptr->forceAtFairlead.fz, NULL); /* @todo: check error */
  bdestroy(alias);

  alias = bformat("fx_a");                          
  success = set_vartype("[N]", alias, 0, &element_ptr->forceAtAnchor.fx, NULL); /* @todo: check error */
  bdestroy(alias);

  alias = bformat("fy_a");                          
  success = set_vartype("[N]", alias, 0, &element_ptr->forceAtAnchor.fy, NULL); /* @todo: check error */
  bdestroy(alias);

  alias = bformat("fz_a");                          
  success = set_vartype("[N]", alias, 0, &element_ptr->forceAtAnchor.fz, NULL); /* @todo: check error */
  bdestroy(alias);

  alias = bformat("V_a");                          
  success = set_vartype("[N]", alias, 0, &element_ptr->VAtAnchor, NULL); /* @todo: check error */
  bdestroy(alias);

  alias = bformat("T");                          
  success = set_vartype("[N]", alias, 0, &element_ptr->T, NULL); /* @todo: check error */
  bdestroy(alias);

  alias = bformat("T_a");                          
  success = set_vartype("[N]", alias, 0, &element_ptr->TAtAnchor, NULL); /* @todo: check error */
  bdestroy(alias);

  return MAP_SAFE;
};


// MAP_ERROR_CODE initialize_external_applied_force(char* unit, char* alias, const int num, VarType* type, char const* property)
// {
//   type->name = malloc(sizeof(char)*strlen(alias)+1);
//   type->units = malloc(sizeof(char)*strlen(unit)+1);
//   
//   strcpy(type->name, alias);
//   strcpy(type->units, unit);
// 
//   type->referenceCounter = 0;
//   type->id = num;  
//   type->value = 0;
//   type->isFixed = true;
// 
//   if (property[0]=='#') { 
//     type->value = 0.0;
//   } else { 
//     /* this variable is constant */    
//     if (is_numeric(property)) { 
//       /* cannot convert to numeric value */
//       return MAP_FATAL;
//     } else { 
//       /* converted to numeric value */
//       type->value = (MapReal)atof(property);      
//     };
//   };
//   return MAP_SAFE;
// };


// /**
//  *
//  */
// MAP_ERROR_CODE initialize_node_sum_force_ptr(char* unit, char* alias, const int num, VarTypePtr* type)
// {
//   type->name = malloc(sizeof(char)*strlen(alias)+1);
//   type->units = malloc(sizeof(char)*strlen(unit)+1);  
//   strcpy(type->name, alias);
//   strcpy(type->units, unit);
// 
//   type->referenceCounter = 0;
//   type->id = num;  
//   type->isFixed = true;
//   *(type->value) = 0.0;
//   return MAP_SAFE;
// };


MAP_ERROR_CODE reset_node(Node* node_ptr)
{
  node_ptr->positionPtr.x.name = NULL;
  node_ptr->positionPtr.x.units = NULL;
  node_ptr->positionPtr.x.value = NULL;
  node_ptr->positionPtr.y.name = NULL;
  node_ptr->positionPtr.y.units = NULL;
  node_ptr->positionPtr.y.value = NULL;
  node_ptr->positionPtr.z.name = NULL;
  node_ptr->positionPtr.z.units = NULL;
  node_ptr->positionPtr.z.value = NULL;
  node_ptr->MApplied.name = NULL;
  node_ptr->MApplied.units = NULL;
  node_ptr->BApplied.name = NULL;
  node_ptr->BApplied.units = NULL;
  node_ptr->sumForcePtr.fx.name = NULL;
  node_ptr->sumForcePtr.fx.units = NULL;
  node_ptr->sumForcePtr.fx.value = NULL;
  node_ptr->sumForcePtr.fy.name = NULL;
  node_ptr->sumForcePtr.fy.units = NULL;
  node_ptr->sumForcePtr.fy.value = NULL;
  node_ptr->sumForcePtr.fz.name = NULL;
  node_ptr->sumForcePtr.fz.units = NULL;
  node_ptr->sumForcePtr.fz.value = NULL;
 
  node_ptr->externalForce.fx.name = NULL;
  node_ptr->externalForce.fx.units = NULL;
  node_ptr->externalForce.fy.name = NULL;
  node_ptr->externalForce.fy.units = NULL;
  node_ptr->externalForce.fz.name = NULL;
  node_ptr->externalForce.fz.units = NULL;
};


MAP_ERROR_CODE associate_element_with_cable_property(Element* element_ptr, ModelData* model_data, const char* word, char* map_msg, MAP_ERROR_CODE* ierr)
{
  MAP_ERROR_CODE success = MAP_SAFE;
  CableLibrary* library_iterator = NULL;
  bstring user_msg = NULL;

  library_iterator = NULL;
  element_ptr->lineProperty = NULL;

  list_iterator_start(&model_data->cableLibrary); /* starting an iteration session */
  while (list_iterator_hasnext(&model_data->cableLibrary)) { /* tell whether more values available */
    library_iterator = (CableLibrary*)list_iterator_next(&model_data->cableLibrary);
    if (biseqcstrcaseless(library_iterator->label, word)) {      
      element_ptr->lineProperty = library_iterator;
      list_iterator_stop(&model_data->cableLibrary); /* ending the iteration session */  
      break;
    }; 
  };
  list_iterator_stop(&model_data->cableLibrary); /* ending the iteration session */  
  if (element_ptr->lineProperty==NULL) {        
    user_msg = bformat("No libraries match <%s>.", word);
    *ierr = map_set_universal_error(user_msg, map_msg, *ierr, MAP_FATAL_27);
    success = bdestroy(user_msg); 
    user_msg = NULL;
    return MAP_FATAL;
  };
  return MAP_SAFE;
};


MAP_ERROR_CODE associate_element_with_anchor_node(Element* element_ptr, ModelData* model_data, const int element_num, const char* word, char* map_msg, MAP_ERROR_CODE* ierr)
{
  MAP_ERROR_CODE success = MAP_SAFE;
  Node* node_iter = NULL;
  int node_num = 0;
  bstring user_msg = NULL;

  element_ptr->anchor = NULL;
  
  if (is_numeric(word)) {
    node_num = (int)atoi(word); 
    node_iter = (Node*)list_get_at(&model_data->node, node_num-1);
    element_ptr->anchor = node_iter; /* create the associate with anchor here */
    if (!node_iter) {
      user_msg = bformat("Element %d.", element_num);
      *ierr = map_set_universal_error(user_msg, map_msg, *ierr, MAP_FATAL_30);        
      success = bdestroy(user_msg); /* clear user_msg to not inavertendly picked up elsewhere */
      user_msg = NULL;
      return MAP_FATAL;
    };
  } else {
    user_msg = bformat("Element %d.", element_num);
    *ierr = map_set_universal_error(user_msg, map_msg, *ierr, MAP_FATAL_28);        
    success = bdestroy(user_msg); 
    user_msg = NULL;
    return MAP_FATAL;
  };
  return MAP_SAFE;
};


MAP_ERROR_CODE associate_element_with_fairlead_node(Element* element_ptr, ModelData* model_data, const int element_num, const char* word, char* map_msg, MAP_ERROR_CODE* ierr)
{
  Node* node_iter = NULL;
  int node_num = 0;
  MAP_ERROR_CODE success = MAP_SAFE;
  bstring user_msg = NULL;

  element_ptr->fairlead = NULL;

  if (is_numeric(word)) {
    node_num = (int)atoi(word); 
    node_iter = (Node*)list_get_at(&model_data->node, node_num-1);
    element_ptr->fairlead = node_iter; /* create the associate with anchor here */
    // printf("element type %s\n", element_ptr->fairlead->positionPtr.x.name->data);
    if (!node_iter) {
      user_msg = bformat("Element %d.", element_num);
      *ierr = map_set_universal_error(user_msg, map_msg, *ierr, MAP_FATAL_31);        
      success = bdestroy(user_msg); /* clear user_msg to not inavertendly picked up elsewhere */
      user_msg = NULL;
      return MAP_FATAL;
    };
  } else {
    user_msg = bformat("Element %d.", element_num);
    *ierr = map_set_universal_error(user_msg, map_msg, *ierr, MAP_FATAL_29);        
    success = bdestroy(user_msg); 
    user_msg = NULL;
    return MAP_FATAL;
  };
  return MAP_SAFE;

};


MAP_ERROR_CODE is_numeric(const char* string)
{
  char* p = NULL;
  if (string==NULL || *string=='\0' || isspace(*string)) {
    return MAP_SAFE;
  };
  strtod (string, &p);
  if (*p=='\0') {
    return MAP_FATAL;
  } else {
    return MAP_SAFE;
  };
};


MAP_ERROR_CODE associate_vartype_ptr(VarTypePtr* type, double* arr, int index)
{
  type->value = &arr[index-1];
  return MAP_SAFE;
};


void copy_target_string(char* target, unsigned char* source)
{
  while (*source) {
    *target = *source;
    source++;
    target++;
  };
  *target = '\0';
};


MAP_ERROR_CODE map_get_version(MAP_InitOutputType_t* io_type)
{
  bstring out_string = NULL;
  int ret = 0;

  /* first set the program version defined in the mapsys.h header file 
   * @todo: program version should be tied to the gi revision numner
   */
  out_string = bformat("<%s>",PROGVERSION);
  if (out_string->slen>MAX_INIT_VERSION_STRING_LENGTH) { /* overflow */
    return MAP_FATAL; /* @todo: give proper error code */
  };
  copy_target_string(io_type->version, out_string->data);
  ret = bdestroy(out_string);

  /* the set the compiling date. This is #defined in the mapsys.h header */
  out_string = bformat("<%c%c%c-%c%c-%c%c%c%c>",BUILD_MONTH_CH0,BUILD_MONTH_CH1,BUILD_MONTH_CH2,BUILD_DAY_CH0,BUILD_DAY_CH1,BUILD_YEAR_CH0,BUILD_YEAR_CH1,BUILD_YEAR_CH2,BUILD_YEAR_CH3);
  if (out_string->slen>MAX_INIT_COMPILING_DATA_STRING_LENGTH) { /* overflow */
    return MAP_FATAL; /* @todo: give proper error code */
  };
  copy_target_string(io_type->compilingData, out_string->data);
  ret = bdestroy(out_string);
  return MAP_SAFE;
};


void print_machine_name_to_screen( ) {
  // __get_machine_name(name);
  printf( "%s Ver. %s ", PROGNAME, PROGVERSION); 
  printf( "%c",BUILD_MONTH_CH0 );// build month
  printf( "%c",BUILD_MONTH_CH1 );
  printf( "%c",BUILD_MONTH_CH2 );
  printf( "-" );
  printf( "%c",BUILD_DAY_CH0 );// build day
  printf( "%c",BUILD_DAY_CH1 );
  printf( "-" );
  printf( "%c",BUILD_YEAR_CH0 ); // build year 
  printf( "%c",BUILD_YEAR_CH1 );
  printf( "%c",BUILD_YEAR_CH2 );
  printf( "%c\n",BUILD_YEAR_CH3 );
}


// /**
//  * @see: http://stackoverflow.com/questions/504810/how-do-i-find-the-current-machines-full-hostname-in-c-hostname-and-domain-info
//  */
// void __get_machine_name(char* machineName)
// {
//   char name[150];  
//   
// #if defined(_WIN32) || defined(_WIN64)
//   int i = 0;
//   TCHAR infoBuf[150];
//   DWORD bufCharCount = 150;
//   memset(name, 0, 150);
//   if (GetComputerName(infoBuf, &bufCharCount)) {
//     for (i=0 ; i<150 ; i++) {
//       name[i] = infoBuf[i];
//     };
//   } else {
//     strcpy(name, "Unknown_Host_Name");
//   };
// #else
//   memset(name, 0, 150);
//   gethostname(name, 150);
// #endif
//   strncpy(machineName ,name, 150);
// };


MAP_ERROR_CODE print_help_to_screen()
{
  print_machine_name_to_screen( );

  printf("MAP Input file section defitions:\n");
  printf("  Line dictionary definitions:\n");   
  printf("    -LineType, --User-defined name of line [-]  \n");   
  printf("    -Diam,     --Line diamater, used to calulate area and line displacement per unit length [m]  \n");   
  printf("    -MassDen,  --Mass (in air) per unit length [kg/m]  \n");   
  printf("    -EA,       --Axial stiffness [N] \n");   
  printf("    -CB,       --Cable/seabed Coulumb friction coefficient [-]  \n");   
  printf("    -CIntDamp, --Internal structural damping coefficient [Pa-s]  \n");   
  printf("    -Ca,       --Cross-flow added-mass coeficient [-]\n");   
  printf("    -Cdn,      --Cross-flow drag coefficient [-]\n");   
  printf("    -Cdt,      --Tangent (skin) drag coefficient[-]\n");   
  printf("  Node property definitions:\n");
  printf("    -Node,     --Node number; first starts at 1 [-]\n");
  printf("    -Type,     --Type of node. Must be one of: VESSEL, FIX, CONNECT [-]\n");
  printf("    -X,        --Node X position. '#' must prefix CONNECT nodes; constitutes user initial guess [m]\n");
  printf("    -Y,        --Node Y position. '#' must prefix CONNECT nodes; constitutes user initial guess [m]\n");
  printf("    -Z,        --Node Z position. '#' must prefix CONNECT nodes; constitutes user initial guess [m]\n");
  printf("    -M,        --Applied point mass at node [kg]\n");
  printf("    -B,        --Applied point buoyancy module at node [m^3]\n");  
  printf("    -FX,       --Applied X external force at node. '#' must prefix VESSEL and FIX nodes [N]\n");
  printf("    -FY,       --Applied Y external force at node. '#' must prefix VESSEL and FIX nodes [N]\n");
  printf("    -FZ,       --Applied Z external force at node. '#' must prefix VESSEL and FIX nodes [N]\n");
  printf("  Element property definitions:\n");
  printf("    -Element,  --Element number; first starts at 1 [-]\n");
  printf("    -LineType, --Must match property defined in 'Line Dictions'[-]\n");
  printf("    -UnstrLen, --Unstretched line length [m]\n");
  printf("    -NodeAnch, --Anchor node number corresponding to 'Node Property Definitions' section [-]\n");
  printf("    -NodeFair, --Fairlead node number corresponding to 'Node Property Definitions' section [-]\n");
  printf("    -Flags,    --User run-time flag; see below [-]\n");

  printf("    \n");
  printf("  Element run-time options definitions\n");
  printf("    Outputs:\n");
  printf("      -gx_pos,       --Fairlead posiion in global X [m]\n");
  printf("      -gy_pos,       --Fairlead posiion in global Y [m]\n");
  printf("      -gx_pos,       --Fairlead posiion in global Z [m]\n");
  printf("      -gx_a_pos,     --Anchor posiion in global X [m]\n");
  printf("      -gy_a_pos,     --Anchor posiion in global Y [m]\n");
  printf("      -gz_a_pos,     --Anchor posiion in global Z [m]\n");
  printf("      -gx_force,     --Fairlead force in global X (include applied forces) [N]\n");
  printf("      -gy_force,     --Fairlead force in global Y (include applied forces) [N]\n");
  printf("      -gz_force,     --Fairlead force in global Z (include applied forces) [N]\n");
  printf("      -h_fair,       --Horizontal force at fairlead (does NOT include applied forces) [N]\n");
  printf("      -v_fair,       --Vertical force at fairlead (does NOT include applied forces) [N]\n");
  printf("      -h_anch,       --Horizontal force at anchor (does NOT include applied forces) [N]\n");
  printf("      -v_anch,       --Vertical force at anchor (does NOT include applied forces) [N]\n");
  printf("      -tension_fair, --Line force-magnitude at fairlead (include applied loads) [N]\n");
  printf("      -tension_anch, --Line force-magnitude at anchor (include applied loads) [N]\n");
  printf("      -azimuth,      --Line lateral offset angle global X axis [deg]\n");
  printf("      -altitude,     --Line inclination angle relative to global XY plane at fiarlead [deg]\n");
  printf("      -lay_length,   --Length of line on seabed [m]\n");
  printf("      -line_tension, -- \n");
  printf("    Model features:\n");
  printf("      -omit_contact,       --Ignore cable/seabed contact\n");
  printf("      -seg_size <10>,      --Number of discrete elements in line\n");
  printf("      -damage_time <NULL>, --Line breakage occurs at specified time [s]\n");
  printf("      -diagnostic,         --Run line solver diagnostics until specified time [s]\n");
  printf("\n");
  printf("  Model option definitions\n");
  printf("    General model features:\n");
  printf("      -ref_position <0.0> <0.0> <0.0>\n");
  printf("      -repeat <NULL> ... <NULL>\n");
  printf("    MSQS solver options:\n");
  printf("      -inner_ftol <float>,\n");
  printf("      -inner_gtol <float>,\n");
  printf("      -inner_xtol <float>,\n");
  printf("      -inner_max_its <int>,\n");
  printf("      -outer_tol <float>,\n");
  printf("      -outer_max_its <int>,\n");
  printf("      -outer_epsilon <float>,\n");
  printf("      -outer_bd,\n");
  printf("      -outer_cd,\n");
  printf("      -outer_fd,\n");
  printf("      -pg_cooked <1000.0> <1.0>,\n");
  printf("      -integration_dt <0.01>,\n");
  printf("    LM model feature (not suported yet):\n");
  printf("      -kb_default      --Seabed stiffness parameter\n");
  printf("      -cb_default      --Seabed damping parameter\n");
  printf("      -wave_kinematics --Enables wave kinematics to drag interaction from surface waves\n");
  printf( "\nMAP++ Copyright (C) 2014 and GNU GPL by Marco Masciola and others\n" );
  printf( "SimCList Copyright (C) 2010 by Mij <http://mij.oltrelinux.com/devel/simclist/>\n" );
  printf( "MinPack Copyright (C) 1999 by the University of Chicago\n" );
  printf( "Modifications to MinPack by Frederic Devernay <http://devernay.free.fr/hacks/cminpack/>\n" );
  printf( "\nMAP++ is free software; see the source for copying conditions.\n" );
  printf( "This software is distributed on an \"AS IS\" BASIS, WITHOUT WARRANTIES\n" );
  printf( "OR CONDITIONS OF ANY KIND, either express or implied. See\n" );
  printf( "<https://gnu.org/licenses/gpl.html> for more details.\n" );
  printf("    \n");
  return MAP_SAFE;
};



