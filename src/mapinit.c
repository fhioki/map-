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
#include "cminpack/minpack.h"
#include "pyprotos.h"



/**
 * This sets the pointers to NULL for the vessel object and gives it default properties. Only 
 * to be used in the python glue code. 
 *
 * @todo: need to associate the node with inputs
 * @acceses: set_vartype_float( )
 * @calledby: mapcall_msqs_init( )
 */
MAP_ERROR_CODE set_vessel(Vessel* floater, const MAP_InputType_t* uType, char* map_msg, MAP_ERROR_CODE* ierr)
{
  MAP_ERROR_CODE success = MAP_SAFE;
  int i = 0;
  int N = uType->x_Len;

  do {
    /* vessel displacement */
    success = set_vartype_float("[m]", "Vessel_X", -999.9, &floater->displacement.x, 0.0); CHECKERRQ(MAP_FATAL_68);
    success = set_vartype_float("[m]", "Vessel_Y", -999.9, &floater->displacement.y, 0.0); CHECKERRQ(MAP_FATAL_68);
    success = set_vartype_float("[m]", "Vessel_Z", -999.9, &floater->displacement.z, 0.0); CHECKERRQ(MAP_FATAL_68);
  
    /* vessel reference origin. When ==[0.0, 0.0, 0.0], then the reference origin is aligned with the SWL */
    success = set_vartype_float("[m]", "Vessel_Xref", -999.9, &floater->refOrigin.x, 0.0); CHECKERRQ(MAP_FATAL_68);
    success = set_vartype_float("[m]", "Vessel_Yref", -999.9, &floater->refOrigin.y, 0.0); CHECKERRQ(MAP_FATAL_68);
    success = set_vartype_float("[m]", "Vessel_Zref", -999.9, &floater->refOrigin.z, 0.0); CHECKERRQ(MAP_FATAL_68);
  
    /* sum force of all fairleads connecte to the vessel */
    success = set_vartype_float("[N]", "Vessel_fx", -999.9, &floater->lineSumForce.fx, 0.0); CHECKERRQ(MAP_FATAL_68);
    success = set_vartype_float("[N]", "Vessel_fy", -999.9, &floater->lineSumForce.fy, 0.0); CHECKERRQ(MAP_FATAL_68);
    success = set_vartype_float("[N]", "Vessel_fz", -999.9, &floater->lineSumForce.fz, 0.0); CHECKERRQ(MAP_FATAL_68);
  
    /* orientation of the vessel. This is used as input from the user */
    success = set_vartype_float("[deg]", "Vessel_phi", -999.9, &floater->orientation.phi, 0.0); CHECKERRQ(MAP_FATAL_68);
    success = set_vartype_float("[deg]", "Vessel_the", -999.9, &floater->orientation.the, 0.0); CHECKERRQ(MAP_FATAL_68);
    success = set_vartype_float("[deg]", "Vessel_psi", -999.9, &floater->orientation.psi, 0.0); CHECKERRQ(MAP_FATAL_68);
  } while(0);

  floater->xi = (double*)malloc(N*sizeof(double));  
  floater->yi = (double*)malloc(N*sizeof(double));  
  floater->zi = (double*)malloc(N*sizeof(double));  

  if (floater->xi==NULL || floater->yi==NULL || floater->zi==NULL) {
    return MAP_FATAL;
  };
  
  for (i=0 ; i<N ; i++) {
    floater->xi[i] = uType->x[i];
    floater->yi[i] = uType->y[i];
    floater->zi[i] = uType->z[i];
  };
    
  return MAP_SAFE;
};


/**
 *  OBTAINED FROM THE CMINPACK SOURCE. The following defintions of the inputs for lmder(...) are provided in the cminpack docmentation.  
 *  
 *  info = __cminpack_func__(lmder)( inner_function_evals, elementIter, m, n, x, fvec, fjac, ldfjac, ftol, xtol, gtol, 
 *                                   maxfev, diag, mode, factor, nprint, &nfev, &njev, ipvt, qtf, wa1, wa2, wa3, wa4);
 *     
 *      - ftol          : a nonnegative input variable. Termination occurs when both the actual and predicted relative reductions in the sum of squares are at most ftol. Therefore, ftol measures the relative error desired in the sum of squares.
 *      - xtol          : a nonnegative input variable. Termination occurs when the relative error between two consecutive iterates is at most xtol. Therefore, xtol measures the relative error desired in the approximate solution.
 *      - gtol          : a nonnegative input variable. Termination occurs when the cosine of the angle between fvec and any column of the Jacobian is at most gtol in absolute value. Therefore, gtol measures the orthogonality desired between the function vector and the columns of the Jacobian.
 *      - maxfev        : a positive integer input variable. Termination occurs when the number of calls to fcn is at least maxfev by the end of an iteration.
 *      - diag          : an array of length n. If mode = 1 (see below), diag is internally set. If mode = 2, diag must contain positive entries that serve as multiplicative scale factors for the variables.
 *      - mode          : an integer input variable. If mode = 1, the variables will be scaled internally. If mode = 2, the scaling is specified by the input diag. Other values of mode are equivalent to mode = 1.
 *      - factor        : a positive input variable used in determining the initial step bound. This bound is set to the product of factor and the euclidean norm of diag*x if the latter is nonzero, or else to factor itself. In most cases factor should lie in the interval (.1,100.). 100. is a generally recommended value.
 *      - nprint        : an integer input variable that enables controlled printing of iterates if it is positive. In this case, fcn is called with iflag = 0 at the beginning of the first iteration and every nprint iterations thereafter and immediately prior to return, with x and fvec available for printing. If nprint is not positive, no special calls of fcn with iflag = 0 are made.
 *      - info          : an integer output variable. If the user has terminated execution, info is set to the (negative) value of iflag. See description of fcn. Otherwise, info is set as follows.
 *      - info=0        : improper input parameters.
 *      - info=1        : both actual and predicted relative reductions in the sum of squares are at most ftol.
 *      - info=2        : relative error between two consecutive iterates is at most xtol.
 *      - info=3        : conditions for info = 1 and info = 2 both hold.
 *      - info=4        : the cosine of the angle between fvec and any column of the Jacobian is at most gtol in absolute value.
 *      - info=5        : number of calls to fcn has reached or exceeded maxfev.
 *      - info=6        : ftol is too small. No further reduction in the sum of squares is possible.
 *      - info=7        : xtol is too small. No further improvement in the approximate solution x is possible.
 *      - info=8        : gtol is too small. fvec is orthogonal to the columns of the Jacobian to machine precision.
 *      - nfev          : an integer output variable set to the number of calls to fcn with iflag = 1.
 *      - njev          : an integer output variable set to the number of calls to fcn with iflag = 2.
 *      - ipvt          : an integer output array of length n. ipvt defines a permutation matrix p such that jac*p = q*r, where jac is the final calculated Jacobian, q is orthogonal (not stored), and r is upper triangular with diagonal elements of nonincreasing magnitude. Column j of p is column ipvt(j) of the identity matrix.
 *      - qtf           : an output array of length n which contains the first n elements of the vector (q transpose)*fvec.
 *      - wa1, wa2, wa3 : are work arrays of length n.
 *      - wa4           : a work array of length m.  
 */
MAP_ERROR_CODE first_solve(ModelData* data, MAP_InputType_t* uType, MAP_ConstraintStateType_t* zType,  MAP_OtherStateType_t* otherType, MAP_OutputType_t* yType, char* map_msg, MAP_ERROR_CODE* ierr)
{
  MAP_ERROR_CODE success = MAP_SAFE;

  if (data->MAP_SOLVE_TYPE==MONOLITHIC) {
    success = line_solve_sequence(data, 0.0, map_msg, ierr);
  } else {
    success = node_solve_sequence(data, uType, zType, otherType, map_msg, ierr); // @todo CHECKERRQ()
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


/**
 * @see mapcall_msqs_init( )
 * @see mapcall_msqs_end( )
 * @see free_outer_solve_data( ) for where this data is deallocated
 */
MAP_ERROR_CODE allocate_outer_solve_data(OuterSolveAttributes* ns, const int size, char* map_msg, MAP_ERROR_CODE* ierr)
{
  int i = 0;
  const int SIZE = 3*size;

  ns->jac = (double**)malloc(SIZE*sizeof(double*));
  ns->l = (double**)malloc(SIZE*sizeof(double*));  
  ns->u = (double**)malloc(SIZE*sizeof(double*));  
  ns->x = (double*)malloc(SIZE*sizeof(double));
  ns->b = (double*)malloc(SIZE*sizeof(double));
  ns->y = (double*)malloc(SIZE*sizeof(double*));  
  
  if (ns->jac==NULL) {
    *ierr = map_set_universal_error("", map_msg, ierr, MAP_FATAL_8);    
    return MAP_FATAL;
  };

  if (ns->x==NULL) {
    *ierr = map_set_universal_error("", map_msg, ierr, MAP_FATAL_8);    
    return MAP_FATAL;
  };

  if (ns->b==NULL) {
    *ierr = map_set_universal_error("", map_msg, ierr, MAP_FATAL_8);    
    return MAP_FATAL;
  };

  if (ns->l==NULL) {
    *ierr = map_set_universal_error("", map_msg, ierr, MAP_FATAL_8);    
    return MAP_FATAL;
  };

  if (ns->u==NULL) {
    *ierr = map_set_universal_error("", map_msg, ierr, MAP_FATAL_8);    
    return MAP_FATAL;
  };

  if (ns->y==NULL) {
    *ierr = map_set_universal_error("", map_msg, ierr, MAP_FATAL_8);    
    return MAP_FATAL;
  };

  for(i=0 ; i<SIZE ; i++) {
    ns->jac[i] = (double*)malloc(SIZE*sizeof(double));    
    ns->l[i] = (double*)malloc(SIZE*sizeof(double));    
    ns->u[i] = (double*)malloc(SIZE*sizeof(double));    
  };

  return MAP_SAFE;
};


/**
 *
 */
MAP_ERROR_CODE free_outer_solve_data(OuterSolveAttributes* ns, const int size, char* map_msg, MAP_ERROR_CODE* ierr)
{
  const int SIZE = 3*size;
  int i = 0;

  if (ns->jac) { /* is it allocated? */
    for(i=0 ; i<SIZE ; i++) {
      MAPFREE(ns->jac[i]);
    };
  };
 
  if (ns->l) { /* is it allocated? */
    for(i=0 ; i<SIZE ; i++) {
     MAPFREE(ns->l[i]);
    };
  };

  if (ns->u) { /* is it allocated? */
    for(i=0 ; i<SIZE ; i++) {
      MAPFREE(ns->u[i]);
   };  
  };

  MAPFREE(ns->jac);
  MAPFREE(ns->l);
  MAPFREE(ns->u);
  MAPFREE(ns->b);
  MAPFREE(ns->x);  
  MAPFREE(ns->y);
  return MAP_SAFE;
};


/**
 * MAP_InputType_t* uType,
 * MAP_ConstraintStateType_t* zType,
 * MAP_OtherStateType_t* otherType,
 * MAP_OutputType_t* yType,
 */
MAP_ERROR_CODE node_solve_sequence(ModelData* data, MAP_InputType_t* uType, MAP_ConstraintStateType_t* zType, MAP_OtherStateType_t* otherType, char* map_msg, MAP_ERROR_CODE* ierr)
{

  OuterSolveAttributes* ns = &data->outer_loop;
  MAP_ERROR_CODE success = MAP_SAFE;
  Element* elementIter = NULL;
  const int THREE = 3;
  const int N = zType->z_Len;
  const int m = THREE*(otherType->Fz_connect_Len); // rows
  const int n = THREE*(zType->z_Len);              // columns
  double error = 0.0;
  int SIZE = THREE*N;
  int col = 0;
  int row = 0;
  int i = 0;
  int j = 0;
  int lineCounter = 0;

  ns->iterationCount = 1;
  do {
    error = 0.0;
    success = line_solve_sequence(data, 0.0, map_msg, ierr); CHECKERRQ(MAP_FATAL_79);
    switch (ns->fd) {
    case BACKWARD_DIFFERENCE :
      success = backward_difference_jacobian(otherType, zType, data, map_msg, ierr); CHECKERRQ(MAP_FATAL_75);
      break;
    case CENTRAL_DIFFERENCE :
      success = central_difference_jacobian(otherType, zType, data, map_msg, ierr); CHECKERRQ(MAP_FATAL_76);
      break;
    case FORWARD_DIFFERENCE :
      success = forward_difference_jacobian(otherType, zType, data, map_msg, ierr); CHECKERRQ(MAP_FATAL_77);
      break;
    }

    success = line_solve_sequence(data, 0.0, map_msg, ierr);
    success = lu(ns, SIZE, map_msg, ierr); CHECKERRQ(MAP_FATAL_74);
    success = lu_back_substitution(ns, SIZE, map_msg, ierr); CHECKERRQ(MAP_FATAL_74);
  
    /* 
       Note that: ns->x = J^(-1) * F
       [x,y,z]_i+1 =  [x,y,z]_i - J^(-1) * F        
    */   
    for (i=0 ; i<N ; i++) { 
      zType->x[i] -= ns->x[THREE*i];
      zType->y[i] -= ns->x[THREE*i+1];
      zType->z[i] -= ns->x[THREE*i+2];
      error += (pow(otherType->Fx_connect[i],2)+ pow(otherType->Fy_connect[i],2) + pow(otherType->Fz_connect[i],2));
    };
    ns->iterationCount++;
    if (ns->iterationCount>ns->maxIts) {
      *ierr = map_set_universal_error("", map_msg, ierr, MAP_FATAL_80);      
      break;
    };
    // printf("Error: %f, tol=%f\n",sqrt(error),ns->tol);
    /* 
       @todo: end when iterations is exceeded. need some way to indicate that simulation did not suuficiently 
       meet termination criteria
    */
  } while (sqrt(error)>ns->tol);

  if (success==MAP_SAFE) {
    return MAP_SAFE;
  } else if (success==MAP_ERROR) {
    return MAP_ERROR;
  } else {
    return MAP_FATAL;
  };
};


/**
 *
 */
MAP_ERROR_CODE solve_line(ModelData* data, double time, char* map_msg, MAP_ERROR_CODE* ierr)
{
  MAP_ERROR_CODE success = MAP_SAFE;
  Element* elementIter = NULL;
  int lineCounter = 1;
  int cx = 0;
  char buffer[64] = "";
      
  list_iterator_start(&data->element);            /* starting an iteration "session" */
  while (list_iterator_hasnext(&data->element)) { /* tell whether more values available */ 
    elementIter = (Element*)list_iterator_next(&data->element);
 
    /* 
       check for fatal errors 
       @todo: this should handle the exception of a perfectly vertical cable. Added to the priority list
    */    
    if (elementIter->l.value<=0.0) {
      cx = map_snprintf(buffer, 64, "Line segment %d, l = %f [m].", lineCounter, elementIter->l.value); assert(cx>=0);
      *ierr = map_set_universal_error(buffer, map_msg, ierr, MAP_FATAL_54);      
      list_iterator_stop(&data->element); /* ending the iteration session */   
      break; 
    } else if (elementIter->h.value<=-MACHINE_EPSILON) {
      cx = map_snprintf(buffer, 64, "Line segment %d, h = %f [m].", lineCounter, elementIter->h.value); assert(cx>=0);
      *ierr = map_set_universal_error(buffer, map_msg, ierr, MAP_FATAL_55);      
      list_iterator_stop(&data->element); /* ending the iteration session */   
      break; 
    } else if (elementIter->lineProperty->omega>0.0) {
      success = check_maximum_line_length(elementIter, elementIter->options.omitContact, map_msg, ierr);
      if (success) {        
        cx = map_snprintf(buffer, 64, "Line segment %d.", lineCounter); assert(cx>=0);
        *ierr = map_set_universal_error(buffer, map_msg, ierr, MAP_FATAL_59);
        break;
      };
    };    
    success = call_minpack_lmder(elementIter, &data->inner_loop, &data->modelOptions, lineCounter, time, map_msg, ierr); CHECKERRQ(MAP_FATAL_79);
    lineCounter++;
  };
  list_iterator_stop(&data->element); /* ending the iteration "session" */    

  if (*ierr==MAP_SAFE) {
    return MAP_SAFE;
  } else if (*ierr==MAP_ERROR) {
    return MAP_ERROR;
  } else {
    return MAP_FATAL;
  };
};


/**
 *
 */
MAP_ERROR_CODE set_model_options_list(ModelData* data, InitializationData* initObj, char* map_msg, MAP_ERROR_CODE* ierr)
{
  int i = 0;
  int index = 0;
  int cx = 0;
  int success = 0;
  int sizeOfString = 0;
  char* word = NULL;
  char* line = NULL;
  MapReal* moreAngles = NULL;
  char buffer[64] = "";  

  // initialize_options(data);
  for (i=0 ; i<=initObj->solverOptionsSize-1 ; i++) { 
    sizeOfString = strlen(initObj->solverOptionsString[i]);
    line = (char*)realloc(line, (sizeOfString+1)*sizeof(char));    
    strcpy(line, initObj->solverOptionsString[i]);    
    word = strtok( line, " ,\n\t\r\0");
    while (word!=NULL) {      
      if (!strcicmp("HELP",word)) {
        success = print_help_to_screen(map_msg,ierr);
      } else if (!strcicmp("INNER_FTOL",word)) {
        word = strtok(NULL, " ,\n\t\r\0");
        success = is_numeric(word); CHECKERRK(MAP_ERROR_2);
        if (!success) {
          data->inner_loop.f_tol = (MapReal)atof(word);        
        };
      } else if (!strcicmp("INNER_GTOL",word)) {
        word = strtok(NULL, " ,\n\t\r\0");
        success = is_numeric(word); CHECKERRK(MAP_ERROR_9);
        if (!success) {
          data->inner_loop.g_tol = (MapReal)atof(word);        
        };
      } else if (!strcicmp("INNER_XTOL",word)) {
        word = strtok( NULL, " ,\n\t\r\0" );
        success = is_numeric(word); CHECKERRK(MAP_ERROR_10);
        if (!success) {
          data->inner_loop.x_tol = (MapReal)atof(word);        
        };
      } else if (!strcicmp("OUTER_TOL",word)) {
        word = strtok(NULL, " ,\n\t\r\0");
        success = is_numeric(word); CHECKERRK(MAP_ERROR_3);
        if (!success) {
          data->outer_loop.tol = (MapReal)atof(word);        
        };
      } else if (!strcicmp("OUTER_MAX_ITS",word)) {
        word = strtok( NULL, " ,\n\t\r\0" );
        success = is_numeric(word); CHECKERRK(MAP_ERROR_3);
        if (!success) {
          data->outer_loop.maxIts = (MapReal)atof(word);        
        };
      } else if (!strcicmp("OUTER_EPSILON",word)) {
        word = strtok( NULL, " ,\n\t\r\0" );
        success = is_numeric(word); CHECKERRK(MAP_ERROR_3);
        if (!success) {
          data->outer_loop.epsilon = (MapReal)atof(word);        
        };
      } else if (!strcicmp("INNER_MAX_ITS",word)) {
        word = strtok(NULL, " ,\n\t\r\0");
        success = is_numeric(word); CHECKERRK(MAP_ERROR_4);
        if (!success) {
          data->inner_loop.max_its = (MapReal)atof(word);        
        };
      } else if (!strcicmp("INTEGRATION_DT",word)) {
        word = strtok(NULL, " ,\n\t\r\0");
        success = is_numeric(word); CHECKERRK(MAP_ERROR_15); 
        if (!success) {
          data->modelOptions.integrationDt = (MapReal)atof(word);        
        };
      } else if (!strcicmp("KB_DEFAULT",word)) {
        word = strtok(NULL, " ,\n\t\r\0");
        success = is_numeric(word); CHECKERRK(MAP_ERROR_16); 
        if (!success) {
          data->modelOptions.integrationDt = (MapReal)atof(word);        
        };
      } else if (!strcicmp("CB_DEFAULT",word)) {
        word = strtok(NULL, " ,\n\t\r\0");
        success = is_numeric(word); CHECKERRK(MAP_ERROR_17); 
        if (!success) {
          data->modelOptions.integrationDt = (MapReal)atof(word);        
        };
      } else if (!strcicmp("WAVE_KINEMATICS",word)) {
        *ierr = map_set_universal_error("", map_msg, ierr, MAP_WARNING_10);
        data->modelOptions.waveKinematics = false;
      } else if (!strcicmp("PG_COOKED",word)) {
        word = strtok(NULL, " ,\t\0");        
        if (word!=NULL) {          
          success = is_numeric(word); CHECKERRK(MAP_WARNING_8);
          if (!success) {
            data->outer_loop.ds = (MapReal)atof(word);        
          } else {
            data->outer_loop.ds = 1.0;
          };
        };
        word = strtok(NULL, " ,\n\t\r\0");        
        if (word!=NULL) {          
          success = is_numeric(word); CHECKERRK(MAP_WARNING_8);
          if (!success) {
            data->outer_loop.d = (MapReal)atof(word);        
          } else {
            data->outer_loop.d = 0.0;
          };
        };
        data->outer_loop.pg = true;
      } else if (!strcicmp("OUTER_BD",word)) {
        data->outer_loop.fd = BACKWARD_DIFFERENCE;
      } else if (!strcicmp("OUTER_CD",word)) {
        data->outer_loop.fd = CENTRAL_DIFFERENCE;
      } else if (!strcicmp("OUTER_FD",word)) {
        data->outer_loop.fd = FORWARD_DIFFERENCE;
      } else if (!strcicmp("REF_POSITION",word)) {
        word = strtok( NULL, " ,\n\t\r\0" );
        success = is_numeric(word); CHECKERRQ(MAP_FATAL_36);
        data->vessel.refOrigin.x.value = (MapReal)atof(word);        
        word = strtok( NULL, " ,\n\t\r\0" );
        success = is_numeric(word); CHECKERRQ(MAP_FATAL_36);
        data->vessel.refOrigin.y.value = (MapReal)atof(word);        
        word = strtok( NULL, " ,\n\t\r\0" );
        success = is_numeric(word); CHECKERRQ(MAP_FATAL_36);
        data->vessel.refOrigin.z.value = (MapReal)atof(word);
      } else if (!strcicmp("REPEAT",word)) {
        word = strtok(NULL, " ,\t\n");
        while (word!=NULL) {      
          index = data->modelOptions.sizeOfRepeatAngles;
          moreAngles = (MapReal*)realloc(data->modelOptions.repeatAngles, (index+1)*sizeof(MapReal));
          if (moreAngles!=NULL) {
            data->modelOptions.repeatAngles = moreAngles;
            success = is_numeric(word); CHECKERRQ(MAP_FATAL_34);
            data->modelOptions.repeatAngles[index] = (MapReal)atof(word);
            data->modelOptions.sizeOfRepeatAngles++;
          } else {
            MAPFREE(moreAngles);
            success = 1; CHECKERRQ(MAP_FATAL_35);
          };
          word = strtok(NULL, " ,\t\n");
        };
      } else {
        cx = map_snprintf(buffer, 64, "Ignoring value '%s'.", word); assert(cx>=0);
        *ierr = map_set_universal_error(buffer, map_msg, ierr, MAP_WARNING_4);
      };
      word = strtok(NULL, " ,\n\t\r\0");
    };
  };
  MAPFREE(line);
  return MAP_SAFE;
};



/**
 *
 */
int set_cable_library_list(ModelData* data, InitializationData* init, char* map_msg, MAP_ERROR_CODE* ierr)
{
  int i = 0;
  int cx = 0;
  int cnt = 0;
  int success = 0;
  int len = 0;
  int sizeOfString = 0;
  char* line = NULL;
  char* word = NULL;
  char buffer[64] = "";  
  CableLibrary in_cableLibrary;

  for(i=0 ; i<=init->librarySize-1 ; i++) { 
    sizeOfString = strlen(init->libraryInputString[i]);
    line = (char*)realloc(line, (sizeOfString+1)*sizeof(char));    
    strcpy(line, init->libraryInputString[i]);        
    word=strtok(line, " ,\n\t\r\0");        
    while (word!=NULL) {            
      if (cnt==0) {
        len = strlen(word);
        in_cableLibrary.label = malloc(sizeof(char*)*(len+1));
        strcpy(in_cableLibrary.label, word);
      } else if (cnt==1) {
        success=is_numeric(word);CHECKERRQ(MAP_FATAL_12);
        in_cableLibrary.diam=(MapReal)atof(word);
      } else if (cnt==2) {
        success=is_numeric(word);CHECKERRQ(MAP_FATAL_13);
        in_cableLibrary.massDensityInAir=(MapReal)atof(word);
      } else if (cnt==3) {
        success=is_numeric(word);CHECKERRQ(MAP_FATAL_14);
        in_cableLibrary.ea=(MapReal)atof(word);
      } else if (cnt==4) {
        success=is_numeric(word); CHECKERRQ(MAP_FATAL_15);
        in_cableLibrary.cb=(MapReal)atof(word);
      } else if (cnt==5) {
        success=is_numeric(word); CHECKERRQ(MAP_FATAL_81);
        in_cableLibrary.cIntDamp=(MapReal)atof(word);
      } else if (cnt==6) {
        success=is_numeric(word); CHECKERRQ(MAP_FATAL_82);
        in_cableLibrary.cAdded=(MapReal)atof(word);
      } else if (cnt==7) {
        success=is_numeric(word); CHECKERRQ(MAP_FATAL_83);
        in_cableLibrary.cDragNormal=(MapReal)atof(word);
      } else if (cnt==8) {
        success=is_numeric(word); CHECKERRQ(MAP_FATAL_84);
        in_cableLibrary.cDragTangent=(MapReal)atof(word);
      } else { /* too many options are listed in the MAP input file. Let's warn users of uncommented information */
        if (strcmp (" ",word)) {
          cx = map_snprintf(buffer, 64, "Value '%s'.", word);assert(cx>=0);
          *ierr=map_set_universal_error(buffer, map_msg, ierr, MAP_WARNING_1);
        };
      };
      cnt++;
      word = strtok(NULL, " ,\n\t\r\0");
    };  
    list_append(&data->cableLibrary, &in_cableLibrary);
    cnt=0;
  };
  
  MAPFREE(line);

  if (*ierr==MAP_FATAL) {
    return MAP_FATAL;
  };
  return MAP_SAFE;  
};


/**
 * Initialized omega (weight per unit length) and cross-section area of a cable. The formula is 
 *
 *   A=\pi*\frac{radius^{2}}{4}
 *   \omega=g*(\mu-A*\rho_{seawater})
 *
 * @acceses: none
 * @calledby: mapcall_msqs_init( )
 */
MAP_ERROR_CODE initialize_cable_library_variables(ModelData* data, MAP_ParameterType_t* paramFortType, char* map_msg, MAP_ERROR_CODE* ierr)
{
  MapReal radius = 0.0;
  MapReal area = 0.0;
  MapReal muCable = 0.0;
  MapReal rhoFluid = 0.0; 
  const MapReal g = paramFortType->g;
  const MapReal PI = 3.14159264;
  CableLibrary* iterCableLibrary=NULL;

  list_iterator_start(&data->cableLibrary); /* starting an iteration "session" */
  while ( list_iterator_hasnext(&data->cableLibrary)) { /* tell whether more values available */ 
    iterCableLibrary = (CableLibrary*)list_iterator_next(&data->cableLibrary);
    radius = iterCableLibrary->diam/2;
    area = PI*pow(radius,2);
    muCable = iterCableLibrary->massDensityInAir;
    rhoFluid = paramFortType->rhoSea;
    iterCableLibrary->omega = g*(muCable-area*rhoFluid);

    iterCableLibrary->a = area;
    if (fabs(iterCableLibrary->omega)<=1) {
      *ierr = map_set_universal_error("", map_msg, ierr, MAP_WARNING_5);
    };
  };
  list_iterator_stop(&data->cableLibrary); /* ending the iteration "session" */    
  
  if (fabs(iterCableLibrary->omega)<=1e-3) {
    return MAP_FATAL;
  }
  return MAP_SAFE;
};


/**
 *  char** list; 
 *  
 *  list = malloc(sizeof(char*)*number_of_row);
 *  for(i=0;i<number_of_row; i++) list[i] = malloc(sizeof(char)*number_of_col);  
 * 
 *  ....
 *  
 *  for(i=0;i<number_of_row; i++) MAPFREE(list[i] );
 *  MAPFREE(list);  
 * 
 *
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 *
 *
 *  char **array = malloc(sizeof(*array) * array_size)
 *  ...
 *  while (storing strings in array)
 *  {
 *      if (amount of strings is going to exceed array_size)
 *      {
 *          new_array_size = array_size * growth_factor
 *          char **temp_array = realloc(array, sizeof(*array) * new_array_size)
 *          if (temp_array)
 *          {
 *              array = temp_array
 *              array_size = new_array_size
 *          }
 *          else
 *              error handling
 *      // rest of processing
 *  }
 */
MAP_ERROR_CODE repeat_nodes(ModelData* dataObj, InitializationData* init, char* map_msg, MAP_ERROR_CODE* ierr)
{
  int i = 0;
  int j = 0;
  int cnt = 0;
  int sizeOfNodeBase = 0;
  int sizeOfString = 0;
  int numRepeat = 0;
  char buffer[256] = "";
  char* word = NULL;
  char lines[256] = "";
  char** tempArray = NULL;
  char tempString[64] = "";
  MapReal angle = -999.9;
  Vector force;     
  Vector forceOriginal;
  Vector position;
  Vector positionOriginal;
  Node newNode;

  force.x = 0.0;
  force.y = 0.0;
  force.z = 0.0;
  position.x = 0.0;
  position.y = 0.0;
  position.z = 0.0;
  forceOriginal.x = 0.0;
  forceOriginal.y = 0.0;
  forceOriginal.z = 0.0;
  positionOriginal.x = 0.0;
  positionOriginal.y = 0.0;
  positionOriginal.z = 0.0;

  sizeOfNodeBase = dataObj->sizeOfNodes;
  
  numRepeat = dataObj->modelOptions.sizeOfRepeatAngles;
  init->sizeOfFullNodeString = (sizeOfNodeBase)*(numRepeat+1); 
  tempArray = malloc(sizeof(char*)*(init->sizeOfFullNodeString));

  for(j=0 ; j<sizeOfNodeBase ; j++) { 
    sizeOfString = strlen(init->nodeInputString[j]);
    tempArray[j] = malloc(sizeof(char)*(sizeOfString+1));
    strcpy(tempArray[j], init->nodeInputString[j]);
  };
  
  for(i=0 ; i<numRepeat ; i++) { 
    for(j=0 ; j<sizeOfNodeBase ; j++) { 
      set_node_null(&newNode);  
      sizeOfString = strlen(init->nodeInputString[j]);
      
      strcpy(lines, init->nodeInputString[j]);
      word = strtok(lines, " ,\n\t\r\0");    
      
      cnt = 0;
      strcpy(buffer, "");
      angle = dataObj->modelOptions.repeatAngles[i]*(DEG2RAD);  // @todo : this needs to change index[0] for angles
      while (word!=NULL) {      
        if (cnt==0) {
          map_snprintf(tempString, 64, "%d  ", ((i+1)*sizeOfNodeBase+j+1));
          map_strcat(buffer, 256, tempString);
        } else if (cnt==1) { /* compare strings to determine what type of node it is connected to */
          if (!strcicmp("FIX",word)) { 
            map_strcat( buffer, 256, "fix  " );
          } else if (!strcicmp("CONNECT",word)) {
            map_strcat(buffer, 256, "connect  ");
          } else {
            map_strcat(buffer, 256, "vessel  ");
          };
        } else if (cnt==2) {           
          if (word[0]=='#') { 
            positionOriginal.x = (MapReal)atof(remove_first_character(word));
          } else {
            positionOriginal.x = (MapReal)atof(word);
          };
        } else if (cnt==3) { 
          if (word[0]=='#') { 
            positionOriginal.y = (MapReal)atof(remove_first_character(word));
          } else {
            positionOriginal.y = (MapReal)atof(word);
          };
        } else if (cnt==4) { 
          position.x =  positionOriginal.x*cos(angle) + positionOriginal.y*sin(angle);
          position.y = -positionOriginal.x*sin(angle) + positionOriginal.y*cos(angle);                    
          if (word[0]=='#') { 
            position.z = (MapReal)atof(remove_first_character(word));
            map_snprintf(tempString, 64, "#%1.4f  #%1.4f  #%1.4f  ", position.x, position.y, position.z);
          } else {
            if (!is_numeric(word)) {
              position.z = (MapReal)atof(word);
              map_snprintf(tempString, 64, "%1.4f  %1.4f  %1.4f  ", position.x, position.y, position.z);
            } else {
              map_snprintf(tempString, 64, "%1.4f  %1.4f  %s  ", position.x, position.y, word);
            };
          };
          map_strcat(buffer, 256, tempString);        
        } else if (cnt==5) { 
          map_snprintf(tempString, 64, "%s  ", word);
          map_strcat(buffer, 256, tempString);
        } else if (cnt==6 ) { 
          map_snprintf(tempString, 64, "%s  ", word);
          map_strcat(buffer, 256, tempString);
        } else if (cnt==7) { 
          if (word[0]=='#') { 
            forceOriginal.x = (MapReal)atof(remove_first_character(word));
          } else {
            forceOriginal.x = (MapReal)atof(word);
          };
        } else if (cnt==8) { 
          if (word[0]=='#') { 
            forceOriginal.y = (MapReal)atof(remove_first_character(word));
          } else {
            forceOriginal.y = (MapReal)atof(word);
          };
        } else if (cnt==9) { 
          force.x =  forceOriginal.x*cos(angle) + forceOriginal.y*sin(angle);
          force.y = -forceOriginal.x*sin(angle) + forceOriginal.y*cos(angle);
          if (word[0]=='#') { 
            force.z = (MapReal)atof(remove_first_character(word));
            map_snprintf(tempString, 64, "#%1.4f  #%1.4f  #%1.4f\n", force.x, force.y, force.z);
          } else {
            force.z = (MapReal)atof(word);
            map_snprintf(tempString, 64, "%1.4f  %1.4f  %1.4f\n", force.x, force.y, force.z);
          };
          map_strcat(buffer, 256, tempString);        
        };
        cnt++;
        word = strtok(NULL, " ,\n\t\r\0");
      };  
      sizeOfString = strlen(buffer);
      tempArray[(i+1)*sizeOfNodeBase+j] = malloc(sizeof(char)*(sizeOfString+1));
      strcpy(tempArray[(i+1)*sizeOfNodeBase+j], buffer);
    };  
  };

  dataObj->sizeOfNodes = init->sizeOfFullNodeString;
  init->expandedNodeInputString = malloc(sizeof(char*)*(init->sizeOfFullNodeString));
  for( i=0 ; i<init->sizeOfFullNodeString ; i++ ) {    
    init->expandedNodeInputString[i] = malloc(sizeof(char)*(strlen(tempArray[i])+1));  
    strcpy(init->expandedNodeInputString[i], tempArray[i]);
    MAPFREE(tempArray[i]);
  };
  MAPFREE(tempArray);

  return MAP_SAFE;
};


/**
 *
 */
MAP_ERROR_CODE repeat_elements(ModelData* dataObj, InitializationData* init, char* map_msg, MAP_ERROR_CODE* ierr)
{
  int i = 0;
  int j = 0;
  int cnt = 0;
  int success = 0;
  int sizeOfElementBase = 0;
  int sizeOfNodeBase = 0;
  int sizeOfString = 0;
  int numRepeat = 0;
  int anchNode = 0.0;
  int fairNode = 0.0;
  char buffer[256] = "";
  char* word = NULL;
  char lines[256] = "";
  char** temp_array  = NULL;
  char temp_string[64] = "";
  Element newElement;

  sizeOfElementBase = dataObj->sizeOfElements;
  sizeOfNodeBase = init->nodeSize;
  
  numRepeat = dataObj->modelOptions.sizeOfRepeatAngles;
  init->sizeOfFullElementString = (sizeOfElementBase)*(numRepeat+1); 
  temp_array = malloc( (init->sizeOfFullElementString)*sizeof(char*));

  for(j=0 ; j<sizeOfElementBase ; j++) { 
    sizeOfString = strlen(init->elementInputString[j]);
    temp_array[j] = malloc(sizeof(char)*(sizeOfString+1));
    strcpy(temp_array[j], init->elementInputString[j]);
  };

  for(i=0 ; i<numRepeat ; i++) { 
    for(j=0 ; j<sizeOfElementBase ; j++) { 
      set_element_null(&newElement);  
      sizeOfString = strlen(init->elementInputString[j]);      
      strcpy(lines, init->elementInputString[j]);
      word = strtok(lines, " ,\n\t\r\0" );    
      cnt = 0;
      buffer[0] = 0;
      temp_string[0] = 0;
      while (word!=NULL) {      
        if (cnt==0) {
          map_snprintf(temp_string, 64, "%d  ", ((i+1)*sizeOfElementBase+j+1));
          map_strcat(buffer, 256, temp_string);
        } else if (cnt==1 ) { 
          map_snprintf(temp_string, 64, "%s  ", word);
          map_strcat(buffer, 256, temp_string);
        } else if (cnt==2) { 
          map_snprintf(temp_string, 64, "%s  ", word);
          map_strcat(buffer, 256, temp_string);
        } else if (cnt==3) { 
          success = is_numeric(word); CHECKERRQ(MAP_FATAL_28);
          if ( !success ) {
            anchNode = (MapReal)atof(word);        
          };
          map_snprintf(temp_string, 64, "%d  ", (i+1)*sizeOfNodeBase+anchNode);
          map_strcat(buffer, 256, temp_string);
        } else if (cnt==4) { 
          success = is_numeric(word); CHECKERRQ(MAP_FATAL_29);
          if ( !success ) {
            fairNode = (MapReal)atof(word);  
          };
          map_snprintf(temp_string, 64, "%d  ", (i+1)*sizeOfNodeBase+fairNode);
          map_strcat(buffer, 256, temp_string);
        } else {
          map_snprintf(temp_string, 64, "%s  ", word);
          map_strcat(buffer, 256, temp_string);
        };
        cnt++;
        word = strtok( NULL, " ,\n\t\r\0" );
      };  
      map_strcat( buffer, 256, "\n" );
      sizeOfString = strlen( buffer );
      temp_array[(i+1)*sizeOfElementBase+j] = malloc(sizeof(char)*(sizeOfString+1));
      strcpy(temp_array[(i+1)*sizeOfElementBase+j], buffer);
    };  
  };

  dataObj->sizeOfElements = init->sizeOfFullElementString;
  init->expandedElementInputString = malloc(sizeof(char*)*(init->sizeOfFullElementString));
  for(i=0 ; i<init->sizeOfFullElementString ; i++) {    
    init->expandedElementInputString[i] = malloc(sizeof(char)*(strlen(temp_array[i])+1));  
    strcpy(init->expandedElementInputString[i], temp_array[i]);
    MAPFREE(temp_array[i]);
  };
  MAPFREE(temp_array);

  return MAP_SAFE;

};


/**
 * Sets init data to NULL or -9999
 *
 * @acceses: none
 * @calledby: mapcall_msqs_init( )
 */
MAP_ERROR_CODE initialize_init_input(InitializationData* init)
{
  if (init==NULL) {
    return MAP_FATAL;
  };    
  init->expandedNodeInputString = NULL;
  init->expandedElementInputString = NULL;
  init->sizeOfFullNodeString = -9999;
  init->sizeOfFullElementString = -9999;  
  return MAP_SAFE;
};


/**
 *
 */
MAP_ERROR_CODE allocate_types_for_nodes(MAP_InputType_t* uType, MAP_ConstraintStateType_t* zType, MAP_OtherStateType_t* otherType, MAP_OutputType_t* yType, ModelData *data, char** const nodeInputString, char *map_msg, MAP_ERROR_CODE *ierr )
{
  int i = 0;
  int cx = 0;
  int cnt = 0;
  int fixNum = 0;
  int vesselNum = 0;
  int connectNum = 0;
  int sizeOfString = 0;
  char* word = NULL;
  char* line = NULL;
  char buffer[64] = "";

  /* 
     First need to run loop to find the number of inputs, outputs, constraints.
     Basically we are just counting the number of VESSEL, FIX, and CONNECT nodes     
  */
  for(i=0 ; i<=data->sizeOfNodes-1 ; i++) {         
    sizeOfString = strlen(nodeInputString[i]);
    line = (char*)realloc(line, (sizeOfString+1)*sizeof(char));    
    strcpy(line, nodeInputString[i]);
    word = strtok(line, " ,\n\t\r\0");
    while(word!=NULL) {            
      if (cnt==0) {
        /* this is ignored because it's not needed for the model assembly */
      } else if (cnt==1) { /* compare strings to determine what type of node it is connected to */        
        if (!strcicmp("FIX",word)) { 
          fixNum++;
        } else if (!strcicmp("CONNECT",word)) {
          connectNum++;
        } else if (!strcicmp("VESSEL",word)) {
          vesselNum++;
        } else {
          /* the node is given an invalid name */
          cx = map_snprintf(buffer, 64, "Value '%s'.", word); assert(cx>=0);
          *ierr = map_set_universal_error(buffer, map_msg, ierr, MAP_FATAL_25);          
        };
      };
      cnt++;
      word = strtok(NULL, " ,\n\t\r\0");
    };
    cnt=0;
  };
  MAPFREE(line);

  otherType->x_Len = fixNum;
  otherType->y_Len = fixNum;
  otherType->z_Len = fixNum;
  otherType->x = (double*)malloc(otherType->x_Len*sizeof(double));
  otherType->y = (double*)malloc(otherType->y_Len*sizeof(double));
  otherType->z = (double*)malloc(otherType->z_Len*sizeof(double));

  /* 
     if the node is VESSEL, then the applied force is an output state. Otherwise, it 
     has to be an other state because it can't be associated with any other type. 
     This is what is done below.
  */
  otherType->Fx_connect_Len = connectNum;
  otherType->Fy_connect_Len = connectNum;
  otherType->Fz_connect_Len = connectNum;
  otherType->Fx_connect = (double*)malloc(otherType->Fx_connect_Len*sizeof(double));
  otherType->Fy_connect = (double*)malloc(otherType->Fy_connect_Len*sizeof(double));
  otherType->Fz_connect = (double*)malloc(otherType->Fz_connect_Len*sizeof(double));

  otherType->Fx_anchor_Len = fixNum;
  otherType->Fy_anchor_Len = fixNum;
  otherType->Fz_anchor_Len = fixNum;
  otherType->Fx_anchor = (double*)malloc(otherType->Fx_anchor_Len*sizeof(double));
  otherType->Fy_anchor = (double*)malloc(otherType->Fy_anchor_Len*sizeof(double));
  otherType->Fz_anchor = (double*)malloc(otherType->Fz_anchor_Len*sizeof(double));
  
  zType->x_Len = connectNum;          
  zType->y_Len = connectNum;          
  zType->z_Len = connectNum;          
  zType->x = (double*)malloc(zType->x_Len*sizeof(double));
  zType->y = (double*)malloc(zType->y_Len*sizeof(double));
  zType->z = (double*)malloc(zType->z_Len*sizeof(double));
  
  uType->x_Len = vesselNum;
  uType->y_Len = vesselNum;
  uType->z_Len = vesselNum;
  yType->Fx_Len = vesselNum;
  yType->Fy_Len = vesselNum;
  yType->Fz_Len = vesselNum;
  uType->x = (double*)malloc(uType->x_Len*sizeof(double));
  uType->y = (double*)malloc(uType->y_Len*sizeof(double));
  uType->z = (double*)malloc(uType->z_Len*sizeof(double));
  yType->Fx = (double*)malloc(yType->Fx_Len*sizeof(double));
  yType->Fy = (double*)malloc(yType->Fy_Len*sizeof(double));
  yType->Fz = (double*)malloc(yType->Fz_Len*sizeof(double));
  return MAP_SAFE;
};


/**
 *
 */
MAP_ERROR_CODE compare_length(int a, int b)
{
  if (a!=b) {
    return MAP_FATAL;
  }; 
  return MAP_SAFE;
};


/**
 * @see: map_free_types to view where fortran types are free'd
 */
int set_node_list(const MAP_ParameterType_t* pType,  MAP_InputType_t* uType, MAP_ConstraintStateType_t* zType, MAP_OtherStateType_t* otherType, MAP_OutputType_t* yType, ModelData* data, char** const nodeInputString, char* map_msg, MAP_ERROR_CODE* ierr)
{
  int i = 0;
  int cx = 0;
  int cnt = 0;
  int success = 0;
  int sizeOfString = 0;
  int fixNum = 0;
  int vesselNum = 0;
  int connectNum = 0;
  Node newNode;
  char buffer[64] = "";
  char* word = NULL;
  char* line = NULL;
  const MapReal depth = pType->depth;

  success = allocate_types_for_nodes(uType, zType, otherType, yType, data, nodeInputString, map_msg, ierr);
  
  for(i=0 ; i<=data->sizeOfNodes-1 ; i++) {         
    set_node_null(&newNode);
    sizeOfString = strlen(nodeInputString[i]);
    line = (char*)realloc(line, (sizeOfString+1)*sizeof(char));    
    strcpy(line, nodeInputString[i]);
    word = strtok(line, " ,\n\t\r\0");    
    while (word!=NULL) {      
      if (cnt==0) {  
        /* this is ignored because it's not needed for the model assembly */
      } else if (cnt==1) { /* compare strings to determine what type of node it is connected to */
        if (!strcicmp("FIX",word)) { 
          newNode.type = FIX;
          fixNum++;
                                         /* VarTypePtr            FAST derived  array index */
          success = associate_vartype_ptr(&newNode.positionPtr.x, otherType->x, fixNum);
          success = associate_vartype_ptr(&newNode.positionPtr.y, otherType->y, fixNum);
          success = associate_vartype_ptr(&newNode.positionPtr.z, otherType->z, fixNum);
          success = associate_vartype_ptr(&newNode.sumForcePtr.fx, otherType->Fx_anchor, fixNum);
          success = associate_vartype_ptr(&newNode.sumForcePtr.fy, otherType->Fy_anchor, fixNum);
          success = associate_vartype_ptr(&newNode.sumForcePtr.fz, otherType->Fz_anchor, fixNum);
        } else if (!strcicmp("CONNECT",word)) {
          newNode.type = CONNECT;
          connectNum++;
          success = associate_vartype_ptr(&newNode.positionPtr.x, zType->x, connectNum);
          success = associate_vartype_ptr(&newNode.positionPtr.y, zType->y, connectNum);
          success = associate_vartype_ptr(&newNode.positionPtr.z, zType->z, connectNum);
          success = associate_vartype_ptr(&newNode.sumForcePtr.fx, otherType->Fx_connect, connectNum);
          success = associate_vartype_ptr(&newNode.sumForcePtr.fy, otherType->Fy_connect, connectNum);
          success = associate_vartype_ptr(&newNode.sumForcePtr.fz, otherType->Fz_connect, connectNum);
        } else if (!strcicmp("VESSEL",word)) {
          newNode.type = VESSEL;
          vesselNum++;
          success = associate_vartype_ptr(&newNode.positionPtr.x, uType->x, vesselNum);
          success = associate_vartype_ptr(&newNode.positionPtr.y, uType->y, vesselNum);
          success = associate_vartype_ptr(&newNode.positionPtr.z, uType->z, vesselNum);
          success = associate_vartype_ptr(&newNode.sumForcePtr.fx, yType->Fx, vesselNum);
          success = associate_vartype_ptr(&newNode.sumForcePtr.fy, yType->Fy, vesselNum);
          success = associate_vartype_ptr(&newNode.sumForcePtr.fz, yType->Fz, vesselNum);
        } else {
          /* the node is given an invalid name */
          cx = map_snprintf(buffer, 64, "Value '%s'.", word); assert(cx>=0);
          *ierr = map_set_universal_error(buffer, map_msg, ierr, MAP_FATAL_25);          
        };
        cx = map_snprintf(buffer, 64, "s(fx)[%d]", i+1); assert(cx>=0);        
        success = initialize_node_sum_force_ptr("[N]", buffer, i, &newNode.sumForcePtr.fx); CHECKERRQ(MAP_FATAL_22); 
        cx = map_snprintf(buffer, 64, "s(fy)[%d]", i+1); assert(cx>=0);        
        success = initialize_node_sum_force_ptr("[N]", buffer, i, &newNode.sumForcePtr.fy); CHECKERRQ(MAP_FATAL_23); 
        cx = map_snprintf(buffer, 64, "s(fz)[%d]", i+1); assert(cx>=0);        
        success = initialize_node_sum_force_ptr("[N]", buffer, i, &newNode.sumForcePtr.fz); CHECKERRQ(MAP_FATAL_24); 
      } else if (cnt==2) { 
        /* set initial X node position values */
        cx = map_snprintf(buffer, 64, "X[%d]", i+1); assert(cx>=0);        
        success = set_vartype_ptr("[m]", buffer, i, &newNode.positionPtr.x, word); CHECKERRQ(MAP_FATAL_17);
      } else if (cnt==3) { 
        /* set initial Y node position values */
        cx = map_snprintf(buffer, 64, "Y[%d]", i+1); assert(cx>=0);        
        success = set_vartype_ptr("[m]", buffer, i, &newNode.positionPtr.y, word); CHECKERRQ(MAP_FATAL_18);
      } else if (cnt==4) { 
        /* set initial Z node position values */
        cx = map_snprintf(buffer, 64, "Z[%d]", i+1); assert(cx>=0);        
        if (!strcicmp("DEPTH", word)) {          
          if (newNode.type!=FIX) { /* can only use 'DEPTH' flag in input file for FIX (anchor) nodes */
            *ierr = map_set_universal_error(buffer, map_msg, ierr, MAP_FATAL_71);            
          } else {
            success = set_vartype_ptr_float("[m]", buffer, i, &newNode.positionPtr.z, depth); CHECKERRQ(MAP_FATAL_19);
          };
        } else { /* all other nodes not using the 'DEPTH' flag */
          success = set_vartype_ptr("[m]", buffer, i, &newNode.positionPtr.z, word); CHECKERRQ(MAP_FATAL_19);
        };        
      } else if (cnt==5) { 
        /* set the node mass */
        cx = map_snprintf(buffer, 64, "M[%d]", i+1); assert(cx>=0);        
        success = set_vartype("[kg]", buffer, i, &newNode.MApplied, word); CHECKERRQ(MAP_FATAL_20);
      } else if (cnt==6) { 
        /* set the node buoyancy */
        cx = map_snprintf(buffer, 64, "B[%d]", i+1); assert(cx>=0);        
        success = set_vartype("[m^3]", buffer, i, &newNode.BApplied, word); CHECKERRQ(MAP_FATAL_21);
      } else if (cnt==7) { 
        /* set applied X external force (or user guess) of the node */        
        cx = map_snprintf(buffer, 64, "FX[%d]", i+1); assert(cx>=0);                
        success = initialize_external_applied_force("[N]", buffer, i, &newNode.externalForce.fx, word); CHECKERRQ(MAP_FATAL_22);
      } else if (cnt==8) { 
        /* set applied Y external force (or user guess) of the node */
        cx = map_snprintf(buffer, 64, "FY[%d]", i+1); assert(cx>=0);        
        success = initialize_external_applied_force("[N]", buffer, i, &newNode.externalForce.fy, word); CHECKERRQ(MAP_FATAL_23);
      } else if (cnt==9) { 
        /* set applied Z external force (or user guess) of the node */
        cx = map_snprintf(buffer, 64, "FZ[%d]", i+1); assert(cx>=0);        
        success = initialize_external_applied_force("[N]", buffer, i, &newNode.externalForce.fz, word); CHECKERRQ(MAP_FATAL_24);
      } else { /* too many user inputs are givenin the input file. Let users know in case was not intended */
        if (strcmp(" ",word)) {
          cx = map_snprintf(buffer, 64, "Value '%s'.", word); assert(cx>=0);
          *ierr = map_set_universal_error(buffer, map_msg, ierr, MAP_WARNING_2);
        };
      };
      cnt++;
      word = strtok(NULL, " ,\n\t\r\0");
    };  
    list_append(&data->node, &newNode);
    cnt = 0;
  };
  
  /* 
     check to make sure the number of allocated array spaces for fortran derived types matches 
     what was actually set in the node initialization front end.
  */     
  do {
    success = compare_length(otherType->Fx_connect_Len, connectNum); CHECKERRQ(MAP_FATAL_49);
    success = compare_length(otherType->Fy_connect_Len, connectNum); CHECKERRQ(MAP_FATAL_49);
    success = compare_length(otherType->Fz_connect_Len, connectNum); CHECKERRQ(MAP_FATAL_49);

    success = compare_length(otherType->Fx_anchor_Len, fixNum); CHECKERRQ(MAP_FATAL_49); // @todo: change error code
    success = compare_length(otherType->Fy_anchor_Len, fixNum); CHECKERRQ(MAP_FATAL_49); // @todo: change error code
    success = compare_length(otherType->Fz_anchor_Len, fixNum); CHECKERRQ(MAP_FATAL_49); // @todo: change error code

    success = compare_length(otherType->x_Len, fixNum); CHECKERRQ(MAP_FATAL_49);
    success = compare_length(otherType->y_Len, fixNum); CHECKERRQ(MAP_FATAL_49);
    success = compare_length(otherType->z_Len, fixNum); CHECKERRQ(MAP_FATAL_49);
              
    success = compare_length(uType->x_Len, vesselNum); CHECKERRQ(MAP_FATAL_50);
    success = compare_length(uType->y_Len, vesselNum); CHECKERRQ(MAP_FATAL_50);
    success = compare_length(uType->z_Len, vesselNum); CHECKERRQ(MAP_FATAL_50);    
              
    success = compare_length(yType->Fx_Len, vesselNum); CHECKERRQ(MAP_FATAL_51);
    success = compare_length(yType->Fy_Len, vesselNum); CHECKERRQ(MAP_FATAL_51);
    success = compare_length(yType->Fz_Len, vesselNum); CHECKERRQ(MAP_FATAL_51);    
              
    success = compare_length(zType->x_Len, connectNum); CHECKERRQ(MAP_FATAL_52);
    success = compare_length(zType->y_Len, connectNum); CHECKERRQ(MAP_FATAL_52);
    success = compare_length(zType->z_Len, connectNum); CHECKERRQ(MAP_FATAL_52);    
  } while (0);  
  
  MAPFREE(line);

  if (*ierr==MAP_FATAL) {
    return 1;
  };
  return 0;  
};


/**
 *
 */
MAP_ERROR_CODE map_set_option_flags( char *parsed_word, Element *elem, char *map_msg, MAP_ERROR_CODE *ierr )
{
  int cx = 0 ;
  int success = 0;
  char buffer[64] = "";
  
  if (!strcicmp("PLOT", parsed_word) ) elem->options.plotFlag = true;
  else if (!strcicmp("GX_POS", parsed_word)) {
    elem->options.gxPosFlag = true;
  } else if (!strcicmp("GY_POS", parsed_word)) {
    elem->options.gyPosFlag = true;
  } else if (!strcicmp("GZ_POS", parsed_word)) {
    elem->options.gzPosFlag = true;
  } else if (!strcicmp("GX_A_POS", parsed_word)) {
    elem->options.gxAnchorPosFlag = true;
  } else if (!strcicmp("GY_A_POS", parsed_word)) {
    elem->options.gyAnchorPosFlag = true;
  } else if (!strcicmp("GZ_A_POS", parsed_word)) {
    elem->options.gzAnchorPosFlag = true;
  } else if (!strcicmp("GX_FORCE", parsed_word)) {
    elem->options.gxForceFlag= true;
  } else if (!strcicmp("GY_FORCE", parsed_word)) {
    elem->options.gyForceFlag = true;
  } else if (!strcicmp("GZ_FORCE", parsed_word)) {
    elem->options.gzForceFlag = true;
  } else if (!strcicmp("H_FAIR", parsed_word)) {
    elem->options.HFlag = true;
  } else if (!strcicmp("H_ANCH", parsed_word)) {
    elem->options.HAnchorFlag = true;
  } else if (!strcicmp("V_FAIR", parsed_word)) {
    elem->options.VFlag = true;
  } else if (!strcicmp("V_ANCH", parsed_word)) {
    elem->options.VAnchorFlag = true;
  } else if (!strcicmp("TENSION_FAIR", parsed_word)) {
    elem->options.fairleadTensionFlag = true;
  } else if (!strcicmp("TENSION_ANCH", parsed_word)) {
    elem->options.anchorTensionFlag = true;
  } else if (!strcicmp("X_EXCURSION", parsed_word)) {
    elem->options.horizontalExcursionFlag = true;
  } else if (!strcicmp("Z_EXCURSION", parsed_word)) {
    elem->options.verticalExcursionFlag = true;
  } else if (!strcicmp("AZIMUTH", parsed_word)) {
    elem->options.azimuthFlag = true;
  } else if (!strcicmp("ALTITUDE", parsed_word)) {
    elem->options.altitudeFlag = true;
  } else if (!strcicmp("ALTITUDE_A", parsed_word)) {
    elem->options.altitudeAnchorFlag = true;
  } else if (!strcicmp("LINE_TENSION", parsed_word)) {
    elem->options.lineTensionFlag = true;
  } else if (!strcicmp("OMIT_CONTACT", parsed_word)) {
    elem->options.omitContact = true;
  } else if (!strcicmp("SEG_SIZE", parsed_word)) {
    parsed_word = strtok(NULL, " ,\n\t\r\0");
    success = is_numeric(parsed_word);
    if (success==MAP_SAFE) {
      elem->segmentSize = (MapReal)atof(parsed_word);
    } else { /* should not cancel the simulation; simply ignore it */      
      *ierr=map_set_universal_error(buffer, map_msg, ierr, MAP_ERROR_18); 
    };          
  } else if (!strcicmp("LAY_LENGTH", parsed_word)) {
    elem->options.layLengthFlag = true;
  } else if (!strcicmp("DAMAGE_TIME", parsed_word)) {
    parsed_word = strtok(NULL, " ,\n\t\r\0");
    success = is_numeric(parsed_word);
    if (success==MAP_SAFE) {
      elem->options.damageTimeFlag = true;
      elem->damageTime = (MapReal)atof(parsed_word);
    } else { /* should not cancel the simulation; simply ignore it */     
  *ierr=map_set_universal_error(buffer, map_msg, ierr, MAP_ERROR_1);
    };          
  } else if (!strcicmp("DIAGNOSTIC", parsed_word)) {
    parsed_word = strtok(NULL, " ,\n\t\r\0");
    success = is_numeric(parsed_word);
    if (success==MAP_SAFE) {
      elem->options.diagnosticsFlag = true;
      elem->diagnosticType = (int)atoi(parsed_word);
    } else {
      /* should not cancel the simulation; simply ignore it */
      *ierr = map_set_universal_error(buffer, map_msg, ierr, MAP_ERROR_14);
      elem->options.diagnosticsFlag = true;
      elem->diagnosticType = 0;
    };          
  } else {
    /* should not cancel the simulation; simply ignore it */
    cx = map_snprintf(buffer, 64, "Option '%s'.", parsed_word); assert(cx>=0);
    *ierr = map_set_universal_error(buffer, map_msg, ierr, MAP_WARNING_3);
  };
  return MAP_SAFE;
};



/**
 *
 */
int set_element_list(MAP_ConstraintStateType_t* zType, ModelData* data, char** const elementInputString, char* map_msg, MAP_ERROR_CODE* ierr)
{
  int i = 0;
  int cx = 0;
  int cnt = 0;
  int success = 0;
  int sizeOfString = 0;
  char buffer[64] = "";
  char* line = NULL;
  char* word = NULL;
  Element newElement;
  Element* elementIter = NULL;
  int lineCounter = 0;

  zType->H_Len = data->sizeOfElements;          
  zType->V_Len = data->sizeOfElements;          
  zType->H = (double*)malloc(zType->H_Len*sizeof(double));
  zType->V = (double*)malloc(zType->V_Len*sizeof(double));

  if (zType->H==NULL) {
    *ierr = map_set_universal_error(buffer, map_msg, ierr, MAP_FATAL_53);
    return MAP_FATAL;
  };

  if (zType->V==NULL) {
    *ierr = map_set_universal_error(buffer, map_msg, ierr, MAP_FATAL_53);
    return MAP_FATAL;
  };

  for(i=0 ; i<=data->sizeOfElements-1 ; i++) {
    sizeOfString = strlen(elementInputString[i]);
    line = (char*)realloc(line, (sizeOfString+1)*sizeof(char));    
    strcpy(line, elementInputString[i]);    
    word = strtok(line, " ,\n\t\r\0");
    while (word!=NULL) {
      do {
        if (cnt==0) {
          /* use this first option as an opportunity to set the run-time flags to false */
          initialize_element(&newElement, i+1);
          success = associate_vartype_ptr(&newElement.H, zType->H, i+1);
          success = associate_vartype_ptr(&newElement.V, zType->V, i+1);          
          cx = map_snprintf(buffer, 64, "V[%d]", i+1);assert(cx>=0);                  
          success = set_vartype_ptr("[N]", buffer, i, &newElement.V, word); CHECKERRQ(MAP_FATAL_32);
          cx = map_snprintf(buffer, 64, "H[%d]", i+1); assert(cx>=0);                            
          success = set_vartype_ptr("[N]", buffer, i, &newElement.H, word); CHECKERRQ(MAP_FATAL_32);
        } else if (cnt==1) {
          success = associate_element_with_cable_property(&newElement, data, word, map_msg, ierr); CHECKERRQ(MAP_FATAL_32);
        } else if (cnt==2) { 
          cx = map_snprintf(buffer, 64, "Lu[%d]", i+1); assert(cx>=0);                
          success = set_vartype("[m]", buffer, i, &newElement.Lu, word); CHECKERRQ(MAP_FATAL_26);
        } else if (cnt==3) {
          success = associate_element_with_anchor_node(data, &newElement, i+1, word, map_msg, ierr); CHECKERRQ(MAP_FATAL_32);        
        } else if (cnt==4) {
          success = associate_element_with_fairlead_node(data, &newElement, i+1, word, map_msg, ierr); CHECKERRQ(MAP_FATAL_32); 
        } else { 
          success = map_set_option_flags (word, &newElement, map_msg, ierr);
        };
      } while (0);
      cnt++;
      word = strtok(NULL, " ,\n\t\r\0");
    };  
    list_append(&data->element, &newElement);
    cnt=0;
  };

  /* check EA, Lu errors. These values cannot be negative. */
  list_iterator_start(&data->element);            /* starting an iteration "session" */
  while (list_iterator_hasnext(&data->element)) { /* tell whether more values available */ 
    elementIter = (Element*)list_iterator_next(&data->element);
    if(elementIter->Lu.value<=0.0) {
      success = MAP_FATAL;
      cx = map_snprintf(buffer, 64, "Line segment %d, Lu = %f [m].", lineCounter, elementIter->Lu.value); assert(cx>=0);
      *ierr = map_set_universal_error(buffer, map_msg, ierr, MAP_FATAL_56);      
    } else if(elementIter->lineProperty->ea<=0.0) {
      success = MAP_FATAL;
      cx = map_snprintf(buffer, 64, "Line segment %d, EA = %f [N].", lineCounter, elementIter->lineProperty->ea); assert(cx>=0);
      *ierr = map_set_universal_error(buffer, map_msg, ierr, MAP_FATAL_57);      
    };
    lineCounter++;
  };
  list_iterator_stop(&data->element); /* ending the iteration "session" */    

  MAPFREE(line);

  if (*ierr==MAP_FATAL) {
    return MAP_FATAL;
  };
  return MAP_SAFE;
};




/**
 * this function returns the size of elements 
 */
size_t cable_element_meter(const void *el) 
{
  /* every element has the constant size of a rectangle structure */
  return sizeof(Element);
};


/**
 * this function returns the size of elements 
 */
size_t cable_library_meter(const void *el) 
{
  /* every element has the constant size of a rectangle structure */
  return sizeof(CableLibrary);
};


/**
 * this function returns the size of elements 
 */
size_t node_meter(const void *el) 
{
  /* every element has the constant size of a rectangle structure */
  return sizeof(Node);
};


/**
 * this function returns the size of elements 
 */
size_t vartype_meter(const void *el) 
{
  /* every element has the constant size of a rectangle structure */
  return sizeof(VarType);
};


/**
 * this function returns the size of elements 
 */
size_t vartype_ptr_meter(const void *el) 
{
  /* every element has the constant size of a rectangle structure */
  return sizeof(VarTypePtr);
};


void set_element_null(Element* ptr)
{
  ptr->lineProperty = NULL;      
  ptr->label = NULL;
  ptr->lineTension  = NULL;
  ptr->anchor = NULL;             /* Anchor node */
  ptr->fairlead = NULL;           /* Fairlead node */

  ptr->psi.name = NULL;
  ptr->psi.units = NULL;
  ptr->alpha.name = NULL;
  ptr->alpha.units = NULL;
  ptr->alphaAtAnchor.name = NULL;
  ptr->alphaAtAnchor.units = NULL;
  ptr->l.name = NULL;
  ptr->l.units = NULL;
  ptr->lb.name = NULL;
  ptr->lb.units = NULL;
  ptr->h.name = NULL;
  ptr->h.units = NULL;
  ptr->H.name = NULL;
  ptr->H.units = NULL;
  ptr->V.name = NULL;
  ptr->V.units = NULL;
  ptr->HAtAnchor.name = NULL;
  ptr->HAtAnchor.units = NULL;
  ptr->VAtAnchor.name = NULL;
  ptr->VAtAnchor.units = NULL;
  ptr->forceAtFairlead.fx.name = NULL;
  ptr->forceAtFairlead.fx.units = NULL;
  ptr->forceAtFairlead.fy.name = NULL;
  ptr->forceAtFairlead.fy.units = NULL;
  ptr->forceAtFairlead.fz.name = NULL;
  ptr->forceAtFairlead.fz.units = NULL;
  ptr->forceAtAnchor.fx.name = NULL;
  ptr->forceAtAnchor.fx.units = NULL;
  ptr->forceAtAnchor.fy.name = NULL;
  ptr->forceAtAnchor.fy.units = NULL;
  ptr->forceAtAnchor.fz.name = NULL;
  ptr->forceAtAnchor.fz.units = NULL;
  ptr->T.name = NULL;
  ptr->T.units = NULL;
  ptr->TAtAnchor.name = NULL;
  ptr->TAtAnchor.units = NULL;
};


/**
 *
 */
void initialize_element(Element* ptr, const int elementNum)
{
  ptr->event.LuRestore = 0.0;
  ptr->event.dLu = 0.0;
  ptr->event.LuMax = 0.0;
  ptr->event.payinFlag = false;
  
  ptr->options.plotFlag = false;
  ptr->options.gxPosFlag = false;
  ptr->options.gyPosFlag = false;
  ptr->options.gzPosFlag = false;
  ptr->options.gxAnchorPosFlag = false;
  ptr->options.gyAnchorPosFlag = false;
  ptr->options.gzAnchorPosFlag = false;
  ptr->options.gxForceFlag = false;
  ptr->options.gyForceFlag = false;
  ptr->options.gzForceFlag = false;
  ptr->options.HFlag = false;
  ptr->options.HAnchorFlag = false;
  ptr->options.VFlag = false;
  ptr->options.VAnchorFlag = false;
  ptr->options.fairleadTensionFlag = false;
  ptr->options.anchorTensionFlag = false;
  ptr->options.horizontalExcursionFlag = false;
  ptr->options.verticalExcursionFlag = false;
  ptr->options.azimuthFlag = false;
  ptr->options.altitudeFlag = false;
  ptr->options.altitudeAnchorFlag = false;
  ptr->options.lineTensionFlag = false;
  ptr->options.omitContact = false;
  ptr->options.layLengthFlag = false;
  ptr->options.damageTimeFlag = false;
  ptr->options.diagnosticsFlag = false;

  initialize_vartype("[deg]", "psi", &ptr->psi, elementNum); 
  initialize_vartype("[deg]", "alpha", &ptr->alpha, elementNum); 
  initialize_vartype("[deg]", "alpha_a", &ptr->alphaAtAnchor, elementNum); 
  initialize_vartype("[m]", "l", &ptr->l, elementNum); 
  initialize_vartype("[m]", "lb", &ptr->lb, elementNum);
  initialize_vartype("[m]", "h", &ptr->h, elementNum);  
  initialize_vartype("[N]", "H_a", &ptr->HAtAnchor, elementNum); 
  initialize_vartype("[N]", "fx", &ptr->forceAtFairlead.fx, elementNum); 
  initialize_vartype("[N]", "fy", &ptr->forceAtFairlead.fy, elementNum);
  initialize_vartype("[N]", "fz", &ptr->forceAtFairlead.fz, elementNum);  
  initialize_vartype("[N]", "fx_a", &ptr->forceAtAnchor.fx, elementNum); 
  initialize_vartype("[N]", "fy_a", &ptr->forceAtAnchor.fy, elementNum); 
  initialize_vartype("[N]", "fz_a", &ptr->forceAtAnchor.fz, elementNum); 
  initialize_vartype("[N]", "V_a" , &ptr->VAtAnchor, elementNum); 
  initialize_vartype("[m]", "T", &ptr->T, elementNum); 
  initialize_vartype("[m]", "T_a", &ptr->TAtAnchor, elementNum); 

  ptr->H.name = NULL;
  ptr->V.name = NULL;
  ptr->H.units = NULL;
  ptr->V.units = NULL;
  
  /* */
  ptr->damageTime = -999.9;
  ptr->diagnosticType = -9999;
  ptr->segmentSize = 10;
};


/**
 *
 */
MAP_ERROR_CODE initialize_external_applied_force(char* unit, char* alias, const int num, VarType* type, char const* property)
{
  type->name = malloc(sizeof(char)*strlen(alias)+1);
  type->units = malloc(sizeof(char)*strlen(unit)+1);
  
  strcpy(type->name, alias);
  strcpy(type->units, unit);

  type->referenceCounter = 0;
  type->id = num;  
  type->value = 0;
  type->isFixed = true;

  if (property[0]=='#') { 
    type->value = 0.0;
  } else { 
    /* this variable is constant */    
    if (is_numeric(property)) { 
      /* cannot convert to numeric value */
      return MAP_FATAL;
    } else { 
      /* converted to numeric value */
      type->value = (MapReal)atof(property);      
    };
  };
  return MAP_SAFE;
};


/**
 *
 */
MAP_ERROR_CODE initialize_node_sum_force_ptr(char* unit, char* alias, const int num, VarTypePtr* type)
{
  type->name = malloc(sizeof(char)*strlen(alias)+1);
  type->units = malloc(sizeof(char)*strlen(unit)+1);  
  strcpy(type->name, alias);
  strcpy(type->units, unit);

  type->referenceCounter = 0;
  type->id = num;  
  type->isFixed = true;
  *(type->value) = 0.0;
  return MAP_SAFE;
};


/**
 *
 */
void set_node_null(Node* ptr)
{
  ptr->positionPtr.x.name = NULL;
  ptr->positionPtr.x.units = NULL;
  ptr->positionPtr.x.value = NULL;
  ptr->positionPtr.y.name = NULL;
  ptr->positionPtr.y.units = NULL;
  ptr->positionPtr.y.value = NULL;
  ptr->positionPtr.z.name = NULL;
  ptr->positionPtr.z.units = NULL;
  ptr->positionPtr.z.value = NULL;
  ptr->MApplied.name = NULL;
  ptr->MApplied.units = NULL;
  ptr->BApplied.name = NULL;
  ptr->BApplied.units = NULL;
  ptr->sumForcePtr.fx.name = NULL;
  ptr->sumForcePtr.fx.units = NULL;
  ptr->sumForcePtr.fx.value = NULL;
  ptr->sumForcePtr.fy.name = NULL;
  ptr->sumForcePtr.fy.units = NULL;
  ptr->sumForcePtr.fy.value = NULL;
  ptr->sumForcePtr.fz.name = NULL;
  ptr->sumForcePtr.fz.units = NULL;
  ptr->sumForcePtr.fz.value = NULL;

  ptr->externalForce.fx.name = NULL;
  ptr->externalForce.fx.units = NULL;
  ptr->externalForce.fy.name = NULL;
  ptr->externalForce.fy.units = NULL;
  ptr->externalForce.fz.name = NULL;
  ptr->externalForce.fz.units = NULL;
};


/**
 *
 */
MAP_ERROR_CODE associate_element_with_cable_property(Element* newElement, ModelData* data, char* word, char* map_msg, MAP_ERROR_CODE* ierr)
{
  CableLibrary* iterCableLibrary = NULL;
  iterCableLibrary = NULL;
  newElement->lineProperty = NULL;

  list_iterator_start(&data->cableLibrary); /* starting an iteration session */
  while (list_iterator_hasnext(&data->cableLibrary)) { /* tell whether more values available */
    iterCableLibrary = (CableLibrary*)list_iterator_next(&data->cableLibrary);
    if (!strcicmp(iterCableLibrary->label, word)) {
      newElement->lineProperty = iterCableLibrary;
      list_iterator_stop(&data->cableLibrary); /* ending the iteration session */  
      break;
    }; 
  };
  list_iterator_stop(&data->cableLibrary); /* ending the iteration session */  
  if (newElement->lineProperty==NULL) {        
    *ierr = map_set_universal_error(word, map_msg, ierr, MAP_FATAL_27);
    return MAP_FATAL;
  };
  return MAP_SAFE;
};


/**
 *
 */
MAP_ERROR_CODE associate_element_with_anchor_node(ModelData* data, Element* new_element, const int element_num, char* word, char* map_msg, MAP_ERROR_CODE* ierr)
{
  Node* iter_node = NULL;
  int node_num = 0;
  char buffer[64] = "";
  int cx = 0;
  int success = 0;

  new_element->anchor = NULL;
  
  success = is_numeric(word); 
  if (success) {
    cx = map_snprintf(buffer, 64, "Element '%d'.", element_num+1);assert(cx>=0);
    *ierr = map_set_universal_error(buffer, map_msg, ierr, MAP_FATAL_28);
    return MAP_FATAL;
  };
  
  node_num = (int)atoi(word); 
  iter_node = (Node*)list_get_at(&data->node, node_num-1);
  if (!iter_node) {
    cx = map_snprintf( buffer, 64, "Element '%d'.", element_num+1); assert(cx>=0);
    *ierr = map_set_universal_error(buffer, map_msg, ierr, MAP_FATAL_30);
    return MAP_FATAL;
  };

  new_element->anchor = iter_node;
  return MAP_SAFE;
};


/**
 *
 */
MAP_ERROR_CODE associate_element_with_fairlead_node(ModelData* data, Element* new_element, const int element_num, char* word, char* map_msg, MAP_ERROR_CODE* ierr)
{
  Node *iter_node = NULL;
  int node_num = 0;
  char buffer[64] = "";
  int cx = 0;
  int success = 0;

  new_element->fairlead = NULL;
  
  success = is_numeric(word); 
  if (success) {
    cx = map_snprintf(buffer, 64, "Element '%d'.", element_num+1); assert(cx>=0);
    *ierr = map_set_universal_error(buffer, map_msg, ierr, MAP_FATAL_29);
    return MAP_FATAL;
  };
  
  node_num = (int)atoi(word); 
  iter_node = (Node*)list_get_at(&data->node, node_num-1);
  if (!iter_node) {
    cx = map_snprintf(buffer, 64, "Element '%d'.", element_num );assert(cx>=0);
    *ierr = map_set_universal_error(buffer, map_msg, ierr, MAP_FATAL_31);
    return MAP_FATAL;
  };
  
  new_element->fairlead = iter_node;

  return MAP_SAFE;
};


/**
 * Create internal state data structures
 *
 * @todo: delete data->z, data->u, data->yList
 * @acceses: none
 * @calledby: mapcall_msqs_init( )
 */
MAP_ERROR_CODE allocate_outlist(ModelData* data, char* map_msg, MAP_ERROR_CODE* ierr)
{
  data->yList = (OutputList*)malloc(sizeof(OutputList)); 
  if (data->yList==NULL) {
    *ierr = map_set_universal_error("", map_msg, ierr, MAP_FATAL_46);    
    return MAP_FATAL;
  };
  return MAP_SAFE;
};


/**
 * Frees internal state data allcoated in the mapcall_msqs_init( ) function
 *
 * @todo: delete additional dependancies in data->z, data->yList, data->u
 * @acceses: none
 * @calledby: mapcall_msqs_end( )
 * @see: allocate_outlist( )
 */
MAP_ERROR_CODE free_outlist(ModelData* data, char* map_msg, MAP_ERROR_CODE* ierr)
{
  MAPFREE(data->yList);
  return MAP_SAFE;
};
