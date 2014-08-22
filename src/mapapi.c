/****************************************************************
 *   Copyright (C) 2014 mdm                                     *
 *   marco[dot]masciola[at]gmail                                *
 *                                                              *
 * Licensed to the Apache Software Foundation (ASF) under one   *
 * or more contributor license agreements.  See the NOTICE file *
 * distributed with this work for additional information        *
 * regarding copyright ownership.  The ASF licenses this file   *
 * to you under the Apache License, Version 2.0 (the            *
 * "License"); you may not use this file except in compliance   *
 * with the License.  You may obtain a copy of the License at   *
 *                                                              *
 *   http://www.apache.org/licenses/LICENSE-2.0                 *
 *                                                              *
 * Unless required by applicable law or agreed to in writing,   *
 * software distributed under the License is distributed on an  *
 * "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY       *
 * KIND, either express or implied.  See the License for the    *
 * specific language governing permissions and limitations      *      
 * under the License.                                           *  
 ****************************************************************/


#include "map.h"
#include "maperror.h"
#include "MAP_Types.h"
#include "mapapi.h"
#include "lineroutines.h"


extern const char MAP_ERROR_STRING[][1024];


MAP_EXTERNCALL void map_offset_vessel(MAP_OtherStateType_t* other_type, MAP_InputType_t* u_type, double x, double y, double z, double phi, double the, double psi, char* map_msg, MAP_ERROR_CODE* ierr)
{
  ModelData* data = other_type->object;
  Vessel* vessel = &data->vessel;
  int i = 0;
  const int u_size = u_type->x_Len;
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

  /* define transformation matrix */
  R[0][0] = cpsi*cthe;    R[0][1] = cpsi*sthe*sphi - spsi*cphi;   R[0][2] = cpsi*sthe*cphi + spsi*sphi;
  R[1][0] = sphi*cthe;    R[1][1] = sphi*sthe*sphi + cpsi*cphi;   R[1][2] = spsi*sthe*cphi - cpsi*sphi;
  R[2][0] = -sthe;        R[2][1] = cthe*sphi;                    R[2][2] = cthe*cphi;

  for (i=0 ; i<u_size ; i++) { 
    /* @todo: need to include the reference position for non-zero reference origins , i.e. r = (xi-ref) 
     *        xi, yi, and zi are the original node position. We are adding the new displacement to it.
     */
    rx = vessel->xi[i];
    ry = vessel->yi[i];
    rz = vessel->zi[i];

    /* matrix-vector product */
    u_type->x[i] = x + rx*R[0][0] + ry*R[0][1] + rz*R[0][2];
    u_type->y[i] = y + rx*R[1][0] + ry*R[1][1] + rz*R[1][2];
    u_type->z[i] = z + rx*R[2][0] + ry*R[2][1] + rz*R[2][2];
  };
};



MAP_EXTERNCALL double** map_linearize_matrix(MAP_InputType_t* u_type, MAP_OtherStateType_t* other_type, MAP_OutputType_t* y_type, MAP_ConstraintStateType_t* z_type, double epsilon, MAP_ERROR_CODE* ierr, char* map_msg)
{
  checkpoint();
  MapReal* x_original = NULL;
  MapReal* y_original = NULL;
  MapReal* z_original = NULL;
  ModelData* data = other_type->object;
  MAP_ERROR_CODE success = MAP_SAFE;
  const int n = u_type->x_Len;
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

   force.fx = (double*)malloc(n*sizeof(double));
   force.fy = (double*)malloc(n*sizeof(double));
   force.fz = (double*)malloc(n*sizeof(double));
   force.mx = (double*)malloc(n*sizeof(double));
   force.my = (double*)malloc(n*sizeof(double));
   force.mz = (double*)malloc(n*sizeof(double));  
   x_original = (double*)malloc(n*sizeof(double));
   y_original = (double*)malloc(n*sizeof(double));
   z_original = (double*)malloc(n*sizeof(double));
   
   /* initialize stuff allocated above to zero */
   for (i=0 ; i<n ; i++) {
    force.fx[i] = 0.0;
    force.fy[i] = 0.0;
    force.fz[i] = 0.0;
    force.mx[i] = 0.0;
    force.my[i] = 0.0;
    force.mz[i] = 0.0;
    x_original[i] = 0.0;
    y_original[i] = 0.0;
    z_original[i] = 0.0;
   };
    
   do {    
     /* first get the original values for the displacements */
     for (k=0 ; k<n ; k++) {
       x_original[k] = u_type->x[k];
       y_original[k] = u_type->y[k];
       z_original[k] = u_type->z[k];      
     };
    
     for (i=0 ; i<SIX ; i++) { /* down, force direction changes */
       success = reset_force_to_zero(force.fx, force.fy, force.fz, force.mx, force.my, force.mz, n);
       if (i==0) {        
         success = fd_x_sequence(other_type, u_type, y_type, z_type, &force, epsilon, n, x_original, map_msg, ierr); CHECKERRQ(MAP_FATAL_62);
         success = calculate_stiffness_2(K[0], &force, epsilon, n); CHECKERRQ(MAP_FATAL_62);
//         // printf("This is the stiffness: %12.2f, %12.2f, %12.2f, %12.2f, %12.2f, %12.2f\n", K[0][0], K[0][1], K[0][2], K[0][3], K[0][4], K[0][5]);
       } else if (i==1) {
         success = fd_y_sequence(other_type, u_type, y_type, z_type, &force, epsilon, n, y_original, map_msg, ierr); CHECKERRQ(MAP_FATAL_63);
         success = calculate_stiffness_2(K[1], &force, epsilon, n); CHECKERRQ(MAP_FATAL_63);
//         // printf("This is the stiffness: %12.2f, %12.2f, %12.2f, %12.2f, %12.2f, %12.2f\n", K[1][0], K[1][1], K[1][2], K[1][3], K[1][4], K[1][5]);
       } else if (i==2) {
         success = fd_z_sequence(other_type, u_type, y_type, z_type, &force, epsilon, n, z_original, map_msg, ierr); CHECKERRQ(MAP_FATAL_64);
         success = calculate_stiffness_2(K[2], &force, epsilon, n); CHECKERRQ(MAP_FATAL_64);
//         // printf("This is the stiffness: %12.2f, %12.2f, %12.2f, %12.2f, %12.2f, %12.2f\n", K[2][0], K[2][1], K[2][2], K[2][3], K[2][4], K[2][5]);
       } else if (i==3) {
         success = fd_phi_sequence(other_type, u_type, y_type, z_type, &force, epsilon, n, x_original, y_original, z_original, map_msg, ierr); //CHECKERRQ(MAP_FATAL_65);
         success = calculate_stiffness_2(K[3], &force, epsilon, n); CHECKERRQ(MAP_FATAL_65);
//         // printf("This is the stiffness: %12.2f, %12.2f, %12.2f, %12.2f, %12.2f, %12.2f\n", K[3][0], K[3][1], K[3][2], K[3][3], K[3][4], K[3][5]);
       } else if (i==4) {
         success = fd_the_sequence(other_type, u_type, y_type, z_type, &force, epsilon, n, x_original, y_original, z_original, map_msg, ierr); //CHECKERRQ(MAP_FATAL_66);
         success = calculate_stiffness_2(K[4], &force, epsilon, n); CHECKERRQ(MAP_FATAL_66);
//         // printf("This is the stiffness: %12.2f, %12.2f, %12.2f, %12.2f, %12.2f, %12.2f\n", K[4][0], K[4][1], K[4][2], K[4][3], K[4][4], K[4][5]);
       } else if (i==5) {
         success = fd_psi_sequence(other_type, u_type, y_type, z_type, &force, epsilon, n, x_original, y_original, z_original, map_msg, ierr); //CHECKERRQ(MAP_FATAL_67);
         success = calculate_stiffness_2(K[5], &force, epsilon, n); CHECKERRQ(MAP_FATAL_67);
//         // printf("This is the stiffness: %12.2f, %12.2f, %12.2f, %12.2f, %12.2f, %12.2f\n", K[5][0], K[5][1], K[5][2], K[5][3], K[5][4], K[5][5]);
       };
     };
   } while (0);  
 
   success = reset_force_to_zero(force.fx, force.fy, force.fz, force.mx, force.my, force.mz, n);
   success = restore_original_displacement(u_type->x, x_original, n);
   success = restore_original_displacement(u_type->y, y_original, n);
   success = restore_original_displacement(u_type->z, z_original, n);
   success = line_solve_sequence(data, 0.0, map_msg, ierr); 
   
  MAPFREE(force.fx);
  MAPFREE(force.fy);
  MAPFREE(force.fz);
  MAPFREE(force.mx);
  MAPFREE(force.my);
  MAPFREE(force.mz);
  MAPFREE(x_original);
  MAPFREE(y_original);
  MAPFREE(z_original);  
  return K;
};


MAP_EXTERNCALL void map_free_linearize_matrix(double** array)
{
  int i = 0;
  for(i=0 ; i<6 ; i++) {
    MAPFREE(array[i]);
  };
  MAPFREE(array);
};


MAP_EXTERNCALL double* map_plot_x_array(MAP_OtherStateType_t* otherType, int i, int numPlotPoints, char *map_msg, MAP_ERROR_CODE *ierr)
{
  ModelData* data = otherType->object;
  const Element* element = (Element*)list_get_at(&data->element, i);  
  MapReal H = 0.0;
  MapReal V = 0.0;
  MapReal EA = 0.0;
  MapReal Lu = 0.0;
  MapReal Lb = 0.0;
  MapReal w = 0.0;
  MapReal S = 0.0;
  MapReal dS = 0.0;
  MapReal fairleadX=0.0;
  MapReal anchorX=0.0;
  MapReal cb = 0.0;
  MapReal lambda = 0.0;
  double* arrayX = NULL;
  char buffer[64] = "";
  int ret = 0;
  int cx = 0;
  int s = 0;
  bstring user_msg = NULL;

  *ierr = MAP_SAFE;
  map_msg[0] = 0;
  // element = (Element*)list_get_at(&data->element, i);

  if (element==NULL) {    
    user_msg = bformat("Element out of range: <%d>.", i);
    *ierr = map_set_universal_error(user_msg, map_msg, *ierr, MAP_FATAL_42);        
    ret = bdestroy(user_msg); 
    user_msg = NULL;
  } else {
    fairleadX = *(element->fairlead->positionPtr.x.value);
    anchorX = *(element->anchor->positionPtr.x.value);
    arrayX = (double*)malloc(numPlotPoints*sizeof(double));
    H = *(element->H.value);
    V = *(element->V.value);  
    EA = element->lineProperty->ea;
    Lu = element->Lu.value;   
    w = element->lineProperty->omega;
    cb = element->lineProperty->cb;
    dS = Lu/(MapReal)(numPlotPoints-1) ;
    
    /* If the cable is not resting on the seabed, we use the classic catenary equation
     * for a hanging chain to plot the mooring line profile. Otherwise if it is, we 
     * the modified version as done in the FAST wind turbine program. 
     *     
     * @ref : J. Jonkman, November 2007. "Dynamic Modeling and Loads Analysis of an 
     *        Offshore Floating Wind Turbine." NREL Technical Report NREL/TP-500-41958.
     */        
    if (element->options.omitContact==true || w<0.0 || (V-w*Lu)>0.0) { /* true when no portion of the line rests on the seabed */
      for (s=0 ; s<numPlotPoints ; s++) {
        arrayX[s] = fairleadX - ( (H/w)*ARCSINH( V/H ) - (H/w)*ARCSINH( (V-S*w)/H ) + (H*S)/(EA) )*cos(element->psi.value);
        S += dS;
      };
    } else {
      Lb = Lu - (V/w);      
      lambda = (Lb-H/(cb*w))>0 ? (Lb-H/(cb*w)) : 0; 
      for (s=0 ; s<numPlotPoints ; s++) {
        if (0<=S && S<=(Lb-H/(cb*w))) { /* for 0 <= s <= Lb - H/(Cb*w) */
          arrayX[s] = S*cos(element->psi.value) + anchorX; 
        } else if ((Lb-H/(cb/w))<S && S<=Lb) { /* for Lb - H/(Cb*w) < s <= Lb */
          arrayX[s] = (S + ((cb*w)/(2*EA)) * (S*S - 2*(Lb-H/(cb*w))*S + (Lb- H/(cb*w))*lambda))*cos(element->psi.value) + anchorX; 
        } else { /* for Lb < s <= L */
          arrayX[s] = (Lb + (H/w)*ARCSINH((w*(S-Lb))/H))*cos(element->psi.value)    
            - ( ((H*S)/(EA)) + ((cb*w)/(2*EA))*(-Lb*Lb + (Lb-H/(cb*w))*lambda))*cos(element->psi.value) + anchorX; 
        };
        S += dS;
      };
    };
  };
  return arrayX;
};


MAP_EXTERNCALL double* map_plot_y_array(MAP_OtherStateType_t* otherType, int i, int numPlotPoints, char* map_msg, MAP_ERROR_CODE* ierr)
{
  ModelData* data = otherType->object;
  Element* element = NULL;
  MapReal H = 0.0;
  MapReal V = 0.0;
  MapReal EA = 0.0;
  MapReal Lu = 0.0;
  MapReal Lb = 0.0;
  MapReal w = 0.0;
  MapReal S = 0.0;
  MapReal dS = 0.0;
  MapReal fairleadY=0.0;
  MapReal anchorY=0.0;
  MapReal cb = 0.0;
  MapReal lambda = 0.0;
  double* arrayY = NULL;
  char buffer[64] = "";
  int cx = 0;
  int s = 0;
  int ret = 0;
  bstring user_msg = NULL;
  
  *ierr = MAP_SAFE;
  map_msg[0] = 0;
  element = (Element*)list_get_at(&data->element, i);
  
  if (element==NULL) {
    user_msg = bformat("Element out of range: <%d>.", i);
    *ierr = map_set_universal_error(user_msg, map_msg, *ierr, MAP_FATAL_42);        
    ret = bdestroy(user_msg); 
    user_msg = NULL;
    //cx = map_snprintf(buffer, 64, "Element out of range: %d.", i); assert(cx>=0);
    //*ierr = map_set_universal_error(buffer, map_msg, ierr, MAP_FATAL_42);
  } else {    
    fairleadY = *(element->fairlead->positionPtr.y.value);
    anchorY = *(element->anchor->positionPtr.y.value);
    arrayY = (double*)malloc(numPlotPoints*sizeof(double));
    H = *(element->H.value);
    V = *(element->V.value);  
    EA = element->lineProperty->ea;
    Lu = element->Lu.value;
    w = element->lineProperty->omega;
    cb = element->lineProperty->cb;
    dS = Lu/(MapReal)(numPlotPoints-1) ;
    
    /* If the cable is not resting on the seabed, we use the classic catenary equation
     * for a hanging chain to plot the mooring line profile. Otherwise if it is, we 
     * the modified version as done in the FAST wind turbine program. 
     *     
     * @ref : J. Jonkman, November 2007. "Dynamic Modeling and Loads Analysis of an 
     *        Offshore Floating Wind Turbine." NREL Technical Report NREL/TP-500-41958.
     */        
    if (element->options.omitContact==true || w<0.0 || (V-w*Lu)>0.0) { /* true when no portion of the line rests on the seabed */
      for (s=0 ; s<numPlotPoints ; s++) {
        arrayY[s] = fairleadY - ( (H/w)*ARCSINH( V/H ) - (H/w)*ARCSINH( (V-S*w)/H ) + (H*S)/(EA) )*sin(element->psi.value);
        S += dS; 
      };
    } else {
      Lb = Lu - (V/w);      
      lambda = (Lb-H/(cb*w))>0 ? (Lb-H/(cb*w)) : 0; 
      for (s=0 ; s<numPlotPoints ; s++) {
        if (0<=S && S<=(Lb-H/(cb*w))) { /* for 0 <= s <= Lb - H/(Cb*w) */
          arrayY[s] = S*sin(element->psi.value) + anchorY; 
        } else if ((Lb-H/(cb/w))<S && S<=Lb) { /* for Lb - H/(Cb*w) < s <= Lb */
          arrayY[s] = (S + ((cb*w)/(2*EA))*(S*S - 2*(Lb-H/(cb*w))*S + (Lb- H/(cb*w))*lambda))*sin(element->psi.value) + anchorY; 
        } else { /* for Lb < s <= L */
          arrayY[s] = (Lb + (H/w)*ARCSINH((w*(S-Lb))/H))*sin(element->psi.value)    
            - (((H*S)/(EA)) + ((cb*w)/(2*EA))*(-Lb*Lb + (Lb-H/(cb*w))*lambda))*sin(element->psi.value) + anchorY; 
        };
        S += dS;
      };
    };
  };
  return arrayY;
};


MAP_EXTERNCALL double* map_plot_z_array(MAP_OtherStateType_t* otherType, int i, int numPlotPoints, char* map_msg, MAP_ERROR_CODE* ierr)
{
  ModelData* data = otherType->object;
  Element* element = NULL;
  MapReal H = 0.0;
  MapReal V = 0.0;
  MapReal Va = 0.0;
  MapReal EA = 0.0;
  MapReal Lu = 0.0;
  MapReal Lb = 0.0;
  MapReal w = 0.0;
  MapReal S = 0.0;
  MapReal dS = 0.0;
  MapReal fairleadZ = 0.0;
  MapReal anchorZ = 0.0;
  double* arrayZ = NULL;
  char buffer[64] = "";
  int ret = 0;
  int cx = 0;
  int s = 0;
  bstring user_msg = NULL;

  *ierr = MAP_SAFE;
  map_msg[0] = 0;
  element = (Element*)list_get_at(&data->element, i);
  
  if (element==NULL){
    user_msg = bformat("Element out of range: <%d>.", i);
    *ierr = map_set_universal_error(user_msg, map_msg, *ierr, MAP_FATAL_42);        
    ret = bdestroy(user_msg); 
    user_msg = NULL;
  } else {
    fairleadZ = *(element->fairlead->positionPtr.z.value);
    anchorZ = *(element->anchor->positionPtr.z.value);    
    arrayZ = (double*)malloc(numPlotPoints*sizeof(double));
    H = *(element->H.value);
    V = *(element->V.value);  
    EA = element->lineProperty->ea;
    Lu = element->Lu.value;
    w = element->lineProperty->omega;
    dS = Lu/(MapReal)(numPlotPoints-1) ;

    /* If the cable is not resting on the seabed, we use the classic catenary equation
     * for a hanging chain to plot the mooring line profile. Otherwise if it is, we 
     * the modified version as done in the FAST wind turbine program. 
     *     
     * @ref : J. Jonkman, November 2007. "Dynamic Modeling and Loads Analysis of an 
     *        Offshore Floating Wind Turbine." NREL Technical Report NREL/TP-500-41958.
     */        
    if (element->options.omitContact==true || w<0.0 || (V-w*Lu)>0.0) { /* true when no portion of the line rests on the seabed */
      for (s=0 ; s<numPlotPoints ; s++) {
        arrayZ[s] =  fairleadZ - ((H/w)*(sqrt(1+pow(V/H,2)) - sqrt(1+pow((V-w*S)/H,2))) + (1/EA)*(V*S+w*S*S/2)); /* Z position of element in global coordinates */
        S += dS;
      };
    } else {
      Lb = Lu - (V/w);      
      for (s=0 ; s<numPlotPoints ; s++) {        
        if (0<=S && S<=Lb) {           
          arrayZ[s] = anchorZ; 
        } else {        
          /* @todo: verify this equation before someone uses this program to  design something that matters */           
          arrayZ[s] = ((H/w)*(sqrt(1 + pow((w*(S-Lb)/H),2)) - 1) + ((w*pow((S-Lb),2))/(2*EA))) + anchorZ; 

        };
        S += dS;
      };      
    };
  };
  return arrayZ;
};


MAP_EXTERNCALL void map_plot_array_free(MapReal* array) 
{
  MAPFREE(array);
}


// MAP_EXTERNCALL MapReal map_residual_function_length(MAP_OtherStateType_t* otherType, int i, char* map_msg, MAP_ERROR_CODE* ierr)
// { 
//   ModelData* data = otherType->object;
//   Element* element = NULL;
//   MapReal Fh = 0.0;
//   MapReal Fv = 0.0;
//   MapReal EA = 0.0;
//   MapReal Lu = 0.0;
//   MapReal length = 0.0;
//   MapReal omega = 0.0;
//   MapReal cb = 0.0;
//   bool contactFlag = false;
//   char buffer[64] = "";
//   int cx = 0;
//   
//   *ierr = MAP_SAFE;
//   map_msg[0] = 0;
//   element = (Element*)list_get_at(&data->element, i);
// 
//   if (element==NULL) {
//     cx = map_snprintf(buffer, 64, "Element out of range: %d.", i); assert(cx>=0);
//     *ierr = map_set_universal_error(buffer, map_msg, ierr, MAP_FATAL_42);
//     return -999.9;
//   };
// 
//   Fh = *(element->H.value);
//   Fv = *(element->V.value);  
//   EA = element->lineProperty->ea;
//   Lu = element->Lu.value;
//   length = element->l.value;
//   omega = element->lineProperty->omega;
//   contactFlag = element->options.omitContact;
//   cb = element->lineProperty->cb;
// 
//   if (contactFlag==true || omega<0.0 || (Fv-omega*Lu)>0.0) { /* true when no portion of the line rests on the seabed */
//     return residual_function_length_no_contact(Fv, Fh, omega, Lu, EA, length);
//   } else { /* true when a portion of the line rests on the seabed and the anchor tension is nonzero */
//     return residual_function_length_contact(Fv, Fh, omega, Lu, EA, length, cb);
//   };
// };


// MAP_EXTERNCALL MapReal map_residual_function_height(MAP_OtherStateType_t* otherType, int i, char* map_msg, MAP_ERROR_CODE* ierr)
// {
//   ModelData* data = otherType->object;
//   Element* element = NULL;
//   MapReal Fh = 0.0;
//   MapReal Fv = 0.0;
//   MapReal EA = 0.0;
//   MapReal Lu = 0.0;
//   MapReal height = 0.0;
//   MapReal omega = 0.0;
//   MapReal cb = 0.0;
//   bool contactFlag = false;
//   char buffer[64] = "";
//   int cx = 0;
// 
//   *ierr = MAP_SAFE;
//   map_msg[0] = 0;
//   element = (Element*)list_get_at( &data->element, i );
// 
//   if (element==NULL) {    
//     cx = map_snprintf(buffer, 64, "Element out of range: %d.", i); assert(cx>=0);
//     *ierr = map_set_universal_error( buffer, map_msg, ierr, MAP_FATAL_42);
//     return -999.9;
//   };
// 
//   Fh = *(element->H.value);
//   Fv = *(element->V.value);  
//   EA = element->lineProperty->ea;
//   Lu = element->Lu.value;
//   height = element->h.value;
//   omega = element->lineProperty->omega;
//   contactFlag = element->options.omitContact;
//   cb = element->lineProperty->cb;
// 
//   if (contactFlag==true || omega<0.0 || (Fv-omega*Lu)>0.0) { /* true when no portion of the line rests on the seabed */
//     return residual_function_height_no_contact(Fv, Fh, omega, Lu, EA, height);
//   } else { /* true when a portion of the line rests on the seabed and the anchor tension is nonzero */
//     return residual_function_height_contact(Fv, Fh, omega, Lu, EA, height, cb);
//   };
// };


// MAP_EXTERNCALL MapReal map_jacobian_dxdh(MAP_OtherStateType_t* otherType, int i, char* map_msg, MAP_ERROR_CODE* ierr)
// {
//   ModelData* data = otherType->object;
//   Element* element = NULL;
//   MapReal Fh = 0.0;
//   MapReal Fv = 0.0;
//   MapReal EA = 0.0;
//   MapReal Lu = 0.0;
//   MapReal omega = 0.0;
//   MapReal cb = 0.0;
//   bool contactFlag = false;
//   char buffer[64] = "";
//   int cx = 0;
// 
//   *ierr = MAP_SAFE;
//   map_msg[0] = 0;
//   element = (Element*)list_get_at(&data->element, i);
// 
//   if (element==NULL) {    
//     cx = map_snprintf(buffer, 64, "Element out of range: %d.", i); assert(cx>=0);
//     *ierr = map_set_universal_error(buffer, map_msg, ierr, MAP_FATAL_42);
//     return -999.9;
//   };
// 
//   Fh = *(element->H.value);
//   Fv = *(element->V.value);  
//   EA = element->lineProperty->ea;
//   Lu = element->Lu.value;
//   omega = element->lineProperty->omega;
//   contactFlag = element->options.omitContact;
//   cb = element->lineProperty->cb;
// 
//   if (contactFlag==true || omega<0.0 || (Fv-omega*Lu)>0.0) { /* true when no portion of the line rests on the seabed */
//     return jacobian_dxdh_no_contact(Fv, Fh, omega, Lu, EA);
//   } else { /* true when a portion of the line rests on the seabed and the anchor tension is nonzero */
//     return jacobian_dxdh_contact(Fv, Fh, omega, Lu, EA, cb);
//   };
// };


// MAP_EXTERNCALL MapReal map_jacobian_dxdv(MAP_OtherStateType_t* otherType, int i, char* map_msg, MAP_ERROR_CODE* ierr)
// {
//   ModelData* data = otherType->object;
//   Element* element = NULL;
//   MapReal Fh = 0.0;
//   MapReal Fv = 0.0;
//   MapReal EA = 0.0;
//   MapReal Lu = 0.0;
//   MapReal omega = 0.0;
//   MapReal cb = 0.0;
//   bool contactFlag = false;
//   char buffer[64] = "";
//   int cx = 0;
// 
//   *ierr = MAP_SAFE;
//   map_msg[0] = 0;
//   element = (Element*)list_get_at(&data->element, i);
// 
//   if (element==NULL) {    
//     cx = map_snprintf(buffer, 64, "Element out of range: %d.", i); assert(cx>=0);
//     *ierr = map_set_universal_error(buffer, map_msg, ierr, MAP_FATAL_42);
//     return -999.9;
//   };
// 
//   Fh = *(element->H.value);
//   Fv = *(element->V.value);  
//   EA = element->lineProperty->ea;
//   Lu = element->Lu.value;
//   omega = element->lineProperty->omega;
//   contactFlag = element->options.omitContact;
//   cb = element->lineProperty->cb;
// 
//   if (contactFlag==true || omega<0.0 || (Fv-omega*Lu)>0.0) { /* true when no portion of the line rests on the seabed */
//     return jacobian_dxdv_no_contact(Fv, Fh, omega, Lu, EA);
//   } else { /* true when a portion of the line rests on the seabed and the anchor tension is nonzero */
//     return jacobian_dxdv_contact(Fv, Fh, omega, Lu, EA, cb);
//   };
// };


// MAP_EXTERNCALL MapReal map_jacobian_dzdh(MAP_OtherStateType_t* otherFortType, int i, char* map_msg, MAP_ERROR_CODE* ierr)
// {
//   ModelData* data = otherFortType->object;
//   Element* element = NULL;
//   MapReal Fh = 0.0;
//   MapReal Fv = 0.0;
//   MapReal EA = 0.0;
//   MapReal Lu = 0.0;
//   MapReal omega = 0.0;
//   MapReal cb = 0.0;
//   bool contactFlag = false;
//   char buffer[64] = "";
//   int cx = 0;
// 
//   *ierr = MAP_SAFE;
//   map_msg[0] = 0;
//   element = (Element*)list_get_at(&data->element, i);
// 
//   if (element==NULL) {    
//     cx = map_snprintf(buffer, 64, "Element out of range: %d.", i ); assert(cx>=0);
//     *ierr = map_set_universal_error(buffer, map_msg, ierr, MAP_FATAL_42);
//     return -999.9;
//   };
// 
//   Fh = *(element->H.value);
//   Fv = *(element->V.value);  
//   EA = element->lineProperty->ea;
//   Lu = element->Lu.value;
//   omega = element->lineProperty->omega;
//   contactFlag = element->options.omitContact;
//   cb = element->lineProperty->cb;
//   
//   if (contactFlag==true || omega<0.0 || (Fv-omega*Lu)>0.0) { /* true when no portion of the line rests on the seabed */
//     return jacobian_dzdh_no_contact(Fv, Fh, omega, Lu, EA);
//   } else { /* true when a portion of the line rests on the seabed and the anchor tension is nonzero */
//     return jacobian_dzdh_contact(Fv, Fh, omega, Lu, EA, cb);
//   };
// };


// MAP_EXTERNCALL MapReal map_jacobian_dzdv(MAP_OtherStateType_t* otherFortType, int i, char* map_msg, MAP_ERROR_CODE* ierr)
// {
//   ModelData* data = otherFortType->object;
//   Element* element = NULL;
//   MapReal Fh = 0.0;
//   MapReal Fv = 0.0;
//   MapReal EA = 0.0;
//   MapReal Lu = 0.0;
//   MapReal omega  = 0.0;
//   MapReal cb = 0.0;
//   char buffer[64] = "";
//   bool contactFlag = false;
//   int cx = 0;
//  
//   *ierr = MAP_SAFE;
//   map_msg[0] = 0;
//   element = (Element*)list_get_at(&data->element, i);
// 
//   if (element==NULL) {    
//     cx = map_snprintf(buffer, 64, "Element out of range: %d.", i); assert(cx>=0);
//     *ierr = map_set_universal_error(buffer, map_msg, ierr, MAP_FATAL_42);
//     return -999.9;
//   };
// 
//   Fh = *(element->H.value);
//   Fv = *(element->V.value);  
//   EA = element->lineProperty->ea;
//   Lu = element->Lu.value;
//   omega  = element->lineProperty->omega;
//   contactFlag = element->options.omitContact;
//   cb = element->lineProperty->cb;
// 
//   if (contactFlag==true || omega<0.0 || (Fv-omega*Lu)>0.0) { /* true when no portion of the line rests on the seabed */
//     return jacobian_dzdv_no_contact(Fv, Fh, omega, Lu, EA);
//   } else { /* true when a portion of the line rests on the seabed and the anchor tension is nonzero */
//     return jacobian_dzdv_contact(Fv, Fh, omega, Lu, EA, cb);
//   };
// };


// MAP_EXTERNCALL void map_get_fairlead_force_2d(double* H, double* V, MAP_OtherStateType_t* other_type, int index, char* map_msg, MAP_ERROR_CODE* ierr)
// {
//   Element* iter_element = NULL;
//   ModelData* model_data = other_type->object;
//   if (index<=list_size(&model_data->element)-1) {
//     iter_element = (Element*)list_get_at(&model_data->element, index);
//     *H = *(iter_element->H.value);
//     *V = *(iter_element->V.value);
//   } else {
//     /* @todo: throw error: element out of range */
//   };
// }


// MAP_EXTERNCALL void map_get_fairlead_force_3d(double* fx, double* fy, double* fz, MAP_OtherStateType_t* other_type, int index, char* map_msg, MAP_ERROR_CODE* ierr)
// {
//   Element* iter_element = NULL;
//   ModelData* model_data = other_type->object;
//   double psi = 0.0;
//   if (index<=list_size(&model_data->element)-1) {
//     iter_element = (Element*)list_get_at(&model_data->element, index);
//     psi = iter_element->psi.value;
//     *fx = *(iter_element->H.value)*cos(psi);
//     *fy = *(iter_element->H.value)*sin(psi);
//     *fz = *(iter_element->V.value);
//   } else {
//     /* @todo: throw error: element out of range */
//   };
// }


MAP_EXTERNCALL int map_size_elements(MAP_OtherStateType_t* other_type, MAP_ERROR_CODE* ierr, char* map_msg)
{
  ModelData* model_data = other_type->object;
  return list_size(&model_data->element);
}
