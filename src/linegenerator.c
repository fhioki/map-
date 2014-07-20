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
#include "pyprotos.h"


/**
 *
 */
MAP_EXTERNCALL double* pyget_plot_x(MAP_OtherStateType_t* otherType, int i, int numPlotPoints, char *map_msg, MAP_ERROR_CODE *ierr)
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
  MapReal fairleadX=0.0;
  MapReal anchorX=0.0;
  MapReal cb = 0.0;
  MapReal lambda = 0.0;
  double* arrayX = NULL;
  char buffer[64] = "";
  int cx = 0;
  int s = 0;

  *ierr = MAP_SAFE;
  map_msg[0] = 0;
  element = (Element*)list_get_at(&data->element, i);

  if (element==NULL) {    
    cx = map_snprintf(buffer, 64, "Element out of range: %d.", i); assert(cx>=0);
    *ierr=map_set_universal_error(buffer, map_msg, ierr, MAP_FATAL_42);
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

    /*
      If the cable is not resting on the seabed, we use the classic catenary equation
      for a hanging chain to plot the mooring line profile. Otherwise if it is, we 
      the modified version as done in the FAST wind turbine program. 
          
      @ref : J. Jonkman, November 2007. "Dynamic Modeling and Loads Analysis of an 
             Offshore Floating Wind Turbine." NREL Technical Report NREL/TP-500-41958.
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


/**
 *
 */
MAP_EXTERNCALL double* pyget_plot_y(MAP_OtherStateType_t* otherType, int i, int numPlotPoints, char* map_msg, MAP_ERROR_CODE* ierr)
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

  *ierr = MAP_SAFE;
  map_msg[0] = 0;
  element = (Element*)list_get_at(&data->element, i);
  
  if (element==NULL) {
    cx = map_snprintf(buffer, 64, "Element out of range: %d.", i); assert(cx>=0);
    *ierr = map_set_universal_error(buffer, map_msg, ierr, MAP_FATAL_42);
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

    /*
      If the cable is not resting on the seabed, we use the classic catenary equation
      for a hanging chain to plot the mooring line profile. Otherwise if it is, we 
      the modified version as done in the FAST wind turbine program. 
          
      @ref : J. Jonkman, November 2007. "Dynamic Modeling and Loads Analysis of an 
             Offshore Floating Wind Turbine." NREL Technical Report NREL/TP-500-41958.
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


/**
 *
 */
MAP_EXTERNCALL double* pyget_plot_z(MAP_OtherStateType_t* otherType, int i, int numPlotPoints, char* map_msg, MAP_ERROR_CODE* ierr)
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
  int cx = 0;
  int s = 0;

  *ierr = MAP_SAFE;
  map_msg[0] = 0;
  element = (Element*)list_get_at(&data->element, i);
  
  if (element==NULL){
    cx = map_snprintf(buffer, 64, "Element out of range: %d.", i); assert(cx>=0);
    *ierr = map_set_universal_error(buffer, map_msg, ierr, MAP_FATAL_42);
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

    /*
      If the cable is not resting on the seabed, we use the classic catenary equation
      for a hanging chain to plot the mooring line profile. Otherwise if it is, we 
      the modified version as done in the FAST wind turbine program. 
          
      @ref : J. Jonkman, November 2007. "Dynamic Modeling and Loads Analysis of an 
             Offshore Floating Wind Turbine." NREL Technical Report NREL/TP-500-41958.
    */        
    if (element->options.omitContact==true || w<0.0 || (V-w*Lu)>0.0) { /* true when no portion of the line rests on the seabed */
      for (s=0 ; s<numPlotPoints ; s++) {
        arrayZ[s] =  fairleadZ - ((H/w)*(sqrt(1+pow(V/H,2)) - sqrt(1+pow((V-w*S)/H,2))) + (1/EA)*(V*S+w*S*S/2)); // Z position of element in global coordinates 
        S += dS;
      };
    } else {
      Lb = Lu - (V/w);      
      for (s=0 ; s<numPlotPoints ; s++) {        
        if (0<=S && S<=Lb) {           
          arrayZ[s] = anchorZ; 
        } else {        
          arrayZ[s] = ((H/w)*(sqrt(1 + pow((w*(S-Lb)/H),2)) - 1) + ((w*pow((S-Lb),2))/(2*EA))) + anchorZ; // @todo: verify this equation
        };
        S += dS;
      };      
    };
  };
  return arrayZ;
};


/**
 *
 */
MAP_EXTERNCALL void pyget_plot_array_free(MapReal* array) 
{
  MAPFREE(array);
}
