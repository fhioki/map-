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
#include "cminpack/minpack.h"


int inner_function_evals(void* elementPtr, int m, int n, const __cminpack_real__* x, __cminpack_real__* fvec, __cminpack_real__* fjac, int ldfjac, int iflag) 
{
  Element* element = (Element*)elementPtr;
  const MapReal Fh = x[0];
  const MapReal Fv = x[1];  
  const MapReal EA = element->lineProperty->ea;
  const MapReal Lu = element->Lu.value;
  const MapReal height = element->h.value;
  const MapReal length = element->l.value;
  const MapReal omega = element->lineProperty->omega;
  const MapReal cb = element->lineProperty->cb;
  const bool contactFlag = element->options.omitContact;
  
  if (iflag==0) {
    return 0;
  };

  if (iflag!=2) {
    if (contactFlag==true || omega<0.0 || (Fv-omega*Lu)>0.0) { /* true when no portion of the line rests on the seabed */
      fvec[0] = residual_function_length_no_contact(Fv, Fh, omega, Lu, EA, length);
      fvec[1] = residual_function_height_no_contact(Fv, Fh, omega, Lu, EA, height); 
    } else { /* true when a portion of the line rests on the seabed and the anchor tension is nonzero */
      fvec[0] = residual_function_length_contact(Fv, Fh, omega, Lu, EA, length, cb);
      fvec[1] = residual_function_height_contact(Fv, Fh, omega, Lu, EA, height, cb); 
    };
  } else {
    if (contactFlag==true || omega<0.0 || (Fv-omega*Lu)>0.0) { /* true when no portion of the line rests on the seabed */
      fjac[0] = jacobian_dxdh_no_contact(Fv, Fh, omega, Lu, EA);
      fjac[1] = jacobian_dxdv_no_contact(Fv, Fh, omega, Lu, EA);
      fjac[2] = jacobian_dzdh_no_contact(Fv, Fh, omega, Lu, EA);
      fjac[3] = jacobian_dzdv_no_contact(Fv, Fh, omega, Lu, EA);
    } else { /* true when a portion of the line rests on the seabed and the anchor tension is nonzero */
      fjac[0] = jacobian_dxdh_contact(Fv, Fh, omega, Lu, EA, cb);
      fjac[1] = jacobian_dxdv_contact(Fv, Fh, omega, Lu, EA, cb);
      fjac[2] = jacobian_dzdh_contact(Fv, Fh, omega, Lu, EA, cb);
      fjac[3] = jacobian_dzdv_contact(Fv, Fh, omega, Lu, EA, cb);
    };
  };
  return 0;
};


/**
 *
 */
MapReal residual_function_length_no_contact(const MapReal V, const MapReal H, const MapReal w, const MapReal Lu, const MapReal EA, const MapReal l)
{  
  return (H/w)*asinh(V/H) - (H/w)*asinh( (V-w*Lu)/H ) + ((H*Lu)/(EA)) - l;
};


/**
 *
 */
MapReal residual_function_height_no_contact(const MapReal V, const MapReal H, const MapReal w, const MapReal Lu, const MapReal EA, const MapReal h)
{
  return (H/w)* sqrt(1 + pow((V/H), 2)) - (H/w)*sqrt(1 + pow(((V-w*Lu)/H), 2)) + 1/(EA)*(V*Lu - (w*Lu*Lu)/2) - h;
};


/**
 *
 */
MapReal jacobian_dxdh_no_contact(const MapReal V, const MapReal H, const MapReal w, const MapReal Lu, const MapReal EA)
{
  return (asinh(V/H) - asinh((V-w*Lu)/H))/w - ((V/H + pow(V/H, 2)/sqrt(1.0 + pow(V/H, 2)))/(V/H + sqrt(1.0 + pow(V/H, 2))) 
                                               - ((V-w*Lu)/H + pow((V-w*Lu)/H, 2)/sqrt(1.0 + pow((V-w*Lu)/H, 2)))
                                               /((V-w*Lu)/H + sqrt(1.0 + pow((V-w*Lu)/H, 2))))/w + (Lu/(EA));      
};


/**
 *
 */
MapReal jacobian_dxdv_no_contact(const MapReal V, const MapReal H, const MapReal w, const MapReal Lu, const MapReal EA)
{
  return ((1.0 + V/H /sqrt(1.0 + pow(V/H, 2)))/(V/H + sqrt(1.0 + pow(V/H, 2))) 
          - (1.0 + (V-w*Lu)/H /sqrt(1.0 + pow( (V-w*Lu)/H , 2)))
          /((V-w*Lu)/H + sqrt(1.0 + pow((V-w*Lu)/H, 2))))/w;
};


/**
 *
 */
MapReal jacobian_dzdh_no_contact(const MapReal V, const MapReal H, const MapReal w, const MapReal Lu, const MapReal EA)
{
  return ( sqrt( 1.0 + pow( V/H , 2) ) - sqrt( 1.0 + pow( (V-w*Lu)/H , 2) ) )/w  
    - ( pow( V/H , 2 )/sqrt( 1.0 + pow( V/H , 2) ) - pow( (V-w*Lu)/H , 2)/sqrt( 1.0 + pow( (V-w*Lu)/H , 2) ) )/w;    
};



/**
 *
 */
MapReal jacobian_dzdv_no_contact(const MapReal V, const MapReal H, const MapReal w, const MapReal Lu, const MapReal EA)
{
  return ( V/H/sqrt( 1.0 + pow( V/H , 2) ) - (V-w*Lu)/H /sqrt( 1.0 + pow( (V-w*Lu)/H , 2) ) )/w + (Lu/(EA));
};


/**
 *
 */
MapReal residual_function_length_contact(const MapReal V, const MapReal H, const MapReal w, const MapReal Lu, const MapReal EA, const MapReal l, const MapReal cb)
{  
  if (-cb*(V-w*Lu)<H) { /* true when a portion of the line rests on the seabed and the anchor tension is nonzero */
    return log( (V/H) + sqrt(1.0 + pow(V/H,2)))*(H/w) - 0.5*(cb/EA)*w*(Lu-V/w)*(Lu-V/w) + (Lu/EA)*H + (Lu-V/w) - l;
  } else { /* 0.0<H<=-CB*(V-w*Lu), A  portion of the line must rest on the seabed and the anchor tension is zero */
    return log((V/H) + sqrt(1.0 + pow(V/H,2)))*(H/w) - 0.5*(cb/EA)*w*((Lu-V/w)*(Lu-V/w) - ((Lu-V/w) - (H/w)/cb)*((Lu-V/w) - (H/w)/cb)) + (Lu/EA)*H + (Lu-V/w) - l;
  };
};


/**
 *
 */
MapReal residual_function_height_contact(const MapReal V, const MapReal H, const MapReal w, const MapReal Lu, const MapReal EA, const MapReal h, const MapReal cb)
{
  // @todo remove this conditional statement because the equations are the same regardless of the outcome
  if (-cb*(V-w*Lu)<H) { /* true when a portion of the line rests on the seabed and the anchor tension is nonzero */
    return (sqrt(1.0 + pow(V/H,2)) - 1.0)*(H/w) + 0.5*V*(V/(w*EA)) - h;
  } else { /* 0.0<H<=-CB*(V-w*Lu), A  portion of the line must rest on the seabed and the anchor tension is zero */
    return (sqrt(1.0 + pow(V/H,2)) - 1.0)*(H/w) + 0.5*V*(V/(w*EA)) - h;
  };
};


/**
 *
 */
MapReal jacobian_dxdh_contact(const MapReal V, const MapReal H, const MapReal w, const MapReal Lu, const MapReal EA, const MapReal cb)
{
  if (-cb*(V-w*Lu)<H) { /* true when a portion of the line rests on the seabed and the anchor tension is nonzero */
    return log((V/H) + sqrt(1.0 + pow(V/H,2)))/w - (((V/H) + (V/H)*(V/H)/sqrt(1.0 + pow(V/H,2)))/((V/H) + sqrt(1.0 + pow(V/H,2))))/w + (Lu/EA);
  } else { /* 0.0<H<=-CB*(V-w*Lu), A  portion of the line must rest on the seabed and the anchor tension is zero */
    return log((V/H) + sqrt(1.0 + pow(V/H,2)))/w - (((V/H) + (V/H)*(V/H)/sqrt(1.0 + pow(V/H,2)))/((V/H) + sqrt(1.0 + pow(V/H,2))))/w + (Lu/EA) - ((Lu-V/w) - (H/w)/cb)/EA;
  };
};


/**
 *
 */
MapReal jacobian_dxdv_contact(const MapReal V, const MapReal H, const MapReal w, const MapReal Lu, const MapReal EA, const MapReal cb)
{
  if (-cb*(V-w*Lu)<H) { /* true when a portion of the line rests on the seabed and the anchor tension is nonzero */
    return ((1.0 + (V/H)/sqrt(1.0 + pow(V/H,2)))/((V/H) + sqrt(1.0 + pow(V/H,2))))/w + (cb/EA)*(Lu-V/w) - 1.0/w;
  } else { /* 0.0<H<=-CB*(V-w*Lu), A  portion of the line must rest on the seabed and the anchor tension is zero */
    return ((1.0 + (V/H)/sqrt(1.0 + pow(V/H,2)))/((V/H) + sqrt(1.0 + pow(V/H,2))))/w + (H/(w*EA)) - 1.0/w;
  };
};


/**
 *
 */
MapReal jacobian_dzdh_contact(const MapReal V, const MapReal H, const MapReal w, const MapReal Lu, const MapReal EA, const MapReal cb)
{
  // @todo remove this conditional statement because the equations are the same regardless of the outcome
  if (-cb*(V-w*Lu)<H) { /* true when a portion of the line rests on the seabed and the anchor tension is nonzero */
    return (sqrt(1.0 + pow(V/H,2)) - 1.0 - (V/H)*(V/H)/sqrt(1.0 + pow(V/H,2)))/w;
  } else { /* 0.0<H<=-CB*(V-w*Lu), A  portion of the line must rest on the seabed and the anchor tension is zero */
    return (sqrt(1.0 + pow(V/H,2)) - 1.0 - (V/H)*(V/H)/sqrt(1.0 + pow(V/H,2)))/w;
  };
};



/**
 *
 */
MapReal jacobian_dzdv_contact(const MapReal V, const MapReal H, const MapReal w, const MapReal Lu, const MapReal EA, const MapReal cb)
{
  // @todo remove this conditional statement because the equations are the same regardless of the outcome
  if (-cb*(V-w*Lu)<H) { /* true when a portion of the line rests on the seabed and the anchor tension is nonzero */
    return ((V/H)/sqrt(1.0 + pow(V/H,2)))/w + (V/(w*EA));
  } else { /* 0.0<H<=-CB*(V-w*Lu), A  portion of the line must rest on the seabed and the anchor tension is zero */
    return ((V/H)/sqrt(1.0 + pow(V/H,2)))/w + (V/(w*EA));
  };  
};


/**
 *
 */
MAP_EXTERNCALL MapReal pyget_residual_function_length(MAP_OtherStateType_t* otherType, int i, char* map_msg, MAP_ERROR_CODE* ierr)
{ 
  ModelData* data = otherType->object;
  Element* element = NULL;
  MapReal Fh = 0.0;
  MapReal Fv = 0.0;
  MapReal EA = 0.0;
  MapReal Lu = 0.0;
  MapReal length = 0.0;
  MapReal omega = 0.0;
  MapReal cb = 0.0;
  bool contactFlag = false;
  char buffer[64] = "";
  int cx = 0;
  
  *ierr = MAP_SAFE;
  map_msg[0] = 0;
  element = (Element*)list_get_at(&data->element, i);

  if (element==NULL) {
    cx = map_snprintf(buffer, 64, "Element out of range: %d.", i); assert(cx>=0);
    *ierr = map_set_universal_error(buffer, map_msg, ierr, MAP_FATAL_42);
    return -999.9;
  };

  Fh = *(element->H.value);
  Fv = *(element->V.value);  
  EA = element->lineProperty->ea;
  Lu = element->Lu.value;
  length = element->l.value;
  omega = element->lineProperty->omega;
  contactFlag = element->options.omitContact;
  cb = element->lineProperty->cb;

  if (contactFlag==true || omega<0.0 || (Fv-omega*Lu)>0.0) { /* true when no portion of the line rests on the seabed */
    return residual_function_length_no_contact(Fv, Fh, omega, Lu, EA, length);
  } else { /* true when a portion of the line rests on the seabed and the anchor tension is nonzero */
    return residual_function_length_contact(Fv, Fh, omega, Lu, EA, length, cb);
  };
};


/**
 *
 */
MAP_EXTERNCALL MapReal pyget_residual_function_height(MAP_OtherStateType_t* otherType, int i, char* map_msg, MAP_ERROR_CODE* ierr)
{
  ModelData* data = otherType->object;
  Element* element = NULL;
  MapReal Fh = 0.0;
  MapReal Fv = 0.0;
  MapReal EA = 0.0;
  MapReal Lu = 0.0;
  MapReal height = 0.0;
  MapReal omega = 0.0;
  MapReal cb = 0.0;
  bool contactFlag = false;
  char buffer[64] = "";
  int cx = 0;

  *ierr = MAP_SAFE;
  map_msg[0] = 0;
  element = (Element*)list_get_at( &data->element, i );

  if (element==NULL) {    
    cx = map_snprintf(buffer, 64, "Element out of range: %d.", i); assert(cx>=0);
    *ierr = map_set_universal_error( buffer, map_msg, ierr, MAP_FATAL_42);
    return -999.9;
  };

  Fh = *(element->H.value);
  Fv = *(element->V.value);  
  EA = element->lineProperty->ea;
  Lu = element->Lu.value;
  height = element->h.value;
  omega = element->lineProperty->omega;
  contactFlag = element->options.omitContact;
  cb = element->lineProperty->cb;

  if (contactFlag==true || omega<0.0 || (Fv-omega*Lu)>0.0) { /* true when no portion of the line rests on the seabed */
    return residual_function_height_no_contact(Fv, Fh, omega, Lu, EA, height);
  } else { /* true when a portion of the line rests on the seabed and the anchor tension is nonzero */
    return residual_function_height_contact(Fv, Fh, omega, Lu, EA, height, cb);
  };
};


/**
 *
 */
MAP_EXTERNCALL MapReal pyget_jacobian_dxdh(MAP_OtherStateType_t* otherType, int i, char* map_msg, MAP_ERROR_CODE* ierr)
{
  ModelData* data = otherType->object;
  Element* element = NULL;
  MapReal Fh = 0.0;
  MapReal Fv = 0.0;
  MapReal EA = 0.0;
  MapReal Lu = 0.0;
  MapReal omega = 0.0;
  MapReal cb = 0.0;
  bool contactFlag = false;
  char buffer[64] = "";
  int cx = 0;

  *ierr = MAP_SAFE;
  map_msg[0] = 0;
  element = (Element*)list_get_at(&data->element, i);

  if (element==NULL) {    
    cx = map_snprintf(buffer, 64, "Element out of range: %d.", i); assert(cx>=0);
    *ierr = map_set_universal_error(buffer, map_msg, ierr, MAP_FATAL_42);
    return -999.9;
  };

  Fh = *(element->H.value);
  Fv = *(element->V.value);  
  EA = element->lineProperty->ea;
  Lu = element->Lu.value;
  omega = element->lineProperty->omega;
  contactFlag = element->options.omitContact;
  cb = element->lineProperty->cb;

  if (contactFlag==true || omega<0.0 || (Fv-omega*Lu)>0.0) { /* true when no portion of the line rests on the seabed */
    return jacobian_dxdh_no_contact(Fv, Fh, omega, Lu, EA);
  } else { /* true when a portion of the line rests on the seabed and the anchor tension is nonzero */
    return jacobian_dxdh_contact(Fv, Fh, omega, Lu, EA, cb);
  };
};


/**
 *
 */
MAP_EXTERNCALL MapReal pyget_jacobian_dxdv(MAP_OtherStateType_t* otherType, int i, char* map_msg, MAP_ERROR_CODE* ierr)
{
  ModelData* data = otherType->object;
  Element* element = NULL;
  MapReal Fh = 0.0;
  MapReal Fv = 0.0;
  MapReal EA = 0.0;
  MapReal Lu = 0.0;
  MapReal omega = 0.0;
  MapReal cb = 0.0;
  bool contactFlag = false;
  char buffer[64] = "";
  int cx = 0;

  *ierr = MAP_SAFE;
  map_msg[0] = 0;
  element = (Element*)list_get_at(&data->element, i);

  if (element==NULL) {    
    cx = map_snprintf(buffer, 64, "Element out of range: %d.", i); assert(cx>=0);
    *ierr = map_set_universal_error(buffer, map_msg, ierr, MAP_FATAL_42);
    return -999.9;
  };

  Fh = *(element->H.value);
  Fv = *(element->V.value);  
  EA = element->lineProperty->ea;
  Lu = element->Lu.value;
  omega = element->lineProperty->omega;
  contactFlag = element->options.omitContact;
  cb = element->lineProperty->cb;

  if (contactFlag==true || omega<0.0 || (Fv-omega*Lu)>0.0) { /* true when no portion of the line rests on the seabed */
    return jacobian_dxdv_no_contact(Fv, Fh, omega, Lu, EA);
  } else { /* true when a portion of the line rests on the seabed and the anchor tension is nonzero */
    return jacobian_dxdv_contact(Fv, Fh, omega, Lu, EA, cb);
  };
};


/**
 *
 */
MAP_EXTERNCALL MapReal pyget_jacobian_dzdh(MAP_OtherStateType_t* otherFortType, int i, char* map_msg, MAP_ERROR_CODE* ierr)
{
  ModelData* data = otherFortType->object;
  Element* element = NULL;
  MapReal Fh = 0.0;
  MapReal Fv = 0.0;
  MapReal EA = 0.0;
  MapReal Lu = 0.0;
  MapReal omega = 0.0;
  MapReal cb = 0.0;
  bool contactFlag = false;
  char buffer[64] = "";
  int cx = 0;

  *ierr = MAP_SAFE;
  map_msg[0] = 0;
  element = (Element*)list_get_at(&data->element, i);

  if (element==NULL) {    
    cx = map_snprintf(buffer, 64, "Element out of range: %d.", i ); assert(cx>=0);
    *ierr = map_set_universal_error(buffer, map_msg, ierr, MAP_FATAL_42);
    return -999.9;
  };

  Fh = *(element->H.value);
  Fv = *(element->V.value);  
  EA = element->lineProperty->ea;
  Lu = element->Lu.value;
  omega = element->lineProperty->omega;
  contactFlag = element->options.omitContact;
  cb = element->lineProperty->cb;
  
  if (contactFlag==true || omega<0.0 || (Fv-omega*Lu)>0.0) { /* true when no portion of the line rests on the seabed */
    return jacobian_dzdh_no_contact(Fv, Fh, omega, Lu, EA);
  } else { /* true when a portion of the line rests on the seabed and the anchor tension is nonzero */
    return jacobian_dzdh_contact(Fv, Fh, omega, Lu, EA, cb);
  };
};


/**
 *
 */
MAP_EXTERNCALL MapReal pyget_jacobian_dzdv(MAP_OtherStateType_t* otherFortType, int i, char* map_msg, MAP_ERROR_CODE* ierr)
{
  ModelData* data = otherFortType->object;
  Element* element = NULL;
  MapReal Fh = 0.0;
  MapReal Fv = 0.0;
  MapReal EA = 0.0;
  MapReal Lu = 0.0;
  MapReal omega  = 0.0;
  MapReal cb = 0.0;
  char buffer[64] = "";
  bool contactFlag = false;
  int cx = 0;
 
  *ierr = MAP_SAFE;
  map_msg[0] = 0;
  element = (Element*)list_get_at(&data->element, i);

  if (element==NULL) {    
    cx = map_snprintf(buffer, 64, "Element out of range: %d.", i); assert(cx>=0);
    *ierr = map_set_universal_error(buffer, map_msg, ierr, MAP_FATAL_42);
    return -999.9;
  };

  Fh = *(element->H.value);
  Fv = *(element->V.value);  
  EA = element->lineProperty->ea;
  Lu = element->Lu.value;
  omega  = element->lineProperty->omega;
  contactFlag = element->options.omitContact;
  cb = element->lineProperty->cb;

  if (contactFlag==true || omega<0.0 || (Fv-omega*Lu)>0.0) { /* true when no portion of the line rests on the seabed */
    return jacobian_dzdv_no_contact(Fv, Fh, omega, Lu, EA);
  } else { /* true when a portion of the line rests on the seabed and the anchor tension is nonzero */
    return jacobian_dzdv_contact(Fv, Fh, omega, Lu, EA, cb);
  };
};


/**
 * check to see if the line is double backing
 ELSEIF ( W  >  0.0_DbKi )  THEN   ! .TRUE. when the line will sink in fluid

         LMax      = XF - EA/W + SQRT( (EA/W)*(EA/W) + 2.0_DbKi*ZF*EA/W )  ! Compute the maximum stretched length of the line with seabed interaction beyond which the line would have to double-back on itself; here the line forms an "L" between the anchor and fairlead (i.e. it is horizontal along the seabed from the anchor, then vertical to the fairlead)

         IF ( ( L  >=  LMax   ) .AND. ( CB >= 0.0_DbKi ) )  &  ! .TRUE. if the line is as long or longer than its maximum possible value with seabed interaction
            CALL ProgAbort ( ' Unstretched mooring line length too large. '// &
                         ' Routine Catenary() cannot solve quasi-static mooring line solution.' )


      ENDIF
 */
MAP_ERROR_CODE check_maximum_line_length(Element* element, const bool contactFlag, char *map_msg, MAP_ERROR_CODE *ierr)
{
  int cx = 0;
  char buffer[64]="";
  MapReal LMax = 0.0;
  MapReal l = element->l.value;
  MapReal h = element->h.value;
  MapReal EA = element->lineProperty->ea;
  MapReal W = element->lineProperty->omega;
  MapReal Lu = element->Lu.value;

  LMax = l - EA/W + sqrt(pow((EA/W),2) + 2.0*h*EA/W);
  if (Lu>=LMax && contactFlag==false) {
    cx = map_snprintf(buffer, 64, "LMax = %f meters.", LMax); assert(cx>=0);
    *ierr = map_set_universal_error(buffer, map_msg, ierr, MAP_FATAL_58);      
    return MAP_FATAL;
  };
  return MAP_SAFE;
};


/**
 *
 */
MapReal get_maximum_line_length(Element* element)
{
  MapReal l = element->l.value;
  MapReal h = element->h.value;
  MapReal EA = element->lineProperty->ea;
  MapReal w = element->lineProperty->omega;
  MapReal Lu = element->Lu.value;

  return (l - EA/w + sqrt(pow((EA/w),2) + 2.0*h*EA/w));
};


/**
 * Ax = b -> LUx = b. Then y is defined to be Ux
 */
MAP_ERROR_CODE lu_back_substitution(OuterSolveAttributes* ns, const int n, char* map_msg, MAP_ERROR_CODE* ierr)
{
  int i = 0;
  int j = 0;

  /* Forward solve Ly = b */
  for (i=0 ; i<n ; i++) {
    ns->y[i] = ns->b[i];
    for (j=0 ; j<i ; j++) {
      ns->y[i] -= (ns->l[i][j])*(ns->y[j]);
    };
    if (fabs(ns->l[i][i])<MACHINE_EPSILON) {
      return MAP_FATAL;
    };
    ns->y[i] /= ns->l[i][i];    
  };

  /* Backward solve Ux = y */
  for (i=n-1 ; i>=0 ; i--) {
    ns->x[i] = ns->y[i];
    for (j=i+1 ; j<n ; j++) {
      ns->x[i] -= (ns->u[i][j])*(ns->x[j]);
    };    
    if (fabs(ns->u[i][i])<MACHINE_EPSILON) {
      return MAP_FATAL;
    };
    ns->x[i] /= ns->u[i][i];
  };
  return MAP_SAFE;
};


/**
 *
 */
MAP_ERROR_CODE lu(OuterSolveAttributes* ns, const int n, char* map_msg, MAP_ERROR_CODE* ierr)
{
  int i = 0;
  int j = 0;
  int k = 0;

  for (i=0 ; i<n ; i++) {
    ns->l[i][i] = 1.0;
    for (j=i+1 ; j<n ; j++) {
      if (fabs(ns->jac[i][i])<MACHINE_EPSILON) {
        return MAP_FATAL;
      };
      ns->l[j][i] = (ns->jac[j][i])/(ns->jac[i][i]);
      /* ns->jac[j][j] = ns->l[j][j] */
      for (k=i+1 ; k<n ; k++) {
        ns->jac[j][k] = ns->jac[j][k] - (ns->l[j][i])*(ns->jac[i][k]);
      };
    };
    
    for (k=i ; k<n ; k++) {
      ns->u[i][k] = ns->jac[i][k] ;
    };
  };

  return MAP_SAFE;
};


/**
 * Forward difference jacobian
 */
MAP_ERROR_CODE forward_difference_jacobian(MAP_OtherStateType_t* otherType, MAP_ConstraintStateType_t* zType, ModelData* data, char* map_msg, MAP_ERROR_CODE* ierr)
{
  OuterSolveAttributes* ns = &data->outer_loop;
  MAP_ERROR_CODE success = MAP_SAFE;
  double originalDisplacement = 0.0;
  const int THREE = 3;
  const int N = zType->z_Len;
  const int m = THREE*(otherType->Fz_connect_Len); // rows
  const int n = THREE*(zType->z_Len);              // columns
  char buffer[64] = "";
  int cx = 0;
  int i = 0;
  int j = 0;
  
  for (i=0 ; i<N ; i++) {
    ns->b[THREE*i] = otherType->Fx_connect[i];
    ns->b[THREE*i+1] = otherType->Fy_connect[i];
    ns->b[THREE*i+2] = otherType->Fz_connect[i];      
  }

  /* 
     First store the connect node force before applying ns->epsilon displacement      
  */
  for (i=0 ; i<n ; i++) {            
    for (j=0 ; j<N ; j++) {            
      ns->jac[THREE*j][i] = -otherType->Fx_connect[j];
      ns->jac[THREE*j+1][i] = -otherType->Fy_connect[j];
      ns->jac[THREE*j+2][i] = -otherType->Fz_connect[j];
    };
  };
    
  for (j=0 ; j<N ; j++) {   // column               
    for (i=0 ; i<N ; i++) { // rows           
      originalDisplacement = zType->x[j];
      zType->x[j] += ns->epsilon;
      success = line_solve_sequence(data, 0.0, map_msg, ierr);
      if (success) {
        cx = map_snprintf(buffer, 64, "Forward difference, x[%d]+delta, row %d, col %d.", j+1, THREE*i, THREE*j); assert(cx>=0);
        *ierr = map_set_universal_error(buffer, map_msg, ierr, MAP_FATAL_78);      
        return MAP_FATAL;
      };
      ns->jac[THREE*i][THREE*j] += otherType->Fx_connect[i];
      ns->jac[THREE*i][THREE*j] /= ns->epsilon;
      ns->jac[THREE*i+1][THREE*j] += otherType->Fy_connect[i];
      ns->jac[THREE*i+1][THREE*j] /= ns->epsilon;
      ns->jac[THREE*i+2][THREE*j] += otherType->Fz_connect[i];
      ns->jac[THREE*i+2][THREE*j] /= ns->epsilon;
      zType->x[j] = originalDisplacement;
        
      originalDisplacement = zType->y[j];
      zType->y[j] += ns->epsilon;
      success = line_solve_sequence(data, 0.0, map_msg, ierr);
      if (success) {
        cx = map_snprintf(buffer, 64, "Forward difference, y[%d]+delta, row %d, col %d.", j+1, THREE*i, THREE*j); assert(cx>=0);
        *ierr = map_set_universal_error(buffer, map_msg, ierr, MAP_FATAL_78);      
        return MAP_FATAL;
      };
      ns->jac[THREE*i][THREE*j+1] += otherType->Fx_connect[i];
      ns->jac[THREE*i][THREE*j+1] /= ns->epsilon;
      ns->jac[THREE*i+1][THREE*j+1] += otherType->Fy_connect[i];
      ns->jac[THREE*i+1][THREE*j+1] /= ns->epsilon;
      ns->jac[THREE*i+2][THREE*j+1] += otherType->Fz_connect[i];
      ns->jac[THREE*i+2][THREE*j+1] /= ns->epsilon;
      zType->y[j] = originalDisplacement;
    
      originalDisplacement = zType->z[j];
      zType->z[j] += ns->epsilon;
      success = line_solve_sequence(data, 0.0, map_msg, ierr);
      if (success) {
        cx = map_snprintf(buffer, 64, "Forward difference, z[%d]+delta, row %d, col %d.", j+1, THREE*i, THREE*j); assert(cx>=0);
        *ierr = map_set_universal_error(buffer, map_msg, ierr, MAP_FATAL_78);      
        return MAP_FATAL;
      };
      ns->jac[THREE*i][THREE*j+2] += otherType->Fx_connect[i];
      ns->jac[THREE*i][THREE*j+2] /= ns->epsilon;
      ns->jac[THREE*i+1][THREE*j+2] += otherType->Fy_connect[i];
      ns->jac[THREE*i+1][THREE*j+2] /= ns->epsilon;
      ns->jac[THREE*i+2][THREE*j+2] += otherType->Fz_connect[i];
      ns->jac[THREE*i+2][THREE*j+2] /= ns->epsilon;
      zType->z[j] = originalDisplacement;
    };
  };

  /* read flag to set scaling parameter */
  if (ns->pg) {
    for (i=0 ; i<THREE*N ; i++) { 
      ns->jac[i][i] += (ns->ds/pow(ns->iterationCount,1.5)+ns->d);
    };
  };

  return MAP_SAFE;
};


/**
 * Backward difference jacobian
 */
MAP_ERROR_CODE backward_difference_jacobian(MAP_OtherStateType_t* otherType, MAP_ConstraintStateType_t* zType, ModelData* data, char* map_msg, MAP_ERROR_CODE* ierr)
{
  OuterSolveAttributes* ns = &data->outer_loop;
  MAP_ERROR_CODE success = MAP_SAFE;
  double originalDisplacement = 0.0;
  const int THREE = 3;
  const int N = zType->z_Len;
  const int m = THREE*(otherType->Fz_connect_Len); // rows
  const int n = THREE*(zType->z_Len);              // columns
  char buffer[64] = "";
  int cx = 0;
  int i = 0;
  int j = 0;
  
  for (i=0 ; i<N ; i++) {
    ns->b[THREE*i] = otherType->Fx_connect[i];
    ns->b[THREE*i+1] = otherType->Fy_connect[i];
    ns->b[THREE*i+2] = otherType->Fz_connect[i];      
  }

  /* 
     First store the connect node force before applying ns->epsilon displacement      
  */
  for (i=0 ; i<n ; i++) {            
    for (j=0 ; j<N ; j++) {            
      ns->jac[THREE*j][i] = otherType->Fx_connect[j];
      ns->jac[THREE*j+1][i] = otherType->Fy_connect[j];
      ns->jac[THREE*j+2][i] = otherType->Fz_connect[j];
    };
  };

  for (j=0 ; j<N ; j++) {   // column               
    for (i=0 ; i<N ; i++) { // rows           
      originalDisplacement = zType->x[j];
      zType->x[j] -= ns->epsilon;
      success = line_solve_sequence(data, 0.0, map_msg, ierr);
      if (success) {
        cx = map_snprintf(buffer, 64, "Backward difference, x[%d]-delta, row %d, col %d.", j+1, THREE*i, THREE*j); assert(cx>=0);
        *ierr = map_set_universal_error(buffer, map_msg, ierr, MAP_FATAL_78);      
        return MAP_FATAL;
      };
      ns->jac[THREE*i][THREE*j] -= otherType->Fx_connect[i];
      ns->jac[THREE*i][THREE*j] /= ns->epsilon;
      ns->jac[THREE*i+1][THREE*j] -= otherType->Fy_connect[i];
      ns->jac[THREE*i+1][THREE*j] /= ns->epsilon;
      ns->jac[THREE*i+2][THREE*j] -= otherType->Fz_connect[i];
      ns->jac[THREE*i+2][THREE*j] /= ns->epsilon;
      zType->x[j] = originalDisplacement;
        
      originalDisplacement = zType->y[j];
      zType->y[j] -= ns->epsilon;
      success = line_solve_sequence(data, 0.0, map_msg, ierr);
      if (success) {
        cx = map_snprintf(buffer, 64, "Backward difference, y[%d]-delta, row %d, col %d.", j+1, THREE*i, THREE*j); assert(cx>=0);
        *ierr = map_set_universal_error(buffer, map_msg, ierr, MAP_FATAL_78);      
        return MAP_FATAL;
      };
      ns->jac[THREE*i][THREE*j+1] -= otherType->Fx_connect[i];
      ns->jac[THREE*i][THREE*j+1] /= ns->epsilon;
      ns->jac[THREE*i+1][THREE*j+1] -= otherType->Fy_connect[i];
      ns->jac[THREE*i+1][THREE*j+1] /= ns->epsilon;
      ns->jac[THREE*i+2][THREE*j+1] -= otherType->Fz_connect[i];
      ns->jac[THREE*i+2][THREE*j+1] /= ns->epsilon;
      zType->y[j] = originalDisplacement;
    
      originalDisplacement = zType->z[j];
      zType->z[j] -= ns->epsilon;
      success = line_solve_sequence(data, 0.0, map_msg, ierr);
      if (success) {
        cx = map_snprintf(buffer, 64, "Backward difference, z[%d]-delta, row %d, col %d.", j+1, THREE*i, THREE*j); assert(cx>=0);
        *ierr = map_set_universal_error(buffer, map_msg, ierr, MAP_FATAL_78);      
        return MAP_FATAL;
      };
      ns->jac[THREE*i][THREE*j+2] -= otherType->Fx_connect[i];
      ns->jac[THREE*i][THREE*j+2] /= ns->epsilon;
      ns->jac[THREE*i+1][THREE*j+2] -= otherType->Fy_connect[i];
      ns->jac[THREE*i+1][THREE*j+2] /= ns->epsilon;
      ns->jac[THREE*i+2][THREE*j+2] -= otherType->Fz_connect[i];
      ns->jac[THREE*i+2][THREE*j+2] /= ns->epsilon;
      zType->z[j] = originalDisplacement;
    };
  };

  /* read flag to set scaling parameter */
  if (ns->pg) {
    for (i=0 ; i<THREE*N ; i++) { 
      ns->jac[i][i] += (ns->ds/pow(ns->iterationCount,1.5)+ns->d);
    };
  };

  return MAP_SAFE;
};


/**
 * central difference jacobian
 */
MAP_ERROR_CODE central_difference_jacobian(MAP_OtherStateType_t* otherType, MAP_ConstraintStateType_t* zType, ModelData* data, char* map_msg, MAP_ERROR_CODE* ierr)
{
  OuterSolveAttributes* ns = &data->outer_loop;
  MAP_ERROR_CODE success = MAP_SAFE;
  double originalDisplacement = 0.0;
  const int THREE = 3;
  const int N = zType->z_Len;
  const int m = THREE*(otherType->Fz_connect_Len); // rows
  const int n = THREE*(zType->z_Len);              // columns
  char buffer[64] = "";
  int cx = 0;
  int i = 0;
  int j = 0;
  
  for (i=0 ; i<N ; i++) {
    ns->b[THREE*i] = otherType->Fx_connect[i];
    ns->b[THREE*i+1] = otherType->Fy_connect[i];
    ns->b[THREE*i+2] = otherType->Fz_connect[i];      
  }

  /* 
     First store the connect node force before applying ns->epsilon displacement      
  */
  for (j=0 ; j<N ; j++) {   // column               
    for (i=0 ; i<N ; i++) { // rows           
      originalDisplacement = zType->x[j];
      zType->x[j] += ns->epsilon;
      success = line_solve_sequence(data, 0.0, map_msg, ierr);
      if (success) {
        cx = map_snprintf(buffer, 64, "Central difference, x[%d]+delta, row %d, col %d.", j+1, THREE*i, THREE*j); assert(cx>=0);
        *ierr = map_set_universal_error(buffer, map_msg, ierr, MAP_FATAL_78);      
        return MAP_FATAL;
      };
      ns->jac[THREE*i][THREE*j] = otherType->Fx_connect[i];      
      ns->jac[THREE*i+1][THREE*j] = otherType->Fy_connect[i];
      ns->jac[THREE*i+2][THREE*j] = otherType->Fz_connect[i];
      zType->x[j] = originalDisplacement;
        
      originalDisplacement = zType->y[j];
      zType->y[j] += ns->epsilon;
      if (success) {
        cx = map_snprintf(buffer, 64, "Central difference, y[%d]+delta, row %d, col %d.", j+1, THREE*i, THREE*j); assert(cx>=0);
        *ierr = map_set_universal_error(buffer, map_msg, ierr, MAP_FATAL_78);      
        return MAP_FATAL;
      };
      success = line_solve_sequence(data, 0.0, map_msg, ierr);
      ns->jac[THREE*i][THREE*j+1] = otherType->Fx_connect[i];
      ns->jac[THREE*i+1][THREE*j+1] = otherType->Fy_connect[i];
      ns->jac[THREE*i+2][THREE*j+1] = otherType->Fz_connect[i];
      zType->y[j] = originalDisplacement;
    
      originalDisplacement = zType->z[j];
      zType->z[j] += ns->epsilon;
      if (success) {
        cx = map_snprintf(buffer, 64, "Central difference, z[%d]+delta, row %d, col %d.", j+1, THREE*i, THREE*j); assert(cx>=0);
        *ierr = map_set_universal_error(buffer, map_msg, ierr, MAP_FATAL_78);      
        return MAP_FATAL;
      };
      success = line_solve_sequence(data, 0.0, map_msg, ierr);
      ns->jac[THREE*i][THREE*j+2] = otherType->Fx_connect[i];
      ns->jac[THREE*i+1][THREE*j+2] = otherType->Fy_connect[i];
      ns->jac[THREE*i+2][THREE*j+2] = otherType->Fz_connect[i];
      zType->z[j] = originalDisplacement;
    };
  };
    
  for (j=0 ; j<N ; j++) {   // column               
    for (i=0 ; i<N ; i++) { // rows           
      originalDisplacement = zType->x[j];
      zType->x[j] -= ns->epsilon;
      if (success) {
        cx = map_snprintf(buffer, 64, "Central difference, x[%d]-delta, row %d, col %d.", j+1, THREE*i, THREE*j); assert(cx>=0);
        *ierr = map_set_universal_error(buffer, map_msg, ierr, MAP_FATAL_78);      
        return MAP_FATAL;
      };
      success = line_solve_sequence(data, 0.0, map_msg, ierr);
      ns->jac[THREE*i][THREE*j] -= otherType->Fx_connect[i];
      ns->jac[THREE*i][THREE*j] /= (2*ns->epsilon);
      ns->jac[THREE*i+1][THREE*j] -= otherType->Fy_connect[i];
      ns->jac[THREE*i+1][THREE*j] /= (2*ns->epsilon);
      ns->jac[THREE*i+2][THREE*j] -= otherType->Fz_connect[i];
      ns->jac[THREE*i+2][THREE*j] /= (2*ns->epsilon);
      zType->x[j] = originalDisplacement;
        
      originalDisplacement = zType->y[j];
      zType->y[j] -= ns->epsilon;
      success = line_solve_sequence(data, 0.0, map_msg, ierr);
      if (success) {
        cx = map_snprintf(buffer, 64, "Central difference, y[%d]-delta, row %d, col %d.", j+1, THREE*i, THREE*j); assert(cx>=0);
        *ierr = map_set_universal_error(buffer, map_msg, ierr, MAP_FATAL_78);      
        return MAP_FATAL;
      };
      ns->jac[THREE*i][THREE*j+1] -= otherType->Fx_connect[i];
      ns->jac[THREE*i][THREE*j+1] /= (2*ns->epsilon);
      ns->jac[THREE*i+1][THREE*j+1] -= otherType->Fy_connect[i];
      ns->jac[THREE*i+1][THREE*j+1] /= (2*ns->epsilon);
      ns->jac[THREE*i+2][THREE*j+1] -= otherType->Fz_connect[i];
      ns->jac[THREE*i+2][THREE*j+1] /= (2*ns->epsilon);
      zType->y[j] = originalDisplacement;
    
      originalDisplacement = zType->z[j];
      zType->z[j] -= ns->epsilon;
      success = line_solve_sequence(data, 0.0, map_msg, ierr);
      if (success) {
        cx = map_snprintf(buffer, 64, "Central difference, z[%d]-delta, row %d, col %d.", j+1, THREE*i, THREE*j); assert(cx>=0);
        *ierr = map_set_universal_error(buffer, map_msg, ierr, MAP_FATAL_78);      
        return MAP_FATAL;
      };
      ns->jac[THREE*i][THREE*j+2] -= otherType->Fx_connect[i];
      ns->jac[THREE*i][THREE*j+2] /= (2*ns->epsilon);
      ns->jac[THREE*i+1][THREE*j+2] -= otherType->Fy_connect[i];
      ns->jac[THREE*i+1][THREE*j+2] /= (2*ns->epsilon);
      ns->jac[THREE*i+2][THREE*j+2] -= otherType->Fz_connect[i];
      ns->jac[THREE*i+2][THREE*j+2] /= (2*ns->epsilon);
      zType->z[j] = originalDisplacement;
    };
  };

  /* read flag to set scaling parameter */
  if (ns->pg) {
    ns->coef = pow(ns->iterationCount,1.5);
    for (i=0 ; i<THREE*N ; i++) { 
      ns->jac[i][i] += (ns->ds/ns->coef + ns->d);
    };
  };

  return MAP_SAFE;
};



MAP_ERROR_CODE call_minpack_lmder(Element* element, InnerSolveAttributes* inner_opt, ModelOptions* opt, const int lineNum, const double time, char* map_msg, MAP_ERROR_CODE* ierr)
{
  char buffer[64] = "";
  int cx = 0;
  MAP_ERROR_CODE success = MAP_SAFE;

  /* initial guess vector is set in set_element_initial_guess(..); otherwise, the previous solution is used as the initial guess */
  inner_opt->x[0] = *(element->H.value);
  inner_opt->x[1] = *(element->V.value);

  element->numFuncEvals = 0;
  element->numJacEvals = 0;

  inner_opt->info = __cminpack_func__(lmder)(inner_function_evals, 
                                      element, 
                                      inner_opt->m, 
                                      inner_opt->n, 
                                      inner_opt->x, 
                                      inner_opt->fvec, 
                                      inner_opt->fjac, 
                                      inner_opt->ldfjac, 
                                      inner_opt->f_tol, 
                                      inner_opt->x_tol, 
                                      inner_opt->g_tol, 
                                      inner_opt->max_its, 
                                      inner_opt->diag,
                                      inner_opt->mode, 
                                      inner_opt->factor, 
                                      inner_opt->nprint, 
                                      &element->numFuncEvals, 
                                      &element->numJacEvals, 
                                      inner_opt->ipvt, 
                                      inner_opt->qtf, 
                                      inner_opt->wa1 ,
                                      inner_opt->wa2 ,
                                      inner_opt->wa3 , 
                                      inner_opt->wa4);
  
  element->residualNorm = (MapReal)__minpack_func__(enorm)(&inner_opt->m, inner_opt->fvec);
  
  if (element->options.diagnosticsFlag && (double)element->diagnosticType>time ) { 
    printf("\n      %4.3f [sec]  Element %d\n",time, lineNum+1);
    printf("      ----------------------------------------------------\n");
    printf("      Residual l2 norm at solution:  %15.7g\n", element->residualNorm);
    printf("      Function evaluations:         %10i\n", element->numFuncEvals);
    printf("      Jacobian evaluations:         %10i\n", element->numJacEvals);
    printf("      Exit parameter                %10i\n\n", inner_opt->info);
  };
  *(element->H.value) = inner_opt->x[0];
  *(element->V.value) = inner_opt->x[1];
  element->convergeReason = inner_opt->info;
  
  switch (inner_opt->info) {
  case 0 :
    success = MAP_FATAL;
    cx = map_snprintf(buffer, 64, "Line segment %d.", lineNum); assert(cx>=0);
    *ierr = map_set_universal_error(buffer, map_msg, ierr, MAP_FATAL_39);
    break;
  case 1 :
    success = MAP_SAFE;
    break;
  case 2 :
    success = MAP_SAFE;
    break;
  case 3 :
    success = MAP_SAFE;
    break;
  case 4 :
    success = MAP_SAFE;
    break;
  case 5 :
    success = MAP_FATAL;
    cx = map_snprintf(buffer, 64, "Line segment %d.", lineNum); assert(cx>=0);
    *ierr = map_set_universal_error(buffer, map_msg, ierr, MAP_FATAL_40);
    break;
  case 6 :
    success = MAP_FATAL;
    cx = map_snprintf(buffer, 64, "Line segment %d.", lineNum); assert(cx>=0);
    *ierr = map_set_universal_error(buffer, map_msg, ierr, MAP_ERROR_11);
    break;
  case 7 :
    success = MAP_FATAL;
    cx = map_snprintf(buffer, 64, "Line segment %d.", lineNum); assert(cx>=0);
    *ierr = map_set_universal_error(buffer, map_msg, ierr, MAP_ERROR_13);
    break;
  case 8 :
    success = MAP_FATAL;
    cx = map_snprintf(buffer, 64, "Line segment %d.", lineNum); assert(cx>=0);
    *ierr = map_set_universal_error(buffer, map_msg, ierr, MAP_ERROR_12);
    break;
  default :
    success = MAP_SAFE;
    break;
  };
  return MAP_SAFE;
};
