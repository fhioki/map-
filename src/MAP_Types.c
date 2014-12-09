//!STARTOFREGISTRYGENERATEDFILE './MAP_Types.c'
//!
//! WARNING This file is generated automatically by the FAST registry
//! Do not edit.  Your changes to this file will be lost.
//!
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "MAP_Types.h"


#ifdef _MSC_VER //define something for Windows MS compiler
#  include "stdbool.h"
#  define CALL __declspec( dllexport )
#else
#  include <stdbool.h>
#  define CALL 
#endif


//#define CALL __attribute__((dllexport) )
//MAP_InitInputType_t* CALL MAP_InitInput_Create() { return ((MAP_InitInputType_t*) malloc( sizeof(MAP_InitInputType_t()))) ; } ;
//void CALL MAP_InitInput_Delete(MAP_InitInputType_t *This) { free(This) ; } ;
void CALL MAP_InitOutput_Create() { } ;
void CALL MAP_InitOutput_Delete(void* none) { } ;
void CALL MAP_ContState_Create() { } ;
void CALL MAP_ContState_Delete(void* none) { } ;
void CALL MAP_DiscState_Create() { } ;
void CALL MAP_DiscState_Delete(void* none) { } ;
//MAP_OtherStateType_t* CALL MAP_OtherState_Create() { return ((MAP_OtherStateType_t*) malloc( sizeof(MAP_OtherStateType_t()))) ; } ;
//void CALL MAP_OtherState_Delete(MAP_OtherStateType_t *This) { free(This) ; } ;
void CALL MAP_ConstrState_Create() { } ;
void CALL MAP_ConstrState_Delete(void* none) { } ;
void CALL MAP_Param_Create() { } ;
void CALL MAP_Param_Delete(void* none) { } ;
void CALL MAP_Input_Create() { } ;
void CALL MAP_Input_Delete(void* none) { } ;
void CALL MAP_Output_Create() { } ;
void CALL MAP_Output_Delete(void* none) { } ;

// int
// C_MAP_PackInitInput( float * ReKiBuf,  int * Re_BufSz ,
//                  double * DbKiBuf, int * Db_BufSz ,
//                  int * IntKiBuf,   int * Int_BufSz ,
//                  MAP_InitInputType_t *InData, char * ErrMsg, int *SizeOnly )
// {
//   int ErrStat ;
//   int OnlySize ;
//   int Re_BufSz2 ;
//   int Db_BufSz2 ;
//   int Int_BufSz2 ;
//   int Re_Xferred = 0 ;
//   int Db_Xferred = 0 ;
//   int Int_Xferred = 0 ;
//   int one         = 1 ;
//   int i,i1,i2,i3,i4,i5 ;
//  // buffers to store meshes and subtypes, if any
// 
//   OnlySize = *SizeOnly ;
// 
//   *Re_BufSz = 0 ;
//   *Db_BufSz = 0 ;
//   *Int_BufSz = 0 ;
//   ReKiBuf = NULL ;
//   DbKiBuf = NULL ;
//   IntKiBuf = NULL ;
//   *Db_BufSz   += 1  ; // gravity
//   *Db_BufSz   += 1  ; // sea_density
//   *Db_BufSz   += 1  ; // depth
//   if ( ! OnlySize ) {
//     if ( *Re_BufSz > 0  ) ReKiBuf  = (float  *)malloc(*Re_BufSz*sizeof(float) ) ;
//     if ( *Db_BufSz > 0  ) DbKiBuf  = (double *)malloc(*Db_BufSz*sizeof(double) ) ;
//     if ( *Int_BufSz > 0 ) IntKiBuf = (int    *)malloc(*Int_BufSz*sizeof(int) ) ;
//     DbKiBuf[Db_Xferred++] = InData->gravity ;
//     DbKiBuf[Db_Xferred++] = InData->sea_density ;
//     DbKiBuf[Db_Xferred++] = InData->depth ;
//   }
//   return(ErrStat) ;
// }
// 
// int
// C_MAP_UnpackInitInput( float * ReKiBuf,  
//                  double * DbKiBuf, 
//                  int * IntKiBuf,   
//                  MAP_InitInputType_t *OutData, char * ErrMsg )
// {
//   int ErrStat ;
//   int Re_BufSz2 = 0 ;
//   int Db_BufSz2 = 0 ;
//   int Int_BufSz2 = 0 ;
//   int Re_Xferred = 0 ;
//   int Db_Xferred = 0 ;
//   int Int_Xferred = 0 ;
//   int Re_CurrSz = 0 ;
//   int Db_CurrSz = 0 ;
//   int Int_CurrSz = 0 ;
//   int one        = 1 ;
//   int i,i1,i2,i3,i4,i5 ;
//  // buffers to store meshes, if any
//   ReKiBuf = NULL ;
//   DbKiBuf = NULL ;
//   IntKiBuf = NULL ;
//   OutData->gravity = DbKiBuf [ Db_Xferred ] ; 
//   Db_Xferred   = Db_Xferred   + 1 ; 
//   OutData->sea_density = DbKiBuf [ Db_Xferred ] ; 
//   Db_Xferred   = Db_Xferred   + 1 ; 
//   OutData->depth = DbKiBuf [ Db_Xferred ] ; 
//   Db_Xferred   = Db_Xferred   + 1 ; 
//   if ( ReKiBuf != NULL )  free(ReKiBuf) ;
//   if ( DbKiBuf != NULL )  free(DbKiBuf) ;
//   if ( IntKiBuf != NULL ) free(IntKiBuf) ;
//   return(ErrStat) ;
// }
// 
// int
// C_MAP_PackInitOutput( float * ReKiBuf,  int * Re_BufSz ,
//                  double * DbKiBuf, int * Db_BufSz ,
//                  int * IntKiBuf,   int * Int_BufSz ,
//                  MAP_InitOutputType_t *InData, char * ErrMsg, int *SizeOnly )
// {
//   int ErrStat ;
//   int OnlySize ;
//   int Re_BufSz2 ;
//   int Db_BufSz2 ;
//   int Int_BufSz2 ;
//   int Re_Xferred = 0 ;
//   int Db_Xferred = 0 ;
//   int Int_Xferred = 0 ;
//   int one         = 1 ;
//   int i,i1,i2,i3,i4,i5 ;
//  // buffers to store meshes and subtypes, if any
// 
//   OnlySize = *SizeOnly ;
// 
//   *Re_BufSz = 0 ;
//   *Db_BufSz = 0 ;
//   *Int_BufSz = 0 ;
//   ReKiBuf = NULL ;
//   DbKiBuf = NULL ;
//   IntKiBuf = NULL ;
//   if ( ! OnlySize ) {
//     if ( *Re_BufSz > 0  ) ReKiBuf  = (float  *)malloc(*Re_BufSz*sizeof(float) ) ;
//     if ( *Db_BufSz > 0  ) DbKiBuf  = (double *)malloc(*Db_BufSz*sizeof(double) ) ;
//     if ( *Int_BufSz > 0 ) IntKiBuf = (int    *)malloc(*Int_BufSz*sizeof(int) ) ;
//   }
//   return(ErrStat) ;
// }
// 
// int
// C_MAP_UnpackInitOutput( float * ReKiBuf,  
//                  double * DbKiBuf, 
//                  int * IntKiBuf,   
//                  MAP_InitOutputType_t *OutData, char * ErrMsg )
// {
//   int ErrStat ;
//   int Re_BufSz2 = 0 ;
//   int Db_BufSz2 = 0 ;
//   int Int_BufSz2 = 0 ;
//   int Re_Xferred = 0 ;
//   int Db_Xferred = 0 ;
//   int Int_Xferred = 0 ;
//   int Re_CurrSz = 0 ;
//   int Db_CurrSz = 0 ;
//   int Int_CurrSz = 0 ;
//   int one        = 1 ;
//   int i,i1,i2,i3,i4,i5 ;
//  // buffers to store meshes, if any
//   ReKiBuf = NULL ;
//   DbKiBuf = NULL ;
//   IntKiBuf = NULL ;
//   if ( ReKiBuf != NULL )  free(ReKiBuf) ;
//   if ( DbKiBuf != NULL )  free(DbKiBuf) ;
//   if ( IntKiBuf != NULL ) free(IntKiBuf) ;
//   return(ErrStat) ;
// }
// 
// CALL void MAP_F2C_InitOutput_writeOutputHdr_C ( MAP_InitOutputType_t *type, char *arr, int len )
// {
//   int i = 0;
//   for( i=0 ; i<=len-1 ; i++ ) type->writeOutputHdr[i] = arr[i];
// }
// 
// CALL void MAP_F2C_InitOutput_writeOutputUnt_C ( MAP_InitOutputType_t *type, char *arr, int len )
// {
//   int i = 0;
//   for( i=0 ; i<=len-1 ; i++ ) type->writeOutputUnt[i] = arr[i];
// }
// 
// int
// C_MAP_PackContState( float * ReKiBuf,  int * Re_BufSz ,
//                  double * DbKiBuf, int * Db_BufSz ,
//                  int * IntKiBuf,   int * Int_BufSz ,
//                  MAP_ContinuousStateType_t *InData, char * ErrMsg, int *SizeOnly )
// {
//   int ErrStat ;
//   int OnlySize ;
//   int Re_BufSz2 ;
//   int Db_BufSz2 ;
//   int Int_BufSz2 ;
//   int Re_Xferred = 0 ;
//   int Db_Xferred = 0 ;
//   int Int_Xferred = 0 ;
//   int one         = 1 ;
//   int i,i1,i2,i3,i4,i5 ;
//  // buffers to store meshes and subtypes, if any
// 
//   OnlySize = *SizeOnly ;
// 
//   *Re_BufSz = 0 ;
//   *Db_BufSz = 0 ;
//   *Int_BufSz = 0 ;
//   ReKiBuf = NULL ;
//   DbKiBuf = NULL ;
//   IntKiBuf = NULL ;
//   *Db_BufSz   += 1  ; // dummy
//   if ( ! OnlySize ) {
//     if ( *Re_BufSz > 0  ) ReKiBuf  = (float  *)malloc(*Re_BufSz*sizeof(float) ) ;
//     if ( *Db_BufSz > 0  ) DbKiBuf  = (double *)malloc(*Db_BufSz*sizeof(double) ) ;
//     if ( *Int_BufSz > 0 ) IntKiBuf = (int    *)malloc(*Int_BufSz*sizeof(int) ) ;
//     DbKiBuf[Db_Xferred++] = InData->dummy ;
//   }
//   return(ErrStat) ;
// }
// 
// int
// C_MAP_UnpackContState( float * ReKiBuf,  
//                  double * DbKiBuf, 
//                  int * IntKiBuf,   
//                  MAP_ContinuousStateType_t *OutData, char * ErrMsg )
// {
//   int ErrStat ;
//   int Re_BufSz2 = 0 ;
//   int Db_BufSz2 = 0 ;
//   int Int_BufSz2 = 0 ;
//   int Re_Xferred = 0 ;
//   int Db_Xferred = 0 ;
//   int Int_Xferred = 0 ;
//   int Re_CurrSz = 0 ;
//   int Db_CurrSz = 0 ;
//   int Int_CurrSz = 0 ;
//   int one        = 1 ;
//   int i,i1,i2,i3,i4,i5 ;
//  // buffers to store meshes, if any
//   ReKiBuf = NULL ;
//   DbKiBuf = NULL ;
//   IntKiBuf = NULL ;
//   OutData->dummy = DbKiBuf [ Db_Xferred ] ; 
//   Db_Xferred   = Db_Xferred   + 1 ; 
//   if ( ReKiBuf != NULL )  free(ReKiBuf) ;
//   if ( DbKiBuf != NULL )  free(DbKiBuf) ;
//   if ( IntKiBuf != NULL ) free(IntKiBuf) ;
//   return(ErrStat) ;
// }
// 
// int
// C_MAP_PackDiscState( float * ReKiBuf,  int * Re_BufSz ,
//                  double * DbKiBuf, int * Db_BufSz ,
//                  int * IntKiBuf,   int * Int_BufSz ,
//                  MAP_DiscreteStateType_t *InData, char * ErrMsg, int *SizeOnly )
// {
//   int ErrStat ;
//   int OnlySize ;
//   int Re_BufSz2 ;
//   int Db_BufSz2 ;
//   int Int_BufSz2 ;
//   int Re_Xferred = 0 ;
//   int Db_Xferred = 0 ;
//   int Int_Xferred = 0 ;
//   int one         = 1 ;
//   int i,i1,i2,i3,i4,i5 ;
//  // buffers to store meshes and subtypes, if any
// 
//   OnlySize = *SizeOnly ;
// 
//   *Re_BufSz = 0 ;
//   *Db_BufSz = 0 ;
//   *Int_BufSz = 0 ;
//   ReKiBuf = NULL ;
//   DbKiBuf = NULL ;
//   IntKiBuf = NULL ;
//   *Db_BufSz   += 1  ; // dummy
//   if ( ! OnlySize ) {
//     if ( *Re_BufSz > 0  ) ReKiBuf  = (float  *)malloc(*Re_BufSz*sizeof(float) ) ;
//     if ( *Db_BufSz > 0  ) DbKiBuf  = (double *)malloc(*Db_BufSz*sizeof(double) ) ;
//     if ( *Int_BufSz > 0 ) IntKiBuf = (int    *)malloc(*Int_BufSz*sizeof(int) ) ;
//     DbKiBuf[Db_Xferred++] = InData->dummy ;
//   }
//   return(ErrStat) ;
// }
// 
// int
// C_MAP_UnpackDiscState( float * ReKiBuf,  
//                  double * DbKiBuf, 
//                  int * IntKiBuf,   
//                  MAP_DiscreteStateType_t *OutData, char * ErrMsg )
// {
//   int ErrStat ;
//   int Re_BufSz2 = 0 ;
//   int Db_BufSz2 = 0 ;
//   int Int_BufSz2 = 0 ;
//   int Re_Xferred = 0 ;
//   int Db_Xferred = 0 ;
//   int Int_Xferred = 0 ;
//   int Re_CurrSz = 0 ;
//   int Db_CurrSz = 0 ;
//   int Int_CurrSz = 0 ;
//   int one        = 1 ;
//   int i,i1,i2,i3,i4,i5 ;
//  // buffers to store meshes, if any
//   ReKiBuf = NULL ;
//   DbKiBuf = NULL ;
//   IntKiBuf = NULL ;
//   OutData->dummy = DbKiBuf [ Db_Xferred ] ; 
//   Db_Xferred   = Db_Xferred   + 1 ; 
//   if ( ReKiBuf != NULL )  free(ReKiBuf) ;
//   if ( DbKiBuf != NULL )  free(DbKiBuf) ;
//   if ( IntKiBuf != NULL ) free(IntKiBuf) ;
//   return(ErrStat) ;
// }
// 
// int
// C_MAP_PackOtherState( float * ReKiBuf,  int * Re_BufSz ,
//                  double * DbKiBuf, int * Db_BufSz ,
//                  int * IntKiBuf,   int * Int_BufSz ,
//                  MAP_OtherStateType_t *InData, char * ErrMsg, int *SizeOnly )
// {
//   int ErrStat ;
//   int OnlySize ;
//   int Re_BufSz2 ;
//   int Db_BufSz2 ;
//   int Int_BufSz2 ;
//   int Re_Xferred = 0 ;
//   int Db_Xferred = 0 ;
//   int Int_Xferred = 0 ;
//   int one         = 1 ;
//   int i,i1,i2,i3,i4,i5 ;
//  // buffers to store meshes and subtypes, if any
// 
//   OnlySize = *SizeOnly ;
// 
//   *Re_BufSz = 0 ;
//   *Db_BufSz = 0 ;
//   *Int_BufSz = 0 ;
//   ReKiBuf = NULL ;
//   DbKiBuf = NULL ;
//   IntKiBuf = NULL ;
//   *Db_BufSz   += InData->H_Len ; // H 
//   *Db_BufSz   += InData->V_Len ; // V 
//   *Db_BufSz   += InData->Ha_Len ; // Ha 
//   *Db_BufSz   += InData->Va_Len ; // Va 
//   *Db_BufSz   += InData->x_Len ; // x 
//   *Db_BufSz   += InData->y_Len ; // y 
//   *Db_BufSz   += InData->z_Len ; // z 
//   *Db_BufSz   += InData->xa_Len ; // xa 
//   *Db_BufSz   += InData->ya_Len ; // ya 
//   *Db_BufSz   += InData->za_Len ; // za 
//   *Db_BufSz   += InData->Fx_connect_Len ; // Fx_connect 
//   *Db_BufSz   += InData->Fy_connect_Len ; // Fy_connect 
//   *Db_BufSz   += InData->Fz_connect_Len ; // Fz_connect 
//   *Db_BufSz   += InData->Fx_anchor_Len ; // Fx_anchor 
//   *Db_BufSz   += InData->Fy_anchor_Len ; // Fy_anchor 
//   *Db_BufSz   += InData->Fz_anchor_Len ; // Fz_anchor 
//   if ( ! OnlySize ) {
//     if ( *Re_BufSz > 0  ) ReKiBuf  = (float  *)malloc(*Re_BufSz*sizeof(float) ) ;
//     if ( *Db_BufSz > 0  ) DbKiBuf  = (double *)malloc(*Db_BufSz*sizeof(double) ) ;
//     if ( *Int_BufSz > 0 ) IntKiBuf = (int    *)malloc(*Int_BufSz*sizeof(int) ) ;
//     for ( i = 0 ; i < InData->H_Len ; i++ ) {
//       if ( !OnlySize ) memcpy( &(DbKiBuf[Db_Xferred+i]), &(InData->H[i]), sizeof(double)) ;
//       Db_Xferred++ ;
//     }
//     for ( i = 0 ; i < InData->V_Len ; i++ ) {
//       if ( !OnlySize ) memcpy( &(DbKiBuf[Db_Xferred+i]), &(InData->V[i]), sizeof(double)) ;
//       Db_Xferred++ ;
//     }
//     for ( i = 0 ; i < InData->Ha_Len ; i++ ) {
//       if ( !OnlySize ) memcpy( &(DbKiBuf[Db_Xferred+i]), &(InData->Ha[i]), sizeof(double)) ;
//       Db_Xferred++ ;
//     }
//     for ( i = 0 ; i < InData->Va_Len ; i++ ) {
//       if ( !OnlySize ) memcpy( &(DbKiBuf[Db_Xferred+i]), &(InData->Va[i]), sizeof(double)) ;
//       Db_Xferred++ ;
//     }
//     for ( i = 0 ; i < InData->x_Len ; i++ ) {
//       if ( !OnlySize ) memcpy( &(DbKiBuf[Db_Xferred+i]), &(InData->x[i]), sizeof(double)) ;
//       Db_Xferred++ ;
//     }
//     for ( i = 0 ; i < InData->y_Len ; i++ ) {
//       if ( !OnlySize ) memcpy( &(DbKiBuf[Db_Xferred+i]), &(InData->y[i]), sizeof(double)) ;
//       Db_Xferred++ ;
//     }
//     for ( i = 0 ; i < InData->z_Len ; i++ ) {
//       if ( !OnlySize ) memcpy( &(DbKiBuf[Db_Xferred+i]), &(InData->z[i]), sizeof(double)) ;
//       Db_Xferred++ ;
//     }
//     for ( i = 0 ; i < InData->xa_Len ; i++ ) {
//       if ( !OnlySize ) memcpy( &(DbKiBuf[Db_Xferred+i]), &(InData->xa[i]), sizeof(double)) ;
//       Db_Xferred++ ;
//     }
//     for ( i = 0 ; i < InData->ya_Len ; i++ ) {
//       if ( !OnlySize ) memcpy( &(DbKiBuf[Db_Xferred+i]), &(InData->ya[i]), sizeof(double)) ;
//       Db_Xferred++ ;
//     }
//     for ( i = 0 ; i < InData->za_Len ; i++ ) {
//       if ( !OnlySize ) memcpy( &(DbKiBuf[Db_Xferred+i]), &(InData->za[i]), sizeof(double)) ;
//       Db_Xferred++ ;
//     }
//     for ( i = 0 ; i < InData->Fx_connect_Len ; i++ ) {
//       if ( !OnlySize ) memcpy( &(DbKiBuf[Db_Xferred+i]), &(InData->Fx_connect[i]), sizeof(double)) ;
//       Db_Xferred++ ;
//     }
//     for ( i = 0 ; i < InData->Fy_connect_Len ; i++ ) {
//       if ( !OnlySize ) memcpy( &(DbKiBuf[Db_Xferred+i]), &(InData->Fy_connect[i]), sizeof(double)) ;
//       Db_Xferred++ ;
//     }
//     for ( i = 0 ; i < InData->Fz_connect_Len ; i++ ) {
//       if ( !OnlySize ) memcpy( &(DbKiBuf[Db_Xferred+i]), &(InData->Fz_connect[i]), sizeof(double)) ;
//       Db_Xferred++ ;
//     }
//     for ( i = 0 ; i < InData->Fx_anchor_Len ; i++ ) {
//       if ( !OnlySize ) memcpy( &(DbKiBuf[Db_Xferred+i]), &(InData->Fx_anchor[i]), sizeof(double)) ;
//       Db_Xferred++ ;
//     }
//     for ( i = 0 ; i < InData->Fy_anchor_Len ; i++ ) {
//       if ( !OnlySize ) memcpy( &(DbKiBuf[Db_Xferred+i]), &(InData->Fy_anchor[i]), sizeof(double)) ;
//       Db_Xferred++ ;
//     }
//     for ( i = 0 ; i < InData->Fz_anchor_Len ; i++ ) {
//       if ( !OnlySize ) memcpy( &(DbKiBuf[Db_Xferred+i]), &(InData->Fz_anchor[i]), sizeof(double)) ;
//       Db_Xferred++ ;
//     }
//   }
//   return(ErrStat) ;
// }
// 
// int
// C_MAP_UnpackOtherState( float * ReKiBuf,  
//                  double * DbKiBuf, 
//                  int * IntKiBuf,   
//                  MAP_OtherStateType_t *OutData, char * ErrMsg )
// {
//   int ErrStat ;
//   int Re_BufSz2 = 0 ;
//   int Db_BufSz2 = 0 ;
//   int Int_BufSz2 = 0 ;
//   int Re_Xferred = 0 ;
//   int Db_Xferred = 0 ;
//   int Int_Xferred = 0 ;
//   int Re_CurrSz = 0 ;
//   int Db_CurrSz = 0 ;
//   int Int_CurrSz = 0 ;
//   int one        = 1 ;
//   int i,i1,i2,i3,i4,i5 ;
//  // buffers to store meshes, if any
//   ReKiBuf = NULL ;
//   DbKiBuf = NULL ;
//   IntKiBuf = NULL ;
//   if ( OutData->H != NULL ) {
//     memcpy( OutData->H,&(DbKiBuf[ Db_Xferred ]),OutData->H_Len) ;
//     Db_Xferred   = Db_Xferred   + OutData->H_Len ; 
//   }
//   if ( OutData->V != NULL ) {
//     memcpy( OutData->V,&(DbKiBuf[ Db_Xferred ]),OutData->V_Len) ;
//     Db_Xferred   = Db_Xferred   + OutData->V_Len ; 
//   }
//   if ( OutData->Ha != NULL ) {
//     memcpy( OutData->Ha,&(DbKiBuf[ Db_Xferred ]),OutData->Ha_Len) ;
//     Db_Xferred   = Db_Xferred   + OutData->Ha_Len ; 
//   }
//   if ( OutData->Va != NULL ) {
//     memcpy( OutData->Va,&(DbKiBuf[ Db_Xferred ]),OutData->Va_Len) ;
//     Db_Xferred   = Db_Xferred   + OutData->Va_Len ; 
//   }
//   if ( OutData->x != NULL ) {
//     memcpy( OutData->x,&(DbKiBuf[ Db_Xferred ]),OutData->x_Len) ;
//     Db_Xferred   = Db_Xferred   + OutData->x_Len ; 
//   }
//   if ( OutData->y != NULL ) {
//     memcpy( OutData->y,&(DbKiBuf[ Db_Xferred ]),OutData->y_Len) ;
//     Db_Xferred   = Db_Xferred   + OutData->y_Len ; 
//   }
//   if ( OutData->z != NULL ) {
//     memcpy( OutData->z,&(DbKiBuf[ Db_Xferred ]),OutData->z_Len) ;
//     Db_Xferred   = Db_Xferred   + OutData->z_Len ; 
//   }
//   if ( OutData->xa != NULL ) {
//     memcpy( OutData->xa,&(DbKiBuf[ Db_Xferred ]),OutData->xa_Len) ;
//     Db_Xferred   = Db_Xferred   + OutData->xa_Len ; 
//   }
//   if ( OutData->ya != NULL ) {
//     memcpy( OutData->ya,&(DbKiBuf[ Db_Xferred ]),OutData->ya_Len) ;
//     Db_Xferred   = Db_Xferred   + OutData->ya_Len ; 
//   }
//   if ( OutData->za != NULL ) {
//     memcpy( OutData->za,&(DbKiBuf[ Db_Xferred ]),OutData->za_Len) ;
//     Db_Xferred   = Db_Xferred   + OutData->za_Len ; 
//   }
//   if ( OutData->Fx_connect != NULL ) {
//     memcpy( OutData->Fx_connect,&(DbKiBuf[ Db_Xferred ]),OutData->Fx_connect_Len) ;
//     Db_Xferred   = Db_Xferred   + OutData->Fx_connect_Len ; 
//   }
//   if ( OutData->Fy_connect != NULL ) {
//     memcpy( OutData->Fy_connect,&(DbKiBuf[ Db_Xferred ]),OutData->Fy_connect_Len) ;
//     Db_Xferred   = Db_Xferred   + OutData->Fy_connect_Len ; 
//   }
//   if ( OutData->Fz_connect != NULL ) {
//     memcpy( OutData->Fz_connect,&(DbKiBuf[ Db_Xferred ]),OutData->Fz_connect_Len) ;
//     Db_Xferred   = Db_Xferred   + OutData->Fz_connect_Len ; 
//   }
//   if ( OutData->Fx_anchor != NULL ) {
//     memcpy( OutData->Fx_anchor,&(DbKiBuf[ Db_Xferred ]),OutData->Fx_anchor_Len) ;
//     Db_Xferred   = Db_Xferred   + OutData->Fx_anchor_Len ; 
//   }
//   if ( OutData->Fy_anchor != NULL ) {
//     memcpy( OutData->Fy_anchor,&(DbKiBuf[ Db_Xferred ]),OutData->Fy_anchor_Len) ;
//     Db_Xferred   = Db_Xferred   + OutData->Fy_anchor_Len ; 
//   }
//   if ( OutData->Fz_anchor != NULL ) {
//     memcpy( OutData->Fz_anchor,&(DbKiBuf[ Db_Xferred ]),OutData->Fz_anchor_Len) ;
//     Db_Xferred   = Db_Xferred   + OutData->Fz_anchor_Len ; 
//   }
//   if ( ReKiBuf != NULL )  free(ReKiBuf) ;
//   if ( DbKiBuf != NULL )  free(DbKiBuf) ;
//   if ( IntKiBuf != NULL ) free(IntKiBuf) ;
//   return(ErrStat) ;
// }
// 
// CALL void MAP_F2C_OtherState_H_C ( MAP_OtherStateType_t *type, double *arr, int len )
// {
//   int i = 0;
//   for( i=0 ; i<=len-1 ; i++ ) type->H[i] = arr[i];
// }
// 
// CALL void MAP_F2C_OtherState_V_C ( MAP_OtherStateType_t *type, double *arr, int len )
// {
//   int i = 0;
//   for( i=0 ; i<=len-1 ; i++ ) type->V[i] = arr[i];
// }
// 
// CALL void MAP_F2C_OtherState_Ha_C ( MAP_OtherStateType_t *type, double *arr, int len )
// {
//   int i = 0;
//   for( i=0 ; i<=len-1 ; i++ ) type->Ha[i] = arr[i];
// }
// 
// CALL void MAP_F2C_OtherState_Va_C ( MAP_OtherStateType_t *type, double *arr, int len )
// {
//   int i = 0;
//   for( i=0 ; i<=len-1 ; i++ ) type->Va[i] = arr[i];
// }
// 
// CALL void MAP_F2C_OtherState_x_C ( MAP_OtherStateType_t *type, double *arr, int len )
// {
//   int i = 0;
//   for( i=0 ; i<=len-1 ; i++ ) type->x[i] = arr[i];
// }
// 
// CALL void MAP_F2C_OtherState_y_C ( MAP_OtherStateType_t *type, double *arr, int len )
// {
//   int i = 0;
//   for( i=0 ; i<=len-1 ; i++ ) type->y[i] = arr[i];
// }
// 
// CALL void MAP_F2C_OtherState_z_C ( MAP_OtherStateType_t *type, double *arr, int len )
// {
//   int i = 0;
//   for( i=0 ; i<=len-1 ; i++ ) type->z[i] = arr[i];
// }
// 
// CALL void MAP_F2C_OtherState_xa_C ( MAP_OtherStateType_t *type, double *arr, int len )
// {
//   int i = 0;
//   for( i=0 ; i<=len-1 ; i++ ) type->xa[i] = arr[i];
// }
// 
// CALL void MAP_F2C_OtherState_ya_C ( MAP_OtherStateType_t *type, double *arr, int len )
// {
//   int i = 0;
//   for( i=0 ; i<=len-1 ; i++ ) type->ya[i] = arr[i];
// }
// 
// CALL void MAP_F2C_OtherState_za_C ( MAP_OtherStateType_t *type, double *arr, int len )
// {
//   int i = 0;
//   for( i=0 ; i<=len-1 ; i++ ) type->za[i] = arr[i];
// }
// 
// CALL void MAP_F2C_OtherState_Fx_connect_C ( MAP_OtherStateType_t *type, double *arr, int len )
// {
//   int i = 0;
//   for( i=0 ; i<=len-1 ; i++ ) type->Fx_connect[i] = arr[i];
// }
// 
// CALL void MAP_F2C_OtherState_Fy_connect_C ( MAP_OtherStateType_t *type, double *arr, int len )
// {
//   int i = 0;
//   for( i=0 ; i<=len-1 ; i++ ) type->Fy_connect[i] = arr[i];
// }
// 
// CALL void MAP_F2C_OtherState_Fz_connect_C ( MAP_OtherStateType_t *type, double *arr, int len )
// {
//   int i = 0;
//   for( i=0 ; i<=len-1 ; i++ ) type->Fz_connect[i] = arr[i];
// }
// 
// CALL void MAP_F2C_OtherState_Fx_anchor_C ( MAP_OtherStateType_t *type, double *arr, int len )
// {
//   int i = 0;
//   for( i=0 ; i<=len-1 ; i++ ) type->Fx_anchor[i] = arr[i];
// }
// 
// CALL void MAP_F2C_OtherState_Fy_anchor_C ( MAP_OtherStateType_t *type, double *arr, int len )
// {
//   int i = 0;
//   for( i=0 ; i<=len-1 ; i++ ) type->Fy_anchor[i] = arr[i];
// }
// 
// CALL void MAP_F2C_OtherState_Fz_anchor_C ( MAP_OtherStateType_t *type, double *arr, int len )
// {
//   int i = 0;
//   for( i=0 ; i<=len-1 ; i++ ) type->Fz_anchor[i] = arr[i];
// }
// 
// int
// C_MAP_PackConstrState( float * ReKiBuf,  int * Re_BufSz ,
//                  double * DbKiBuf, int * Db_BufSz ,
//                  int * IntKiBuf,   int * Int_BufSz ,
//                  MAP_ConstraintStateType_t *InData, char * ErrMsg, int *SizeOnly )
// {
//   int ErrStat ;
//   int OnlySize ;
//   int Re_BufSz2 ;
//   int Db_BufSz2 ;
//   int Int_BufSz2 ;
//   int Re_Xferred = 0 ;
//   int Db_Xferred = 0 ;
//   int Int_Xferred = 0 ;
//   int one         = 1 ;
//   int i,i1,i2,i3,i4,i5 ;
//  // buffers to store meshes and subtypes, if any
// 
//   OnlySize = *SizeOnly ;
// 
//   *Re_BufSz = 0 ;
//   *Db_BufSz = 0 ;
//   *Int_BufSz = 0 ;
//   ReKiBuf = NULL ;
//   DbKiBuf = NULL ;
//   IntKiBuf = NULL ;
//   *Db_BufSz   += InData->H_Len ; // H 
//   *Db_BufSz   += InData->V_Len ; // V 
//   *Db_BufSz   += InData->x_Len ; // x 
//   *Db_BufSz   += InData->y_Len ; // y 
//   *Db_BufSz   += InData->z_Len ; // z 
//   if ( ! OnlySize ) {
//     if ( *Re_BufSz > 0  ) ReKiBuf  = (float  *)malloc(*Re_BufSz*sizeof(float) ) ;
//     if ( *Db_BufSz > 0  ) DbKiBuf  = (double *)malloc(*Db_BufSz*sizeof(double) ) ;
//     if ( *Int_BufSz > 0 ) IntKiBuf = (int    *)malloc(*Int_BufSz*sizeof(int) ) ;
//     for ( i = 0 ; i < InData->H_Len ; i++ ) {
//       if ( !OnlySize ) memcpy( &(DbKiBuf[Db_Xferred+i]), &(InData->H[i]), sizeof(double)) ;
//       Db_Xferred++ ;
//     }
//     for ( i = 0 ; i < InData->V_Len ; i++ ) {
//       if ( !OnlySize ) memcpy( &(DbKiBuf[Db_Xferred+i]), &(InData->V[i]), sizeof(double)) ;
//       Db_Xferred++ ;
//     }
//     for ( i = 0 ; i < InData->x_Len ; i++ ) {
//       if ( !OnlySize ) memcpy( &(DbKiBuf[Db_Xferred+i]), &(InData->x[i]), sizeof(double)) ;
//       Db_Xferred++ ;
//     }
//     for ( i = 0 ; i < InData->y_Len ; i++ ) {
//       if ( !OnlySize ) memcpy( &(DbKiBuf[Db_Xferred+i]), &(InData->y[i]), sizeof(double)) ;
//       Db_Xferred++ ;
//     }
//     for ( i = 0 ; i < InData->z_Len ; i++ ) {
//       if ( !OnlySize ) memcpy( &(DbKiBuf[Db_Xferred+i]), &(InData->z[i]), sizeof(double)) ;
//       Db_Xferred++ ;
//     }
//   }
//   return(ErrStat) ;
// }
// 
// int
// C_MAP_UnpackConstrState( float * ReKiBuf,  
//                  double * DbKiBuf, 
//                  int * IntKiBuf,   
//                  MAP_ConstraintStateType_t *OutData, char * ErrMsg )
// {
//   int ErrStat ;
//   int Re_BufSz2 = 0 ;
//   int Db_BufSz2 = 0 ;
//   int Int_BufSz2 = 0 ;
//   int Re_Xferred = 0 ;
//   int Db_Xferred = 0 ;
//   int Int_Xferred = 0 ;
//   int Re_CurrSz = 0 ;
//   int Db_CurrSz = 0 ;
//   int Int_CurrSz = 0 ;
//   int one        = 1 ;
//   int i,i1,i2,i3,i4,i5 ;
//  // buffers to store meshes, if any
//   ReKiBuf = NULL ;
//   DbKiBuf = NULL ;
//   IntKiBuf = NULL ;
//   if ( OutData->H != NULL ) {
//     memcpy( OutData->H,&(DbKiBuf[ Db_Xferred ]),OutData->H_Len) ;
//     Db_Xferred   = Db_Xferred   + OutData->H_Len ; 
//   }
//   if ( OutData->V != NULL ) {
//     memcpy( OutData->V,&(DbKiBuf[ Db_Xferred ]),OutData->V_Len) ;
//     Db_Xferred   = Db_Xferred   + OutData->V_Len ; 
//   }
//   if ( OutData->x != NULL ) {
//     memcpy( OutData->x,&(DbKiBuf[ Db_Xferred ]),OutData->x_Len) ;
//     Db_Xferred   = Db_Xferred   + OutData->x_Len ; 
//   }
//   if ( OutData->y != NULL ) {
//     memcpy( OutData->y,&(DbKiBuf[ Db_Xferred ]),OutData->y_Len) ;
//     Db_Xferred   = Db_Xferred   + OutData->y_Len ; 
//   }
//   if ( OutData->z != NULL ) {
//     memcpy( OutData->z,&(DbKiBuf[ Db_Xferred ]),OutData->z_Len) ;
//     Db_Xferred   = Db_Xferred   + OutData->z_Len ; 
//   }
//   if ( ReKiBuf != NULL )  free(ReKiBuf) ;
//   if ( DbKiBuf != NULL )  free(DbKiBuf) ;
//   if ( IntKiBuf != NULL ) free(IntKiBuf) ;
//   return(ErrStat) ;
// }
// 
// CALL void MAP_F2C_ConstrState_H_C ( MAP_ConstraintStateType_t *type, double *arr, int len )
// {
//   int i = 0;
//   for( i=0 ; i<=len-1 ; i++ ) type->H[i] = arr[i];
// }
// 
// CALL void MAP_F2C_ConstrState_V_C ( MAP_ConstraintStateType_t *type, double *arr, int len )
// {
//   int i = 0;
//   for( i=0 ; i<=len-1 ; i++ ) type->V[i] = arr[i];
// }
// 
// CALL void MAP_F2C_ConstrState_x_C ( MAP_ConstraintStateType_t *type, double *arr, int len )
// {
//   int i = 0;
//   for( i=0 ; i<=len-1 ; i++ ) type->x[i] = arr[i];
// }
// 
// CALL void MAP_F2C_ConstrState_y_C ( MAP_ConstraintStateType_t *type, double *arr, int len )
// {
//   int i = 0;
//   for( i=0 ; i<=len-1 ; i++ ) type->y[i] = arr[i];
// }
// 
// CALL void MAP_F2C_ConstrState_z_C ( MAP_ConstraintStateType_t *type, double *arr, int len )
// {
//   int i = 0;
//   for( i=0 ; i<=len-1 ; i++ ) type->z[i] = arr[i];
// }
// 
// int
// C_MAP_PackParam( float * ReKiBuf,  int * Re_BufSz ,
//                  double * DbKiBuf, int * Db_BufSz ,
//                  int * IntKiBuf,   int * Int_BufSz ,
//                  MAP_ParameterType_t *InData, char * ErrMsg, int *SizeOnly )
// {
//   int ErrStat ;
//   int OnlySize ;
//   int Re_BufSz2 ;
//   int Db_BufSz2 ;
//   int Int_BufSz2 ;
//   int Re_Xferred = 0 ;
//   int Db_Xferred = 0 ;
//   int Int_Xferred = 0 ;
//   int one         = 1 ;
//   int i,i1,i2,i3,i4,i5 ;
//  // buffers to store meshes and subtypes, if any
// 
//   OnlySize = *SizeOnly ;
// 
//   *Re_BufSz = 0 ;
//   *Db_BufSz = 0 ;
//   *Int_BufSz = 0 ;
//   ReKiBuf = NULL ;
//   DbKiBuf = NULL ;
//   IntKiBuf = NULL ;
//   *Db_BufSz   += 1  ; // g
//   *Db_BufSz   += 1  ; // depth
//   *Db_BufSz   += 1  ; // rho_sea
//   if ( ! OnlySize ) {
//     if ( *Re_BufSz > 0  ) ReKiBuf  = (float  *)malloc(*Re_BufSz*sizeof(float) ) ;
//     if ( *Db_BufSz > 0  ) DbKiBuf  = (double *)malloc(*Db_BufSz*sizeof(double) ) ;
//     if ( *Int_BufSz > 0 ) IntKiBuf = (int    *)malloc(*Int_BufSz*sizeof(int) ) ;
//     DbKiBuf[Db_Xferred++] = InData->g ;
//     DbKiBuf[Db_Xferred++] = InData->depth ;
//     DbKiBuf[Db_Xferred++] = InData->rho_sea ;
//   }
//   return(ErrStat) ;
// }
// 
// int
// C_MAP_UnpackParam( float * ReKiBuf,  
//                  double * DbKiBuf, 
//                  int * IntKiBuf,   
//                  MAP_ParameterType_t *OutData, char * ErrMsg )
// {
//   int ErrStat ;
//   int Re_BufSz2 = 0 ;
//   int Db_BufSz2 = 0 ;
//   int Int_BufSz2 = 0 ;
//   int Re_Xferred = 0 ;
//   int Db_Xferred = 0 ;
//   int Int_Xferred = 0 ;
//   int Re_CurrSz = 0 ;
//   int Db_CurrSz = 0 ;
//   int Int_CurrSz = 0 ;
//   int one        = 1 ;
//   int i,i1,i2,i3,i4,i5 ;
//  // buffers to store meshes, if any
//   ReKiBuf = NULL ;
//   DbKiBuf = NULL ;
//   IntKiBuf = NULL ;
//   OutData->g = DbKiBuf [ Db_Xferred ] ; 
//   Db_Xferred   = Db_Xferred   + 1 ; 
//   OutData->depth = DbKiBuf [ Db_Xferred ] ; 
//   Db_Xferred   = Db_Xferred   + 1 ; 
//   OutData->rho_sea = DbKiBuf [ Db_Xferred ] ; 
//   Db_Xferred   = Db_Xferred   + 1 ; 
//   if ( ReKiBuf != NULL )  free(ReKiBuf) ;
//   if ( DbKiBuf != NULL )  free(DbKiBuf) ;
//   if ( IntKiBuf != NULL ) free(IntKiBuf) ;
//   return(ErrStat) ;
// }
// 
// int
// C_MAP_PackInput( float * ReKiBuf,  int * Re_BufSz ,
//                  double * DbKiBuf, int * Db_BufSz ,
//                  int * IntKiBuf,   int * Int_BufSz ,
//                  MAP_InputType_t *InData, char * ErrMsg, int *SizeOnly )
// {
//   int ErrStat ;
//   int OnlySize ;
//   int Re_BufSz2 ;
//   int Db_BufSz2 ;
//   int Int_BufSz2 ;
//   int Re_Xferred = 0 ;
//   int Db_Xferred = 0 ;
//   int Int_Xferred = 0 ;
//   int one         = 1 ;
//   int i,i1,i2,i3,i4,i5 ;
//  // buffers to store meshes and subtypes, if any
//   float   * Re_PtFairDisplacement_Buf ;
//   double  * Db_PtFairDisplacement_Buf ;
//   int     * Int_PtFairDisplacement_Buf ;
// 
//   OnlySize = *SizeOnly ;
// 
//   *Re_BufSz = 0 ;
//   *Db_BufSz = 0 ;
//   *Int_BufSz = 0 ;
//   ReKiBuf = NULL ;
//   DbKiBuf = NULL ;
//   IntKiBuf = NULL ;
//   *Db_BufSz   += InData->x_Len ; // x 
//   *Db_BufSz   += InData->y_Len ; // y 
//   *Db_BufSz   += InData->z_Len ; // z 
//   if ( ! OnlySize ) {
//     if ( *Re_BufSz > 0  ) ReKiBuf  = (float  *)malloc(*Re_BufSz*sizeof(float) ) ;
//     if ( *Db_BufSz > 0  ) DbKiBuf  = (double *)malloc(*Db_BufSz*sizeof(double) ) ;
//     if ( *Int_BufSz > 0 ) IntKiBuf = (int    *)malloc(*Int_BufSz*sizeof(int) ) ;
//     for ( i = 0 ; i < InData->x_Len ; i++ ) {
//       if ( !OnlySize ) memcpy( &(DbKiBuf[Db_Xferred+i]), &(InData->x[i]), sizeof(double)) ;
//       Db_Xferred++ ;
//     }
//     for ( i = 0 ; i < InData->y_Len ; i++ ) {
//       if ( !OnlySize ) memcpy( &(DbKiBuf[Db_Xferred+i]), &(InData->y[i]), sizeof(double)) ;
//       Db_Xferred++ ;
//     }
//     for ( i = 0 ; i < InData->z_Len ; i++ ) {
//       if ( !OnlySize ) memcpy( &(DbKiBuf[Db_Xferred+i]), &(InData->z[i]), sizeof(double)) ;
//       Db_Xferred++ ;
//     }
//   }
//   return(ErrStat) ;
// }
// 
// int
// C_MAP_UnpackInput( float * ReKiBuf,  
//                  double * DbKiBuf, 
//                  int * IntKiBuf,   
//                  MAP_InputType_t *OutData, char * ErrMsg )
// {
//   int ErrStat ;
//   int Re_BufSz2 = 0 ;
//   int Db_BufSz2 = 0 ;
//   int Int_BufSz2 = 0 ;
//   int Re_Xferred = 0 ;
//   int Db_Xferred = 0 ;
//   int Int_Xferred = 0 ;
//   int Re_CurrSz = 0 ;
//   int Db_CurrSz = 0 ;
//   int Int_CurrSz = 0 ;
//   int one        = 1 ;
//   int i,i1,i2,i3,i4,i5 ;
//  // buffers to store meshes, if any
//   float   * Re_PtFairDisplacement_Buf ;
//   double  * Db_PtFairDisplacement_Buf ;
//   int     * Int_PtFairDisplacement_Buf ;
//   ReKiBuf = NULL ;
//   DbKiBuf = NULL ;
//   IntKiBuf = NULL ;
//   if ( OutData->x != NULL ) {
//     memcpy( OutData->x,&(DbKiBuf[ Db_Xferred ]),OutData->x_Len) ;
//     Db_Xferred   = Db_Xferred   + OutData->x_Len ; 
//   }
//   if ( OutData->y != NULL ) {
//     memcpy( OutData->y,&(DbKiBuf[ Db_Xferred ]),OutData->y_Len) ;
//     Db_Xferred   = Db_Xferred   + OutData->y_Len ; 
//   }
//   if ( OutData->z != NULL ) {
//     memcpy( OutData->z,&(DbKiBuf[ Db_Xferred ]),OutData->z_Len) ;
//     Db_Xferred   = Db_Xferred   + OutData->z_Len ; 
//   }
//   if ( ReKiBuf != NULL )  free(ReKiBuf) ;
//   if ( DbKiBuf != NULL )  free(DbKiBuf) ;
//   if ( IntKiBuf != NULL ) free(IntKiBuf) ;
//   return(ErrStat) ;
// }
// 
// CALL void MAP_F2C_Input_x_C ( MAP_InputType_t *type, double *arr, int len )
// {
//   int i = 0;
//   for( i=0 ; i<=len-1 ; i++ ) type->x[i] = arr[i];
// }
// 
// CALL void MAP_F2C_Input_y_C ( MAP_InputType_t *type, double *arr, int len )
// {
//   int i = 0;
//   for( i=0 ; i<=len-1 ; i++ ) type->y[i] = arr[i];
// }
// 
// CALL void MAP_F2C_Input_z_C ( MAP_InputType_t *type, double *arr, int len )
// {
//   int i = 0;
//   for( i=0 ; i<=len-1 ; i++ ) type->z[i] = arr[i];
// }
// 
// int
// C_MAP_PackOutput( float * ReKiBuf,  int * Re_BufSz ,
//                  double * DbKiBuf, int * Db_BufSz ,
//                  int * IntKiBuf,   int * Int_BufSz ,
//                  MAP_OutputType_t *InData, char * ErrMsg, int *SizeOnly )
// {
//   int ErrStat ;
//   int OnlySize ;
//   int Re_BufSz2 ;
//   int Db_BufSz2 ;
//   int Int_BufSz2 ;
//   int Re_Xferred = 0 ;
//   int Db_Xferred = 0 ;
//   int Int_Xferred = 0 ;
//   int one         = 1 ;
//   int i,i1,i2,i3,i4,i5 ;
//  // buffers to store meshes and subtypes, if any
//   float   * Re_ptFairleadLoad_Buf ;
//   double  * Db_ptFairleadLoad_Buf ;
//   int     * Int_ptFairleadLoad_Buf ;
// 
//   OnlySize = *SizeOnly ;
// 
//   *Re_BufSz = 0 ;
//   *Db_BufSz = 0 ;
//   *Int_BufSz = 0 ;
//   ReKiBuf = NULL ;
//   DbKiBuf = NULL ;
//   IntKiBuf = NULL ;
//   *Db_BufSz   += InData->Fx_Len ; // Fx 
//   *Db_BufSz   += InData->Fy_Len ; // Fy 
//   *Db_BufSz   += InData->Fz_Len ; // Fz 
//   *Db_BufSz   += InData->wrtOutput_Len ; // wrtOutput 
//   if ( ! OnlySize ) {
//     if ( *Re_BufSz > 0  ) ReKiBuf  = (float  *)malloc(*Re_BufSz*sizeof(float) ) ;
//     if ( *Db_BufSz > 0  ) DbKiBuf  = (double *)malloc(*Db_BufSz*sizeof(double) ) ;
//     if ( *Int_BufSz > 0 ) IntKiBuf = (int    *)malloc(*Int_BufSz*sizeof(int) ) ;
//     for ( i = 0 ; i < InData->Fx_Len ; i++ ) {
//       if ( !OnlySize ) memcpy( &(DbKiBuf[Db_Xferred+i]), &(InData->Fx[i]), sizeof(double)) ;
//       Db_Xferred++ ;
//     }
//     for ( i = 0 ; i < InData->Fy_Len ; i++ ) {
//       if ( !OnlySize ) memcpy( &(DbKiBuf[Db_Xferred+i]), &(InData->Fy[i]), sizeof(double)) ;
//       Db_Xferred++ ;
//     }
//     for ( i = 0 ; i < InData->Fz_Len ; i++ ) {
//       if ( !OnlySize ) memcpy( &(DbKiBuf[Db_Xferred+i]), &(InData->Fz[i]), sizeof(double)) ;
//       Db_Xferred++ ;
//     }
//     for ( i = 0 ; i < InData->wrtOutput_Len ; i++ ) {
//       if ( !OnlySize ) memcpy( &(DbKiBuf[Db_Xferred+i]), &(InData->wrtOutput[i]), sizeof(double)) ;
//       Db_Xferred++ ;
//     }
//   }
//   return(ErrStat) ;
// }
// 
// int
// C_MAP_UnpackOutput( float * ReKiBuf,  
//                  double * DbKiBuf, 
//                  int * IntKiBuf,   
//                  MAP_OutputType_t *OutData, char * ErrMsg )
// {
//   int ErrStat ;
//   int Re_BufSz2 = 0 ;
//   int Db_BufSz2 = 0 ;
//   int Int_BufSz2 = 0 ;
//   int Re_Xferred = 0 ;
//   int Db_Xferred = 0 ;
//   int Int_Xferred = 0 ;
//   int Re_CurrSz = 0 ;
//   int Db_CurrSz = 0 ;
//   int Int_CurrSz = 0 ;
//   int one        = 1 ;
//   int i,i1,i2,i3,i4,i5 ;
//  // buffers to store meshes, if any
//   float   * Re_ptFairleadLoad_Buf ;
//   double  * Db_ptFairleadLoad_Buf ;
//   int     * Int_ptFairleadLoad_Buf ;
//   ReKiBuf = NULL ;
//   DbKiBuf = NULL ;
//   IntKiBuf = NULL ;
//   if ( OutData->Fx != NULL ) {
//     memcpy( OutData->Fx,&(DbKiBuf[ Db_Xferred ]),OutData->Fx_Len) ;
//     Db_Xferred   = Db_Xferred   + OutData->Fx_Len ; 
//   }
//   if ( OutData->Fy != NULL ) {
//     memcpy( OutData->Fy,&(DbKiBuf[ Db_Xferred ]),OutData->Fy_Len) ;
//     Db_Xferred   = Db_Xferred   + OutData->Fy_Len ; 
//   }
//   if ( OutData->Fz != NULL ) {
//     memcpy( OutData->Fz,&(DbKiBuf[ Db_Xferred ]),OutData->Fz_Len) ;
//     Db_Xferred   = Db_Xferred   + OutData->Fz_Len ; 
//   }
//   if ( OutData->wrtOutput != NULL ) {
//     memcpy( OutData->wrtOutput,&(DbKiBuf[ Db_Xferred ]),OutData->wrtOutput_Len) ;
//     Db_Xferred   = Db_Xferred   + OutData->wrtOutput_Len ; 
//   }
//   if ( ReKiBuf != NULL )  free(ReKiBuf) ;
//   if ( DbKiBuf != NULL )  free(DbKiBuf) ;
//   if ( IntKiBuf != NULL ) free(IntKiBuf) ;
//   return(ErrStat) ;
// }
// 
// CALL void MAP_F2C_Output_Fx_C ( MAP_OutputType_t *type, double *arr, int len )
// {
//   int i = 0;
//   for( i=0 ; i<=len-1 ; i++ ) type->Fx[i] = arr[i];
// }
// 
// CALL void MAP_F2C_Output_Fy_C ( MAP_OutputType_t *type, double *arr, int len )
// {
//   int i = 0;
//   for( i=0 ; i<=len-1 ; i++ ) type->Fy[i] = arr[i];
// }
// 
// CALL void MAP_F2C_Output_Fz_C ( MAP_OutputType_t *type, double *arr, int len )
// {
//   int i = 0;
//   for( i=0 ; i<=len-1 ; i++ ) type->Fz[i] = arr[i];
// }
// 
// CALL void MAP_F2C_Output_wrtOutput_C ( MAP_OutputType_t *type, double *arr, int len )
// {
//   int i = 0;
//   for( i=0 ; i<=len-1 ; i++ ) type->wrtOutput[i] = arr[i];
// }
//!ENDOFREGISTRYGENERATEDFILE
