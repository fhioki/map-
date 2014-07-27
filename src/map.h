/***************************************************************************
 *   Copyright (C) 2014 mdm                                                *
 *   marco.masciola@gmail.com                                              *
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


#ifndef _MAP_H
#define _MAP_H


#include "mapsys.h"
#include "simclist.h"


typedef enum {
  NONE,
  FIX,
  CONNECT,
  VESSEL
} NodeType;


typedef enum {
  BACKWARD_DIFFERENCE,
  CENTRAL_DIFFERENCE,
  FORWARD_DIFFERENCE  
} FdType;


typedef enum {  
  MONOLITHIC,   /* for MSQS, elements only (no connect nodes) */
  PARTITIONED,  /* for MSQS system with connect nodes */
} SolveType;


typedef struct Fd {
  double* fx;
  double* fy;
  double* fz;
  double* mx;
  double* my;
  double* mz;
} Fd;


typedef struct {
  MapReal value;   /* The value */
  char* units;   /* units for printing information to a summary file or output buffer */
  char* name;    /* name of the variable. This is used for identifying it in the output buffer */
  bool isFixed; /* If isFixed = true, then we are not solving for this var */
  int referenceCounter;  /* For ensuring the variable is assigned to one of: input, param, or constraint */
  int id;      /* node or element this value is attached to */
} VarType;


typedef struct {
  MapReal* value;   /* The value */
  char* units;   /* units for printing information to a summary file or output buffer */
  char* name;    /* name of the variable. This is used for identifying it in the output buffer */
  bool isFixed; /* If isFixed = true, then we are not solving for this var */
  int referenceCounter;  /* For ensuring the variable is assigned to one of: input, param, or constraint */
  int id;      /* node or element this value is attached to */
} VarTypePtr;


struct Vector_t {
  MapReal x;
  MapReal y;
  MapReal z;
}; typedef struct Vector_t Vector;
  

struct Point_t {
  VarType x;
  VarType y;
  VarType z;  
}; typedef struct Point_t Point;


struct PointPtr_t {
  VarTypePtr x;
  VarTypePtr y;
  VarTypePtr z;  
}; typedef struct PointPtr_t PointPtr;


struct EulerAngle_t {
  VarType phi;
  VarType the;
  VarType psi;  
}; typedef struct EulerAngle_t EulerAngle;


struct Force_t {
  VarType fx;
  VarType fy;
  VarType fz;  
}; typedef struct Force_t Force;


struct ForcePtr_t {
  VarTypePtr fx;
  VarTypePtr fy;
  VarTypePtr fz;  
}; typedef struct ForcePtr_t ForcePtr;


/**
 * @brief Central point where all 'VESSEL' nodes can be displaced. Instead of displacing all nodes individually, the vessel can be displaced,
 *        then helper functions can be called to displacement the nodes. The vessel only reference 'input' nodes.
 */
struct Vessel_t {
  Point displacement;     /**< User-specified vessel displacement. This is the [m]*/
  Point refOrigin;        /**< Center of rotation origin. The moments are taken about this is point. The reference point is with respect to the FAST reference origin (equal to the
                           * SWL at zero vessel dispalcements) [m]*/
  EulerAngle orientation; /**< Vessel orientation [deg]*/
  Force lineSumForce;     /**< Sum force of all nodes connecting to the vessel [N] */
  MapReal* xi;            /**< initial node connection point in body frame.This is equal to uType->x at initialization [m] */
  MapReal* yi;            /**< initial node connection point in body frame This is equal to uType->y at initialization [m] */
  MapReal* zi;            /**< initial node connection point in body frame This is equal to uType->z at initialization [m] */
}; typedef struct Vessel_t Vessel;


struct OutputList_t{
  list_t outList;
  list_t outListPtr;
}; typedef struct OutputList_t OutputList;


typedef struct {
  bool plotFlag;  // @rm
  bool gxPosFlag;
  bool gyPosFlag;
  bool gzPosFlag;
  bool gxAnchorPosFlag;
  bool gyAnchorPosFlag;
  bool gzAnchorPosFlag;
  bool gxForceFlag;
  bool gyForceFlag;
  bool gzForceFlag;
  bool HFlag;
  bool VFlag;
  bool VAnchorFlag;
  bool HAnchorFlag;
  bool fairleadTensionFlag;
  bool anchorTensionFlag;
  bool horizontalExcursionFlag;
  bool verticalExcursionFlag;
  bool azimuthFlag;
  bool altitudeFlag;
  bool altitudeAnchorFlag;
  bool lineTensionFlag;
  bool omitContact;
  bool layLengthFlag;
  bool damageTimeFlag;
  bool diagnosticsFlag;
} LineOptions;


/**
 * @struct CableLibrary_t
 * @brief Defines cable properties for a line. These values are fixed with time and connot change.
 */
struct CableLibrary_t {
  MapReal diam;             /**< Cable diameter, [m] */
  MapReal massDensityInAir; /**< Cable density in air [kg/m] */
  MapReal ea;               /**< Element stiffness [N] */
  MapReal omega;            /**< cable weight per length in seawater [N/m] */
  MapReal a;                /**< cross-sectional area [m^2] */
  MapReal cb;               /**< Cable/seabed friction coefficient [non-dimensional] */
  MapReal cIntDamp;         /**< Internal (structural) damping coefficient [non-dimensional] */
  MapReal cAdded;           /**< Added mass coefficient [non-dimensional] */
  MapReal cDragNormal;      /**< Quadtradice drag coefficient in the cable cross-flow direction [non-dimensional] */
  MapReal cDragTangent;     /**< Tangential drag oefficient [non-dimensional] */
  char* label;              /**< Provides the string a recognizable name (such as 'nylon' or 'steel') */
}; typedef struct CableLibrary_t CableLibrary;


struct Events_t {
  MapReal LuRestore;
  MapReal LuMax;
  MapReal dLu;
  bool payinFlag;  
}; typedef struct Events_t Events;


struct Node_t {
  NodeType type;
  VarType MApplied;
  VarType BApplied;
  PointPtr positionPtr; /* this is a Ptr because it points to a fortran type */
  ForcePtr sumForcePtr; /* this is a Ptr because it points to a fortran type */ 
  Force externalForce;    
}; typedef struct Node_t Node;


struct Element_t {
  Events event;
  CableLibrary* lineProperty; /* line properties */
  LineOptions options;        /* run-time options flag */
  VarTypePtr H;               /* Horizontal fairlead force in the local cable elemenet frame */
  VarTypePtr V;               /* Vertical fairlead force in the local cable elemenet frame */  
  MapReal* lineTension;       /* array of line tension along 's' [kN] */ 
  MapReal damageTime;         /* time to damage this element and return zero force to the glue code */
  MapReal residualNorm;
  VarType psi;                /* angle of roation between global X and local x axis [deg] */
  VarType alpha;              /* angle of inclication [deg] */
  VarType alphaAtAnchor;      /* angle of inclication at anchor [deg] */
  VarType l;                  /* horizontal cable excursion [m] */
  VarType lb;                 /* length of element touching the seabed [m] */
  VarType Lu;                 /* unstretched cable length [m] */
  VarType h;                  /* vertical cable excursion [m] */
  VarType HAtAnchor;          /* Horizontal anchor force in the local cable elemenet frame */
  VarType VAtAnchor;          /* Vertical anchor force in the local cable elemenet frame */
  VarType T;                  /* Tension magnitude [kN] */
  VarType TAtAnchor;          /* Tension magnitude at anchor [kN] */
  Force forceAtFairlead;      // @rm is this even necessary? I don't think so. Element should not store node forces. They only can contribute force in sumForcePrt
  Force forceAtAnchor;        // @rm is this even necessary? I don't think so. Element should not store node forces. They only can contribute force in sumForcePrt
  char* label;                /* reference a pre-defined property in the line dictionary */
  Node* anchor;               /* Anchor node */
  Node* fairlead;             /* Fairlead node */
  int segmentSize;
  int diagnosticType;         /* none=0, first iteration only=2, all iterations otherwise */
  int numFuncEvals;           /* number of function evaluations */ 
  int numJacEvals;            /* number of function evaluations */      
  int convergeReason;         /*   - info=0 : improper input parameters.
                               *   - info=1 : both actual and predicted relative reductions in the sum of squares are at most ftol.
                               *   - info=2 : relative error between two consecutive iterates is at most xtol.
                               *   - info=3 : conditions for info = 1 and info = 2 both hold.
                               *   - info=4 : the cosine of the angle between fvec and any column of the Jacobian is at most gtol in absolute value.
                               *   - info=5 : number of calls to fcn has reached or exceeded maxfev.
                               *   - info=6 : ftol is too small. No further reduction in the sum of squares is possible.
                               *   - info=7 : xtol is too small. No further improvement in the approximate solution x is possible.
                               *   - info=8 : gtol is too small. fvec is orthogonal to the columns of the Jacobian to machine precision.
                               */
}; typedef struct Element_t Element;


struct ModelOptions_t {
  MapReal* repeatAngles;
  MapReal innerFTol; /* @todo: this should be moved in MinPackDataInner */
  MapReal innerGTol; /* @todo: this should be moved in MinPackDataInner */
  MapReal innerXTol; /* @todo: this should be moved in MinPackDataInner */
  int innerMaxIts; /* @todo: this should be moved in MinPackDataInner */
  MapReal integrationDt; /* LM model specific */
  MapReal kbLm; /* LM model specific */
  MapReal cbLm; /* LM model specific */
  bool waveKinematics; /* LM model specific */
  int sizeOfRepeatAngles;
}; typedef struct ModelOptions_t ModelOptions;


struct MinPackDataInner_t {
  MapReal** nodeJac;
  MapReal x[2];               /* array of variables the length of n */
  MapReal fvec[2];            /* function evaluations (residual) */
  MapReal fjac[4];            /* jacobian. This is a little convoluted because the jacobian is not an array */
  MapReal wa1[2];             /* work array of length n */
  MapReal wa2[2];             /* work array of length n */
  MapReal wa3[2];             /* work array of length n */
  MapReal wa4[2];             /* work array of length m */  
  double diag[2];            
  double qtf[2];             
  double factor;             
  int ldfjac; /* number of columns in fjac */
  int mode;             
  int nprint;           
  int info;             
  int ipvt[2];
  int m; /* number of functions */ 
  int n; /* number of variables */
}; typedef struct MinPackDataInner_t MinPackDataInner;


struct MinPackDataOuter_t {
  FdType fd;
  MapReal** jac;
  MapReal** l;
  MapReal** u;
  MapReal* b;
  MapReal* x;
  MapReal* y;
  MapReal ds;
  MapReal d;
  MapReal epsilon;
  MapReal normError;
  MapReal tol;
  MapReal coef;
  bool pg;
  int maxIts;
  int iterationCount;
}; typedef struct MinPackDataOuter_t MinPackDataOuter;


struct ModelData_t {
  SolveType MAP_SOLVE_TYPE;
  ModelOptions modelOptions;
  MinPackDataInner solverData;
  MinPackDataOuter outerSolveData;
  OutputList* yList; /* this should be where the internal C data is stored */
  Vessel vessel;
  list_t cableLibrary; /**< Cable library link list; stores cable properties, e.g., @see CableLibrary_t */
  list_t element;
  list_t node;
  int sizeOfCableLibrary;
  int sizeOfElements;
  int sizeOfNodes;
}; typedef struct ModelData_t ModelData;


struct InitializationData_t{
  char** libraryInputString;
  char** nodeInputString;
  char** elementInputString;
  char** solverOptionsString;
  char** expandedNodeInputString;    /**< this is the full string duplicating information in nodeInputString when the 'repeat' flag is used */
  char** expandedElementInputString; /**< this is the full string duplicating information in nodeElementString when the 'repeat' flag is used */
  char* summaryFileName;
  int sizeOfFullNodeString; 
  int sizeOfFullElementString; 
  int librarySize;
  int nodeSize;
  int elementSize;
  int solverOptionsSize;
}; typedef struct InitializationData_t InitializationData;


#endif /* _MAP_H */
