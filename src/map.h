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


#ifndef _MAP_H
#define _MAP_H


#include "mapsys.h"
#include "simclist/simclist.h"
#include "bstring/bstrlib.h"


/**
 * @brief Associates the node with a particular type. Fix nodes are anchor points
 *        (ussually) and cannot move with time. Connect nodes are intermediaries
 *        connecting two lines. When a connect node is defined, the force-balance
 *        equation is minimized to obtain the cable profile X, Y, and Z node position
 *        at connect nodes are associated with constraint types. Vessel nodes attach
 *        to the platform. X, Y, and Z, node positions are associated with input 
 *        types, and the corresponding fx, fy, and fz forces are the outputs.  
 *        Moments are calculated by the calling program. 
 */
typedef enum NodeType_enum {
  NONE,     /**< */
  FIX,      /**< */
  CONNECT,  /**< */
  VESSEL    /**< */
} NodeType;


/**
 * @brief Finite different routine when solving the outer-loop iterations. 
 *        This is also used when the linearization routine is called. 
 */
typedef enum FdType_enum {
  BACKWARD_DIFFERENCE, /**< */
  CENTRAL_DIFFERENCE,  /**< */
  FORWARD_DIFFERENCE   /**< */
} FdType;


/**
 * @struct SolveType
 * @brief Set internally. If a connect node is present in the MAP input
 *        file, then the solve typeis partitioned. Partitioned describes the
 *        process of breaking the problem into two blocks: an inner-loop solve
 *        (where the catenary equations are minimized) and an outer-loop solve
 *        (where the force-balance equations are minimized). Likewise, when no
 *        connect nodes are present, the solver type is monolithic and only the
 *        non-linear catenary equations need to be solved. 
 */
typedef enum SolveType_enum {  
  MONOLITHIC,   /**< for MSQS, elements only (no connect nodes) */
  PARTITIONED,  /**< for MSQS system with connect nodes */
} SolveType;


/**
 * @brief Finite-difference structure. This is used locally withint FD routines to conveniently store the
 *        force and moment when vessel is displaced by epsilon. 
 */
struct Fd_t {
  double* fx;
  double* fy;
  double* fz;
  double* mx;
  double* my;
  double* mz;
}; typedef struct Fd_t Fd;


/**
 * @brief Fundamental MAP type. Unless it is used locally, every variable should be darclared as a VarType
 *        (or VarTypePtr, described below). The VarType provides convenience for printing variable 
 *        information, such as the units, it's name (such as 'FX[1]' for X-direction node 1 force to the 
 *        output file), and counting references.
 */
struct VarType_t {
  MapReal value;         /**< the value */
  bstring units;           /*< units for printing information to a summary file or output buffer */
  bstring name;            /*< name of the variable. This is used for identifying it in the output buffer */
  // char* units;           /*< units for printing information to a summary file or output buffer */
  // char* name;            /*< name of the variable. This is used for identifying it in the output buffer */
  bool isFixed;          /**< if isFixed = true, then we are not solving for this variable */
  int referenceCounter;  /**< for ensuring the variable is assigned to one of: input, param, or constraint */
  int id;                /**< node or element this value is attached to */
}; typedef struct VarType_t VarType;


/**
 * @brief Serves the same function as VarType, but treats value as a pointer. This preserves the FAST integration
 *        to native Fortran derivived types. Instead of the variable residing in C, the value parameter points
 *        to a variable allocated in Fortran. This feature is also preserved with Python binding. 
 */
struct VarTypePtr_t {
  MapReal* value;        /**< the value */
  bstring units;           /*< units for printing information to a summary file or output buffer */
  bstring name;            /*< name of the variable. This is used for identifying it in the output buffer */
  // char* units;           /**< units for printing information to a summary file or output buffer */
  // char* name;            /**< name of the variable. This is used for identifying it in the output buffer */
  bool isFixed;          /**< If isFixed = true, then we are not solving for this variable */
  int referenceCounter;  /**< For ensuring the variable is assigned to one of: input, param, or constraint */
  int id;                /**< node or element this value is attached to */
}; typedef struct VarTypePtr_t VarTypePtr;


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



/**
 * @brief References a list of VarType's (in the case of out_list) and VarTypePtr's
 *        (in the case of out_list_ptr) for output file data streaming. The idea is
 *        to loop through the list and print the 'value' and 'name[id]' members of
 *        VarType and VarTypePtr.
 * @see   
 * @todo  This really should be made polymorphic so only one variable is needed. 
 */
struct OutputList_t{
  list_t out_list;     /**< Outputs associated with VarType */
  list_t out_list_ptr; /**< Outputs associated with VarTypePtr */
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
  bstring label;            /**< Provides the string a recognizable name (such as 'nylon' or 'steel') */
  // char* label;              /**< Provides the string a recognizable name (such as 'nylon' or 'steel') */
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
  CableLibrary* lineProperty; /**< line properties */
  LineOptions options;        /**< run-time options flag */
  VarTypePtr H;               /**< Horizontal fairlead force in the local cable elemenet frame */
  VarTypePtr V;               /**< Vertical fairlead force in the local cable elemenet frame */  
  MapReal* lineTension;       /**< array of line tension along 's' [N] */ 
  MapReal damageTime;         /**< time to damage this element and return zero force to the glue code */
  MapReal residualNorm;       
  VarType psi;                /**< angle of roation between global X and local x axis [deg] */
  VarType alpha;              /**< angle of inclication [deg] */
  VarType alphaAtAnchor;      /**< angle of inclication at anchor [deg] */
  VarType l;                  /**< horizontal cable excursion [m] */
  VarType lb;                 /**< length of element touching the seabed [m] */
  VarType Lu;                 /**< unstretched cable length [m] */
  VarType h;                  /**< vertical cable excursion [m] */
  VarType HAtAnchor;          /**< Horizontal anchor force in the local cable elemenet frame */
  VarType VAtAnchor;          /**< Vertical anchor force in the local cable elemenet frame */
  VarType T;                  /**< Tension magnitude [N] */
  VarType TAtAnchor;          /**< Tension magnitude at anchor [N] */
  Force forceAtFairlead;      //*< @rm is this even necessary? I don't think so. Element should not store node forces. They only can contribute force in sumForcePrt
  Force forceAtAnchor;        //*< @rm is this even necessary? I don't think so. Element should not store node forces. They only can contribute force in sumForcePrt
  // char* label;                /**< reference a pre-defined property in the line dictionary */
  bstring label;              /**< reference a pre-defined property in the line dictionary */
  Node* anchor;               /**< Anchor node */
  Node* fairlead;             /**< Fairlead node */
  int segmentSize;
  int diagnosticType;         /**< none=0, first iteration only=2, all iterations otherwise */
  int numFuncEvals;           /**< number of function evaluations */ 
  int numJacEvals;            /**< number of function evaluations */      
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
  MapReal* repeat_angle;
  MapReal integrationDt; /**< Integration time step [sec]. LM model specific */
  MapReal kbLm;          /**< Seabed stiffeness coefficient [N/m]. LM model specific */
  MapReal cbLm;          /**< Seabed damping parameter [N-s/m]. LM model specific */
  bool waveKinematics;   /**< Enable wave kinematics o calculated relative flui velcity. LM model specific */
  int repeat_angle_size;
}; typedef struct ModelOptions_t ModelOptions;


struct InnerSolveAttributes_t {
  MapReal f_tol; 
  MapReal g_tol; 
  MapReal x_tol; 
  int max_its;   
  MapReal** node_jac;
  MapReal x[2];               /**< array of variables the length of n */
  MapReal fvec[2];            /**< function evaluations (residual) */
  MapReal fjac[4];            /**< jacobian. This is a little convoluted because the jacobian is not an array */
  MapReal wa1[2];             /**< work array of length n */
  MapReal wa2[2];             /**< work array of length n */
  MapReal wa3[2];             /**< work array of length n */
  MapReal wa4[2];             /**< work array of length m */  
  double diag[2];            
  double qtf[2];             
  double factor;             
  int ldfjac;                 /**< number of columns in fjac */
  int mode;             
  int nprint;           
  int info;             
  int ipvt[2];
  int m;                      /**< number of functions */ 
  int n;                      /**< number of variables */
}; typedef struct InnerSolveAttributes_t InnerSolveAttributes;


struct OuterSolveAttributes_t {
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
}; typedef struct OuterSolveAttributes_t OuterSolveAttributes;


struct ModelData_t {
  SolveType MAP_SOLVE_TYPE;        /**< Identifies the solver type: single line, partitioned (multisegmented), and lumped-mass/FEA. Initialized in {@link initialize_model_data_to_null} */
  ModelOptions modelOptions;       /**< Contains global model options. Default setting in {@link initialize_model_options_to_defaults} */
  InnerSolveAttributes inner_loop; /**< Inner-loop (line level) solver options. Default settings in {@link initialize_solver_data_to_null} */
  OuterSolveAttributes outer_loop; /**< Outer-loop (node level) solver options. Default settings in {@link initialize_model_options_to_defaults} */
  OutputList* yList;               /**< Output stream. Set to null at initialization */
  Vessel vessel;                   /**< Vessel for the mooring instance. Initialized in {@link initialize_vessel_to_null}. Associated VarType's are set in {@link set_vessel} */
  list_t cableLibrary;             /**< Cable library link list; stores cable properties, e.g., @see CableLibrary_t */
  list_t element;                  /**< Line link list */
  list_t node;                     /**< Node link list */
  int sizeOfCableLibrary;          /**< Number of cable types defined in the cable library section of the input file */ 
  int sizeOfElements;              /**< Number of lines */
  int sizeOfNodes;                 /**< Number of nodes */
}; typedef struct ModelData_t ModelData;


/**
 * @details MAP options from parsed input file. Note that MAP does not readon the input file. This is done by the calling program.
 *          The calling program simply sets library_input_string, node_input_string, element_input_string, and solver_options_string.
 *          MAP then parses this string and expands them if necessary depending on the '{@link ModelOptions_t}' repeat_angle flag.
 */
struct InitializationData_t {
  struct bstrList* libraryInputString;         /**< library property string from input file. MAP does not read contents from input string; must be done by calling program */
  struct bstrList* nodeInputString;            /**< raw (non-expanded) node input string. MAP does not read contents from input string; must be done by calling program */
  struct bstrList* elementInputString;         /**< raw (non-expanded) element input string(MAP does not read contents from input string; must be done by calling program */
  struct bstrList* solverOptionsString;        /**< model poptions input string */
  struct bstrList* expandedNodeInputString;    /**< full node input string duplicating information in nodeInputString when the 'repeat' flag is used */
  struct bstrList* expandedElementInputString; /**< full element input string duplicating information in nodeElementString when the 'repeat' flag is used */
  bstring summaryFileName;                     /**< summary file name. Can be set through {@link map_set_summary_file_name()} */
  // char** libraryInputString;         /**< library property string from input file. MAP does not read contents from input string; must be done by calling program */
  // char** nodeInputString;            /**< raw (non-expanded) node input string. MAP does not read contents from input string; must be done by calling program */
  // char** elementInputString;         /**< raw (non-expanded) element input string(MAP does not read contents from input string; must be done by calling program */
  // char** solverOptionsString;        /**< model poptions input string */
  // char** expandedNodeInputString;    /**< full node input string duplicating information in nodeInputString when the 'repeat' flag is used */
  // char** expandedElementInputString; /**< full element input string duplicating information in nodeElementString when the 'repeat' flag is used */
  // char* summaryFileName;             /**< summary file name. Can be set through {@link map_set_summary_file_name()} */
  // int sizeOfFullNodeString;          /**< number of node entries after expansion, i.e., after repeats */
  // int sizeOfFullElementString;       /**< number of element entries after expansion, i.e., after repeats */
  // int librarySize;                   /**< number of cable types defined in the cable library section of the input string */
  // int nodeSize;                      /**< number of nodes after expansion */
  // int elementSize;                   /**< number of elements after expansion */
  int solverOptionsSize;             /**< number of model options. This does not necessarily correspond to solver option, i.e., 'help' */
}; typedef struct InitializationData_t InitializationData;


#endif /* _MAP_H */
