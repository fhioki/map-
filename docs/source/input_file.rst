Input File
==========

A MAP input file initializes a mooring mode to provide essential information regarding the mooring line properties, declaring the number of nodes and lines in the system, association between nodes and lines, boundary conditions, and initial states.
We use the extension ``<*.map>`` to identify the MAP input file.
The sample MAP input deck given below::

	--------------- LINE DICTIONARY ----------------------------------------------
	LineType  Diam      MassDenInAir   EA        CB   CIntDamp  Ca    Cdn    Cdt
	(-)       (m)       (kg/m)        (N)        (-)   (Pa-s)   (-)   (-)    (-)
	steel     0.25       319.07     9817477042   1.0   -999.9 -999.9 -999.9 -999.9
	poly      0.30       98.96      989601685    1.0   -999.9 -999.9 -999.9 -999.9
	--------------- NODE PROPERTIES ----------------------------------------------
	Node Type       X       Y       Z      M     B     FX    FY    FZ
	(-)  (-)       (m)     (m)     (m)    (kg)  (mˆ3)  (N)   (N)   (N)
	1    fix     400        0     depth    0     0      #     #     #
	2    connect #90       #0    #-80      0     0      0     0     0   
	3    vessel   20        20    -10      0     0      #     #     #
	4    vessel   20       -20    -10      0     0      #     #     #
	--------------- LINE PROPERTIES ----------------------------------------------
	Line    LineType  UnstrLen  NodeAnch  NodeFair  Flags
	(-)      (-)       (m)       (-)       (-)       (-)
	1       steel      450        1         2  
	2       poly      #90         2         3  
	3       poly      #90         2         4  
	--------------- SOLVER OPTIONS------------------------------------------------
	Option
	(-)

acts as a template for defining parameters in a mooring system demonstrated here: 

.. figure:: nstatic/input1a_image.png
    :align: center
    :width: 400pt

.. Note::
   Environmental properties like water depth, sea density, and gavity constant are set by the calling program. 
   They are purposely absent in the MAP++ input file to prevent force imbalances from coefficient mismatches. 

The MAP input file is divided into four sections:

* **LINE DICTIONARY:** Defines the material properties of the line.

* **NODE PROPERTIES:** Defines boundary constraints and extensional limits.

* **LINE PROPERTIES:** Associates a line with material properties and connectivity between nodes. 

* **SOLVER OPTIONS:** Run-time options to engage different solve strategies.

Line Dictionary
---------------

   **LineType**     -- user-defined name [-]

   **Diam**         -- material diameter [m]

   **MassDenInAir** -- mass density in air [kg/m^3]

   **EA**           -- element axial stiffness [N/m]

   **CB**           -- cable/seabed friction coefficient [-]

   **CIntDamp**     -- unused

   **Ca**           -- unused

   **Cdn**          -- unused

   **Cdt**          -- unused

Node Properties
---------------

   **Node** -- node number (sequential).

   **Type** --  type of node, which can be one of ``FIX``, ``CONNECT``, or ``VESSEL``. Vessel implied the node motion is prescribed.

   **X**    -- X connection point relative to the global origin (if ``FIX`` or ``CONNECT``) or local vessel reference frame (if ``VESSEL``).  ``Connect`` nodes must be preceeded by a ``#`` is indicate this is as an initial guess.

   **Y**    -- Y connection point relative to the global origin (if ``FIX`` or ``CONNECT``) or local vessel reference frame (if ``VESSEL``).  ``Connect`` nodes must be preceeded by a ``#`` is indicate this is as an initial guess.

   **Z**    -- Z connection point relative to the global origin (if ``FIX`` or ``CONNECT``) or local vessel reference frame (if ``VESSEL``).  ``Connect`` nodes must be preceeded by a ``#`` is indicate this is as an initial guess.

   **M**    -- point mass applied to the node. The force appled to the node is :math:`M\times g` applied in the direction of gravity.

   **B**    -- displaced volume applied to node. The  force applied to the node is :math:`B\times \rho \times g` applied in the direction opposite of gravity. 

   **FX**   -- X direction external force applied to node if ``CONNECT``. In the node is a ``VESSEL`` or ``FIX`` type, the ``#`` key must be used to indicate this is as an initial guess/iterated value.

   **FY**   -- Y direction external force applied to node if ``CONNECT``. In the node is a ``VESSEL`` or ``FIX`` type, the ``#`` key must be used to indicate this is as an initial guess/iterated value.

   **FZ**   -- Z direction external force applied to node if ``CONNECT``. In the node is a ``VESSEL`` or ``FIX`` type, the ``#`` key must be used to indicate this is as an initial guess/iterated value.

Line Properties
---------------

   **Line**      -- line number (sequential).

   **LineType**  -- line type. Must be one type defined in ``LineType`` from dictionary.

   **UnstrLen**  -- unstretched line length.

   **NodeAnch**  -- anchor node number.

   **NodeFair**  -- fairlead node number.

   **Flags**     -- line flag (defined below).

Flags
----- 
Flags are applied to individual lines. 
These flags control the output text stream:

 * ``GX_POS`` - global X fairlead position [m]     
 * ``GY_POS`` - global Y fairlead position [m]  
 * ``GZ_POS`` - global Z fairlead position [m]  
 * ``GX_A_POS`` - global X position of anchor [m]
 * ``GY_A_POS`` - global Y position of anchor [m]  
 * ``GZ_A_POS`` - global Z position of anchor [m]
 * ``GX_FORCE`` - global X fairlead force [N]
 * ``GY_FORCE`` - global Y fairlead force [N]
 * ``GZ_FORCE`` - global Z fairlead force [N]
 * ``H_FAIR`` - horizontal (XY plane) fairlead force [N] 
 * ``H_ANCH`` - horizontal (XY plane) anchor force [N]  
 * ``V_FAIR`` - vertical (Z axis) fairlead force [N]
 * ``V_ANCH`` - vertical (Z axis) anchor force [N]  
 * ``TENSION_FAIR`` - fairlead force magnitude, [N] 
 * ``TENSION_ANCH`` - anchor force magnitude, [N]  
 * ``X_EXCURSION`` - line horizontal excursion [m]
 * ``Z_EXCURSION`` - line veritical excursion [m]
 * ``AZIMUTH`` - line azimuth angle with respect to the inertial reference frame [deg]
 * ``ALTITUDE`` - angle of declination at the fairlead [deg]
 * ``ALTITUDE_ANCH`` - line lift-off angle at the anchor [deg]

The follow flags enable/disable features for each line they are applied to:

 * ``LINE_TENSION`` - line tension force magnitude at fairlead [N]
 * ``OMIT_CONTACT`` - ignore seabed boundary and treat line as freely hanging
 * ``LINEAR_SPRING`` - model the line as a linear spring. Intended for taut lines. 
 * ``LAY_LENGTH`` - amount of line laying on the seabed [m]
 * ``DAMAGE_TIME`` - time [sec] to disconnect fairlead from node
 * ``DIAGNOSTIC`` - run diagonostics on line for each solve iteration

Solver Options
--------------
Solver options are applied to the entire model domain.

 * ``HELP`` - prints a list of options on the command line when MAP++ initializes
 * ``INNER_FTOL`` - inner loop function tolerance
 * ``INNER_GTOL`` - desired orthogonality between the function evaluations and Jacobian column
 * ``INNER_XTOL`` - inner loop consecutive iterate tolerance
 * ``INNER_MAX_ITS`` - maximum inner loop iterations 
 * ``OUTER_MAX_ITS`` - maximum outer loop iterations
 * ``OUTER_TOL`` - outer loop tolerance
 * ``OUTER_EPSILON`` - Not used
 * ``INTEGRATION_DT`` - Not used
 * ``KB_DEFAULT`` - Not used
 * ``CB_DEFAULT`` - Not used
 * ``OUTER_CD`` - central difference Jacobian (outer loop solve only)
 * ``OUTER_BD`` - backward difference Jacobian (outer loop solve only)
 * ``OUTER_FD`` - forward difference Jacobian (outer loop solve only)
 * ``LM_MODEL`` - Not used
 * ``PG_COOKED`` - use the relaxation algorithm developed in :cite:`peyrot1979`
 * ``KRYLOV_ACCELERATOR`` - use the Krylov accelerator algorithm developed in :cite:`scott2009`
 * ``REPEAT`` - repeat the element/nodes defined in the input file by mirroring the mooring pattern with a rotation about the Z-axis
 * ``REF_POSITION`` - reference position

Default Options
---------------

.. include-comment:: ../../src/mapinit.c
