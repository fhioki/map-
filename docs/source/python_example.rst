Python Example
==============

.. literalinclude:: ../../python_driver/sphinx_example.py
   :language: python

Output
------

.. code-block:: bash

   MAP++ (Mooring Analysis Program++) Ver. 1.20.00 Nov-29-2015
   MAP environment properties (set externally)...
       Gravity constant          [m/s^2]  : 9.81
       Sea density               [kg/m^3] : 1025.00
       Water depth               [m]      : 350.00
       Vessel reference position [m]      : 0.00 , 0.00 , 0.00
   
   Here is the linearized stiffness matrix with 0.0 vessel displacement:

   [[  1.99e+04   1.54e-03   1.27e-02  -2.81e-01  -2.00e+05   1.03e-01]
    [  6.15e-03   1.99e+04   1.25e-02   2.00e+05  -1.57e-02  -4.96e-01]
    [ -3.12e-03  -6.42e-03   2.27e+04   1.41e-01   2.71e-01   9.51e-02]
    [  7.40e-03   2.00e+05   3.45e-02   1.70e+08  -1.15e-01  -3.00e+01]
    [ -2.00e+05   2.63e-03   4.88e-01  -4.74e-02   1.70e+08   6.39e-02]
    [  9.96e-04  -5.59e-01   9.52e-04  -5.42e+01  -7.34e-02   1.21e+08]]
   
   Here is the linearized stiffness matrix with 5.00 surge vessel displacement:

   [[  1.96e+04  -1.29e-04   1.17e+03  -1.21e-01  -2.15e+05  -6.88e-02]
    [  4.98e-03   2.07e+04   2.99e-04   1.81e+05  -6.10e-02   1.72e+03]
    [  1.17e+03  -3.82e-03   2.32e+04   2.48e-02  -1.19e+04   1.01e-01]
    [  3.94e-03   2.00e+05   1.84e-03   1.70e+08  -2.45e-01  -2.99e+01]
    [ -2.00e+05   2.96e-03   4.90e-01  -8.84e-02   1.70e+08  -4.57e-02]
    [ -4.52e-03  -5.56e-01   1.52e-03  -5.44e+01  -1.77e-01   1.21e+08]]

   Line 0: H = 651461.23 [N]  V = 1178548.28 [N]
   Line 0: Fx = -651461.22 [N]  Fy = 117.69 [N]  Fz = 1178548.28 [N]
   
