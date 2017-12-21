.. _act_wind_turbine_aerodynamics:

============================================
 Actuator Wind Turbine Aerodynamics Modeling
============================================

Theory
======

Wind turbine rotor, tower, and nacelle aerodynamic effects can be
modeled using actuator representations. Compared to resolving the
geometry of the turbine, actuator modeling alleviates the need for a
complex body-fitted meshes, can relax time step restrictions, and
eliminates the need for turbulence modeling at the turbine surfaces.
This comes at the expense of a loss of fine-scale detail, for example,
the boundary layers of the wind turbine surfaces are not resolved.
However, actuator methods well represent wind turbine wakes in the mid
to far downstream regions where wake interactions are important.

Actuator methods usually fall within the classes of disks, lines,
surface, or some blend between the disk and line (i.e., the swept
actuator line). Most commonly, the force over the actuator is computed,
and then applied as a body-force source term, :math:`f_i`, to the
Favre-filtered momentum equation 

.. math::
   :label: fav-mom-bodyforce

     \int \frac{\partial \bar{\rho} \tilde{u}_i} {\partial t} {\rm d}V
     + \int \bar{\rho} \tilde{u}_i \tilde{u}_j n_j {\rm d}S 
     + \int \bar{P} n_i {\rm d}S = \int \bar{\tau}_{ij} n_j {\rm d}S 
     + \int \tau_{u_i u_j} n_j {\rm d}S  
     + \int \left(\bar{\rho} - \rho_{\circ} \right) g_i {\rm d}V \\
     + \int f_i {\rm d}V,

The body-force term :math:`f_i` is volumetric and is a force per unit
volume. The actuator forces, :math:`F'_i`, are not volumetric. They
exist along lines or on surfaces and are force per unit length or area.
Therefore, a projection function, :math:`g`, is used to project the
actuator forces into the fluid volume as volumetric forces. A simple and
commonly used projection function is a uniform Gaussian as proposed by
S{\o}rensen and Shen :cite:`Sorensen:2002`,

.. math:: g(\vec{r}) = \frac{1}{\pi^{3/2} \epsilon^3} e^{-\left( \left| \vec{r} \right|/\epsilon \right)^2},

where :math:`\vec{r}` is the position vector between the fluid point of
interest to a particular point on the actuator, and :math:`\epsilon` is
the width of the Gaussian, which determines how diluted the body force
become. As an example, for an actuator line extending from :math:`l=0`
to :math:`L`, the body force at point :math:`(x,y,z)` due to the line is
given by

.. math::
   :label: force-integral
           
   f_i(x,y,z) = \int_0^L g\left(\vec{r}\left(l\right)\right) F'_i\left(l\right) \: \textrm{d} l.
   

Here, the projection function’s position vector is a function of
position on the actuator line. The part of the line nearest to the point in
the fluid at :math:`(x,y,z)` has most weight.

The force along an actuator line or over an actuator disk is often
computed using blade element theory, where it is convenient to discretize
the actuator into a set of elements. For example, with the actuator line,
the line is broken into discrete line segments, and the force at the center
of each element, :math:`F_i^k`, is computed. Here, :math:`k` is the actuator
element index. These actuator points are independent of the fluid mesh.
This set of point forces is then projected onto the fluid mesh using any
desired projection function, :math:`g(\vec{r})`, as described above.
This is convenient because the integral given in Equation
:eq:`force-integral` can become the summation

.. math::
   :label: force-summation
           
   f_i(x,y,z) = \sum_{k=0}^N g(\vec{r}^k) F_i^k.
   

This summation well approximates the integral given in Equation
:eq:`force-integral` so long as the ratio of actuator element size to
projection function width :math:`\epsilon` does not exceed a certain threshold.

Design
======

The initial actuator capability implemented in Nalu is focused on the actuator line algorithm. However, the class hierarchy is designed with the potential to add other actuator source terms such as actuator disk, swept actuator line and actuator surface capability in the future. The ``ActuatorLineFAST`` class couples Nalu with the third party library OpenFAST for actuator line simulations of wind turbines. OpenFAST (https://nwtc.nrel.gov/FAST), available from https://github.com/OpenFAST/openfast, is a aero-hydro-servo-elastic tool to model wind turbine developed by the National Renewable Energy Laboratory (NREL). The ``ActuatorLineFAST`` class will help Nalu effectively act as an inflow module to OpenFAST by supplying the velocity field information.

We have tested  actuator line implementation to be reasonably scalable. Actuators require searches and parallel communication of blade element velocities and forces, so our implementation should be scalable. Scalability is affected by the number of actuator turbines, the actuator element  density, and the resolution of the mesh surrounding the actuators (i.e., the number of mesh elements that will receive body force). Further testing on scalability is underway with the demonstration of this capability to simulate the OWEZ wind farm.

The actuator line implementation allows for flexible blades that are not necessarily straight (prebend and sweep). The current implementation requires a fixed time step when coupled to OpenFAST, but allows the time step in Nalu to be an integral multiple of the OpenFAST time step. Initially, a simple time lagged FSI model is used to interface Nalu with the turbine model in OpenFAST:

  + The velocity at time step at time step 'n' is sampled at the actuator points and sent 
    to OpenFAST,
  + OpenFAST advances the turbines upto the next Nalu time step 'n+1',
  + The body forces at the actuator points are converted to the source terms of the momentum 
    equation to advance Nalu to the next time step 'n+1'.
    
We are currently working on advanced FSI algorithms along with verification using an MMS approach.
 
The actuator implementation is flexible enough to incorporate a variety of future wind turbine technology capabilities. For example, it is possible that the nacelle may actively tilt for wake steering. The actuator capability is also able to handle a variety of turbines types within one simulation. The current capability allows the modeling of not only the rotor with actuators, but also the tower. However, an aerodynamic model still needs to be implemented for the nacelle.
