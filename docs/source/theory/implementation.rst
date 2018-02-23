Implementation
--------------

.. warning::

   All information here could be factually incorrect. If you find anything wrong, feel free to put in a pull request.

This page contains my attempt to document the implementation of the algorithms previsouly described in this manual in Nalu.

Nodal Gradient
++++++++++++++

There are two methods to calculate the gradient of a variable at a node, viz., using the Gauss theorem and using a Projected Nodal Gradient (PNG). 

Gauss Theorem
~~~~~~~~~~~~~

.. _cvfem-nalu-element-dualvolume:

.. figure:: images/nalu_element_dualvolume.pdf
   :width: 500px
   
   A control volume centered about a finite-element node in a collection of 2-D quadrilateral elements.

The `AssembleNodalGradElemAlgorithm` class computes the nodal gradient by converting a volume integral to a surface integral using Gauss theorem. The value of :math:`\phi` at the sub-control surface integration points is calculated using shape functions.
   
.. math::

   \int \int \int_V \frac{\partial \phi}{\partial x_j} {\rm d}V = \left . \frac{\partial \phi}{\partial x_j} \right |_{node} \Delta V = \int \int_S \phi n_j {\rm d}A = \displaystyle \sum_{scs} \phi_{ip-scs} dA_j^{scs}


The code loops over every element and adds to the value of the gradient at the corresponding "left" and "right" nodes through `gradQL` and `gradQR` respectively.

.. code-block:: c++

    for ( int ip = 0; ip < numScsIp_; ++ip ) {

      // left and right nodes for this ip
      const int il = lrscv[2*ip];
      const int ir = lrscv[2*ip+1];

      stk::mesh::Entity nodeL = node_rels[il];
      stk::mesh::Entity nodeR = node_rels[ir];

      // pointer to fields to assemble
      double *gradQL = stk::mesh::field_data(dqdx_, nodeL );
      double *gradQR = stk::mesh::field_data(dqdx_, nodeR );

      // interpolate to scs point; operate on saved off ws_field
      double qIp = 0.0;
      const int offSet = ip*nodesPerElement_;
      for ( int ic = 0; ic < nodesPerElement_; ++ic ) {
        qIp += p_shape_function_[offSet+ic]*p_scalarQ[ic];
      }

      // left and right volume
      double inv_volL = 1.0/p_dualVolume[il];
      double inv_volR = 1.0/p_dualVolume[ir];

      // assemble to il/ir
      for ( int j = 0; j < nDim_; ++j ) {
        double fac = qIp*p_scs_areav[ip*nDim_+j];
        gradQL[j] += fac*inv_volL;
        gradQR[j] -= fac*inv_volR;
      }
    }


Projected Nodal Gradient
~~~~~~~~~~~~~~~~~~~~~~~~

The projected nodal gradient method tries to minimize the difference between the gradient evaluated using the Gauss theorem and that using sub control volumes.

.. math::
   \left . \frac{\partial \phi}{\partial x_j} \right |_{{\rm node - gauss}} &= \frac{1}{\Delta V} \int \int_S \phi n_j {\rm d}A =  \frac{1}{\Delta V} \displaystyle \sum_{scs} \phi_{ip-scs} dA_j^{scs} \\
   \left . \frac{\partial \phi}{\partial x_j} \right |_{{\rm node - scv}} &= \frac{1}{\Delta V} \displaystyle \sum_{scv} \Delta V_{scv} \; \left ( \frac{\partial \phi}{\partial x_j} \right )_{ip-scv}

The `ProjectedNodalGradientEquationSystem` class uses the `AssemblePNGElemSolverAlgorithm` as it's `solverAlgDriver_` to assemble a system of equations to solve for :math:`(\partial \phi / \partial x_j)_{node-scv}^{**}` given a predicted value of :math:`(\partial \phi / \partial x_j)_{node-scv}^{*}` as:

.. math::
   \left . \frac{\partial \phi}{\partial x_j} \right |_{{\rm node - scv}}^{**} - \left . \frac{\partial \phi}{\partial x_j} \right |_{{\rm node - scv}}^{*} &= \left . \frac{\partial \phi}{\partial x_j} \right |_{{\rm node - scv}}^{*}  + \left . \frac{\partial \phi}{\partial x_j} \right |_{{\rm node - gauss}} \\
   \displaystyle \sum_{scv} \Delta V_{scv} \; \left ( \frac{\partial \phi}{\partial x_j}^{**} -  \frac{\partial \phi}{\partial x_j}^{*} \right )_{ip-scv} &= - \displaystyle \sum_{scv} \Delta V_{scv} \; \left ( \frac{\partial \phi}{\partial x_j}^{*} \right )_{ip-scv} + \displaystyle \sum_{scs} \phi_{ip-scs} dA_j^{scs}

The code first loops over every sub-control surface integration points of an element and adds the second term to the RHS first. Then, the code loops over the sub-control volume integration points of an element and adds the RHS and LHS terms.

.. code-block:: c++

      // handle scs first; all RHS as if it is a source term
      for ( int ip = 0; ip < numScsIp; ++ip ) {

        const int ipNdim = ip*nDim;

        const int offSetSF = ip*nodesPerElement;

        // left and right nodes for this ip
        const int il = lrscv[2*ip];
        const int ir = lrscv[2*ip+1];

        // save off some offsets
        const int ilNdim = il*nDim;
        const int irNdim = ir*nDim;

        // compute scs point values
        double scalarQIp = 0.0;
        for ( int ic = 0; ic < nodesPerElement; ++ic ) {
          const double r = p_shape_function_scs[offSetSF+ic];
          scalarQIp += r*p_scalarQ[ic];
        }
      
        // add residual for each component i
        for ( int i = 0; i < nDim; ++i ) {  
          const int indexL = ilNdim + i;
          const int indexR = irNdim + i;

          const double axi = p_scs_areav[ipNdim+i];

          // right hand side; L and R
          const double rhsFac = -scalarQIp*axi;
          p_rhs[indexL] -= rhsFac;
          p_rhs[indexR] += rhsFac;  
        }
      }

      // handle scv LHS second
      for ( int ip = 0; ip < numScvIp; ++ip ) {

        // nearest node to ip
        const int nearestNode = ipNodeMap[ip];

        // save off some offsets and sc_volume at this ip
        const int nnNdim = nearestNode*nDim;
        const int offSetSF = ip*nodesPerElement;
        const double scV = p_scv_volume[ip];

        // zero out scv
        for ( int j = 0; j < nDim; ++j )
          p_dqdxScv[j] = 0.0;
        
        for ( int ic = 0; ic < nodesPerElement; ++ic ) {
          const double r = p_shape_function_scv[offSetSF+ic];
          for ( int j = 0; j < nDim; ++j ) {
            p_dqdxScv[j] += r*p_dqdx[ic*nDim+j];
          }
        }
      
        // assemble rhs
        for ( int i = 0; i < nDim; ++i ) {
          p_rhs[nnNdim+i] -= p_dqdxScv[i]*scV;
        }

        // manage LHS
        for ( int ic = 0; ic < nodesPerElement; ++ic ) {
          
          const int icNdim = ic*nDim;
      
          // save off shape function
          const double r = p_shape_function_scv[offSetSF+ic];
      
          const double lhsfac = r*scV;
      
          for ( int i = 0; i < nDim; ++i ) {
            const int indexNN = nnNdim + i;
            const int rowNN = indexNN*nodesPerElement*nDim;
            const int rNNiC_i = rowNN+icNdim+i;
            p_lhs[rNNiC_i] += lhsfac;
          }
        }
      }

      apply_coeff(connected_nodes, scratchIds, scratchVals, rhs, lhs, __FILE__);
   
Time stepping method in Nalu
++++++++++++++++++++++++++++

The time stepping method in Nalu is described in the Fuego theory manual
:cite:`FuegoTheoryManual:2016` for the backward Euler time discretization. The
implementation details of the BDF2 time stepping scheme used in Nalu is
described here. The Navier-Stokes equations are written as 

.. math::
   :label: fav-mom-nalu
           
   {\bf F}_i (\rho^{n+1}, u_i^{n+1}, P^{n+1}) - \int \left . \frac{\partial \rho u_i}{\partial t} \right |^{n+1} {\rm d}V &= 0, \\
   {\bf F}_i (\rho^{n+1}, u_i^{n+1}, P^{n+1}) - \frac{ (\gamma_1 \rho^{n+1} {u_i}^{n+1} + \gamma_2 \rho^n {u_i}^{n} + \gamma_3 \rho^n {u_i}^{n-1})}{\Delta t} \Delta V &=0,

where

.. math::
           
   {\bf F}_i (\rho^{n+1} u_i^{n+1}) &= - \int \rho^{n+1} u_i^{n+1} u_j^{n+1} n_j {\rm d}S  + \int \tau_{ij}^{n+1} n_j {\rm d}S - \int P^{n+1} n_i {\rm d}S - \int \left(\rho^{n+1} - \rho_{\circ} \right) g_i {\rm d}V, \\
   &= - \int u_i^{n+1} \dot{m}^{n+1}  + \int \tau_{ij}^{n+1} n_j {\rm d}S  - \int P^{n+1} n_i {\rm d}S - \int \left(\rho^{n+1} - \rho_{\circ} \right) g_i {\rm d}V. \\
   
   
and :math:`\gamma_i` are factors for BDF2 time discretization scheme (see
:numref:`theory_time_discretization`). As is typical of incompressible flow
solvers, the mass flow rate through the sub-control surfaces is tracked
independent of the velocity to maintain conservation of mass. The following
conventions are used:

.. math::

   \phi^* &= \textrm{ Predicted value of } \phi \textrm{ at } n+1 \textrm{ time step before linear solve} \\
   \widehat{\phi} = \phi^{**} &= \textrm{ Predicted value of } \phi \textrm{ at } n+1 \textrm{ time step after linear solve}


The Newton's method is used along with a linearization procedure to predict a
solution to the Navier-Stokes equations at time step :math:`n+1` as

.. math::
   :label: fav-mom-nalu-newton
           
   \mathbf{A}_{ij} \; \delta u_{j} &= {\bf F}_i^{*} - \frac{ (\gamma_1 \rho^{*} {u_i}^{*} + \gamma_2 \rho^n {u_i}^{n} + \gamma_3 \rho^n {u_i}^{n-1})}{\Delta t} \Delta V, \\
   \textrm{where } \delta u_{j} &= u_i^{**} - u_i^*, \\
   \mathbf{A}_{ij} &= \left ( \frac{ \gamma_1 \rho^{*}}{\Delta t} \Delta V \delta_{ij} - \left . \frac{\partial F_i}{\partial u_j} \right |^{*} \right ), \\
   \textrm{and } {\bf F}_i^{*} &= - \int u_i^* \dot{m}^*  + \int \tau_{ij}^* n_j {\rm d}S  - \int P^* n_i {\rm d}S - \int \left(\rho^* - \rho_{\circ} \right) g_i {\rm d}V.


After each Newton or *outer* iteration, :math:`\phi^{**}` is a better approximation to :math:`\phi^{n+1}` compared to :math:`\phi^*`. :math:`\rho*` and :math:`\dot{m}^*` are retained constant through each outer iteration. :math:`{\bf F} (\rho^{*} u_i^{**})` is linear in :math:`u_i^*` and hence

.. math::
   :label: linearize-f-phi-star
           
   {\bf F}_i^* = \left . \frac{\partial F_i}{\partial u_j} \right |^{*} u_j^{*} - \int P^{*} n_i {\rm d}S - \int \left(\rho^{*} - \rho_{\circ} \right) g_i {\rm d}V 
   
Applying Eq. :eq:`linearize-f-phi-star` to Eq. :eq:`fav-mom-nalu-newton`, we get the
linearized momentum predictor equation solved in Nalu.

.. math::   
   :label: fav-mom-nalu-linearize-f
      
   {\bf A}_{ij} \; \delta u_j &= \left . \frac{\partial F_i}{\partial u_j} \right |^{*} u_j^{*} - \int P^{*} n_i {\rm d}S - \int \left(\rho^{*} - \rho_{\circ} \right) g_i {\rm d}V  \\
   & \quad \quad  - \frac{ (\gamma_1 \rho^{*} {u_i}^{*} + \gamma_2 \rho^{n} {u_i}^{n} + \gamma_3 \rho^{n-1} {u_i}^{n-1})}{\Delta t} \Delta V \\
   {\bf A}_{ij} \; \delta u_j &= \left (\frac{ \gamma_1 \rho^{*}}{\Delta t} \Delta V \delta_{ij} - \left . \frac{\partial F_i}{\partial u_j} \right |^{*} \right ) {u_j}^{*} - \int P^{*} n_i {\rm d}S - \int \left(\rho^{*} - \rho_{\circ} \right) g_i {\rm d}V \\
   & \quad - \frac{ (\gamma_2 \rho^{n} {u_i}^{n} + \gamma_3 \rho^{n-1} {u_i}^{n-1})}{\Delta t} \Delta V  \\
   {\bf A}_{ij} \; \delta u_j & = {\bf A}_{ij} \; u_j^{*} - \int P^{*} n_i {\rm d}S - \int \left(\rho^{*} - \rho_{\circ} \right) g_i {\rm d}V \\
   & \quad - \frac{ (\gamma_2 \rho^{n} {u_i}^{n} + \gamma_3 \rho^{n-1} {u_i}^{n-1})}{\Delta t} \Delta V

:math:`u_i^{**}` will not satisfy the continuity equation. A correction step is
performed at the end of each outer iteration to make :math:`u_i^{**}`
satisfy the continuity equation as 

.. math::

   u_i^{n+1} &= u_i^{**} - \frac{\tau_3}{\rho} {\bf G} \Delta P^{**} \\
   \textrm{where } \Delta P^{**} &= P^{**} - P^*


As described in :numref:`theory_errors_splitting_stabilization`, the continuity
equation to be satisfied along with the splitting and stabilization errors is

.. math::
   :label: eq-continuity
           
   {\bf D } \rho u^{**} = b + \left ({\bf L_1} - {\bf D} \tau_3 {\bf G} \right ) \Delta P^{**} + \left ({\bf L_2} - {\bf D} \tau_2 {\bf G} \right ) P^{*}

where :math:`b` contains any source terms when the velocity field is not
divergence free and the other terms are the errors due to pressure stabilization
as shown by Domino :cite:`Domino:2006`. The final pressure Poisson equation
solved to enforce continuity at each outer iteration is
   
.. math::
   :label: eq-pressure

   u^{n+1} &= u^{**} - \frac{\tau_3}{\rho} {\bf G} \Delta P^{**} \\
   b + \left ({\bf L_1} - {\bf D} \tau_3 {\bf G} \right ) \Delta P^{**} &+ \left ({\bf L_2} - {\bf D} \tau_2 {\bf G} \right ) P^{*} \\
   &= {\bf D}(\rho u^{n+1}) = {\bf D} ( \rho \widehat{u}) - {\bf D }( \tau_3 {\bf G} \Delta P^{**} ) \\
   b + {\bf L_1} \Delta P^{**} &= {\bf D} (\rho \widehat{u}) - \left ({\bf L_2} - {\bf D} \tau_2 {\bf G} \right ) P^{*} \\
   -{\bf L_1} \Delta P^{**} &= {\bf D} \rho \widehat{u} + {\bf D} \tau_2 {\bf G} P^{*} - {\bf L_2} P^{*} \\
   -{\bf L_1} \Delta P^{**} &= - {\bf D} \rho \widehat{u} - {\bf D} \tau_2 {\bf G} P^{*} + {\bf L_2} P^{*} + b
  
Thus, the final set of equations solved at each outer iteration is

.. math::

   {\bf A}_{ij} \; \delta u_j & = {\bf A}_{ij} \; u_j^{*} - \int P^{*} n_i {\rm d}S - \int \left(\rho^{*} - \rho_{\circ} \right) g_i {\rm d}V \\
   & \quad - \frac{ (\gamma_2 \rho^{n} {u_i}^{n} + \gamma_3 \rho^{n-1} {u_i}^{n-1})}{\Delta t} \Delta V \\
   -{\bf L_1} \Delta P^{**} &= - {\bf D} \rho \widehat{u} - {\bf D} \tau_2 {\bf G} P^{*} + {\bf L_2} P^{*} + b \\
   u_i^{n+1} &= u_i^{**} - \frac{\tau_3}{\rho} {\bf G} \Delta P^{**}
  
Implementation in code
++++++++++++++++++++++

The core CFD algorithm is the ``LowMachEquationSystem`` as shown in fig. :num:`low-mach-eq-system`.

.. _low-mach-eq-system:

.. figure:: images/lowMachEqSystem.pdf
   :width: 100%

   ``LowMachEqSystem`` class in Nalu.
   
The implementation of the momentum predictor in Nalu uses :eq:`fav-mom-nalu-linearize-f`.

.. math::
   :label: mom-predictor-nalu
           
   \left ( \frac{ \gamma_1 \rho^{*}}{\Delta t} \Delta V \delta_{ij} - \left . \frac{\partial F_i}{\partial u_j} \right |^{u_i^{*}} \right ) \left ( u_j^{**} - u_j^{*} \right ) &= \left . \frac{\partial F_i}{\partial u_j} \right |^{u_i^{*}} u_j^{*} - \int \bar{P}^{*} n_i {\rm d}S + \int \left(\rho^{*} - \rho_{\circ} \right) g_i {\rm d}V \\
   & \quad - \frac{ (\gamma_1 \rho^{*} {u_i}^{*} + \gamma_2 \rho {u_i}^{n} + \gamma_3 \rho {u_i}^{n-1})}{\Delta t} \Delta V


When the element based solver is used for the momentum predictor, Eq. :eq:`mom-predictor-nalu` is split into the three parts. The ``MomentumMassBDF2NodeSuppAlg::node_execute`` function handles the following part

.. math::
   
   \frac{ \gamma_1 \rho^{*}}{\Delta t} \Delta V \delta_{ij} \left ( u_j^{**} - u_j^{*} \right ) + \cdots = - \int \bar{P}^{*} n_i {\rm d}S - \frac{ (\gamma_1 \rho^{*} {u_i}^{*} + \gamma_2 \rho {u_i}^{n} + \gamma_3 \rho {u_i}^{n-1})}{\Delta t} \Delta V + \cdots
   

.. code-block :: c++

   void MomentumMassBDF2NodeSuppAlg::node_execute( double *lhs, double *rhs,  stk::mesh::Entity node)
   {
    // deal with lumped mass matrix (diagonal matrix)
    const double *uNm1      =  stk::mesh::field_data(*velocityNm1_, node);
    const double *uN        =  stk::mesh::field_data(*velocityN_, node);
    const double *uNp1      =  stk::mesh::field_data(*velocityNp1_, node);
    const double rhoNm1     = *stk::mesh::field_data(*densityNm1_, node);
    const double rhoN       = *stk::mesh::field_data(*densityN_, node);
    const double rhoNp1     = *stk::mesh::field_data(*densityNp1_, node);
    const double dualVolume = *stk::mesh::field_data(*dualNodalVolume_, node);
    const double *dpdx = stk::mesh::field_data(*dpdx_, node);
   
    const double lhsfac = gamma1_*rhoNp1*dualVolume/dt_;
    const int nDim = nDim_;
    for ( int i = 0; i < nDim; ++i ) {
      rhs[i] += -(gamma1_*rhoNp1*uNp1[i] + gamma2_*rhoN*uN[i] + gamma3_*rhoNm1*uNm1[i])*dualVolume/dt_
      - dpdx[i]*dualVolume;
      const int row = i*nDim;
      lhs[row+i] += lhsfac;
    }
   }


``AssembleMomentumElemSolverAlgorithm::execute`` function handles the following part

.. math::

   \left ( - \left . \frac{\partial F_i}{\partial u_j} \right |^{u_i^{*}} \right ) \left ( u_j^{**} - u_j^{*} \right ) + \cdots = \left . \frac{\partial F_i}{\partial u_j} \right |^{u_i^{*}} u_j^{*} + \cdots

From Sec. :ref:`advStab`, the advection term of the Favre averaged momentum equation integrated over a dual control volume surrounding a node is (:math:`u_i` is simply referred to as :math:`u_i`)

.. math::
   :label: adv-term-favre-mom
           
   - {\bf F_i}^{u_i^{*}}_{adv} = \int \rho u_m u_{i_{ip}} A_m = \sum \dot{m} u_{i_{ip}}.


where ``ip`` is a sub control surface integration point. If the velocity at the integration point is expressed as a linear combination of the velocity at the nodes surrounding it,

.. math::

   \dot{m} u_{i_{ip}} = \dot{m} \sum \zeta_k u_{i}^k


Hence, the contribution of the combination of node ``k`` and integration point ``ip`` to the advection part of the Jacobian matrix :math:`\partial F_i/ \partial u_j` is

.. math::
   
   - \left . \frac{\partial F_i}{\partial u_j} \right |^{u_i^{*}}_{adv-ip-k} = \dot{m} \zeta_k \delta_{ij}. 

From Sec. :ref:`advStab` the interpolated velocity at any subcontrol surface integration point ``ip`` using cell Peclet blending :math:`\eta`

.. math::
    :label: ui-ip-cell-peclet
            
    u_{i_{ip}} = \eta u_{i_{upw}} + (1-\eta) u_{i_{gcds}},

where the upwind operator :math:`u_{i_{upw}}` is written as 

.. math::
   :label: ui-upw
           
   \phi_{upw} &= \alpha_{upw}\tilde \phi^L_{upw} + \left(1-\alpha_{upw}\right)\phi_{cds}; \dot m > 0, \nonumber \\
                & \quad \alpha_{upw}\tilde\phi^R_{upw} + \left(1-\alpha_{upw}\right)\phi_{cds}; \dot m < 0. \\
                &= 0.5 \frac{\dot{m} + \left | \dot{m} \right |}{\dot{m}} \alpha_{upw}\tilde \phi^L_{upw} + 0.5 \frac{\dot{m} - \left | \dot{m} \right |}{\dot{m}} \alpha_{upw}\tilde\phi^R_{upw} + \left(1-\alpha_{upw}\right)\phi_{cds}

for :math:`\phi = u_i`. From Sec. :ref:`advStab`, the generalized central difference operator :math:`u_{i_{gcds}}` is

.. math::
   :label: ui-gcds
   
   \phi_{gcds} &= \frac{1}{2} \left(  \hat\phi^L_{upw} + \hat\phi^R_{upw} \right) \\
               &= \frac{1}{2} \left(  \alpha \tilde \phi^L_{upw} + \left(1-\alpha\right) \phi_{cds} + \alpha \tilde \phi^R_{upw} + \left(1-\alpha\right) \phi_{cds} \right) \\
               &= \frac{\alpha}{2} \left( \tilde \phi^L_{upw} + \tilde \phi^R_{upw} \right) + \left(1-\alpha\right) \phi_{cds}
   
for :math:`\phi = u_i`. Plugging eqs. :eq:`ui-gcds` and :eq:`ui-upw` into eq. :eq:`ui-ip-cell-peclet`,

.. math::

   u_{i_{ip}} &= \eta \left ( 0.5 \frac{\dot{m} + \left | \dot{m} \right |}{\dot{m}} \alpha_{upw}\tilde \phi^L_{upw} + 0.5 \frac{\dot{m} + \left | \dot{m} \right |}{\dot{m}} \alpha_{upw}\tilde\phi^R_{upw} \right . \\
   & \left . + \left( 1-\alpha_{upw} \right) \phi_{cds} \right )  \\
   & + (1-\eta) \left ( \frac{\alpha}{2} \left( \tilde \phi^L_{upw} + \tilde \phi^R_{upw} \right) + \left(1-\alpha\right) \phi_{cds} \right ) ; \\
   \dot{m} u_{i_{ip}} & = \eta \left ( 0.5 (\dot{m} + \left | \dot{m} \right |) \alpha_{upw} + \dot{m} (1-\eta) \frac{\alpha}{2} \right ) \tilde \phi^L_{upw} \\
   & + \eta \left ( 0.5 (\dot{m} - \left | \dot{m} \right |) \alpha_{upw} + \dot{m} (1-\eta) \frac{\alpha}{2} \right ) \tilde \phi^R_{upw} \\
   & + \dot{m} \left ( \left(1-\alpha_{upw}\right) + (1-\eta) \left(1-\alpha\right) \right ) \phi_{cds}


for :math:`\phi = u_i`. Although the above equations are written for a dual control volume around a single node, the implementation in the code loops over the elements as


.. code-block:: c++

   for ( stk::mesh::Bucket::size_type k = 0 ; k < length ; ++k ) {
      // get elem
      stk::mesh::Entity elem = b[k];
   }

The implementation of the advection term is as follows.
               
.. code-block:: c++

   for ( int ip = 0; ip < numScsIp; ++ip ) {
     for ( int i = 0; i < nDim; ++i ) {

        // 2nd order central
        const double uiIp = p_uIp[i];
        
        // upwind
        const double uiUpwind = (tmdot > 0) ? alphaUpw*p_uIpL[i] + (om_alphaUpw)*uiIp
        : alphaUpw*p_uIpR[i] + (om_alphaUpw)*uiIp;
        
        // generalized central (2nd and 4th order)
        const double uiHatL = alpha*p_uIpL[i] + om_alpha*uiIp;
        const double uiHatR = alpha*p_uIpR[i] + om_alpha*uiIp;
        const double uiCds = 0.5*(uiHatL + uiHatR);
        
        // total advection; pressure contribution in time term
        const double aflux = tmdot*(pecfac*uiUpwind + om_pecfac*uiCds);
     
        const int indexL = ilNdim + i;
        const int indexR = irNdim + i;
        const int rowL = indexL*nodesPerElement*nDim;
        const int rowR = indexR*nodesPerElement*nDim;
        const int rLiL_i = rowL+ilNdim+i;
        const int rLiR_i = rowL+irNdim+i;
        const int rRiR_i = rowR+irNdim+i;
        const int rRiL_i = rowR+ilNdim+i;
        // right hand side; L and R
        p_rhs[indexL] -= aflux;
        p_rhs[indexR] += aflux;
        // upwind advection (includes 4th); left node
        const double alhsfacL = 0.5*(tmdot+std::abs(tmdot))*pecfac*alphaUpw
        + 0.5*alpha*om_pecfac*tmdot;
          p_lhs[rLiL_i] += alhsfacL;
          p_lhs[rRiL_i] -= alhsfacL;
        // upwind advection (includes 4th); right node
        const double alhsfacR = 0.5*(tmdot-std::abs(tmdot))*pecfac*alphaUpw
        + 0.5*alpha*om_pecfac*tmdot;
          p_lhs[rRiR_i] -= alhsfacR;
          p_lhs[rLiR_i] += alhsfacR;
     }

     for ( int ic = 0; ic < nodesPerElement; ++ic ) {
        const int icNdim = ic*nDim;
        
        // shape function
        const double r = p_shape_function[offSetSF+ic];
        const double lhsfacAdv = r*tmdot*(pecfac*om_alphaUpw + om_pecfac*om_alpha);
        for ( int i = 0; i < nDim; ++i ) {

          const int indexL = ilNdim + i;
          const int indexR = irNdim + i;
          
          const int rowL = indexL*nodesPerElement*nDim;
          const int rowR = indexR*nodesPerElement*nDim;
          
          const int rLiC_i = rowL+icNdim+i;
          const int rRiC_i = rowR+icNdim+i;
          
          // advection operator  lhs; rhs handled above
          // lhs; il then ir
          p_lhs[rLiC_i] += lhsfacAdv;
          p_lhs[rRiC_i] -= lhsfacAdv;
        }
     }
   }

From Sec. :ref:`advStab`, the viscous term term of the Favre averaged momentum equation integrated over a dual control volume surrounding a node is

.. math::
   :label: viscous-term-nalu

   {\bf F_i}^{u_i^{*}}_{visc} &= \int 2 \mu_{eff} \left( \tilde{S}_{ij} - \frac{1}{3} \tilde{S}_{kk} \delta_{ij} \right) n_j {\rm d}S \\
           &= \int \mu_{eff} \left( \frac{\partial u_i}{\partial x_j} + \frac{\partial u_j}{\partial x_i}  - \frac{2}{3} \tilde{S}_{kk} \delta_{ij} \right) n_j {\rm d}S \\
           &= \sum_{ip} \mu_{eff} \left (  \sum_m \frac{\partial N_m}{\partial x_j} {\rm d}A_j u_{i_m} + \sum_m \frac{\partial N_m}{\partial x_i} {\rm d}A_j u_{j_m} \right ) \\
           & \quad - \sum_{ip} \left ( \mu_{eff} \frac{2}{3} \tilde{S}_{kk}  {\rm d}A_i \right )

.. code-block:: c++

   for ( int ip = 0; ip < numScsIp; ++ip ) {
     for ( int ic = 0; ic < nodesPerElement; ++ic ) {
       for ( int i = 0; i < nDim; ++i ) {

          const int indexL = ilNdim + i;
          const int indexR = irNdim + i;
          
          const int rowL = indexL*nodesPerElement*nDim;
          const int rowR = indexR*nodesPerElement*nDim;
          
          const int rLiC_i = rowL+icNdim+i;
          const int rRiC_i = rowR+icNdim+i;
          
          // viscous stress
          const int offSetDnDx = nDim*nodesPerElement*ip + icNdim;
          double lhs_riC_i = 0.0;
          for ( int j = 0; j < nDim; ++j ) {
             const double axj = p_scs_areav[ipNdim+j];
             const double uj = p_velocityNp1[icNdim+j];
             // -mu*dui/dxj*A_j; fixed i over j loop; see below..
             const double lhsfacDiff_i = -muIp*p_dndx[offSetDnDx+j]*axj;
              // lhs; il then ir
              lhs_riC_i += lhsfacDiff_i;

              // -mu*duj/dxi*A_j
              const double lhsfacDiff_j = -muIp*p_dndx[offSetDnDx+i]*axj;
              // lhs; il then ir
              p_lhs[rowL+icNdim+j] += lhsfacDiff_j;
              p_lhs[rowR+icNdim+j] -= lhsfacDiff_j;
              // rhs; il then ir
              p_rhs[indexL] -= lhsfacDiff_j*uj;
              p_rhs[indexR] += lhsfacDiff_j*uj;
           }

           // deal with accumulated lhs and flux for -mu*dui/dxj*Aj
           p_lhs[rLiC_i] += lhs_riC_i;
           p_lhs[rRiC_i] -= lhs_riC_i;
           const double ui = p_velocityNp1[icNdim+i];
           p_rhs[indexL] -= lhs_riC_i*ui;
           p_rhs[indexR] += lhs_riC_i*ui;
          }
       }
     }     
   }
   

The implementation of the continuity equation in Nalu solves equation :eq:`eq-pressure` for the difference in pressure :math:`\Delta p^{n+1}` in the integral form.

.. math::
   :label: eq-continuity-implementation
      
   \int & \left ( -{\bf L_1} \Delta p^{n+1} = -  {\bf D} \rho \widehat{u} - {\bf D} \tau_2 {\bf G} p^{n} + {\bf L_2} p^{n} + b \right ) \; d\textrm{V} \\
   - \int & {\bf L_1} \Delta p^{n+1} \; d\textrm{V} = - \int {\bf D} \rho \widehat{u} \; d\textrm{V} - \int {\bf D} \tau_2 {\bf G} p^n \; d\textrm{V} + \int {\bf L_2} p^{n} \; d\textrm{V} + \int b \; d\textrm{V} \\
   - \int & \tau_1 \frac{\partial^2 \Delta p^{n+1}}{\partial x_j^2}  \; d\textrm{V} = - \int \frac{\partial \rho \widehat{u}_j}{\partial x_j}  \; d\textrm{V} - \int  \frac{\partial}{\partial x_j} \left ( \tau_2  \frac{\partial p^n}{\partial x_j} \right ) \; d\textrm{V} + \int \tau_2 \frac{\partial^2 p^{n}}{\partial x_j^2} \; d\textrm{V} + \int b \; d\textrm{V} \\
   - \int_S & \tau_1 \frac{\partial \Delta p^{n+1}}{\partial x_j}  \; d\textrm{S}_j = - \int_S \rho \widehat{u}_j \; d\textrm{S}_j - \int_S \left ( \tau_2 \frac{\partial p^n}{\partial x_j} \right ) \; d\textrm{S}_j + \int_S \tau_2 \frac{\partial p^{n}}{\partial x_j} \; d\textrm{S}_j + \int b \; d\textrm{V} \\
   - \int_S & \tau_1 \left ( \sum_{k \in Nodes} \frac{\partial N_k}{\partial x_j} \Delta p^{n+1}_k \right ) \; d\textrm{S}_j = - \int_S \left ( \sum_{k \in Nodes} N_k \rho \widehat{u}_{j,k} \right ) \; d\textrm{S}_j - \int_S \left ( \tau_2 \left ( \sum_{k \in Nodes} N_k \left . \frac{\partial p^n}{\partial x_j} \right |_k \right ) \right ) \; d\textrm{S}_j + \int_S \tau_2 \left ( \sum_{k \in Nodes} \frac{\partial N_k}{\partial x_j} p^{n}_k \right ) \; d\textrm{S}_j + \int b \; d\textrm{V} 
     
     
The last line in Eq. :eq:`eq-continuity-implementation` shows how :math:`{\bf D} \tau_2 {\bf G} p^{n}` and :math:`{\bf L_2} p^{n}` are implemented differently. In the implementation below :math:`{\bf L_2} p^{n}` is `p_dpdx` and :math:`{\bf D} \tau_2 {\bf G} p^{n}` is `p_Gpdx`. The code provides two options to evaluate the :math:`{\bf D} \rho \widehat{u}` using the blending factor `interpTogether` as

.. math::
   :label: eq-rho-u-interptogether

   (\rho \widehat{u})_{scs} = interpTogether \left ( \sum_{k \in Nodes} N_k \rho_k \widehat{u}_{j,k} \right ) + (1.0 - interpTogether) \left ( \sum_{k \in Nodes} N_k \rho_k \right ) \left ( \sum_{k \in Nodes} N_k \widehat{u}_{j,k} \right )

The implementation of the assembly of the continuity equation described above can be found in `AssembleContinuityElemSolverAlgorithm.C`.
   
.. code-block:: c++   

      for ( int ip = 0; ip < numScsIp; ++ip ) {

        // left and right nodes for this ip
        const int il = lrscv[2*ip];
        const int ir = lrscv[2*ip+1];

        // corresponding matrix rows
        int rowL = il*nodesPerElement;
        int rowR = ir*nodesPerElement;

        // setup for ip values; sneak in geometry for possible reduced sens
        for ( int j = 0; j < nDim; ++j ) {
          p_uIp[j] = 0.0;
          p_rho_uIp[j] = 0.0;
          p_GpdxIp[j] = 0.0;
          p_dpdxIp[j] = 0.0;
        }
        double rhoIp = 0.0;

        const int offSet = ip*nodesPerElement;
        for ( int ic = 0; ic < nodesPerElement; ++ic ) {

          const double r = p_shape_function[offSet+ic];
          const double nodalPressure = p_pressure[ic];
          const double nodalRho = p_density[ic];

          rhoIp += r*nodalRho;

          double lhsfac = 0.0;
          const int offSetDnDx = nDim*nodesPerElement*ip + ic*nDim;
          for ( int j = 0; j < nDim; ++j ) {
            p_GpdxIp[j] += r*p_Gpdx[nDim*ic+j];
            p_uIp[j] += r*p_vrtm[nDim*ic+j];
            p_rho_uIp[j] += r*nodalRho*p_vrtm[nDim*ic+j];
            p_dpdxIp[j] += p_dndx[offSetDnDx+j]*nodalPressure;
            lhsfac += -p_dndx_lhs[offSetDnDx+j]*p_scs_areav[ip*nDim+j];
          }

          // assemble to lhs; left
          p_lhs[rowL+ic] += lhsfac;

          // assemble to lhs; right
          p_lhs[rowR+ic] -= lhsfac;

        }

        // assemble mdot
        double mdot = 0.0;
        for ( int j = 0; j < nDim; ++j ) {
          mdot += (interpTogether*p_rho_uIp[j] + om_interpTogether*rhoIp*p_uIp[j] 
                   - projTimeScale*(p_dpdxIp[j] - p_GpdxIp[j]))*p_scs_areav[ip*nDim+j];
        }

        // residual; left and right
        p_rhs[il] -= mdot/projTimeScale;
        p_rhs[ir] += mdot/projTimeScale;
      }



