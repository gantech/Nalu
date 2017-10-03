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

The time stepping method in Nalu is described in the Fuego theory manual :cite:`FuegoTheoryManual:2016` for the backward Euler time discretization. I present a version of this time stepping scheme adapated to the BDF2 time discretization of the momentum equation described earlier. The Favre-averaged momentum equation is written in integral form as

.. math::
   :label: fav-mom-nalu
           
   \int \left . \frac{\partial \rho u_i}{\partial t} \right |^{n+1} {\rm d}V &=  \frac{ (\gamma_1 \rho^{n+1} {u_i}^{n+1} + \gamma_2 \rho^n {u_i}^{n} + \gamma_3 \rho^n {u_i}^{n-1})}{\Delta t} \Delta V \\
   &= {\bf F}_i (\rho^{n+1} u_i^{n+1})  - \int \bar{P}^{n+1} n_i {\rm d}S - \int \left(\rho^{n+1} - \rho_{\circ} \right) g_i {\rm d}V

where

.. math::
   :label: fav-mom-nalu-f
           
   {\bf F}_i (\rho^{n+1} u_i^{n+1}) &= \int \rho^{n+1} u_i^{n+1} u_j^{n+1} n_j {\rm d}S  + \int \bar{\tau}_{ij}^{n+1} n_j {\rm d}S + \int \tau_{u_i u_j}^{n+1} n_j {\rm d}S \\
   &= \int \rho^{n+1} u_i^{n+1} \dot{m}^{n+1}  + \int \bar{\tau}_{ij}^{n+1} n_j {\rm d}S + \int \tau_{u_i u_j}^{n+1} n_j {\rm d}S \\
   
   
The following conventions are used:

.. math::

   \phi^* &= \textrm{ Predicted value of } \phi \textrm{ at } n+1 \textrm{ time step before linear solve} \\
   \widehat{\phi} = \phi^{**} &= \textrm{ Predicted value of } \phi \textrm{ at } n+1 \textrm{ time step after linear solve}

Nalu uses a predictor for the density :math:`\rho^{n+1} = \rho^*` and the mass flow rate through the sub-control surfaces :math:`\dot{m}^{n+1} = \dot{m}^*`. Nalu then corrects for these quantities through outer iterations and hence retains :math:`\rho` and :math:`\dot{m}` constant through each outer iteration. Hence Nalu uses

.. math::
   
   {\bf F}_i (\rho^{n+1} u_i^{n+1}) \approx {\bf F}_i (\rho^{*} u_i^{n+1}) = \int \rho^{*} u_i^{n+1} \dot{m}^{*}  + \int \bar{\tau}_{ij}^{n+1} n_j {\rm d}S + \int \tau_{u_i u_j}^{n+1} n_j {\rm d}S

and solves the following linearized momentum equation.

.. math::
   
   \int \left . \frac{\partial \rho u_i}{\partial t} \right |^{n+1} {\rm d}V \approx {\bf F}_i (\rho^{*} u_i^{n+1}) - \int \bar{P}^{n+1} n_i {\rm d}S - \int \left(\rho^{*} - \rho_{\circ} \right) g_i {\rm d}V

Nalu uses a predictor-corrector method to calculate :math:`u_i^{n+1}` and :math:`P^{n+1}`. First, a momentum predictor step is used to estimate :math:`u_i^{**}` by solving

.. math::
   
   &\frac{ (\gamma_1 \rho^{*} {u_i}^{**} + \gamma_2 \rho^n {u_i}^{n} + \gamma_3 \rho^n {u_i}^{n-1})}{\Delta t} \Delta V \\
   &= {\bf F}_i (\rho^{*} u_i^{**}) - \int P^{*} n_i {\rm d}S - \int \left(\rho^{*} - \rho_{\circ} \right) g_i {\rm d}V + - \int (P^{**} - P^{*}) n_i {\rm d}S, \\
   &= {\bf F}_i (\rho^{*} u_i^{**}) - \int P^{*} n_i {\rm d}S - \int \left(\rho^{*} - \rho_{\circ} \right) g_i {\rm d}V + \epsilon,

where :math:`\epsilon` is an error that reduces with increasing number of outer iterations. :math:`u_i^{**}` will not satisfy the continuity equation. A correction step is performed later to make :math:`u_i^{n+1}` satisfy the continuity equation. :math:`{\bf F} (\rho^{*} u_i^{**})` is linear in :math:`u_i` and hence

.. math::
   :label: linearize-f-phi-star
           
   {\bf F}_i (\rho^{*} u_i^{**}) = \frac{\partial F_i}{\partial u_j} u_j^{**}


Applying Eq. :eq:`linearize-f-phi-star` to Eq. :eq:`fav-mom-nalu`, we get the linearized momentum equation solved in Nalu.
   
.. math::   
   :label: fav-mom-nalu-linearize-f
      
   & \frac{ (\gamma_1 \rho^{*} {u_i}^{**} + \gamma_2 \rho^n {u_i}^{n} + \gamma_3 \rho^{n-1} {u_i}^{n-1})}{\Delta t} \Delta V = \\
   & \quad \quad \frac{\partial F_i}{\partial u_j} u_j^{**} - \int P^{*} n_i {\rm d}S - \int \left(\rho^{*} - \rho_{\circ} \right) g_i {\rm d}V  \\
   & \frac{ (\gamma_1 \rho^{*} {u_i}^{**} - \gamma_1 \rho^{*} {u_i}^{*}) }{\Delta t} \Delta V - \frac{\partial F_i}{\partial u_j} \left ( u_j^{**} - u_j^{*} \right ) = \\
   & \quad \quad \frac{\partial F_i}{\partial u_j} u_j^{*} - \int P^{*} n_i {\rm d}S - \int \left(\rho^{*} - \rho_{\circ} \right) g_i {\rm d}V  \\
   & \quad \quad  - \frac{ (\gamma_1 \rho^{*} {u_i}^{*} + \gamma_2 \rho^{n} {u_i}^{n} + \gamma_3 \rho^{n-1} {u_i}^{n-1})}{\Delta t} \Delta V \\
   & \left ( \frac{ \gamma_1 \rho^{*}}{\Delta t} \Delta V \delta_{ij} - \left . \frac{\partial F_i}{\partial u_j} \right |^{u_i^{*}} \right ) \left ( u_j^{**} - u_j^{*} \right ) = \left (\frac{ \gamma_1 \rho^{*}}{\Delta t} \Delta V \delta_{ij} - \frac{\partial F_i}{\partial u_j} \right ) {u_j}^{*}  \\
   & \quad \quad - \int P^{*} n_i {\rm d}S - \int \left(\rho^{*} - \rho_{\circ} \right) g_i {\rm d}V - \frac{ (\gamma_2 \rho^{n} {u_i}^{n} + \gamma_3 \rho^{n-1} {u_i}^{n-1})}{\Delta t} \Delta V \\
   & \quad {\bf A}_{ij} \; \delta u_j = {\bf A}_{ij} \; u_j^{*} - \int P^{*} n_i {\rm d}S - \int \left(\rho^{*} - \rho_{\circ} \right) g_i {\rm d}V \\
   & \quad \quad \quad - \frac{ (\gamma_2 \rho^{n} {u_i}^{n} + \gamma_3 \rho^{n-1} {u_i}^{n-1})}{\Delta t} \Delta V


The correction to make :math:`u^{**}` satisfy the continuity equation is

.. math::

   u_i^{n+1} = u_i^{**} - \tau {\bf G} \Delta P^{**}
   
.. math::
   :label: 

   u^{n+1} &= u^{**} - \tau_3 {\bf G} \Delta p^{n+1} \\
   0 = {\bf D}(u^{n+1}) &= {\bf D} (\widehat{u}) - {\bf D }( \tau_3 {\bf G} \Delta p^{n+1} ) \\
   0 &= {\bf D} (\widehat{u}) - {\bf D }( \tau_3 {\bf G} p^{n+1} ) + {\bf D }( \tau_3 {\bf G} p^{n} ) \\
   -{\bf L_1} \Delta p^{n+1} &= {\bf D} (\widehat{u} + \tau_3 {\bf G} p^{n}) - {\bf D }( \tau_3 {\bf G} p^{n+1} ) - {\bf L_1} \Delta p^{n+1} \\
   -{\bf L_1} \Delta p^{n+1} &= {\bf D} (\widehat{u} + \tau_3 {\bf G} p^{n}) - {\bf D }( \tau_3 {\bf G} p^{n+1} ) - {\bf L_1} p^{n+1} + {\bf L_1} p^{n}  
   
Domino's :cite:`Domino:2006` generalized sequence for the incremental pressure projection method with stabilization is (with the change from :math:`p^{n+1/2}` and :math:`p^{n-1/2}` to :math:`p^{n+1}` and :math:`p^{n}`

.. math::
   :label: pres-proj-sequence

   {\bf A} \Delta \widehat{u} &= f - {\bf G} p^{n} - {\bf A} u^n \\
   -{\bf L_1} \Delta p^{n+1} &= D \left ( \widehat{u} + \tau_2 {\bf G} p^{n} \right ) + {\bf L_2} p^{n} + b \\
   u^{n+1} &= \widehat{u} - \tau_3 {\bf G} \Delta p^{n+1}
           
   
where

.. math::

   \Delta \widehat{u} &= \widehat{u} - u^n \\
   \Delta p^{n+1} &= p^{n+1} - p^{n}

and :math:`L_1` and :math:`L_2` are Laplacian operators such that

.. math::

   L_1 \phi &= \tau_1 \nabla \phi \cdot d {\bf A} \\
   L_2 \phi &= \tau_2 \nabla \phi \cdot d {\bf A}   


Expanding Eq. :eq:`pres-proj-sequence`,

.. math::

   
   {\bf A} (\widehat{u} - u^n) &= f - {\bf G} p^{n} - {\bf A} u^n \\
   -{\bf L_1} (p^{n+1}-p^{n}) &= -{\bf D} \left ( \widehat{u} + \tau_2 {\bf G} p^{n} \right ) + {\bf L_2} p^{n} + b \\
   {\bf A} u^{n+1} &= {\bf A} \widehat{u} - \tau_3 {\bf G} \Delta p^{n+1} \\
   & \cdots \\
   {\bf A} \left (u^{n+1} + \tau_3 {\bf G} \Delta p^{n+1} \right ) &= f - {\bf G} p^{n} \\
   -{\bf L_1} \Delta p^{n+1} &= -{\bf D} \left ( u^{n+1} + \tau_3 {\bf G} \Delta p^{n+1} + \tau_2 {\bf G} p^{n} \right ) + {\bf L_2} p^{n} + b \\
   & \cdots \\
   {\bf A} u^{n+1} + {\bf G} p^{n+1} + {\bf A} \tau_3 {\bf G} \Delta p^{n+1} &= f - {\bf G} p^{n} + {\bf G} p^{n+1}  \\
   {\bf D} u^{n+1} &= {\bf L_1} \Delta p^{n+1} - {\bf D} \tau_3 {\bf G} \Delta p^{n+1} - {\bf D} \tau_2 {\bf G} p^{n} + {\bf L_2} p^{n} + b \\
   & \cdots \\
   {\bf A} u^{n+1} + {\bf G} p^{n+1} &= ({\bf I}- \tau_3 {\bf A} )  {\bf G} \Delta p^{n+1} \\
   {\bf D} u^{n+1} &= \left ({\bf L_1} - {\bf D} \tau_3 {\bf G} \right ) \Delta p^{n+1} + \left ({\bf L_2} - {\bf D} \tau_2 {\bf G} \right ) \Delta p^{n+1}


Hence the discrete momentum and continuity equations in matrix form with errors becomes

   
.. math::
   :label: mom-continuity-nalu

   \left[
       \begin{array}{lr}
         {\bf A}  &  {\bf G}  \\
         {\bf D}  &  {\bf 0}
       \end{array}
     \right]
   %
     \left[
       \begin{array}{l}
         {\bf u}^{n+1}  \\
         p^{n+1} 
       \end{array}
     \right] =
   %
     \left[
       \begin{array}{l}
         {\bf f}  \\
         0
       \end{array}
     \right]    + 
      \left[
       \begin{array}{l}
         ({\bf I}- \tau_3 {\bf A } ){\bf G}(p^{n+1}-p^{n}) \\ 
         \epsilon({\bf L_i},\tau_i, {\bf D}, {\bf G})
     \end{array}
     \right] 


where the error term that appears for the discrete continuity solve is given by,

.. math::
   :label: contErrorDefined

   \epsilon ({\bf L_i}, \tau_i, {\bf D}, {\bf G}) &= \left ( ({\bf L_1}-{\bf D}\tau_3{\bf G}) \right .\\
   &- \left . ({\bf L_2}-{\bf D}\tau_2{\bf G}) \right ) (p^{n+1}-p^{n})


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
   















