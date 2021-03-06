/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#include <master_element/MasterElement.h>
#include <master_element/MasterElementFunctions.h>
#include <master_element/Quad42DCVFEM.h>

#include <master_element/MasterElementHO.h>
#include <master_element/MasterElementUtils.h>

#include <element_promotion/LagrangeBasis.h>
#include <element_promotion/TensorProductQuadratureRule.h>
#include <element_promotion/QuadratureRule.h>
#include <AlgTraits.h>

#include <NaluEnv.h>
#include <FORTRAN_Proto.h>

#include <stk_util/environment/ReportHandler.hpp>
#include <stk_topology/topology.hpp>

#include <iostream>

#include <cmath>
#include <limits>
#include <array>
#include <map>
#include <memory>

namespace sierra{
namespace nalu{

//-------- quad_derivative -----------------------------------------------------
void quad_derivative(const std::vector<double> &par_coord, 
                     SharedMemView<DoubleType***>& deriv) {
  const double half = 0.5;
  const size_t npts = deriv.dimension(0);

  for(size_t j=0; j<npts; ++j) {
    const DoubleType s1 = par_coord[2*j+0];
    const DoubleType s2 = par_coord[2*j+1];
// shape function derivative in the s1 direction -
    deriv(j,0,0) = - half + s2;
    deriv(j,1,0) =   half - s2;
    deriv(j,2,0) =   half + s2;
    deriv(j,3,0) = - half - s2;

// shape function derivative in the s2 direction -
    deriv(j,0,1) = - half + s1;
    deriv(j,1,1) = - half - s1;
    deriv(j,2,1) =   half + s1;
    deriv(j,3,1) =   half - s1;
  }
}

//-------- quad_gradient_operator -----------------------------------------------------
template<int nint, int npe>
void quad_gradient_operator(const SharedMemView<DoubleType***>& deriv,
                            const SharedMemView<DoubleType**>&  coords,
                            SharedMemView<DoubleType***>& gradop) {

  for (size_t ki=0; ki<nint; ++ki) {
    DoubleType dx_ds1 = 0.0;
    DoubleType dx_ds2 = 0.0;
    DoubleType dy_ds1 = 0.0;
    DoubleType dy_ds2 = 0.0;
 
// calculate the jacobian at the integration station -
    for (size_t kn=0; kn<npe; ++kn) {
      dx_ds1 += deriv(ki,kn,0)*coords(kn,0);
      dx_ds2 += deriv(ki,kn,1)*coords(kn,0);
      dy_ds1 += deriv(ki,kn,0)*coords(kn,1);
      dy_ds2 += deriv(ki,kn,1)*coords(kn,1);
    }
// calculate the determinate of the jacobian at the integration station -
    const DoubleType det_j = dx_ds1*dy_ds2 - dy_ds1*dx_ds2;

// protect against a negative or small value for the determinate of the 
// jacobian. The value of real_min (set in precision.par) represents 
// the smallest Real value (based upon the precision set for this 
// compilation) which the machine can represent - 
    const DoubleType test = stk::math::if_then_else(det_j > 1.e+6*MEconstants::realmin, det_j, 1.0);
    const DoubleType denom = 1.0/test;

// compute the gradient operators at the integration station -

    const DoubleType ds1_dx =  denom*dy_ds2;
    const DoubleType ds2_dx = -denom*dy_ds1;
    const DoubleType ds1_dy = -denom*dx_ds2;
    const DoubleType ds2_dy =  denom*dx_ds1;

    for (size_t kn=0; kn<npe; ++kn) {
      gradop(ki,kn,0) = deriv(ki,kn,0)*ds1_dx + deriv(ki,kn,1)*ds2_dx;
      gradop(ki,kn,1) = deriv(ki,kn,0)*ds1_dy + deriv(ki,kn,1)*ds2_dy;
    }
  }
}

//--------------------------------------------------------------------------
//-------- constructor -----------------------------------------------------
//--------------------------------------------------------------------------
Quad42DSCV::Quad42DSCV()
  : MasterElement()
{
  nDim_ = 2;
  nodesPerElement_ = 4;
  numIntPoints_ = 4;

  // define ip node mappings
  ipNodeMap_.resize(4);
  ipNodeMap_[0] = 0; ipNodeMap_[1] = 1; ipNodeMap_[2] = 2; ipNodeMap_[3] = 3;

  // standard integration location
  intgLoc_.resize(8);    
  intgLoc_[0]  = -0.25; intgLoc_[1]  = -0.25; 
  intgLoc_[2]  = +0.25; intgLoc_[3]  = -0.25; 
  intgLoc_[4]  = +0.25; intgLoc_[5]  = +0.25; 
  intgLoc_[6]  = -0.25; intgLoc_[7]  = +0.25; 

  // shifted integration location
  intgLocShift_.resize(8);    
  intgLocShift_[0]  = -0.50; intgLocShift_[1]  = -0.50; 
  intgLocShift_[2]  = +0.50; intgLocShift_[3]  = -0.50; 
  intgLocShift_[4]  = +0.50; intgLocShift_[5]  = +0.50; 
  intgLocShift_[6]  = -0.50; intgLocShift_[7]  = +0.50; 
}

//--------------------------------------------------------------------------
//-------- destructor ------------------------------------------------------
//--------------------------------------------------------------------------
Quad42DSCV::~Quad42DSCV()
{
  // does nothing
}

//--------------------------------------------------------------------------
//-------- ipNodeMap -------------------------------------------------------
//--------------------------------------------------------------------------
const int *
Quad42DSCV::ipNodeMap(
  int /*ordinal*/)
{
  // define scv->node mappings
  return &ipNodeMap_[0];
}

//--------------------------------------------------------------------------
//-------- determinant -----------------------------------------------------
//--------------------------------------------------------------------------
void Quad42DSCV::determinant(
  SharedMemView<DoubleType**> &coords,
  SharedMemView<DoubleType*> &vol) {

  const int npe  = nodesPerElement_;
  const int nint = numIntPoints_;

// Gaussian quadrature points within an interval [-.25,+.25]
  const double gpp =  0.144337567;
  const double gpm = -0.144337567;
  const double cvm = -0.25;
  const double cvp =  0.25;

  const double half    = 0.5;  
  const double zero    = 0.0;  
  const double one16th = 0.0625;  

  DoubleType deriv[2][4];
  DoubleType shape_fcn[4];

//   store sub-volume centroids
  const double xi[2][4]  ={{cvm,cvp,cvp,cvm},
                           {cvm,cvm,cvp,cvp}};
  const double xigp[2][4]={{gpm,gpp,gpp,gpm},
                           {gpm,gpm,gpp,gpp}};
  DoubleType dx_ds1 = zero;
  DoubleType dx_ds2 = zero;
  DoubleType dy_ds1 = zero;
  DoubleType dy_ds2 = zero;
// 2d cartesian, no cross-section area
  for (int ki=0; ki<nint; ++ki) {
    vol(ki) = zero;

    for (int kq=0; kq<nint; ++kq) {
      dx_ds1 = zero;
      dx_ds2 = zero;
      dy_ds1 = zero;
      dy_ds2 = zero;

      const double ximod  = xi[0][ki] + xigp[0][kq];
      const double etamod = xi[1][ki] + xigp[1][kq];

      deriv[0][0] = -(half - etamod);
      deriv[0][1] =  (half - etamod);
      deriv[0][2] =  (half + etamod);
      deriv[0][3] = -(half + etamod);

      deriv[1][0] = -(half - ximod);
      deriv[1][1] = -(half + ximod);
      deriv[1][2] =  (half + ximod);
      deriv[1][3] =  (half - ximod);

      shape_fcn[0] = (half - ximod)*(half - etamod);
      shape_fcn[1] = (half + ximod)*(half - etamod);
      shape_fcn[2] = (half + ximod)*(half + etamod);
      shape_fcn[3] = (half - ximod)*(half + etamod);

// calculate the jacobian at the integration station -
      for (int kn=0; kn<npe; ++kn) {
        dx_ds1 += deriv[0][kn]*coords(kn,0);
        dx_ds2 += deriv[1][kn]*coords(kn,0);
        dy_ds1 += deriv[0][kn]*coords(kn,1);
        dy_ds2 += deriv[1][kn]*coords(kn,1);
      }
// calculate the determinate of the jacobian at the integration station -
      const DoubleType det_j = (dx_ds1*dy_ds2 - dy_ds1*dx_ds2);

      vol(ki) +=  det_j*one16th;
    } 
  } 
}

//--------------------------------------------------------------------------
//-------- grad_op ---------------------------------------------------------
//--------------------------------------------------------------------------
void Quad42DSCV::grad_op(
  SharedMemView<DoubleType**>& coords,
  SharedMemView<DoubleType***>& gradop,
  SharedMemView<DoubleType***>& deriv) {

  quad_derivative(intgLoc_, deriv);
  quad_gradient_operator<Traits::numScsIp_, Traits::nodesPerElement_>(deriv, coords, gradop);
}

void Quad42DSCV::determinant(
  const int nelem,
  const double *coords,
  double *volume,
  double *error)
{
  int lerr = 0;

  SIERRA_FORTRAN(quad_scv_det)
    ( &nelem, &nodesPerElement_, &numIntPoints_, coords,
      volume, error, &lerr );
}

//--------------------------------------------------------------------------
//-------- shifted_grad_op -------------------------------------------------
//--------------------------------------------------------------------------
void Quad42DSCV::shifted_grad_op(
  SharedMemView<DoubleType**>& coords,
  SharedMemView<DoubleType***>& gradop,
  SharedMemView<DoubleType***>& deriv) {

  quad_derivative(intgLocShift_, deriv);
  quad_gradient_operator<Traits::numScsIp_, Traits::nodesPerElement_>(deriv, coords, gradop);
}

//--------------------------------------------------------------------------
//-------- shape_fcn -------------------------------------------------------
//--------------------------------------------------------------------------
void
Quad42DSCV::shape_fcn(double *shpfc)
{
  quad_shape_fcn(numIntPoints_, &intgLoc_[0], shpfc);
}

//--------------------------------------------------------------------------
//-------- shifted_shape_fcn -----------------------------------------------
//--------------------------------------------------------------------------
void
Quad42DSCV::shifted_shape_fcn(double *shpfc)
{
  quad_shape_fcn(numIntPoints_, &intgLocShift_[0], shpfc);
}

//--------------------------------------------------------------------------
//-------- quad_shape_fcn ---------------------------------------------------
//--------------------------------------------------------------------------
void
Quad42DSCV::quad_shape_fcn(
  const int  &npts,
  const double *isoParCoord, 
  double *shape_fcn)
{
  for (int j = 0; j < npts; ++j ) {
    const int fourj = 4*j;
    const int k = 2*j;
    const double s1 = isoParCoord[k];
    const double s2 = isoParCoord[k+1];
    shape_fcn[    fourj] = 1.0/4.0 + 0.5*(-s1 - s2 ) + s1*s2;
    shape_fcn[1 + fourj] = 1.0/4.0 + 0.5*( s1 - s2 ) - s1*s2;
    shape_fcn[2 + fourj] = 1.0/4.0 + 0.5*( s1 + s2 ) + s1*s2;
    shape_fcn[3 + fourj] = 1.0/4.0 + 0.5*(-s1 + s2 ) - s1*s2;
  }
}

//--------------------------------------------------------------------------
//-------- constructor -----------------------------------------------------
//--------------------------------------------------------------------------
Quad42DSCS::Quad42DSCS()
  : MasterElement()
{
  nDim_ = 2;
  nodesPerElement_ = 4;
  numIntPoints_ = 4;
  scaleToStandardIsoFac_ = 2.0;

  // define L/R mappings
  lrscv_.resize(8);
  lrscv_[0]  = 0; lrscv_[1]  = 1;
  lrscv_[2]  = 1; lrscv_[3]  = 2;
  lrscv_[4]  = 2; lrscv_[5]  = 3;
  lrscv_[6]  = 0; lrscv_[7]  = 3;
  
  // define opposing node
  oppNode_.resize(8);
  // face 0; nodes 0,1
  oppNode_[0] = 3; oppNode_[1] = 2;
  // face 1; nodes 1,2
  oppNode_[2] = 0; oppNode_[3] = 3;
  // face 2; nodes 2,3
  oppNode_[4] = 1; oppNode_[5] = 0;
  // face 3; nodes 3,0
  oppNode_[6] = 2; oppNode_[7] = 1;

  // define opposing face
  oppFace_.resize(8);
  // face 0
  oppFace_[0]  = 3; oppFace_[1] = 1;
  // face 1
  oppFace_[2]  = 0; oppFace_[3] = 2;
  // face 2
  oppFace_[4]  = 1; oppFace_[5] = 3;  
  // face 3
  oppFace_[6]  = 2; oppFace_[7] = 0;

  // standard integration location
  intgLoc_.resize(8);    
  intgLoc_[0] =  0.00; intgLoc_[1] = -0.25; // surf 1; 1->2
  intgLoc_[2] =  0.25; intgLoc_[3] =  0.00; // surf 2; 2->3
  intgLoc_[4] =  0.00; intgLoc_[5] =  0.25; // surf 3; 3->4
  intgLoc_[6] = -0.25; intgLoc_[7] =  0.00; // surf 3; 1->5

  // shifted
  intgLocShift_.resize(8);
  intgLocShift_[0] =  0.00; intgLocShift_[1] = -0.50;
  intgLocShift_[2] =  0.50; intgLocShift_[3] =  0.00;
  intgLocShift_[4] =  0.00; intgLocShift_[5] =  0.50;
  intgLocShift_[6] = -0.50; intgLocShift_[7] =  0.00;

  // exposed face
  intgExpFace_.resize(16);
  // face 0; scs 0, 1; nodes 0,1
  intgExpFace_[0]  = -0.25; intgExpFace_[1]  = -0.50; 
  intgExpFace_[2]  =  0.25; intgExpFace_[3]  = -0.50;
  // face 1; scs 0, 1; nodes 1,2
  intgExpFace_[4]  =  0.50; intgExpFace_[5]  = -0.25;
  intgExpFace_[6]  =  0.50; intgExpFace_[7]  =  0.25;
  // face 2, surf 0, 1; nodes 2,3
  intgExpFace_[8]  =  0.25; intgExpFace_[9]  =  0.50;
  intgExpFace_[10] = -0.25; intgExpFace_[11] =  0.50;
  // face 3, surf 0, 1; nodes 3,0
  intgExpFace_[12] = -0.50; intgExpFace_[13] =  0.25;
  intgExpFace_[14] = -0.50; intgExpFace_[15] = -0.25; 

  // boundary integration point ip node mapping (ip on an ordinal to local node number)
  ipNodeMap_.resize(8); // 2 ips * 4 faces
  // face 0;
  ipNodeMap_[0] = 0;  ipNodeMap_[1] = 1;  
  // face 1;
  ipNodeMap_[2] = 1;  ipNodeMap_[3] = 2; 
  // face 2; 
  ipNodeMap_[4] = 2;  ipNodeMap_[5] = 3;  
  // face 3;
  ipNodeMap_[6] = 3;  ipNodeMap_[7] = 0; 

  sideNodeOrdinals_ = {
      0, 1,
      1, 2,
      2, 3,
      3, 0
  };

  std::vector<std::vector<double>> nodeLocations =
  {
      {-0.5,-0.5}, {+0.5,-0.5},
      {+0.5,+0.5}, {-0.5,+0.5}
  };
  intgExpFaceShift_.resize(16);
  int index = 0;
  stk::topology topo = stk::topology::QUADRILATERAL_4_2D;
  for (unsigned k = 0; k < topo.num_sides(); ++k) {
    stk::topology side_topo = topo.side_topology(k);
    const int* ordinals = side_node_ordinals(k);
    for (unsigned n = 0; n < side_topo.num_nodes(); ++n) {
      intgExpFaceShift_[2*index + 0] = nodeLocations[ordinals[n]][0];
      intgExpFaceShift_[2*index + 1] = nodeLocations[ordinals[n]][1];
      ++index;
    }
  }
}

//--------------------------------------------------------------------------
//-------- destructor ------------------------------------------------------
//--------------------------------------------------------------------------
Quad42DSCS::~Quad42DSCS()
{
  // does nothing
}

//--------------------------------------------------------------------------
//-------- ipNodeMap -------------------------------------------------------
//--------------------------------------------------------------------------
const int *
Quad42DSCS::ipNodeMap(
  int ordinal)
{
  // define ip->node mappings for each face (ordinal); 
  return &ipNodeMap_[ordinal*2];
}


//--------------------------------------------------------------------------
//-------- side_node_ordinals ----------------------------------------------
//--------------------------------------------------------------------------
const int *
Quad42DSCS::side_node_ordinals(
  int ordinal)
{
  // define face_ordinal->node_ordinal mappings for each face (ordinal);
  return &sideNodeOrdinals_[ordinal*2];
}

//--------------------------------------------------------------------------
//-------- determinant -----------------------------------------------------
//--------------------------------------------------------------------------
void Quad42DSCS::determinant(
  SharedMemView<DoubleType**>& coords,
  SharedMemView<DoubleType**>& areav) {
  const double zero   = 0.0;
  const double one    = 1.0;
  const double half   = 0.5;
  const double one4th = 0.25;
  const int kx = 0;
  const int ky = 1;
  // Cartesian
  const double a1 = one;
  const double a2 = zero;
  const double a3 = zero;

  DoubleType coord_mid_face[2][4];

  // calculate element mid-point coordinates
  const DoubleType x1 = (coords(0,kx) + coords(1,kx) + coords(2,kx) + coords(3,kx)) * one4th;
  const DoubleType y1 = (coords(0,ky) + coords(1,ky) + coords(2,ky) + coords(3,ky)) * one4th;
  // calculate element mid-face coordinates
  coord_mid_face[kx][0] = ( coords(0,kx)+coords(1,kx) )*half;
  coord_mid_face[kx][1] = ( coords(1,kx)+coords(2,kx) )*half;
  coord_mid_face[kx][2] = ( coords(2,kx)+coords(3,kx) )*half;
  coord_mid_face[kx][3] = ( coords(3,kx)+coords(0,kx) )*half;

  coord_mid_face[ky][0] = ( coords(0,ky)+coords(1,ky) )*half;
  coord_mid_face[ky][1] = ( coords(1,ky)+coords(2,ky) )*half;
  coord_mid_face[ky][2] = ( coords(2,ky)+coords(3,ky) )*half;
  coord_mid_face[ky][3] = ( coords(3,ky)+coords(0,ky) )*half;
  // Control surface 1
  {
    const DoubleType x2 = coord_mid_face[kx][0];
    const DoubleType y2 = coord_mid_face[ky][0];
    const DoubleType rr = a1 + a2*(x1+x2) + a3*(y1+y2);
    areav(0,kx) = -(y2 - y1)*rr;
    areav(0,ky) =  (x2 - x1)*rr;
  }
  // Control surface 2
  {
    const DoubleType x2 = coord_mid_face[kx][1];
    const DoubleType y2 = coord_mid_face[ky][1];
    const DoubleType rr = a1 + a2*(x1+x2) + a3*(y1+y2);
    areav(1,kx) = -(y2 - y1)*rr;
    areav(1,ky) =  (x2 - x1)*rr;
  }
  // Control surface 3
  {
    const DoubleType x2 = coord_mid_face[kx][2];
    const DoubleType y2 = coord_mid_face[ky][2];
    const DoubleType rr = a1 + a2*(x1+x2) + a3*(y1+y2);
    areav(2,kx) = -(y2 - y1)*rr;
    areav(2,ky) =  (x2 - x1)*rr;
  }
  // Control surface 4
  {
    const DoubleType x2 = coord_mid_face[kx][3];
    const DoubleType y2 = coord_mid_face[ky][3];
    const DoubleType rr = a1 + a2*(x1+x2) + a3*(y1+y2);
    areav(3,kx) =  (y2 - y1)*rr;
    areav(3,ky) = -(x2 - x1)*rr;
  }
}

void Quad42DSCS::determinant(
  const int nelem,
  const double *coords,
  double *areav,
  double *error)
{
  SIERRA_FORTRAN(quad_scs_det)
    ( &nelem, &nodesPerElement_, &numIntPoints_, coords, areav );

  // all is always well; no error checking
  *error = 0;
}

//--------------------------------------------------------------------------
//-------- grad_op ---------------------------------------------------------
//--------------------------------------------------------------------------
void Quad42DSCS::grad_op(
  SharedMemView<DoubleType**>& coords,
  SharedMemView<DoubleType***>& gradop,
  SharedMemView<DoubleType***>& deriv) {

  quad_derivative(intgLoc_, deriv);
  quad_gradient_operator<Traits::numScsIp_, Traits::nodesPerElement_>(deriv, coords, gradop);
}

void Quad42DSCS::grad_op(
  const int nelem,
  const double *coords,
  double *gradop,
  double *deriv,
  double *det_j,
  double *error)
{
  int lerr = 0;

  SIERRA_FORTRAN(quad_derivative)
    ( &numIntPoints_, &intgLoc_[0], deriv );
  
  SIERRA_FORTRAN(quad_gradient_operator)
    ( &nelem,
      &nodesPerElement_,
      &numIntPoints_,
      deriv,
      coords, gradop, det_j, error, &lerr );
  
  if ( lerr )
    NaluEnv::self().naluOutput() << "sorry, negative Quad42DSCS volume.." << std::endl;
}

//--------------------------------------------------------------------------
//-------- shifted_grad_op -------------------------------------------------
//--------------------------------------------------------------------------
void Quad42DSCS::shifted_grad_op(
  SharedMemView<DoubleType**>& coords,
  SharedMemView<DoubleType***>& gradop,
  SharedMemView<DoubleType***>& deriv) {
  quad_derivative(intgLocShift_, deriv);
  quad_gradient_operator<Traits::numScsIp_, Traits::nodesPerElement_>(deriv, coords, gradop);
}

void Quad42DSCS::shifted_grad_op(
  const int nelem,
  const double *coords,
  double *gradop,
  double *deriv,
  double *det_j,
  double *error)
{
  int lerr = 0;

  SIERRA_FORTRAN(quad_derivative)
    ( &numIntPoints_, &intgLocShift_[0], deriv );

  SIERRA_FORTRAN(quad_gradient_operator)
    ( &nelem,
      &nodesPerElement_,
      &numIntPoints_,
      deriv,
      coords, gradop, det_j, error, &lerr );

  if ( lerr )
    NaluEnv::self().naluOutput() << "sorry, negative Quad42DSCS volume.." << std::endl;
}

//--------------------------------------------------------------------------
//-------- face_grad_op ----------------------------------------------------
//--------------------------------------------------------------------------
void Quad42DSCS::face_grad_op(
  const int nelem,
  const int face_ordinal,
  const double *coords,
  double *gradop,
  double *det_j,
  double *error)
{
  int lerr = 0;
  int npf = 2;
  int ndim = 2;

  const int nface = 1;
  double dpsi[8];

  for ( int n=0; n<nelem; n++ ) {
    
    for ( int k=0; k<npf; k++ ) {
      
      const int row = 4*face_ordinal + k*ndim;
      SIERRA_FORTRAN(quad_derivative)
        ( &nface, &intgExpFace_[row], dpsi );
      
      SIERRA_FORTRAN(quad_gradient_operator)
        ( &nface,
          &nodesPerElement_,
          &nface,
          dpsi,
          &coords[8*n], &gradop[k*nelem*8+n*8], &det_j[npf*n+k], error, &lerr );
      
      if ( lerr )
        NaluEnv::self().naluOutput() << "sorry, issue with face_grad_op.." << std::endl;
      
    }
  }
}

//--------------------------------------------------------------------------
//-------- shifted_face_grad_op --------------------------------------------
//--------------------------------------------------------------------------
void Quad42DSCS::shifted_face_grad_op(
  const int nelem,
  const int face_ordinal,
  const double *coords,
  double *gradop,
  double *det_j,
  double *error)
{
  int lerr = 0;
  int npf = 2;
  int ndim = 2;

  const int nface = 1;
  double dpsi[8];

  for ( int n=0; n<nelem; n++ ) {

    for ( int k=0; k<npf; k++ ) {

      const int row = 4*face_ordinal + k*ndim;
      SIERRA_FORTRAN(quad_derivative)
        ( &nface, &intgExpFaceShift_[row], dpsi );

      SIERRA_FORTRAN(quad_gradient_operator)
        ( &nface,
          &nodesPerElement_,
          &nface,
          dpsi,
          &coords[8*n], &gradop[k*nelem*8+n*8], &det_j[npf*n+k], error, &lerr );

      if ( lerr )
        NaluEnv::self().naluOutput() << "sorry, issue with face_grad_op.." << std::endl;

    }
  }
}

//--------------------------------------------------------------------------
//-------- gij -------------------------------------------------------------
//--------------------------------------------------------------------------
void Quad42DSCS::gij(
  SharedMemView<DoubleType**>& coords,
  SharedMemView<DoubleType***>& gupper,
  SharedMemView<DoubleType***>& glower,
  SharedMemView<DoubleType***>& deriv) {
     
  const int npe  = nodesPerElement_;
  const int nint = numIntPoints_;

  DoubleType dx_ds[2][2], ds_dx[2][2];

  for (int ki=0; ki<nint; ++ki) {
    dx_ds[0][0] = 0.0;
    dx_ds[0][1] = 0.0;
    dx_ds[1][0] = 0.0;
    dx_ds[1][1] = 0.0;
 
// calculate the jacobian at the integration station -
    for (int kn=0; kn<npe; ++kn) {
      dx_ds[0][0] += deriv(ki,kn,0)*coords(kn,0);
      dx_ds[0][1] += deriv(ki,kn,1)*coords(kn,0);
      dx_ds[1][0] += deriv(ki,kn,0)*coords(kn,1);
      dx_ds[1][1] += deriv(ki,kn,1)*coords(kn,1);
    }
// calculate the determinate of the jacobian at the integration station -
    const DoubleType det_j = dx_ds[0][0]*dx_ds[1][1] - dx_ds[1][0]*dx_ds[0][1];

// clip
    const double test = stk::simd::get_data(det_j,0);
    const DoubleType denom = (test <= 1.e6*MEconstants::realmin) ? 1.0 : 1.0/det_j;

// compute the inverse jacobian
    ds_dx[0][0] =  dx_ds[1][1]*denom;
    ds_dx[0][1] = -dx_ds[0][1]*denom;
    ds_dx[1][0] = -dx_ds[1][0]*denom;
    ds_dx[1][1] =  dx_ds[0][0]*denom;
         
    for (int i=0; i<2; ++i) {
      for (int j=0; j<2; ++j) {
        gupper(ki,j,i) = dx_ds[i][0]*dx_ds[j][0]+dx_ds[i][1]*dx_ds[j][1];
        glower(ki,j,i) = ds_dx[0][i]*ds_dx[0][j]+ds_dx[1][i]*ds_dx[1][j];
      } 
    }
  }
}

void Quad42DSCS::gij(
  const double *coords,
  double *gupperij,
  double *glowerij,
  double *deriv)
{
  SIERRA_FORTRAN(twod_gij)
    ( &nodesPerElement_,
      &numIntPoints_,
      deriv,
      coords, gupperij, glowerij);
}

//--------------------------------------------------------------------------
//-------- adjacentNodes ---------------------------------------------------
//--------------------------------------------------------------------------
const int *
Quad42DSCS::adjacentNodes()
{
  // define L/R mappings
  return &lrscv_[0];
}

//--------------------------------------------------------------------------
//-------- opposingNodes --------------------------------------------------
//--------------------------------------------------------------------------
int
Quad42DSCS::opposingNodes(
  const int ordinal,
  const int node)
{
  return oppNode_[ordinal*2+node];
}

//--------------------------------------------------------------------------
//-------- opposingFace --------------------------------------------------
//--------------------------------------------------------------------------
int
Quad42DSCS::opposingFace(
  const int ordinal,
  const int node)
{
  return oppFace_[ordinal*2+node];
}

//--------------------------------------------------------------------------
//-------- shape_fcn -------------------------------------------------------
//--------------------------------------------------------------------------
void
Quad42DSCS::shape_fcn(double *shpfc)
{
  quad_shape_fcn(numIntPoints_, &intgLoc_[0], shpfc);
}

//--------------------------------------------------------------------------
//-------- shifted_shape_fcn -----------------------------------------------
//--------------------------------------------------------------------------
void
Quad42DSCS::shifted_shape_fcn(double *shpfc)
{
  quad_shape_fcn(numIntPoints_, &intgLocShift_[0], shpfc);
}

//--------------------------------------------------------------------------
//-------- quad_shape_fcn ---------------------------------------------------
//--------------------------------------------------------------------------
void
Quad42DSCS::quad_shape_fcn(
  const int  &npts,
  const double *isoParCoord, 
  double *shape_fcn)
{
  for (int j = 0; j < npts; ++j ) {
    const int fourj = 4*j;
    const int k = 2*j;
    const double s1 = isoParCoord[k];
    const double s2 = isoParCoord[k+1];
    shape_fcn[    fourj] = 1.0/4.0 + 0.5*(-s1 - s2 ) + s1*s2;
    shape_fcn[1 + fourj] = 1.0/4.0 + 0.5*( s1 - s2 ) - s1*s2;
    shape_fcn[2 + fourj] = 1.0/4.0 + 0.5*( s1 + s2 ) + s1*s2;
    shape_fcn[3 + fourj] = 1.0/4.0 + 0.5*(-s1 + s2 ) - s1*s2;
  }
}

//--------------------------------------------------------------------------
//-------- isInElement -----------------------------------------------------
//--------------------------------------------------------------------------
double
Quad42DSCS::isInElement(
  const double *elemNodalCoord,
  const double *pointCoord,
  double *isoParCoord )
{
  // square of the desired norm, 1.0e-8
  const double isInElemConverged = 1.0e-16;
  const int maxNonlinearIter = 10;

  // -1:1 isoparametric range
  
  // Translate element so that (x,y) coordinates of the first node are (0,0)
  double x[4] = {0.,
                 elemNodalCoord[1] - elemNodalCoord[0],
                 elemNodalCoord[2] - elemNodalCoord[0],
                 elemNodalCoord[3] - elemNodalCoord[0] };
  double y[4] = {0.,
                 elemNodalCoord[5] - elemNodalCoord[4],
                 elemNodalCoord[6] - elemNodalCoord[4],
                 elemNodalCoord[7] - elemNodalCoord[4] };
  
  // (xp,yp) is the point at which we're searching for (xi,eta)
  // (must translate this also)
  
  double xp = pointCoord[0] - elemNodalCoord[0];
  double yp = pointCoord[1] - elemNodalCoord[4];
  
  // Newton-Raphson iteration for (xi,eta)
  double j[4];
  double f[2];
  double shapefct[4];
  
  double xinew = 0.5;     // initial guess
  double etanew = 0.5;
  
  double xicur = 0.5;
  double etacur = 0.5;
  
  double xidiff[2] = { 1.0, 1.0};
  int i = 0;
  
  bool converged = false;
  
  do {
    xicur = xinew;
    etacur = etanew;

    j[0]=  0.25*(1.00-etacur)*x[1]
	  +0.25*(1.00+etacur)*x[2]
	  -0.25*(1.00+etacur)*x[3];

    j[1]= -0.25*(1.00+xicur)*x[1]
	  +0.25*(1.00+xicur)*x[2]
	  +0.25*(1.00-xicur)*x[3];

    j[2]=  0.25*(1.00-etacur)*y[1]
	  +0.25*(1.00+etacur)*y[2]
	  -0.25*(1.00+etacur)*y[3];

    j[3]= -0.25*(1.00+xicur)*y[1]
	  +0.25*(1.00+xicur)*y[2]
	  +0.25*(1.00-xicur)*y[3];

    double jdet = j[0]*j[3] - j[1]*j[2];

    shapefct[0]=0.25*(1.00-etacur)*(1.00-xicur);
    shapefct[1]=0.25*(1.00-etacur)*(1.00+xicur);
    shapefct[2]=0.25*(1.00+etacur)*(1.00+xicur);
    shapefct[3]=0.25*(1.00+etacur)*(1.00-xicur);


    f[0] = (shapefct[1]*x[1]+shapefct[2]*x[2]+shapefct[3]*x[3]) - xp;
    f[1] = (shapefct[1]*y[1]+shapefct[2]*y[2]+shapefct[3]*y[3]) - yp;

    xinew  = xicur  - ( f[0]*j[3] - f[1]*j[1])/jdet;
    etanew = etacur - (-f[0]*j[2] + f[1]*j[0])/jdet;

    xidiff[0] = xinew  - xicur;
    xidiff[1] = etanew - etacur;
    
    double vectorNorm = xidiff[0]*xidiff[0] + xidiff[1]*xidiff[1];
    converged = (vectorNorm < isInElemConverged);
  }  while ( !converged && (++i < maxNonlinearIter) );

  // set a bad value
  isoParCoord[0] = isoParCoord[1] = 1.0e6;
  double dist = 1.0e6;
  if ( i < maxNonlinearIter ) {
    isoParCoord[0] = xinew;
    isoParCoord[1] = etanew;
    dist = (std::abs(xinew) > std::abs(etanew))
      ? std::abs(xinew) : std::abs(etanew);
  }
  return dist;
}

//--------------------------------------------------------------------------
//-------- interpolatePoint ------------------------------------------------
//--------------------------------------------------------------------------
void
Quad42DSCS::interpolatePoint(
  const int &nComp,
  const double *isoParCoord,
  const double *field,
  double *result )
{
  // -1:1 isoparametric range
  const double xi   = isoParCoord[0];
  const double eta  = isoParCoord[1];

  for ( int i = 0; i < nComp; i++ )
  {
    // Base 'field array' index for ith component
    int b = 4*i;

    result[i] = 0.25 * (
      (1-eta) * (1-xi) * field[b+0] +
      (1-eta) * (1+xi) * field[b+1] +
      (1+eta) * (1+xi) * field[b+2] +
      (1+eta) * (1-xi) * field[b+3] ) ;
  }  
}

//--------------------------------------------------------------------------
//-------- general_shape_fcn -----------------------------------------------
//--------------------------------------------------------------------------
void
Quad42DSCS::general_shape_fcn(
  const int numIp,
  const double *isoParCoord,
  double *shpfc)
{
  // -1:1 isoparametric range
  const double npe = nodesPerElement_;
  for ( int ip = 0; ip < numIp; ++ip ) {
    
    const int rowIpc = 2*ip;
    const int rowSfc = npe*ip;
    
    const double s1 = isoParCoord[rowIpc];
    const double s2 = isoParCoord[rowIpc+1];
    shpfc[rowSfc  ] = 0.25*(1.0-s1)*(1.0-s2);
    shpfc[rowSfc+1] = 0.25*(1.0+s1)*(1.0-s2);
    shpfc[rowSfc+2] = 0.25*(1.0+s1)*(1.0+s2);
    shpfc[rowSfc+3] = 0.25*(1.0-s1)*(1.0+s2);
    
  }
}

//--------------------------------------------------------------------------
//-------- general_face_grad_op --------------------------------------------
//--------------------------------------------------------------------------
void 
Quad42DSCS::general_face_grad_op(
  const int face_ordinal,
  const double *isoParCoord,
  const double *coords,
  double *gradop,
  double *det_j,
  double *error)
{
  int lerr = 0;
  const int nface = 1;

  double dpsi[8];

  SIERRA_FORTRAN(quad_derivative)
    ( &nface, &isoParCoord[0], dpsi );
      
  SIERRA_FORTRAN(quad_gradient_operator)
    ( &nface,
      &nodesPerElement_,
      &nface,
      dpsi,
      &coords[0], &gradop[0], &det_j[0], error, &lerr );
  
  if ( lerr )
    NaluEnv::self().naluOutput() << "Quad42DSCS::general_face_grad_op: issue.." << std::endl;
  
}

//--------------------------------------------------------------------------
//-------- sidePcoords_to_elemPcoords --------------------------------------
//--------------------------------------------------------------------------
void 
Quad42DSCS::sidePcoords_to_elemPcoords(
  const int & side_ordinal,
  const int & npoints,
  const double *side_pcoords,
  double *elem_pcoords)
{
  switch (side_ordinal) {
  case 0:
    for (int i=0; i<npoints; i++) {
      elem_pcoords[i*2+0] = 0.5*side_pcoords[i];
      elem_pcoords[i*2+1] = -0.5;
    }
    break;
  case 1:
    for (int i=0; i<npoints; i++) {
      elem_pcoords[i*2+0] = 0.5;
      elem_pcoords[i*2+1] = 0.5*side_pcoords[i];
    }
    break;
  case 2:
    for (int i=0; i<npoints; i++) {
      elem_pcoords[i*2+0] = -0.5*side_pcoords[i];
      elem_pcoords[i*2+1] = 0.5;
    }
    break;
  case 3:
    for (int i=0; i<npoints; i++) {
      elem_pcoords[i*2+0] = -0.5;
      elem_pcoords[i*2+1] = -0.5*side_pcoords[i];
    }
    break;
  default:
    throw std::runtime_error("Quad42DSCS::sideMap invalid ordinal");
  }
}

} // namespace nalu
} // namespace sierra
