//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//
// CAK O-star wind simulation in spherical polar coodinates.
//
// TODO:
//   - Finish initialising magnetic fields -> include B-fiels in x1 BC energy.
//   - Sort out polar BCs:
//     - the idexing is off for the poler BCs and r direction, it maybe that the 
//       code taken from the field_loop_polar.cpp problem is specific to that problem.
//       -> Ask on github site.
//
//======================================================================================== 

// Athena++ headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../parameter_input.hpp"
#include "../coordinates/coordinates.hpp"
#include "../mesh/mesh.hpp"
#include "../eos/eos.hpp"
#include "../hydro/hydro.hpp"
#include "../globals.hpp"
#include "../mesh/mesh_refinement.hpp"
#include <iostream>
#include <sstream>   // stringstream
#include <stdexcept> // runtime_error
#include <string>    // c_str()
#include <algorithm> // std::min
#include <cmath>
#include "../bvals/bvals.hpp"
#include "../field/field.hpp"

int RefinementCondition(MeshBlock *pmb);
void MyBoundary_ix1(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim, FaceField &b,
                    Real time, Real dt, int is, int ie, int js, int je, int ks, int ke);
void AccelRotation(MeshBlock *pmb, const Real time, const Real dt, const AthenaArray<Real> &prim, 
              const AthenaArray<Real> &bcc, AthenaArray<Real> &cons);
void AccelGravity(MeshBlock *pmb, const Real time, const Real dt, const AthenaArray<Real> &prim, 
              const AthenaArray<Real> &bcc, AthenaArray<Real> &cons);
void AccelCAK(MeshBlock *pmb, const Real time, const Real dt, const AthenaArray<Real> &prim, 
              const AthenaArray<Real> &bcc, AthenaArray<Real> &cons);
void Cooling(MeshBlock *pmb, const Real time, const Real dt, const AthenaArray<Real> &prim, 
              const AthenaArray<Real> &bcc, AthenaArray<Real> &cons);
Real dveldr(Real *x, Real *fn);

//========================================================================================
//! \fn void Mesh::InitUserMeshData(ParameterInput *pin)
//  \brief Function to initialize problem-specific data in Mesh class.  Can also be used
//  to initialize variables which are global to (and therefore can be passed to) other
//  functions in this file.  Called in Mesh constructor.
//========================================================================================

void Mesh::InitUserMeshData(ParameterInput *pin)
{
  if(adaptive==true){
    EnrollUserRefinementCondition(RefinementCondition);
  }
  EnrollUserExplicitSourceFunction(AccelGravity);
  EnrollUserExplicitSourceFunction(AccelRotation);
  EnrollUserExplicitSourceFunction(AccelCAK);
  EnrollUserBoundaryFunction(INNER_X1, MyBoundary_ix1);
  EnrollUserExplicitSourceFunction(Cooling);
  return;
}

//========================================================================================
//! \fn void MeshBlock::ProblemGenerator(ParameterInput *pin)
//  \brief Should be used to set initial conditions.
//========================================================================================

void MeshBlock::ProblemGenerator(ParameterInput *pin)
{

  Real gmma  = peos->GetGamma();
  Real gm1  = gmma - 1;

  Real pi = 3.1415926536;
  Real kb = 1.38064852e-16; //[erg K-1]
  Real mp = 1.67262171e-24; //[g]
  Real Big_G = 6.67259e-8;  //[cm3 g-1 s-2]
  Real Msun = 1.989e+33;    //[g]
  Real Rsun = 6.955e+10;    //[cm]
  Real Lsun = 3.846e+33;    //[erg s-1]

  Real Mratio = pin->GetReal("problem","M_RATIO");
  Real Lratio = pin->GetReal("problem","L_RATIO");
  Real Bcgs = pin->GetReal("problem","B_CGS");
  Real T = pin->GetReal("problem","TT");
  Real mu = pin->GetReal("problem","MU");
  Real a = pin->GetReal("problem","AA");
  Real b = pin->GetReal("problem","b_law");
  Real Q = pin->GetReal("problem","QQ");
  Real a_eff = pin->GetReal("problem","aa_eff");
  Real beta = pin->GetReal("problem","BB");
  Real Cs_p = pin->GetReal("problem","Cs_P");
  Real Rratio = 9.0;

  Real UNIT_DENSITY = 1.0e-12;
  Real UNIT_LENGTH = Rratio*6.955e+10;
  Real UNIT_VELOCITY = 1.0e+5;
  Real UNIT_TIME = UNIT_LENGTH/UNIT_VELOCITY;
  Real UNIT_MASS = UNIT_DENSITY*pow(UNIT_LENGTH,3);
  Real UNIT_G = Big_G*UNIT_DENSITY*SQR(UNIT_TIME);
  Real UNIT_L = pow(UNIT_TIME,-3)*(UNIT_DENSITY*pow(UNIT_LENGTH,5));
  Real UNIT_PRESSURE = 1.0e-02;
  Real UNIT_kB = (kb*pow(UNIT_TIME,2))/(UNIT_DENSITY*pow(UNIT_LENGTH,5));
  Real UNIT_B = sqrt((4.0*pi*UNIT_DENSITY*UNIT_VELOCITY*UNIT_VELOCITY));

  Real M_star = Mratio*Msun/UNIT_MASS;
  Real Rs = Rratio*Rsun/UNIT_LENGTH;
  Real Edd = (2.6e-5*(Lratio)*(1.0/Mratio));
  Real L = (Lratio*Lsun/UNIT_L);
  Real c = 3.0e+5;
  Real M_dot = pow(1.0+a_eff,-(1.0/a_eff)) * a_eff * pow(1.0-a_eff,-1)*
             (L/(c*c))*pow(((Q*Edd)*pow(1.0-Edd,-1)),pow(a_eff,-1)-1.0);
  Real cs = sqrt(2.0*kb*T/mp)/UNIT_VELOCITY;
  Real Bq = Bcgs/UNIT_B;
  Real v_esc = sqrt(2.0*UNIT_G*M_star*(1.0-Edd));
  Real v_inf = v_esc*sqrt((a/(1.0-a)));
  Real vv, rho, prs, gs, r, theta, phi;

  AthenaArray<Real> a1,a2,a3;
  int nx1 = (ie-is)+1 + 2*(NGHOST);
  int nx2 = (je-js)+1 + 2*(NGHOST);
  int nx3 = (ke-ks)+1 + 2*(NGHOST);
  a1.NewAthenaArray(nx3,nx2,nx1);
  a2.NewAthenaArray(nx3,nx2,nx1);
  a3.NewAthenaArray(nx3,nx2,nx1);

  for (int k=ks; k<=ke; k++) {
  for (int j=js; j<=je; j++) {
  for (int i=is; i<=ie; i++) {

    r = pcoord->x1v(i);
    theta = pcoord->x2v(j);
    vv = v_inf*pow(1.0 - (1.0/r),b);

    rho = M_dot/(4.0*pi*vv*r*r);
    prs = cs*cs*rho/gmma;

    phydro->u(IDN,k,j,i) = rho;
    phydro->u(IM1,k,j,i) = rho*vv;
    phydro->u(IM2,k,j,i) = 0.0;
    phydro->u(IM3,k,j,i) = 0.0;
    phydro->u(IEN,k,j,i) = prs/(rho*gm1) + 0.5*rho*SQR(vv);

    if (MAGNETIC_FIELDS_ENABLED){ 
      a1(k,j,i) = 0.0;
      a2(k,j,i) = 0.0;
      a3(k,j,i) = Bq*Rs*pow(r,-2)*sin(theta);
    }

  }}}

  if (MAGNETIC_FIELDS_ENABLED){

    // Initialize interface fields
    AthenaArray<Real> area,len,len_p1;
    area.NewAthenaArray(nx1);
    len.NewAthenaArray(nx1);
    len_p1.NewAthenaArray(nx1);
/*
    // for 1,2,3-D
    for (int k=ks; k<=ke; ++k) {
      // reset loop limits for polar boundary
      int jl=js; int ju=je+1;
      if (block_bcs[INNER_X2] == 5) jl=js+1;
      if (block_bcs[OUTER_X2] == 5) ju=je;
      for (int j=jl; j<=ju; ++j) {
        pcoord->Face2Area(k,j,is,ie,area);
        pcoord->Edge3Length(k,j,is,ie+1,len);
        for (int i=is; i<=ie; ++i) {
          pfield->b.x2f(k,j,i) = -1.0*(len(i+1)*a3(k,j,i+1) - len(i)*a3(k,j,i))/area(i);
    }}}
    for (int k=ks; k<=ke+1; ++k) {
      for (int j=js; j<=je; ++j) {
        pcoord->Face3Area(k,j,is,ie,area);
        pcoord->Edge2Length(k,j,is,ie+1,len);
        for (int i=is; i<=ie; ++i) {
          pfield->b.x3f(k,j,i) = (len(i+1)*a2(k,j,i+1) - len(i)*a2(k,j,i))/area(i);
    }}}
    // for 2D and 3D
    if (block_size.nx2 > 1) {
      for (int k=ks; k<=ke; ++k) {
        for (int j=js; j<=je; ++j) {
          pcoord->Face1Area(k,j,is,ie+1,area);
          pcoord->Edge3Length(k,j  ,is,ie+1,len);
          pcoord->Edge3Length(k,j+1,is,ie+1,len_p1);
          for (int i=is; i<=ie+1; ++i) {
            pfield->b.x1f(k,j,i) = (len_p1(i)*a3(k,j+1,i) - len(i)*a3(k,j,i))/area(i);
      }}}
      for (int k=ks; k<=ke+1; ++k) {
        for (int j=js; j<=je; ++j) {
          pcoord->Face3Area(k,j,is,ie,area);
          pcoord->Edge1Length(k,j  ,is,ie,len);
          pcoord->Edge1Length(k,j+1,is,ie,len_p1);
          for (int i=is; i<=ie; ++i) {
            pfield->b.x3f(k,j,i) -= (len_p1(i)*a1(k,j+1,i) - len(i)*a1(k,j,i))/area(i);
      }}}
    }
    // for 3D only
    if (block_size.nx3 > 1) {
      for (int k=ks; k<=ke; ++k) {
        for (int j=js; j<=je; ++j) {
          pcoord->Face1Area(k,j,is,ie+1,area);
          pcoord->Edge2Length(k  ,j,is,ie+1,len);
          pcoord->Edge2Length(k+1,j,is,ie+1,len_p1);
          for (int i=is; i<=ie+1; ++i) {
            pfield->b.x1f(k,j,i) -= (len_p1(i)*a2(k+1,j,i) - len(i)*a2(k,j,i))/area(i);
      }}}
      for (int k=ks; k<=ke; ++k) {
        // reset loop limits for polar boundary
        int jl=js; int ju=je+1;
        if (block_bcs[INNER_X2] == 5) jl=js+1;
        if (block_bcs[OUTER_X2] == 5) ju=je;
        for (int j=jl; j<=ju; ++j) {
          pcoord->Face2Area(k,j,is,ie,area);
          pcoord->Edge1Length(k  ,j,is,ie,len);
          pcoord->Edge1Length(k+1,j,is,ie,len_p1);
          for (int i=is; i<=ie; ++i) {
            pfield->b.x2f(k,j,i) += (len_p1(i)*a1(k+1,j,i) - len(i)*a1(k,j,i))/area(i);
      }}}
    }
*/
    // initialize interface B
    for (int k=ks; k<=ke; k++) {
    for (int j=js; j<=je; j++) {
    for (int i=is; i<=ie+1; i++) {
      pfield->b.x1f(k,j,i) = (a3(k,j+1,i) - a3(k,j,i))/pcoord->dx2f(j) -
                          (a2(k+1,j,i) - a2(k,j,i))/pcoord->dx3f(k);
    }}}
    for (int k=ks; k<=ke; k++) {
    for (int j=js; j<=je+1; j++) {
    for (int i=is; i<=ie; i++) {
      pfield->b.x2f(k,j,i) = (a1(k+1,j,i) - a1(k,j,i))/pcoord->dx3f(k) -
                          (a3(k,j,i+1) - a3(k,j,i))/pcoord->dx1f(i);
    }}}
    for (int k=ks; k<=ke+1; k++) {
    for (int j=js; j<=je; j++) {
    for (int i=is; i<=ie; i++) {
      pfield->b.x3f(k,j,i) = (a2(k,j,i+1) - a2(k,j,i))/pcoord->dx1f(i) -
                          (a1(k,j+1,i) - a1(k,j,i))/pcoord->dx2f(j);
    }}}

    a1.DeleteAthenaArray();
    a2.DeleteAthenaArray();
    a3.DeleteAthenaArray();
    area.DeleteAthenaArray();
    len.DeleteAthenaArray();
    len_p1.DeleteAthenaArray();

    for (int k=ks; k<=ke; k++) {
    for (int j=js; j<=je; j++) {
    for (int i=is; i<=ie; i++) {
      phydro->u(IEN,k,j,i) +=
              0.5*(SQR(0.5*(pfield->b.x1f(k,j,i+1) + pfield->b.x1f(k,j,i)))
                 + SQR(0.5*(pfield->b.x2f(k,j+1,i) + pfield->b.x2f(k,j,i)))
                 + SQR(0.5*(pfield->b.x3f(k+1,j,i) + pfield->b.x3f(k,j,i))));
    }}}

  }

  return;
}

//========================================================================================
//! \fn void MeshBlock::AccelGravity()
//  \brief Used to calculate acceleration due to gravity.
//========================================================================================

void AccelGravity(MeshBlock *pmb, const Real time, const Real dt, const AthenaArray<Real> &prim,
             const AthenaArray<Real> &bcc, AthenaArray<Real> &cons)
{

  Real Big_G = 6.67259e-8;  //[cm3 g-1 s-2]
  Real Msun = 1.989e+33;    //[g]
  Real Rsun = 6.955e+10;    //[cm]

  Real Mratio = 26.6;
  Real Lratio = 1.15e+5;
  Real Rratio = 9.0;

  Real UNIT_DENSITY = 1.0e-12;
  Real UNIT_LENGTH = Rratio*Rsun;
  Real UNIT_VELOCITY = 1.0e+5;
  Real UNIT_TIME = UNIT_LENGTH/UNIT_VELOCITY;
  Real UNIT_MASS = UNIT_DENSITY*pow(UNIT_LENGTH,3);
  Real UNIT_G = Big_G*UNIT_DENSITY*SQR(UNIT_TIME);

  Real M_star = Mratio*Msun/UNIT_MASS;
  Real Edd = (2.6e-5*(Lratio)*(1.0/Mratio));
  Real r = 0.0, gs = 0.0;

  for (int k=pmb->ks; k<=pmb->ke; ++k) {
  for (int j=pmb->js; j<=pmb->je; ++j) {
  for (int i=pmb->is; i<=pmb->ie; ++i) {

    r = pmb->pcoord->x1v(i);
    gs = -UNIT_G*(M_star - Edd)/r/r;;
    cons(IM1,k,j,i) += prim(IDN,k,j,i)*gs*dt;

  }}}
  return;
}

//========================================================================================
//! \fn void MeshBlock::AccelRotation()
//  \brief used to calculate the acceleration due to the rotating frame of the star.
//========================================================================================

void AccelRotation(MeshBlock *pmb, const Real time, const Real dt, const AthenaArray<Real> &prim,
             const AthenaArray<Real> &bcc, AthenaArray<Real> &cons)
{

  Real Big_G = 6.67259e-8;  //[cm3 g-1 s-2]
  Real Msun = 1.989e+33;    //[g]
  Real Rsun = 6.955e+10;    //[cm]

  Real Mratio = 26.6;
  Real omega = 0.05; 
  Real Rratio = 9.0;

  Real UNIT_DENSITY = 1.0e-12;
  Real UNIT_LENGTH = Rratio*Rsun;
  Real UNIT_VELOCITY = 1.0e+5;
  Real UNIT_MASS = UNIT_DENSITY*pow(UNIT_LENGTH,3);
  Real UNIT_TIME = UNIT_LENGTH/UNIT_VELOCITY;
  Real UNIT_G = Big_G*UNIT_DENSITY*SQR(UNIT_TIME);
  Real M_star = Mratio*Msun/UNIT_MASS;
  Real Rs = Rratio*Rsun/UNIT_LENGTH;
  Real omega_fr = pow(omega,2)*(8.0*UNIT_G*M_star/(27.0*Rs));

  Real r, theta, phi;
  Real v_r, v_theta, v_phi;
  Real F[3];

  Real x, y, z;
  Real vx, vy, vz;
  Real x_ce, y_ce, z_ce;

  for (int k=pmb->ks; k<=pmb->ke; ++k) {
  for (int j=pmb->js; j<=pmb->je; ++j) {
  for (int i=pmb->is; i<=pmb->ie; ++i) {

    r = pmb->pcoord->x1v(i);
    theta = pmb->pcoord->x2v(j);
    phi = pmb->pcoord->x3v(k); 
    v_r = prim(IVX,k,j,i);
    v_theta = prim(IVY,k,j,i);
    v_phi = prim(IVZ,k,j,i);

    // transform to cartesian coordinates.
    x = r*sin(theta)*cos(phi);
    y = r*sin(theta)*sin(phi);
    z = r*cos(theta);

    // get cartesian velocity components.
    vx = v_r*sin(v_theta)*cos(v_phi);
    vy = v_r*sin(v_theta)*sin(v_phi);
    vz = v_r*cos(v_theta);

    // calculate coriolis and centrifugal accelerations.
    x_ce = omega_fr*omega_fr*x + 2.0*omega_fr*vy;
    y_ce = omega_fr*omega_fr*y - 2.0*omega_fr*vx;
    z_ce = 0.0;

    // transform accelerations back to spherical polar coordinates.
    F[0] = x_ce*sin(theta)*cos(phi) + y_ce*sin(theta)*sin(phi) + z_ce*cos(theta);
    F[1] = x_ce*cos(theta)*cos(phi) + y_ce*cos(theta)*sin(phi) - z_ce*sin(theta);
    F[2] = -x_ce*sin(phi) + y_ce*cos(phi);

    cons(IM1,k,j,i) += prim(IDN,k,j,i)*F[0]*dt;
    cons(IM2,k,j,i) += prim(IDN,k,j,i)*F[1]*dt;
    cons(IM3,k,j,i) += prim(IDN,k,j,i)*F[2]*dt;

  }}}
  return;

}

//========================================================================================
//! \fn void MeshBlock::AccelCAK()
//  \brief Used to calculate the accelration due to line driving accourding to mCAK
//         theory. This can be improved by including the non-radial acceleration 
//         i.e. Pitard 2009.
//========================================================================================

void AccelCAK(MeshBlock *pmb, const Real time, const Real dt, const AthenaArray<Real> &prim,
             const AthenaArray<Real> &bcc, AthenaArray<Real> &cons)
{

  Real gmma  = pmb->peos->GetGamma();
  Real gm1  = gmma - 1;

  Real pi = 3.1415926536;
  Real kb = 1.38064852e-16; //[erg K-1]
  Real mp = 1.67262171e-24; //[g]
  Real Big_G = 6.67259e-8;  //[cm3 g-1 s-2]
  Real Msun = 1.989e+33;    //[g]
  Real Rsun = 6.955e+10;    //[cm]
  Real Lsun = 3.846e+33;    //[erg s-1]
  Real c = 3.0e+5;

  Real Cs_p = 4.0;
  Real Mratio = 26.6;
  Real Lratio = 1.15e+5;
  Real Bcgs = 300.0;
  Real T = 36.3e+3;
  Real mu = 1.09;
  Real a = 0.6;
  Real b = 0.8;
  Real Q = 700.0;
  Real a_eff = 0.55;
  Real beta = 0.0;
  Real omega = 0.05; 
  Real Rratio = 9.0;

  Real UNIT_DENSITY = 1.0e-12;
  Real UNIT_LENGTH = Rratio*6.955e+10;
  Real UNIT_VELOCITY = 1.0e+5;
  Real UNIT_TIME = UNIT_LENGTH/UNIT_VELOCITY;
  Real UNIT_MASS = UNIT_DENSITY*pow(UNIT_LENGTH,3);
  Real UNIT_G = Big_G*UNIT_DENSITY*SQR(UNIT_TIME);
  Real UNIT_L = pow(UNIT_TIME,-3)*(UNIT_DENSITY*pow(UNIT_LENGTH,5));
  Real UNIT_PRESSURE = 1.0e-02;
  Real UNIT_kB = (kb*pow(UNIT_TIME,2))/(UNIT_DENSITY*pow(UNIT_LENGTH,5));
  Real UNIT_B = sqrt((4.0*pi*UNIT_DENSITY*UNIT_VELOCITY*UNIT_VELOCITY));

  Real M_star = Mratio*Msun/UNIT_MASS;
  Real Rs = Rratio*Rsun/UNIT_LENGTH;
  Real Edd = (2.6e-5*(Lratio)*(1.0/Mratio));
  Real L = (Lratio*Lsun/UNIT_L);
  Real M_dot = pow(1.0+a_eff,-(1.0/a_eff)) * a_eff * pow(1.0-a_eff,-1)*
             (L/(c*c))*pow(((Q*Edd)*pow(1.0-Edd,-1)),pow(a_eff,-1)-1.0);
  Real cs = sqrt(2.0*kb*T/mp)/UNIT_VELOCITY;
  Real Bq = Bcgs/UNIT_B;
  Real v_esc = sqrt(2.0*UNIT_G*M_star*(1.0-Edd));
  Real v_inf = v_esc*sqrt((a/(1.0-a)));
  Real ke = ((4.0*pi*UNIT_G*M_star*c*Edd)/L); 
  Real A = ((1.0/(1.0-a))*((ke*L*Q)/(4.0*pi*c)));;

  Real r;
  Real dvdr, ddr, gL;
  Real nu2_c, B, sigma, f;
  Real v_r, Temp;

  for (int k=pmb->ks; k<=pmb->ke; ++k) {
  for (int j=pmb->js; j<=pmb->je; ++j) {
  for (int i=pmb->is; i<=pmb->ie; ++i) {

    r = pmb->pcoord->x1v(i);
    v_r = prim(IVX,k,j,i);

    ddr = pmb->pcoord->x1v(i) - pmb->pcoord->x1v(i-1);
    dvdr = fabs((prim(IVX,k,j,i+1) - prim(IVX,k,j,i-1))/(2.0*ddr));

    nu2_c = (1.0-(1.0/(r*r)));
    B = (prim(IDN,k,j,i)*Q*c*ke);
    sigma = (r/fabs(v_r))*(dvdr)-1.0;

    f = ((pow(1.0+sigma,1.0+a)-pow(1.0+sigma*nu2_c,1.0+a))/
        ((1.0+a)*(1.0-nu2_c)*sigma*pow(1.0+sigma,a)));

    Temp = prim(IPR,k,j,i)*mu*(mp/UNIT_MASS)/(prim(IDN,k,j,i)*UNIT_kB);

    if (Temp < 1.0e+6){
      gL = (f*A*pow(r,-2)*pow(dvdr/B,a));
    } else {
      gL = 0.0;
    }

    cons(IM1,k,j,i) += prim(IDN,k,j,i)*gL*dt;

  }}}
  return;
}

//========================================================================================
//! \fn void MeshBlock::MyBoundary_ix1()
//  \brief Bounday fot the stellar surface.
//========================================================================================

void MyBoundary_ix1(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim, FaceField &b,
                    Real time, Real dt, int is, int ie, int js, int je, int ks, int ke)
{

  Real gmma  = pmb->peos->GetGamma();
  Real gm1  = gmma - 1;

  Real pi = 3.1415926536;
  Real kb = 1.38064852e-16; //[erg K-1]
  Real mp = 1.67262171e-24; //[g]
  Real Big_G = 6.67259e-8;  //[cm3 g-1 s-2]
  Real Msun = 1.989e+33;    //[g]
  Real Rsun = 6.955e+10;    //[cm]
  Real Lsun = 3.846e+33;    //[erg s-1]
  Real c = 3.0e+5;

  Real Cs_p = 4.0;
  Real Mratio = 26.6;
  Real Lratio = 1.15e+5;
  Real Bcgs = 300.0;
  Real T = 36.3e+3;
  Real mu = 1.09;
  Real a = 0.6;
  Real B = 0.8;
  Real Q = 700.0;
  Real a_eff = 0.55;
  Real beta = 0.0;
  Real Rratio = 9.0;

  Real UNIT_DENSITY = 1.0e-12;
  Real UNIT_LENGTH = Rratio*6.955e+10;
  Real UNIT_VELOCITY = 1.0e+5;
  Real UNIT_TIME = UNIT_LENGTH/UNIT_VELOCITY;
  Real UNIT_MASS = UNIT_DENSITY*pow(UNIT_LENGTH,3);
  Real UNIT_G = Big_G*UNIT_DENSITY*SQR(UNIT_TIME);
  Real UNIT_L = pow(UNIT_TIME,-3)*(UNIT_DENSITY*pow(UNIT_LENGTH,5));
  Real UNIT_PRESSURE = 1.0e-02;
  Real UNIT_kB = (kb*pow(UNIT_TIME,2))/(UNIT_DENSITY*pow(UNIT_LENGTH,5));
  Real UNIT_B = sqrt((4.0*pi*UNIT_DENSITY*UNIT_VELOCITY*UNIT_VELOCITY));

  Real M_star = Mratio*Msun/UNIT_MASS;
  Real Rs = Rratio*Rsun/UNIT_LENGTH;
  Real Edd = (2.6e-5*(Lratio)*(1.0/Mratio));
  Real L = (Lratio*Lsun/UNIT_L);
  Real M_dot = pow(1.0+a_eff,-(1.0/a_eff)) * a_eff * pow(1.0-a_eff,-1)*
             (L/(c*c))*pow(((Q*Edd)*pow(1.0-Edd,-1)),pow(a_eff,-1)-1.0);
  Real cs = sqrt(2.0*kb*T/mp)/UNIT_VELOCITY;

  Real rho = 0.0;
  Real prs = 0.0;

  for (int k=ks; k<=ke; ++k) {
  for (int j=js; j<=je; ++j) {
  for (int i=0; i<=(NGHOST); ++i) {

    rho = M_dot/(4.0*pi*(cs/Cs_p));
    prs = cs*cs*rho/gmma;

    prim(IDN,k,j,i) = rho;
    prim(IM1,k,j,i) = cs/Cs_p;
    prim(IM2,k,j,i) = 0.0;
    prim(IM3,k,j,i) = 0.0;
    prim(IEN,k,j,i) = prs/(rho*gm1) + 0.5*rho*SQR(cs/Cs_p);

  }}}

  if (MAGNETIC_FIELDS_ENABLED) {
    for (int k=ks; k<=ke; ++k) {
    for (int j=js; j<=je; ++j) {
      for (int i=0; i<=(NGHOST); ++i) {
        b.x1f(k,j,(is-i)) = b.x1f(k,j,is);
    }}}
    for (int k=ks; k<=ke; ++k) {
    for (int j=js; j<=je+1; ++j) {
      for (int i=0; i<=(NGHOST); ++i) {
        b.x2f(k,j,(is-i)) = b.x2f(k,j,is);
    }}}

    for (int k=ks; k<=ke+1; ++k) {
    for (int j=js; j<=je; ++j) {
      for (int i=0; i<=(NGHOST); ++i) {
        b.x3f(k,j,(is-i)) = b.x3f(k,j,is);
    }}}
  }

}

//========================================================================================
//! \fn void MeshBlock::Cooling(MeshBlock *pmb, const Real time, const Real dt, 
//                      const AthenaArray<Real> &prim, const AthenaArray<Real> &bcc, 
//                      AthenaArray<Real> &cons)
//  \brief Cooling source term. This is a simple cooling function that is simply a 
//         place holder for a proper one to be coded up later.
//========================================================================================

void Cooling(MeshBlock *pmb, const Real time, const Real dt, const AthenaArray<Real> &prim,
             const AthenaArray<Real> &bcc, AthenaArray<Real> &cons)
{
  Real g = pmb->peos->GetGamma();
  Real tau = 0.01;
  for (int k=pmb->ks; k<=pmb->ke; ++k) {
    for (int j=pmb->js; j<=pmb->je; ++j) {
      for (int i=pmb->is; i<=pmb->ie; ++i) {
        Real temp = (g-1.0)*prim(IEN,k,j,i)/prim(IDN,k,j,i);
        cons(IEN,k,j,i) -= dt*prim(IDN,k,j,i)*(temp - 10.0)/tau/(g-1.0);
      }
    }
  }
  return;
}

//========================================================================================
//! \fn void MeshBlock::RefinementCondition(MeshBlock *pmb)
//  \brief Used to calculate the refinment criteria for AMR.
//========================================================================================

int RefinementCondition(MeshBlock *pmb)
{
  AthenaArray<Real> &w = pmb->phydro->w;
  Real maxeps=0.0;
  int k=pmb->ks;
  for(int j=pmb->js; j<=pmb->je; j++) {
    for(int j=pmb->js; j<=pmb->je; j++) {
      for(int i=pmb->is; i<=pmb->ie; i++) {
        Real epsr= (std::abs(w(IDN,k,j,i+1)-2.0*w(IDN,k,j,i)+w(IDN,k,j,i-1))
                   +std::abs(w(IDN,k,j+1,i)-2.0*w(IDN,k,j,i)+w(IDN,k,j-1,i))
                   +std::abs(w(IDN,k+1,j,i)-2.0*w(IDN,k,j,i)+w(IDN,k-1,j,i)))/w(IDN,k,j,i);
        Real epsp= (std::abs(w(IEN,k,j,i+1)-2.0*w(IEN,k,j,i)+w(IEN,k,j,i-1))
                   +std::abs(w(IEN,k,j+1,i)-2.0*w(IEN,k,j,i)+w(IEN,k,j-1,i))
                   +std::abs(w(IEN,k+1,j,i)-2.0*w(IEN,k,j,i)+w(IEN,k-1,j,i)))/w(IEN,k,j,i);
        Real eps = std::max(epsr, epsp);
        maxeps = std::max(maxeps, eps);
  }}}
  if(maxeps > 0.075) return 1;
  //if(maxeps < 2.8) return -1;
  return 0;
} 

//========================================================================================
//! \fn void MeshBlock::dveldr(Real *x, Real *fn)
//  \brief Function for calculating the Forberg's algorithm:
//         http://faculty.washington.edu/rjl/fdmbook/matlab/fdcoeffF.m
//         This is a matlab script but has been translated to c++.
//========================================================================================

Real dveldr(Real *x, Real *fn)
{

  int n = 5, m = 1, mn = 0;
  int i, j, s, i1, j1, s1;

  Real c0, c1, c2, c3, c4;
  Real xbar = x[2];
  Real C[n-1][m+1];
  Real dvdr = 0.0;

  for (i=0; i<n; i++){
    for (j=0; j<m+1; j++){
      C[i][j] = 0.0;
  }}

  C[0][0] = 1.0;
  c0 = 1;
  c3 = x[0] - xbar;

  for (i=0; i<n; i++){
    i1 = i+1;
    mn = std::min(i,m);
    c1 = 1.0;
    c4 = c3;
    c3 = x[i1] - xbar;
    for (j=0; j<=i-1; j++){
      j1 = j+1;
      c2 = x[i1] - x[j1];
      c1 = c1*c2;
      if (j==i-1) {
        for (s=mn; s>=1; s--){
          s1 = s+1;
          C[i1][s1] = c0*((s+1)*C[i1-1][s1-1] - c4*C[i1-1][s1])/c1;
        }
        C[i1][0] = -c0*c4*C[i1-1][0]/c1;
      }
      for (s=mn; s>=0; s--){
        s1 = s+1;
        C[j1][s1] = (c3*C[j1][s1] - s*C[j1][s1-1])/c2;
      }
      C[j1][0] = c3*C[j1][0]/c2;
    }
    c0 = c1;
  }
  for (i=0; i<n; i++){
    dvdr = C[i][m+1]*fn[i];
  }
  return dvdr;

}

