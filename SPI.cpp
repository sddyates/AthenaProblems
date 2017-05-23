//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
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
#include <sstream>    // stringstream
#include <stdexcept>  // runtime_error
#include <string>     // c_str()
#include <cmath>
#include "../bvals/bvals.hpp"
#include "../field/field.hpp"

int RefinementCondition(MeshBlock *pmb);

void MySource(MeshBlock *pmb, const Real time, const Real dt, const AthenaArray<Real> &prim,
              const AthenaArray<Real> &bcc, AthenaArray<Real> &cons);

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

  EnrollUserExplicitSourceFunction(MySource);

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
  Real Mj = 0.0009543*Msun; //[g]
  Real Rj = 0.10045*Rsun;   //[cm]

  Real UNIT_DENSITY = 1.0e-15;
  Real UNIT_LENGTH = 6.955e+10;
  Real UNIT_VELOCITY = 1.0e+5;
  Real UNIT_TIME = UNIT_LENGTH/UNIT_VELOCITY;
  Real UNIT_MASS = UNIT_DENSITY*pow(UNIT_LENGTH,3);
  Real UNIT_G = Big_G*UNIT_DENSITY*SQR(UNIT_TIME);

  Real a = pin->GetReal("problem","sep")*1.49597892e+11/UNIT_LENGTH;

  Real Ms = pin->GetReal("problem","M_star")*Msun/UNIT_MASS;
  Real Rs = pin->GetReal("problem","R_star")*Rsun/UNIT_LENGTH;
  Real Ts = pin->GetReal("problem","T_star");
  Real RHs = pin->GetReal("problem","RH_star")/UNIT_DENSITY;
  Real B0s = pin->GetReal("problem","B0_star");

  Real Mp = pin->GetReal("problem","M_planet")*Mj/UNIT_MASS;
  Real Rp = pin->GetReal("problem","R_planet")*Rj/UNIT_LENGTH;
  Real Tp = pin->GetReal("problem","T_planet");
  Real RHp = pin->GetReal("problem","RH_planet")/UNIT_DENSITY;
  Real B0p = pin->GetReal("problem","B0_planet");

  Real v_escp = sqrt(2.0*UNIT_G*Mp/Rp);
  Real v_escs = sqrt(2.0*UNIT_G*Ms/Rs);
  Real csp = sqrt((2.0*kb*Tp)/mp)/UNIT_VELOCITY;
  Real css = sqrt((2.0*kb*Ts)/mp)/UNIT_VELOCITY;
  Real omega_orb = sqrt(UNIT_G*Ms/pow(a,3));
  Real omegas = omega_orb;      
  Real omegap = omega_orb;
  Real omega_fr = omega_orb;
  Real Pp = csp*csp*RHp/gmma;                                                           
  Real Ps = css*css*RHs/gmma;
  Real rcs = UNIT_G*Ms/(2.0*css*css);
  Real rcp = UNIT_G*Mp/(2.0*csp*csp);
  Real gs = 0.0;
  Real gp = 0.0;
  Real gs_in = 0.0;
  Real gp_in = 0.0;

  Real PRSs = 0.0;
  Real RHOs = 0.0;
  Real PRSp = 0.0;
  Real RHOp = 0.0;

  Real rs, thetas, phis;
  Real rp, thetap, phip;

  Real lambdas = 0.0;
  Real lambdap = 0.0;
  Real etas = 0.0;
  Real etap = 0.0;
  Real b = 0.0;
  Real o = 0.0;
  Real psi = 0.0; 
  Real function = 0.0;
  Real dfunction = 0.0;
  Real step = 0.0;

  Real v_ws = 0.0;
  Real v_ws_x1 = 0.0;
  Real v_ws_x2 = 0.0;
  Real v_ws_x3 = 0.0;
  Real v_wp = 0.0;
  Real v_wp_x1 = 0.0;
  Real v_wp_x2 = 0.0;
  Real v_wp_x3 = 0.0;

  for (int k=ks; k<=ke; k++) {
  for (int j=js; j<=je; j++) {
  for (int i=is; i<=ie; i++) {

    Real x = pcoord->x1v(i);
    Real y = pcoord->x2v(j);
    Real z = pcoord->x3v(k);

    rs = sqrt(SQR(x) + SQR(y) + SQR(z));
    rp = sqrt(SQR(x - a) + SQR(y) + SQR(z));

    thetap = acos(z/rp);
    phip = atan2(y,(x-a));
    thetas = acos(z/rs);
    phis = atan2(y,x);

    // Stellar interior.
    if (rs < 0.5*Rs) {    
      phydro->u(IDN,k,j,i) = RHs;
      phydro->u(IM1,k,j,i) = 0.0;
      phydro->u(IM2,k,j,i) = 0.0;
      phydro->u(IM3,k,j,i) = 0.0;
      phydro->u(IEN,k,j,i) = Ps/(RHs*gm1);
    } else if (rs >= 0.5*Rs && rs < Rs) {
      phydro->u(IDN,k,j,i) = RHs;
      phydro->u(IM1,k,j,i) = 0.0;
      phydro->u(IM2,k,j,i) = 0.0;
      phydro->u(IM3,k,j,i) = 0.0;
      phydro->u(IEN,k,j,i) = Ps/(RHs*gm1);
    } else if (rs >= Rs && rp > 1.5*Rp) {
      // Newton raphson method Star.
      lambdas = 0.5*pow(v_escs/css,2);
      etas = rs/Rs;
      b = -3.0 - 4.0*log(lambdas/2.0) + 4.0*log(etas) + 2.0*(lambdas/etas);
      o = 0;
      if (rs <= rcs) {
        psi = 2.e-8;
      } else {
        psi = 2.5;
      }
      do {
        function  = -psi + log(psi) + b;
        dfunction = -1.0 + (1.0/psi);
        step = function/dfunction;
        psi = psi - step;
        o++;
      } while (fabs(step) > 1.0e-8 && o < 1000);
      v_ws = css*sqrt(psi);
      if(isnan(v_ws)){
        v_ws = 0.0;
      }
      v_ws_x1 = sin(thetas)*(v_ws*cos(phis)+sin(phis)*rs*(omega_fr+omegas));
      v_ws_x2 = sin(thetas)*(v_ws*sin(phis)-cos(phis)*rs*(omega_fr+omegas));
      v_ws_x3 = v_ws*cos(thetas);
      PRSs = Ps*exp(lambdas*(Rs/rs-1.0)-0.5*pow(v_ws/css,2));
      RHOs = (RHs/Ps)*PRSs;

      phydro->u(IDN,k,j,i) = RHOs;
      phydro->u(IM1,k,j,i) = RHOs*v_ws_x1;
      phydro->u(IM2,k,j,i) = RHOs*v_ws_x2;
      phydro->u(IM3,k,j,i) = RHOs*v_ws_x3;
      phydro->u(IEN,k,j,i) = PRSs/(RHOs*gm1) + 0.5*RHOs*(SQR(v_ws_x1) + SQR(v_ws_x2) + SQR(v_ws_x3));


    }

    // Planetary interior.
    if (rp < 0.5*Rp) {
      phydro->u(IDN,k,j,i) = RHp;
      phydro->u(IM1,k,j,i) = 0.0;
      phydro->u(IM2,k,j,i) = 0.0;
      phydro->u(IM3,k,j,i) = 0.0;
      phydro->u(IEN,k,j,i) = Pp/(RHp*gm1);
    } else if (rp >= 0.5*Rp && rp < Rp) {
      phydro->u(IDN,k,j,i) = RHp;
      phydro->u(IM1,k,j,i) = 0.0;
      phydro->u(IM2,k,j,i) = 0.0;
      phydro->u(IM3,k,j,i) = 0.0;
      phydro->u(IEN,k,j,i) = Pp/(RHp*gm1);
    } else if (rp >= Rp && rp <= 1.5*Rp) {
      lambdap = 0.5*pow(v_escp/csp,2);
      etap = rp/Rp;
      b = -3.0 - 4.0*log(lambdap/2.0) + 4.0*log(etap)+2.0*(lambdap/etap);
      o = 0;
      if (rp <= rcp) {
        psi = 2.e-8;
      } else {
        psi = 2.5;
      }
      do {
        function  = -psi + log(psi) + b;
        dfunction = -1.0 + (1.0/psi);
        step=function/dfunction;
        psi=psi-step;
        o++;
      } while (fabs(step) > 1.0e-8 && o < 1000);
      v_wp = csp*sqrt(psi);
      if(isnan(v_wp)){
        v_wp = 0.0;
      }
      v_wp_x1 = sin(thetap)*(v_wp*cos(phip)+sin(phis)*rs*omega_fr-sin(phip)*rp*omegap);
      v_wp_x2 = sin(thetap)*(v_wp*sin(phip)-cos(phis)*rs*omega_fr+a*omega_orb + cos(phip)*rp*omegap);
      v_wp_x3 = v_wp*cos(thetap);
      PRSp = Pp*exp(lambdap*(Rp/rp-1.0)-0.5*pow(v_wp/csp,2));
      RHOp = (RHp/Pp)*PRSp;

      phydro->u(IDN,k,j,i) = RHOp;
      phydro->u(IM1,k,j,i) = RHOp*v_wp_x1;
      phydro->u(IM2,k,j,i) = RHOp*v_wp_x2;
      phydro->u(IM3,k,j,i) = RHOp*v_wp_x3;
      phydro->u(IEN,k,j,i) = PRSp/(RHOp*gm1) + 0.5*RHOp*(SQR(v_wp_x1) + SQR(v_wp_x2) + SQR(v_wp_x3));
    }

  }}}


  return;
}

//========================================================================================
//! \fn void MeshBlock::UserWorkInLoop(void)
//  \brief Function called once every time step for user-defined work.
//========================================================================================

void MeshBlock::UserWorkInLoop(void)
{
  // do nothing
  return;
}


void MySource(MeshBlock *pmb, const Real time, const Real dt, const AthenaArray<Real> &prim,
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
  Real Mj = 0.0009543*Msun; //[g]
  Real Rj = 0.10045*Rsun;   //[cm]

  Real UNIT_DENSITY = 1.0e-15;
  Real UNIT_LENGTH = 6.955e+10;
  Real UNIT_VELOCITY = 1.0e+5;
  Real UNIT_TIME = UNIT_LENGTH/UNIT_VELOCITY;
  Real UNIT_MASS = UNIT_DENSITY*pow(UNIT_LENGTH,3);
  Real UNIT_G = Big_G*UNIT_DENSITY*SQR(UNIT_TIME);

  Real a = 4.7*1.49597892e+11/UNIT_LENGTH;

  Real Ms = 1.0*Msun/UNIT_MASS;
  Real Rs = 1.0*Rsun/UNIT_LENGTH;
  Real Ts = 1.0e+6;
  Real RHs = 5.0e-15/UNIT_DENSITY;
  Real B0s = 2.0;

  Real Mp = 0.5*Mj/UNIT_MASS;
  Real Rp = 1.5*Rj/UNIT_LENGTH;
  Real Tp = 1.0e+4;
  Real RHp = 7.0e-16/UNIT_DENSITY;
  Real B0p = 1.0;

  Real v_escp = sqrt(2.0*UNIT_G*Mp/Rp);
  Real v_escs = sqrt(2.0*UNIT_G*Ms/Rs);
  Real csp = sqrt((2.0*kb*Tp)/mp)/UNIT_VELOCITY;
  Real css = sqrt((2.0*kb*Ts)/mp)/UNIT_VELOCITY;
  Real omega_orb = sqrt(UNIT_G*Ms/pow(a,3));
  Real omegas = omega_orb;      
  Real omegap = omega_orb;
  Real omega_fr = omega_orb;
  Real Pp = csp*csp*RHp/gmma;                                                           
  Real Ps = css*css*RHs/gmma;
  Real rcs = UNIT_G*Ms/(2.0*css*css);
  Real rcp = UNIT_G*Mp/(2.0*csp*csp);
  Real gs = 0.0;
  Real gp = 0.0;
  Real gs_in = 0.0;
  Real gp_in = 0.0;

  Real PRSs = 0.0;
  Real RHOs = 0.0;
  Real PRSp = 0.0;
  Real RHOp = 0.0;

  Real rs, thetas, phis;
  Real rp, thetap, phip;

  Real lambdas = 0.0;
  Real lambdap = 0.0;
  Real etas = 0.0;
  Real etap = 0.0;
  Real b = 0.0;
  Real o = 0.0;
  Real psi = 0.0; 
  Real function = 0.0;
  Real dfunction = 0.0;
  Real step = 0.0;

  Real v_ws = 0.0;
  Real v_ws_x1 = 0.0;
  Real v_ws_x2 = 0.0;
  Real v_ws_x3 = 0.0;
  Real v_wp = 0.0;
  Real v_wp_x1 = 0.0;
  Real v_wp_x2 = 0.0;
  Real v_wp_x3 = 0.0;

  Real x = 0.0;
  Real y = 0.0;
  Real z = 0.0;

  for (int k=pmb->ks; k<=pmb->ke; ++k) {
  for (int j=pmb->js; j<=pmb->je; ++j) {
  for (int i=pmb->is; i<=pmb->ie; ++i) {

    x = pmb->pcoord->x1v(i);
    y = pmb->pcoord->x2v(j);
    z = pmb->pcoord->x3v(k);

    rs = sqrt(SQR(x) + SQR(y) + SQR(z));
    rp = sqrt(SQR(x - a) + SQR(y) + SQR(z));

    thetap = acos(z/rp);
    phip = atan2(y,(x-a));
    thetas = acos(z/rs);
    phis = atan2(y,x);

    // Stellar interior.
    if (rs < 0.5*Rs) {
      cons(IDN,k,j,i) = RHs;
      cons(IM1,k,j,i) = 0.0;
      cons(IM2,k,j,i) = 0.0;
      cons(IM3,k,j,i) = 0.0;
      cons(IEN,k,j,i) = Ps/(RHs*gm1);
    } else if (rs >= 0.5*Rs && rs < Rs) {
      cons(IDN,k,j,i) = RHs;
      cons(IM1,k,j,i) = 0.0;
      cons(IM2,k,j,i) = 0.0;
      cons(IM3,k,j,i) = 0.0;
      cons(IEN,k,j,i) = Ps/(RHs*gm1);
    } else if (rs >= Rs && rs <= 1.5*Rs) {
      // Newton raphson method Star.
      lambdas = 0.5*pow(v_escs/css,2);
      etas = rs/Rs;
      b = -3.0 - 4.0*log(lambdas/2.0) + 4.0*log(etas) + 2.0*(lambdas/etas);
      o = 0;
      if (rs <= rcs) {
        psi = 2.e-8;
      } else {
        psi = 2.5;
      }
      do {
        function  = -psi + log(psi) + b;
        dfunction = -1.0 + (1.0/psi);
        step = function/dfunction;
        psi = psi - step;
        o++;
      } while (fabs(step) > 1.0e-8 && o < 1000);
      v_ws = css*sqrt(psi);
      if(isnan(v_ws)){
        v_ws = 0.0;
      }
      v_ws_x1 = sin(thetas)*(v_ws*cos(phis)+sin(phis)*rs*(omega_fr+omegas));
      v_ws_x2 = sin(thetas)*(v_ws*sin(phis)-cos(phis)*rs*(omega_fr+omegas));
      v_ws_x3 = v_ws*cos(thetas);
      PRSs = Ps*exp(lambdas*(Rs/rs-1.0)-0.5*pow(v_ws/css,2));
      RHOs = (RHs/Ps)*PRSs;

      cons(IDN,k,j,i) = RHOs;
      cons(IM1,k,j,i) = RHOs*v_ws_x1;
      cons(IM2,k,j,i) = RHOs*v_ws_x2;
      cons(IM3,k,j,i) = RHOs*v_ws_x3;
      cons(IEN,k,j,i) = PRSs/(RHOs*gm1) + 0.5*RHOs*(SQR(v_ws_x1) + SQR(v_ws_x2) + SQR(v_ws_x3));
    }

    // Planetary interior.
    if (rp < 0.5*Rp) {
      cons(IDN,k,j,i) = RHp;
      cons(IM1,k,j,i) = 0.0;
      cons(IM2,k,j,i) = 0.0;
      cons(IM3,k,j,i) = 0.0;
      cons(IEN,k,j,i) = Pp/(RHp*gm1);
    } else if (rp >= 0.5*Rp && rp < Rp) {
      cons(IDN,k,j,i) = RHp;
      cons(IM1,k,j,i) = 0.0;
      cons(IM2,k,j,i) = 0.0;
      cons(IM3,k,j,i) = 0.0;
      cons(IEN,k,j,i) = Pp/(RHp*gm1);
    } else if (rp >= Rp && rp <= 1.5*Rp) {
      lambdap = 0.5*pow(v_escp/csp,2);
      etap = rp/Rp;
      b = -3.0 - 4.0*log(lambdap/2.0) + 4.0*log(etap)+2.0*(lambdap/etap);
      o = 0;
      if (rp <= rcp) {
        psi = 2.e-8;
      } else {
        psi = 2.5;
      }
      do {
        function = -psi + log(psi) + b;
        dfunction = -1.0 + (1.0/psi);
        step=function/dfunction;
        psi=psi-step;
        o++;
      } while (fabs(step) > 1.0e-8 && o < 1000);
      v_wp = csp*sqrt(psi);
      if(isnan(v_wp)){
        v_wp = 0.0;
      }
      v_wp_x1 = sin(thetap)*(v_wp*cos(phip)+sin(phis)*rs*omega_fr-sin(phip)*rp*omegap);
      v_wp_x2 = sin(thetap)*(v_wp*sin(phip)-cos(phis)*rs*omega_fr+a*omega_orb + cos(phip)*rp*omegap);
      v_wp_x3 = v_wp*cos(thetap);
      PRSp = Pp*exp(lambdap*(Rp/rp-1.0)-0.5*pow(v_wp/csp,2));
      RHOp = (RHp/Pp)*PRSp;

      cons(IDN,k,j,i) = RHOp;
      cons(IM1,k,j,i) = RHOp*v_wp_x1;
      cons(IM2,k,j,i) = RHOp*v_wp_x2;
      cons(IM3,k,j,i) = RHOp*v_wp_x3;
      cons(IEN,k,j,i) = PRSp/(RHOp*gm1) + 0.5*RHOp*(SQR(v_wp_x1) + SQR(v_wp_x2) + SQR(v_wp_x3));
    }

    Real vx = cons(IM1,k,j,i)/cons(IDN,k,j,i);
    Real vy = cons(IM2,k,j,i)/cons(IDN,k,j,i);

    // - Gravity outside bodies 
    Real gs = -UNIT_G*Ms/rs/rs;
    Real gp = -UNIT_G*Mp/rp/rp;
    // - Gravity inside bodies 
    Real gs_in = -(4.0/3.0)*pi*UNIT_G*RHs;     
    Real gp_in = -(4.0/3.0)*pi*UNIT_G*RHp;
    // - Coriolis and centrifugal forces 
    Real Fin_x = omega_fr*omega_fr*x + 2.0*omega_fr*vy;
    Real Fin_y = omega_fr*omega_fr*y - 2.0*omega_fr*vx;

    Real g_x = 0.0;
    Real g_y = 0.0;
    Real g_z = 0.0;

    if (rs > 1.5*Rs && rp > 1.5*Rp){ 
        g_x = Fin_x + gs*x/rs + gp*(x - a)/rp;
        g_y = Fin_y + gs*y/rs + gp*y/rp;
        g_z = gs*z/rs + gp*z/rp;
    } else if (rs < Rs) { // - Star interal gravity 
        g_x = gs_in*x;
        g_y = gs_in*y;
        g_z = gs_in*z;
    } else if (rp < Rp) { // - Planet interal gravity 
        g_x = gp_in*(x - a);
        g_y = gp_in*y;
        g_z = gp_in*z;
    }

    //cons(IM1,k,j,i) += cons(IDN,k,j,i)*g_x*dt;
    //cons(IM2,k,j,i) += cons(IDN,k,j,i)*g_y*dt;
    //cons(IM3,k,j,i) += cons(IDN,k,j,i)*g_z*dt;

  }}}
  return;
}


int RefinementCondition(MeshBlock *pmb)
{
  AthenaArray<Real> &w = pmb->phydro->w;
  Real maxeps=0.0;
  int k=pmb->ks;
  for(int j=pmb->js; j<=pmb->je; j++) {
    for(int i=pmb->is; i<=pmb->ie; i++) {
      Real epsr= (std::abs(w(IDN,k,j,i+1)-2.0*w(IDN,k,j,i)+w(IDN,k,j,i-1))
                 +std::abs(w(IDN,k,j+1,i)-2.0*w(IDN,k,j,i)+w(IDN,k,j-1,i)))/w(IDN,k,j,i);
      Real epsp= (std::abs(w(IEN,k,j,i+1)-2.0*w(IEN,k,j,i)+w(IEN,k,j,i-1))
                 +std::abs(w(IEN,k,j+1,i)-2.0*w(IEN,k,j,i)+w(IEN,k,j-1,i)))/w(IEN,k,j,i);
      Real eps = std::max(epsr, epsp);
      maxeps = std::max(maxeps, eps);
    }
  }
  if(maxeps > 1.0) return 1;
  //if(maxeps < 2.8) return -1;
  return 0;
}  
