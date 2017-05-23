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

Real TwoDimensionalInterp(MeshBlock *pmb, int &i, int &j, int &k, Real &x, Real &y, 
                            Real &r, Real &ddr, const AthenaArray<Real> &prim);
Real ThreeDimensionalInterp(MeshBlock *pmb, int &i, int &j, int &k, Real &x, Real &y, Real &z, 
                              Real &r, Real &ddr, const AthenaArray<Real> &prim);
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

  Real rho = 0.0;
  Real prs = 0.0;
  Real gs = 0.0;

  Real r = 0.0;
  Real theta = 0.0;
  Real phi = 0.0;

  Real x = 0.0;
  Real y = 0.0;
  Real z = 0.0;

  Real vv = 0.0;
  Real vx = 0.0;
  Real vy = 0.0;
  Real vz = 0.0;

  for (int k=ks; k<=ke; k++) {
  for (int j=js; j<=je; j++) {
  for (int i=is; i<=ie; i++) {

    Real x = pcoord->x1v(i);
    Real y = pcoord->x2v(j);
    Real z = pcoord->x3v(k);

    r = sqrt(SQR(x) + SQR(y) + SQR(z));

    theta = acos(z/r);
    phi = atan2(y,x);

    vv = v_inf*pow(1.0 - (1.0/r),b);
    vx = vv*x/r;
    vy = vv*y/r;
    vz = vv*z/r;

    // Stellar interior.
    if (r < 0.5) {

      rho = M_dot/(4.0*pi*(cs/Cs_p));
      prs = cs*cs*rho/gmma;

      phydro->u(IDN,k,j,i) = rho;
      phydro->u(IM1,k,j,i) = 0.0;
      phydro->u(IM2,k,j,i) = 0.0;
      phydro->u(IM3,k,j,i) = 0.0;
      phydro->u(IEN,k,j,i) = prs/(rho*gm1);

    } else if (r >= 0.5 && r < 1.0) {

      rho = M_dot/(4.0*pi*(cs/Cs_p));
      prs = cs*cs*rho/gmma;

      phydro->u(IDN,k,j,i) = rho;
      phydro->u(IM1,k,j,i) = 0.0;
      phydro->u(IM2,k,j,i) = 0.0;
      phydro->u(IM3,k,j,i) = 0.0;
      phydro->u(IEN,k,j,i) = prs/(rho*gm1);

    } else if (r >= 1.0) {

      rho = M_dot/(4.0*pi*vv*r*r);
      prs = cs*cs*rho/gmma;

      phydro->u(IDN,k,j,i) = rho;
      phydro->u(IM1,k,j,i) = rho*vx;
      phydro->u(IM2,k,j,i) = rho*vy;
      phydro->u(IM3,k,j,i) = rho*vz;
      phydro->u(IEN,k,j,i) = prs/(rho*gm1) + 0.5*rho*(SQR(vx) + SQR(vy) + SQR(vz));

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
  Real omega = 0.5; 
  Real shell = 1.5; 

  //Real Mratio = pmb->pin->GetReal("problem","Mratio");
  //Real Lratio = pin->GetReal("problem","Lratio");
  //Real Bcgs = pin->GetReal("problem","B_CGS");
  //Real T = pin->GetReal("problem","TT");
  //Real mu = pin->GetReal("problem","MU");
  //Real a = pin->GetReal("problem","AA");
  //Real b = pin->GetReal("problem","b_law");
  //Real Q = pin->GetReal("problem","QQ");
  //Real a_eff = pin->GetReal("problem","aa_eff");
  //Real beta = pin->GetReal("problem","BB");
  //Real Cs_p = pin->GetReal("problem","Cs_P");
  //Real omega_fr = pin->GetReal("problem","OMEGA");
  //Real shell = pin->GetReal("problem","SHELL");
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
  Real omega_fr = pow(omega,2)*(8.0/27.0)*UNIT_G*M_star;

  Real rho = 0.0;
  Real prs = 0.0;
  Real gs = 0.0;

  Real r = 0.0;
  Real theta = 0.0;
  Real phi = 0.0;

  Real x = 0.0;
  Real y = 0.0;
  Real z = 0.0;

  Real vv = 0.0;
  Real vx = 0.0;
  Real vy = 0.0;
  Real vz = 0.0;

  Real vrI[2];
  Real ddr = 0.0;
  Real gL = 0.0;

  int p;

  for (int k=pmb->ks; k<=pmb->ke; ++k) {
  for (int j=pmb->js; j<=pmb->je; ++j) {
  for (int i=pmb->is; i<=pmb->ie; ++i) {

    x = pmb->pcoord->x1v(i);
    y = pmb->pcoord->x2v(j);
    z = pmb->pcoord->x3v(k);

    r = sqrt(SQR(x) + SQR(y) + SQR(z));

    theta = acos(z/r);
    phi = atan2(y,x);

    vv = v_inf*pow(1.0 - (1.0/r),b);
    vx = vv*x/r;
    vy = vv*y/r;
    vz = vv*z/r;

    // Stellar interior.
    if (r < 0.5) {
      rho = M_dot/(4.0*pi*(cs/Cs_p));
      prs = cs*cs*rho/gmma;
      cons(IDN,k,j,i) = rho;
      cons(IM1,k,j,i) = 0.0;
      cons(IM2,k,j,i) = 0.0;
      cons(IM3,k,j,i) = 0.0;
      cons(IEN,k,j,i) = prs/(rho*gm1);
    } else if (r >= 0.5 && r < 1.0) {
      rho = M_dot/(4.0*pi*(cs/Cs_p));
      prs = cs*cs*rho/gmma;
      cons(IDN,k,j,i) = rho;
      cons(IM1,k,j,i) = 0.0;
      cons(IM2,k,j,i) = 0.0;
      cons(IM3,k,j,i) = 0.0;
      cons(IEN,k,j,i) = prs/(rho*gm1);
    } else if (r >= 1.0 && r < shell) {
      rho = M_dot/(4.0*pi*vv*r*r);
      prs = cs*cs*rho/gmma;
      cons(IDN,k,j,i) = rho;
      cons(IM1,k,j,i) = rho*vx;
      cons(IM2,k,j,i) = rho*vy;
      cons(IM3,k,j,i) = rho*vz;
      cons(IEN,k,j,i) = prs/(rho*gm1) + 0.5*rho*(SQR(vx) + SQR(vy) + SQR(vz));
    } else if (r >= shell) {

      for (p = 0; p < 2; p++) {
        Real dx = pmb->pcoord->x1v(i+1) - x;
        Real dy = pmb->pcoord->x2v(j+1) - y;
        Real dz = pmb->pcoord->x3v(k+1) - z;
        Real dr = sqrt(SQR(dx) + SQR(dy) + SQR(dz));
        if (p == 0) {
          ddr = -0.9*dr;
        } else {
          ddr = 0.9*dr;
        }
        //vrI[p] = TwoDimensionalInterp(pmb, i, j, k, x, y, r, ddr, prim);
        vrI[p] = ThreeDimensionalInterp(pmb, i, j, k, x, y, z, r, ddr, prim);
      }

      Real vr    = prim(IVX,k,j,i)*x/r + prim(IVY,k,j,i)*y/r + prim(IVZ,k,j,i)*z/r;
      Real dvdr = fabs((vrI[1] - vrI[0])/(2.0*ddr));
      Real nu2_c = (1.0-(1.0/(r*r)));
      Real B     = (prim(IDN,k,j,i)*Q*c*ke);
      Real sigma = (r/fabs(vr))*(dvdr)-1.0;
      Real f     = ((pow(1.0+sigma,1.0+a)-pow(1.0+sigma*nu2_c,1.0+a))/((1.0+a)*(1.0-nu2_c)*sigma*pow(1.0+sigma,a)));
      gL = (f*A*pow(r,-2)*pow(dvdr/B,a)); 

      // Set acceleration from gravity, rotation and mCAK.
      Real gs = -UNIT_G*(M_star - Edd)/r/r;
      Real Fin_x = 0.0;//omega_fr*omega_fr*x + 2.0*omega_fr*prim(IVY,k,j,i);
      Real Fin_y = 0.0;//omega_fr*omega_fr*y - 2.0*omega_fr*prim(IVX,k,j,i);

      Real g_x = Fin_x + gs*x/r + gL*x/r;
      Real g_y = Fin_y + gs*y/r + gL*y/r;
      Real g_z = gs*z/r + gL*z/r;

      cons(IM1,k,j,i) += prim(IDN,k,j,i)*g_x*dt;
      cons(IM2,k,j,i) += prim(IDN,k,j,i)*g_y*dt;
      cons(IM3,k,j,i) += prim(IDN,k,j,i)*g_z*dt;

/*
      if (cons(IDN,k,j,i) < 1.0e-10 || prim(IDN,k,j,i) < 1.0e-10) {
        std::cout << "vrI[0] =" << vrI[0] << "\n";
        std::cout << "vrI[1] =" << vrI[1] << "\n";
        std::cout << "nu2_c =" << nu2_c << "\n";
        std::cout << "B =" << B << "\n";
        std::cout << "sigma =" << sigma << "\n";
        std::cout << "f =" << f << "\n";
        std::cout << "dvdr =" << dvdr << "\n";
        std::cout << "gL =" << gL << "\n";
      }
*/

    }

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
                 +std::abs(w(IDN,k,j+1,i)-2.0*w(IDN,k,j,i)+w(IDN,k,j-1,i))
                 +std::abs(w(IDN,k+1,j,i)-2.0*w(IDN,k,j,i)+w(IDN,k-1,j,i)))/w(IDN,k,j,i);
      Real epsp= (std::abs(w(IEN,k,j,i+1)-2.0*w(IEN,k,j,i)+w(IEN,k,j,i-1))
                 +std::abs(w(IEN,k,j+1,i)-2.0*w(IEN,k,j,i)+w(IEN,k,j-1,i))
                 +std::abs(w(IEN,k+1,j,i)-2.0*w(IEN,k,j,i)+w(IEN,k-1,j,i)))/w(IEN,k,j,i);
      Real eps = std::max(epsr, epsp);
      maxeps = std::max(maxeps, eps);
    }
  }
  if(maxeps > 0.075) return 1;
  //if(maxeps < 2.8) return -1;
  return 0;
} 

Real TwoDimensionalInterp(MeshBlock *pmb, int &i, int &j, int &k, Real &x, Real &y, 
                            Real &r, Real &ddr, const AthenaArray<Real> &prim)
{
/*
       _____________________________
  j+1 |              |              | 
      |              |              |
      |              |     *        |
      |              |              |
      |              |              |
    j |______________|______________|
      |              |              | 
      |              |              |
      |              |              |
      |              |              |
      |              |              |
  j-1 |______________|______________|

     i - 1          i            i + 1


  yb 3               4
      |```|`````````|
      |   |         |
   yI |___|_________|
      |   |         |
      |   |         |
      |___|_________|
  ya 1    xI         2
      xa            xb
*/

  int u, s;
  Real Ntot;
  Real vrI;
  Real vI[2];
  Real N[4];
  Real V[2][4];
  Real xa, xb;
  Real ya, yb;
  Real theta = atan2(x,y);
  Real xI = (r+ddr)*sin(theta);
  Real yI = (r+ddr)*cos(theta);

  if (xI > x && yI > y){
    xa = pmb->pcoord->x1v(i); xb = pmb->pcoord->x1v(i+1);
    ya = pmb->pcoord->x2v(j); yb = pmb->pcoord->x2v(j+1);
    for (u = IVX; u < IVY+1; u++) {
      V[u-1][0] = prim(u,k,j,i);
      V[u-1][1] = prim(u,k,j,i+1);
      V[u-1][2] = prim(u,k,j+1,i);
      V[u-1][3] = prim(u,k,j+1,i+1);      
    }
  } else if (xI < x && yI > y){
    xa = pmb->pcoord->x1v(i-1); xb = pmb->pcoord->x1v(i);
    ya = pmb->pcoord->x2v(j); yb = pmb->pcoord->x2v(j+1);
    for (u = IVX; u < IVY+1; u++) {
      V[u-1][0] = prim(u,k,j,i-1);
      V[u-1][1] = prim(u,k,j,i);
      V[u-1][2] = prim(u,k,j+1,i-1);
      V[u-1][3] = prim(u,k,j+1,i);
    }
  } else if (xI < x && yI < y){
    xa = pmb->pcoord->x1v(i-1); xb = pmb->pcoord->x1v(i);
    ya = pmb->pcoord->x2v(j-1); yb = pmb->pcoord->x2v(j);
    for (u = IVX; u < IVY+1; u++) {
      V[u-1][0] = prim(u,k,j-1,i-1);
      V[u-1][1] = prim(u,k,j-1,i);
      V[u-1][2] = prim(u,k,j,i-1);
      V[u-1][3] = prim(u,k,j,i);
    }
  } else if (xI > x && yI < y){
    xa = pmb->pcoord->x1v(i); xb = pmb->pcoord->x1v(i+1);
    ya = pmb->pcoord->x2v(j-1); yb = pmb->pcoord->x2v(j);
    for (u = IVX; u < IVY+1; u++) {
      V[u-1][0] = prim(u,k,j-1,i);
      V[u-1][1] = prim(u,k,j-1,i+1);
      V[u-1][2] = prim(u,k,j,i);
      V[u-1][3] = prim(u,k,j,i+1);
    }
  }

  // Find total volume.
  N[0] = (xb - xI)*(yb - yI);
  N[1] = (xI - xa)*(yb - yI);
  N[2] = (xb - xI)*(yI - ya);
  N[3] = (xI - xa)*(yI - ya);
  Ntot = N[0] + N[1] + N[2] + N[3];
  // Normalise volumes by total.
  N[0] /= Ntot; 
  N[1] /= Ntot; 
  N[2] /= Ntot; 
  N[3] /= Ntot;

  // ==========================================
  // vI contains the interpolated velovities 
  // 
  // vI[0] = x velocity 
  // vI[1] = y "
  //
  // V contains the velocity componants at 
  // each of the sample points.
  // 
  // V[0, :] = x velocity at each sample point.
  // V[1, :] = y velocity at each sample point.
  //
  // N contains the waighting volumes for the 
  // the corner velocities.
  // ==========================================

  vI[0] = 0.0; vI[1] = 0.0;
  for (s = 0; s < 4; s++) {
    vI[0] += V[0][s]*N[s];
    vI[1] += V[1][s]*N[s];
  }
  vrI = (vI[0]*x + vI[1]*y)/r;
  return vrI;

}


Real ThreeDimensionalInterp(MeshBlock *pmb, int &i, int &j, int &k, Real &x, Real &y, 
                              Real &z, Real &r, Real &ddr, const AthenaArray<Real> &prim)
{
  int u, s;
  Real Ntot;
  Real vrI;
  Real vI[3];
  Real N[8];
  Real V[3][8];
  Real xa, xb;
  Real ya, yb;
  Real za, zb;
  Real phi = atan2(y,x);
  Real theta = acos(z/r);
  Real xI = (r+ddr)*sin(theta)*cos(phi);
  Real yI = (r+ddr)*sin(theta)*sin(phi);
  Real zI = (r+ddr)*cos(theta);

  if (xI > x && yI > y && zI > z){
    xa = pmb->pcoord->x1v(i); xb = pmb->pcoord->x1v(i+1);
    ya = pmb->pcoord->x2v(j); yb = pmb->pcoord->x2v(j+1);
    za = pmb->pcoord->x3v(k); zb = pmb->pcoord->x3v(k+1);
    for (u = IVX; u < IVZ+1; u++) {
      V[u-1][0] = prim(u,k,j,i);
      V[u-1][1] = prim(u,k,j,i+1);
      V[u-1][2] = prim(u,k,j+1,i);
      V[u-1][3] = prim(u,k,j+1,i+1);
      V[u-1][4] = prim(u,k+1,j,i);
      V[u-1][5] = prim(u,k+1,j,i+1);
      V[u-1][6] = prim(u,k+1,j+1,i);
      V[u-1][7] = prim(u,k+1,j+1,i+1);
    }
  }
  if (xI < x && yI > y && zI > z){
    xa = pmb->pcoord->x1v(i-1); xb = pmb->pcoord->x1v(i);
    ya = pmb->pcoord->x2v(j); yb = pmb->pcoord->x2v(j+1);
    za = pmb->pcoord->x3v(k); zb = pmb->pcoord->x3v(k+1);
    for (u = IVX; u < IVZ+1; u++) {
      V[u-1][0] = prim(u,k,j,i-1);
      V[u-1][1] = prim(u,k,j,i);
      V[u-1][2] = prim(u,k,j+1,i-1);
      V[u-1][3] = prim(u,k,j+1,i);
      V[u-1][4] = prim(u,k+1,j,i-1);
      V[u-1][5] = prim(u,k+1,j,i);
      V[u-1][6] = prim(u,k+1,j+1,i-1);
      V[u-1][7] = prim(u,k+1,j+1,i);
    }
  }
  if (xI > x && yI < y && zI > z){
    xa = pmb->pcoord->x1v(i); xb = pmb->pcoord->x1v(i+1);
    ya = pmb->pcoord->x2v(j-1); yb = pmb->pcoord->x2v(j);
    za = pmb->pcoord->x3v(k); zb = pmb->pcoord->x3v(k+1);
    for (u = IVX; u < IVZ+1; u++) {
      V[u-1][0] = prim(u,k,j-1,i);
      V[u-1][1] = prim(u,k,j-1,i+1);
      V[u-1][2] = prim(u,k,j,i);
      V[u-1][3] = prim(u,k,j,i+1);
      V[u-1][4] = prim(u,k+1,j-1,i);
      V[u-1][5] = prim(u,k+1,j-1,i+1);
      V[u-1][6] = prim(u,k+1,j,i);
      V[u-1][7] = prim(u,k+1,j,i+1);
    }
  }
  if (xI < x && yI < y && zI > z){
    xa = pmb->pcoord->x1v(i-1); xb = pmb->pcoord->x1v(i);
    ya = pmb->pcoord->x2v(j-1); yb = pmb->pcoord->x2v(j);
    za = pmb->pcoord->x3v(k); zb = pmb->pcoord->x3v(k+1);
    for (u = IVX; u < IVZ+1; u++) {
      V[u-1][0] = prim(u,k,j-1,i-1);
      V[u-1][1] = prim(u,k,j-1,i);
      V[u-1][2] = prim(u,k,j,i-1);
      V[u-1][3] = prim(u,k,j,i);
      V[u-1][4] = prim(u,k+1,j-1,i-1);
      V[u-1][5] = prim(u,k+1,j-1,i);
      V[u-1][6] = prim(u,k+1,j,i-1);
      V[u-1][7] = prim(u,k+1,j,i);
    }
  }

  // zI < zP 
  if (xI > x && yI > y && zI < z){
    xa = pmb->pcoord->x1v(i); xb = pmb->pcoord->x1v(i+1);
    ya = pmb->pcoord->x2v(j); yb = pmb->pcoord->x2v(j+1);
    za = pmb->pcoord->x3v(k-1); zb = pmb->pcoord->x3v(k);
    for (u = IVX; u < IVZ+1; u++) {
      V[u-1][0] = prim(u,k-1,j,i);
      V[u-1][1] = prim(u,k-1,j,i+1);
      V[u-1][2] = prim(u,k-1,j+1,i);
      V[u-1][3] = prim(u,k-1,j+1,i+1);
      V[u-1][4] = prim(u,k,j,i);
      V[u-1][5] = prim(u,k,j,i+1);
      V[u-1][6] = prim(u,k,j+1,i);
      V[u-1][7] = prim(u,k,j+1,i+1);
    }
  }
  if (xI < x && yI > y && zI < z){
    xa = pmb->pcoord->x1v(i-1); xb = pmb->pcoord->x1v(i);
    ya = pmb->pcoord->x2v(j); yb = pmb->pcoord->x2v(j+1);
    za = pmb->pcoord->x3v(k-1); zb = pmb->pcoord->x3v(k);
    for (u = IVX; u < IVZ+1; u++) {
      V[u-1][0] = prim(u,k-1,j,i-1);
      V[u-1][1] = prim(u,k-1,j,i);
      V[u-1][2] = prim(u,k-1,j+1,i-1);
      V[u-1][3] = prim(u,k-1,j+1,i);
      V[u-1][4] = prim(u,k,j,i-1);
      V[u-1][5] = prim(u,k,j,i);
      V[u-1][6] = prim(u,k,j+1,i-1);
      V[u-1][7] = prim(u,k,j+1,i);
    }
  }
  if (xI > x && yI < y && zI < z){
    xa = pmb->pcoord->x1v(i); xb = pmb->pcoord->x1v(i+1);
    ya = pmb->pcoord->x2v(j-1); yb = pmb->pcoord->x2v(j);
    za = pmb->pcoord->x3v(k-1); zb = pmb->pcoord->x3v(k);
    for (u = IVX; u < IVZ+1; u++) {
      V[u-1][0] = prim(u,k-1,j-1,i);
      V[u-1][1] = prim(u,k-1,j-1,i+1);
      V[u-1][2] = prim(u,k-1,j,i);
      V[u-1][3] = prim(u,k-1,j,i+1);
      V[u-1][4] = prim(u,k,j-1,i);
      V[u-1][5] = prim(u,k,j-1,i+1);
      V[u-1][6] = prim(u,k,j,i);
      V[u-1][7] = prim(u,k,j,i+1);
    }
  }
  if (xI < x && yI < y && zI < z){
    xa = pmb->pcoord->x1v(i-1); xb = pmb->pcoord->x1v(i);
    ya = pmb->pcoord->x2v(j-1); yb = pmb->pcoord->x2v(j);
    za = pmb->pcoord->x3v(k-1); zb = pmb->pcoord->x3v(k);
    for (u = IVX; u < IVZ+1; u++) {
      V[u-1][0] = prim(u,k-1,j-1,i-1);
      V[u-1][1] = prim(u,k-1,j-1,i);
      V[u-1][2] = prim(u,k-1,j,i-1);
      V[u-1][3] = prim(u,k-1,j,i);
      V[u-1][4] = prim(u,k,j-1,i-1);
      V[u-1][5] = prim(u,k,j-1,i);
      V[u-1][6] = prim(u,k,j,i-1);
      V[u-1][7] = prim(u,k,j,i);
    }
  }

  // Find total volume.
  N[0] = (xb - xI)*(yb - yI)*(zb - zI);
  N[1] = (xI - xa)*(yb - yI)*(zb - zI);
  N[2] = (xb - xI)*(yI - ya)*(zb - zI);
  N[3] = (xI - xa)*(yI - ya)*(zb - zI);
  N[4] = (xb - xI)*(yb - yI)*(zI - za);
  N[5] = (xI - xa)*(yb - yI)*(zI - za);
  N[6] = (xb - xI)*(yI - ya)*(zI - za);
  N[7] = (xI - xa)*(yI - ya)*(zI - za);
  Ntot = N[0] + N[1] + N[2] + N[3] + N[4] + N[5] + N[6] + N[7];
  // Normalise volumes by total.
  N[0] /= Ntot; 
  N[1] /= Ntot; 
  N[2] /= Ntot; 
  N[3] /= Ntot;
  N[4] /= Ntot; 
  N[5] /= Ntot; 
  N[6] /= Ntot; 
  N[7] /= Ntot;
         
  // ==========================================
  // vI contains the interpolated velovities 
  // 
  // vI[0] = x velocity 
  // vI[1] = y "
  // vI[2] = z "
  //
  // V contains the velocity componants at 
  // each of the sample points.
  // 
  // V[0][:] = x velocity at each sample point.
  // V[1][:] = y velocity at each sample point.
  // V[2][:] = z velocity at each sample point.
  //
  // N contains the waighting volumes for the 
  // the corner velocities.
  // ==========================================

  vI[0] = 0.0; vI[1] = 0.0; vI[2] = 0.0;
  for (s = 0; s < 8; s++) {
    vI[0] += V[0][s]*N[s];
    vI[1] += V[1][s]*N[s];
    vI[2] += V[2][s]*N[s];
  }
  vrI = (vI[0]*x + vI[1]*y + vI[2]*z)/r;

 
  if (isnan(vrI)) {
    std::cout << "r =" << r << ", ddr = " << ddr << ", r+ddr =" << r+ddr << "\n";
    std::cout << "x  = " << x << ", y  =" << y << " z =" << z << "\n";
    std::cout << "xa = " << xa << ", ya = " << ya << " za =" << za << "\n";
    std::cout << "xI = " << xI << ", yI = " << yI << " zI =" << zI << "\n";
    std::cout << "xb = " << xb << ", yb = " << yb << " zb =" << zb << "\n";
    std::cout << "N[0] = " << N[0] << "\n";
    std::cout << "N[1] = " << N[1] << "\n";
    std::cout << "N[2] = " << N[2] << "\n";
    std::cout << "N[3] = " << N[3] << "\n";
    std::cout << "N[4] = " << N[4] << "\n";
    std::cout << "N[5] = " << N[5] << "\n";
    std::cout << "N[6] = " << N[6] << "\n";
    std::cout << "N[7] = " << N[7] << "\n";
    std::cout << "Ntot = " << N[0]+N[1]+N[2]+N[3]+N[4]+N[5]+N[6]+N[7] << "\n";
    std::cout << "V[0][0] = " << V[0][0] << ", V[1][0] = " << V[1][0] << ", V[2][0] = " << V[2][0] << "\n";
    std::cout << "V[0][1] = " << V[0][1] << ", V[1][1] = " << V[1][1] << ", V[2][1] = " << V[2][1] << "\n";
    std::cout << "V[0][2] = " << V[0][2] << ", V[1][2] = " << V[1][2] << ", V[2][2] = " << V[2][2] << "\n";
    std::cout << "V[0][3] = " << V[0][3] << ", V[1][3] = " << V[1][3] << ", V[2][3] = " << V[2][3] << "\n";
    std::cout << "V[0][4] = " << V[0][4] << ", V[1][4] = " << V[1][4] << ", V[2][4] = " << V[2][4] << "\n";
    std::cout << "V[0][5] = " << V[0][5] << ", V[1][5] = " << V[1][5] << ", V[2][5] = " << V[2][5] << "\n";
    std::cout << "V[0][6] = " << V[0][6] << ", V[1][6] = " << V[1][6] << ", V[2][6] = " << V[2][6] << "\n";
    std::cout << "V[0][7] = " << V[0][7] << ", V[1][7] = " << V[1][7] << ", V[2][7] = " << V[2][7] << "\n";
    std::cout << "vI[0] = " << vI[0] << "\n";
    std::cout << "vI[1] = " << vI[1] << "\n";
    std::cout << "vI[2] = " << vI[2] << "\n";
    std::cout << "vrI = " << vrI << "\n";
    std::cout << "\n";
  }


  return vrI;
}


