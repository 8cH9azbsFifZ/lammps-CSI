/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under 
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing author: Eric Simon (Cray)
------------------------------------------------------------------------- */

#include "math.h"
#include "stdlib.h"
#include "dihedral_class2.h"
#include "atom.h"
#include "neighbor.h"
#include "update.h"
#include "domain.h"
#include "comm.h"
#include "force.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;

#define MIN(A,B) ((A) < (B)) ? (A) : (B)
#define MAX(A,B) ((A) > (B)) ? (A) : (B)

#define TOLERANCE 0.05
#define SMALL     0.0000001

/* ---------------------------------------------------------------------- */

DihedralClass2::DihedralClass2(LAMMPS *lmp) : Dihedral(lmp)
{
  PI = 4.0*atan(1.0);
}

/* ---------------------------------------------------------------------- */

DihedralClass2::~DihedralClass2()
{
  if (allocated) {
    memory->sfree(setflag);
    memory->sfree(setflag_d);
    memory->sfree(setflag_mbt);
    memory->sfree(setflag_ebt);
    memory->sfree(setflag_at);
    memory->sfree(setflag_aat);
    memory->sfree(setflag_bb13t);

    memory->sfree(k1);
    memory->sfree(k2);
    memory->sfree(k3);
    memory->sfree(phi1);
    memory->sfree(phi2);
    memory->sfree(phi3);

    memory->sfree(mbt_f1);
    memory->sfree(mbt_f2);
    memory->sfree(mbt_f3);
    memory->sfree(mbt_r0);

    memory->sfree(ebt_f1_1);
    memory->sfree(ebt_f2_1);
    memory->sfree(ebt_f3_1);
    memory->sfree(ebt_r0_1);

    memory->sfree(ebt_f1_2);
    memory->sfree(ebt_f2_2);
    memory->sfree(ebt_f3_2);
    memory->sfree(ebt_r0_2);

    memory->sfree(at_f1_1);
    memory->sfree(at_f2_1);
    memory->sfree(at_f3_1);
    memory->sfree(at_theta0_1);

    memory->sfree(at_f1_2);
    memory->sfree(at_f2_2);
    memory->sfree(at_f3_2);
    memory->sfree(at_theta0_2);

    memory->sfree(aat_k);
    memory->sfree(aat_theta0_1);
    memory->sfree(aat_theta0_2);

    memory->sfree(bb13t_k);
    memory->sfree(bb13t_r10);
    memory->sfree(bb13t_r30);
  }
}

/* ---------------------------------------------------------------------- */

void DihedralClass2::compute(int eflag, int vflag)
{
  int i1,i2,i3,i4,i,j,k,n,type,factor;
  double rfactor;
  double delx1,dely1,delz1,dely2,delz2,delx2m,dely2m,delz2m;
  double delx2,dely3,delz3,r1mag2,r1,r2mag2,r2,r3mag2,r3;
  double sb1,rb1,sb2,rb2,sb3,rb3,c0,r12c1;
  double r12c2,costh12,costh13,costh23,sc1,sc2,s1,s2,c;
  double cosphi,phi,sinphi,a11,a22,a33,a12,a13,a23,sx1,sx2;
  double sx12,sy1,sy2,sy12,sz1,sz2,sz12,dphi1,dphi2,dphi3;
  double de_dihedral,t1,t2,t3,t4,cos2phi,cos3phi,bt1,bt2;
  double bt3,sumbte,db,sumbtf,at1,at2,at3,da,da1,da2,r1_0;
  double r3_0,dr1,dr2,tk1,tk2,vx1,vx2,vx3,vy1,vy2,vy3,vz1;
  double vz2,vz3,delx3,s12,sin2;
  double dcosphidr[4][3],dphidr[4][3],dbonddr[3][4][3],dthetadr[2][4][3];
  double fabcd[4][3];

  energy = 0.0;
  if (vflag) for (n = 0; n < 6; n++) virial[n] = 0.0;

  double **x = atom->x;
  double **f = atom->f;
  int **dihedrallist = neighbor->dihedrallist;
  int ndihedrallist = neighbor->ndihedrallist;
  int nlocal = atom->nlocal;
  int newton_bond = force->newton_bond;

  for (n = 0; n < ndihedrallist; n++) {

    i1 = dihedrallist[n][0];
    i2 = dihedrallist[n][1];
    i3 = dihedrallist[n][2];
    i4 = dihedrallist[n][3];
    type = dihedrallist[n][4];

    if (newton_bond) factor = 4;
    else {
      factor = 0;
      if (i1 < nlocal) factor++;
      if (i2 < nlocal) factor++;
      if (i3 < nlocal) factor++;
      if (i4 < nlocal) factor++;
      }
    rfactor = 0.25 * factor;

    // 1st bond

    delx1 = x[i1][0] - x[i2][0];
    dely1 = x[i1][1] - x[i2][1];
    delz1 = x[i1][2] - x[i2][2];
    domain->minimum_image(delx1,dely1,delz1);

    // 2nd bond

    delx2 = x[i3][0] - x[i2][0];
    dely2 = x[i3][1] - x[i2][1];
    delz2 = x[i3][2] - x[i2][2];
    domain->minimum_image(delx2,dely2,delz2);

    delx2m = -delx2;
    dely2m = -dely2;
    delz2m = -delz2;
    domain->minimum_image(delx2m,dely2m,delz2m);

    // 3rd bond

    delx3 = x[i4][0] - x[i3][0];
    dely3 = x[i4][1] - x[i3][1];
    delz3 = x[i4][2] - x[i3][2];
    domain->minimum_image(delx3,dely3,delz3);

    // distances

    r1mag2 = delx1*delx1 + dely1*dely1 + delz1*delz1;
    r1 = sqrt(r1mag2);
    r2mag2 = delx2*delx2 + dely2*dely2 + delz2*delz2;
    r2 = sqrt(r2mag2);
    r3mag2 = delx3*delx3 + dely3*dely3 + delz3*delz3;
    r3 = sqrt(r3mag2);

    sb1 = 1.0/r1mag2;
    rb1 = 1.0/r1;
          
    sb2 = 1.0/r2mag2;
    rb2 = 1.0/r2;

    sb3 = 1.0/r3mag2;
    rb3 = 1.0/r3;

    c0 = (delx1*delx3 + dely1*dely3 + delz1*delz3) * rb1*rb3;

    // angles

    r12c1 = rb1*rb2;
    r12c2 = rb2*rb3;
    costh12 = (delx1*delx2 + dely1*dely2 + delz1*delz2) * r12c1;
    costh13 = c0;
    costh23 = (delx2m*delx3 + dely2m*dely3 + delz2m*delz3) * r12c2;
          
    // cos and sin of 2 angles and final c

    sin2 = MAX(1.0 - costh12*costh12,0.0);
    sc1 = sqrt(sin2);
    if (sc1 < SMALL) sc1 = SMALL;
    sc1 = 1.0/sc1;
          
    sin2 = MAX(1.0 - costh23*costh23,0.0);
    sc2 = sqrt(sin2);
    if (sc2 < SMALL) sc2 = SMALL;
    sc2 = 1.0/sc2;
          
    s1 = sc1 * sc1;
    s2 = sc2 * sc2;
    s12 = sc1 * sc2;
    c = (c0 + costh12*costh23) * s12;

    // error check

    if (c > 1.0 + TOLERANCE || c < (-1.0 - TOLERANCE)) {
      int me;
      MPI_Comm_rank(world,&me);
      if (screen) {
	fprintf(screen,"Dihedral problem: %d %d %d %d %d %d\n",
		me,update->ntimestep,
		atom->tag[i1],atom->tag[i2],atom->tag[i3],atom->tag[i4]);
	fprintf(screen,"  1st atom: %d %g %g %g\n",
		me,x[i1][0],x[i1][1],x[i1][2]);
	fprintf(screen,"  2nd atom: %d %g %g %g\n",
		me,x[i2][0],x[i2][1],x[i2][2]);
	fprintf(screen,"  3rd atom: %d %g %g %g\n",
		me,x[i3][0],x[i3][1],x[i3][2]);
	fprintf(screen,"  4th atom: %d %g %g %g\n",
		me,x[i4][0],x[i4][1],x[i4][2]);
      }
    }

    if (c > 1.0) c = 1.0;
    if (c < -1.0) c = -1.0;
    cosphi = c;
    phi = acos(c);

    sinphi = sqrt(1.0 - c*c);
    sinphi = MAX(sinphi,SMALL);

    a11 = -c*sb1*s1;
    a22 = sb2 * (2.0*costh13*s12 - c*(s1+s2));
    a33 = -c*sb3*s2;
    a12 = r12c1 * (costh12*c*s1 + costh23*s12);
    a13 = rb1*rb3*s12;
    a23 = r12c2 * (-costh23*c*s2 - costh12*s12);
          
    sx1  = a11*delx1 + a12*delx2 + a13*delx3;
    sx2  = a12*delx1 + a22*delx2 + a23*delx3;
    sx12 = a13*delx1 + a23*delx2 + a33*delx3;
    sy1  = a11*dely1 + a12*dely2 + a13*dely3;
    sy2  = a12*dely1 + a22*dely2 + a23*dely3;
    sy12 = a13*dely1 + a23*dely2 + a33*dely3;
    sz1  = a11*delz1 + a12*delz2 + a13*delz3;
    sz2  = a12*delz1 + a22*delz2 + a23*delz3;
    sz12 = a13*delz1 + a23*delz2 + a33*delz3;

    // set up d(cos(phi))/d(r) and dphi/dr arrays

    dcosphidr[0][0] = -sx1;
    dcosphidr[0][1] = -sy1;
    dcosphidr[0][2] = -sz1;
    dcosphidr[1][0] = sx2 + sx1;
    dcosphidr[1][1] = sy2 + sy1;
    dcosphidr[1][2] = sz2 + sz1;
    dcosphidr[2][0] = sx12 - sx2;
    dcosphidr[2][1] = sy12 - sy2;
    dcosphidr[2][2] = sz12 - sz2;
    dcosphidr[3][0] = -sx12;
    dcosphidr[3][1] = -sy12;
    dcosphidr[3][2] = -sz12;

    for (i = 0; i < 4; i++)
      for (j = 0; j < 3; j++)
	dphidr[i][j] = -dcosphidr[i][j] / sinphi;

    // energy

    dphi1 = phi - phi1[type];
    dphi2 = 2.0*phi - phi2[type];
    dphi3 = 3.0*phi - phi3[type];
    
    if (eflag) energy += rfactor * (k1[type]*(1.0 - cos(dphi1)) +
				    k2[type]*(1.0 - cos(dphi2)) +
				    k3[type]*(1.0 - cos(dphi3)));
          
    de_dihedral = k1[type]*sin(dphi1) + 2.0*k2[type]*sin(dphi2) +
      3.0*k3[type]*sin(dphi3);

    // torsion forces on all 4 atoms

    for (i = 0; i < 4; i++)
      for (j = 0; j < 3; j++)
	fabcd[i][j] = de_dihedral*dphidr[i][j];

    // set up d(bond)/d(r) array
    // dbonddr(i,j,k) = bond i, atom j, coordinate k

    for (i = 0; i < 3; i++)
      for (j = 0; j < 4; j++)
	for (k = 0; k < 3; k++)
	  dbonddr[i][j][k] = 0.0;
    
    // bond1
    
    dbonddr[0][0][0] = delx1 / r1;
    dbonddr[0][0][1] = dely1 / r1;
    dbonddr[0][0][2] = delz1 / r1;
    dbonddr[0][1][0] = -delx1 / r1;
    dbonddr[0][1][1] = -dely1 / r1;
    dbonddr[0][1][2] = -delz1 / r1;

    // bond2

    dbonddr[1][1][0] = delx2 / r2;
    dbonddr[1][1][1] = dely2 / r2;
    dbonddr[1][1][2] = delz2 / r2;
    dbonddr[1][2][0] = -delx2 / r2;
    dbonddr[1][2][1] = -dely2 / r2;
    dbonddr[1][2][2] = -delz2 / r2;

    // bond3
    
    dbonddr[2][2][0] = delx3 / r3;
    dbonddr[2][2][1] = dely3 / r3;
    dbonddr[2][2][2] = delz3 / r3;
    dbonddr[2][3][0] = -delx3 / r3;
    dbonddr[2][3][1] = -dely3 / r3;
    dbonddr[2][3][2] = -delz3 / r3;

    // set up d(theta)/d(r) array
    // dthetadr(i,j,k) = angle i, atom j, coordinate k

    for (i = 0; i < 2; i++)
      for (j = 0; j < 4; j++)
	for (k = 0; k < 3; k++)
	  dthetadr[i][j][k] = 0.0;
    
    t1 = costh12 / r1mag2;
    t2 = costh23 / r2mag2;
    t3 = costh12 / r2mag2;
    t4 = costh23 / r3mag2;
    
    // angle12
    
    dthetadr[0][0][0] = sc1 * ((t1 * delx1) - (delx2 * r12c1));
    dthetadr[0][0][1] = sc1 * ((t1 * dely1) - (dely2 * r12c1));
    dthetadr[0][0][2] = sc1 * ((t1 * delz1) - (delz2 * r12c1));
    
    dthetadr[0][1][0] = sc1 * ((-t1 * delx1) + (delx2 * r12c1) +
			       (-t3 * delx2) + (delx1 * r12c1));
    dthetadr[0][1][1] = sc1 * ((-t1 * dely1) + (dely2 * r12c1) +
			       (-t3 * dely2) + (dely1 * r12c1));
    dthetadr[0][1][2] = sc1 * ((-t1 * delz1) + (delz2 * r12c1) +
			       (-t3 * delz2) + (delz1 * r12c1));
    
    dthetadr[0][2][0] = sc1 * ((t3 * delx2) - (delx1 * r12c1)); 
    dthetadr[0][2][1] = sc1 * ((t3 * dely2) - (dely1 * r12c1));
    dthetadr[0][2][2] = sc1 * ((t3 * delz2) - (delz1 * r12c1));
    
    // angle23
    
    dthetadr[1][1][0] = sc2 * ((t2 * delx2) + (delx3 * r12c2));
    dthetadr[1][1][1] = sc2 * ((t2 * dely2) + (dely3 * r12c2));
    dthetadr[1][1][2] = sc2 * ((t2 * delz2) + (delz3 * r12c2));
    
    dthetadr[1][2][0] = sc2 * ((-t2 * delx2) - (delx3 * r12c2) +
			       (t4 * delx3) + (delx2 * r12c2));
    dthetadr[1][2][1] = sc2 * ((-t2 * dely2) - (dely3 * r12c2) +
			       (t4 * dely3) + (dely2 * r12c2));
    dthetadr[1][2][2] = sc2 * ((-t2 * delz2) - (delz3 * r12c2) +
			       (t4 * delz3) + (delz2 * r12c2));
    
    dthetadr[1][3][0] = -sc2 * ((t4 * delx3) + (delx2 * r12c2));
    dthetadr[1][3][1] = -sc2 * ((t4 * dely3) + (dely2 * r12c2));
    dthetadr[1][3][2] = -sc2 * ((t4 * delz3) + (delz2 * r12c2));
    
    // mid-bond/torsion coupling
    // energy on bond2 (middle bond)
    
    cos2phi = cos(2.0*phi);
    cos3phi = cos(3.0*phi);
    
    bt1 = mbt_f1[type] * cosphi;
    bt2 = mbt_f2[type] * cos2phi;
    bt3 = mbt_f3[type] * cos3phi;
    sumbte = bt1 + bt2 + bt3;
    db = r2 - mbt_r0[type];
    if (eflag) energy += rfactor * db * sumbte;
    
    // force on bond2
    
    bt1 = -mbt_f1[type] * sinphi;
    bt2 = -2.0 * mbt_f2[type] * sin(2.0*phi);
    bt3 = -3.0 * mbt_f3[type] * sin(3.0*phi);
    sumbtf = bt1 + bt2 + bt3;
    
    for (i = 0; i < 4; i++)
      for (j = 0; j < 3; j++)
	fabcd[i][j] += db*sumbtf*dphidr[i][j] + sumbte*dbonddr[1][i][j];

    // end-bond/torsion coupling
    // energy on bond1 (first bond)

    bt1 = ebt_f1_1[type] * cosphi;
    bt2 = ebt_f2_1[type] * cos2phi;
    bt3 = ebt_f3_1[type] * cos3phi;
    sumbte = bt1 + bt2 + bt3;

    db = r1 - ebt_r0_1[type];
    if (eflag) energy += rfactor * db * (bt1+bt2+bt3);

    // force on bond1

    bt1 = ebt_f1_1[type] * sinphi;
    bt2 = 2.0 * ebt_f2_1[type] * sin(2.0*phi);
    bt3 = 3.0 * ebt_f3_1[type] * sin(3.0*phi);
    sumbtf = bt1 + bt2 + bt3;

    for (i = 0; i < 4; i++)
      for (j = 0; j < 3; j++)
	fabcd[i][j] -= db*sumbtf*dphidr[i][j] + sumbte*dbonddr[0][i][j];

    // end-bond/torsion coupling
    // energy on bond3 (last bond)

    bt1 = ebt_f1_2[type] * cosphi;
    bt2 = ebt_f2_2[type] * cos2phi;
    bt3 = ebt_f3_2[type] * cos3phi;
    sumbte = bt1 + bt2 + bt3;

    db = r3 - ebt_r0_2[type];
    if (eflag) energy += rfactor * db * (bt1+bt2+bt3);

    // force on bond3

    bt1 = -ebt_f1_2[type] * sinphi;
    bt2 = -2.0 * ebt_f2_2[type] * sin(2.0*phi);
    bt3 = -3.0 * ebt_f3_2[type] * sin(3.0*phi);
    sumbtf = bt1 + bt2 + bt3;

    for (i = 0; i < 4; i++)
      for (j = 0; j < 3; j++)
	fabcd[i][j] += db*sumbtf*dphidr[i][j] + sumbte*dbonddr[2][i][j];

    // angle/torsion coupling
    // energy on angle1

    at1 = at_f1_1[type] * cosphi;
    at2 = at_f2_1[type] * cos2phi;
    at3 = at_f3_1[type] * cos3phi;
    sumbte = at1 + at2 + at3;

    da = acos(costh12) - at_theta0_1[type];
    if (eflag) energy += rfactor * da * (at1+at2+at3);

    // force on angle1

    bt1 = at_f1_1[type] * sinphi;
    bt2 = 2.0 * at_f2_1[type] * sin(2.0*phi);
    bt3 = 3.0 * at_f3_1[type] * sin(3.0*phi);
    sumbtf = bt1 + bt2 + bt3;

    for (i = 0; i < 4; i++)
      for (j = 0; j < 3; j++)
	fabcd[i][j] -= da*sumbtf*dphidr[i][j] + sumbte*dthetadr[0][i][j];

    // energy on angle2

    at1 = at_f1_2[type] * cosphi;
    at2 = at_f2_2[type] * cos2phi;
    at3 = at_f3_2[type] * cos3phi;
    sumbte = at1 + at2 + at3;

    da = acos(costh23) - at_theta0_2[type];
    if (eflag) energy += rfactor *da * (at1+at2+at3);

    // force on angle2

    bt1 = -at_f1_2[type] * sinphi;
    bt2 = -2.0 * at_f2_2[type] * sin(2.0*phi);
    bt3 = -3.0 * at_f3_2[type] * sin(3.0*phi);
    sumbtf = bt1 + bt2 + bt3;

    for (i = 0; i < 4; i++)
      for (j = 0; j < 3; j++)
	fabcd[i][j] += da*sumbtf*dphidr[i][j] + sumbte*dthetadr[1][i][j];

    // angle/angle/torsion coupling

    da1 = acos(costh12) - aat_theta0_1[type];
    da2 = acos(costh23) - aat_theta0_2[type];
          
    // energy

    if (eflag) energy += rfactor * aat_k[type]*da1*da2*cosphi;

    // force

    for (i = 0; i < 4; i++)
      for (j = 0; j < 3; j++)
	fabcd[i][j] -= aat_k[type] * 
	  (cosphi * (da2*dthetadr[0][i][j] - da1*dthetadr[1][i][j]) +
	   sinphi * da1*da2*dphidr[i][j]);

    // bond1/bond3 coupling

    if (fabs(bb13t_k[type]) > SMALL) {

      r1_0 = bb13t_r10[type];
      r3_0 = bb13t_r30[type];
      dr1 = r1 - r1_0;
      dr2 = r3 - r3_0;
      tk1 = -bb13t_k[type] * dr1 / r3;
      tk2 = -bb13t_k[type] * dr2 / r1;

      if (eflag) energy += rfactor * bb13t_k[type]*dr1*dr2;
        
      fabcd[0][0] += tk2 * delx1;
      fabcd[0][1] += tk2 * dely1;
      fabcd[0][2] += tk2 * delz1;

      fabcd[1][0] -= tk2 * delx1;
      fabcd[1][1] -= tk2 * dely1;
      fabcd[1][2] -= tk2 * delz1;
        
      fabcd[2][0] -= tk1 * delx3;
      fabcd[2][1] -= tk1 * dely3;
      fabcd[2][2] -= tk1 * delz3;

      fabcd[3][0] += tk1 * delx3;
      fabcd[3][1] += tk1 * dely3;
      fabcd[3][2] += tk1 * delz3;
    }

    // apply force to each of 4 atoms

    if (newton_bond || i1 < nlocal) {
      f[i1][0] += fabcd[0][0];
      f[i1][1] += fabcd[0][1];
      f[i1][2] += fabcd[0][2];
    }

    if (newton_bond || i2 < nlocal) {
      f[i2][0] += fabcd[1][0];
      f[i2][1] += fabcd[1][1];
      f[i2][2] += fabcd[1][2];
    }

    if (newton_bond || i3 < nlocal) {
      f[i3][0] += fabcd[2][0];
      f[i3][1] += fabcd[2][1];
      f[i3][2] += fabcd[2][2];
    }

    if (newton_bond || i4 < nlocal) {
      f[i4][0] += fabcd[3][0];
      f[i4][1] += fabcd[3][1];
      f[i4][2] += fabcd[3][2];
    }

    // virial contribution

    if (vflag) {
      vx1 = fabcd[0][0];
      vx2 = fabcd[2][0] + fabcd[3][0];
      vx3 = fabcd[3][0];

      vy1 = fabcd[0][1];
      vy2 = fabcd[2][1] + fabcd[3][1];
      vy3 = fabcd[3][1];

      vz1 = fabcd[0][2];
      vz2 = fabcd[2][2] + fabcd[3][2];
      vz3 = fabcd[3][2];

      virial[0] += rfactor * (delx1*vx1 + delx2*vx2 + delx3*vx3);
      virial[1] += rfactor * (dely1*vy1 + dely2*vy2 + dely3*vy3);
      virial[2] += rfactor * (delz1*vz1 + delz2*vz2 + delz3*vz3);
      virial[3] += rfactor * (delx1*vy1 + delx2*vy2 + delx3*vy3);
      virial[4] += rfactor * (delx1*vz1 + delx2*vz2 + delx3*vz3);
      virial[5] += rfactor * (dely1*vz1 + dely2*vz2 + dely3*vz3);
    }
  }
}

/* ---------------------------------------------------------------------- */

void DihedralClass2::allocate()
{
  allocated = 1;
  int n = atom->ndihedraltypes;

  k1 = (double *) memory->smalloc((n+1)*sizeof(double),"dihedral:k1");
  k2 = (double *) memory->smalloc((n+1)*sizeof(double),"dihedral:k2");
  k3 = (double *) memory->smalloc((n+1)*sizeof(double),"dihedral:k3");
  phi1 = (double *) memory->smalloc((n+1)*sizeof(double),"dihedral:phi1");
  phi2 = (double *) memory->smalloc((n+1)*sizeof(double),"dihedral:phi2");
  phi3 = (double *) memory->smalloc((n+1)*sizeof(double),"dihedral:phi3");

  mbt_f1 = (double *) memory->smalloc((n+1)*sizeof(double),"dihedral:mbt_f1");
  mbt_f2 = (double *) memory->smalloc((n+1)*sizeof(double),"dihedral:mbt_f2");
  mbt_f3 = (double *) memory->smalloc((n+1)*sizeof(double),"dihedral:mbt_f3");
  mbt_r0 = (double *) memory->smalloc((n+1)*sizeof(double),"dihedral:mbt_r0");

  ebt_f1_1 = (double *) 
    memory->smalloc((n+1)*sizeof(double),"dihedral:ebt_f1_1");
  ebt_f2_1 = (double *) 
    memory->smalloc((n+1)*sizeof(double),"dihedral:ebt_f2_1");
  ebt_f3_1 = (double *) 
    memory->smalloc((n+1)*sizeof(double),"dihedral:ebt_f3_1");
  ebt_r0_1 = (double *) 
    memory->smalloc((n+1)*sizeof(double),"dihedral:ebt_r0_1");

  ebt_f1_2 = (double *) 
    memory->smalloc((n+1)*sizeof(double),"dihedral:ebt_f1_2");
  ebt_f2_2 = (double *) 
    memory->smalloc((n+1)*sizeof(double),"dihedral:ebt_f2_2");
  ebt_f3_2 = (double *) 
    memory->smalloc((n+1)*sizeof(double),"dihedral:ebt_f3_2");
  ebt_r0_2 = (double *) 
    memory->smalloc((n+1)*sizeof(double),"dihedral:ebt_r0_2");

  at_f1_1 = (double *) 
    memory->smalloc((n+1)*sizeof(double),"dihedral:at_f1_1");
  at_f2_1 = (double *) 
    memory->smalloc((n+1)*sizeof(double),"dihedral:at_f2_1");
  at_f3_1 = (double *) 
    memory->smalloc((n+1)*sizeof(double),"dihedral:at_f3_1");
  at_theta0_1 = (double *) 
    memory->smalloc((n+1)*sizeof(double),"dihedral:at_theta0_1");

  at_f1_2 = (double *) 
    memory->smalloc((n+1)*sizeof(double),"dihedral:at_f1_2");
  at_f2_2 = (double *) 
    memory->smalloc((n+1)*sizeof(double),"dihedral:at_f2_2");
  at_f3_2 = (double *) 
    memory->smalloc((n+1)*sizeof(double),"dihedral:at_f3_2");
  at_theta0_2 = (double *) 
    memory->smalloc((n+1)*sizeof(double),"dihedral:at_theta0_2");

  aat_k = (double *) memory->smalloc((n+1)*sizeof(double),"dihedral:aat_k");
  aat_theta0_1 = (double *) 
    memory->smalloc((n+1)*sizeof(double),"dihedral:aat_theta0_1");
  aat_theta0_2 = (double *) 
    memory->smalloc((n+1)*sizeof(double),"dihedral:aat_theta0_2");

  bb13t_k = (double *) 
    memory->smalloc((n+1)*sizeof(double),"dihedral:bb13t_k");
  bb13t_r10 = (double *) 
    memory->smalloc((n+1)*sizeof(double),"dihedral:bb13t_r10");
  bb13t_r30 = (double *) 
    memory->smalloc((n+1)*sizeof(double),"dihedral:bb13t_r30");

  setflag = (int *) memory->smalloc((n+1)*sizeof(int),"dihedral:setflag");
  setflag_d = (int *) 
    memory->smalloc((n+1)*sizeof(int),"dihedral:setflag_d");
  setflag_mbt = (int *) 
    memory->smalloc((n+1)*sizeof(int),"dihedral:setflag_mbt");
  setflag_ebt = (int *) 
    memory->smalloc((n+1)*sizeof(int),"dihedral:setflag_ebt");
  setflag_at = (int *) 
    memory->smalloc((n+1)*sizeof(int),"dihedral:setflag_at");
  setflag_aat = (int *) 
    memory->smalloc((n+1)*sizeof(int),"dihedral:setflag_aat");
  setflag_bb13t = (int *) 
    memory->smalloc((n+1)*sizeof(int),"dihedral:setflag_bb13t");
  for (int i = 1; i <= n; i++)
    setflag[i] = setflag_d[i] = setflag_mbt[i] = setflag_ebt[i] = 
      setflag_at[i] = setflag_aat[i] = setflag_bb13t[i] = 0;
}

/* ----------------------------------------------------------------------
   set coeffs for one or more types
   which = 0 -> Dihedral coeffs
   which = 1 -> MiddleBondTorsion coeffs
   which = 2 -> EndBondTorsion coeffs
   which = 3 -> AngleTorsion coeffs
   which = 4 -> AngleAngleTorsion coeffs
   which = 5 -> BondBond13Torsion coeffs
------------------------------------------------------------------------- */

void DihedralClass2::coeff(int which, int narg, char **arg)
{
  if (which < 0 || which > 5)
    error->all("Invalid coeffs for this dihedral style");
  if (!allocated) allocate();

  int ilo,ihi;
  force->bounds(arg[0],atom->ndihedraltypes,ilo,ihi);

  int count = 0;

  if (which == 0) {
    if (narg != 7) error->all("Incorrect args for dihedral coefficients");

    double k1_one = atof(arg[1]);
    double phi1_one = atof(arg[2]);
    double k2_one = atof(arg[3]);
    double phi2_one = atof(arg[4]);
    double k3_one = atof(arg[5]);
    double phi3_one = atof(arg[6]);
    
    // convert phi's from degrees to radians

    for (int i = ilo; i <= ihi; i++) {
      k1[i] = k1_one;
      phi1[i] = phi1_one/180.0 * PI;
      k2[i] = k2_one;
      phi2[i] = phi2_one/180.0 * PI;
      k3[i] = k3_one;
      phi3[i] = phi3_one/180.0 * PI;
      setflag_d[i] = 1;
      count++;
    }
  }

  if (which == 1) {
    if (narg != 5) error->all("Incorrect args for dihedral coefficients");

    double f1_one = atof(arg[1]);
    double f2_one = atof(arg[2]);
    double f3_one = atof(arg[3]);
    double r0_one = atof(arg[4]);
    
    for (int i = ilo; i <= ihi; i++) {
      mbt_f1[i] = f1_one;
      mbt_f2[i] = f2_one;
      mbt_f3[i] = f3_one;
      mbt_r0[i] = r0_one;
      setflag_mbt[i] = 1;
      count++;
    }
  }

  if (which == 2) {
    if (narg != 9) error->all("Incorrect args for dihedral coefficients");

    double f1_1_one = atof(arg[1]);
    double f2_1_one = atof(arg[2]);
    double f3_1_one = atof(arg[3]);
    double f1_2_one = atof(arg[4]);
    double f2_2_one = atof(arg[5]);
    double f3_2_one = atof(arg[6]);
    double r0_1_one = atof(arg[7]);
    double r0_2_one = atof(arg[8]);
    
    for (int i = ilo; i <= ihi; i++) {
      ebt_f1_1[i] = f1_1_one;
      ebt_f2_1[i] = f2_1_one;
      ebt_f3_1[i] = f3_1_one;
      ebt_f1_2[i] = f1_2_one;
      ebt_f2_2[i] = f2_2_one;
      ebt_f3_2[i] = f3_2_one;
      ebt_r0_1[i] = r0_1_one;
      ebt_r0_2[i] = r0_2_one;
      setflag_ebt[i] = 1;
      count++;
    }
  }

  if (which == 3) {
    if (narg != 9) error->all("Incorrect args for dihedral coefficients");

    double f1_1_one = atof(arg[1]);
    double f2_1_one = atof(arg[2]);
    double f3_1_one = atof(arg[3]);
    double f1_2_one = atof(arg[4]);
    double f2_2_one = atof(arg[5]);
    double f3_2_one = atof(arg[6]);
    double theta0_1_one = atof(arg[7]);
    double theta0_2_one = atof(arg[8]);

    // convert theta0's from degrees to radians
    
    for (int i = ilo; i <= ihi; i++) {
      at_f1_1[i] = f1_1_one;
      at_f2_1[i] = f2_1_one;
      at_f3_1[i] = f3_1_one;
      at_f1_2[i] = f1_2_one;
      at_f2_2[i] = f2_2_one;
      at_f3_2[i] = f3_2_one;
      at_theta0_1[i] = theta0_1_one/180.0 * PI;
      at_theta0_2[i] = theta0_2_one/180.0 * PI;
      setflag_at[i] = 1;
      count++;
    }
  }

  if (which == 4) {
    if (narg != 4) error->all("Incorrect args for dihedral coefficients");

    double k_one = atof(arg[1]);
    double theta0_1_one = atof(arg[2]);
    double theta0_2_one = atof(arg[3]);

    // convert theta0's from degrees to radians
    
    for (int i = ilo; i <= ihi; i++) {
      aat_k[i] = k_one;
      aat_theta0_1[i] = theta0_1_one/180.0 * PI;
      aat_theta0_2[i] = theta0_2_one/180.0 * PI;
      setflag_aat[i] = 1;
      count++;
    }
  }

  if (which == 5) {
    if (narg != 4) error->all("Incorrect args for dihedral coefficients");

    double k_one = atof(arg[1]);
    double r10_one = atof(arg[2]);
    double r30_one = atof(arg[3]);
    
    for (int i = ilo; i <= ihi; i++) {
      bb13t_k[i] = k_one;
      bb13t_r10[i] = r10_one;
      bb13t_r30[i] = r30_one;
      setflag_bb13t[i] = 1;
      count++;
    }
  }

  if (count == 0) error->all("Incorrect args for dihedral coefficients");

  for (int i = ilo; i <= ihi; i++)
    if (setflag_d[i] == 1 && setflag_mbt[i] == 1 && setflag_ebt[i] == 1 &&
	setflag_at[i] == 1 && setflag_aat[i] == 1 && setflag_bb13t[i] == 1)
      setflag[i] = 1;
}

/* ----------------------------------------------------------------------
   proc 0 writes out coeffs to restart file 
------------------------------------------------------------------------- */

void DihedralClass2::write_restart(FILE *fp)
{
  fwrite(&k1[1],sizeof(double),atom->ndihedraltypes,fp);
  fwrite(&k2[1],sizeof(double),atom->ndihedraltypes,fp);
  fwrite(&k3[1],sizeof(double),atom->ndihedraltypes,fp);
  fwrite(&phi1[1],sizeof(double),atom->ndihedraltypes,fp);
  fwrite(&phi2[1],sizeof(double),atom->ndihedraltypes,fp);
  fwrite(&phi3[1],sizeof(double),atom->ndihedraltypes,fp);

  fwrite(&mbt_f1[1],sizeof(double),atom->ndihedraltypes,fp);
  fwrite(&mbt_f2[1],sizeof(double),atom->ndihedraltypes,fp);
  fwrite(&mbt_f3[1],sizeof(double),atom->ndihedraltypes,fp);
  fwrite(&mbt_r0[1],sizeof(double),atom->ndihedraltypes,fp);

  fwrite(&ebt_f1_1[1],sizeof(double),atom->ndihedraltypes,fp);
  fwrite(&ebt_f2_1[1],sizeof(double),atom->ndihedraltypes,fp);
  fwrite(&ebt_f3_1[1],sizeof(double),atom->ndihedraltypes,fp);
  fwrite(&ebt_r0_1[1],sizeof(double),atom->ndihedraltypes,fp);

  fwrite(&ebt_f1_2[1],sizeof(double),atom->ndihedraltypes,fp);
  fwrite(&ebt_f2_2[1],sizeof(double),atom->ndihedraltypes,fp);
  fwrite(&ebt_f3_2[1],sizeof(double),atom->ndihedraltypes,fp);
  fwrite(&ebt_r0_2[1],sizeof(double),atom->ndihedraltypes,fp);

  fwrite(&at_f1_1[1],sizeof(double),atom->ndihedraltypes,fp);
  fwrite(&at_f2_1[1],sizeof(double),atom->ndihedraltypes,fp);
  fwrite(&at_f3_1[1],sizeof(double),atom->ndihedraltypes,fp);
  fwrite(&at_theta0_1[1],sizeof(double),atom->ndihedraltypes,fp);

  fwrite(&at_f1_2[1],sizeof(double),atom->ndihedraltypes,fp);
  fwrite(&at_f2_2[1],sizeof(double),atom->ndihedraltypes,fp);
  fwrite(&at_f3_2[1],sizeof(double),atom->ndihedraltypes,fp);
  fwrite(&at_theta0_2[1],sizeof(double),atom->ndihedraltypes,fp);

  fwrite(&aat_k[1],sizeof(double),atom->ndihedraltypes,fp);
  fwrite(&aat_theta0_1[1],sizeof(double),atom->ndihedraltypes,fp);
  fwrite(&aat_theta0_2[1],sizeof(double),atom->ndihedraltypes,fp);

  fwrite(&bb13t_k[1],sizeof(double),atom->ndihedraltypes,fp);
  fwrite(&bb13t_r10[1],sizeof(double),atom->ndihedraltypes,fp);
  fwrite(&bb13t_r30[1],sizeof(double),atom->ndihedraltypes,fp);
}

/* ----------------------------------------------------------------------
   proc 0 reads coeffs from restart file, bcasts them 
------------------------------------------------------------------------- */

void DihedralClass2::read_restart(FILE *fp)
{
  allocate();

  if (comm->me == 0) {
    fread(&k1[1],sizeof(double),atom->ndihedraltypes,fp);
    fread(&k2[1],sizeof(double),atom->ndihedraltypes,fp);
    fread(&k3[1],sizeof(double),atom->ndihedraltypes,fp);
    fread(&phi1[1],sizeof(double),atom->ndihedraltypes,fp);
    fread(&phi2[1],sizeof(double),atom->ndihedraltypes,fp);
    fread(&phi3[1],sizeof(double),atom->ndihedraltypes,fp);

    fread(&mbt_f1[1],sizeof(double),atom->ndihedraltypes,fp);
    fread(&mbt_f2[1],sizeof(double),atom->ndihedraltypes,fp);
    fread(&mbt_f3[1],sizeof(double),atom->ndihedraltypes,fp);
    fread(&mbt_r0[1],sizeof(double),atom->ndihedraltypes,fp);

    fread(&ebt_f1_1[1],sizeof(double),atom->ndihedraltypes,fp);
    fread(&ebt_f2_1[1],sizeof(double),atom->ndihedraltypes,fp);
    fread(&ebt_f3_1[1],sizeof(double),atom->ndihedraltypes,fp);
    fread(&ebt_r0_1[1],sizeof(double),atom->ndihedraltypes,fp);

    fread(&ebt_f1_2[1],sizeof(double),atom->ndihedraltypes,fp);
    fread(&ebt_f2_2[1],sizeof(double),atom->ndihedraltypes,fp);
    fread(&ebt_f3_2[1],sizeof(double),atom->ndihedraltypes,fp);
    fread(&ebt_r0_2[1],sizeof(double),atom->ndihedraltypes,fp);

    fread(&at_f1_1[1],sizeof(double),atom->ndihedraltypes,fp);
    fread(&at_f2_1[1],sizeof(double),atom->ndihedraltypes,fp);
    fread(&at_f3_1[1],sizeof(double),atom->ndihedraltypes,fp);
    fread(&at_theta0_1[1],sizeof(double),atom->ndihedraltypes,fp);

    fread(&at_f1_2[1],sizeof(double),atom->ndihedraltypes,fp);
    fread(&at_f2_2[1],sizeof(double),atom->ndihedraltypes,fp);
    fread(&at_f3_2[1],sizeof(double),atom->ndihedraltypes,fp);
    fread(&at_theta0_2[1],sizeof(double),atom->ndihedraltypes,fp);

    fread(&aat_k[1],sizeof(double),atom->ndihedraltypes,fp);
    fread(&aat_theta0_1[1],sizeof(double),atom->ndihedraltypes,fp);
    fread(&aat_theta0_2[1],sizeof(double),atom->ndihedraltypes,fp);

    fread(&bb13t_k[1],sizeof(double),atom->ndihedraltypes,fp);
    fread(&bb13t_r10[1],sizeof(double),atom->ndihedraltypes,fp);
    fread(&bb13t_r30[1],sizeof(double),atom->ndihedraltypes,fp);
  }

  MPI_Bcast(&k1[1],atom->ndihedraltypes,MPI_DOUBLE,0,world);
  MPI_Bcast(&k2[1],atom->ndihedraltypes,MPI_DOUBLE,0,world);
  MPI_Bcast(&k3[1],atom->ndihedraltypes,MPI_DOUBLE,0,world);
  MPI_Bcast(&phi1[1],atom->ndihedraltypes,MPI_DOUBLE,0,world);
  MPI_Bcast(&phi2[1],atom->ndihedraltypes,MPI_DOUBLE,0,world);
  MPI_Bcast(&phi3[1],atom->ndihedraltypes,MPI_DOUBLE,0,world);

  MPI_Bcast(&mbt_f1[1],atom->ndihedraltypes,MPI_DOUBLE,0,world);
  MPI_Bcast(&mbt_f2[1],atom->ndihedraltypes,MPI_DOUBLE,0,world);
  MPI_Bcast(&mbt_f3[1],atom->ndihedraltypes,MPI_DOUBLE,0,world);
  MPI_Bcast(&mbt_r0[1],atom->ndihedraltypes,MPI_DOUBLE,0,world);

  MPI_Bcast(&ebt_f1_1[1],atom->ndihedraltypes,MPI_DOUBLE,0,world);
  MPI_Bcast(&ebt_f2_1[1],atom->ndihedraltypes,MPI_DOUBLE,0,world);
  MPI_Bcast(&ebt_f3_1[1],atom->ndihedraltypes,MPI_DOUBLE,0,world);
  MPI_Bcast(&ebt_r0_1[1],atom->ndihedraltypes,MPI_DOUBLE,0,world);

  MPI_Bcast(&ebt_f1_2[1],atom->ndihedraltypes,MPI_DOUBLE,0,world);
  MPI_Bcast(&ebt_f2_2[1],atom->ndihedraltypes,MPI_DOUBLE,0,world);
  MPI_Bcast(&ebt_f3_2[1],atom->ndihedraltypes,MPI_DOUBLE,0,world);
  MPI_Bcast(&ebt_r0_2[1],atom->ndihedraltypes,MPI_DOUBLE,0,world);

  MPI_Bcast(&at_f1_1[1],atom->ndihedraltypes,MPI_DOUBLE,0,world);
  MPI_Bcast(&at_f2_1[1],atom->ndihedraltypes,MPI_DOUBLE,0,world);
  MPI_Bcast(&at_f3_1[1],atom->ndihedraltypes,MPI_DOUBLE,0,world);
  MPI_Bcast(&at_theta0_1[1],atom->ndihedraltypes,MPI_DOUBLE,0,world);

  MPI_Bcast(&at_f1_2[1],atom->ndihedraltypes,MPI_DOUBLE,0,world);
  MPI_Bcast(&at_f2_2[1],atom->ndihedraltypes,MPI_DOUBLE,0,world);
  MPI_Bcast(&at_f3_2[1],atom->ndihedraltypes,MPI_DOUBLE,0,world);
  MPI_Bcast(&at_theta0_2[1],atom->ndihedraltypes,MPI_DOUBLE,0,world);

  MPI_Bcast(&aat_k[1],atom->ndihedraltypes,MPI_DOUBLE,0,world);
  MPI_Bcast(&aat_theta0_1[1],atom->ndihedraltypes,MPI_DOUBLE,0,world);
  MPI_Bcast(&aat_theta0_2[1],atom->ndihedraltypes,MPI_DOUBLE,0,world);

  MPI_Bcast(&bb13t_k[1],atom->ndihedraltypes,MPI_DOUBLE,0,world);
  MPI_Bcast(&bb13t_r10[1],atom->ndihedraltypes,MPI_DOUBLE,0,world);
  MPI_Bcast(&bb13t_r30[1],atom->ndihedraltypes,MPI_DOUBLE,0,world);

  for (int i = 1; i <= atom->ndihedraltypes; i++) setflag[i] = 1;
}
