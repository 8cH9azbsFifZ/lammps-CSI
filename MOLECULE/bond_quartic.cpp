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
   Contributing authors: Chris Lorenz and Mark Stevens (SNL)
------------------------------------------------------------------------- */

#include "math.h"
#include "stdlib.h"
#include "bond_quartic.h"
#include "atom.h"
#include "neighbor.h"
#include "domain.h"
#include "comm.h"
#include "update.h"
#include "force.h"
#include "pair.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

BondQuartic::BondQuartic(LAMMPS *lmp) : Bond(lmp)
{
  TWO_1_3 = pow(2.0,(1.0/3.0));
}

/* ---------------------------------------------------------------------- */

BondQuartic::~BondQuartic()
{
  if (allocated) {
    memory->sfree(setflag);
    memory->sfree(k);
    memory->sfree(b1);
    memory->sfree(b2);
    memory->sfree(rc);
    memory->sfree(u0);
  }
}

/* ---------------------------------------------------------------------- */

void BondQuartic::compute(int eflag, int vflag)
{
  int i1,i2,n,m,type,factor,itype,jtype;
  double delx,dely,delz,r,rsq,dr,r2,ra,rb,fforce,sr2,sr6,rfactor;
  Pair::One one;

  energy = 0.0;
  eng_vdwl = 0.0;
  if (vflag) for (n = 0; n < 6; n++) virial[n] = 0.0;

  double **cutsq = force->pair->cutsq;
  double **x = atom->x;
  double **f = atom->f;
  int **bondlist = neighbor->bondlist;
  int nbondlist = neighbor->nbondlist;
  int nlocal = atom->nlocal;
  int newton_bond = force->newton_bond;

  for (n = 0; n < nbondlist; n++) {

    // skip bond if already broken

    if (bondlist[n][2] <= 0) continue;

    i1 = bondlist[n][0];
    i2 = bondlist[n][1];
    type = bondlist[n][2];

    if (newton_bond) factor = 2;
    else {
      factor = 0;
      if (i1 < nlocal) factor++;
      if (i2 < nlocal) factor++;
    }
    rfactor = 0.5*factor;

    delx = x[i1][0] - x[i2][0];
    dely = x[i1][1] - x[i2][1];
    delz = x[i1][2] - x[i2][2];
    domain->minimum_image(delx,dely,delz);

    rsq = delx*delx + dely*dely + delz*delz;

    // if bond breaks, set type to 0
    //   both in temporary bondlist and permanent bond_type
    // if this proc owns both atoms,
    //   negate bond_type twice if other atom stores it
    // if other proc owns 2nd atom, other proc will also break bond

    if (rsq > rc[type]*rc[type]) {
      bondlist[n][2] = 0;
      for (m = 0; m < atom->num_bond[i1]; m++)
	if (atom->bond_atom[i1][m] == atom->tag[i2])
	  atom->bond_type[i1][m] = 0;
      if (i2 < atom->nlocal)
	for (m = 0; m < atom->num_bond[i2]; m++)
	  if (atom->bond_atom[i2][m] == atom->tag[i1])
	    atom->bond_type[i2][m] = 0;
      continue;
    }

    // subtract out pairwise contribution from 2 atoms via pair->single()
    // required since special_bond = 1,1,1

    itype = atom->type[i1];
    jtype = atom->type[i2];

    if (rsq < cutsq[itype][jtype]) {
      force->pair->single(i1,i2,itype,jtype,rsq,1.0,1.0,eflag,one);
      fforce = -one.fforce;
      if (eflag) eng_vdwl -= one.eng_vdwl + one.eng_coul;
    } else fforce = 0.0;

    // quartic bond
    // 1st portion is from quartic term
    // 2nd portion is from LJ term cut at 2^(1/6) with eps = sigma = 1.0

    r = sqrt(rsq);
    dr = r - rc[type];
    r2 = dr*dr;
    ra = dr - b1[type];
    rb = dr - b2[type];
    fforce += -k[type]/r * (r2*(ra+rb) + 2.0*dr*ra*rb);
    
    if (rsq < TWO_1_3) {
      sr2 = 1.0/rsq;
      sr6 = sr2*sr2*sr2;
      fforce += 48.0*sr6*(sr6-0.5)/rsq;
    }

    if (eflag) {
      energy += rfactor*(k[type]*r2*ra*rb + u0[type]);
      if (rsq < TWO_1_3) energy += rfactor * (4.0*sr6*(sr6-1.0) + 1.0);
    }

    // apply force to each of 2 atoms

    if (newton_bond || i1 < nlocal) {
      f[i1][0] += delx*fforce;
      f[i1][1] += dely*fforce;
      f[i1][2] += delz*fforce;
    }

    if (newton_bond || i2 < nlocal) {
      f[i2][0] -= delx*fforce;
      f[i2][1] -= dely*fforce;
      f[i2][2] -= delz*fforce;
    }

    // virial contribution

    if (vflag) {
      virial[0] += rfactor * delx*delx*fforce;
      virial[1] += rfactor * dely*dely*fforce;
      virial[2] += rfactor * delz*delz*fforce;
      virial[3] += rfactor * delx*dely*fforce;
      virial[4] += rfactor * delx*delz*fforce;
      virial[5] += rfactor * dely*delz*fforce;
    }
  }
}

/* ---------------------------------------------------------------------- */

void BondQuartic::allocate()
{
  allocated = 1;
  int n = atom->nbondtypes;

  k = (double *) memory->smalloc((n+1)*sizeof(double),"bond:k");
  b1 = (double *) memory->smalloc((n+1)*sizeof(double),"bond:b1");
  b2 = (double *) memory->smalloc((n+1)*sizeof(double),"bond:b2");
  rc = (double *) memory->smalloc((n+1)*sizeof(double),"bond:rc");
  u0 = (double *) memory->smalloc((n+1)*sizeof(double),"bond:u0");

  setflag = (int *) memory->smalloc((n+1)*sizeof(int),"bond:setflag");
  for (int i = 1; i <= n; i++) setflag[i] = 0;
}

/* ----------------------------------------------------------------------
   set coeffs for one or more types
------------------------------------------------------------------------- */

void BondQuartic::coeff(int narg, char **arg)
{
  if (narg != 6) error->all("Incorrect args for bond coefficients");
  if (!allocated) allocate();

  int ilo,ihi;
  force->bounds(arg[0],atom->nbondtypes,ilo,ihi);

  double k_one = atof(arg[1]);
  double b1_one = atof(arg[2]);
  double b2_one = atof(arg[3]);
  double rc_one = atof(arg[4]);
  double u0_one = atof(arg[5]);

  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
    k[i] = k_one;
    b1[i] = b1_one;
    b2[i] = b2_one;
    rc[i] = rc_one;
    u0[i] = u0_one;
    setflag[i] = 1;
    count++;
  }

  if (count == 0) error->all("Incorrect args for bond coefficients");
}

/* ----------------------------------------------------------------------
   check if pair defined and special_bond settings are valid
------------------------------------------------------------------------- */

void BondQuartic::init_style()
{
  if (force->pair == NULL || force->pair->single_enable == 0)
    error->all("Pair style does not support bond_style quartic");
  if (force->angle)
    error->all("Bond style quartic cannot be used with 3,4-body interactions");
  if (force->dihedral)
    error->all("Bond style quartic cannot be used with 3,4-body interactions");
  if (force->improper)
    error->all("Bond style quartic cannot be used with 3,4-body interactions");

  // special bonds must be 1 1 1

  if (force->special_lj[1] != 1.0 || force->special_lj[2] != 1.0 ||
      force->special_lj[3] != 1.0)
    error->all("Must use special bonds = 1,1,1 with bond style quartic");
}

/* ----------------------------------------------------------------------
   return an equilbrium bond length 
------------------------------------------------------------------------- */

double BondQuartic::equilibrium_distance(int i)
{
  return 0.97;
}

/* ----------------------------------------------------------------------
   proc 0 writes out coeffs to restart file 
------------------------------------------------------------------------- */

void BondQuartic::write_restart(FILE *fp)
{
  fwrite(&k[1],sizeof(double),atom->nbondtypes,fp);
  fwrite(&b1[1],sizeof(double),atom->nbondtypes,fp);
  fwrite(&b2[1],sizeof(double),atom->nbondtypes,fp);
  fwrite(&rc[1],sizeof(double),atom->nbondtypes,fp);
  fwrite(&u0[1],sizeof(double),atom->nbondtypes,fp);
}

/* ----------------------------------------------------------------------
   proc 0 reads coeffs from restart file, bcasts them 
------------------------------------------------------------------------- */

void BondQuartic::read_restart(FILE *fp)
{
  allocate();

  if (comm->me == 0) {
    fread(&k[1],sizeof(double),atom->nbondtypes,fp);
    fread(&b1[1],sizeof(double),atom->nbondtypes,fp);
    fread(&b2[1],sizeof(double),atom->nbondtypes,fp);
    fread(&rc[1],sizeof(double),atom->nbondtypes,fp);
    fread(&u0[1],sizeof(double),atom->nbondtypes,fp);
  }
  MPI_Bcast(&k[1],atom->nbondtypes,MPI_DOUBLE,0,world);
  MPI_Bcast(&b1[1],atom->nbondtypes,MPI_DOUBLE,0,world);
  MPI_Bcast(&b2[1],atom->nbondtypes,MPI_DOUBLE,0,world);
  MPI_Bcast(&rc[1],atom->nbondtypes,MPI_DOUBLE,0,world);
  MPI_Bcast(&u0[1],atom->nbondtypes,MPI_DOUBLE,0,world);

  for (int i = 1; i <= atom->nbondtypes; i++) setflag[i] = 1;
}


/* ---------------------------------------------------------------------- */

void BondQuartic::single(int type, double rsq, int i, int j, double rfactor,
			 int eflag, double &fforce, double &eng)
{
  double r,dr,r2,ra,rb,sr2,sr6;

  fforce = eng = 0.0;
  if (type <= 0) return;

  // subtract out pairwise contribution from 2 atoms via pair->single()
  // required since special_bond = 1,1,1

  int itype = atom->type[i];
  int jtype = atom->type[j];
  
  if (rsq < force->pair->cutsq[itype][jtype]) {
    Pair::One one;
    force->pair->single(i,j,itype,jtype,rsq,1.0,1.0,eflag,one);
    fforce = -one.fforce;
    if (eflag) eng = -one.eng_coul - one.eng_vdwl;
  }

  // quartic bond
  // 1st portion is from quartic term
  // 2nd portion is from LJ term cut at 2^(1/6) with eps = sigma = 1.0

  r = sqrt(rsq);
  dr = r - rc[type];
  r2 = dr*dr;
  ra = dr - b1[type];
  rb = dr - b2[type];
  fforce += -k[type]/r * (r2*(ra+rb) + 2.0*dr*ra*rb);
  
  if (rsq < TWO_1_3) {
    sr2 = 1.0/rsq;
    sr6 = sr2*sr2*sr2;
    fforce += 48.0*sr6*(sr6-0.5)/rsq;
  }
    
  if (eflag) {
    eng += rfactor*(k[type]*r2*ra*rb + u0[type]);
    if (rsq < TWO_1_3) eng += rfactor * (4.0*sr6*(sr6-1.0) + 1.0);
  }
}
