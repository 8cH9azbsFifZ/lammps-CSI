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
   Contributing author: Paul Crozier (SNL)
------------------------------------------------------------------------- */

#include "math.h"
#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "pair_lj_charmm_coul_long.h"
#include "atom.h"
#include "comm.h"
#include "force.h"
#include "kspace.h"
#include "update.h"
#include "integrate.h"
#include "respa.h"
#include "memory.h"
#include "neighbor.h"
#include "error.h"

using namespace LAMMPS_NS;

#define MIN(a,b) ((a) < (b) ? (a) : (b))
#define MAX(a,b) ((a) > (b) ? (a) : (b))

#define EWALD_F   1.12837917
#define EWALD_P   0.3275911
#define A1        0.254829592
#define A2       -0.284496736
#define A3        1.421413741
#define A4       -1.453152027
#define A5        1.061405429

/* ---------------------------------------------------------------------- */

PairLJCharmmCoulLong::PairLJCharmmCoulLong(LAMMPS *lmp) : Pair(lmp)
{
  respa_enable = 1;
  ftable = NULL;
}

/* ---------------------------------------------------------------------- */

PairLJCharmmCoulLong::~PairLJCharmmCoulLong()
{
  if (allocated) {
    memory->destroy_2d_int_array(setflag);
    memory->destroy_2d_double_array(cutsq);

    memory->destroy_2d_double_array(epsilon);
    memory->destroy_2d_double_array(sigma);
    memory->destroy_2d_double_array(eps14);
    memory->destroy_2d_double_array(sigma14);
    memory->destroy_2d_double_array(lj1);
    memory->destroy_2d_double_array(lj2);
    memory->destroy_2d_double_array(lj3);
    memory->destroy_2d_double_array(lj4);
    memory->destroy_2d_double_array(lj14_1);
    memory->destroy_2d_double_array(lj14_2);
    memory->destroy_2d_double_array(lj14_3);
    memory->destroy_2d_double_array(lj14_4);
  }
  if (ftable) free_tables();
}

/* ---------------------------------------------------------------------- */

void PairLJCharmmCoulLong::compute(int eflag, int vflag)
{
  int i,j,k,numneigh,itype,jtype,itable;
  double qtmp,xtmp,ytmp,ztmp,delx,dely,delz,fraction,table;
  double r,r2inv,r6inv,forcecoul,forcelj,fforce,factor_coul,factor_lj;
  double grij,expm2,prefactor,t,erfc;
  double factor,phicoul,philj,switch1,switch2;
  int *neighs;
  double **f;
  float rsq;
  int *int_rsq = (int *) &rsq;

  eng_vdwl = eng_coul = 0.0;
  if (vflag) for (i = 0; i < 6; i++) virial[i] = 0.0;

  if (vflag == 2) f = update->f_pair;
  else f = atom->f;
  double **x = atom->x;
  double *q = atom->q;
  int *type = atom->type;
  int nlocal = atom->nlocal;
  int nall = atom->nlocal + atom->nghost;
  double *special_coul = force->special_coul;
  double *special_lj = force->special_lj;
  int newton_pair = force->newton_pair;
  double qqrd2e = force->qqrd2e;

  // loop over neighbors of my atoms

  for (i = 0; i < nlocal; i++) {
    qtmp = q[i];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    itype = type[i];
    neighs = neighbor->firstneigh[i];
    numneigh = neighbor->numneigh[i];

    for (k = 0; k < numneigh; k++) {
      j = neighs[k];

      if (j < nall) factor_coul = factor_lj = 1.0;
      else {
	factor_coul = special_coul[j/nall];
	factor_lj = special_lj[j/nall];
	j %= nall;
      }

      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx*delx + dely*dely + delz*delz;

      if (rsq < cut_bothsq) {
	r2inv = 1.0/rsq;

	if (rsq < cut_coulsq) {
	  if (!ncoultablebits || rsq <= tabinnersq) {
	    r = sqrtf(rsq);
	    grij = g_ewald * r;
	    expm2 = exp(-grij*grij);
	    t = 1.0 / (1.0 + EWALD_P*grij);
	    erfc = t * (A1+t*(A2+t*(A3+t*(A4+t*A5)))) * expm2;
	    prefactor = qqrd2e * qtmp*q[j]/r;
	    forcecoul = prefactor * (erfc + EWALD_F*grij*expm2);
	    if (factor_coul < 1.0) forcecoul -= (1.0-factor_coul)*prefactor;
	  } else {
	    itable = *int_rsq & ncoulmask;
	    itable >>= ncoulshiftbits;
	    fraction = (rsq - rtable[itable]) * drtable[itable];
	    table = ftable[itable] + fraction*dftable[itable];
	    forcecoul = qtmp*q[j] * table;
	    if (factor_coul < 1.0) {
	      table = ctable[itable] + fraction*dctable[itable];
	      prefactor = qtmp*q[j] * table;
	      forcecoul -= (1.0-factor_coul)*prefactor;
	    }
	  }
	} else forcecoul = 0.0;

	if (rsq < cut_ljsq) {
	  r6inv = r2inv*r2inv*r2inv;
	  jtype = type[j];
	  forcelj = r6inv * (lj1[itype][jtype]*r6inv - lj2[itype][jtype]);
	  if (rsq > cut_lj_innersq) {
	    switch1 = (cut_ljsq-rsq) * (cut_ljsq-rsq) *
	      (cut_ljsq + 2.0*rsq - 3.0*cut_lj_innersq) / denom_lj;
	    switch2 = 12.0*rsq * (cut_ljsq-rsq) * 
	      (rsq-cut_lj_innersq) / denom_lj;
	    philj = r6inv * (lj3[itype][jtype]*r6inv - lj4[itype][jtype]);
	    forcelj = forcelj*switch1 + philj*switch2;
	  }
	} else forcelj = 0.0;

	fforce = (forcecoul + factor_lj*forcelj) * r2inv;

	f[i][0] += delx*fforce;
	f[i][1] += dely*fforce;
	f[i][2] += delz*fforce;
	if (newton_pair || j < nlocal) {
	  f[j][0] -= delx*fforce;
	  f[j][1] -= dely*fforce;
	  f[j][2] -= delz*fforce;
	}

	if (eflag) {
	  if (newton_pair || j < nlocal) factor = 1.0;
	  else factor = 0.5;
	  if (rsq < cut_coulsq) {
	    if (!ncoultablebits || rsq <= tabinnersq)
	      phicoul = prefactor*erfc;
	    else {
	      table = etable[itable] + fraction*detable[itable];
	      phicoul = qtmp*q[j] * table;
	    }
	    if (factor_coul < 1.0) phicoul -= (1.0-factor_coul)*prefactor;
	    eng_coul += factor*phicoul;
	  }
	  if (rsq < cut_ljsq) {
	    philj = r6inv*(lj3[itype][jtype]*r6inv-lj4[itype][jtype]);
	    if (rsq > cut_lj_innersq) {
	      switch1 = (cut_ljsq-rsq) * (cut_ljsq-rsq) *
		(cut_ljsq + 2.0*rsq - 3.0*cut_lj_innersq) / denom_lj;
	      philj *= switch1;
	    }
	    eng_vdwl += factor*factor_lj*philj;
	  }
	}

	if (vflag == 1) {
	  if (newton_pair == 0 && j >= nlocal) fforce *= 0.5;
	  virial[0] += delx*delx*fforce;
	  virial[1] += dely*dely*fforce;
	  virial[2] += delz*delz*fforce;
	  virial[3] += delx*dely*fforce;
	  virial[4] += delx*delz*fforce;
	  virial[5] += dely*delz*fforce;
	}
      }
    }
  }
  if (vflag == 2) virial_compute();
}

/* ---------------------------------------------------------------------- */

void PairLJCharmmCoulLong::compute_inner()
{
  int i,j,k,numneigh,itype,jtype;
  double qtmp,xtmp,ytmp,ztmp,delx,dely,delz;
  double rsq,r2inv,r6inv,forcecoul,forcelj,fforce,factor_coul,factor_lj;
  double rsw;
  int *neighs;

  double **f = atom->f;
  double **x = atom->x;
  double *q = atom->q;
  int *type = atom->type;
  int nlocal = atom->nlocal;
  int nall = atom->nlocal + atom->nghost;
  double *special_coul = force->special_coul;
  double *special_lj = force->special_lj;
  int newton_pair = force->newton_pair;
  double qqrd2e = force->qqrd2e;

  double cut_out_on = cut_respa[0];
  double cut_out_off = cut_respa[1];
  
  double cut_out_diff = cut_out_off - cut_out_on;
  double cut_out_on_sq = cut_out_on*cut_out_on;
  double cut_out_off_sq = cut_out_off*cut_out_off;

  // loop over neighbors of my atoms

  for (i = 0; i < nlocal; i++) {
    qtmp = q[i];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    itype = type[i];
    neighs = neighbor->firstneigh_inner[i];
    numneigh = neighbor->numneigh_inner[i];

    for (k = 0; k < numneigh; k++) {
      j = neighs[k];

      if (j < nall) factor_coul = factor_lj = 1.0;
      else {
	factor_coul = special_coul[j/nall];
	factor_lj = special_lj[j/nall];
	j %= nall;
      }

      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx*delx + dely*dely + delz*delz;

      if (rsq < cut_out_off_sq) {
        r2inv = 1.0/rsq;
        forcecoul = qqrd2e * qtmp*q[j]*sqrt(r2inv);
        if (factor_coul < 1.0) forcecoul -= (1.0-factor_coul)*forcecoul;

	r6inv = r2inv*r2inv*r2inv;
	jtype = type[j];
	forcelj = r6inv * (lj1[itype][jtype]*r6inv - lj2[itype][jtype]);

	fforce = (forcecoul + factor_lj*forcelj) * r2inv;

        if (rsq > cut_out_on_sq) {
          rsw = (sqrt(rsq) - cut_out_on)/cut_out_diff; 
	  fforce  *= 1.0 + rsw*rsw*(2.0*rsw-3.0);
        }

        f[i][0] += delx*fforce;
        f[i][1] += dely*fforce;
        f[i][2] += delz*fforce;
        if (newton_pair || j < nlocal) {
          f[j][0] -= delx*fforce;
          f[j][1] -= dely*fforce;
          f[j][2] -= delz*fforce;
        }
      }
    }
  }
}

/* ---------------------------------------------------------------------- */

void PairLJCharmmCoulLong::compute_middle()
{
  int i,j,k,numneigh,itype,jtype;
  double qtmp,xtmp,ytmp,ztmp,delx,dely,delz;
  double rsq,r2inv,r6inv,forcecoul,forcelj,fforce,factor_coul,factor_lj;
  double philj,switch1,switch2;
  double rsw;
  int *neighs;

  double **f = atom->f;
  double **x = atom->x;
  double *q = atom->q;
  int *type = atom->type;
  int nlocal = atom->nlocal;
  int nall = atom->nlocal + atom->nghost;
  double *special_coul = force->special_coul;
  double *special_lj = force->special_lj;
  int newton_pair = force->newton_pair;
  double qqrd2e = force->qqrd2e;

  double cut_in_off = cut_respa[0];
  double cut_in_on = cut_respa[1];
  double cut_out_on = cut_respa[2];
  double cut_out_off = cut_respa[3];

  double cut_in_diff = cut_in_on - cut_in_off;
  double cut_out_diff = cut_out_off - cut_out_on;
  double cut_in_off_sq = cut_in_off*cut_in_off;
  double cut_in_on_sq = cut_in_on*cut_in_on;
  double cut_out_on_sq = cut_out_on*cut_out_on;
  double cut_out_off_sq = cut_out_off*cut_out_off;

  // loop over neighbors of my atoms

  for (i = 0; i < nlocal; i++) {

    qtmp = q[i];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    itype = type[i];
    neighs = neighbor->firstneigh_middle[i];
    numneigh = neighbor->numneigh_middle[i];

    for (k = 0; k < numneigh; k++) {
      j = neighs[k];

      if (j < nall) factor_coul = factor_lj = 1.0;
      else {
	factor_coul = special_coul[j/nall];
	factor_lj = special_lj[j/nall];
	j %= nall;
      }
      
      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx*delx + dely*dely + delz*delz;

      if (rsq < cut_out_off_sq && rsq > cut_in_off_sq) {
	r2inv = 1.0/rsq;
        forcecoul = qqrd2e * qtmp*q[j]*sqrt(r2inv);
        if (factor_coul < 1.0) forcecoul -= (1.0-factor_coul)*forcecoul;

	r6inv = r2inv*r2inv*r2inv;
	jtype = type[j];
	forcelj = r6inv * (lj1[itype][jtype]*r6inv - lj2[itype][jtype]);
	if (rsq > cut_lj_innersq) {
	  switch1 = (cut_ljsq-rsq) * (cut_ljsq-rsq) *
	    (cut_ljsq + 2.0*rsq - 3.0*cut_lj_innersq) / denom_lj;
	  switch2 = 12.0*rsq * (cut_ljsq-rsq) * 
	    (rsq-cut_lj_innersq) / denom_lj;
	  philj = r6inv * (lj3[itype][jtype]*r6inv - lj4[itype][jtype]);
	  forcelj = forcelj*switch1 + philj*switch2;
	}

	fforce = (forcecoul + factor_lj*forcelj) * r2inv;
        if (rsq < cut_in_on_sq) {
	  rsw = (sqrt(rsq) - cut_in_off)/cut_in_diff; 
	  fforce *= rsw*rsw*(3.0 - 2.0*rsw);
	}
        if (rsq > cut_out_on_sq) {
	  rsw = (sqrt(rsq) - cut_out_on)/cut_out_diff; 
	  fforce *= 1.0 + rsw*rsw*(2.0*rsw - 3.0);
	}

	f[i][0] += delx*fforce;
	f[i][1] += dely*fforce;
	f[i][2] += delz*fforce;
	if (newton_pair || j < nlocal) {
	  f[j][0] -= delx*fforce;
	  f[j][1] -= dely*fforce;
	  f[j][2] -= delz*fforce;
	}
      }
    }
  }
}

/* ---------------------------------------------------------------------- */

void PairLJCharmmCoulLong::compute_outer(int eflag, int vflag)
{
  int i,j,k,numneigh,itype,jtype,itable;
  double qtmp,xtmp,ytmp,ztmp,delx,dely,delz,fraction,table;
  double r,r2inv,r6inv,forcecoul,forcelj,fforce,factor_coul,factor_lj;
  double grij,expm2,prefactor,t,erfc;
  double factor,phicoul,philj,switch1,switch2;
  double rsw;
  int *neighs;
  float rsq;
  int *int_rsq = (int *) &rsq;

  eng_vdwl = eng_coul = 0.0;
  if (vflag) for (i = 0; i < 6; i++) virial[i] = 0.0;

  double **f = atom->f;
  double **x = atom->x;
  double *q = atom->q;
  int *type = atom->type;
  int nlocal = atom->nlocal;
  int nall = atom->nlocal + atom->nghost;
  double *special_coul = force->special_coul;
  double *special_lj = force->special_lj;
  int newton_pair = force->newton_pair;
  double qqrd2e = force->qqrd2e;
  
  double cut_in_off = cut_respa[2];
  double cut_in_on = cut_respa[3];

  double cut_in_diff = cut_in_on - cut_in_off;
  double cut_in_off_sq = cut_in_off*cut_in_off;
  double cut_in_on_sq = cut_in_on*cut_in_on;

  // loop over neighbors of my atoms

  for (i = 0; i < nlocal; i++) {
    
    qtmp = q[i];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    itype = type[i];
    neighs = neighbor->firstneigh[i];
    numneigh = neighbor->numneigh[i];
    
    for (k = 0; k < numneigh; k++) {
      j = neighs[k];
      
      if (j < nall) factor_coul = factor_lj = 1.0;
      else {
	factor_coul = special_coul[j/nall];
	factor_lj = special_lj[j/nall];
	j %= nall;
      }
      
      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx*delx + dely*dely + delz*delz;
      jtype = type[j];
      
      if (rsq < cut_bothsq) {
	r2inv = 1.0/rsq;
	
	if (rsq < cut_coulsq) {
	  if (!ncoultablebits || rsq <= tabinnersq) {
	    r = sqrtf(rsq);
	    grij = g_ewald * r;
	    expm2 = exp(-grij*grij);
	    t = 1.0 / (1.0 + EWALD_P*grij);
	    erfc = t * (A1+t*(A2+t*(A3+t*(A4+t*A5)))) * expm2;
	    prefactor = qqrd2e * qtmp*q[j]/r;
	    forcecoul = prefactor * (erfc + EWALD_F*grij*expm2 - 1.0);
	    if (rsq > cut_in_off_sq) {
	      if (rsq < cut_in_on_sq) {
		rsw = (r - cut_in_off)/cut_in_diff; 
		forcecoul += prefactor*rsw*rsw*(3.0 - 2.0*rsw);
		if (factor_coul < 1.0)
		  forcecoul -= 
		    (1.0-factor_coul)*prefactor*rsw*rsw*(3.0 - 2.0*rsw);
	      } else {
		forcecoul += prefactor;
		if (factor_coul < 1.0)
		  forcecoul -= (1.0-factor_coul)*prefactor;
	      }
	    }
	  } else {
	    itable = *int_rsq & ncoulmask;
	    itable >>= ncoulshiftbits;
	    fraction = (rsq - rtable[itable]) * drtable[itable];
	    table = ftable[itable] + fraction*dftable[itable];
	    forcecoul = qtmp*q[j] * table;
	    if (factor_coul < 1.0) {
	      table = ctable[itable] + fraction*dctable[itable];
	      prefactor = qtmp*q[j] * table;
	      forcecoul -= (1.0-factor_coul)*prefactor;
	    }
	  }
	} else forcecoul = 0.0;
	
	if (rsq < cut_ljsq && rsq > cut_in_off_sq) {
	  r6inv = r2inv*r2inv*r2inv;
	  forcelj = r6inv * (lj1[itype][jtype]*r6inv - lj2[itype][jtype]);
	  if (rsq > cut_lj_innersq) {
	    switch1 = (cut_ljsq-rsq) * (cut_ljsq-rsq) *
	      (cut_ljsq + 2.0*rsq - 3.0*cut_lj_innersq) / denom_lj;
	    switch2 = 12.0*rsq * (cut_ljsq-rsq) * 
	      (rsq-cut_lj_innersq) / denom_lj;
	    philj = r6inv * (lj3[itype][jtype]*r6inv - lj4[itype][jtype]);
	    forcelj = forcelj*switch1 + philj*switch2;
	  }
	  if (rsq < cut_in_on_sq) {
	    rsw = (sqrtf(rsq) - cut_in_off)/cut_in_diff; 
	    forcelj *= rsw*rsw*(3.0 - 2.0*rsw);
	  }
	} else forcelj = 0.0;
	
	fforce = (forcecoul + forcelj) * r2inv;
	
	f[i][0] += delx*fforce;
	f[i][1] += dely*fforce;
	f[i][2] += delz*fforce;
	if (newton_pair || j < nlocal) {
	  f[j][0] -= delx*fforce;
	  f[j][1] -= dely*fforce;
	  f[j][2] -= delz*fforce;
	}
	
	if (eflag) {
	  if (newton_pair || j < nlocal) factor = 1.0;
	  else factor = 0.5;
	  if (rsq < cut_coulsq) {
	    if (!ncoultablebits || rsq <= tabinnersq) {
	      phicoul = prefactor*erfc;
	      if (factor_coul < 1.0) phicoul -= (1.0-factor_coul)*prefactor;
	    } else {
	      table = etable[itable] + fraction*detable[itable];
	      phicoul = qtmp*q[j] * table;
	      if (factor_coul < 1.0) {
		table = ptable[itable] + fraction*dptable[itable];
		prefactor = qtmp*q[j] * table;
		phicoul -= (1.0-factor_coul)*prefactor;
	      }
	    }
	    eng_coul += factor*phicoul;
	  }
	  if (rsq < cut_ljsq) {
	    r6inv = r2inv*r2inv*r2inv;
	    philj = r6inv*(lj3[itype][jtype]*r6inv-lj4[itype][jtype]);
	    if (rsq > cut_lj_innersq) {
	      switch1 = (cut_ljsq-rsq) * (cut_ljsq-rsq) *
		(cut_ljsq + 2.0*rsq - 3.0*cut_lj_innersq) / denom_lj;
	      philj *= switch1;
	    }
	    eng_vdwl += factor*factor_lj*philj;
	  }
	}
	
	if (vflag) {
	  if (rsq < cut_coulsq) {
	    if (!ncoultablebits || rsq <= tabinnersq) {
	      forcecoul = prefactor * (erfc + EWALD_F*grij*expm2);
	      if (factor_coul < 1.0) forcecoul -= (1.0-factor_coul)*prefactor;
	    } else {
	      table = vtable[itable] + fraction*dvtable[itable];
	      forcecoul = qtmp*q[j] * table;
	      if (factor_coul < 1.0) {
		table = ptable[itable] + fraction*dptable[itable];
		prefactor = qtmp*q[j] * table;
		forcecoul -= (1.0-factor_coul)*prefactor;
	      }
	    }
	  } else forcecoul = 0.0;

	  if (rsq <= cut_in_off_sq) {
	    r6inv = r2inv*r2inv*r2inv;
	    forcelj = r6inv * (lj1[itype][jtype]*r6inv - lj2[itype][jtype]);
	    if (rsq > cut_lj_innersq) {
	      switch1 = (cut_ljsq-rsq) * (cut_ljsq-rsq) *
		(cut_ljsq + 2.0*rsq - 3.0*cut_lj_innersq) / denom_lj;
	      switch2 = 12.0*rsq * (cut_ljsq-rsq) * 
		(rsq-cut_lj_innersq) / denom_lj;
	      philj = r6inv * (lj3[itype][jtype]*r6inv - lj4[itype][jtype]);
	      forcelj = forcelj*switch1 + philj*switch2;
	    }
	  } else if (rsq <= cut_in_on_sq) {
	    forcelj = r6inv * (lj1[itype][jtype]*r6inv - lj2[itype][jtype]);
	    if (rsq > cut_lj_innersq) {
	      switch1 = (cut_ljsq-rsq) * (cut_ljsq-rsq) *
		(cut_ljsq + 2.0*rsq - 3.0*cut_lj_innersq) / denom_lj;
	      switch2 = 12.0*rsq * (cut_ljsq-rsq) * 
		(rsq-cut_lj_innersq) / denom_lj;
	      philj = r6inv * (lj3[itype][jtype]*r6inv - lj4[itype][jtype]);
	      forcelj = forcelj*switch1 + philj*switch2;
	    }
	  }
	  
	  fforce = (forcecoul + factor_lj*forcelj) * r2inv;

	  if (newton_pair == 0 && j >= nlocal) fforce *= 0.5;
	  virial[0] += delx*delx*fforce;
	  virial[1] += dely*dely*fforce;
	  virial[2] += delz*delz*fforce;
	  virial[3] += delx*dely*fforce;
	  virial[4] += delx*delz*fforce;
	  virial[5] += dely*delz*fforce;
	}
      }
    }
  }
}

/* ----------------------------------------------------------------------
   allocate all arrays
------------------------------------------------------------------------- */

void PairLJCharmmCoulLong::allocate()
{
  allocated = 1;
  int n = atom->ntypes;

  setflag = memory->create_2d_int_array(n+1,n+1,"pair:setflag");
  for (int i = 1; i <= n; i++)
    for (int j = i; j <= n; j++)
      setflag[i][j] = 0;

  cutsq = memory->create_2d_double_array(n+1,n+1,"pair:cutsq");

  epsilon = memory->create_2d_double_array(n+1,n+1,"pair:epsilon");
  sigma = memory->create_2d_double_array(n+1,n+1,"pair:sigma");
  eps14 = memory->create_2d_double_array(n+1,n+1,"pair:eps14");
  sigma14 = memory->create_2d_double_array(n+1,n+1,"pair:sigma14");
  lj1 = memory->create_2d_double_array(n+1,n+1,"pair:lj1");
  lj2 = memory->create_2d_double_array(n+1,n+1,"pair:lj2");
  lj3 = memory->create_2d_double_array(n+1,n+1,"pair:lj3");
  lj4 = memory->create_2d_double_array(n+1,n+1,"pair:lj4");
  lj14_1 = memory->create_2d_double_array(n+1,n+1,"pair:lj14_1");
  lj14_2 = memory->create_2d_double_array(n+1,n+1,"pair:lj14_2");
  lj14_3 = memory->create_2d_double_array(n+1,n+1,"pair:lj14_3");
  lj14_4 = memory->create_2d_double_array(n+1,n+1,"pair:lj14_4");
}

/* ----------------------------------------------------------------------
   global settings
   unlike other pair styles,
     there are no individual pair settings that these override
------------------------------------------------------------------------- */

void PairLJCharmmCoulLong::settings(int narg, char **arg)
{
  if (narg != 2 && narg != 3) error->all("Illegal pair_style command");

  cut_lj_inner = atof(arg[0]);
  cut_lj = atof(arg[1]);
  if (narg == 2) cut_coul = cut_lj;
  else cut_coul = atof(arg[2]);
}

/* ----------------------------------------------------------------------
   set coeffs for one or more type pairs
------------------------------------------------------------------------- */

void PairLJCharmmCoulLong::coeff(int narg, char **arg)
{
  if (narg != 4 && narg != 6) error->all("Illegal pair_coeff command");
  if (!allocated) allocate();

  int ilo,ihi,jlo,jhi;
  force->bounds(arg[0],atom->ntypes,ilo,ihi);
  force->bounds(arg[1],atom->ntypes,jlo,jhi);

  double epsilon_one = atof(arg[2]);
  double sigma_one = atof(arg[3]);
  double eps14_one = epsilon_one;
  double sigma14_one = sigma_one;
  if (narg == 6) {
    eps14_one = atof(arg[4]);
    sigma14_one = atof(arg[5]);
  }

  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
    for (int j = MAX(jlo,i); j <= jhi; j++) {
      epsilon[i][j] = epsilon_one;
      sigma[i][j] = sigma_one;
      eps14[i][j] = eps14_one;
      sigma14[i][j] = sigma14_one;
      setflag[i][j] = 1;
      count++;
    }
  }

  if (count == 0) error->all("Incorrect args for pair coefficients");
}

/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

double PairLJCharmmCoulLong::init_one(int i, int j)
{
  // always mix arithmetically

  if (setflag[i][j] == 0) {
    epsilon[i][j] = sqrt(epsilon[i][i]*epsilon[j][j]);
    sigma[i][j] = 0.5 * (sigma[i][i] + sigma[j][j]);
    eps14[i][j] = sqrt(eps14[i][i]*eps14[j][j]);
    sigma14[i][j] = 0.5 * (sigma14[i][i] + sigma14[j][j]);
  }

  double cut = MAX(cut_lj,cut_coul);

  lj1[i][j] = 48.0 * epsilon[i][j] * pow(sigma[i][j],12.0);
  lj2[i][j] = 24.0 * epsilon[i][j] * pow(sigma[i][j],6.0);
  lj3[i][j] = 4.0 * epsilon[i][j] * pow(sigma[i][j],12.0);
  lj4[i][j] = 4.0 * epsilon[i][j] * pow(sigma[i][j],6.0);
  lj14_1[i][j] = 48.0 * eps14[i][j] * pow(sigma14[i][j],12.0);
  lj14_2[i][j] = 24.0 * eps14[i][j] * pow(sigma14[i][j],6.0);
  lj14_3[i][j] = 4.0 * eps14[i][j] * pow(sigma14[i][j],12.0);
  lj14_4[i][j] = 4.0 * eps14[i][j] * pow(sigma14[i][j],6.0);
     
  lj1[j][i] = lj1[i][j];
  lj2[j][i] = lj2[i][j];
  lj3[j][i] = lj3[i][j];
  lj4[j][i] = lj4[i][j];
  lj14_1[j][i] = lj14_1[i][j];
  lj14_2[j][i] = lj14_2[i][j];
  lj14_3[j][i] = lj14_3[i][j];
  lj14_4[j][i] = lj14_4[i][j];

  return cut;
}

/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */

void PairLJCharmmCoulLong::init_style()
{
  if (!atom->q_flag)
    error->all("Pair style lj/charmm/coul/long requires atom attribute q");

  // require cut_lj_inner < cut_lj

  if (cut_lj_inner >= cut_lj) 
    error->all("Pair inner cutoff >= Pair outer cutoff");

  cut_lj_innersq = cut_lj_inner * cut_lj_inner;
  cut_ljsq = cut_lj * cut_lj;
  cut_coulsq = cut_coul * cut_coul;
  cut_bothsq = MAX(cut_ljsq,cut_coulsq);

  denom_lj = (cut_ljsq-cut_lj_innersq) * (cut_ljsq-cut_lj_innersq) * 
    (cut_ljsq-cut_lj_innersq);

  // set & error check interior rRESPA cutoffs

  cut_respa = NULL;
  if (strcmp(update->integrate_style,"respa") == 0) {
    if (((Respa *) update->integrate)->level_inner >= 0) {
      cut_respa = ((Respa *) update->integrate)->cutoff;
      if (MIN(cut_lj,cut_coul) < cut_respa[3])
	error->all("Pair cutoff < Respa interior cutoff");
      if (cut_lj_inner < cut_respa[1])
	error->all("Pair inner cutoff < Respa interior cutoff");
    }
  } else cut_respa = NULL;

  // insure use of KSpace long-range solver, set g_ewald

  if (force->kspace == NULL) 
    error->all("Pair style is incompatible with KSpace style");
  else if (strcmp(force->kspace_style,"ewald") == 0)
    g_ewald = force->kspace->g_ewald;
  else if (strcmp(force->kspace_style,"pppm") == 0)
    g_ewald = force->kspace->g_ewald;
  else if (strcmp(force->kspace_style,"pppm2") == 0)
    g_ewald = force->kspace->g_ewald;
  else error->all("Pair style is incompatible with KSpace style");

  // setup force tables

  if (ncoultablebits) init_tables();
}

/* ----------------------------------------------------------------------
   setup force tables used in compute routines
------------------------------------------------------------------------- */

void PairLJCharmmCoulLong::init_tables()
{
  int masklo,maskhi;
  double r,grij,expm2,derfc,rsw;
  double qqrd2e = force->qqrd2e;

  tabinnersq = tabinner*tabinner;
  init_bitmap(tabinner,cut_coul,ncoultablebits,
	      masklo,maskhi,ncoulmask,ncoulshiftbits);
  
  int ntable = 1;
  for (int i = 0; i < ncoultablebits; i++) ntable *= 2;
  
  // linear lookup tables of length N = 2^ncoultablebits
  // stored value = value at lower edge of bin
  // d values = delta from lower edge to upper edge of bin

  if (ftable) free_tables();
  
  rtable = (double *) memory->smalloc(ntable*sizeof(double),"pair:rtable");
  ftable = (double *) memory->smalloc(ntable*sizeof(double),"pair:ftable");
  ctable = (double *) memory->smalloc(ntable*sizeof(double),"pair:ctable");
  etable = (double *) memory->smalloc(ntable*sizeof(double),"pair:etable");
  drtable = (double *) memory->smalloc(ntable*sizeof(double),"pair:drtable");
  dftable = (double *) memory->smalloc(ntable*sizeof(double),"pair:dftable");
  dctable = (double *) memory->smalloc(ntable*sizeof(double),"pair:dctable");
  detable = (double *) memory->smalloc(ntable*sizeof(double),"pair:detable");

  if (cut_respa == NULL) {
    vtable = ptable = dvtable = dptable = NULL;
  } else {
    vtable = (double *) memory->smalloc(ntable*sizeof(double),"pair:vtable");
    ptable = (double *) memory->smalloc(ntable*sizeof(double),"pair:ptable");
    dvtable = (double *) memory->smalloc(ntable*sizeof(double),"pair:dvtable");
    dptable = (double *) memory->smalloc(ntable*sizeof(double),"pair:dptable");
  }

  float rsq;
  int *int_rsq = (int *) &rsq;  
  float minrsq;
  int *int_minrsq = (int *) &minrsq;
  int itablemin;
  *int_minrsq = 0 << ncoulshiftbits;
  *int_minrsq = *int_minrsq | maskhi;
    
  for (int i = 0; i < ntable; i++) {
    *int_rsq = i << ncoulshiftbits;
    *int_rsq = *int_rsq | masklo;
    if (rsq < tabinnersq) {
      *int_rsq = i << ncoulshiftbits;
      *int_rsq = *int_rsq | maskhi;
    }
    r = sqrtf(rsq);
    grij = g_ewald * r;
    expm2 = exp(-grij*grij);
    derfc = erfc(grij);
    if (cut_respa == NULL) {
      rtable[i] = rsq;
      ftable[i] = qqrd2e/r * (derfc + EWALD_F*grij*expm2);
      ctable[i] = qqrd2e/r;
      etable[i] = qqrd2e/r * derfc;
    } else {
      rtable[i] = rsq;
      ftable[i] = qqrd2e/r * (derfc + EWALD_F*grij*expm2 - 1.0);
      ctable[i] = 0.0;
      etable[i] = qqrd2e/r * derfc;
      ptable[i] = qqrd2e/r;
      vtable[i] = qqrd2e/r * (derfc + EWALD_F*grij*expm2);
      if (rsq > cut_respa[2]*cut_respa[2]) {
	if (rsq < cut_respa[3]*cut_respa[3]) {
	  rsw = (r - cut_respa[2])/(cut_respa[3] - cut_respa[2]); 
	  ftable[i] += qqrd2e/r * rsw*rsw*(3.0 - 2.0*rsw);
	  ctable[i] = qqrd2e/r * rsw*rsw*(3.0 - 2.0*rsw);
	} else {
	  ftable[i] = qqrd2e/r * (derfc + EWALD_F*grij*expm2);
	  ctable[i] = qqrd2e/r;
	}
      }
    }
    minrsq = MIN(minrsq,rsq);
  }

  tabinnersq = minrsq;
  
  int ntablem1 = ntable - 1;
  
  for (int i = 0; i < ntablem1; i++) {
    drtable[i] = 1.0/(rtable[i+1] - rtable[i]);
    dftable[i] = ftable[i+1] - ftable[i];
    dctable[i] = ctable[i+1] - ctable[i];
    detable[i] = etable[i+1] - etable[i];
  }

  if (cut_respa) {
    for (int i = 0; i < ntablem1; i++) {
      dvtable[i] = vtable[i+1] - vtable[i];
      dptable[i] = ptable[i+1] - ptable[i];
    }
  }
  
  // get the delta values for the last table entries 
  // tables are connected periodically between 0 and ntablem1
    
  drtable[ntablem1] = 1.0/(rtable[0] - rtable[ntablem1]);
  dftable[ntablem1] = ftable[0] - ftable[ntablem1];
  dctable[ntablem1] = ctable[0] - ctable[ntablem1];
  detable[ntablem1] = etable[0] - etable[ntablem1];
  if (cut_respa) {
    dvtable[ntablem1] = vtable[0] - vtable[ntablem1];
    dptable[ntablem1] = ptable[0] - ptable[ntablem1];
  }

  // get the correct delta values at itablemax    
  // smallest r is in bin itablemin
  // largest r is in bin itablemax, which is itablemin-1,
  //   or ntablem1 if itablemin=0
  // deltas at itablemax only needed if corresponding rsq < cut*cut
  // if so, compute deltas between rsq and cut*cut 
	
  double f_tmp,c_tmp,e_tmp,p_tmp,v_tmp;
  itablemin = *int_minrsq & ncoulmask;
  itablemin >>= ncoulshiftbits;  
  int itablemax = itablemin - 1; 
  if (itablemin == 0) itablemax = ntablem1;     
  *int_rsq = itablemax << ncoulshiftbits;
  *int_rsq = *int_rsq | maskhi;          

  if (rsq < cut_coulsq) {
    rsq = cut_coulsq;  
    r = sqrtf(rsq);
    grij = g_ewald * r;
    expm2 = exp(-grij*grij);
    derfc = erfc(grij);

    if (cut_respa == NULL) {
      f_tmp = qqrd2e/r * (derfc + EWALD_F*grij*expm2);
      c_tmp = qqrd2e/r;
      e_tmp = qqrd2e/r * derfc;
    } else {
      f_tmp = qqrd2e/r * (derfc + EWALD_F*grij*expm2 - 1.0);
      c_tmp = 0.0;
      e_tmp = qqrd2e/r * derfc;
      p_tmp = qqrd2e/r;
      v_tmp = qqrd2e/r * (derfc + EWALD_F*grij*expm2);
      if (rsq > cut_respa[2]*cut_respa[2]) {
        if (rsq < cut_respa[3]*cut_respa[3]) {
          rsw = (r - cut_respa[2])/(cut_respa[3] - cut_respa[2]); 
          f_tmp += qqrd2e/r * rsw*rsw*(3.0 - 2.0*rsw);
          c_tmp = qqrd2e/r * rsw*rsw*(3.0 - 2.0*rsw);
        } else {
          f_tmp = qqrd2e/r * (derfc + EWALD_F*grij*expm2);
          c_tmp = qqrd2e/r;
        }
      }
    }

    drtable[itablemax] = 1.0/(rsq - rtable[itablemax]);   
    dftable[itablemax] = f_tmp - ftable[itablemax];
    dctable[itablemax] = c_tmp - ctable[itablemax];
    detable[itablemax] = e_tmp - etable[itablemax];
    if (cut_respa) {
      dvtable[itablemax] = v_tmp - vtable[itablemax];
      dptable[itablemax] = p_tmp - ptable[itablemax];
    }   
  }
}

/* ----------------------------------------------------------------------
  proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairLJCharmmCoulLong::write_restart(FILE *fp)
{
  write_restart_settings(fp);

  int i,j;
  for (i = 1; i <= atom->ntypes; i++)
    for (j = i; j <= atom->ntypes; j++) {
      fwrite(&setflag[i][j],sizeof(int),1,fp);
      if (setflag[i][j]) {
	fwrite(&epsilon[i][j],sizeof(double),1,fp);
	fwrite(&sigma[i][j],sizeof(double),1,fp);
	fwrite(&eps14[i][j],sizeof(double),1,fp);
	fwrite(&sigma14[i][j],sizeof(double),1,fp);
      }
    }
}

/* ----------------------------------------------------------------------
  proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairLJCharmmCoulLong::read_restart(FILE *fp)
{
  read_restart_settings(fp);

  allocate();

  int i,j;
  int me = comm->me;
  for (i = 1; i <= atom->ntypes; i++)
    for (j = i; j <= atom->ntypes; j++) {
      if (me == 0) fread(&setflag[i][j],sizeof(int),1,fp);
      MPI_Bcast(&setflag[i][j],1,MPI_INT,0,world);
      if (setflag[i][j]) {
	if (me == 0) {
	  fread(&epsilon[i][j],sizeof(double),1,fp);
	  fread(&sigma[i][j],sizeof(double),1,fp);
	  fread(&eps14[i][j],sizeof(double),1,fp);
	  fread(&sigma14[i][j],sizeof(double),1,fp);
	}
	MPI_Bcast(&epsilon[i][j],1,MPI_DOUBLE,0,world);
	MPI_Bcast(&sigma[i][j],1,MPI_DOUBLE,0,world);
	MPI_Bcast(&eps14[i][j],1,MPI_DOUBLE,0,world);
	MPI_Bcast(&sigma14[i][j],1,MPI_DOUBLE,0,world);
      }
    }
}

/* ----------------------------------------------------------------------
  proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairLJCharmmCoulLong::write_restart_settings(FILE *fp)
{
  fwrite(&cut_lj_inner,sizeof(double),1,fp);
  fwrite(&cut_lj,sizeof(double),1,fp);
  fwrite(&cut_coul,sizeof(double),1,fp);
  fwrite(&offset_flag,sizeof(int),1,fp);
  fwrite(&mix_flag,sizeof(int),1,fp);
}

/* ----------------------------------------------------------------------
  proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairLJCharmmCoulLong::read_restart_settings(FILE *fp)
{
  if (comm->me == 0) {
    fread(&cut_lj_inner,sizeof(double),1,fp);
    fread(&cut_lj,sizeof(double),1,fp);
    fread(&cut_coul,sizeof(double),1,fp);
    fread(&offset_flag,sizeof(int),1,fp);
    fread(&mix_flag,sizeof(int),1,fp);
  }
  MPI_Bcast(&cut_lj_inner,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&cut_lj,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&cut_coul,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&offset_flag,1,MPI_INT,0,world);
  MPI_Bcast(&mix_flag,1,MPI_INT,0,world);
}

/* ----------------------------------------------------------------------
   free memory for tables used in pair computations
------------------------------------------------------------------------- */

void PairLJCharmmCoulLong::free_tables()
{
  memory->sfree(rtable);
  memory->sfree(drtable);
  memory->sfree(ftable);
  memory->sfree(dftable);
  memory->sfree(ctable);
  memory->sfree(dctable);
  memory->sfree(etable);
  memory->sfree(detable);
  memory->sfree(vtable);
  memory->sfree(dvtable);
  memory->sfree(ptable);
  memory->sfree(dptable);
}

/* ---------------------------------------------------------------------- */

void PairLJCharmmCoulLong::single(int i, int j, int itype, int jtype,
				  double rsq,
				  double factor_coul, double factor_lj,
				  int eflag, One &one)
{
  double r2inv,r6inv,r,grij,expm2,t,erfc,prefactor;
  double switch1,switch2,fraction,table,forcecoul,forcelj,phicoul,philj;
  int itable;

  r2inv = 1.0/rsq;
  if (rsq < cut_coulsq) {
    if (!ncoultablebits || rsq <= tabinnersq) {
      r = sqrt(rsq);
      grij = g_ewald * r;
      expm2 = exp(-grij*grij);
      t = 1.0 / (1.0 + EWALD_P*grij);
      erfc = t * (A1+t*(A2+t*(A3+t*(A4+t*A5)))) * expm2;
      prefactor = force->qqrd2e * atom->q[i]*atom->q[j]/r;
      forcecoul = prefactor * (erfc + EWALD_F*grij*expm2);
      if (factor_coul < 1.0) forcecoul -= (1.0-factor_coul)*prefactor;
    } else {
      float rsq_single = rsq;
      int *int_rsq = (int *) &rsq_single;
      itable = *int_rsq & ncoulmask;
      itable >>= ncoulshiftbits;
      fraction = (rsq_single - rtable[itable]) * drtable[itable];
      table = ftable[itable] + fraction*dftable[itable];
      forcecoul = atom->q[i]*atom->q[j] * table;
      if (factor_coul < 1.0) {
	table = ctable[itable] + fraction*dctable[itable];
	prefactor = atom->q[i]*atom->q[j] * table;
	forcecoul -= (1.0-factor_coul)*prefactor;
      }
    }
  } else forcecoul = 0.0;
  if (rsq < cut_ljsq) {
    r6inv = r2inv*r2inv*r2inv;
    forcelj = r6inv * (lj1[itype][jtype]*r6inv - lj2[itype][jtype]);
    if (rsq > cut_lj_innersq) {
      switch1 = (cut_ljsq-rsq) * (cut_ljsq-rsq) *
	(cut_ljsq + 2.0*rsq - 3.0*cut_lj_innersq) / denom_lj;
      switch2 = 12.0*rsq * (cut_ljsq-rsq) * 
	(rsq-cut_lj_innersq) / denom_lj;
      philj = r6inv * (lj3[itype][jtype]*r6inv - lj4[itype][jtype]);
      forcelj = forcelj*switch1 + philj*switch2;
    }
  } else forcelj = 0.0;
  one.fforce = (forcecoul + factor_lj*forcelj) * r2inv;
  
  if (eflag) {
    if (rsq < cut_coulsq) {
      if (!ncoultablebits || rsq <= tabinnersq)
	phicoul = prefactor*erfc;
      else {
	table = etable[itable] + fraction*detable[itable];
	phicoul = atom->q[i]*atom->q[j] * table;
      }
      if (factor_coul < 1.0) phicoul -= (1.0-factor_coul)*prefactor;
      one.eng_coul = phicoul;
    } else one.eng_coul = 0.0;
    if (rsq < cut_ljsq) {
      philj = r6inv*(lj3[itype][jtype]*r6inv-lj4[itype][jtype]);
      if (rsq > cut_lj_innersq) {
	switch1 = (cut_ljsq-rsq) * (cut_ljsq-rsq) *
	  (cut_ljsq + 2.0*rsq - 3.0*cut_lj_innersq) / denom_lj;
	philj *= switch1;
      }
      one.eng_vdwl = factor_lj*philj;
    } else one.eng_vdwl = 0.0;
  }
}

/* ---------------------------------------------------------------------- */

void PairLJCharmmCoulLong::extract_charmm(double ***p_lj14_1, 
					  double ***p_lj14_2,
					  double ***p_lj14_3,
					  double ***p_lj14_4,
					  int *p_implicit_flag)
{
  *p_lj14_1 = lj14_1;
  *p_lj14_2 = lj14_2;
  *p_lj14_3 = lj14_3;
  *p_lj14_4 = lj14_4;
  *p_implicit_flag = 0;
}

/* ---------------------------------------------------------------------- */

void PairLJCharmmCoulLong::extract_long(double *p_cut_coul)
{
  *p_cut_coul = cut_coul;
}
