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
   Contributing authors:  
     James Fischer, High Performance Technologies, Inc.
     David Richie, Stone Ridge Technology
     Vincent Natol, Stone Ridge Technology
------------------------------------------------------------------------- */

#ifndef PAIR_LJ_CHARMM_COUL_LONG_OPT_H
#define PAIR_LJ_CHARMM_COUL_LONG_OPT_H

#include "pair_lj_charmm_coul_long.h"
#include "math.h"
#include "stdio.h"
#include "stdlib.h"
#include "string.h"
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

#define MIN(a,b) ((a) < (b) ? (a) : (b))
#define MAX(a,b) ((a) > (b) ? (a) : (b))

#define EWALD_F   1.12837917
#define EWALD_P   0.3275911
#define A1        0.254829592
#define A2       -0.284496736
#define A3        1.421413741
#define A4       -1.453152027
#define A5        1.061405429

namespace LAMMPS_NS {

class PairLJCharmmCoulLongOpt : public PairLJCharmmCoulLong {
 public:
  PairLJCharmmCoulLongOpt(class LAMMPS *);
  void compute(int, int);

 private:
  template < int EFLAG, int VFLAG, int NEWTON_PAIR > void eval();
};

template < int EFLAG, int VFLAG, int NEWTON_PAIR >
void PairLJCharmmCoulLongOpt::eval()
{
  typedef struct { double x,y,z; } vec3_t;
  
  typedef struct {
    double cutsq,lj1,lj2,lj3,lj4,offset;
    double _pad[2];
  } fast_alpha_t;
  
  int i,j,k,itype,jtype,itable;
  double fraction,table;
  double r,r2inv,r6inv,forcecoul,forcelj,fforce,factor_coul,factor_lj;
  double grij,expm2,prefactor,t,erfc;
  double factor,phicoul,philj,switch1,switch2;
  
  double** __restrict__ f;
  
  float rsq;
  int *int_rsq = (int *) &rsq;
  
  eng_vdwl = eng_coul = 0.0;
  if (VFLAG) for (i = 0; i < 6; i++) virial[i] = 0.0;
  
  if (VFLAG == 2) f = update->f_pair;
  else f = atom->f;
  
  double** __restrict__ x = atom->x;
  double* __restrict__ q = atom->q;
  int* __restrict__ type = atom->type;
  int nlocal = atom->nlocal;
  int nall = atom->nlocal + atom->nghost;
  double* __restrict__ special_coul = force->special_coul;
  double* __restrict__ special_lj = force->special_lj;
  double qqrd2e = force->qqrd2e;
  int** __restrict__ firstneigh = neighbor->firstneigh;
  int* __restrict__ num = neighbor->numneigh;
  
  vec3_t* __restrict__ xx = (vec3_t*)x[0];
  vec3_t* __restrict__ ff = (vec3_t*)f[0];
  
  int ntypes = atom->ntypes;
  int ntypes2 = ntypes*ntypes;
  
  double tmp_coef1 = 1.0/denom_lj;
  double tmp_coef2 = cut_ljsq - 3.0*cut_lj_innersq;
  
  fast_alpha_t* __restrict__ fast_alpha = 
    (fast_alpha_t*)malloc(ntypes2*sizeof(fast_alpha_t));
  for( int i = 0; i < ntypes; i++) for( int j = 0; j < ntypes; j++) {
    fast_alpha_t& a = fast_alpha[i*ntypes+j];
    a.cutsq = cutsq[i+1][j+1];
    a.lj1 = lj1[i+1][j+1];
    a.lj2 = lj2[i+1][j+1];
    a.lj3 = lj3[i+1][j+1];
    a.lj4 = lj4[i+1][j+1];
  }
  fast_alpha_t* __restrict__ tabsix = fast_alpha;
  
  // loop over neighbors of my atoms
  
  for (i = 0; i < nlocal; i++) {
    double qtmp = q[i];
    double xtmp = xx[i].x;
    double ytmp = xx[i].y;
    double ztmp = xx[i].z;
    itype = type[i] - 1;
    int* __restrict__ neighs = firstneigh[i];
    int numneigh = num[i];
    
    double tmpfx = 0.0;
    double tmpfy = 0.0;
    double tmpfz = 0.0;
    
    fast_alpha_t* __restrict__ tabsixi = (fast_alpha_t*) &tabsix[itype*ntypes];
    for (k = 0; k < numneigh; k++) {
      j = neighs[k];
      
      if (j < nall) {
	double delx = xtmp - xx[j].x;
	double dely = ytmp - xx[j].y;
	double delz = ztmp - xx[j].z;
	rsq = delx*delx + dely*dely + delz*delz;
	double tmp_coef3 = qtmp*q[j];
	
	if (rsq < cut_bothsq) {
	  r2inv = 1.0/rsq;
	  
	  forcecoul = 0.0;
	  if (rsq < cut_coulsq) {
	    if (!ncoultablebits || rsq <= tabinnersq) {
	      r = sqrtf(rsq);
	      grij = g_ewald * r;
	      expm2 = exp(-grij*grij);
	      t = 1.0 / (1.0 + EWALD_P*grij);
	      erfc = t * (A1+t*(A2+t*(A3+t*(A4+t*A5)))) * expm2;
	      prefactor = qqrd2e * tmp_coef3/r;
	      forcecoul = prefactor * (erfc + EWALD_F*grij*expm2);
	    } else {
	      itable = *int_rsq & ncoulmask;
	      itable >>= ncoulshiftbits;
	      fraction = (rsq - rtable[itable]) * drtable[itable];
	      table = ftable[itable] + fraction*dftable[itable];
	      forcecoul = tmp_coef3 * table;
	    }
	  }  
	  
	  forcelj = 0.0;
	  if (rsq < cut_ljsq) {
	    r6inv = r2inv*r2inv*r2inv;
	    jtype = type[j] - 1;
	    fast_alpha_t& a = tabsixi[jtype];
	    forcelj = r6inv * (a.lj1*r6inv - a.lj2);
	    if (rsq > cut_lj_innersq) {
	      switch1 = (cut_ljsq-rsq) * (cut_ljsq-rsq) *
		(tmp_coef2 + 2.0*rsq) * tmp_coef1;
	      switch2 = 12.0*rsq * (cut_ljsq-rsq) * 
		(rsq-cut_lj_innersq) * tmp_coef1;
	      philj = r6inv * (a.lj3*r6inv - a.lj4);
	      forcelj = forcelj*switch1 + philj*switch2;
	    }
	  } 
	  
	  fforce = (forcecoul + forcelj) * r2inv;
	  
	  tmpfx += delx*fforce;
	  tmpfy += dely*fforce;
	  tmpfz += delz*fforce;
	  if (NEWTON_PAIR || j < nlocal) {
	    ff[j].x -= delx*fforce;
	    ff[j].y -= dely*fforce;
	    ff[j].z -= delz*fforce;
	  }
	  
	  if (EFLAG) {
	    if (NEWTON_PAIR || j < nlocal) factor = 1.0;
	    else factor = 0.5;
	    if (rsq < cut_coulsq) {
	      if (!ncoultablebits || rsq <= tabinnersq)
		phicoul = prefactor*erfc;
	      else {
		table = etable[itable] + fraction*detable[itable];
		phicoul = tmp_coef3 * table;
	      }
	      eng_coul += factor*phicoul;
	    }
	    if (rsq < cut_ljsq) {
	      fast_alpha_t& a = tabsixi[jtype];
	      philj = r6inv*(a.lj3*r6inv-a.lj4);
	      if (rsq > cut_lj_innersq) {
		switch1 = (cut_ljsq-rsq) * (cut_ljsq-rsq) *
		  (tmp_coef2 + 2.0*rsq) * tmp_coef1;
		philj *= switch1;
	      }
	      eng_vdwl += factor*philj;
	    }
	  }
	  
	  if (VFLAG == 1) {
	    if (NEWTON_PAIR == 0 && j >= nlocal) fforce *= 0.5;
	    virial[0] += delx*delx*fforce;
	    virial[1] += dely*dely*fforce;
	    virial[2] += delz*delz*fforce;
	    virial[3] += delx*dely*fforce;
	    virial[4] += delx*delz*fforce;
	    virial[5] += dely*delz*fforce;
	  }
	}
      } else {
	factor_coul = special_coul[j/nall];
	factor_lj = special_lj[j/nall];
	j %= nall;
	double delx = xtmp - xx[j].x;
	double dely = ytmp - xx[j].y;
	double delz = ztmp - xx[j].z;
	rsq = delx*delx + dely*dely + delz*delz;
	double tmp_coef3 = qtmp*q[j];
	
	if (rsq < cut_bothsq) {
	  r2inv = 1.0/rsq;
	  
	  forcecoul = 0.0;
	  if (rsq < cut_coulsq) {
	    if (!ncoultablebits || rsq <= tabinnersq) {
	      r = sqrtf(rsq);
	      grij = g_ewald * r;
	      expm2 = exp(-grij*grij);
	      t = 1.0 / (1.0 + EWALD_P*grij);
	      erfc = t * (A1+t*(A2+t*(A3+t*(A4+t*A5)))) * expm2;
	      prefactor = qqrd2e * tmp_coef3/r;
	      forcecoul = prefactor * (erfc + EWALD_F*grij*expm2);
	      if (factor_coul < 1.0) {
		forcecoul -= (1.0-factor_coul)*prefactor;
	      }
	    } else {
	      itable = *int_rsq & ncoulmask;
	      itable >>= ncoulshiftbits;
	      fraction = (rsq - rtable[itable]) * drtable[itable];
	      table = ftable[itable] + fraction*dftable[itable];
	      forcecoul = tmp_coef3 * table;
	      if (factor_coul < 1.0) {
		table = ctable[itable] + fraction*dctable[itable];
		prefactor = tmp_coef3 * table;
		forcecoul -= (1.0-factor_coul)*prefactor;
	      }
	    }
	  }  
	  
	  forcelj = 0.0;
	  if (rsq < cut_ljsq) {
	    r6inv = r2inv*r2inv*r2inv;
	    jtype = type[j] - 1;
	    fast_alpha_t& a = tabsixi[jtype];
	    forcelj = r6inv * (a.lj1*r6inv - a.lj2);
	    if (rsq > cut_lj_innersq) {
	      switch1 = (cut_ljsq-rsq) * (cut_ljsq-rsq) *
		(tmp_coef2 + 2.0*rsq) * tmp_coef1;
	      switch2 = 12.0*rsq * (cut_ljsq-rsq) * 
		(rsq-cut_lj_innersq) * tmp_coef1;
	      fast_alpha_t& a = tabsixi[jtype];
	      philj = r6inv * (a.lj3*r6inv - a.lj4);
	      forcelj = forcelj*switch1 + philj*switch2;
	    }
	  }  
	  
	  fforce = (forcecoul + factor_lj*forcelj) * r2inv;
	  
	  tmpfx += delx*fforce;
	  tmpfy += dely*fforce;
	  tmpfz += delz*fforce;
	  if (NEWTON_PAIR || j < nlocal) {
	    ff[j].x -= delx*fforce;
	    ff[j].y -= dely*fforce;
	    ff[j].z -= delz*fforce;
	  }
	  
	  if (EFLAG) {
	    if (NEWTON_PAIR || j < nlocal) factor = 1.0;
	    else factor = 0.5;
	    if (rsq < cut_coulsq) {
	      if (!ncoultablebits || rsq <= tabinnersq)
		phicoul = prefactor*erfc;
	      else {
		table = etable[itable] + fraction*detable[itable];
		phicoul = tmp_coef3 * table;
	      }
	      if (factor_coul < 1.0)
		phicoul -= (1.0-factor_coul)*prefactor;
	      eng_coul += factor*phicoul;
	    }
	    if (rsq < cut_ljsq) {
	      fast_alpha_t& a = tabsixi[jtype];
	      philj = r6inv*(a.lj3*r6inv-a.lj4);
	      if (rsq > cut_lj_innersq) {
		switch1 = (cut_ljsq-rsq) * (cut_ljsq-rsq) *
		  (tmp_coef2 + 2.0*rsq) * tmp_coef1;
		philj *= switch1;
	      }
	      eng_vdwl += factor*factor_lj*philj;
	    }
	  }
	  
	  if (VFLAG == 1) {
	    if (NEWTON_PAIR == 0 && j >= nlocal) fforce *= 0.5;
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
    ff[i].x += tmpfx;
    ff[i].y += tmpfy;
    ff[i].z += tmpfz;
  }
  if (VFLAG == 2) virial_compute();
  free(fast_alpha); fast_alpha = 0;
}
 
}

#endif
