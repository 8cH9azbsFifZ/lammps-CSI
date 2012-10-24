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

#ifndef PAIR_LJ_CUT_OPT_H
#define PAIR_LJ_CUT_OPT_H

#include "stdlib.h"
#include "pair_lj_cut.h"
#include "atom.h"
#include "update.h"
#include "force.h"
#include "neighbor.h"

namespace LAMMPS_NS {

class PairLJCutOpt : public PairLJCut {
 public:
  PairLJCutOpt(class LAMMPS *);
  void compute(int, int);

 private:
  template < int EFLAG, int VFLAG, int NEWTON_PAIR > void eval();
};

template < int EFLAG, int VFLAG, int NEWTON_PAIR >
void PairLJCutOpt::eval()
{
  typedef struct { double x,y,z; } vec3_t;
  
  typedef struct {
    double cutsq,lj1,lj2,lj3,lj4,offset;
    double _pad[2];
  } fast_alpha_t;
  
  double** __restrict__ f;
  
  eng_vdwl = 0.0;
  if (VFLAG) for (int i = 0; i < 6; i++) virial[i] = 0.0;
  
  if (VFLAG == 2) f = update->f_pair;
  else f = atom->f;
  
  double** __restrict__ x = atom->x;
  int* __restrict__ type = atom->type;
  int nlocal = atom->nlocal;
  int nall = atom->nlocal + atom->nghost;
  double* __restrict__ special_lj = force->special_lj;
  int** __restrict__ firstneigh = neighbor->firstneigh;
  int* __restrict__ num = neighbor->numneigh;
  
  vec3_t* __restrict__ xx = (vec3_t*)x[0];
  vec3_t* __restrict__ ff = (vec3_t*)f[0];
  
  int ntypes = atom->ntypes;
  int ntypes2 = ntypes*ntypes;
  
  fast_alpha_t* __restrict__ fast_alpha = 
    (fast_alpha_t*) malloc(ntypes2*sizeof(fast_alpha_t));
  for( int i = 0; i < ntypes; i++) for( int j = 0; j < ntypes; j++) {
    fast_alpha_t& a = fast_alpha[i*ntypes+j];
    a.cutsq = cutsq[i+1][j+1];
    a.lj1 = lj1[i+1][j+1];
    a.lj2 = lj2[i+1][j+1];
    a.lj3 = lj3[i+1][j+1];
    a.lj4 = lj4[i+1][j+1];
    a.offset = offset[i+1][j+1];
  }
  fast_alpha_t* __restrict__ tabsix = fast_alpha;
  
  // loop over neighbors of my atoms
  
  for (int i = 0; i < nlocal; i++) {
    double xtmp = xx[i].x;
    double ytmp = xx[i].y;
    double ztmp = xx[i].z;
    int itype = type[i] - 1;
    int* __restrict__ neighs = firstneigh[i];
    int numneigh = num[i];
    
    double tmpfx = 0.0;
    double tmpfy = 0.0;
    double tmpfz = 0.0;
    
    fast_alpha_t* __restrict__ tabsixi = (fast_alpha_t*)&tabsix[itype*ntypes];
    
    for (int k = 0; k < numneigh; k++) {
      int j = neighs[k];
      double factor_lj;
      
      if (j < nall) {
	double delx = xtmp - xx[j].x;
	double dely = ytmp - xx[j].y;
	double delz = ztmp - xx[j].z;
	double rsq = delx*delx + dely*dely + delz*delz;
	
	int jtype = type[j] - 1;
	
	fast_alpha_t& a = tabsixi[jtype];
	
	if (rsq < a.cutsq) {
	  double r2inv = 1.0/rsq;
	  double r6inv = r2inv*r2inv*r2inv;
	  double forcelj = r6inv * (a.lj1*r6inv - a.lj2);
	  if (EFLAG) {
	    double philj = r6inv*(a.lj3*r6inv-a.lj4) - a.offset;
	    if (NEWTON_PAIR || j < nlocal) eng_vdwl += philj;
	    else eng_vdwl += 0.5*philj;
	  }
	  double fforce = forcelj*r2inv;
	  
	  tmpfx += delx*fforce;
	  tmpfy += dely*fforce;
	  tmpfz += delz*fforce;
	  if (NEWTON_PAIR || j < nlocal) {
	    ff[j].x -= delx*fforce;
	    ff[j].y -= dely*fforce;
	    ff[j].z -= delz*fforce;
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
	
	factor_lj = special_lj[j/nall];
	j = j % nall;
	double delx = xtmp - xx[j].x;
	double dely = ytmp - xx[j].y;
	double delz = ztmp - xx[j].z;
	double rsq = delx*delx + dely*dely + delz*delz;
	
	int jtype1 = type[j];
	int jtype = jtype1 - 1;
	
	fast_alpha_t& a = tabsixi[jtype];
	if (rsq < a.cutsq) {
	  
	  double r2inv = 1.0/rsq;
	  double r6inv = r2inv*r2inv*r2inv;
	  
	  fast_alpha_t& a = tabsixi[jtype];
	  double forcelj = r6inv * (a.lj1*r6inv - a.lj2);
	  if (EFLAG) {
	    double philj = r6inv*(a.lj3*r6inv-a.lj4) - a.offset;
	    if (NEWTON_PAIR || j < nlocal) eng_vdwl += factor_lj*philj;
	    else eng_vdwl += 0.5*factor_lj*philj;
	  }
	  
	  double fforce = factor_lj*forcelj*r2inv;
	  
	  tmpfx += delx*fforce;
	  tmpfy += dely*fforce;
	  tmpfz += delz*fforce;
	  if (NEWTON_PAIR || j < nlocal) {
	    ff[j].x -= delx*fforce;
	    ff[j].y -= dely*fforce;
	    ff[j].z -= delz*fforce;
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
