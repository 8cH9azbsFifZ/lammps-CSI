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

#include "math.h"
#include "stdlib.h"
#include "pair_lj_cut_coul_debye.h"
#include "atom.h"
#include "neighbor.h"
#include "force.h"
#include "comm.h"
#include "update.h"
#include "error.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

PairLJCutCoulDebye::PairLJCutCoulDebye(LAMMPS *lmp) : PairLJCutCoulCut(lmp) {}

/* ---------------------------------------------------------------------- */

void PairLJCutCoulDebye::compute(int eflag, int vflag)
{
  int i,j,k,numneigh,itype,jtype;
  double qtmp,xtmp,ytmp,ztmp,delx,dely,delz;
  double rsq,r2inv,r6inv,forcecoul,forcelj,fforce,factor_coul,factor_lj;
  double factor,phicoul,philj;
  double r,rinv,screening;
  int *neighs;
  double **f;

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
      jtype = type[j];

      if (rsq < cutsq[itype][jtype]) {
	r2inv = 1.0/rsq;

	if (rsq < cut_coulsq[itype][jtype]) {
	  r = sqrt(rsq);
	  rinv = 1.0/r;
	  screening = exp(-kappa*r);
	  forcecoul = qqrd2e * qtmp*q[j] * screening * (kappa + rinv);
	} else forcecoul = 0.0;

	if (rsq < cut_ljsq[itype][jtype]) {
	  r6inv = r2inv*r2inv*r2inv;
	  forcelj = r6inv * (lj1[itype][jtype]*r6inv - lj2[itype][jtype]);
	} else forcelj = 0.0;

	fforce = (factor_coul*forcecoul + factor_lj*forcelj) * r2inv;

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
	  if (rsq < cut_coulsq[itype][jtype]) {
	    phicoul = qqrd2e * qtmp*q[j] * rinv * screening;
	    eng_coul += factor*factor_coul*phicoul;
	  }
	  if (rsq < cut_ljsq[itype][jtype]) {
	    philj = r6inv*(lj3[itype][jtype]*r6inv-lj4[itype][jtype]) -
	      offset[itype][jtype];
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

/* ----------------------------------------------------------------------
   global settings 
------------------------------------------------------------------------- */

void PairLJCutCoulDebye::settings(int narg, char **arg)
{
  if (narg < 2 || narg > 3) error->all("Illegal pair_style command");

  kappa = atof(arg[0]);
  cut_lj_global = atof(arg[1]);
  if (narg == 2) cut_coul_global = cut_lj_global;
  else cut_coul_global = atof(arg[2]);

  // reset cutoffs that were previously set from data file

  if (allocated) {
    int i,j;
    for (i = 1; i <= atom->ntypes; i++)
      for (j = i+1; j <= atom->ntypes; j++)
	if (setflag[i][j] == 1) {
	  cut_lj[i][j] = cut_lj_global;
	  cut_coul[i][j] = cut_coul_global;
	}
  }
}

/* ----------------------------------------------------------------------
   proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairLJCutCoulDebye::write_restart(FILE *fp)
{
  write_restart_settings(fp);

  int i,j;
  for (i = 1; i <= atom->ntypes; i++)
    for (j = i; j <= atom->ntypes; j++) {
      fwrite(&setflag[i][j],sizeof(int),1,fp);
      if (setflag[i][j]) {
	fwrite(&epsilon[i][j],sizeof(double),1,fp);
	fwrite(&sigma[i][j],sizeof(double),1,fp);
	fwrite(&cut_lj[i][j],sizeof(double),1,fp);
	fwrite(&cut_coul[i][j],sizeof(double),1,fp);
      }
    }
}

/* ----------------------------------------------------------------------
   proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairLJCutCoulDebye::read_restart(FILE *fp)
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
	  fread(&cut_lj[i][j],sizeof(double),1,fp);
	  fread(&cut_coul[i][j],sizeof(double),1,fp);
	}
	MPI_Bcast(&epsilon[i][j],1,MPI_DOUBLE,0,world);
	MPI_Bcast(&sigma[i][j],1,MPI_DOUBLE,0,world);
	MPI_Bcast(&cut_lj[i][j],1,MPI_DOUBLE,0,world);
	MPI_Bcast(&cut_coul[i][j],1,MPI_DOUBLE,0,world);
      }
    }
}

/* ----------------------------------------------------------------------
   proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairLJCutCoulDebye::write_restart_settings(FILE *fp)
{
  fwrite(&cut_lj_global,sizeof(double),1,fp);
  fwrite(&cut_coul_global,sizeof(double),1,fp);
  fwrite(&kappa,sizeof(double),1,fp);
  fwrite(&offset_flag,sizeof(int),1,fp);
  fwrite(&mix_flag,sizeof(int),1,fp);
}

/* ----------------------------------------------------------------------
   proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairLJCutCoulDebye::read_restart_settings(FILE *fp)
{
  if (comm->me == 0) {
    fread(&cut_lj_global,sizeof(double),1,fp);
    fread(&cut_coul_global,sizeof(double),1,fp);
    fread(&kappa,sizeof(double),1,fp);
    fread(&offset_flag,sizeof(int),1,fp);
    fread(&mix_flag,sizeof(int),1,fp);
  }
  MPI_Bcast(&cut_lj_global,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&cut_coul_global,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&kappa,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&offset_flag,1,MPI_INT,0,world);
  MPI_Bcast(&mix_flag,1,MPI_INT,0,world);
}

/* ---------------------------------------------------------------------- */

void PairLJCutCoulDebye::single(int i, int j, int itype, int jtype, double rsq,
				double factor_coul, double factor_lj,
				int eflag, One &one)
{
  double r2inv,r6inv,r,rinv,screening,forcecoul,forcelj,phicoul,philj;

  r2inv = 1.0/rsq;
  if (rsq < cut_coulsq[itype][jtype]) {
    r = sqrt(rsq);
    rinv = 1.0/r;
    screening = exp(-kappa*r);
    forcecoul = force->qqrd2e * atom->q[i]*atom->q[j] *
      screening * (kappa + rinv);
  } else forcecoul = 0.0;
  if (rsq < cut_ljsq[itype][jtype]) {
    r6inv = r2inv*r2inv*r2inv;
    forcelj = r6inv * (lj1[itype][jtype]*r6inv - lj2[itype][jtype]);
  } else forcelj = 0.0;
  one.fforce = (factor_coul*forcecoul + factor_lj*forcelj) * r2inv;
  
  if (eflag) {
    if (rsq < cut_coulsq[itype][jtype]) {
      phicoul = force->qqrd2e * atom->q[i]*atom->q[j] * rinv * screening;
      one.eng_coul = factor_coul*phicoul;
    } else one.eng_coul = 0.0;
    if (rsq < cut_ljsq[itype][jtype]) {
      philj = r6inv*(lj3[itype][jtype]*r6inv-lj4[itype][jtype]) -
	offset[itype][jtype];
      one.eng_vdwl = factor_lj*philj;
    } else one.eng_vdwl = 0.0;
  }
}
