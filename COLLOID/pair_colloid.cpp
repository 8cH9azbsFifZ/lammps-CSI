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
   Contributing author: Pieter in 't Veld (SNL)
------------------------------------------------------------------------- */

#include "math.h"
#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "pair_colloid.h"
#include "atom.h"
#include "comm.h"
#include "force.h"
#include "neighbor.h"
#include "update.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;

#define MIN(a,b) ((a) < (b) ? (a) : (b))
#define MAX(a,b) ((a) > (b) ? (a) : (b))

enum{SMALL_SMALL,SMALL_LARGE,LARGE_LARGE};

/* ---------------------------------------------------------------------- */

PairColloid::PairColloid(LAMMPS *lmp) : Pair(lmp) {}

/* ---------------------------------------------------------------------- */

PairColloid::~PairColloid()
{
  if (allocated) {
    memory->destroy_2d_int_array(setflag);
    memory->destroy_2d_double_array(cutsq);

    memory->destroy_2d_int_array(form);
    memory->destroy_2d_double_array(a12);
    memory->destroy_2d_double_array(sigma);
    memory->destroy_2d_double_array(d1);
    memory->destroy_2d_double_array(d2);
    memory->destroy_2d_double_array(a1);
    memory->destroy_2d_double_array(a2);
    memory->destroy_2d_double_array(cut);		
    memory->destroy_2d_double_array(offset);
    memory->destroy_2d_double_array(sigma3);
    memory->destroy_2d_double_array(sigma6);
    memory->destroy_2d_double_array(lj1);
    memory->destroy_2d_double_array(lj2);
    memory->destroy_2d_double_array(lj3);
    memory->destroy_2d_double_array(lj4);
  }
}

/* ---------------------------------------------------------------------- */

void PairColloid::compute(int eflag, int vflag)
{
  int i,j,k,numneigh,itype,jtype;
  double xtmp,ytmp,ztmp,delx,dely,delz;
  double rsq,r,fforce,forcelj,factor_lj,phi;
  double r2inv,r6inv,c1,c2,aij,s2,s6,fR,dUR,dUA;
  double K[9],h[4],g[4];
  int *neighs;
  double **f;

  eng_vdwl = 0.0;
  if (vflag) for (i = 0; i < 6; i++) virial[i] = 0.0;

  if (vflag == 2) f = update->f_pair;
  else f = atom->f;
  double **x = atom->x;
  int *type = atom->type;
  int nlocal = atom->nlocal;
  int nall = atom->nlocal + atom->nghost;
  double *special_lj = force->special_lj;
  int newton_pair = force->newton_pair;

  // loop over neighbors of my atoms

  for (i = 0; i < nlocal; ++i) {
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    itype = type[i];
    neighs = neighbor->firstneigh[i];
    numneigh = neighbor->numneigh[i];
    
    for (k = 0; k < numneigh; k++) {
      j = neighs[k];

      if (j < nall) factor_lj = 1.0;
      else {
	factor_lj = special_lj[j/nall];
	j %= nall;
      }
     
      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx*delx + dely*dely + delz*delz;
      jtype = type[j];

      if (rsq >= cutsq[itype][jtype]) continue;

      switch (form[itype][jtype]) {
      case SMALL_SMALL: 
	r2inv = 1.0/rsq;
	r6inv = r2inv*r2inv*r2inv;
	forcelj = r6inv * (lj1[itype][jtype]*r6inv - lj2[itype][jtype]);
	fforce = factor_lj*forcelj*r2inv;
	if (eflag) phi = r6inv*(r6inv*lj3[itype][jtype]-lj4[itype][jtype]) -
		     offset[itype][jtype];
	break;
	
      case SMALL_LARGE:
	c2 = a2[itype][jtype];
	K[1] = c2*c2;
	K[2] = rsq;
	K[0] = K[1] - rsq;
	K[4] = rsq*rsq;
	K[3] = K[1] - K[2];
	K[3] *= K[3]*K[3];
	K[6] = K[3]*K[3];
	fR = sigma3[itype][jtype]*a12[itype][jtype]*c2*K[1]/K[3];
	fforce = 4.0/15.0*sqrt(rsq)*fR*factor_lj * 
	  (2.0*(K[1]+K[2]) * (K[1]*(5.0*K[1]+22.0*K[2])+5.0*K[4]) * 
	   sigma6[itype][jtype]/K[6]-5.0) / K[0];
	if (eflag) 
	  phi = 2.0/9.0*fR * 
	    (1.0-(K[1]*(K[1]*(K[1]/3.0+3.0*K[2])+4.2*K[4])+K[2]*K[4]) *
	     sigma6[itype][jtype]/K[6]) - offset[itype][jtype];
	break;

      case LARGE_LARGE:
	r = sqrt(rsq);
	c1 = a1[itype][jtype];
	c2 = a2[itype][jtype];
	K[0] = c1*c2;
	K[1] = c1+c2;
	K[2] = c1-c2;
	K[3] = K[1]+r;
	K[4] = K[1]-r;
	K[5] = K[2]+r;
	K[6] = K[2]-r;
	K[7] = 1.0/(K[3]*K[4]);
	K[8] = 1.0/(K[5]*K[6]);
	g[0] = pow(K[3],-7.0);
	g[1] = pow(K[4],-7.0);
	g[2] = pow(K[5],-7.0);
	g[3] = pow(K[6],-7.0);
	h[0] = ((K[3]+5.0*K[1])*K[3]+30.0*K[0])*g[0];
	h[1] = ((K[4]+5.0*K[1])*K[4]+30.0*K[0])*g[1];
	h[2] = ((K[5]+5.0*K[2])*K[5]-30.0*K[0])*g[2];
	h[3] = ((K[6]+5.0*K[2])*K[6]-30.0*K[0])*g[3];
	g[0] *= 42.0*K[0]/K[3]+6.0*K[1]+K[3];
	g[1] *= 42.0*K[0]/K[4]+6.0*K[1]+K[4];
	g[2] *= -42.0*K[0]/K[5]+6.0*K[2]+K[5];
	g[3] *= -42.0*K[0]/K[6]+6.0*K[2]+K[6];
	
	fR = a12[itype][jtype]*sigma6[itype][jtype]/r/37800.0;
	phi = fR * (h[0]-h[1]-h[2]+h[3]);
	dUR = phi/r + 5.0*fR*(g[0]+g[1]-g[2]-g[3]);
	dUA = -a12[itype][jtype]/3.0*r*((2.0*K[0]*K[7]+1.0)*K[7] + 
					(2.0*K[0]*K[8]-1.0)*K[8]);
	fforce = factor_lj * (dUR+dUA)/r;
	if (eflag)
	  phi += a12[itype][jtype]/6.0*(2.0*K[0]*(K[7]+K[8])-log(K[8]/K[7])) - 
	    offset[itype][jtype];
	break;
      }
      
      if (eflag) {
	if (newton_pair || j < nlocal) eng_vdwl += factor_lj*phi;
	else eng_vdwl += 0.5*factor_lj*phi;
      }

      f[i][0] += delx*fforce;
      f[i][1] += dely*fforce;
      f[i][2] += delz*fforce;
      if (newton_pair || j < nlocal) {
	f[j][0] -= delx*fforce;
	f[j][1] -= dely*fforce;
	f[j][2] -= delz*fforce;
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
  if (vflag == 2) virial_compute();
}

/* ----------------------------------------------------------------------
   allocate all arrays 
------------------------------------------------------------------------- */

void PairColloid::allocate()
{
  allocated = 1;
  int n = atom->ntypes;

  setflag = memory->create_2d_int_array(n+1,n+1,"pair:setflag");
  for (int i = 1; i <= n; i++)
    for (int j = i; j <= n; j++)
      setflag[i][j] = 0;

  cutsq = memory->create_2d_double_array(n+1,n+1,"pair:cutsq");

  form = memory->create_2d_int_array(n+1,n+1,"pair:form");
  a12 = memory->create_2d_double_array(n+1,n+1,"pair:a12");
  sigma = memory->create_2d_double_array(n+1,n+1,"pair:sigma");
  d1 = memory->create_2d_double_array(n+1,n+1,"pair:d1");
  d2 = memory->create_2d_double_array(n+1,n+1,"pair:d2");
  a1 = memory->create_2d_double_array(n+1,n+1,"pair:a1");
  a2 = memory->create_2d_double_array(n+1,n+1,"pair:a2");
  cut = memory->create_2d_double_array(n+1,n+1,"pair:cut");
  offset = memory->create_2d_double_array(n+1,n+1,"pair:offset");
  sigma3 = memory->create_2d_double_array(n+1,n+1,"pair:sigma3");
  sigma6 = memory->create_2d_double_array(n+1,n+1,"pair:sigma6");
  lj1 = memory->create_2d_double_array(n+1,n+1,"pair:lj1");
  lj2 = memory->create_2d_double_array(n+1,n+1,"pair:lj2");
  lj3 = memory->create_2d_double_array(n+1,n+1,"pair:lj3");
  lj4 = memory->create_2d_double_array(n+1,n+1,"pair:lj4");
}

/* ----------------------------------------------------------------------
   global settings 
------------------------------------------------------------------------- */

void PairColloid::settings(int narg, char **arg)
{
  if (narg != 1) error->all("Illegal pair_style command");

  cut_global = atof(arg[0]);

  // reset cutoffs that have been explicitly set

  if (allocated) {
    int i,j;
    for (i = 1; i <= atom->ntypes; i++)
      for (j = i+1; j <= atom->ntypes; j++)
	if (setflag[i][j]) cut[i][j] = cut_global;
  }
}

/* ----------------------------------------------------------------------
   set coeffs for one or more type pairs
------------------------------------------------------------------------- */

void PairColloid::coeff(int narg, char **arg)
{
  if (narg < 6 || narg > 7) error->all("Incorrect args for pair coefficients");
  if (!allocated) allocate();

  int ilo,ihi,jlo,jhi;
  force->bounds(arg[0],atom->ntypes,ilo,ihi);
  force->bounds(arg[1],atom->ntypes,jlo,jhi);

  double a12_one = atof(arg[2]);
  double sigma_one = atof(arg[3]);
  double d1_one = atof(arg[4]);
  double d2_one = atof(arg[5]);

  double cut_one = cut_global;
  if (narg == 7) cut_one = atof(arg[6]);

  if (d1_one < 0.0 || d2_one < 0.0) 
    error->all("Invalid d1 or d2 value for pair colloid coeff");

  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
    for (int j = MAX(jlo,i); j <= jhi; j++) {
      a12[i][j] = a12_one;
      sigma[i][j] = sigma_one;
      if (i == j && d1_one != d2_one)
	error->all("Invalid d1 or d2 value for pair colloid coeff");
      d1[i][j] = d1_one;
      d2[i][j] = d2_one;
      cut[i][j] = cut_one;
      setflag[i][j] = 1;
      count++;
    }
  }

  if (count == 0) error->all("Incorrect args for pair coefficients");
}

/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

double PairColloid::init_one(int i, int j)
{
  if (setflag[i][j] == 0) {
    a12[i][j] = mix_energy(a12[i][i],a12[j][j],sigma[i][i],sigma[j][j]);
    sigma[i][j] = mix_distance(sigma[i][i],sigma[j][j]);
    d1[i][j] = mix_distance(d1[i][i],d1[j][j]);
    d2[i][j] = mix_distance(d2[i][i],d2[j][j]);
    cut[i][j] = mix_distance(cut[i][i],cut[j][j]);
  }

  sigma3[i][j] = sigma[i][j]*sigma[i][j]*sigma[i][j];
  sigma6[i][j] = sigma3[i][j]*sigma3[i][j];

  if (d1[i][j] == 0.0 && d2[i][j] == 0.0) form[i][j] = SMALL_SMALL;
  else if (d1[i][j] == 0.0 || d2[i][j] == 0.0) form[i][j] = SMALL_LARGE;
  else form[i][j] = LARGE_LARGE;

  // for SMALL_SMALL, a1/a2 do not need to be set
  // for SMALL_LARGE, a1 does not need to be set, a2 = radius for i,j and j,i
  // for LARGE_LARGE, a1/a2 are radii, swap them for j,i

  if (form[i][j] == SMALL_LARGE) {
    if (d1[i][j] > 0.0) a2[i][j] = 0.5*d1[i][j];
    else a2[i][j] = 0.5*d2[i][j];
    a2[j][i] = a2[i][j];
  } else if (form[i][j] == LARGE_LARGE) {
    a2[j][i] = a1[i][j] = 0.5*d1[i][j];
    a1[j][i] = a2[i][j] = 0.5*d2[i][j];
  }

  form[j][i] = form[i][j];
  a12[j][i] = a12[i][j];
  sigma[j][i] = sigma[i][j];
  sigma3[j][i] = sigma3[i][j];
  sigma6[j][i] = sigma6[i][j];
  cut[j][i] = cut[i][j];
  cutsq[j][i] = cutsq[i][j] = cut[i][j] * cut[i][j];

  double epsilon = a12[i][j]/144.0;
  lj1[j][i] = lj1[i][j] = 48.0 * epsilon * sigma6[i][j] * sigma6[i][j];
  lj2[j][i] = lj2[i][j] = 24.0 * epsilon * sigma6[i][j];
  lj3[j][i] = lj3[i][j] = 4.0 * epsilon * sigma6[i][j] * sigma6[i][j];
  lj4[j][i] = lj4[i][j] = 4.0 * epsilon * sigma6[i][j];

  offset[j][i] = offset[i][j] = 0.0;
  if (offset_flag) {
    One one;
    single(0,0,i,j,cutsq[i][j],0.0,1.0,1,one);
    offset[j][i] = offset[i][j] = one.eng_vdwl;
  }

  return cut[i][j];
}

/* ----------------------------------------------------------------------
   proc 0 writes to restart file 
------------------------------------------------------------------------- */

void PairColloid::write_restart(FILE *fp)
{
  write_restart_settings(fp);

  int i,j,flag;
  double d;
  for (i = 1; i <= atom->ntypes; i++)
    for (j = i; j <= atom->ntypes; j++) {
      fwrite(&setflag[i][j],sizeof(int),1,fp);
      if (setflag[i][j]) {
	fwrite(&a12[i][j],sizeof(double),1,fp);
	fwrite(&sigma[i][j],sizeof(double),1,fp);
	fwrite(&d1[i][j],sizeof(double),1,fp);
	fwrite(&d2[i][j],sizeof(double),1,fp);
	fwrite(&cut[i][j],sizeof(double),1,fp);
      }
    }
}

/* ----------------------------------------------------------------------
   proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairColloid::read_restart(FILE *fp)
{
  read_restart_settings(fp);
  allocate();

  int i,j;
  int me = comm->me;
  double d;
  for (i = 1; i <= atom->ntypes; i++)
    for (j = i; j <= atom->ntypes; j++) {
      if (me == 0) fread(&setflag[i][j],sizeof(int),1,fp);
      MPI_Bcast(&setflag[i][j],1,MPI_INT,0,world);
      if (setflag[i][j]) {
	if (me == 0) {
	  fread(&a12[i][j],sizeof(double),1,fp);
	  fread(&sigma[i][j],sizeof(double),1,fp);
	  fread(&d1[i][j],sizeof(double),1,fp);
	  fread(&d2[i][j],sizeof(double),1,fp);
	  fread(&cut[i][j],sizeof(double),1,fp);
	}
	MPI_Bcast(&a12[i][j],1,MPI_DOUBLE,0,world);
	MPI_Bcast(&sigma[i][j],1,MPI_DOUBLE,0,world);
	MPI_Bcast(&d1[i][j],1,MPI_DOUBLE,0,world);
	MPI_Bcast(&d2[i][j],1,MPI_DOUBLE,0,world);
	MPI_Bcast(&cut[i][j],1,MPI_DOUBLE,0,world);
      }
    }
}

/* ----------------------------------------------------------------------
   proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairColloid::write_restart_settings(FILE *fp)
{
  fwrite(&cut_global,sizeof(double),1,fp);
  fwrite(&offset_flag,sizeof(int),1,fp);
  fwrite(&mix_flag,sizeof(int),1,fp);
}

/* ----------------------------------------------------------------------
   proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairColloid::read_restart_settings(FILE *fp)
{
  int me = comm->me;
  if (me == 0) {
    fread(&cut_global,sizeof(double),1,fp);
    fread(&offset_flag,sizeof(int),1,fp);
    fread(&mix_flag,sizeof(int),1,fp);
  }
  MPI_Bcast(&cut_global,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&offset_flag,1,MPI_INT,0,world);
  MPI_Bcast(&mix_flag,1,MPI_INT,0,world);
}

/* ---------------------------------------------------------------------- */

void PairColloid::single(int i, int j, int itype, int jtype, double rsq,
		       double factor_coul, double factor_lj, int eflag,
		       One &one)
{
  double K[9],h[4],g[4];
  double r,r2inv,r6inv,forcelj,c1,c2,aij,s6,phi,fR,dUR,dUA;

  switch (form[itype][jtype]) {
  case SMALL_SMALL:
    r2inv = 1.0/rsq;
    r6inv = r2inv*r2inv*r2inv;
    forcelj = r6inv * (lj1[itype][jtype]*r6inv - lj2[itype][jtype]);
    one.fforce = factor_lj*forcelj*r2inv;
    if (eflag) phi = r6inv*(r6inv*lj3[itype][jtype]-lj4[itype][jtype]) -
		     offset[itype][jtype];
    break;

  case SMALL_LARGE:
    r = sqrt(rsq);
    c2 = a2[itype][jtype];
    K[1] = c2*c2;
    K[2] = rsq;
    K[4] = rsq*rsq;
    K[3] = K[1] - K[2];
    K[3] *= K[3]*K[3];
    K[6] = K[3]*K[3];
    fR = sigma3[itype][jtype]*a12[itype][jtype]*c2*K[1]/K[3];
    one.fforce = 4.0/15.0*r*fR*factor_lj * 
      (2.0*(K[1]+K[2])*(K[1]*(5.0*K[1]+22.0*K[2])+5.0*K[4]) * 
       sigma6[itype][jtype]/K[6] - 5.0)/K[0];
    if (eflag) 
      phi = 2.0/9.0*fR * 
	(1.0-(K[1]*(K[1]*(K[1]/3.0+3.0*K[2])+4.2*K[4])+K[2]*K[4]) * 
	 sigma6[itype][jtype]/K[6]) - offset[itype][jtype];
    break;

  case LARGE_LARGE:
    r = sqrt(rsq);
    c1 = a1[itype][jtype];
    c2 = a2[itype][jtype];
    K[0] = c1*c2;
    K[1] = c1+c2;
    K[2] = c1-c2;
    K[3] = K[1]+r;
    K[4] = K[1]-r;
    K[5] = K[2]+r;
    K[6] = K[2]-r;
    K[7] = 1.0/(K[3]*K[4]);
    K[8] = 1.0/(K[5]*K[6]);
    g[0] = pow(K[3],-7.0);
    g[1] = pow(K[4],-7.0);
    g[2] = pow(K[5],-7.0);
    g[3] = pow(K[6],-7.0);
    h[0] = ((K[3]+5.0*K[1])*K[3]+30.0*K[0])*g[0];
    h[1] = ((K[4]+5.0*K[1])*K[4]+30.0*K[0])*g[1];
    h[2] = ((K[5]+5.0*K[2])*K[5]-30.0*K[0])*g[2];
    h[3] = ((K[6]+5.0*K[2])*K[6]-30.0*K[0])*g[3];
    g[0] *= 42.0*K[0]/K[3]+6.0*K[1]+K[3];
    g[1] *= 42.0*K[0]/K[4]+6.0*K[1]+K[4];
    g[2] *= -42.0*K[0]/K[5]+6.0*K[2]+K[5];
    g[3] *= -42.0*K[0]/K[6]+6.0*K[2]+K[6];
    
    fR = a12[itype][jtype]*sigma6[itype][jtype]/r/37800.0;
    phi = fR * (h[0]-h[1]-h[2]+h[3]);
    dUR = phi/r + 5.0*fR*(g[0]+g[1]-g[2]-g[3]);
    dUA = -a12[itype][jtype]/3.0*r*((2.0*K[0]*K[7]+1.0)*K[7] + 
				    (2.0*K[0]*K[8]-1.0)*K[8]);
    one.fforce = factor_lj*(dUR+dUA)/r;
    
    if (eflag)
      phi += a12[itype][jtype]/6.0*(2.0*K[0]*(K[7]+K[8])-log(K[8]/K[7])) - 
	offset[itype][jtype];
    break;
  }

  if (eflag) {
    one.eng_vdwl = factor_lj*phi;
    one.eng_coul = 0.0;
  }
}
