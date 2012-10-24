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
   Contributing author: Mike Brown (SNL)
------------------------------------------------------------------------- */

#include "math.h"
#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "pair_gayberne.h"
#include "math_extra.h"
#include "atom.h"
#include "comm.h"
#include "force.h"
#include "neighbor.h"
#include "update.h"
#include "integrate.h"
#include "respa.h"
#include "memory.h"
#include "error.h"

#define MIN(a,b) ((a) < (b) ? (a) : (b))
#define MAX(a,b) ((a) > (b) ? (a) : (b))

using namespace LAMMPS_NS;

enum{SPHERE_SPHERE,SPHERE_ELLIPSE,ELLIPSE_ELLIPSE};

/* ---------------------------------------------------------------------- */

PairGayBerne::PairGayBerne(LAMMPS *lmp) : Pair(lmp)
{
  single_enable = 0;
}

/* ----------------------------------------------------------------------
   free all arrays
------------------------------------------------------------------------- */

PairGayBerne::~PairGayBerne()
{
  if (allocated) {
    memory->destroy_2d_int_array(setflag);
    memory->destroy_2d_double_array(cutsq);

    memory->destroy_2d_int_array(form);
    memory->destroy_2d_double_array(epsilon);
    memory->destroy_2d_double_array(sigma);
    memory->destroy_2d_double_array(shape);
    memory->destroy_2d_double_array(well);
    memory->destroy_2d_double_array(cut);
    memory->destroy_2d_double_array(lj1);
    memory->destroy_2d_double_array(lj2);
    memory->destroy_2d_double_array(lj3);
    memory->destroy_2d_double_array(lj4);
    memory->destroy_2d_double_array(offset);
    delete [] lshape;
    delete [] setwell;
  }
}

/* ---------------------------------------------------------------------- */

void PairGayBerne::compute(int eflag, int vflag)
{
  int i,j,k,m,numneigh,itype,jtype;
  double one_eng,rsq,r2inv,r6inv,forcelj;
  double fforce[3],ttor[3],rtor[3],r12[3];
  double a1[3][3],b1[3][3],g1[3][3],a2[3][3],b2[3][3],g2[3][3],temp[3][3];
  int *neighs;
  double factor_lj;
  double **f;

  eng_vdwl = 0.0;
  if (vflag) for (i = 0; i < 6; i++) virial[i] = 0.0;

  if (vflag == 2) f = update->f_pair;
  else f = atom->f;
  double **x = atom->x;
  double **quat = atom->quat;
  double **tor = atom->torque;
  int *type = atom->type;
  int nlocal = atom->nlocal;
  int nall = atom->nlocal + atom->nghost;
  double *special_lj = force->special_lj;
  int newton_pair = force->newton_pair;

  // loop over neighbors of my atoms

  for (i = 0; i < nlocal; i++) {
    itype = type[i];
    neighs = neighbor->firstneigh[i];
    numneigh = neighbor->numneigh[i];

    MathExtra::quat_to_mat_trans(quat[i],a1);
    MathExtra::diag_times3(well[itype],a1,temp);
    MathExtra::transpose_times3(a1,temp,b1);
    MathExtra::diag_times3(shape[itype],a1,temp);
    MathExtra::transpose_times3(a1,temp,g1);

    for (k = 0; k < numneigh; k++) {
      j = neighs[k];

      if (j < nall) factor_lj = 1.0;
      else {
        factor_lj = special_lj[j/nall];
        j %= nall;
      }

      // r12 = center to center vector

      r12[0] = x[j][0]-x[i][0];
      r12[1] = x[j][1]-x[i][1];
      r12[2] = x[j][2]-x[i][2];
      rsq = MathExtra::dot3(r12,r12);
      jtype = type[j];

      // compute if less than cutoff

      if (rsq < cutsq[itype][jtype]) {

	switch (form[itype][jtype]) {
	case SPHERE_SPHERE:
	  r2inv = 1.0/rsq;
	  r6inv = r2inv*r2inv*r2inv;
	  forcelj = r6inv * (lj1[itype][jtype]*r6inv - lj2[itype][jtype]);
	  forcelj *= -r2inv;
	  if (eflag) one_eng = 
		       r6inv*(r6inv*lj3[itype][jtype]-lj4[itype][jtype]) -
		       offset[itype][jtype];
	  fforce[0] = r12[0]*forcelj;
	  fforce[1] = r12[1]*forcelj;
	  fforce[2] = r12[2]*forcelj;
	  ttor[0] = ttor[1] = ttor[2] = 0.0;
	  rtor[0] = rtor[1] = rtor[2] = 0.0;
	  break;

	default:
	  MathExtra::quat_to_mat_trans(quat[j],a2);
	  MathExtra::diag_times3(well[jtype],a2,temp);
	  MathExtra::transpose_times3(a2,temp,b2);
	  MathExtra::diag_times3(shape[jtype],a2,temp);
	  MathExtra::transpose_times3(a2,temp,g2);
	  one_eng = gayberne_analytic(i,j,a1,a2,b1,b2,g1,g2,r12,rsq,
				      fforce,ttor,rtor);
	  break;
	}

        fforce[0] *= factor_lj;
	fforce[1] *= factor_lj;
	fforce[2] *= factor_lj;
        ttor[0] *= factor_lj;
	ttor[1] *= factor_lj;
	ttor[2] *= factor_lj;

        f[i][0] += fforce[0];
	f[i][1] += fforce[1];
	f[i][2] += fforce[2];
        tor[i][0] += ttor[0];
	tor[i][1] += ttor[1];
	tor[i][2] += ttor[2];

        if (newton_pair || j < nlocal) {
          rtor[0] *= factor_lj;
	  rtor[1] *= factor_lj;
	  rtor[2] *= factor_lj;
          f[j][0] -= fforce[0];
	  f[j][1] -= fforce[1];
	  f[j][2] -= fforce[2];
          tor[j][0] += rtor[0];
	  tor[j][1] += rtor[1];
	  tor[j][2] += rtor[2];
        }

        if (eflag) {
	  if (newton_pair || j < nlocal) eng_vdwl += factor_lj*one_eng;
	  else eng_vdwl += 0.5*factor_lj*one_eng;
        }

        if (vflag == 1) {
	  if (newton_pair == 0 && j >= nlocal) 
	    for (m = 0; m < 6; m++) fforce[m] *= 0.5;
	  virial[0] += r12[0]*r12[0]*fforce[0];
	  virial[1] += r12[1]*r12[1]*fforce[1];
	  virial[2] += r12[2]*r12[2]*fforce[2];
	  virial[3] += r12[0]*r12[1]*fforce[0];
	  virial[4] += r12[0]*r12[2]*fforce[1];
	  virial[5] += r12[1]*r12[2]*fforce[2];
	}
      }
    }
  }
  if (vflag == 2) virial_compute();
}

/* ----------------------------------------------------------------------
   allocate all arrays
------------------------------------------------------------------------- */

void PairGayBerne::allocate()
{
  allocated = 1;
  int n = atom->ntypes;

  setflag = memory->create_2d_int_array(n+1,n+1,"pair:setflag");
  for (int i = 1; i <= n; i++)
    for (int j = i; j <= n; j++)
      setflag[i][j] = 0;

  cutsq = memory->create_2d_double_array(n+1,n+1,"pair:cutsq");

  form = memory->create_2d_int_array(n+1,n+1,"pair:form");
  epsilon = memory->create_2d_double_array(n+1,n+1,"pair:epsilon");
  sigma = memory->create_2d_double_array(n+1,n+1,"pair:sigma");
  shape = memory->create_2d_double_array(n+1,3,"pair:shape");
  well = memory->create_2d_double_array(n+1,3,"pair:well");
  cut = memory->create_2d_double_array(n+1,n+1,"pair:cut");
  lj1 = memory->create_2d_double_array(n+1,n+1,"pair:lj1");
  lj2 = memory->create_2d_double_array(n+1,n+1,"pair:lj2");
  lj3 = memory->create_2d_double_array(n+1,n+1,"pair:lj3");
  lj4 = memory->create_2d_double_array(n+1,n+1,"pair:lj4");
  offset = memory->create_2d_double_array(n+1,n+1,"pair:offset");

  lshape = new double[n+1];
  setwell = new int[n+1];
  for (int i = 1; i <= n; i++) setwell[i] = 0;
}

/* ----------------------------------------------------------------------
   global settings
------------------------------------------------------------------------- */

void PairGayBerne::settings(int narg, char **arg)
{
  if (narg != 4) error->all("Illegal pair_style command");

  gamma = atof(arg[0]);
  upsilon = atof(arg[1])/2.0;
  mu = atof(arg[2]);
  cut_global = atof(arg[3]);
  
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

void PairGayBerne::coeff(int narg, char **arg)
{
  if (narg < 10 || narg > 11)
    error->all("Incorrect args for pair coefficients");
  if (!allocated) allocate();

  int ilo,ihi,jlo,jhi;
  force->bounds(arg[0],atom->ntypes,ilo,ihi);
  force->bounds(arg[1],atom->ntypes,jlo,jhi);

  double epsilon_one = atof(arg[2]);
  double sigma_one = atof(arg[3]);
  double eia_one = atof(arg[4]);
  double eib_one = atof(arg[5]);
  double eic_one = atof(arg[6]);
  double eja_one = atof(arg[7]);
  double ejb_one = atof(arg[8]);
  double ejc_one = atof(arg[9]);
  
  double cut_one = cut_global;
  if (narg == 11) cut_one = atof(arg[10]);

  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
    for (int j = MAX(jlo,i); j <= jhi; j++) {
      epsilon[i][j] = epsilon_one;
      sigma[i][j] = sigma_one;
      cut[i][j] = cut_one;
      if (eia_one != 0.0 || eib_one != 0.0 || eic_one != 0.0) {
	well[i][0] = pow(eia_one,-1.0/mu);
        well[i][1] = pow(eib_one,-1.0/mu);
	well[i][2] = pow(eic_one,-1.0/mu);
	if (eia_one == 1.0 && eib_one == 1.0 && eic_one == 1.0) setwell[i] = 2;
	else setwell[i] = 1;
      }
      if (eja_one != 0.0 || ejb_one != 0.0 || ejc_one != 0.0) {
	well[j][0] = pow(eja_one,-1.0/mu);
        well[j][1] = pow(ejb_one,-1.0/mu);
	well[j][2] = pow(ejc_one,-1.0/mu);
	if (eja_one == 1.0 && ejb_one == 1.0 && ejc_one == 1.0) setwell[j] = 2;
	else setwell[j] = 1;
      }
      setflag[i][j] = 1;
      count++;
    }
  }

  if (count == 0) error->all("Incorrect args for pair coefficients");
}

/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

double PairGayBerne::init_one(int i, int j)
{
  if (setwell[i] == 0 || setwell[j] == 0)
    error->all("Pair gayberne epsilon a,b,c coeffs are not all set");

  if (setflag[i][j] == 0) {
    epsilon[i][j] = mix_energy(epsilon[i][i],epsilon[j][j],
			       sigma[i][i],sigma[j][j]);
    sigma[i][j] = mix_distance(sigma[i][i],sigma[j][j]);
    cut[i][j] = mix_distance(cut[i][i],cut[j][j]);
  }
  
  lj1[i][j] = 48.0 * epsilon[i][j] * pow(sigma[i][j],12.0);
  lj2[i][j] = 24.0 * epsilon[i][j] * pow(sigma[i][j],6.0);
  lj3[i][j] = 4.0 * epsilon[i][j] * pow(sigma[i][j],12.0);
  lj4[i][j] = 4.0 * epsilon[i][j] * pow(sigma[i][j],6.0);
     
  if (offset_flag) {
    double ratio = sigma[i][j] / cut[i][j];
    offset[i][j] = 4.0 * epsilon[i][j] * (pow(ratio,12.0) - pow(ratio,6.0));
  } else offset[i][j] = 0.0;

  int ishape = 0;
  if (atom->shape[i][0] != atom->shape[i][1] || 
      atom->shape[i][0] != atom->shape[i][2] || 
      atom->shape[i][1] != atom->shape[i][2]) ishape = 1;
  if (setwell[i] == 1) ishape = 1;
  int jshape = 0;
  if (atom->shape[j][0] != atom->shape[j][1] || 
      atom->shape[j][0] != atom->shape[j][2] || 
      atom->shape[j][1] != atom->shape[j][2]) jshape = 1;
  if (setwell[j] == 1) jshape = 1;
  
  if (ishape == 0 && jshape == 0) form[i][j] = SPHERE_SPHERE;
  else if (ishape == 0 || jshape == 0) form[i][j] = SPHERE_ELLIPSE;
  else form[i][j] = ELLIPSE_ELLIPSE;

  form[j][i] = form[i][j];
  epsilon[j][i] = epsilon[i][j];
  sigma[j][i] = sigma[i][j];
  cut[j][i] = cut[i][j];
  lj1[j][i] = lj1[i][j];
  lj2[j][i] = lj2[i][j];
  lj3[j][i] = lj3[i][j];
  lj4[j][i] = lj4[i][j];
  offset[j][i] = offset[i][j];

  return cut[i][j];
}

/* ----------------------------------------------------------------------
   check for necessary atom arrays
------------------------------------------------------------------------- */

void PairGayBerne::init_style()
{
  if (!atom->quat_flag || !atom->torque_flag)
    error->all("Pair gayberne requires atom attributes quat, torque");

  // per-type shape precalculations

  for (int i = 1; i <= atom->ntypes; i++) {
    if (setwell[i]) {
      double *one = atom->shape[i];
      shape[i][0] = one[0]*one[0];
      shape[i][1] = one[1]*one[1];
      shape[i][2] = one[2]*one[2];
      lshape[i] = (one[0]*one[1]+one[2]*one[2])*sqrt(one[0]*one[1]);
    }
  }
}

/* ----------------------------------------------------------------------
   proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairGayBerne::write_restart(FILE *fp)
{
  write_restart_settings(fp);

  int i,j;
  for (i = 1; i <= atom->ntypes; i++) {
    fwrite(&setwell[i],sizeof(int),1,fp);
    if (setwell[i]) fwrite(&well[i][0],sizeof(double),3,fp);
    for (j = i; j <= atom->ntypes; j++) {
      fwrite(&setflag[i][j],sizeof(int),1,fp);
      if (setflag[i][j]) {
        fwrite(&epsilon[i][j],sizeof(double),1,fp);
        fwrite(&sigma[i][j],sizeof(double),1,fp);
        fwrite(&cut[i][j],sizeof(double),1,fp);
      }
    }
  }
}

/* ----------------------------------------------------------------------
   proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairGayBerne::read_restart(FILE *fp)
{
  read_restart_settings(fp);
  allocate();

  int i,j;
  int me = comm->me;
  for (i = 1; i <= atom->ntypes; i++) {
    if (me == 0) fread(&setwell[i],sizeof(int),1,fp);
    MPI_Bcast(&setwell[i],1,MPI_INT,0,world);
    if (setwell[i]) {
      if (me == 0) fread(&well[i][0],sizeof(double),3,fp);
      MPI_Bcast(&well[i][0],3,MPI_DOUBLE,0,world);
    }
    for (j = i; j <= atom->ntypes; j++) {
      if (me == 0) fread(&setflag[i][j],sizeof(int),1,fp);
      MPI_Bcast(&setflag[i][j],1,MPI_INT,0,world);
      if (setflag[i][j]) {
        if (me == 0) {
          fread(&epsilon[i][j],sizeof(double),1,fp);
          fread(&sigma[i][j],sizeof(double),1,fp);
          fread(&cut[i][j],sizeof(double),1,fp);
        }
        MPI_Bcast(&epsilon[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&sigma[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&cut[i][j],1,MPI_DOUBLE,0,world);
      }
    }
  }
}

/* ----------------------------------------------------------------------
   proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairGayBerne::write_restart_settings(FILE *fp)
{
  fwrite(&gamma,sizeof(double),1,fp);
  fwrite(&upsilon,sizeof(double),1,fp);
  fwrite(&mu,sizeof(double),1,fp);
  fwrite(&cut_global,sizeof(double),1,fp);
  fwrite(&mix_flag,sizeof(int),1,fp);
}

/* ----------------------------------------------------------------------
   proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairGayBerne::read_restart_settings(FILE *fp)
{
  int me = comm->me;
  if (me == 0) {
    fread(&gamma,sizeof(double),1,fp);
    fread(&upsilon,sizeof(double),1,fp);
    fread(&mu,sizeof(double),1,fp);
    fread(&cut_global,sizeof(double),1,fp);
    fread(&mix_flag,sizeof(int),1,fp);
  }
  MPI_Bcast(&gamma,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&upsilon,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&mu,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&cut_global,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&mix_flag,1,MPI_INT,0,world);
}

/* ----------------------------------------------------------------------
   compute analytic energy, force (fforce), and torque (ttor & rtor)
   based on rotation matrices a and precomputed matrices b and g
   if newton is off, rtor is not calculated for ghost atoms
------------------------------------------------------------------------- */

double PairGayBerne::gayberne_analytic(const int i,const int j,double a1[3][3],
                                       double a2[3][3], double b1[3][3],
                                       double b2[3][3], double g1[3][3],
                                       double g2[3][3], double *r12,
                                       const double rsq, double *fforce,
                                       double *ttor, double *rtor)
{
  double tempv[3], tempv2[3];
  double temp[3][3];
  double temp1,temp2,temp3;

  int *type = atom->type;
  int newton_pair = force->newton_pair;
  int nlocal = atom->nlocal;

  double r12hat[3];
  MathExtra::normalize3(r12,r12hat);
  double r = sqrt(rsq);

  // compute distance of closest approach

  double g12[3][3];
  MathExtra::plus3(g1,g2,g12);
  double kappa[3];
  MathExtra::mldivide3(g12,r12,kappa,error);

  // tempv = G12^-1*r12hat

  tempv[0] = kappa[0]/r;
  tempv[1] = kappa[1]/r;
  tempv[2] = kappa[2]/r;
  double sigma12 = MathExtra::dot3(r12hat,tempv);
  sigma12 = pow(0.5*sigma12,-0.5);
  double h12 = r-sigma12;

  // energy
  // compute u_r

  double varrho = sigma[type[i]][type[j]]/(h12+gamma*sigma[type[i]][type[j]]);
  double varrho6 = pow(varrho,6.0);
  double varrho12 = varrho6*varrho6;
  double u_r = 4.0*epsilon[type[i]][type[j]]*(varrho12-varrho6);

  // compute eta_12

  double eta = 2.0*lshape[type[i]]*lshape[type[j]];
  double det_g12 = MathExtra::det3(g12);
  eta = pow(eta/det_g12,upsilon);

  // compute chi_12

  double b12[3][3];
  double iota[3];
  MathExtra::plus3(b1,b2,b12);
  MathExtra::mldivide3(b12,r12,iota,error);

  // tempv = G12^-1*r12hat

  tempv[0] = iota[0]/r;
  tempv[1] = iota[1]/r;
  tempv[2] = iota[2]/r;
  double chi = MathExtra::dot3(r12hat,tempv);
  chi = pow(chi*2.0,mu);

  // force
  // compute dUr/dr

  temp1 = (2.0*varrho12*varrho-varrho6*varrho)/sigma[type[i]][type[j]];
  temp1 = temp1*24.0*epsilon[type[i]][type[j]];
  double u_slj = temp1*pow(sigma12,3.0)/2.0;
  double dUr[3];
  temp2 = MathExtra::dot3(kappa,r12hat);
  double uslj_rsq = u_slj/rsq;
  dUr[0] = temp1*r12hat[0]+uslj_rsq*(kappa[0]-temp2*r12hat[0]);
  dUr[1] = temp1*r12hat[1]+uslj_rsq*(kappa[1]-temp2*r12hat[1]);
  dUr[2] = temp1*r12hat[2]+uslj_rsq*(kappa[2]-temp2*r12hat[2]);

  // compute dChi_12/dr

  double dchi[3];
  temp1 = MathExtra::dot3(iota,r12hat);
  temp2 = -4.0/rsq*mu*pow(chi,(mu-1.0)/mu);
  dchi[0] = temp2*(iota[0]-temp1*r12hat[0]);
  dchi[1] = temp2*(iota[1]-temp1*r12hat[1]);
  dchi[2] = temp2*(iota[2]-temp1*r12hat[2]);

  temp1 = -eta*u_r;
  temp2 = eta*chi;
  fforce[0] = temp1*dchi[0]-temp2*dUr[0];
  fforce[1] = temp1*dchi[1]-temp2*dUr[1];
  fforce[2] = temp1*dchi[2]-temp2*dUr[2];

  // torque for particle 1 and 2
  // compute dUr

  tempv[0] = -uslj_rsq*kappa[0];
  tempv[1] = -uslj_rsq*kappa[1];
  tempv[2] = -uslj_rsq*kappa[2];
  MathExtra::row_times3(kappa,g1,tempv2);
  MathExtra::cross3(tempv,tempv2,dUr);
  double dUr2[3];

  if (newton_pair || j < nlocal) {
    MathExtra::row_times3(kappa,g2,tempv2);
    MathExtra::cross3(tempv,tempv2,dUr2);
  }

  // compute d_chi

  MathExtra::row_times3(iota,b1,tempv);
  MathExtra::cross3(tempv,iota,dchi);
  temp1 = -4.0/rsq;
  dchi[0] *= temp1;
  dchi[1] *= temp1;
  dchi[2] *= temp1;
  double dchi2[3];

  if (newton_pair || j < nlocal) {
    MathExtra::row_times3(iota,b2,tempv);
    MathExtra::cross3(tempv,iota,dchi2);
    dchi2[0] *= temp1;
    dchi2[1] *= temp1;
    dchi2[2] *= temp1;
  }

  // compute d_eta

  double deta[3];
  deta[0] = deta[1] = deta[2] = 0.0;
  compute_eta_torque(g12,a1,shape[type[i]],temp);
  temp1 = -eta*upsilon;
  for (int m = 0; m < 3; m++) {
    for (int y = 0; y < 3; y++) tempv[y] = temp1*temp[m][y];
    MathExtra::cross3(a1[m],tempv,tempv2);
    deta[0] += tempv2[0];
    deta[1] += tempv2[1];
    deta[2] += tempv2[2];
  }

  // compute d_eta for particle 2

  double deta2[3];
  if (newton_pair || j < nlocal) {
    deta2[0] = deta2[1] = deta2[2] = 0.0;
    compute_eta_torque(g12,a2,shape[type[j]],temp);
    for (int m = 0; m < 3; m++) {
      for (int y = 0; y < 3; y++) tempv[y] = temp1*temp[m][y];
      MathExtra::cross3(a2[m],tempv,tempv2);
      deta2[0] += tempv2[0];
      deta2[1] += tempv2[1];
      deta2[2] += tempv2[2];
    }
  }

  // torque

  temp1 = u_r*eta;
  temp2 = u_r*chi;
  temp3 = chi*eta;

  ttor[0] = (temp1*dchi[0]+temp2*deta[0]+temp3*dUr[0]) * -1.0;
  ttor[1] = (temp1*dchi[1]+temp2*deta[1]+temp3*dUr[1]) * -1.0;
  ttor[2] = (temp1*dchi[2]+temp2*deta[2]+temp3*dUr[2]) * -1.0;

  if (newton_pair || j < nlocal) {
    rtor[0] = (temp1*dchi2[0]+temp2*deta2[0]+temp3*dUr2[0]) * -1.0;
    rtor[1] = (temp1*dchi2[1]+temp2*deta2[1]+temp3*dUr2[1]) * -1.0;
    rtor[2] = (temp1*dchi2[2]+temp2*deta2[2]+temp3*dUr2[2]) * -1.0;
  }

  return temp1*chi;
}

/* ----------------------------------------------------------------------
   torque contribution from eta
   computes trace in the last doc equation for the torque derivative
   code comes from symbolic solver dump
   m is g12, m2 is a_i, s is the shape for the particle
------------------------------------------------------------------------- */

void PairGayBerne::compute_eta_torque(double m[3][3], double m2[3][3],
                                      double *s, double ans[3][3])
{
  double den = m[1][0]*m[0][2]*m[2][1]-m[0][0]*m[1][2]*m[2][1]-
    m[0][2]*m[2][0]*m[1][1]+m[0][1]*m[2][0]*m[1][2]-
    m[1][0]*m[0][1]*m[2][2]+m[0][0]*m[1][1]*m[2][2];
  
  ans[0][0] = s[0]*(m[1][2]*m[0][1]*m2[0][2]+2.0*m[1][1]*m[2][2]*m2[0][0]-
		    m[1][1]*m2[0][2]*m[0][2]-2.0*m[1][2]*m2[0][0]*m[2][1]+
		    m2[0][1]*m[0][2]*m[2][1]-m2[0][1]*m[0][1]*m[2][2]-
		    m[1][0]*m[2][2]*m2[0][1]+m[2][0]*m[1][2]*m2[0][1]+
		    m[1][0]*m2[0][2]*m[2][1]-m2[0][2]*m[2][0]*m[1][1])/den;
  
  ans[0][1] = s[0]*(m[0][2]*m2[0][0]*m[2][1]-m[2][2]*m2[0][0]*m[0][1]+
		    2.0*m[0][0]*m[2][2]*m2[0][1]-m[0][0]*m2[0][2]*m[1][2]-
		    2.0*m[2][0]*m[0][2]*m2[0][1]+m2[0][2]*m[1][0]*m[0][2]-
		    m[2][2]*m[1][0]*m2[0][0]+m[2][0]*m2[0][0]*m[1][2]+
		    m[2][0]*m2[0][2]*m[0][1]-m2[0][2]*m[0][0]*m[2][1])/den;
  
  ans[0][2] = s[0]*(m[0][1]*m[1][2]*m2[0][0]-m[0][2]*m2[0][0]*m[1][1]-
		    m[0][0]*m[1][2]*m2[0][1]+m[1][0]*m[0][2]*m2[0][1]-
		    m2[0][1]*m[0][0]*m[2][1]-m[2][0]*m[1][1]*m2[0][0]+
		    2.0*m[1][1]*m[0][0]*m2[0][2]-2.0*m[1][0]*m2[0][2]*m[0][1]+
		    m[1][0]*m[2][1]*m2[0][0]+m[2][0]*m2[0][1]*m[0][1])/den;
  
  ans[1][0] = s[1]*(-m[1][1]*m2[1][2]*m[0][2]+2.0*m[1][1]*m[2][2]*m2[1][0]+
		    m[1][2]*m[0][1]*m2[1][2]-2.0*m[1][2]*m2[1][0]*m[2][1]+
		    m2[1][1]*m[0][2]*m[2][1]-m2[1][1]*m[0][1]*m[2][2]-
		    m[1][0]*m[2][2]*m2[1][1]+m[2][0]*m[1][2]*m2[1][1]-
		    m2[1][2]*m[2][0]*m[1][1]+m[1][0]*m2[1][2]*m[2][1])/den;
  
  ans[1][1] = s[1]*(m[0][2]*m2[1][0]*m[2][1]-m[0][1]*m[2][2]*m2[1][0]+
		    2.0*m[2][2]*m[0][0]*m2[1][1]-m2[1][2]*m[0][0]*m[1][2]-
		    2.0*m[2][0]*m2[1][1]*m[0][2]-m[1][0]*m[2][2]*m2[1][0]+
		    m[2][0]*m[1][2]*m2[1][0]+m[1][0]*m2[1][2]*m[0][2]-
		    m[0][0]*m2[1][2]*m[2][1]+m2[1][2]*m[0][1]*m[2][0])/den;
  
  ans[1][2] = s[1]*(m[0][1]*m[1][2]*m2[1][0]-m[0][2]*m2[1][0]*m[1][1]-
		    m[0][0]*m[1][2]*m2[1][1]+m[1][0]*m[0][2]*m2[1][1]+
		    2.0*m[1][1]*m[0][0]*m2[1][2]-m[0][0]*m2[1][1]*m[2][1]+
		    m[0][1]*m[2][0]*m2[1][1]-m2[1][0]*m[2][0]*m[1][1]-
		    2.0*m[1][0]*m[0][1]*m2[1][2]+m[1][0]*m2[1][0]*m[2][1])/den;
  
  ans[2][0] = s[2]*(-m[1][1]*m[0][2]*m2[2][2]+m[0][1]*m[1][2]*m2[2][2]+
		    2.0*m[1][1]*m2[2][0]*m[2][2]-m[0][1]*m2[2][1]*m[2][2]+
		    m[0][2]*m[2][1]*m2[2][1]-2.0*m2[2][0]*m[2][1]*m[1][2]-
		    m[1][0]*m2[2][1]*m[2][2]+m[1][2]*m[2][0]*m2[2][1]-
		    m[1][1]*m[2][0]*m2[2][2]+m[2][1]*m[1][0]*m2[2][2])/den;
  
  ans[2][1] = s[2]*-(m[0][1]*m[2][2]*m2[2][0]-m[0][2]*m2[2][0]*m[2][1]-
		     2.0*m2[2][1]*m[0][0]*m[2][2]+m[1][2]*m2[2][2]*m[0][0]+
		     2.0*m2[2][1]*m[0][2]*m[2][0]+m[1][0]*m2[2][0]*m[2][2]-
		     m[1][0]*m[0][2]*m2[2][2]-m[1][2]*m[2][0]*m2[2][0]+
		     m[0][0]*m2[2][2]*m[2][1]-m2[2][2]*m[0][1]*m[2][0])/den;
  
  ans[2][2] = s[2]*(m[0][1]*m[1][2]*m2[2][0]-m[0][2]*m2[2][0]*m[1][1]-
		    m[0][0]*m[1][2]*m2[2][1]+m[1][0]*m[0][2]*m2[2][1]-
		    m[1][1]*m[2][0]*m2[2][0]-m[2][1]*m2[2][1]*m[0][0]+
		    2.0*m[1][1]*m2[2][2]*m[0][0]+m[2][1]*m[1][0]*m2[2][0]+
		    m[2][0]*m[0][1]*m2[2][1]-2.0*m2[2][2]*m[1][0]*m[0][1])/den;
}
