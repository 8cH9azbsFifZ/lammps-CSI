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

#include "mpi.h"
#include "string.h"
#include "stdlib.h"
#include "fix_ave_force.h"
#include "atom.h"
#include "update.h"
#include "respa.h"
#include "error.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

FixAveForce::FixAveForce(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  if (narg != 6) error->all("Illegal fix aveforce command");

  xflag = yflag = zflag = 1;
  if (strcmp(arg[3],"NULL") == 0) xflag = 0;
  else xvalue = atof(arg[3]);
  if (strcmp(arg[4],"NULL") == 0) yflag = 0;
  else yvalue = atof(arg[4]);
  if (strcmp(arg[5],"NULL") == 0) zflag = 0;
  else zvalue = atof(arg[5]);
}

/* ---------------------------------------------------------------------- */

int FixAveForce::setmask()
{
  int mask = 0;
  mask |= POST_FORCE;
  mask |= POST_FORCE_RESPA;
  mask |= MIN_POST_FORCE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixAveForce::init()
{
  if (strcmp(update->integrate_style,"respa") == 0)
    nlevels_respa = ((Respa *) update->integrate)->nlevels;

  // ncount = total # of atoms in group

  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  int count = 0;
  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) count++;
  MPI_Allreduce(&count,&ncount,1,MPI_INT,MPI_SUM,world);
}

/* ---------------------------------------------------------------------- */

void FixAveForce::setup()
{
  if (strcmp(update->integrate_style,"verlet") == 0)
    post_force(1);
  else
    for (int ilevel = 0; ilevel < nlevels_respa; ilevel++) {
      ((Respa *) update->integrate)->copy_flevel_f(ilevel);
      post_force_respa(1,ilevel,0);
      ((Respa *) update->integrate)->copy_f_flevel(ilevel);
    }
}

/* ---------------------------------------------------------------------- */

void FixAveForce::min_setup()
{
  post_force(1);
}

/* ---------------------------------------------------------------------- */

void FixAveForce::post_force(int vflag)
{
  // sum forces on participating atoms

  double **f = atom->f;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  double sum[3];
  sum[0] = sum[1] = sum[2] = 0.0;

  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) {
      sum[0] += f[i][0];
      sum[1] += f[i][1];
      sum[2] += f[i][2];
    }

  // average the force on participating atoms
  // add in requested amount

  double sumall[3];
  MPI_Allreduce(sum,sumall,3,MPI_DOUBLE,MPI_SUM,world);
  sumall[0] = sumall[0]/ncount + xvalue;
  sumall[1] = sumall[1]/ncount + yvalue;
  sumall[2] = sumall[2]/ncount + zvalue;

  // set force of all participating atoms to same value
  // only for active dimensions

  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) {
      if (xflag) f[i][0] = sumall[0];
      if (yflag) f[i][1] = sumall[1];
      if (zflag) f[i][2] = sumall[2];
    }
}

/* ---------------------------------------------------------------------- */

void FixAveForce::post_force_respa(int vflag, int ilevel, int iloop)
{
  // ave + extra force on outermost level
  // just ave on inner levels

  if (ilevel == nlevels_respa-1) post_force(vflag);
  else {
    double **f = atom->f;
    int *mask = atom->mask;
    int nlocal = atom->nlocal;

    double sum[3];
    sum[0] = sum[1] = sum[2] = 0.0;

    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) {
	sum[0] += f[i][0];
	sum[1] += f[i][1];
	sum[2] += f[i][2];
      }

    double sumall[3];
    MPI_Allreduce(sum,sumall,3,MPI_DOUBLE,MPI_SUM,world);
    sumall[0] = sumall[0]/ncount;
    sumall[1] = sumall[1]/ncount;
    sumall[2] = sumall[2]/ncount;

    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) {
	if (xflag) f[i][0] = sumall[0];
	if (yflag) f[i][1] = sumall[1];
	if (zflag) f[i][2] = sumall[2];
      }
  }
}

/* ---------------------------------------------------------------------- */

void FixAveForce::min_post_force(int vflag)
{
  post_force(vflag);
}
