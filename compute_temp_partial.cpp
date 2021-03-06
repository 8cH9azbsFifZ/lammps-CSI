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
#include "stdlib.h"
#include "compute_temp_partial.h"
#include "atom.h"
#include "force.h"
#include "modify.h"
#include "fix.h"
#include "group.h"
#include "error.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

ComputeTempPartial::ComputeTempPartial(LAMMPS *lmp, int narg, char **arg) : 
  Compute(lmp, narg, arg)
{
  if (narg != 6) error->all("Illegal compute temp/partial command");

  xflag = atoi(arg[3]);
  yflag = atoi(arg[4]);
  zflag = atoi(arg[5]);

  scalar_flag = vector_flag = 1;
  size_vector = 6;
  extensive = 0;
  tempflag = 1;

  vector = new double[6];
}

/* ---------------------------------------------------------------------- */

ComputeTempPartial::~ComputeTempPartial()
{
  delete [] vector;
}

/* ---------------------------------------------------------------------- */

void ComputeTempPartial::init()
{
  fix_dof = 0;
  for (int i = 0; i < modify->nfix; i++)
    fix_dof += modify->fix[i]->dof(igroup);
  recount();
}

/* ---------------------------------------------------------------------- */

void ComputeTempPartial::recount()
{
  double natoms = group->count(igroup);
  dof = (xflag+yflag+zflag) * natoms;
  dof -= extra_dof + fix_dof;
  if (dof > 0) tfactor = force->mvv2e / (dof * force->boltz);
  else tfactor = 0.0;
}

/* ---------------------------------------------------------------------- */

double ComputeTempPartial::compute_scalar()
{
  double **v = atom->v;
  double *mass = atom->mass;
  double *rmass = atom->rmass;
  int *type = atom->type;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  double t = 0.0;

  if (mass) {
    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit)
	t += (xflag*v[i][0]*v[i][0] + yflag*v[i][1]*v[i][1] + 
	      zflag*v[i][2]*v[i][2]) * mass[type[i]];
  } else {
    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit)
	t += (xflag*v[i][0]*v[i][0] + yflag*v[i][1]*v[i][1] + 
	      zflag*v[i][2]*v[i][2]) * rmass[i];
  }
    
  MPI_Allreduce(&t,&scalar,1,MPI_DOUBLE,MPI_SUM,world);
  if (dynamic) recount();
  scalar *= tfactor;
  return scalar;
}

/* ---------------------------------------------------------------------- */

void ComputeTempPartial::compute_vector()
{
  int i;

  double **v = atom->v;
  double *mass = atom->mass;
  double *rmass = atom->rmass;
  int *type = atom->type;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  double massone,t[6];
  for (i = 0; i < 6; i++) t[i] = 0.0;

  for (i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) {
      if (mass) massone = mass[type[i]];
      else massone = rmass[i];
      t[0] += massone * xflag*v[i][0]*v[i][0];
      t[1] += massone * yflag*v[i][1]*v[i][1];
      t[2] += massone * zflag*v[i][2]*v[i][2];
      t[3] += massone * xflag*yflag*v[i][0]*v[i][1];
      t[4] += massone * xflag*zflag*v[i][0]*v[i][2];
      t[5] += massone * yflag*zflag*v[i][1]*v[i][2];
    }

  MPI_Allreduce(t,vector,6,MPI_DOUBLE,MPI_SUM,world);
  for (i = 0; i < 6; i++) vector[i] *= force->mvv2e;
}
