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
#include "string.h"
#include "fix_wall_lj93.h"
#include "atom.h"
#include "domain.h"
#include "update.h"
#include "output.h"
#include "respa.h"
#include "error.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

FixWallLJ93::FixWallLJ93(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  if (narg != 8) error->all("Illegal fix wall/lj93 command");

  if (strcmp(arg[3],"xlo") == 0) {
    dim = 0;
    side = -1;
  } else if (strcmp(arg[3],"xhi") == 0) {
    dim = 0;
    side = 1;
  } else if (strcmp(arg[3],"ylo") == 0) {
    dim = 1;
    side = -1;
  } else if (strcmp(arg[3],"yhi") == 0) {
    dim = 1;
    side = 1;
  } else if (strcmp(arg[3],"zlo") == 0) {
    dim = 2;
    side = -1;
  } else if (strcmp(arg[3],"zhi") == 0) {
    dim = 2;
    side = 1;
  } else error->all("Illegal fix wall/lj93 command");

  coord = atof(arg[4]);
  epsilon = atof(arg[5]);
  sigma = atof(arg[6]);
  cutoff = atof(arg[7]);

  coeff1 = 6.0/5.0 * epsilon * pow(sigma,9.0);
  coeff2 = 3.0 * epsilon * pow(sigma,3.0);
  coeff3 = 2.0/15.0 * epsilon * pow(sigma,9.0);
  coeff4 = epsilon * pow(sigma,3.0);

  double rinv = 1.0/cutoff;
  double r2inv = rinv*rinv;
  double r4inv = r2inv*r2inv;
  offset = coeff3*r4inv*r4inv*rinv - coeff4*r2inv*rinv;

  if (dim == 0 && domain->xperiodic)
    error->all("Cannot use wall in periodic dimension");
  if (dim == 1 && domain->yperiodic)
    error->all("Cannot use wall in periodic dimension");
  if (dim == 2 && domain->zperiodic)
    error->all("Cannot use wall in periodic dimension");
}

/* ---------------------------------------------------------------------- */

int FixWallLJ93::setmask()
{
  int mask = 0;
  mask |= POST_FORCE;
  mask |= THERMO_ENERGY;
  mask |= POST_FORCE_RESPA;
  mask |= MIN_POST_FORCE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixWallLJ93::init()
{
  if (strcmp(update->integrate_style,"respa") == 0)
    nlevels_respa = ((Respa *) update->integrate)->nlevels;
}

/* ---------------------------------------------------------------------- */

void FixWallLJ93::setup()
{
  eflag_enable = 1;
  if (strcmp(update->integrate_style,"verlet") == 0)
    post_force(1);
  else {
    ((Respa *) update->integrate)->copy_flevel_f(nlevels_respa-1);
    post_force_respa(1,nlevels_respa-1,0);
    ((Respa *) update->integrate)->copy_f_flevel(nlevels_respa-1);
  }
  eflag_enable = 0;
}

/* ---------------------------------------------------------------------- */

void FixWallLJ93::min_setup()
{
  eflag_enable = 1;
  post_force(1);
}

/* ---------------------------------------------------------------------- */

void FixWallLJ93::post_force(int vflag)
{
  bool eflag = false;
  if (eflag_enable) eflag = true;
  else if (output->next_thermo == update->ntimestep) eflag = true;

  double **x = atom->x;
  double **f = atom->f;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  double delta,rinv,r2inv,r4inv,r10inv,eng;
  if (eflag) eng = 0.0;

  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) {
      if (side == -1) delta = x[i][dim] - coord;
      else delta = coord - x[i][dim];
      if (delta <= 0.0) continue;
      if (delta > cutoff) continue;
      rinv = 1.0/delta;
      r2inv = rinv*rinv;
      r4inv = r2inv*r2inv;
      r10inv = r4inv*r4inv*r2inv;
      f[i][dim] -= (coeff1*r10inv - coeff2*r4inv) * side;
      if (eflag) eng += coeff3*r4inv*r4inv*rinv - coeff4*r2inv*rinv - offset;
    }

  if (eflag) MPI_Allreduce(&eng,&etotal,1,MPI_DOUBLE,MPI_SUM,world);
}

/* ---------------------------------------------------------------------- */

void FixWallLJ93::post_force_respa(int vflag, int ilevel, int iloop)
{
  if (ilevel == nlevels_respa-1) post_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixWallLJ93::min_post_force(int vflag)
{
  post_force(vflag);
}

/* ---------------------------------------------------------------------- */

double FixWallLJ93::thermo(int n)
{
  if (n == 0) return etotal;
  else return 0.0;
}
