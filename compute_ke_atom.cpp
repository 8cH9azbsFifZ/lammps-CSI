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

#include "string.h"
#include "compute_ke_atom.h"
#include "atom.h"
#include "modify.h"
#include "comm.h"
#include "force.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

ComputeKEAtom::ComputeKEAtom(LAMMPS *lmp, int narg, char **arg) :
  Compute(lmp, narg, arg)
{
  if (narg != 3) error->all("Illegal compute ke/atom command");

  peratom_flag = 1;
  size_peratom = 0;

  nmax = 0;
  ke = NULL;
}

/* ---------------------------------------------------------------------- */

ComputeKEAtom::~ComputeKEAtom()
{
  memory->sfree(ke);
}

/* ---------------------------------------------------------------------- */

void ComputeKEAtom::init()
{
  int count = 0;
  for (int i = 0; i < modify->ncompute; i++)
    if (strcmp(modify->compute[i]->style,"ke/atom") == 0) count++;
  if (count > 1 && comm->me == 0)
    error->warning("More than one compute ke/atom");
}

/* ---------------------------------------------------------------------- */

void ComputeKEAtom::compute_peratom()
{
  // grow ke array if necessary

  if (atom->nlocal > nmax) {
    memory->sfree(ke);
    nmax = atom->nmax;
    ke = (double *) memory->smalloc(nmax*sizeof(double),"compute/ke/atom:ke");
    scalar_atom = ke;
  }

  // compute kinetic energy for each atom in group

  double mvv2e = force->mvv2e;
  double **v = atom->v;
  double *mass = atom->mass;
  double *rmass = atom->rmass;
  int *mask = atom->mask;
  int *type = atom->type;
  int nlocal = atom->nlocal;

  if (mass)
    for (int i = 0; i < nlocal; i++) {
      if (mask[i] & groupbit) {
	ke[i] = 0.5 * mvv2e * mass[type[i]] *
	  (v[i][0]*v[i][0] + v[i][1]*v[i][1] + v[i][2]*v[i][2]);
      }
    }
  else
    for (int i = 0; i < nlocal; i++) {
      if (mask[i] & groupbit) {
	ke[i] = 0.5 * mvv2e * rmass[i] *
	  (v[i][0]*v[i][0] + v[i][1]*v[i][1] + v[i][2]*v[i][2]);
      }
    }
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based array
------------------------------------------------------------------------- */

int ComputeKEAtom::memory_usage()
{
  int bytes = nmax * sizeof(double);
  return bytes;
}