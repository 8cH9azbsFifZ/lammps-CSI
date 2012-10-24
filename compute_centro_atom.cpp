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
#include "compute_centro_atom.h"
#include "atom.h"
#include "modify.h"
#include "update.h"
#include "neighbor.h"
#include "force.h"
#include "pair.h"
#include "comm.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

ComputeCentroAtom::ComputeCentroAtom(LAMMPS *lmp, int narg, char **arg) :
  Compute(lmp, narg, arg)
{
  if (narg != 3) error->all("Illegal compute centro/atom command");

  peratom_flag = 1;
  size_peratom = 0;
  neigh_full_once = 1;

  nmax = 0;
  centro = NULL;
  maxneigh = 0;
  distsq = NULL;
  nearest = NULL;
}

/* ---------------------------------------------------------------------- */

ComputeCentroAtom::~ComputeCentroAtom()
{
  memory->sfree(centro);
  memory->sfree(distsq);
  memory->sfree(nearest);
}

/* ---------------------------------------------------------------------- */

void ComputeCentroAtom::init()
{
  int count = 0;
  for (int i = 0; i < modify->ncompute; i++)
    if (strcmp(modify->compute[i]->style,"centro/atom") == 0) count++;
  if (count > 1 && comm->me == 0)
    error->warning("More than one compute centro/atom");
}

/* ---------------------------------------------------------------------- */

void ComputeCentroAtom::compute_peratom()
{
  int j,k,jj,kk,n,numneigh;
  double xtmp,ytmp,ztmp,delx,dely,delz,rsq,value;
  int *neighs;
  double pairs[66];

  double nIDEAL = 0.,
         nPARTIALDISLOCATION = 0.,
         nSTACKINGFAULT = 0.,
         nSURFACE = 0.;

  // grow centro array if necessary

  if (atom->nlocal > nmax) {
    memory->sfree(centro);
    nmax = atom->nmax;
    centro = (double *) 
      memory->smalloc(nmax*sizeof(double),"compute/centro/atom:centro");
    scalar_atom = centro;
  }

  // if needed, build a full neighbor list

  if (!neighbor->full_every) neighbor->build_full();

  // compute centro-symmetry parameter for each atom in group
  // use full neighbor list

  double **x = atom->x;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  int nall = atom->nlocal + atom->nghost;
  double cutsq = force->pair->cutforce * force->pair->cutforce;

  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) {
      xtmp = x[i][0];
      ytmp = x[i][1];
      ztmp = x[i][2];
      neighs = neighbor->firstneigh_full[i];
      numneigh = neighbor->numneigh_full[i];

      // insure distsq and nearest arrays are long enough

      if (numneigh > maxneigh) {
	memory->sfree(distsq);
	memory->sfree(nearest);
	maxneigh = numneigh;
	distsq = (double *) memory->smalloc(maxneigh*sizeof(double),
					    "compute/centro/atom:distsq");
	nearest = (int *) memory->smalloc(maxneigh*sizeof(int),
					  "compute/centro/atom:nearest");
      }

      // loop over list of all neighbors within force cutoff
      // distsq[] = distance sq to each
      // nearest[] = atom indices of neighbors

      n = 0;
      for (k = 0; k < numneigh; k++) {
	j = neighs[k];
	if (j >= nall) j %= nall;

	delx = xtmp - x[j][0];
	dely = ytmp - x[j][1];
	delz = ztmp - x[j][2];
	rsq = delx*delx + dely*dely + delz*delz;
	if (rsq < cutsq) {
	  distsq[n] = rsq;
	  nearest[n++] = j;
	}
      }

      // if not 12 neighbors, centro = 0.0

      if (n < 12) {
	centro[i] = 0.0;
	continue;
      }

      // store 12 nearest neighs in 1st 12 locations of distsq and nearest

      select2(12,n,distsq,nearest);

      // R = Ri + Rj for each of 66 i,j pairs among 12 neighbors
      // pairs = squared length of each R

      n = 0;
      for (j = 0; j < 12; j++) {
	jj = nearest[j];
	for (k = j+1; k < 12; k++) {
	  kk = nearest[k];
	  delx = x[jj][0] + x[kk][0] - 2.0*xtmp;
	  dely = x[jj][1] + x[kk][1] - 2.0*ytmp;
	  delz = x[jj][2] + x[kk][2] - 2.0*ztmp;
	  pairs[n++] = delx*delx + dely*dely + delz*delz;
	}
      }

      // store 6 smallest pair distances in 1st 6 locations of pairs

      select(6,66,pairs);

      // centrosymmetry = sum of 6 smallest squared values

      value = 0.0;
      for (j = 0; j < 6; j++) value += pairs[j];
      centro[i] = value;
    }
}

/* ----------------------------------------------------------------------
   2 select routines from Numerical Recipes (slightly modified)
   find k smallest values in array of length n
   2nd routine sorts auxiliary array at same time
------------------------------------------------------------------------- */

#define SWAP(a,b)   tmp = a; a = b; b = tmp;
#define ISWAP(a,b) itmp = a; a = b; b = itmp;

void ComputeCentroAtom::select(int k, int n, double *arr)
{
  int i,ir,j,l,mid;
  double a,tmp;

  arr--;
  l = 1;
  ir = n;
  for (;;) {
    if (ir <= l+1) {
      if (ir == l+1 && arr[ir] < arr[l]) {
	SWAP(arr[l],arr[ir])
      }
      return;
    } else {
      mid=(l+ir) >> 1;
      SWAP(arr[mid],arr[l+1])
      if (arr[l] > arr[ir]) {
	SWAP(arr[l],arr[ir])
      }
      if (arr[l+1] > arr[ir]) {
	SWAP(arr[l+1],arr[ir])
      }
      if (arr[l] > arr[l+1]) {
	SWAP(arr[l],arr[l+1])
      }
      i = l+1;
      j = ir;
      a = arr[l+1];
      for (;;) {
	do i++; while (arr[i] < a);
	do j--; while (arr[j] > a);
	if (j < i) break;
	SWAP(arr[i],arr[j])
      }
      arr[l+1] = arr[j];
      arr[j] = a;
      if (j >= k) ir = j-1;
      if (j <= k) l = i;
    }
  }
}

/* ---------------------------------------------------------------------- */

void ComputeCentroAtom::select2(int k, int n, double *arr, int *iarr)
{
  int i,ir,j,l,mid,ia,itmp;
  double a,tmp;

  arr--;
  iarr--;
  l = 1;
  ir = n;
  for (;;) {
    if (ir <= l+1) {
      if (ir == l+1 && arr[ir] < arr[l]) {
	SWAP(arr[l],arr[ir])
	ISWAP(iarr[l],iarr[ir])
      }
      return;
    } else {
      mid=(l+ir) >> 1;
      SWAP(arr[mid],arr[l+1])
      ISWAP(iarr[mid],iarr[l+1])
      if (arr[l] > arr[ir]) {
	SWAP(arr[l],arr[ir])
	ISWAP(iarr[l],iarr[ir])
      }
      if (arr[l+1] > arr[ir]) {
	SWAP(arr[l+1],arr[ir])
	ISWAP(iarr[l+1],iarr[ir])
      }
      if (arr[l] > arr[l+1]) {
	SWAP(arr[l],arr[l+1])
	ISWAP(iarr[l],iarr[l+1])
      }
      i = l+1;
      j = ir;
      a = arr[l+1];
      ia = iarr[l+1];
      for (;;) {
	do i++; while (arr[i] < a);
	do j--; while (arr[j] > a);
	if (j < i) break;
	SWAP(arr[i],arr[j])
	ISWAP(iarr[i],iarr[j])
      }
      arr[l+1] = arr[j];
      arr[j] = a;
      iarr[l+1] = iarr[j];
      iarr[j] = ia;
      if (j >= k) ir = j-1;
      if (j <= k) l = i;
    }
  }
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based array
------------------------------------------------------------------------- */

int ComputeCentroAtom::memory_usage()
{
  int bytes = nmax * sizeof(double);
  return bytes;
}
