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

#include "neighbor.h"
#include "atom.h"
#include "error.h"

using namespace LAMMPS_NS;

/* ----------------------------------------------------------------------
   N^2 search for all neighbors
   every neighbor pair appears in list of both atoms i and j
------------------------------------------------------------------------- */

void Neighbor::full_nsq()
{
  int i,j,n,itype,jtype,which;
  double xtmp,ytmp,ztmp,delx,dely,delz,rsq;
  int *neighptr;

  double **x = atom->x;
  int *type = atom->type;
  int *mask = atom->mask;
  int *molecule = atom->molecule;
  int nlocal = atom->nlocal;
  int nall = atom->nlocal + atom->nghost;
  int molecular = atom->molecular;

  int npage = 0;
  int npnt = 0;

  for (i = 0; i < nlocal; i++) {

    if (pgsize - npnt < oneatom) {
      npnt = 0;
      npage++;
      if (npage == maxpage_full) add_pages_full(npage);
    }

    neighptr = &pages_full[npage][npnt];
    n = 0;

    itype = type[i];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];

    // loop over all atoms, owned and ghost
    // skip i = j

    for (j = 0; j < nall; j++) {
      if (i == j) continue;
      if (exclude && exclusion(i,j,type,mask,molecule)) continue;

      jtype = type[j];
      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx*delx + dely*dely + delz*delz;
      if (rsq <= cutneighsq[itype][jtype]) {
	if (molecular) which = find_special(i,j);
	else which = 0;
	if (which == 0) neighptr[n++] = j;
	else if (which > 0) neighptr[n++] = which*nall + j;
      }
    }

    firstneigh_full[i] = neighptr;
    numneigh_full[i] = n;
    npnt += n;
    if (npnt >= pgsize)
      error->one("Neighbor list overflow, boost neigh_modify one or page");
  }
}

/* ----------------------------------------------------------------------
   binned neighbor list construction for all neighbors
   every neighbor pair appears in list of both atoms i and j
------------------------------------------------------------------------- */

void Neighbor::full_bin()
{
  int i,j,k,n,itype,jtype,ibin,which;
  double xtmp,ytmp,ztmp,delx,dely,delz,rsq;
  int *neighptr;

  // bin local & ghost atoms

  bin_atoms();

  // loop over each atom, storing neighbors

  double **x = atom->x;
  int *type = atom->type;
  int *mask = atom->mask;
  int *molecule = atom->molecule;
  int nlocal = atom->nlocal;
  int nall = atom->nlocal + atom->nghost;
  int molecular = atom->molecular;

  int npage = 0;
  int npnt = 0;

  for (i = 0; i < nlocal; i++) {

    if (pgsize - npnt < oneatom) {
      npnt = 0;
      npage++;
      if (npage == maxpage_full) add_pages_full(npage);
    }

    neighptr = &pages_full[npage][npnt];
    n = 0;

    itype = type[i];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];

    // loop over all atoms in surrounding bins in stencil including self
    // skip i = j

    ibin = coord2bin(x[i]);

    for (k = 0; k < nstencil_full; k++) {
      for (j = binhead[ibin+stencil_full[k]]; j >= 0; j = bins[j]) {
	if (i == j) continue;
	if (exclude && exclusion(i,j,type,mask,molecule)) continue;

	jtype = type[j];
	delx = xtmp - x[j][0];
	dely = ytmp - x[j][1];
	delz = ztmp - x[j][2];
	rsq = delx*delx + dely*dely + delz*delz;

	if (rsq <= cutneighsq[itype][jtype]) {
	  if (molecular) which = find_special(i,j);
	  else which = 0;
	  if (which == 0) neighptr[n++] = j;
	  else if (which > 0) neighptr[n++] = which*nall + j;
	}
      }
    }

    firstneigh_full[i] = neighptr;
    numneigh_full[i] = n;
    npnt += n;
    if (npnt >= pgsize)
      error->one("Neighbor list overflow, boost neigh_modify one or page");
  }
}

/* ----------------------------------------------------------------------
   binned neighbor list construction for all neighbors
   multi-type stencil is itype dependent and is distance checked
   every neighbor pair appears in list of both atoms i and j
------------------------------------------------------------------------- */

void Neighbor::full_bin_multi()
{
  int i,j,k,n,itype,jtype,ibin,which,ns;
  double xtmp,ytmp,ztmp,delx,dely,delz,rsq;
  int *neighptr,*s;
  double *cutsq,*distsq;

  // bin local & ghost atoms

  bin_atoms();

  // loop over each atom, storing neighbors

  double **x = atom->x;
  int *type = atom->type;
  int *mask = atom->mask;
  int *molecule = atom->molecule;
  int nlocal = atom->nlocal;
  int nall = atom->nlocal + atom->nghost;
  int molecular = atom->molecular;

  int npage = 0;
  int npnt = 0;

  for (i = 0; i < nlocal; i++) {

    if (pgsize - npnt < oneatom) {
      npnt = 0;
      npage++;
      if (npage == maxpage_full) add_pages_full(npage);
    }

    neighptr = &pages_full[npage][npnt];
    n = 0;

    itype = type[i];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];

    // loop over all atoms in other bins in stencil, including self
    // skip if i,j neighbor cutoff is less than bin distance
    // skip i = j

    ibin = coord2bin(x[i]);
    s = stencil_full_multi[itype];
    distsq = distsq_full_multi[itype];
    cutsq = cutneighsq[itype];
    ns = nstencil_full_multi[itype];
    for (k = 0; k < ns; k++) {
      for (j = binhead[ibin+s[k]]; j >= 0; j = bins[j]) {
	jtype = type[j];
	if (cutsq[jtype] < distsq[k]) continue;
	if (i == j) continue;
	if (exclude && exclusion(i,j,type,mask,molecule)) continue;

	delx = xtmp - x[j][0];
	dely = ytmp - x[j][1];
	delz = ztmp - x[j][2];
	rsq = delx*delx + dely*dely + delz*delz;

	if (rsq <= cutneighsq[itype][jtype]) {
	  if (molecular) which = find_special(i,j);
	  else which = 0;
	  if (which == 0) neighptr[n++] = j;
	  else if (which > 0) neighptr[n++] = which*nall + j;
	}
      }
    }

    firstneigh_full[i] = neighptr;
    numneigh_full[i] = n;
    npnt += n;
    if (npnt >= pgsize)
      error->one("Neighbor list overflow, boost neigh_modify one or page");
  }
}
