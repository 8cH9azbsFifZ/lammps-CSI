/* -----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under 
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------ */

#include "math.h"
#include "string.h"
#include "bond_hybrid.h"
#include "atom.h"
#include "neighbor.h"
#include "domain.h"
#include "comm.h"
#include "force.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;

#define EXTRA 1000

/* ---------------------------------------------------------------------- */

BondHybrid::BondHybrid(LAMMPS *lmp) : Bond(lmp)
{
  nstyles = 0;
}

/* ---------------------------------------------------------------------- */

BondHybrid::~BondHybrid()
{
  if (nstyles) {
    for (int i = 0; i < nstyles; i++) delete styles[i];
    delete [] styles;
    for (int i = 0; i < nstyles; i++) delete [] keywords[i];
    delete [] keywords;
  }

  if (allocated) {
    memory->sfree(setflag);
    memory->sfree(map);
    delete [] nbondlist;
    delete [] maxbond;
    for (int i = 0; i < nstyles; i++)
      memory->destroy_2d_int_array(bondlist[i]);
    delete [] bondlist;
  }
}

/* ---------------------------------------------------------------------- */

void BondHybrid::compute(int eflag, int vflag)
{
  int i,m,n;

  // save ptrs to original bondlist

  int nbondlist_orig = neighbor->nbondlist;
  int **bondlist_orig = neighbor->bondlist;

  // if this is re-neighbor step, create sub-style bondlists
  // nbondlist[] = length of each sub-style list
  // realloc sub-style bondlist if necessary
  // load sub-style bondlist with 3 values from original bondlist

  if (neighbor->ago == 0) {
    for (m = 0; m < nstyles; m++) nbondlist[m] = 0;
    for (i = 0; i < nbondlist_orig; i++)
      nbondlist[map[bondlist_orig[i][2]]]++;
    for (m = 0; m < nstyles; m++) {
      if (nbondlist[m] > maxbond[m]) {
	memory->destroy_2d_int_array(bondlist[m]);
	maxbond[m] = nbondlist[m] + EXTRA;
	bondlist[m] = (int **)
	  memory->create_2d_int_array(maxbond[m],3,"bond_hybrid:bondlist");
      }
      nbondlist[m] = 0;
    }
    for (i = 0; i < nbondlist_orig; i++) {
      m = map[bondlist_orig[i][2]];
      n = nbondlist[m];
      bondlist[m][n][0] = bondlist_orig[i][0];
      bondlist[m][n][1] = bondlist_orig[i][1];
      bondlist[m][n][2] = bondlist_orig[i][2];
      nbondlist[m]++;
    }
  }
  
  // call each sub-style's compute function
  // must set neighbor->bondlist to sub-style bondlist before call
  // accumulate sub-style energy,virial in hybrid's energy,virial

  energy = 0.0;
  eng_vdwl = 0.0;
  if (vflag) for (n = 0; n < 6; n++) virial[n] = 0.0;

  for (m = 0; m < nstyles; m++) {
    if (styles[m] == NULL) continue;
    neighbor->nbondlist = nbondlist[m];
    neighbor->bondlist = bondlist[m];
    styles[m]->compute(eflag,vflag);
    if (eflag) {
      energy += styles[m]->energy;
      eng_vdwl += styles[m]->eng_vdwl;
    }
    if (vflag) for (n = 0; n < 6; n++) virial[n] += styles[m]->virial[n];
  }

  // restore ptrs to original bondlist

  neighbor->nbondlist = nbondlist_orig;
  neighbor->bondlist = bondlist_orig;
}

/* ---------------------------------------------------------------------- */

void BondHybrid::allocate()
{
  allocated = 1;
  int n = atom->nbondtypes;

  map = (int *) memory->smalloc((n+1)*sizeof(int),"bond:map");
  setflag = (int *) memory->smalloc((n+1)*sizeof(int),"bond:setflag");
  for (int i = 1; i <= n; i++) setflag[i] = 0;

  nbondlist = new int[nstyles];
  maxbond = new int[nstyles];
  bondlist = new int**[nstyles];
  for (int m = 0; m < nstyles; m++) maxbond[m] = 0;
  for (int m = 0; m < nstyles; m++) bondlist[m] = NULL;
}

/* ----------------------------------------------------------------------
   create one bond style for each arg in list
------------------------------------------------------------------------- */

void BondHybrid::settings(int narg, char **arg)
{
  nstyles = narg;
  styles = new Bond*[nstyles];
  keywords = new char*[nstyles];

  for (int m = 0; m < nstyles; m++) {
    for (int i = 0; i < m; i++)
      if (strcmp(arg[m],arg[i]) == 0) 
	error->all("Bond style hybrid cannot use same bond style twice");
    if (strcmp(arg[m],"hybrid") == 0) 
      error->all("Bond style hybrid cannot have hybrid as an argument");
    styles[m] = force->new_bond(arg[m]);
    keywords[m] = new char[strlen(arg[m])+1];
    strcpy(keywords[m],arg[m]);
  }
}

/* ----------------------------------------------------------------------
   set coeffs for one type
---------------------------------------------------------------------- */

void BondHybrid::coeff(int narg, char **arg)
{
  if (!allocated) allocate();

  int ilo,ihi;
  force->bounds(arg[0],atom->nbondtypes,ilo,ihi);

  // 2nd arg = bond style name (harmonic, fene, etc)

  int m;
  for (m = 0; m < nstyles; m++)
    if (strcmp(arg[1],keywords[m]) == 0) break;
  if (m == nstyles) error->all("Bond coeff for hybrid has invalid style");

  // set low-level coefficients for each bondtype
  // replace 2nd arg with i, call coeff() with no 1st arg
  // if sub-style is NULL for "none", still set setflag

  for (int i = ilo; i <= ihi; i++) {
    sprintf(arg[1],"%d",i);
    map[i] = m;
    if (styles[m]) styles[m]->coeff(narg-1,&arg[1]);
    setflag[i] = 1;
  }
}

/* ---------------------------------------------------------------------- */

void BondHybrid::init_style()
{
  for (int m = 0; m < nstyles; m++) 
    if (styles[m]) styles[m]->init_style();
}

/* ----------------------------------------------------------------------
   return an equilbrium bond length 
------------------------------------------------------------------------- */

double BondHybrid::equilibrium_distance(int i)
{
  return styles[map[i]]->equilibrium_distance(i);
}

/* ----------------------------------------------------------------------
   proc 0 writes to restart file
------------------------------------------------------------------------- */

void BondHybrid::write_restart(FILE *fp)
{
  fwrite(&nstyles,sizeof(int),1,fp);

  int n;
  for (int m = 0; m < nstyles; m++) {
    n = strlen(keywords[m]) + 1;
    fwrite(&n,sizeof(int),1,fp);
    fwrite(keywords[m],sizeof(char),n,fp);
  }
}

/* ----------------------------------------------------------------------
   proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void BondHybrid::read_restart(FILE *fp)
{
  allocate();

  int me = comm->me;
  if (me == 0) fread(&nstyles,sizeof(int),1,fp);
  MPI_Bcast(&nstyles,1,MPI_INT,0,world);
  styles = new Bond*[nstyles];
  keywords = new char*[nstyles];
  
  int n;
  for (int m = 0; m < nstyles; m++) {
    if (me == 0) fread(&n,sizeof(int),1,fp);
    MPI_Bcast(&n,1,MPI_INT,0,world);
    keywords[m] = new char[n];
    if (me == 0) fread(keywords[m],sizeof(char),n,fp);
    MPI_Bcast(keywords[m],n,MPI_CHAR,0,world);
    styles[m] = force->new_bond(keywords[m]);
  }
}

/* ---------------------------------------------------------------------- */

void BondHybrid::single(int type, double rsq, int i, int j, double rfactor,
			int eflag, double &fforce, double &eng)
{
  if (styles[map[type]]) 
    styles[map[type]]->single(type,rsq,i,j,rfactor,eflag,fforce,eng);
}

/* ----------------------------------------------------------------------
   memory usage
------------------------------------------------------------------------- */

int BondHybrid::memory_usage()
{
  int bytes = 0;
  for (int m = 0; m < nstyles; m++) bytes += maxbond[m]*3 * sizeof(int);
  for (int m = 0; m < nstyles; m++)
    if (styles[m]) bytes += styles[m]->memory_usage();
  return bytes;
}
