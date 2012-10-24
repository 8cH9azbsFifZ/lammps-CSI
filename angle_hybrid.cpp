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
#include "angle_hybrid.h"
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

AngleHybrid::AngleHybrid(LAMMPS *lmp) : Angle(lmp)
{
  nstyles = 0;
}

/* ---------------------------------------------------------------------- */

AngleHybrid::~AngleHybrid()
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
    delete [] nanglelist;
    delete [] maxangle;
    for (int i = 0; i < nstyles; i++)
      memory->destroy_2d_int_array(anglelist[i]);
    delete [] anglelist;
  }
}

/* ---------------------------------------------------------------------- */

void AngleHybrid::compute(int eflag, int vflag)
{
  int i,m,n;

  // save ptrs to original anglelist

  int nanglelist_orig = neighbor->nanglelist;
  int **anglelist_orig = neighbor->anglelist;

  // if this is re-neighbor step, create sub-style anglelists
  // nanglelist[] = length of each sub-style list
  // realloc sub-style anglelist if necessary
  // load sub-style anglelist with 4 values from original anglelist

  if (neighbor->ago == 0) {
    for (m = 0; m < nstyles; m++) nanglelist[m] = 0;
    for (i = 0; i < nanglelist_orig; i++)
      nanglelist[map[anglelist_orig[i][3]]]++;
    for (m = 0; m < nstyles; m++) {
      if (nanglelist[m] > maxangle[m]) {
	memory->destroy_2d_int_array(anglelist[m]);
	maxangle[m] = nanglelist[m] + EXTRA;
	anglelist[m] = (int **)
	  memory->create_2d_int_array(maxangle[m],4,"angle_hybrid:anglelist");
      }
      nanglelist[m] = 0;
    }
    for (i = 0; i < nanglelist_orig; i++) {
      m = map[anglelist_orig[i][3]];
      n = nanglelist[m];
      anglelist[m][n][0] = anglelist_orig[i][0];
      anglelist[m][n][1] = anglelist_orig[i][1];
      anglelist[m][n][2] = anglelist_orig[i][2];
      anglelist[m][n][3] = anglelist_orig[i][3];
      nanglelist[m]++;
    }
  }
  
  // call each sub-style's compute function
  // must set neighbor->anglelist to sub-style anglelist before call
  // accumulate sub-style energy,virial in hybrid's energy,virial

  energy = 0.0;
  if (vflag) for (n = 0; n < 6; n++) virial[n] = 0.0;

  for (m = 0; m < nstyles; m++) {
    if (styles[m] == NULL) continue;
    neighbor->nanglelist = nanglelist[m];
    neighbor->anglelist = anglelist[m];
    styles[m]->compute(eflag,vflag);
    if (eflag) energy += styles[m]->energy;
    if (vflag) for (n = 0; n < 6; n++) virial[n] += styles[m]->virial[n];
  }

  // restore ptrs to original anglelist

  neighbor->nanglelist = nanglelist_orig;
  neighbor->anglelist = anglelist_orig;
}

/* ---------------------------------------------------------------------- */

void AngleHybrid::allocate()
{
  allocated = 1;
  int n = atom->nangletypes;

  map = (int *) memory->smalloc((n+1)*sizeof(int),"angle:map");
  setflag = (int *) memory->smalloc((n+1)*sizeof(int),"angle:setflag");
  for (int i = 1; i <= n; i++) setflag[i] = 0;

  nanglelist = new int[nstyles];
  maxangle = new int[nstyles];
  anglelist = new int**[nstyles];
  for (int m = 0; m < nstyles; m++) maxangle[m] = 0;
  for (int m = 0; m < nstyles; m++) anglelist[m] = NULL;
}

/* ----------------------------------------------------------------------
   create one angle style for each arg in list
------------------------------------------------------------------------- */

void AngleHybrid::settings(int narg, char **arg)
{
  nstyles = narg;
  styles = new Angle*[nstyles];
  keywords = new char*[nstyles];

  for (int m = 0; m < nstyles; m++) {
    for (int i = 0; i < m; i++)
      if (strcmp(arg[m],arg[i]) == 0) 
	error->all("Angle style hybrid cannot use same angle style twice");
    if (strcmp(arg[m],"hybrid") == 0) 
      error->all("Angle style hybrid cannot have hybrid as an argument");
    styles[m] = force->new_angle(arg[m]);
    keywords[m] = new char[strlen(arg[m])+1];
    strcpy(keywords[m],arg[m]);
  }
}

/* ----------------------------------------------------------------------
   set coeffs for one type
---------------------------------------------------------------------- */

void AngleHybrid::coeff(int which, int narg, char **arg)
{
  if (!allocated) allocate();

  int ilo,ihi;
  force->bounds(arg[0],atom->nangletypes,ilo,ihi);

  // 2nd arg = angle style name (harmonic, etc)

  int m;
  for (m = 0; m < nstyles; m++)
    if (strcmp(arg[1],keywords[m]) == 0) break;
  if (m == nstyles) error->all("Angle coeff for hybrid has invalid style");

  // set low-level coefficients for each angletype
  // replace 2nd arg with i, call coeff() with no 1st arg
  // if sub-style is NULL for "none", still set setflag

  for (int i = ilo; i <= ihi; i++) {
    sprintf(arg[1],"%d",i);
    map[i] = m;
    if (styles[m]) styles[m]->coeff(which,narg-1,&arg[1]);
    if (styles[m] == NULL) setflag[i] = 1;
    else setflag[i] = styles[m]->setflag[i];
  }
}

/* ----------------------------------------------------------------------
   return an equilbrium angle length 
------------------------------------------------------------------------- */

double AngleHybrid::equilibrium_angle(int i)
{
  return styles[map[i]]->equilibrium_angle(i);
}

/* ----------------------------------------------------------------------
   proc 0 writes to restart file
------------------------------------------------------------------------- */

void AngleHybrid::write_restart(FILE *fp)
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

void AngleHybrid::read_restart(FILE *fp)
{
  allocate();

  int me = comm->me;
  if (me == 0) fread(&nstyles,sizeof(int),1,fp);
  MPI_Bcast(&nstyles,1,MPI_INT,0,world);
  styles = new Angle*[nstyles];
  keywords = new char*[nstyles];
  
  int n;
  for (int m = 0; m < nstyles; m++) {
    if (me == 0) fread(&n,sizeof(int),1,fp);
    MPI_Bcast(&n,1,MPI_INT,0,world);
    keywords[m] = new char[n];
    if (me == 0) fread(keywords[m],sizeof(char),n,fp);
    MPI_Bcast(keywords[m],n,MPI_CHAR,0,world);
    styles[m] = force->new_angle(keywords[m]);
  }
}

/* ---------------------------------------------------------------------- */

double AngleHybrid::single(int type, int i1, int i2, int i3, double rfactor)
{
  if (styles[map[type]]) 
    return styles[map[type]]->single(type,i1,i2,i3,rfactor);
  else return 0.0;
}

/* ----------------------------------------------------------------------
   memory usage
------------------------------------------------------------------------- */

int AngleHybrid::memory_usage()
{
  int bytes = 0;
  for (int m = 0; m < nstyles; m++) bytes += maxangle[m]*4 * sizeof(int);
  for (int m = 0; m < nstyles; m++) 
    if (styles[m]) bytes += styles[m]->memory_usage();
  return bytes;
}
