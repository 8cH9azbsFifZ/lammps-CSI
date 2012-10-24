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
   Contributing author: Ravi Agrawal (Northwestern U)
------------------------------------------------------------------------- */

#include "math.h"
#include "string.h"
#include "stdlib.h"
#include "fix_indent.h"
#include "atom.h"
#include "domain.h"
#include "lattice.h"
#include "update.h"
#include "output.h"
#include "respa.h"
#include "error.h"
#include "mpi.h"

#define PI 3.1415 // for contactarea
#define MAX(A,B) ((A) > (B)) ? (A) : (B)
#define MAGIC_MAX -1000.0
#define MAGIC_MIN 1000.0

using namespace LAMMPS_NS;

enum{NONE,SPHERE,CYLINDER};

/* ---------------------------------------------------------------------- */

FixIndent::FixIndent(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  if (narg < 4) error->all("Illegal fix indent command");
  k = atof(arg[3]);

  // set defaults

  istyle = NONE;
  vx = vy = vz = 0.0;
  radflag = 0;
  r0_start = 0.0;
  scaleflag = 1;

  // read options from end of input line

  options(narg-4,&arg[4]);

  // setup scaling

  if (scaleflag && domain->lattice == NULL)
    error->all("Use of fix indent with undefined lattice");

  double xscale,yscale,zscale;
  if (scaleflag) {
    xscale = domain->lattice->xlattice;
    yscale = domain->lattice->ylattice;
    zscale = domain->lattice->zlattice;
  }
  else xscale = yscale = zscale = 1.0;

  // apply scaling factors to force constant, velocity, and geometry

  k /= xscale;
  k3 = k/3.0;
  vx *= xscale;
  vy *= yscale;
  vz *= zscale;

  if (istyle == SPHERE) {
    x0 *= xscale;
    y0 *= yscale;
    z0 *= zscale;
    r0_stop *= xscale;
    r0_start *= xscale;
  } else if (istyle == CYLINDER) {
    if (cdim == 0) {
      c1 *= yscale;
      c2 *= zscale;
      r0_stop *= xscale;
      r0_start *= xscale;
    } else if (cdim == 1) {
      c1 *= xscale;
      c2 *= zscale;
      r0_stop *= yscale;
      r0_start *= yscale;
    } else if (cdim == 2) {
      c1 *= xscale;
      c2 *= yscale;
      r0_stop *= zscale;
      r0_start *= zscale;
    }
  } else error->all("Illegal fix indent command");
}

/* ---------------------------------------------------------------------- */

int FixIndent::setmask()
{
  int mask = 0;
  mask |= POST_FORCE;
  mask |= THERMO_ENERGY;
  mask |= POST_FORCE_RESPA;
  mask |= MIN_POST_FORCE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixIndent::init()
{
  if (strcmp(update->integrate_style,"respa") == 0)
    nlevels_respa = ((Respa *) update->integrate)->nlevels;
}

/* ---------------------------------------------------------------------- */

void FixIndent::setup()
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

void FixIndent::min_setup()
{
  eflag_enable = 1;
  post_force(1);
}

/* ---------------------------------------------------------------------- */

void FixIndent::post_force(int vflag)
{
  bool eflag = false;
  if (eflag_enable) eflag = true;
  else if (output->next_thermo == update->ntimestep) eflag = true;

  // set current r0
  // for minimization, always set to r0_stop

  double r0;
  if (!radflag || update->whichflag) r0 = r0_stop;
  else {
    double delta = update->ntimestep - update->beginstep;
    delta /= update->endstep - update->beginstep;
    r0 = r0_start + delta * (r0_stop-r0_start);
  }

  double eng;
  if (eflag) eng = 0.0;

  // spherical indenter
  
  /***** Indentation Properties *****/
  double FIndXLocal = 0., FIndYLocal = 0., FIndZLocal = 0.;
  double IndXLocal, IndYLocal, IndZLocal;
  double XContactMinLocal = MAGIC_MIN, XContactMaxLocal = MAGIC_MAX,
         YContactMinLocal = MAGIC_MIN, YContactMaxLocal = MAGIC_MAX,
         ZContactMinLocal = MAGIC_MIN, ZContactMaxLocal = MAGIC_MAX,
         XContactMin = MAGIC_MIN, XContactMax = MAGIC_MAX,
         YContactMin = MAGIC_MIN, YContactMax = MAGIC_MAX,
         ZContactMin = MAGIC_MIN, ZContactMax = MAGIC_MAX;
  double AcAtomarProjectedLocal = 0.;
  int NcLocal = 0;
  /***** *****/

  if (istyle == SPHERE) {

    // x1,y1,z1 = current position of indenter from original x0,y0,z0

    double delta = (update->ntimestep - update->beginstep) * update->dt;
    double x1 = x0 + delta*vx;
    double y1 = y0 + delta*vy;
    double z1 = z0 + delta*vz;

    double **x = atom->x;
    double **f = atom->f;
    int *mask = atom->mask;
    int nlocal = atom->nlocal;

    double delx,dely,delz,r,dr,fmag;

    /***** Indentation Properties *****/
    IndXLocal = x1;
    IndYLocal = y1;
    IndZLocal = z1;

    double drx, dry, drz;
    /***** *****/

    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) {
	      delx = x[i][0] - x1;
      	dely = x[i][1] - y1;
      	delz = x[i][2] - z1;
      	r = sqrt(delx*delx + dely*dely + delz*delz);
      	dr = r - r0;
      	if (dr >= 0.0) continue;
      	fmag = k*dr*dr;
         drx = delx/r;
         dry = dely/r;
         drz = delz/r;
      	f[i][0] += drx*fmag;
      	f[i][1] += dry*fmag;
      	f[i][2] += drz*fmag;
      	if (eflag) eng -= k3 * dr*dr*dr;

         /***** Indentation Properties *****/
         FIndXLocal += drx*fmag;
         FIndYLocal += dry*fmag;
         FIndZLocal += drz*fmag;
         
         NcLocal++;
        
         // Extremal contact positions
         if (x[i][0] > XContactMaxLocal) { XContactMaxLocal = x[i][0]; }
         if (x[i][0] < XContactMinLocal) { XContactMinLocal = x[i][0]; }

         if (x[i][1] > YContactMaxLocal) { YContactMaxLocal = x[i][1]; }
         if (x[i][1] < YContactMinLocal) { YContactMinLocal = x[i][1]; }

         if (x[i][2] > ZContactMaxLocal) { ZContactMaxLocal = x[i][2]; }
         if (x[i][2] < ZContactMinLocal) { ZContactMinLocal = x[i][2]; }


         AcAtomarProjectedLocal -= dry;
      }
      /***** *****/


  // cylindrical indenter

  } else {

    // c1new,c2new = current coords of indenter axis from original c1,c2
	      
    double delta = (update->ntimestep - update->beginstep) * update->dt;
    double c1new,c2new;
    if (cdim == 0) {
      c1new = c1 + delta*vy;
      c2new = c2 + delta*vz;
    } else if (cdim == 1) {
      c1new = c1 + delta*vx;
      c2new = c2 + delta*vz;
    } else {
      c1new = c1 + delta*vx;
      c2new = c2 + delta*vy;
    }
    
    double **x = atom->x;
    double **f = atom->f;
    int *mask = atom->mask;
    int nlocal = atom->nlocal;
    
    double delx,dely,delz,r,dr,fmag;
    
    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) {
	if (cdim == 0) {
	  delx = 0;
	  dely = x[i][1] - c1new;
	  delz = x[i][2] - c2new;
	} else if (cdim == 1) {
	  delx = x[i][0] - c1new;
	  dely = 0;
	  delz = x[i][2] - c2new;
	} else {
	  delx = x[i][0] - c1new;
	  dely = x[i][1] - c2new;
	  delz = 0;
	}
	r = sqrt(delx*delx + dely*dely + delz*delz);
	dr = r - r0;
	if (dr >= 0.0) continue;
	fmag = k*dr*dr;
	f[i][0] += delx*fmag/r;
	f[i][1] += dely*fmag/r;
	f[i][2] += delz*fmag/r;
	if (eflag) eng -= k3 * dr*dr*dr;
      }
  }

  if (eflag) MPI_Allreduce(&eng,&etotal,1,MPI_DOUBLE,MPI_SUM,world);

  /***** Calc Global Indentation Properties *****/
  
  // Force and position
# define PATCH_INDENTATION
# ifdef PATCH_INDENTATION  
  FIndX = 0.0;
  MPI_Allreduce(&FIndXLocal,&FIndX,1,MPI_DOUBLE,MPI_SUM,world);
  FIndY = 0.0;
  MPI_Allreduce(&FIndYLocal,&FIndY,1,MPI_DOUBLE,MPI_SUM,world);
  FIndZ = 0.0;
  MPI_Allreduce(&FIndZLocal,&FIndZ,1,MPI_DOUBLE,MPI_SUM,world);
  // FIXME: Check if this is needed!
  //if (IndXLocal >= 0.0) MPI_Allreduce(&IndXLocal,&IndX,1,MPI_DOUBLE,MPI_MAX,world);
  //else MPI_Allreduce(&IndXLocal,&IndX,1,MPI_DOUBLE,MPI_MIN,world);
  //if (IndYLocal >= 0.0) MPI_Allreduce(&IndYLocal,&IndY,1,MPI_DOUBLE,MPI_MAX,world);
  //else MPI_Allreduce(&IndYLocal,&IndY,1,MPI_DOUBLE,MPI_MIN,world);
  //if (IndZLocal >= 0.0) MPI_Allreduce(&IndZLocal,&IndZ,1,MPI_DOUBLE,MPI_MAX,world);
  //else MPI_Allreduce(&IndZLocal,&IndZ,1,MPI_DOUBLE,MPI_MIN,world);
  IndX = IndXLocal;
  IndY = IndYLocal;
  IndZ = IndZLocal;


  // Extremal Positions
   MPI_Allreduce (&XContactMaxLocal, &XContactMax, 1, MPI_DOUBLE, MPI_MAX, world);
   MPI_Allreduce (&YContactMinLocal, &YContactMin, 1, MPI_DOUBLE, MPI_MIN, world);
   MPI_Allreduce (&ZContactMaxLocal, &ZContactMax, 1, MPI_DOUBLE, MPI_MAX, world);
   MPI_Allreduce (&XContactMinLocal, &XContactMin, 1, MPI_DOUBLE, MPI_MIN, world);
   MPI_Allreduce (&YContactMaxLocal, &YContactMax, 1, MPI_DOUBLE, MPI_MAX, world);
   MPI_Allreduce (&ZContactMinLocal, &ZContactMin, 1, MPI_DOUBLE, MPI_MIN, world);

  if (YContactMax != MAGIC_MAX && YContactMin != MAGIC_MIN)
     PlasticDepth = YContactMax-YContactMin;
  else 
     PlasticDepth = 0.0;
# ifdef DEBUG     
  int rnk = 0;
  MPI_Comm_rank(MPI_COMM_WORLD,&rnk);
  if (rnk == 0) {
     printf ("\n->     %f %f %f %f", YContactMin, YContactMax, PlasticDepth, IndY);
     fflush(stdout);
  }
# endif

  // Atomar contact areas
  Nc = 0;
  MPI_Allreduce(&NcLocal,&Nc,1,MPI_INT,MPI_SUM,world);
  AcAtomarMax = (double)Nc * PI;
  AcAtomarProjected = 0.0;
  MPI_Allreduce(&AcAtomarProjectedLocal,&AcAtomarProjected,1,MPI_DOUBLE,MPI_SUM,world);
  AcAtomarProjected *= PI;

  // Continuous contact areas
  double Rind = r0;
  AcMeyer = PI*(2.0*Rind*PlasticDepth-PlasticDepth*PlasticDepth);
  AcBrinell = 2.0*PI*Rind*PlasticDepth;
  if (XContactMax == MAGIC_MAX || YContactMax == MAGIC_MAX || ZContactMax == MAGIC_MAX ||
        XContactMin == MAGIC_MIN || YContactMin == MAGIC_MIN || ZContactMin == MAGIC_MIN)
     AcElliptic = 0.;
  else
      AcElliptic = PI*.25*(XContactMax-XContactMin)*(ZContactMax-ZContactMin);
# endif
  /***** *****/
}

/* ---------------------------------------------------------------------- */

void FixIndent::post_force_respa(int vflag, int ilevel, int iloop)
{
  if (ilevel == nlevels_respa-1) post_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixIndent::min_post_force(int vflag)
{
  post_force(vflag);
}

/* ---------------------------------------------------------------------- */

double FixIndent::thermo(int n)
{
  if (n == 0) return etotal;
  else return 0.0;
}

/* ---------------------------------------------------------------------- */

int FixIndent::thermo_fields(int n, int *flags, char **keywords)
{
   if (n == 0) return 1;
   flags[0] = 3;
   strcpy(keywords[0],"Indent");
   return 1;
}
  
/* ---------------------------------------------------------------------- */


/* ----------------------------------------------------------------------
   parse optional parameters at end of input line 
------------------------------------------------------------------------- */

void FixIndent::options(int narg, char **arg)
{
  if (narg < 0) error->all("Illegal fix indent command");

  int iarg = 0;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"sphere") == 0) {
      if (iarg+5 > narg) error->all("Illegal fix indent command");
      x0 = atof(arg[iarg+1]);
      y0 = atof(arg[iarg+2]);
      z0 = atof(arg[iarg+3]);
      r0_stop = atof(arg[iarg+4]);
      istyle = SPHERE;
      iarg += 5;
    } else if (strcmp(arg[iarg],"cylinder") == 0) {
      if (iarg+5 > narg) error->all("Illegal fix indent command");
      if (strcmp(arg[iarg+1],"x") == 0) cdim = 0;
      else if (strcmp(arg[iarg+1],"y") == 0) cdim = 1;
      else if (strcmp(arg[iarg+1],"z") == 0) cdim = 2;
      else error->all("Illegal fix indent command");
      c1 = atof(arg[iarg+2]);
      c2 = atof(arg[iarg+3]);
      r0_stop = atof(arg[iarg+4]);
      istyle = CYLINDER;
      iarg += 5;
    } else if (strcmp(arg[iarg],"vel") == 0) {
      if (iarg+4 > narg) error->all("Illegal fix indent command");
      vx = atof(arg[iarg+1]);
      vy = atof(arg[iarg+2]);
      vz = atof(arg[iarg+3]);
      iarg += 4;
    } else if (strcmp(arg[iarg],"rstart") == 0) {
      if (iarg+2 > narg) error->all("Illegal fix indent command");
      radflag = 1;
      r0_start = atof(arg[iarg+1]);
      iarg += 2;
    } else if (strcmp(arg[iarg],"units") == 0) {
      if (iarg+2 > narg) error->all("Illegal fix indent command");
      if (strcmp(arg[iarg+1],"box") == 0) scaleflag = 0;
      else if (strcmp(arg[iarg+1],"lattice") == 0) scaleflag = 1;
      else error->all("Illegal fix indent command");
      iarg += 2;
    } else error->all("Illegal fix indent command");
  }
}
