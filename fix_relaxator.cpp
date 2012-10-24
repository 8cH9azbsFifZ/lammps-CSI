/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under 
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.

   Copyright (C) 2007 G. Ziegenhain, gerolf@ziegenhain.com
------------------------------------------------------------------------- */

#include "string.h"
#include "stdlib.h"
#include "fix_relaxator.h"
#include "atom.h"
#include "update.h"
#include "respa.h"
#include "error.h"
#include "compute.h"
#include "modify.h"
#include "domain.h"
#include "neighbor.h"
#include <math.h>

#define INDEX_BAD -1

#define XX 0
#define YY 1
#define ZZ 2
#define XY 3
#define YX 3
#define YZ 4
#define ZY 4
#define ZX 5
#define XZ 5

#define X 0
#define Y 1
#define Z 2

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

FixRelaxator::FixRelaxator (LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
   if (narg != 3+9) error->all ("Illegal fix relaxator command");

   double GigaPascal2eVAngstrom3 = 0.0062415097;

   press_applied[XX] = atof ( arg[3] );
   press_applied[YY] = atof ( arg[4] );
   press_applied[ZZ] = atof ( arg[5] );
   press_applied[XY] = atof ( arg[6] );
   press_applied[YZ] = atof ( arg[7] );
   press_applied[ZX] = atof ( arg[8] );
   damping_factor = atof ( arg[9] );
   bulk_modulus = GigaPascal2eVAngstrom3 * atof ( arg[10] );
   shear_modulus = GigaPascal2eVAngstrom3 * atof ( arg[11] );


   restart_global = 1;
   pressure_every = 1;
   box_change = 1;

   if (domain->triclinic) error->all ("Doesn't work with triclinic box");
   if (domain->xperiodic || domain->yperiodic || domain->zperiodic) error->all ("Doesn't work with periodic BC");
}

/* ---------------------------------------------------------------------- */

int FixRelaxator::setmask ()
{
  int mask = 0;
  mask |= POST_FORCE; //INITIAL_INTEGRATE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixRelaxator::init ()
{
}

/* ---------------------------------------------------------------------- */

void FixRelaxator::setup ()
{
}

/* ---------------------------------------------------------------------- */

void FixRelaxator::min_setup ()
{
//  post_force ();
}

/* ---------------------------------------------------------------------- */

int FixRelaxator::FindIndexPressure ()
{
   IndexPressure = INDEX_BAD;
   int i;
   for (i = 0; i < modify->ncompute; i++) 
      if (strcmp (modify->compute[i]->style, "pressure") == 0)
         IndexPressure = i;
   if (IndexPressure == INDEX_BAD)
      error->warning ("FindIndexPressure: No compute index pressure found");
   return i;
}

/* ---------------------------------------------------------------------- */

int FixRelaxator::FindIndexTemperature ()
{
   IndexTemperature = INDEX_BAD;
   int i;
   for (i = 0; i < modify->ncompute; i++) 
      if (strcmp (modify->compute[i]->style, "temperature") == 0)
         IndexTemperature = i;
   if (IndexTemperature == INDEX_BAD)
      error->warning ("FindIndexPressure: No compute index temperature found");
   return i;
}

/* ---------------------------------------------------------------------- */

void FixRelaxator::DetermineStress ()
{
//   FindIndexTemperature ();
//   modify->compute[IndexTemperature]->compute_vector();
   double volume = domain->xprd * domain->yprd * domain->zprd;

   FindIndexPressure ();
   modify->compute[IndexPressure]->compute_vector();

   press[XX] = modify->compute[IndexPressure]->vector[XX]/volume - press_applied[XX];
   press[YY] = modify->compute[IndexPressure]->vector[YY]/volume - press_applied[YY];
   press[ZZ] = modify->compute[IndexPressure]->vector[ZZ]/volume - press_applied[ZZ];
   press[XY] = modify->compute[IndexPressure]->vector[XY]/volume - press_applied[XY];
   press[YZ] = modify->compute[IndexPressure]->vector[YZ]/volume - press_applied[YZ];
   press[ZX] = modify->compute[IndexPressure]->vector[ZX]/volume - press_applied[ZX];

   hydrostatic_pressure = (press[XX]+press[YY]+press[ZZ]) / 3.;
}

/* ---------------------------------------------------------------------- */

void FixRelaxator::DetermineStrain ()
{
   strain[XX] = hydrostatic_pressure/bulk_modulus + (press[XX]-hydrostatic_pressure)/shear_modulus;
   strain[YY] = hydrostatic_pressure/bulk_modulus + (press[YY]-hydrostatic_pressure)/shear_modulus;
   strain[ZZ] = hydrostatic_pressure/bulk_modulus + (press[ZZ]-hydrostatic_pressure)/shear_modulus;
   strain[XY] = press[XY]/shear_modulus;
   strain[YZ] = press[YZ]/shear_modulus;
   strain[ZX] = press[ZX]/shear_modulus;
   
   //printf ("\n   Pressure  %10.10f %10.10f %10.10f %10.10f %10.10f %10.10f\n", press[XX], press[YY] , press[ZZ],press[XY],press[YZ],press[ZX]);
   //printf ("\n   Strain    %10.10f %10.10f %10.10f %10.10f %10.10f %10.10f\n", strain[XX], strain[YY] , strain[ZZ],strain[XY],strain[YZ],strain[ZX]);
   //fflush (stdout);
}

/* ---------------------------------------------------------------------- */

void FixRelaxator::ResizeBox (double *r) 
{
  /* 
   if (domain->xperiodic) {
      if (r[X] < domain->boxlo[X]) domain->boxlo[X] = r[X];
      if (r[X] > domain->boxhi[X]) domain->boxhi[X] = r[X];
      domain->prd[X] = domain->xprd = domain->boxhi[X] - domain->boxlo[X];
   }
   if (domain->yperiodic) {
      if (r[Y] < domain->boxlo[Y]) domain->boxlo[Y] = r[Y];
      if (r[Y] > domain->boxhi[Y]) domain->boxhi[Y] = r[Y];
      domain->prd[Y] = domain->yprd = domain->boxhi[Y] - domain->boxlo[Y];
   }
   if (domain->zperiodic) {
      if (r[Z] < domain->boxlo[Z]) domain->boxlo[Z] = r[Z];
      if (r[Z] > domain->boxhi[Z]) domain->boxhi[Z] = r[Z];
      domain->prd[Z] = domain->zprd = domain->boxhi[Z] - domain->boxlo[Z];
   }
   */
   
   double oldlo, oldhi, ctr;
   if (domain->xperiodic && strain[XX] > 0.) {
      oldlo = domain->boxlo[X];
      oldhi = domain->boxhi[X];
      ctr = 0.5 * (oldlo + oldhi);
      domain->boxlo[X] = (oldlo-ctr)*(1.+.5*damping_factor*strain[XX]) + ctr;
      domain->boxhi[X] = (oldhi-ctr)*(1.+.5*damping_factor*strain[XX]) + ctr;
      domain->prd[X] = domain->xprd = domain->boxhi[X] - domain->boxlo[X];
   }
   if (domain->yperiodic && strain[YY] > 0.) {
      oldlo = domain->boxlo[Y];
      oldhi = domain->boxhi[Y];
      ctr = 0.5 * (oldlo + oldhi);
      domain->boxlo[Y] = (oldlo-ctr)*(1.+.5*damping_factor*strain[YY]) + ctr;
      domain->boxhi[Y] = (oldhi-ctr)*(1.+.5*damping_factor*strain[YY]) + ctr;
      domain->prd[Y] = domain->xprd = domain->boxhi[Y] - domain->boxlo[Y];
   }
   if (domain->zperiodic && strain[ZZ] > 0.) {
      oldlo = domain->boxlo[Z];
      oldhi = domain->boxhi[Z];
      ctr = 0.5 * (oldlo + oldhi);
      domain->boxlo[Z] = (oldlo-ctr)*(1.+.5*damping_factor*strain[ZZ]) + ctr;
      domain->boxhi[Z] = (oldhi-ctr)*(1.+.5*damping_factor*strain[ZZ]) + ctr;
      domain->prd[Z] = domain->xprd = domain->boxhi[Z] - domain->boxlo[Z];
   }

   /*
   double box[3], center[3];
   double box_[3];
   box[X] = domain->boxhi[X]-domain->boxlo[X];
   box[Y] = domain->boxhi[Y]-domain->boxlo[Y];
   box[Z] = domain->boxhi[Z]-domain->boxlo[Z];

   center[X] = .5*(domain->boxhi[X]+domain->boxlo[X]);
   center[Y] = .5*(domain->boxhi[Y]+domain->boxlo[Y]);
   center[Z] = .5*(domain->boxhi[Z]+domain->boxlo[Z]);
   */
/*
   box_[X] = strain[XX]*r[i][X] + strain[XY]*r[i][Y] + strain[XZ]*r[i][Z];
   box_[Y] = strain[YX]*r[i][X] + strain[YY]*r[i][Y] + strain[YZ]*r[i][Z];
   box_[Z] = strain[ZX]*r[i][X] + strain[ZY]*r[i][Y] + strain[ZZ]*r[i][Z];

   r[i][X] += damping_factor * r_[X];
   r[i][Y] += damping_factor * r_[Y];
   r[i][Z] += damping_factor * r_[Z];
*/
   //domain->print_box (" ");
}

/* ---------------------------------------------------------------------- */

void FixRelaxator::post_force (int vflag)
{
   int *mask = atom->mask;

   DetermineStress ();
   DetermineStrain ();

   double **r = atom->x;
   double r_[3];
   double r_hat[3];

   for (int i = 0; i < atom->nlocal; i++) //+atom->nghost; i++)
      if (mask[i] & groupbit) {
         double Norm = 1. / sqrt ( r[i][X]*r[i][X] + r[i][Y]*r[i][Y] + r[i][Z]*r[i][Z] );
         r_hat[X] = r[i][X]*Norm;
         r_hat[Y] = r[i][Y]*Norm;
         r_hat[Z] = r[i][Z]*Norm;

         r_[X] = strain[XX]*r_hat[X] + strain[XY]*r_hat[Y] + strain[XZ]*r_hat[Z];
         r_[Y] = strain[YX]*r_hat[X] + strain[YY]*r_hat[Y] + strain[YZ]*r_hat[Z];
         r_[Z] = strain[ZX]*r_hat[X] + strain[ZY]*r_hat[Y] + strain[ZZ]*r_hat[Z];
         /*
         r_[X] = strain[XX]*r[i][X] + strain[XY]*r[i][Y] + strain[XZ]*r[i][Z];
         r_[Y] = strain[YX]*r[i][X] + strain[YY]*r[i][Y] + strain[YZ]*r[i][Z];
         r_[Z] = strain[ZX]*r[i][X] + strain[ZY]*r[i][Y] + strain[ZZ]*r[i][Z];
         */
         //printf ("\n   -> %d %f %f %f %f %f %f", i, r[i][X], r[i][Y], r[i][Z], r_[X], r_[Y], r_[Z]);
         //fflush (stdout);

         r[i][X] += damping_factor * r_[X];
         r[i][Y] += damping_factor * r_[Y];
         r[i][Z] += damping_factor * r_[Z];

//         if (i < atom->nlocal)
            ResizeBox (r[i]);
//         domain->remap (atom->x[i], atom->image[i]);
      }

  // neighbor->build_full ();

//   ResizeBox (r[1]);
}

