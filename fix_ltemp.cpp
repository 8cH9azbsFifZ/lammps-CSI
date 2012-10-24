/* ----------------------------------------------------------------------
 * (C)opyright Gerolf Ziegenhain 
 *    mail.gerolf@ziegenhain.com
------------------------------------------------------------------------- */

#include "string.h"
#include "fix_ltemp.h"
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

FixLtemp::FixLtemp(LAMMPS *lmp, int narg, char **arg) : Fix(lmp, narg, arg) {
  if (narg != 3) error->all("Illegal fix ltemp command");

  neigh_full_once = 1;

  ltemp = NULL;
  nmax = 0;
}

/* ---------------------------------------------------------------------- */

FixLtemp::~FixLtemp() {
  memory->sfree (ltemp);
}

/* ---------------------------------------------------------------------- */

int FixLtemp::setmask() {
  int mask = 0;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixLtemp::init() {
   // Check for only one fix
  int count = 0;
  for (int i = 0; i < modify->nfix; i++)
    if (strcmp(modify->fix[i]->style,"LTEMP") == 0) count++;
  if (count > 1 && comm->me == 0)
    error->warning("More than one dump custom with a ltemp attribute");
}

/* ---------------------------------------------------------------------- */

void FixLtemp::dump() {
   int nlocal = atom->nlocal;
   int *mask = atom->mask;

   // Allocate memory
   if (nlocal > nmax) {
      memory->sfree (ltemp);
      nmax = atom->nmax;
      ltemp = (double *) memory->smalloc (nmax*sizeof (double), "ltemp:ltemp");
   }

   // If needed: build a full neighbor list
   if (!neighbor->full_every) neighbor->build_full();

   int i, j;
   int *neighs;
   int numneigh;
   double factor;
   double **v = atom->v;
   double v1, v2, v3;
   for (i = 0; i < nlocal; i++) 
      if (mask[i] & groupbit) {
         ltemp[i] = 0.0;
         for (j = 0; j < neighbor->numneigh_full[i]; j++) {
            if (atom->mass)
               factor = atom->mass[atom->type[j]];
            else
               factor = atom->rmass[j];

            v1 = v[j][0]-v[i][0];
            v1 = v[j][1]-v[i][1];
            v1 = v[j][2]-v[i][2];

            ltemp[i] += ( v1*v1 + v2*v2 + v3*v3 ) * factor;
         }
         ltemp[i] *= force->mvv2e / ( 3.0 * neighbor->numneigh_full[i] * force->boltz );
      }
}

/* ---------------------------------------------------------------------- */

int FixLtemp::memory_usage() {
   int bytes = nmax * sizeof (double);
   return bytes;
}
