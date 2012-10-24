/* ----------------------------------------------------------------------
 * (C)opyright Gerolf Ziegenhain 
 *    mail.gerolf@ziegenhain.com
------------------------------------------------------------------------- */

#include "string.h"
#include "fix_ldens.h"
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

FixLdens::FixLdens(LAMMPS *lmp, int narg, char **arg) : Fix(lmp, narg, arg) {
  if (narg != 3) error->all("Illegal fix ldens command");

  neigh_full_once = 1;

  ldens = NULL;
  nmax = 0;
}

/* ---------------------------------------------------------------------- */

FixLdens::~FixLdens() {
  memory->sfree (ldens);
}

/* ---------------------------------------------------------------------- */

int FixLdens::setmask() {
  int mask = 0;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixLdens::init() {
   // Check for only one fix
  int count = 0;
  for (int i = 0; i < modify->nfix; i++)
    if (strcmp(modify->fix[i]->style,"LDENS") == 0) count++;
  if (count > 1 && comm->me == 0)
    error->warning("More than one dump custom with a ldens attribute");
}

/* ---------------------------------------------------------------------- */

void FixLdens::dump() {
   int nlocal = atom->nlocal;
   int *mask = atom->mask;

   // Allocate memory
   if (nlocal > nmax) {
      memory->sfree (ldens);
      nmax = atom->nmax;
      ldens = (double *) memory->smalloc (nmax*sizeof (double), "ldens:ldens");
   }

   // If needed: build a full neighbor list
   if (!neighbor->full_every) neighbor->build_full();

   double fac = 3.0 / (4.0 * 3.1415);
   int i;
   for (i = 0; i < nlocal; i++) {
      if (mask[i] & groupbit) {
         double rcut = *(neighbor->cuttype+atom->type[i]);
         ldens[i] = fac * (double)neighbor->numneigh_full[i] / (rcut*rcut*rcut);
      }
   }
}

/* ---------------------------------------------------------------------- */

int FixLdens::memory_usage() {
   int bytes = nmax * sizeof (double);
   return bytes;
}
