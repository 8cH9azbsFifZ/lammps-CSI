/* ----------------------------------------------------------------------
 * (C)opyright Gerolf Ziegenhain 2007
 *    mail.gerolf@ziegenhain.com
------------------------------------------------------------------------- */

#ifndef FIX_LDENS_H
#define FIX_LDENS_H

#include "fix.h"

namespace LAMMPS_NS {
   class FixLdens : public Fix {
      friend class DumpCustom;

      public:
         FixLdens(class LAMMPS *, int, char **);
         ~FixLdens();
         int setmask();
         void init();
         void dump();
         int memory_usage();

      private:
         int nmax;
         double *ldens;
   };
}

#endif
