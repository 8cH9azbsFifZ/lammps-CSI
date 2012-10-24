/* ----------------------------------------------------------------------
 * (C)opyright Gerolf Ziegenhain 2007
 *    mail.gerolf@ziegenhain.com
------------------------------------------------------------------------- */

#ifndef FIX_LTEMP_H
#define FIX_LTEMP_H

#include "fix.h"

namespace LAMMPS_NS {
   class FixLtemp : public Fix {
      friend class DumpCustom;

      public:
         FixLtemp(class LAMMPS *, int, char **);
         ~FixLtemp();
         int setmask();
         void init();
         void dump();
         int memory_usage();

      private:
         int nmax;
         double *ltemp;
   };
}

#endif
