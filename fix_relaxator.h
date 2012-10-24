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

#ifndef FIX_RELAXATOR_H
#define FIX_RELAXATOR_H

#include "fix.h"
#include "pointers.h"

namespace LAMMPS_NS {

class FixRelaxator : public Fix {
   public:
      FixRelaxator (class LAMMPS *, int, char **);
      int setmask ();
      void init ();
      void setup ();
      void min_setup ();
      void post_force (int vflag);

   private:
      double bulk_modulus;
      double shear_modulus;
      double damping_factor;
      double press[6], hydrostatic_pressure, press_applied[6];
      int IndexPressure, IndexTemperature;
      double strain[6];

      int FindIndexPressure ();
      int FindIndexTemperature ();
      void DetermineStress ();
      void DetermineStrain ();
      void ResizeBox (double *r);
};

}

#endif
