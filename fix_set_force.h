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

#ifndef FIX_SET_FORCE_H
#define FIX_SET_FORCE_H

#include "fix.h"

namespace LAMMPS_NS {

class FixSetForce : public Fix {
 public:
  FixSetForce(class LAMMPS *, int, char **);
  int setmask();
  void init();
  void setup();
  void min_setup();
  void post_force(int);
  void post_force_respa(int, int, int);
  void min_post_force(int);
  double thermo(int);

 private:
  int flagx,flagy,flagz;
  double xvalue,yvalue,zvalue;
  double foriginal[3];
  int nlevels_respa;
};

}

#endif
