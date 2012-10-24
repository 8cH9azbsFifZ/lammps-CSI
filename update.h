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

#ifndef UPDATE_H
#define UPDATE_H

#include "pointers.h"

namespace LAMMPS_NS {

class Update : protected Pointers {
 public:
  double dt;                      // timestep
  double tolerance;               // minimizer tolerance
  int ntimestep;                  // current step (dynamics or min iter)
  int nsteps;                     // # of steps to run (dynamics or min iter)
  int firststep,laststep;         // 1st & last step of this run
  int beginstep,endstep;          // 1st and last step of multiple runs
  int first_update;               // 0 before initial update, 1 after
  int max_eval;                   // max force evaluations for minimizer
  int whichflag;                  // 0 for time integration, 1 for minimization

  double **f_pair;                // used by pair to compute force & virial
  int maxpair;

  char *unit_style;

  class Integrate *integrate;
  char *integrate_style;

  class Min *minimize;
  char *minimize_style;

  Update(class LAMMPS *);
  ~Update();
  void init();
  void set_units(char *);
  void create_integrate(int, char **);
  void create_minimize(int, char **);
  int memory_usage();
};

}

#endif
