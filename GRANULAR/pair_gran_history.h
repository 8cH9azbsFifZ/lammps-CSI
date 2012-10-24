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

#ifndef PAIR_GRAN_HISTORY_H
#define PAIR_GRAN_HISTORY_H

#include "pair.h"

namespace LAMMPS_NS {

class PairGranHistory : public Pair {
 public:
  PairGranHistory(class LAMMPS *);
  ~PairGranHistory();
  virtual void compute(int, int);
  void settings(int, char **);
  void coeff(int, char **);
  double init_one(int, int);
  void init_style();
  void write_restart(FILE *);
  void read_restart(FILE *);
  void write_restart_settings(FILE *);
  void read_restart_settings(FILE *);

  void extract_gran(double *, double *, double *, int *);

 protected:
  double xkk,xkkt,xmu;
  int dampflag;
  double gamman;
  double dt,gamman_dl,gammas_dl;
  int freeze_group_bit;
  int history;
  class FixShearHistory *fix_history;

  void allocate();
};

}

#endif
