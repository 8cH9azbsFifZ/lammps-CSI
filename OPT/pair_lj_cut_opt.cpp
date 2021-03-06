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
   Contributing authors:  
     James Fischer, High Performance Technologies, Inc.
     David Richie, Stone Ridge Technology
     Vincent Natoli, Stone Ridge Technology
------------------------------------------------------------------------- */

#include "pair_lj_cut_opt.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

PairLJCutOpt::PairLJCutOpt(LAMMPS *lmp) : PairLJCut(lmp) {}

/* ---------------------------------------------------------------------- */

void PairLJCutOpt::compute(int eflag, int vflag)
{
  if (eflag) {
    if (force->newton_pair) {
      switch (vflag) {
      case 0: return eval<1,0,1>();
      case 1: return eval<1,1,1>();
      case 2: return eval<1,2,1>();
      }
    } else {
      switch (vflag) {
      case 0: return eval<1,0,0>();
      case 1: return eval<1,1,0>();
      case 2: return eval<1,2,0>();
      }
    }
    
  } else {
    if (force->newton_pair) {
      switch (vflag) {
      case 0: return eval<0,0,1>();
      case 1: return eval<0,1,1>();
      case 2: return eval<0,2,1>();
      }
    } else {
      switch (vflag) {
      case 0: return eval<0,0,0>();
      case 1: return eval<0,1,0>();
      case 2: return eval<0,2,0>();
      }
    }
  }
}
