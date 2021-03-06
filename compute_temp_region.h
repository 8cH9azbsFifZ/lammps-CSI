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

#ifndef COMPUTE_TEMP_REGION_H
#define COMPUTE_TEMP_REGION_H

#include "compute.h"

namespace LAMMPS_NS {

class ComputeTempRegion : public Compute {
 public:
  ComputeTempRegion(class LAMMPS *, int, char **);
  ~ComputeTempRegion();
  void init();
  double compute_scalar();
  void compute_vector();

 private:
  int iregion;
};

}

#endif
