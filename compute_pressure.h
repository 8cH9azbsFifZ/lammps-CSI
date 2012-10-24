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

#ifndef COMPUTE_PRESSURE_H
#define COMPUTE_PRESSURE_H

#include "compute.h"

namespace LAMMPS_NS {

class ComputePressure : public Compute {
 public:
  ComputePressure(class LAMMPS *, int, char **);
  ~ComputePressure();
  void init();
  double compute_scalar();
  void compute_vector();

 private:
  double boltz,nktv2p,inv_volume;
  int nvirial;
  double **vptr;
  int kspaceflag;
  double *kspace_virial;
  Compute *temperature;
  double virial[6];

  void virial_compute(int);
};

}

#endif
