/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov
                  
   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under 
   the GNU General Public License.
                                        
   See the README file in the top-level LAMMPS directory.
                                               
   This module identifies the local lattice structure using the method
   supposed in PRB(2006)73:054104.
                                                             
   Copyright (C) 2007 G. Ziegenhain, gerolf@ziegenhain.com
                                      
------------------------------------------------------------------------- */


#ifndef COMPUTE_ACKLAND_ATOM_H
#define COMPUTE_ACKLAND_ATOM_H

#include "compute.h"

namespace LAMMPS_NS {

class ComputeAcklandAtom : public Compute {
 public:
  ComputeAcklandAtom(class LAMMPS *, int, char **);
  ~ComputeAcklandAtom();
  void init();
  void compute_peratom();
  int memory_usage();

  // Public concentrations for thermal output
  double ConcentrationUNKNOWN,
      ConcentrationBCC,
      ConcentrationFCC,
      ConcentrationHCP,
      ConcentrationICO;

 private:
  int nmax,maxneigh;
  double *distsq;
  int *nearest, *nearest_n0, *nearest_n1;
  int *chi;
  double *structure;

  void select(int, int, double *);
  void select2(int, int, double *, int *);
};

}

#endif
