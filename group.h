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

#ifndef GROUP_H
#define GROUP_H

#include "stdio.h"
#include "pointers.h"

namespace LAMMPS_NS {

class Group : protected Pointers {
 public:
  int me;
  int ngroup;                  // # of defined groups
  char **names;                // name of each group
  int *bitmask;                // one-bit mask for each group
  int *inversemask;            // inverse mask for each group

  Group(class LAMMPS *);
  ~Group();
  void assign(int, char **);         // assign atoms to a group
  void create(char *, int *);        // add flagged atoms to a group
  int find(char *);                  // lookup name in list of groups
  void write_restart(FILE *);
  void read_restart(FILE *);

  double count(int);                       // count atoms in group
  double mass(int);                        // total mass of atoms in group
  double charge(int);                      // total charge of atoms in group
  void bounds(int, double *);              // bounds of atoms in group
  void xcm(int, double, double *);         // center-of-mass coords of group
  void vcm(int, double, double *);         // center-of-mass velocity of group
  void fcm(int, double *);                 // total force on group
  double ke(int);                          // kinetic energy of group
  double gyration(int, double, double *);  // radius-of-gyration of group
  void angmom(int, double *, double *);    // angular momentum of group
  void inertia(int, double *, double [3][3]);     // inertia tensor
  void omega(double *, double [3][3], double *);  // angular velocity
};

}

#endif
