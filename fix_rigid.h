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

#ifndef FIX_RIGID_H
#define FIX_RIGID_H

#include "fix.h"

namespace LAMMPS_NS {

class FixRigid : public Fix {
 public:
  FixRigid(class LAMMPS *, int, char **);
  ~FixRigid();
  int setmask();
  void init();
  void setup();
  void initial_integrate();
  void final_integrate();
  void initial_integrate_respa(int, int);
  void final_integrate_respa(int);

  int memory_usage();
  void grow_arrays(int);
  void copy_arrays(int, int);
  int pack_exchange(int, double *);
  int unpack_exchange(int, double *);

  int dof(int);
  void dilate(int, double, double, double, double);

 private:
  double dtv,dtf,dtq;
  double *step_respa;
  int pressure_flag;

  int nbody;                // # of rigid bodies
  int *nrigid;              // # of atoms in each rigid body
  double *masstotal;        // total mass of each rigid body
  double **xcm;             // coords of center-of-mass of each rigid body
  double **vcm;             // velocity of center-of-mass of each
  double **fcm;             // force on center-of-mass of each
  double **inertia;         // 3 principal components of inertia of each
  double **ex_space,**ey_space,**ez_space;
                            // principal axes of each in space coords
  double **angmom;          // angular momentum of each in space coords
  double **omega;           // angular velocity of each in space coords
  double **torque;          // torque on each rigid body in space coords
  double **quat;            // quaternion of each rigid body

  int *body;                // which body each atom is part of (-1 if none)
  double **displace;        // displacement of each atom in body coords

  double **sum,**all;       // work vectors for each rigid body

  int jacobi(double **, double *, double **);
  void rotate(double **, int, int, int, int, double, double);
  void qcreate(double **, double *);
  void multiply(double *, double *, double *);
  void normalize(double *);
  void richardson(double *, double *, double *, double *,
		  double *, double *, double *);
  void omega_from_mq(double *, double *, double *,
		     double *, double *, double *);
  void exyz_from_q(double *, double *, double *, double *);
  void set_xv(int);
  void set_v(int);
};

}

#endif
