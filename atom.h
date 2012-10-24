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

#ifndef ATOM_H
#define ATOM_H

#include "pointers.h"

namespace LAMMPS_NS {

class Atom : protected Pointers {
 public:
  char *atom_style;
  class AtomVec *avec;

  // atom counts

  double natoms;                // total # of atoms in system, could be 0
  int nlocal,nghost;            // # of owned and ghost atoms on this proc
  int nmax;                     // max # of owned+ghost in arrays on this proc
  int tag_enable;               // 0/1 if atom ID tags are defined
  int molecular;                // 0 = atomic, 1 = molecular system

  int ntypes,nbondtypes,nangletypes,ndihedraltypes,nimpropertypes;
  int nbonds,nangles,ndihedrals,nimpropers;
  int bond_per_atom,angle_per_atom,dihedral_per_atom,improper_per_atom;

  // per-atom arrays
  // customize by adding new array

  int *tag,*type,*mask,*image;
  double **x,**v,**f;

  int *molecule;
  double *q,**mu;
  double **xorient,**quat,**omega,**angmom,**torque;
  double *radius,*density,*rmass,*vfrac;

  int maxspecial;
  int **nspecial,**special;

  int *num_bond;
  int **bond_type;
  int **bond_atom;

  int *num_angle;
  int **angle_type;
  int **angle_atom1,**angle_atom2,**angle_atom3;

  int *num_dihedral;
  int **dihedral_type;
  int **dihedral_atom1,**dihedral_atom2,**dihedral_atom3,**dihedral_atom4;

  int *num_improper;
  int **improper_type;
  int **improper_atom1,**improper_atom2,**improper_atom3,**improper_atom4;

  // per-atom array existence flags
  // customize by adding new flag

  int molecule_flag;
  int q_flag,mu_flag;
  int xorient_flag,quat_flag,omega_flag,angmom_flag,torque_flag;
  int radius_flag,density_flag,rmass_flag,vfrac_flag;

  // extra peratom info in restart file destined for fix & diag 

  double **extra;

  // per-type arrays

  double *mass,**shape,*dipole;
  int *mass_setflag,*shape_setflag,*dipole_setflag;

  // callback ptrs for atom arrays managed by fix classes

  int nextra_grow,nextra_restart;             // # of callbacks of each type
  int *extra_grow,*extra_restart;             // index of fix to callback to
  int nextra_grow_max,nextra_restart_max;     // size of callback lists
  int nextra_store;

  int map_style;                  // default or user-specified style of map
                                  // 0 = none, 1 = array, 2 = hash

  // functions

  Atom(class LAMMPS *);
  ~Atom();

  void create_avec(char *, int, char **);
  class AtomVec *new_avec(char *, int, char **);
  void init();

  int style_match(char *);
  void modify_params(int, char **);
  void tag_extend();
  int tag_consecutive();

  int parse_data(char *);
  int count_words(char *);

  void data_atoms(int, char *);
  void data_vels(int, char *);
  void data_bonds(int, char *);
  void data_angles(int, char *);
  void data_dihedrals(int, char *);
  void data_impropers(int, char *);

  void allocate_type_arrays();
  void set_mass(char *);
  void set_mass(int, double);
  void set_mass(int, char **);
  void set_mass(double *);
  void check_mass();
  void set_shape(char *);
  void set_shape(int, char **);
  void set_shape(double **);
  void check_shape();
  void set_dipole(char *);
  void set_dipole(int, char **);
  void set_dipole(double *);
  void check_dipole();

  void add_callback(int);
  void delete_callback(char *, int);
  void update_callback(int);

  int memory_usage();
  int memcheck(const char *);

  // functions for global to local ID mapping
  // map lookup function inlined for efficiency
  
  inline int map(int global) {
    if (map_style == 1) return map_array[global];
    else return map_find_hash(global);
  };

  void map_init();
  void map_clear();
  void map_set();
  void map_one(int, int);
  void map_delete();
  int map_find_hash(int);

 private:

  // data for global to local ID mapping

  int map_tag_max;
  int *map_array;

  struct HashElem {
    int global;                   // key to search on = global ID
    int local;                    // value associated with key = local index
    int next;                     // next entry in this bucket, -1 if last
  };
  int map_nhash;                  // # of entries hash table can hold
  int map_nused;                  // # of actual entries in hash table
  int map_free;                   // ptr to 1st unused entry in hash table
  int map_nbucket;                // # of hash buckets
  int *map_bucket;                // ptr to 1st entry in each bucket
  HashElem *map_hash;             // hash table
  int *primes;                    // table of prime #s for hashing
  int nprimes;                    // # of primes

  int memlength;                  // allocated size of memstr
  char *memstr;                   // string of array names already counted
};

}

#endif
