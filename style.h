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

#ifdef AngleInclude
#endif

#ifdef AngleClass
#endif

#ifdef AtomInclude
#include "atom_vec_atomic.h"
#include "atom_vec_charge.h"
#include "atom_vec_hybrid.h"
#endif

#ifdef AtomClass
AtomStyle(atomic,AtomVecAtomic)
AtomStyle(charge,AtomVecCharge)
AtomStyle(hybrid,AtomVecHybrid)
#endif

#ifdef BondInclude
#endif

#ifdef BondClass
#endif

#ifdef CommandInclude
#include "create_atoms.h"
#include "create_box.h"
#include "delete_atoms.h"
#include "delete_bonds.h"
#include "displace_atoms.h"
#include "minimize.h"
#include "read_data.h"
#include "read_restart.h"
#include "replicate.h"
#include "run.h"
#include "set.h"
#include "shell.h"
#include "temper.h"
#include "velocity.h"
#include "write_restart.h"
#endif

#ifdef CommandClass
CommandStyle(create_atoms,CreateAtoms)
CommandStyle(create_box,CreateBox)
CommandStyle(delete_atoms,DeleteAtoms)
CommandStyle(delete_bonds,DeleteBonds)
CommandStyle(displace_atoms,DisplaceAtoms)
CommandStyle(minimize,Minimize)
CommandStyle(read_data,ReadData)
CommandStyle(read_restart,ReadRestart)
CommandStyle(replicate,Replicate)
CommandStyle(run,Run)
CommandStyle(set,Set)
CommandStyle(shell,Shell)
CommandStyle(temper,Temper)
CommandStyle(velocity,Velocity)
CommandStyle(write_restart,WriteRestart)
#endif

#ifdef ComputeInclude
#include "compute_centro_atom.h"
#include "compute_ackland_atom.h"   
#include "compute_coord_atom.h"
#include "compute_epair_atom.h"
#include "compute_etotal_atom.h"
#include "compute_ke_atom.h"
#include "compute_pressure.h"
#include "compute_rotate_dipole.h"
#include "compute_rotate_gran.h"
#include "compute_stress_atom.h"
#include "compute_temp.h"
#include "compute_temp_deform.h"
#include "compute_temp_partial.h"
#include "compute_temp_ramp.h"
#include "compute_temp_region.h"
#include "compute_variable.h"
#include "compute_variable_atom.h"
#endif

#ifdef ComputeClass
ComputeStyle(centro/atom,ComputeCentroAtom)
ComputeStyle(ackland/atom,ComputeAcklandAtom)
ComputeStyle(coord/atom,ComputeCoordAtom)
ComputeStyle(epair/atom,ComputeEpairAtom)
ComputeStyle(etotal/atom,ComputeEtotalAtom)
ComputeStyle(ke/atom,ComputeKEAtom)
ComputeStyle(pressure,ComputePressure)
ComputeStyle(rotate/dipole,ComputeRotateDipole)
ComputeStyle(rotate/gran,ComputeRotateGran)
ComputeStyle(stress/atom,ComputeStressAtom)
ComputeStyle(temp,ComputeTemp)
ComputeStyle(temp/deform,ComputeTempDeform)
ComputeStyle(temp/partial,ComputeTempPartial)
ComputeStyle(temp/ramp,ComputeTempRamp)
ComputeStyle(temp/region,ComputeTempRegion)
ComputeStyle(variable,ComputeVariable)
ComputeStyle(variable/atom,ComputeVariableAtom)
#endif

#ifdef DihedralInclude
#endif

#ifdef DihedralClass
#endif

#ifdef DumpInclude
#include "dump_atom.h"
#include "dump_custom.h"
#include "dump_dcd.h"
#include "dump_xyz.h"
#endif

#ifdef DumpClass
DumpStyle(atom,DumpAtom)
DumpStyle(custom,DumpCustom)
DumpStyle(dcd,DumpDCD)
DumpStyle(xyz,DumpXYZ)
#endif

#ifdef FixInclude
#include "fix_add_force.h"
#include "fix_ave_force.h"
#include "fix_ave_spatial.h"
#include "fix_ave_time.h"
#include "fix_com.h"
#include "fix_drag.h"
#include "fix_deform.h"
#include "fix_deposit.h"
#include "fix_efield.h"
#include "fix_enforce2d.h"
#include "fix_gravity.h"
#include "fix_gyration.h"
#include "fix_heat.h"
#include "fix_indent.h"
#include "fix_langevin.h"
#include "fix_line_force.h"
#include "fix_minimize.h"
#include "fix_msd.h"
#include "fix_momentum.h"
#include "fix_nph.h"
#include "fix_npt.h"
#include "fix_nve.h"
#include "fix_nve_noforce.h"
#include "fix_nvt.h"
#include "fix_nvt_sllod.h"
#include "fix_plane_force.h"
#include "fix_print.h"
#include "fix_orient_fcc.h"
#include "fix_rdf.h"
#include "fix_recenter.h"
#include "fix_respa.h"
#include "fix_rigid.h"
#include "fix_relaxator.h"   
#include "fix_set_force.h"
#include "fix_shake.h"
#include "fix_shear_history.h"
#include "fix_spring.h"
#include "fix_spring_rg.h"
#include "fix_spring_self.h"
#include "fix_temp_rescale.h"
#include "fix_tmd.h"
#include "fix_viscous.h"
#include "fix_wall_lj126.h"
#include "fix_wall_lj93.h"
#include "fix_wall_reflect.h"
#include "fix_wiggle.h"
#endif

#ifdef FixClass
FixStyle(addforce,FixAddForce)
FixStyle(aveforce,FixAveForce)
FixStyle(ave/spatial,FixAveSpatial)
FixStyle(ave/time,FixAveTime)
FixStyle(com,FixCOM)
FixStyle(drag,FixDrag)
FixStyle(deform,FixDeform)
FixStyle(deposit,FixDeposit)
FixStyle(efield,FixEfield)
FixStyle(enforce2d,FixEnforce2D)
FixStyle(gravity,FixGravity)
FixStyle(gyration,FixGyration)
FixStyle(heat,FixHeat)
FixStyle(indent,FixIndent)
FixStyle(langevin,FixLangevin)
FixStyle(lineforce,FixLineForce)
FixStyle(MINIMIZE,FixMinimize)
FixStyle(momentum,FixMomentum)
FixStyle(msd,FixMSD)
FixStyle(nph,FixNPH)
FixStyle(npt,FixNPT)
FixStyle(nve,FixNVE)
FixStyle(nve,FixNVENoforce)
FixStyle(nvt,FixNVT)
FixStyle(nvt/sllod,FixNVTSlodd)
FixStyle(orient/fcc,FixOrientFCC)
FixStyle(print,FixPrint)
FixStyle(planeforce,FixPlaneForce)
FixStyle(rdf,FixRDF)
FixStyle(recenter,FixRecenter)
FixStyle(RESPA,FixRespa)
FixStyle(rigid,FixRigid)
FixStyle(setforce,FixSetForce)
FixStyle(relaxator,FixRelaxator)   
FixStyle(shake,FixShake)
FixStyle(SHEAR_HISTORY,FixShearHistory)
FixStyle(spring,FixSpring)
FixStyle(spring/rg,FixSpringRG)
FixStyle(spring/self,FixSpringSelf)
FixStyle(temp/rescale,FixTempRescale)
FixStyle(tmd,FixTMD)
FixStyle(viscous,FixViscous)
FixStyle(wall/lj126,FixWallLJ126)
FixStyle(wall/lj93,FixWallLJ93)
FixStyle(wall/reflect,FixWallReflect)
FixStyle(wiggle,FixWiggle)
#endif

#ifdef ImproperInclude
#endif

#ifdef ImproperClass
#endif

#ifdef IntegrateInclude
#include "respa.h"
#include "verlet.h"
#endif

#ifdef IntegrateClass
IntegrateStyle(respa,Respa)
IntegrateStyle(verlet,Verlet)
# endif

#ifdef KSpaceInclude
#endif

#ifdef KSpaceClass
#endif

#ifdef MinimizeInclude
#include "min_cg.h"
#include "min_cg_fr.h"
#include "min_sd.h"
#endif

#ifdef MinimizeClass
MinimizeStyle(cg,MinCG)
MinimizeStyle(cg/fr,MinCGFR)
MinimizeStyle(sd,MinSD)
# endif

#ifdef PairInclude
#include "pair_buck.h"
#include "pair_buck_coul_cut.h"
#include "pair_hybrid.h"
#include "pair_lj_cut.h"
#include "pair_lj_cut_coul_cut.h"
#include "pair_lj_cut_coul_debye.h"
#include "pair_lj_expand.h"
#include "pair_lj_smooth.h"
#include "pair_morse.h"
#include "pair_soft.h"
#include "pair_table.h"
#include "pair_yukawa.h"
#endif

#ifdef PairClass
PairStyle(buck,PairBuck)
PairStyle(buck/coul/cut,PairBuckCoulCut)
PairStyle(hybrid,PairHybrid)
PairStyle(lj/cut,PairLJCut)
PairStyle(lj/cut/coul/cut,PairLJCutCoulCut)
PairStyle(lj/cut/coul/debye,PairLJCutCoulDebye)
PairStyle(lj/expand,PairLJExpand)
PairStyle(lj/smooth,PairLJSmooth)
PairStyle(morse,PairMorse)
PairStyle(soft,PairSoft)
PairStyle(table,PairTable)
PairStyle(yukawa,PairYukawa)
#endif

#ifdef RegionInclude
#include "region_block.h"
#include "region_cylinder.h"
#include "region_intersect.h"
#include "region_prism.h"
#include "region_sphere.h"
#include "region_union.h"
#endif

#ifdef RegionClass
RegionStyle(block,RegBlock)
RegionStyle(cylinder,RegCylinder)
RegionStyle(intersect,RegIntersect)
RegionStyle(prism,RegPrism)
RegionStyle(sphere,RegSphere)
RegionStyle(union,RegUnion)
#endif

// style files for optional packages

#include "style_asphere.h"
#include "style_class2.h"
#include "style_colloid.h"
#include "style_dipole.h"
#include "style_dpd.h"
#include "style_granular.h"
#include "style_kspace.h"
#include "style_manybody.h"
#include "style_meam.h"
#include "style_molecule.h"
#include "style_opt.h"
#include "style_poems.h"
#include "style_xtc.h"

// user add-ons

#include "style_user.h"
