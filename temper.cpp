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
   Contributing author: Mark Sears (SNL)
------------------------------------------------------------------------- */

#include "math.h"
#include "stdlib.h"
#include "string.h"
#include "temper.h"
#include "universe.h"
#include "domain.h"
#include "update.h"
#include "integrate.h"
#include "modify.h"
#include "force.h"
#include "output.h"
#include "thermo.h"
#include "fix.h"
#include "fix_nvt.h"
#include "fix_langevin.h"
#include "random_park.h"
#include "finish.h"
#include "timer.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;

#define NVT      1
#define LANGEVIN 2
// #define TEMPER_DEBUG 1
/* ---------------------------------------------------------------------- */

Temper::Temper(LAMMPS *lmp) : Pointers(lmp) {}

/* ---------------------------------------------------------------------- */

Temper::~Temper()
{
  MPI_Comm_free(&roots);
  if (ranswap) delete ranswap;
  delete ranboltz;
  delete [] set_temp;
  delete [] temp2world;
  delete [] world2temp;
  delete [] world2root;
}

/* ----------------------------------------------------------------------
   perform tempering with inter-world swaps
------------------------------------------------------------------------- */

void Temper::command(int narg, char **arg)
{
  if (universe->nworlds == 1) 
    error->all("Must have more than one processor partition to temper");

  if (narg != 6 && narg != 7) error->universe_all("Illegal temper command");

  if (domain->box_exist == 0) 
    error->all("Temper command before simulation box is defined");

  update->nsteps = atoi(arg[0]);
  update->beginstep = update->firststep = update->ntimestep;
  update->endstep = update->laststep = update->firststep + update->nsteps;
  update->whichflag = 0;

  lmp->init();
  
  // grab temper command args

  int nsteps = update->nsteps;
  nevery = atoi(arg[1]);
  double temp = atof(arg[2]);

  for (whichfix = 0; whichfix < modify->nfix; whichfix++)
    if (strcmp(arg[3],modify->fix[whichfix]->id) == 0) break;
  if (whichfix == modify->nfix) 
    error->universe_all("Tempering fix ID is not defined");

  seed_swap = atoi(arg[4]);
  seed_boltz = atoi(arg[5]);

  my_set_temp = universe->iworld;
  if (narg == 7) my_set_temp = atoi(arg[6]);

  // swap frequency must evenly divide total # of timesteps

  nswaps = nsteps/nevery;
  if (nswaps*nevery != nsteps) 
    error->universe_all("Non integer # of swaps in temper command");

  // thermodynamics must be computed on swap steps
  // potential energy must be computed by thermo

  if (nevery % output->thermo_every)
    error->universe_all("Thermodynamics not computed on tempering swap steps");

  if (output->thermo->peflag == 0)
    error->universe_all("Thermodynamics must compute PE for temper");

  // fix style must be appropriate for temperature control

  if (strcmp(modify->fix[whichfix]->style,"nvt") == 0) fixstyle = NVT;
  else if (strcmp(modify->fix[whichfix]->style,"langevin") == 0) 
    fixstyle = LANGEVIN;
  else error->universe_all("Tempering fix is not valid");

  // local storage

  me_universe = universe->me;
  MPI_Comm_rank(world,&me);
  nworlds = universe->nworlds;
  iworld = universe->iworld;
  boltz = force->boltz;

  // create MPI communicator for root proc from each world

  int color;
  if (me == 0) color = 0;
  else color = 1;
  MPI_Comm_split(universe->uworld,color,0,&roots);

  // RNGs for swaps and Boltzmann test
  // warm up Boltzmann RNG

  if (seed_swap) ranswap = new RanPark(lmp,seed_swap);
  else ranswap = NULL;
  ranboltz = new RanPark(lmp,seed_boltz + me_universe);
  for (int i = 0; i < 100; i++) ranboltz->uniform();

  // world2root[i] = global proc that is root proc of world i

  world2root = new int[nworlds];
  if (me == 0) 
    MPI_Allgather(&me_universe,1,MPI_INT,world2root,1,MPI_INT,roots);
  MPI_Bcast(world2root,nworlds,MPI_INT,0,world);

  // create static list of set temperatures
  // allgather tempering arg "temp" across root procs
  // bcast from each root to other procs in world

  set_temp = new double[nworlds];
  if (me == 0) MPI_Allgather(&temp,1,MPI_DOUBLE,set_temp,1,MPI_DOUBLE,roots);
  MPI_Bcast(set_temp,nworlds,MPI_DOUBLE,0,world);

  // create world2temp only on root procs from my_set_temp
  // create temp2world on root procs from world2temp,
  //   then bcast to all procs within world

  world2temp = new int[nworlds];
  temp2world = new int[nworlds];
  if (me == 0) {
    MPI_Allgather(&my_set_temp,1,MPI_INT,world2temp,1,MPI_INT,roots);
    for (int i = 0; i < nworlds; i++) temp2world[world2temp[i]] = i;
  }
  MPI_Bcast(temp2world,nworlds,MPI_INT,0,world);

  // if restarting tempering, reset temp target of Fix to current my_set_temp

  if (narg == 7) {
    double new_temp = set_temp[my_set_temp];
    if (fixstyle == NVT) 
      ((FixNVT *) modify->fix[whichfix])->reset_target(new_temp);
    else if (fixstyle == LANGEVIN)
      ((FixLangevin *) modify->fix[whichfix])->reset_target(new_temp);
  }

  // setup tempering runs

  int i,which,partner,swap,partner_set_temp,partner_world;
  double pe,pe_partner,boltz_factor,new_temp;
  MPI_Status status;

  if (me_universe == 0 && universe->uscreen) 
    fprintf(universe->uscreen,"Setting up tempering ...\n");

  update->integrate->setup();

  timer->barrier_start(TIME_LOOP);

  if (me_universe == 0) {
    if (universe->uscreen) fprintf(universe->uscreen,"Step T1 T2 ...\n");
    if (universe->ulogfile) fprintf(universe->ulogfile,"Step T1 T2 ...\n");
    print_status();
  }

  for (int iswap = 0; iswap < nswaps; iswap++) {

    // run for nevery timesteps

    update->integrate->iterate(nevery);

    // which = which of 2 kinds of swaps to do (0,1)

    if (!ranswap) which = iswap % 2;
    else if (ranswap->uniform() < 0.5) which = 0;
    else which = 1;

    // partner_set_temp = which set temp I am partnering with for this swap

    if (which == 0) {
      if (my_set_temp % 2 == 0) partner_set_temp = my_set_temp + 1;
      else partner_set_temp = my_set_temp - 1;
    } else {
      if (my_set_temp % 2 == 1) partner_set_temp = my_set_temp + 1;
      else partner_set_temp = my_set_temp - 1;
    }

    // partner = proc ID to swap with
    // if partner = -1, then I am not a proc that swaps

    partner = -1;
    if (me == 0 && partner_set_temp >= 0 && partner_set_temp < nworlds) {
      partner_world = temp2world[partner_set_temp];
      partner = world2root[partner_world];
    }

    // swap with a partner, only root procs in each world participate
    // hi proc sends PE to low proc
    // lo proc make Boltzmann decision on whether to swap
    // lo proc communicates decision back to hi proc

    swap = 0;
    if (partner != -1) {
      pe = output->thermo->potential_energy;
      if (me_universe > partner) 
	MPI_Send(&pe,1,MPI_DOUBLE,partner,0,universe->uworld);
      else
	MPI_Recv(&pe_partner,1,MPI_DOUBLE,partner,0,universe->uworld,&status);

      if (me_universe < partner) {
	boltz_factor = (pe - pe_partner) * 
	  (1.0/(boltz*set_temp[my_set_temp]) - 
	   1.0/(boltz*set_temp[partner_set_temp]));
	if (boltz_factor >= 0.0) swap = 1;
	else if (ranboltz->uniform() < exp(boltz_factor)) swap = 1;
      }

      if (me_universe < partner) 
	MPI_Send(&swap,1,MPI_INT,partner,0,universe->uworld);
      else
	MPI_Recv(&swap,1,MPI_INT,partner,0,universe->uworld,&status);

#ifdef TEMPER_DEBUG
      if (me_universe < partner)
	printf("SWAP %d & %d: yes = %d,Ts = %d %d, PEs = %g %g, Bz = %g %g\n",
	       me_universe,partner,swap,my_set_temp,partner_set_temp,
	       pe,pe_partner,boltz_factor,exp(boltz_factor));
#endif

    }

    // bcast swap result to other procs in my world

    MPI_Bcast(&swap,1,MPI_INT,0,world);

    // if my world swapped, all procs in world reset temp target of Fix

    if (swap) {
      new_temp = set_temp[partner_set_temp];
      if (fixstyle == NVT) 
	((FixNVT *) modify->fix[whichfix])->reset_target(new_temp);
      else if (fixstyle == LANGEVIN)
	((FixLangevin *) modify->fix[whichfix])->reset_target(new_temp);
    }

    // update my_set_temp and temp2world on every proc
    // root procs update their value if swap took place
    // allgather across root procs
    // bcast within my world

    if (swap) my_set_temp = partner_set_temp;
    if (me == 0) {
      MPI_Allgather(&my_set_temp,1,MPI_INT,world2temp,1,MPI_INT,roots);
      for (i = 0; i < nworlds; i++) temp2world[world2temp[i]] = i;
    }
    MPI_Bcast(temp2world,nworlds,MPI_INT,0,world);

    // print out current swap status

    if (me_universe == 0) print_status();
  }

  timer->barrier_stop(TIME_LOOP);

  Finish finish(lmp);
  finish.end(1);
  update->whichflag = -1;
}

/* ----------------------------------------------------------------------
   proc 0 prints current tempering status
------------------------------------------------------------------------- */

void Temper::print_status()
{
  if (universe->uscreen) {
    fprintf(universe->uscreen,"%d ",update->ntimestep);
    for (int i = 0; i < nworlds; i++) 
      fprintf(universe->uscreen,"%d ",world2temp[i]);
    fprintf(universe->uscreen,"\n");
  }
  if (universe->ulogfile) {
    fprintf(universe->ulogfile,"%d ",update->ntimestep);
    for (int i = 0; i < nworlds; i++)
      fprintf(universe->ulogfile,"%d ",world2temp[i]);
    fprintf(universe->ulogfile,"\n");
  }
}
