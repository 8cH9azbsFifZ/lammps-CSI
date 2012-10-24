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

#include "mpi.h"
#include "stdlib.h"
#include "string.h"
#include "thermo.h"
#include "atom.h"
#include "comm.h"
#include "domain.h"
#include "update.h"
#include "modify.h"
#include "fix.h"
#include "compute.h"
#include "input.h"
#include "variable.h"
#include "force.h"
#include "pair.h"
#include "bond.h"
#include "angle.h"
#include "dihedral.h"
#include "improper.h"
#include "kspace.h"
#include "timer.h"
#include "memory.h"
#include "error.h"
#include "fix_indent.h" // for indenter properties

using namespace LAMMPS_NS;

// customize a new keyword by adding to this list:

// step, atoms, cpu, temp, press, pe, ke, etotal, enthalpy
// evdwl, ecoul, epair, ebond, eangle, edihed, eimp, emol, elong, etail
// vol, lx, ly, lz, xlo, xhi, ylo, yhi, zlo, zhi
// pxx, pyy, pzz, pxy, pxz, pyz
// drot, grot (rotational KE for dipole and granular particles)
// tave, pave, eave, peave (time-averaged quantities)

// customize a new thermo style by adding a DEFINE to this list

#define ONE "step temp epair emol etotal press"
#define MULTI "etotal ke temp pe ebond eangle edihed eimp evdwl ecoul elong press"
#define GRANULAR "step atoms ke grot"
#define DIPOLE "step temp epair drot etotal press"

enum{IGNORE,WARN,ERROR};           // same as write_restart.cpp
enum{ONELINE,MULTILINE};
enum{INT,FLOAT};

#define MAXLINE 1024
#define DELTA   8

/* ---------------------------------------------------------------------- */

Thermo::Thermo(LAMMPS *lmp, int narg, char **arg) : Pointers(lmp)
{
  MPI_Comm_rank(world,&me);

  int n = strlen(arg[0]) + 1;
  style = new char[n];
  strcpy(style,arg[0]);

  thermoflag = 1;

  // set thermo_modify defaults

  normuserflag = 0;
  lineflag = ONELINE;
  lostflag = ERROR;
  lostbefore = 0;
  flushflag = 0;
  nwindow = 10;
  ncurrent_t = ncurrent_p = ncurrent_e = ncurrent_pe = -1;
  npartial_t = npartial_p = npartial_e = npartial_pe = 0;

  // set style and corresponding lineflag
  // custom style builds its own line of keywords
  // customize a new thermo style by adding to if statement

  line = new char[MAXLINE];

  if (strcmp(style,"one") == 0) {
    strcpy(line,ONE);
  } else if (strcmp(style,"multi") == 0) {
    strcpy(line,MULTI);
    lineflag = MULTILINE;
  } else if (strcmp(style,"granular") == 0) {
    strcpy(line,GRANULAR);
  } else if (strcmp(style,"dipole") == 0) {
    strcpy(line,DIPOLE);

  } else if (strcmp(style,"custom") == 0) {
    if (narg == 1) error->all("Illegal thermo style custom command");
    line[0] = '\0';
    for (int iarg = 1; iarg < narg; iarg++) {
      strcat(line,arg[iarg]);
      strcat(line," ");
    }
    line[strlen(line)-1] = '\0';

  } else error->all("Illegal thermo style command");

  // ptrs, flags, IDs for compute objects thermo may use or create

  temperature = NULL;
  pressure = NULL;
  rotate_dipole = NULL;
  rotate_gran = NULL;

  index_temp = index_press = index_drot = index_grot = -1;
  internal_drot = internal_grot = 0;

  id_temp = "thermo_temp";
  id_press = "thermo_pressure";
  id_drot = "thermo_rotate_dipole";
  id_grot = "thermo_rotate_gran";

  // count fields in line
  // allocate per-field memory
  // process line of keywords

  nfield_initial = atom->count_words(line);
  allocate();
  parse_fields(line);

  // create the requested compute styles
  // temperature and pressure always exist b/c Output class created them

  if (index_drot >= 0) {
    create_compute(id_drot,"rotate/dipole",NULL);
    internal_drot = 1;
  }
  if (index_grot >= 0) {
    create_compute(id_grot,"rotate/gran",NULL);
    internal_grot = 1;
  }

  // format strings

  format_multi = "---------------- Step %8d ----- "
                 "CPU = %11.4f (sec) ----------------";
  format_int_one_def = "%8d";
  format_int_multi_def = "%14d";
  format_g_def = "%12.8g";
  format_f_def = "%14.4f";
  format_int_user = NULL;
  format_float_user = NULL;

  // average quantities

  tsum = psum = esum = pesum = 0.0;
  tpast = new double[nwindow];
  ppast = new double[nwindow];
  epast = new double[nwindow];
  pepast = new double[nwindow];

  // prepare indenter properties
  IndexFixIndentBad = 666;
  FindIndexFixIndent();

  IndexComputeAcklandBad = 666;
  FindIndexComputeAckland();
}

/* ---------------------------------------------------------------------- */

// Find Index of the Fix for the Indenter
int Thermo::FindIndexFixIndent()
{
   IndexFixIndent = IndexFixIndentBad;
   int i = 0;
   for (i = 0; i < modify->nfix; i++) 
      if (strcmp("indent",modify->fix[i]->style) == 0)
         IndexFixIndent = i;
   if (IndexFixIndent == IndexFixIndentBad) 
      error->warning ("FindIndexFixIndent: No fix _indent_ found");
   return i;
}

/* ---------------------------------------------------------------------- */

// Find Index of the Ackland Computation
int Thermo::FindIndexComputeAckland()
{
   IndexComputeAckland = IndexComputeAcklandBad;
   int i = 0;
   for (i = 0; i < modify->ncompute; i++)
      if (strcmp(modify->compute[i]->style,"ackland/atom") == 0)
         IndexComputeAckland = i;
   if (IndexComputeAckland == IndexComputeAcklandBad)
      error->warning ("FindIndexComputeAckland: No fix _ackland/atom_ found");
   return i;
}   

/* ---------------------------------------------------------------------- */

Thermo::~Thermo()
{
  delete [] style;
  delete [] line;

  // delete Compute classes if thermo created them

  if (index_drot >= 0 && internal_drot)
    modify->delete_compute(id_compute[index_drot]);
  if (index_grot >= 0 && internal_grot)
    modify->delete_compute(id_compute[index_grot]);

  deallocate();

  // format strings

  delete [] format_int_user;
  delete [] format_float_user;

  // average arrays

  delete [] tpast;
  delete [] ppast;
  delete [] epast;
  delete [] pepast;
}

/* ---------------------------------------------------------------------- */

void Thermo::init()
{
  int i,n;

  // set normvalue to default setting unless user has specified it

  if (normuserflag) normvalue = normuser;
  else if (strcmp(style,"granular") == 0) normvalue = 0;
  else if (strcmp(update->unit_style,"lj") == 0) normvalue = 1;
  else normvalue = 0;

  // add Volume field if volume changes and not style = custom
  // this check must come after domain init, so box_change is set

  nfield = nfield_initial;
  if (domain->box_change && strcmp(style,"custom") != 0)
    addfield("Volume",&Thermo::compute_vol,FLOAT);

  // set format string for each field
  // include keyword if lineflag = MULTILINE
  // add '/n' every 3 values if lineflag = MULTILINE
  // add trailing '/n' to last value

  char *ptr;
  for (i = 0; i < nfield; i++) {
    format[i][0] = '\0';
    if (lineflag == MULTILINE && i % 3 == 0) strcat(format[i],"\n");

    if (format_user[i]) ptr = format_user[i];
    else if (vtype[i] == INT && format_int_user) ptr = format_int_user;
    else if (vtype[i] == INT && lineflag == ONELINE) ptr = format_int_one_def;
    else if (vtype[i] == INT && lineflag == MULTILINE) ptr = format_int_multi_def;
    else if (vtype[i] == FLOAT && format_float_user) ptr = format_float_user;
    else if (lineflag == ONELINE) ptr = format_g_def;
    else if (lineflag == MULTILINE) ptr = format_f_def;

    n = strlen(format[i]);
    if (lineflag == ONELINE) sprintf(&format[i][n],"%s ",ptr);
    else sprintf(&format[i][n],"%-8s = %s ",keyword[i],ptr);

    if (i == nfield-1) strcat(format[i],"\n");
  }

  // find current ptr for each Compute ID

  int icompute;
  for (i = 0; i < ncompute; i++) {
    icompute = modify->find_compute(id_compute[i]);
    if (icompute < 0) error->all("Could not find thermo compute ID");
    computes[i] = modify->compute[icompute];
  }

  // find current ptr for each Fix ID

  int ifix;
  for (i = 0; i < nfix; i++) {
    ifix = modify->find_fix(id_fix[i]);
    if (ifix < 0) error->all("Could not find thermo fix ID");
    fixes[i] = modify->fix[ifix];
  }

  // set ptrs to keyword-specific Compute objects

  if (index_temp >= 0) temperature = computes[index_temp];
  if (index_press >= 0) pressure = computes[index_press];
  if (index_drot >= 0) rotate_dipole = computes[index_drot];
  if (index_grot >= 0) rotate_gran = computes[index_grot];
}

/* ---------------------------------------------------------------------- */

void Thermo::header()
{
  if (lineflag == MULTILINE) return;

  int loc = 0;
  for (int i = 0; i < nfield; i++)
    loc += sprintf(&line[loc],"%s ",keyword[i]);
  sprintf(&line[loc],"\n");
  
  if (me == 0) {
    if (screen) fprintf(screen,line);
    if (logfile) fprintf(logfile,line);
  }
}

/* ---------------------------------------------------------------------- */

void Thermo::compute(int flag)
{
  int i;
  double tmp;

  firststep = flag;

  // check for lost atoms
  // turn off normflag if natoms = 0 to avoid divide by 0

  natoms = lost_check();
  if (natoms == 0) normflag = 0;
  else normflag = normvalue;

  // invoke Compute methods needed for thermo keywords
  // call compute_scalar() if which = 0, else compute_vector()

  for (i = 0; i < ncompute; i++) {
    if (compute_which[i] % 2 == 0) tmp = computes[i]->compute_scalar();
    if (compute_which[i] > 0) computes[i]->compute_vector();
  }

  // if lineflag = MULTILINE, prepend step/cpu header line

  int loc = 0;
  if (lineflag == MULTILINE) {
    double cpu;
    if (flag) cpu = timer->elapsed(TIME_LOOP);
    else cpu = 0.0;
    loc = sprintf(&line[loc],format_multi,update->ntimestep,cpu);
  }

  // add each thermo value to line with its specific format

  for (ifield = 0; ifield < nfield; ifield++) {
    (this->*vfunc[ifield])();
    if (vtype[ifield] == INT) loc += sprintf(&line[loc],format[ifield],ivalue);
    else loc += sprintf(&line[loc],format[ifield],dvalue);
  }

  // print line to screen and logfile

  if (me == 0) {
    if (screen) fprintf(screen,line);
    if (logfile) {
      fprintf(logfile,line);
      if (flushflag) fflush(logfile);
    }
  }
}

/* ----------------------------------------------------------------------
   check for lost atoms, return current number of atoms
------------------------------------------------------------------------- */

double Thermo::lost_check()
{
  // ntotal = current # of atoms

  double ntotal;
  double rlocal = atom->nlocal;
  MPI_Allreduce(&rlocal,&ntotal,1,MPI_DOUBLE,MPI_SUM,world);
  if (ntotal == atom->natoms) return ntotal;

  // if not checking or already warned, just return

  if (lostflag == IGNORE) return ntotal;
  if (lostflag == WARN && lostbefore == 1) return ntotal;

  // error message

  if (lostflag == ERROR) {
    char str[128];
    sprintf(str,"Lost atoms: original %.15g current %.15g",
	    atom->natoms,ntotal);
    error->all(str);
  }

  // warning message

  char str[128];
  sprintf(str,"Lost atoms: original %.15g current %.15g",atom->natoms,ntotal);
  if (me == 0) error->warning(str);
  lostbefore = 1;
  return ntotal;
}

/* ----------------------------------------------------------------------
   modify thermo parameters
------------------------------------------------------------------------- */

void Thermo::modify_params(int narg, char **arg)
{
  if (narg == 0) error->all("Illegal thermo_modify command");

  int iarg = 0;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"temp") == 0) {
      if (iarg+2 > narg) error->all("Illegal thermo_modify command");
      if (index_temp < 0) error->all("Thermo style does not use temp");
      delete [] id_compute[index_temp];
      int n = strlen(arg[iarg+1]) + 1;
      id_compute[index_temp] = new char[n];
      strcpy(id_compute[index_temp],arg[iarg+1]);

      int icompute = modify->find_compute(id_compute[index_temp]);
      if (icompute < 0) error->all("Could not find thermo_modify temp ID");
      temperature = modify->compute[icompute];

      if (temperature->tempflag == 0)
	error->all("Thermo_modify temp ID does not compute temperature");
      if (temperature->igroup != 0 && comm->me == 0)
	error->warning("Temperature for thermo pressure is not for group all");

      // reset id_pre of pressure to new temp ID
      // either pressure currently being used by thermo or "thermo_pressure"

      if (index_press >= 0) {
	icompute = modify->find_compute(id_compute[index_press]);
	if (icompute < 0) error->all("Press ID for thermo does not exist");
      } else icompute = modify->find_compute("thermo_pressure");

      delete [] modify->compute[icompute]->id_pre;
      modify->compute[icompute]->id_pre = new char[n];
      strcpy(modify->compute[icompute]->id_pre,arg[iarg+1]);

      iarg += 2;

    } else if (strcmp(arg[iarg],"press") == 0) {
      if (iarg+2 > narg) error->all("Illegal thermo_modify command");
      if (index_press < 0) error->all("Thermo style does not use press");
      delete [] id_compute[index_press];
      int n = strlen(arg[iarg+1]) + 1;
      id_compute[index_press] = new char[n];
      strcpy(id_compute[index_press],arg[iarg+1]);

      int icompute = modify->find_compute(id_compute[index_press]);
      if (icompute < 0) error->all("Could not find thermo_modify press ID");
      pressure = modify->compute[icompute];

      if (pressure->pressflag == 0)
	error->all("Thermo_modify press ID does not compute pressure");

      // if id_pre of new pressure is not being computed, add to compute list
      // swap it with pressure in list so id_pre will be computed first
      // OK to call add_compute with "which" acting as index

      int which = compute_which[index_press];
      int ncompute_current = ncompute;
      icompute = add_compute(pressure->id_pre,which);
      if (icompute == ncompute_current) {
	int iswap = compute_which[index_press];
	compute_which[index_press] = compute_which[icompute];
	compute_which[icompute] = iswap;
	char *cswap = id_compute[index_press];
	id_compute[index_press] = id_compute[icompute];
	id_compute[icompute] = cswap;
	index_press = icompute;
      }

      iarg += 2;

    } else if (strcmp(arg[iarg],"drot") == 0) {
      if (iarg+2 > narg) error->all("Illegal thermo_modify command");
      if (index_drot < 0) error->all("Thermo style does not use drot");
      if (internal_drot) {
	modify->delete_compute(id_compute[index_drot]);
	internal_drot = 0;
      }
      delete [] id_compute[index_drot];
      int n = strlen(arg[iarg+1]) + 1;
      id_compute[index_drot] = new char[n];
      strcpy(id_compute[index_drot],arg[iarg+1]);

      int icompute = modify->find_compute(id_compute[index_drot]);
      if (icompute < 0) error->all("Could not find thermo_modify drot ID");
      iarg += 2;

    } else if (strcmp(arg[iarg],"grot") == 0) {
      if (iarg+2 > narg) error->all("Illegal thermo_modify command");
      if (index_grot < 0) error->all("Thermo style does not use grot");
      if (internal_grot) {
	modify->delete_compute(id_compute[index_grot]);
	internal_grot = 0;
      }
      delete [] id_compute[index_grot];
      int n = strlen(arg[iarg+1]) + 1;
      id_compute[index_grot] = new char[n];
      strcpy(id_compute[index_grot],arg[iarg+1]);

      int icompute = modify->find_compute(id_compute[index_grot]);
      if (icompute < 0) error->all("Could not find thermo_modify grot ID");
      iarg += 2;

    } else if (strcmp(arg[iarg],"lost") == 0) {
      if (iarg+2 > narg) error->all("Illegal thermo_modify command");
      if (strcmp(arg[iarg+1],"ignore") == 0) lostflag = IGNORE;
      else if (strcmp(arg[iarg+1],"warn") == 0) lostflag = WARN;
      else if (strcmp(arg[iarg+1],"error") == 0) lostflag = ERROR;
      else error->all("Illegal thermo_modify command");
      iarg += 2;

    } else if (strcmp(arg[iarg],"norm") == 0) {
      if (iarg+2 > narg) error->all("Illegal thermo_modify command");
      normuserflag = 1;
      if (strcmp(arg[iarg+1],"no") == 0) normuser = 0;
      else if (strcmp(arg[iarg+1],"yes") == 0) normuser = 1;
      else error->all("Illegal thermo_modify command");
      iarg += 2;

    } else if (strcmp(arg[iarg],"flush") == 0) {
      if (iarg+2 > narg) error->all("Illegal thermo_modify command");
      if (strcmp(arg[iarg+1],"no") == 0) flushflag = 0;
      else if (strcmp(arg[iarg+1],"yes") == 0) flushflag = 1;
      else error->all("Illegal thermo_modify command");
      iarg += 2;

    } else if (strcmp(arg[iarg],"line") == 0) {
      if (iarg+2 > narg) error->all("Illegal thermo_modify command");
      if (strcmp(arg[iarg+1],"one") == 0) lineflag = ONELINE;
      else if (strcmp(arg[iarg+1],"multi") == 0) lineflag = MULTILINE;
      else error->all("Illegal thermo_modify command");
      iarg += 2;

    } else if (strcmp(arg[iarg],"format") == 0) {
      if (iarg+3 > narg) error->all("Illegal thermo_modify command");
      if (strcmp(arg[iarg+1],"int") == 0) {
	if (format_int_user) delete [] format_int_user;
	int n = strlen(arg[iarg+2]) + 1;
	format_int_user = new char[n];
	strcpy(format_int_user,arg[iarg+2]);
      } else if (strcmp(arg[iarg+1],"float") == 0) {
	if (format_float_user) delete [] format_float_user;
	int n = strlen(arg[iarg+2]) + 1;
	format_float_user = new char[n];
	strcpy(format_float_user,arg[iarg+2]);
      } else {
	int i = atoi(arg[iarg+1]) - 1;
	if (i < 0 || i >= nfield_initial)
	  error->all("Illegal thermo_modify command");
	if (format_user[i]) delete [] format_user[i];
	int n = strlen(arg[iarg+2]) + 1;
	format_user[i] = new char[n];
	strcpy(format_user[i],arg[iarg+2]);
      }
      iarg += 3;

    } else if (strcmp(arg[iarg],"window") == 0) {
      if (iarg+2 > narg) error->all("Illegal thermo_modify command");
      nwindow = atoi(arg[iarg+1]);
      if (nwindow <= 0) error->all("Illegal thermo_modify command");
      ncurrent_t = ncurrent_p = ncurrent_e = ncurrent_pe = -1;
      npartial_t = npartial_p = npartial_e = npartial_pe = 0;
      tsum = psum  = esum = pesum = 0.0;
      delete [] tpast;
      delete [] ppast;
      delete [] epast;
      delete [] pepast;
      tpast = new double[nwindow];
      ppast = new double[nwindow];
      epast = new double[nwindow];
      pepast = new double[nwindow];
      iarg += 2;

    } else error->all("Illegal thermo_modify command");
  }
}

/* ----------------------------------------------------------------------
   allocate all per-field memory
   allow each c_ID to imply 2 Compute objects (if it has id_pre)
------------------------------------------------------------------------- */

void Thermo::allocate()
{
  // n = specified fields + Volume field (added at run time)

  int n = nfield_initial + 1;

  keyword = new char*[n];
  for (int i = 0; i < n; i++) keyword[i] = new char[32];
  vfunc = new FnPtr[n];
  vtype = new int[n];

  format = new char*[n];
  for (int i = 0; i < n; i++) format[i] = new char[32];
  format_user = new char*[n];
  for (int i = 0; i < n; i++) format_user[i] = NULL;

  field2object = new int[n];
  arg_object = new int[n];

  ncompute = 0;
  id_compute = new char*[2*n];
  compute_which = new int[2*n];
  computes = new Compute*[2*n];

  nfix = 0;
  id_fix = new char*[n];
  fixes = new Fix*[n];

  nvariable = 0;
  id_variable = new char*[n];
}

/* ----------------------------------------------------------------------
   deallocate all per-field memory
------------------------------------------------------------------------- */

void Thermo::deallocate()
{
  int n = nfield_initial + 1;

  for (int i = 0; i < n; i++) delete [] keyword[i];
  delete [] keyword;
  delete [] vfunc;
  delete [] vtype;

  for (int i = 0; i < n; i++) delete [] format[i];
  delete [] format;
  for (int i = 0; i < n; i++) delete [] format_user[i];
  delete [] format_user;

  delete [] field2object;
  delete [] arg_object;

  for (int i = 0; i < ncompute; i++) delete [] id_compute[i];
  delete [] id_compute;
  delete [] compute_which;
  delete [] computes;

  for (int i = 0; i < nfix; i++) delete [] id_fix[i];
  delete [] id_fix;
  delete [] fixes;

  for (int i = 0; i < nvariable; i++) delete [] id_variable[i];
  delete [] id_variable;
}

/* ----------------------------------------------------------------------
   parse list of thermo keywords from str
   set compute flags (temp, press, etc)
------------------------------------------------------------------------- */

void Thermo::parse_fields(char *str)
{
  nfield = 0;
  peflag = 0;

  // customize a new keyword by adding to if statement

  char *word = strtok(str," \0");
  while (word) {

    if (strcmp(word,"step") == 0) {
      addfield("Step",&Thermo::compute_step,INT);
    } else if (strcmp(word,"atoms") == 0) {
      addfield("Atoms",&Thermo::compute_atoms,INT);
    } else if (strcmp(word,"cpu") == 0) {
      addfield("CPU",&Thermo::compute_cpu,FLOAT);

    } else if (strcmp(word,"temp") == 0) {
      addfield("Temp",&Thermo::compute_temp,FLOAT);
      index_temp = add_compute(id_temp,0);
    } else if (strcmp(word,"press") == 0) {
      addfield("Press",&Thermo::compute_press,FLOAT);
      index_temp = add_compute(id_temp,0);
      index_press = add_compute(id_press,0);
    } else if (strcmp(word,"pe") == 0) {
      addfield("PotEng",&Thermo::compute_pe,FLOAT);
      peflag = 1;
    } else if (strcmp(word,"ke") == 0) {
      addfield("KinEng",&Thermo::compute_ke,FLOAT);
      index_temp = add_compute(id_temp,0);
    } else if (strcmp(word,"etotal") == 0) {
      addfield("TotEng",&Thermo::compute_etotal,FLOAT);
      peflag = 1;
      index_temp = add_compute(id_temp,0);
    } else if (strcmp(word,"enthalpy") == 0) {
      addfield("Enthalpy",&Thermo::compute_enthalpy,FLOAT);
      index_temp = add_compute(id_temp,0);
      index_press = add_compute(id_press,0);

    } else if (strcmp(word,"evdwl") == 0) {
      addfield("E_vdwl",&Thermo::compute_evdwl,FLOAT);
    } else if (strcmp(word,"ecoul") == 0) {
      addfield("E_coul",&Thermo::compute_ecoul,FLOAT);
    } else if (strcmp(word,"epair") == 0) {
      addfield("E_pair",&Thermo::compute_epair,FLOAT);
    } else if (strcmp(word,"ebond") == 0) {
      addfield("E_bond",&Thermo::compute_ebond,FLOAT);
    } else if (strcmp(word,"eangle") == 0) {
      addfield("E_angle",&Thermo::compute_eangle,FLOAT);
    } else if (strcmp(word,"edihed") == 0) {
      addfield("E_dihed",&Thermo::compute_edihed,FLOAT);
    } else if (strcmp(word,"eimp") == 0) {
      addfield("E_impro",&Thermo::compute_eimp,FLOAT);
    } else if (strcmp(word,"emol") == 0) {
      addfield("E_mol",&Thermo::compute_emol,FLOAT);
    } else if (strcmp(word,"elong") == 0) {
      addfield("E_long",&Thermo::compute_elong,FLOAT);
    } else if (strcmp(word,"etail") == 0) {
      addfield("E_tail",&Thermo::compute_etail,FLOAT);

    } else if (strcmp(word,"vol") == 0) {
      addfield("Volume",&Thermo::compute_vol,FLOAT);
    } else if (strcmp(word,"lx") == 0) {
      addfield("Lx",&Thermo::compute_lx,FLOAT);
    } else if (strcmp(word,"ly") == 0) {
      addfield("Ly",&Thermo::compute_ly,FLOAT);
    } else if (strcmp(word,"lz") == 0) {
      addfield("Lz",&Thermo::compute_lz,FLOAT);

    } else if (strcmp(word,"xlo") == 0) {
      addfield("Xlo",&Thermo::compute_xlo,FLOAT);
    } else if (strcmp(word,"xhi") == 0) {
      addfield("Xhi",&Thermo::compute_xhi,FLOAT);
    } else if (strcmp(word,"ylo") == 0) {
      addfield("Ylo",&Thermo::compute_ylo,FLOAT);
    } else if (strcmp(word,"yhi") == 0) {
      addfield("Yhi",&Thermo::compute_yhi,FLOAT);
    } else if (strcmp(word,"zlo") == 0) {
      addfield("Zlo",&Thermo::compute_zlo,FLOAT);
    } else if (strcmp(word,"zhi") == 0) {
      addfield("Zhi",&Thermo::compute_zhi,FLOAT);

    } else if (strcmp(word,"pxx") == 0) {
      addfield("Pxx",&Thermo::compute_pxx,FLOAT);
      index_temp = add_compute(id_temp,1);
      index_press = add_compute(id_press,1);
    } else if (strcmp(word,"pyy") == 0) {
      addfield("Pyy",&Thermo::compute_pyy,FLOAT);
      index_temp = add_compute(id_temp,1);
      index_press = add_compute(id_press,1);
    } else if (strcmp(word,"pzz") == 0) {
      addfield("Pzz",&Thermo::compute_pzz,FLOAT);
      index_temp = add_compute(id_temp,1);
      index_press = add_compute(id_press,1);
    } else if (strcmp(word,"pxy") == 0) {
      addfield("Pxy",&Thermo::compute_pxy,FLOAT);
      index_temp = add_compute(id_temp,1);
      index_press = add_compute(id_press,1);
    } else if (strcmp(word,"pxz") == 0) {
      addfield("Pxz",&Thermo::compute_pxz,FLOAT);
      index_temp = add_compute(id_temp,1);
      index_press = add_compute(id_press,1);
    } else if (strcmp(word,"pyz") == 0) {
      addfield("Pyz",&Thermo::compute_pyz,FLOAT);
      index_temp = add_compute(id_temp,1);
      index_press = add_compute(id_press,1);

    } else if (strcmp(word,"drot") == 0) {
      addfield("RotKEdip",&Thermo::compute_drot,FLOAT);
      index_drot = add_compute(id_drot,0);
    } else if (strcmp(word,"grot") == 0) {
      addfield("RotKEgrn",&Thermo::compute_grot,FLOAT);
      index_grot = add_compute(id_grot,0);

    } else if (strcmp(word,"tave") == 0) {
      addfield("T_ave",&Thermo::compute_tave,FLOAT);
      index_temp = add_compute(id_temp,0);
    } else if (strcmp(word,"pave") == 0) {
      addfield("P_ave",&Thermo::compute_pave,FLOAT);
      index_temp = add_compute(id_temp,0);
      index_press = add_compute(id_press,0);
    } else if (strcmp(word,"eave") == 0) {
      addfield("E_ave",&Thermo::compute_eave,FLOAT);
      peflag = 1;
      index_temp = add_compute(id_temp,0);
    } else if (strcmp(word,"peave") == 0) {
      addfield("PE_ave",&Thermo::compute_peave,FLOAT);
      peflag = 1;
    //***** Indenter Properties *****//
    } else if (strcmp(word,"indX") == 0) { addfield("indX",&Thermo::compute_indX,FLOAT);
    } else if (strcmp(word,"indY") == 0) { addfield("indY",&Thermo::compute_indY,FLOAT);
    } else if (strcmp(word,"indZ") == 0) { addfield("indZ",&Thermo::compute_indZ,FLOAT);
    } else if (strcmp(word,"indFX") == 0) { addfield("indFX",&Thermo::compute_indFX,FLOAT);
    } else if (strcmp(word,"indFY") == 0) { addfield("indFY",&Thermo::compute_indFY,FLOAT);
    } else if (strcmp(word,"indFZ") == 0) { addfield("indFZ",&Thermo::compute_indFZ,FLOAT);
    } else if (strcmp(word,"indNc") == 0) { addfield("indNc",&Thermo::compute_indNc,FLOAT);
    } else if (strcmp(word,"indAcAtomarMax") == 0) { addfield("indAcAtomarMax",&Thermo::compute_indAcAtomarMax,FLOAT);
    } else if (strcmp(word,"indAcAtomarProjected") == 0) { addfield("indAcAtomarProjected",&Thermo::compute_indAcAtomarProjected,FLOAT);
    } else if (strcmp(word,"indAcBrinell") == 0) { addfield("indAcBrinell",&Thermo::compute_indAcBrinell,FLOAT);
    } else if (strcmp(word,"indAcMeyer") == 0) { addfield("indAcMeyer",&Thermo::compute_indAcMeyer,FLOAT);
    } else if (strcmp(word,"indAcElliptic") == 0) { addfield("indAcElliptic",&Thermo::compute_indAcElliptic,FLOAT);
    } else if (strcmp(word,"indPlasticDepth") == 0) { addfield("indPlasticDepth",&Thermo::compute_indPlasticDepth,FLOAT);

    // Ackland Concentrations
    } else if (strcmp(word,"cAckUNKNOWN") == 0) { addfield("cAckUNKNOWN",&Thermo::compute_cAckUNKNOWN,FLOAT);
    } else if (strcmp(word,"cAckICO") == 0) { addfield("cAckICO",&Thermo::compute_cAckICO,FLOAT);
    } else if (strcmp(word,"cAckBCC") == 0) { addfield("cAckBCC",&Thermo::compute_cAckBCC,FLOAT);
    } else if (strcmp(word,"cAckFCC") == 0) { addfield("cAckFCC",&Thermo::compute_cAckFCC,FLOAT);
    } else if (strcmp(word,"cAckHCP") == 0) { addfield("cAckHCP",&Thermo::compute_cAckHCP,FLOAT);





    // compute value = c_ID, fix value = f_ID, variable value = v_ID
    // if no trailing [], then arg is set to 0, else arg is between []
    // copy = at most 8 chars of ID to pass to addfield
    // for compute, if has pre-compute object, first add it to list

    } else if ((strncmp(word,"c_",2) == 0) || (strncmp(word,"f_",2) == 0) ||
	       (strncmp(word,"v_",2) == 0)) {

      int n = strlen(word);
      char *id = new char[n];
      strcpy(id,&word[2]);
      char copy[9];
      strncpy(copy,id,8);
      copy[8] = '\0';

      char *ptr = strchr(id,'[');
      if (ptr) {
	if (id[strlen(id)-1] != ']')
	  error->all("Invalid keyword in thermo_style custom command");
	arg_object[nfield] = atoi(ptr+1);
	*ptr = '\0';
      } else arg_object[nfield] = 0;

      if (word[0] == 'c') {
	n = modify->find_compute(id);
	if (n < 0) error->all("Could not find thermo custom compute ID");
	if (arg_object[nfield] == 0 && modify->compute[n]->scalar_flag == 0)
	  error->all("Thermo compute ID does not compute scalar info");
	if (arg_object[nfield] > 0 && modify->compute[n]->vector_flag == 0)
	  error->all("Thermo compute ID does not compute vector info");
	if (arg_object[nfield] > 0 && 
	    arg_object[nfield] > modify->compute[n]->size_vector)
	  error->all("Thermo compute ID vector is not large enough");
	if (modify->compute[n]->id_pre)
	  int tmp = add_compute(modify->compute[n]->id_pre,arg_object[nfield]);
	field2object[nfield] = add_compute(id,arg_object[nfield]);
	addfield(copy,&Thermo::compute_compute,FLOAT);

      } else if (word[0] == 'f') {
	n = modify->find_fix(id);
	if (n < 0) error->all("Could not find thermo custom fix ID");
	field2object[nfield] = add_fix(id);
	addfield(copy,&Thermo::compute_fix,FLOAT);

      } else {
	n = input->variable->find(id);
	if (n < 0) error->all("Could not find thermo custom variable ID");
	field2object[nfield] = add_variable(id);
	addfield(copy,&Thermo::compute_variable,FLOAT);
      }

      delete [] id;

    } else error->all("Invalid keyword in thermo_style custom command");

    word = strtok(NULL," \0");
  }
}

/* ----------------------------------------------------------------------
   add field to list of quantities to print
------------------------------------------------------------------------- */

void Thermo::addfield(char *key, FnPtr func, int typeflag)
{
  strcpy(keyword[nfield],key);
  vfunc[nfield] = func;
  vtype[nfield] = typeflag;
  nfield++;
}

/* ----------------------------------------------------------------------
   add compute ID to list of Compute objects to call
   if already in list, do not add, just return location, else add to list
   convert index into which param
     index = 0 -> scalar, index >= 1 -> vector
     which = 1 -> scalar only, which = 2 -> vector only, which = 3 -> both
------------------------------------------------------------------------- */

int Thermo::add_compute(char *id, int index)
{
  int icompute;
  for (icompute = 0; icompute < ncompute; icompute++)
    if (strcmp(id,id_compute[icompute]) == 0) {
      if (index == 0 && compute_which[icompute] == 1)
	compute_which[icompute] = 2;
      if (index > 0 && compute_which[icompute] == 0)
	compute_which[icompute] = 2;
      break;
    }
  if (icompute < ncompute) return icompute;

  int n = strlen(id) + 1;
  id_compute[ncompute] = new char[n];
  strcpy(id_compute[ncompute],id);
  if (index == 0) compute_which[ncompute] = 0;
  else compute_which[ncompute] = 1;
  ncompute++;
  return ncompute-1;
}

/* ----------------------------------------------------------------------
   add fix ID to list of Fix objects to call
------------------------------------------------------------------------- */

int Thermo::add_fix(char *id)
{
  int n = strlen(id) + 1;
  id_fix[nfix] = new char[n];
  strcpy(id_fix[nfix],id);
  nfix++;
  return nfix-1;
}

/* ----------------------------------------------------------------------
   add variable ID to list of Variables to evaluate
------------------------------------------------------------------------- */

int Thermo::add_variable(char *id)
{
  int n = strlen(id) + 1;
  id_variable[nvariable] = new char[n];
  strcpy(id_variable[nvariable],id);
  nvariable++;
  return nvariable-1;
}

/* ----------------------------------------------------------------------
   create a Compute object for group all
------------------------------------------------------------------------- */

void Thermo::create_compute(char *id, char *cstyle, char *extra)
{
  char **newarg = new char*[4];
  newarg[0] = id;
  newarg[1] = "all";
  newarg[2] = cstyle;
  if (extra) newarg[3] = extra;
  if (extra) modify->add_compute(4,newarg);
  else modify->add_compute(3,newarg);
  delete [] newarg;
}

/* ----------------------------------------------------------------------
   compute a single thermodyanmic value matching word to custom list
   called when a variable is evaluated in input script
   return value as double in answer
   return 0 if OK, 1 if str is invalid
   customize a new keyword by adding to if statement
------------------------------------------------------------------------- */

int Thermo::evaluate_keyword(char *word, double *answer)
{
  // could do tests to insure user does not invoke variable at wrong time

  // all keywords:
  //   if (domain->box_exist == 0)
  //   error->all("Variable equal keyword used before simulation box defined");
  // keywords except step,atoms,vol,ly,ly,lz
  //   if (update->first_update == 0)
  //     error->all("Variable equal keyword used before initial run");
  //   if (update->ntimestep != output->next_thermo && me == 0)
  //     error->warning("Variable equal keyword used with non-current thermo");
  // keywords that require Compute objects like T,P
  //   ptrs like temperature will not be set if init() hasn't been called
  //   Compute objects will not exist if thermo style isn't using them

  // toggle thermoflag off/on so inidividual compute routines don't assume
  //   they are being called from compute() which pre-calls Compute objects

  thermoflag = 0;

  if (strcmp(word,"step") == 0) {
    compute_step();
    dvalue = ivalue;
  } else if (strcmp(word,"atoms") == 0) {
    compute_atoms();
    dvalue = ivalue;
  }
  else if (strcmp(word,"cpu") == 0) compute_cpu();

  else if (strcmp(word,"temp") == 0) compute_temp();
  else if (strcmp(word,"press") == 0) compute_press();
  else if (strcmp(word,"pe") == 0) compute_pe();
  else if (strcmp(word,"ke") == 0) compute_ke();
  else if (strcmp(word,"etotal") == 0) compute_etotal();
  else if (strcmp(word,"enthalpy") == 0) compute_enthalpy();

  else if (strcmp(word,"evdwl") == 0) compute_evdwl();
  else if (strcmp(word,"ecoul") == 0) compute_ecoul();
  else if (strcmp(word,"epair") == 0) compute_epair();
  else if (strcmp(word,"ebond") == 0) compute_ebond();
  else if (strcmp(word,"eangle") == 0) compute_eangle();
  else if (strcmp(word,"edihed") == 0) compute_edihed();
  else if (strcmp(word,"eimp") == 0) compute_eimp();
  else if (strcmp(word,"emol") == 0) compute_emol();
  else if (strcmp(word,"elong") == 0) compute_elong();
  else if (strcmp(word,"etail") == 0) compute_etail();

  else if (strcmp(word,"vol") == 0) compute_vol();
  else if (strcmp(word,"lx") == 0) compute_lx();
  else if (strcmp(word,"ly") == 0) compute_ly();
  else if (strcmp(word,"lz") == 0) compute_lz();

  else if (strcmp(word,"xlo") == 0) compute_xlo();
  else if (strcmp(word,"xhi") == 0) compute_xhi();
  else if (strcmp(word,"ylo") == 0) compute_ylo();
  else if (strcmp(word,"yhi") == 0) compute_yhi();
  else if (strcmp(word,"zlo") == 0) compute_zlo();
  else if (strcmp(word,"zhi") == 0) compute_zhi();

  else if (strcmp(word,"pxx") == 0) compute_pxx();
  else if (strcmp(word,"pyy") == 0) compute_pyy();
  else if (strcmp(word,"pzz") == 0) compute_pzz();
  else if (strcmp(word,"pxy") == 0) compute_pxy();
  else if (strcmp(word,"pxz") == 0) compute_pxz();
  else if (strcmp(word,"pyz") == 0) compute_pyz();

  else if (strcmp(word,"drot") == 0) compute_drot();
  else if (strcmp(word,"grot") == 0) compute_grot();

  else if (strcmp(word,"tave") == 0) compute_tave();
  else if (strcmp(word,"pave") == 0) compute_pave();
  else if (strcmp(word,"eave") == 0) compute_eave();
  else if (strcmp(word,"peave") == 0) compute_peave();

  else return 1;

  thermoflag = 1;

  *answer = dvalue;
  return 0;
}

/* ----------------------------------------------------------------------
   extraction of Compute, Fix, Variable results
   ignore thermoflag, since these 3 routines only called by Thermo::compute(),
     not by variable evaluation
   compute value is normalized by atoms if returning extensive quantities
   fix value is normalized (so should return extensive quantity)
   variable value is not normalized (so formula should normalize if desired)
------------------------------------------------------------------------- */

/* ---------------------------------------------------------------------- */

void Thermo::compute_compute()
{
  int index = field2object[ifield];
  if (arg_object[ifield] == 0) dvalue = computes[index]->scalar;
  else dvalue = computes[index]->vector[arg_object[ifield]];
  if (computes[index]->extensive && normflag) dvalue /= natoms;
}

/* ---------------------------------------------------------------------- */

void Thermo::compute_fix()
{
  int index = field2object[ifield];
  dvalue = fixes[index]->thermo(arg_object[ifield]);
  if (normflag) dvalue /= natoms;
}

/* ---------------------------------------------------------------------- */

void Thermo::compute_variable()
{
  int index = field2object[ifield];
  dvalue = atof(input->variable->retrieve(id_variable[index]));
}

/* ----------------------------------------------------------------------
   one routine for every thermo keyword can output
   thermoflag = 1 if called by Thermo::compute()
   thermoflag = 0 if called by evaluate_keyword() via Variable
   set ivalue/dvalue if value is integer/double
   customize a new keyword by adding a method
------------------------------------------------------------------------- */

void Thermo::compute_step()
{
  ivalue = update->ntimestep;
}

/* ---------------------------------------------------------------------- */

void Thermo::compute_atoms()
{
  if (thermoflag) ivalue = static_cast<int> (natoms);
  else ivalue = static_cast<int> (atom->natoms);
}

/* ---------------------------------------------------------------------- */

void Thermo::compute_cpu()
{
  if (firststep == 0) dvalue = 0.0;
  else dvalue = timer->elapsed(TIME_LOOP);
}

/* ---------------------------------------------------------------------- */

void Thermo::compute_temp()
{
  if (thermoflag) dvalue = temperature->scalar;
  else dvalue = temperature->compute_scalar();
}

/* ---------------------------------------------------------------------- */

void Thermo::compute_press()
{
  if (thermoflag) dvalue = pressure->scalar;
  else {
    dvalue = temperature->compute_scalar();
    dvalue = pressure->compute_scalar();
  }
}

/* ---------------------------------------------------------------------- */

void Thermo::compute_pe()
{
  double tmp = 0.0;
  if (force->pair) tmp += force->pair->eng_vdwl + force->pair->eng_coul;
  if (force->bond) tmp += force->bond->eng_vdwl;
  if (force->dihedral) 
    tmp += force->dihedral->eng_vdwl + force->dihedral->eng_coul;

  if (atom->molecular) {
    if (force->bond) tmp += force->bond->energy;
    if (force->angle) tmp += force->angle->energy;
    if (force->dihedral) tmp += force->dihedral->energy;
    if (force->improper) tmp += force->improper->energy;
  }

  MPI_Allreduce(&tmp,&dvalue,1,MPI_DOUBLE,MPI_SUM,world);

  if (force->kspace) dvalue += force->kspace->energy;
  if (force->pair && force->pair->tail_flag) {
    double volume = domain->xprd * domain->yprd * domain->zprd;
    dvalue += force->pair->etail / volume;
  }
  if (modify->n_thermo_energy) dvalue += modify->thermo_energy();

  if (normflag) dvalue /= natoms;

  // also store PE in potential_energy so other classes can grab it

  potential_energy = dvalue;
}

/* ---------------------------------------------------------------------- */

void Thermo::compute_ke()
{
  if (thermoflag) dvalue = temperature->scalar;
  else dvalue = temperature->compute_scalar();
  dvalue *= 0.5 * temperature->dof * force->boltz;
  if (normflag) dvalue /= natoms;
}

/* ---------------------------------------------------------------------- */

void Thermo::compute_etotal()
{
  compute_pe();
  double ke;
  if (thermoflag) ke = temperature->scalar;
  else ke = temperature->compute_scalar();
  ke *= 0.5 * temperature->dof * force->boltz;
  if (normflag) ke /= natoms;
  dvalue += ke;
}

/* ---------------------------------------------------------------------- */

void Thermo::compute_enthalpy()
{
  compute_etotal();
  double etmp = dvalue;

  compute_vol();
  double vtmp = dvalue;
  if (normflag) vtmp /= natoms;

  compute_press();
  double ptmp = dvalue;

  dvalue = etmp + ptmp*vtmp/(force->nktv2p);
}

/* ---------------------------------------------------------------------- */

void Thermo::compute_evdwl()
{
  double tmp = 0.0;
  if (force->pair) tmp += force->pair->eng_vdwl;
  if (force->bond) tmp += force->bond->eng_vdwl;
  if (force->dihedral) tmp += force->dihedral->eng_vdwl;

  MPI_Allreduce(&tmp,&dvalue,1,MPI_DOUBLE,MPI_SUM,world);

  if (force->pair && force->pair->tail_flag) {
    double volume = domain->xprd * domain->yprd * domain->zprd;
    dvalue += force->pair->etail / volume;
  }

  if (normflag) dvalue /= natoms;
}

/* ---------------------------------------------------------------------- */

void Thermo::compute_ecoul()
{
  double tmp = 0.0;
  if (force->pair) tmp += force->pair->eng_coul;
  if (force->dihedral) tmp += force->dihedral->eng_coul;
  MPI_Allreduce(&tmp,&dvalue,1,MPI_DOUBLE,MPI_SUM,world);
  if (normflag) dvalue /= natoms;
}

/* ---------------------------------------------------------------------- */

void Thermo::compute_epair()
{
  double tmp = 0.0;
  if (force->pair) tmp += force->pair->eng_vdwl + force->pair->eng_coul;
  if (force->bond) tmp += force->bond->eng_vdwl;
  if (force->dihedral) 
    tmp += force->dihedral->eng_vdwl + force->dihedral->eng_coul;
  MPI_Allreduce(&tmp,&dvalue,1,MPI_DOUBLE,MPI_SUM,world);

  if (force->kspace) dvalue += force->kspace->energy;
  if (force->pair && force->pair->tail_flag) {
    double volume = domain->xprd * domain->yprd * domain->zprd;
    dvalue += force->pair->etail / volume;
  }

  if (normflag) dvalue /= natoms;
}

/* ---------------------------------------------------------------------- */

void Thermo::compute_ebond()
{
  if (force->bond) {
    double tmp = force->bond->energy;
    MPI_Allreduce(&tmp,&dvalue,1,MPI_DOUBLE,MPI_SUM,world);
    if (normflag) dvalue /= natoms;
  } else dvalue = 0.0;
}

/* ---------------------------------------------------------------------- */

void Thermo::compute_eangle()
{
  if (force->angle) {
    double tmp = force->angle->energy;
    MPI_Allreduce(&tmp,&dvalue,1,MPI_DOUBLE,MPI_SUM,world);
    if (normflag) dvalue /= natoms;
  } else dvalue = 0.0;
}

/* ---------------------------------------------------------------------- */

void Thermo::compute_edihed()
{
  if (force->dihedral) {
    double tmp = force->dihedral->energy;
    MPI_Allreduce(&tmp,&dvalue,1,MPI_DOUBLE,MPI_SUM,world);
    if (normflag) dvalue /= natoms;
  } else dvalue = 0.0;
}

/* ---------------------------------------------------------------------- */

void Thermo::compute_eimp()
{
  if (force->improper) {
    double tmp = force->improper->energy;
    MPI_Allreduce(&tmp,&dvalue,1,MPI_DOUBLE,MPI_SUM,world);
    if (normflag) dvalue /= natoms;
  } else dvalue = 0.0;
}

/* ---------------------------------------------------------------------- */

void Thermo::compute_emol()
{
  double tmp = 0.0;
  if (atom->molecular) {
    if (force->bond) tmp += force->bond->energy;
    if (force->angle) tmp += force->angle->energy;
    if (force->dihedral) tmp += force->dihedral->energy;
    if (force->improper) tmp += force->improper->energy;
    MPI_Allreduce(&tmp,&dvalue,1,MPI_DOUBLE,MPI_SUM,world);
    if (normflag) dvalue /= natoms;
  } else dvalue = 0.0;
}

/* ---------------------------------------------------------------------- */

void Thermo::compute_elong()
{
  if (force->kspace) {
    dvalue = force->kspace->energy;
    if (normflag) dvalue /= natoms;
  } else dvalue = 0.0;
}

/* ---------------------------------------------------------------------- */

void Thermo::compute_etail()
{
  if (force->pair && force->pair->tail_flag) {
    double volume = domain->xprd * domain->yprd * domain->zprd;
    dvalue = force->pair->etail / volume;
    if (normflag) dvalue /= natoms;
  } else dvalue = 0.0;
}

/* ---------------------------------------------------------------------- */

void Thermo::compute_vol()
{
  dvalue = domain->xprd * domain->yprd * domain->zprd;
}

/* ---------------------------------------------------------------------- */

void Thermo::compute_lx()
{
  dvalue = domain->xprd;
}

/* ---------------------------------------------------------------------- */

void Thermo::compute_ly()
{
  dvalue = domain->yprd;
}

/* ---------------------------------------------------------------------- */

void Thermo::compute_lz()
{
  dvalue = domain->zprd;
}

/* ---------------------------------------------------------------------- */

void Thermo::compute_xlo()
{
  dvalue = domain->boxlo[0];
}

/* ---------------------------------------------------------------------- */

void Thermo::compute_xhi()
{
  dvalue = domain->boxhi[0];
}

/* ---------------------------------------------------------------------- */

void Thermo::compute_ylo()
{
  dvalue = domain->boxlo[1];
}

/* ---------------------------------------------------------------------- */

void Thermo::compute_yhi()
{
  dvalue = domain->boxhi[1];
}

/* ---------------------------------------------------------------------- */

void Thermo::compute_zlo()
{
  dvalue = domain->boxlo[2];
}

/* ---------------------------------------------------------------------- */

void Thermo::compute_zhi()
{
  dvalue = domain->boxhi[2];
}

/* ---------------------------------------------------------------------- */

void Thermo::compute_pxx()
{
  if (!thermoflag) {
    temperature->compute_vector();
    pressure->compute_vector();
  }
  dvalue = pressure->vector[0];
}

/* ---------------------------------------------------------------------- */

void Thermo::compute_pyy()
{
  if (!thermoflag) {
    temperature->compute_vector();
    pressure->compute_vector();
  }
  dvalue = pressure->vector[1];
}

/* ---------------------------------------------------------------------- */

void Thermo::compute_pzz()
{
  if (!thermoflag) {
    temperature->compute_vector();
    pressure->compute_vector();
  }
  dvalue = pressure->vector[2];
}

/* ---------------------------------------------------------------------- */

void Thermo::compute_pxy()
{
  if (!thermoflag) {
    temperature->compute_vector();
    pressure->compute_vector();
  }
  dvalue = pressure->vector[3];
}

/* ---------------------------------------------------------------------- */

void Thermo::compute_pxz()
{
  if (!thermoflag) {
    temperature->compute_vector();
    pressure->compute_vector();
  }
  dvalue = pressure->vector[4];
}

/* ---------------------------------------------------------------------- */

void Thermo::compute_pyz()
{
  if (!thermoflag) {
    temperature->compute_vector();
    pressure->compute_vector();
  }
  dvalue = pressure->vector[5];
}

/* ---------------------------------------------------------------------- */

void Thermo::compute_drot()
{
  if (thermoflag) dvalue = rotate_dipole->scalar;
  else dvalue = rotate_dipole->compute_scalar();
  if (normflag) dvalue /= natoms;
}

/* ---------------------------------------------------------------------- */

void Thermo::compute_grot()
{
  if (thermoflag) dvalue = rotate_gran->scalar;
  else dvalue = rotate_gran->compute_scalar();
  if (normflag) dvalue /= natoms;
}

/* ---------------------------------------------------------------------- */

void Thermo::compute_tave()
{
  compute_temp();

  if (firststep == 0) {
    if (npartial_t == 0) {
      ncurrent_t = 0;
      npartial_t = 1;
      tsum = tpast[0] = dvalue;
    }
    dvalue = tsum/npartial_t;
  } else {
    ncurrent_t++;
    if (ncurrent_t == nwindow) ncurrent_t = 0;
    if (npartial_t == nwindow) tsum -= tpast[ncurrent_t];
    else npartial_t++;
    tpast[ncurrent_t] = dvalue;
    tsum += tpast[ncurrent_t];
    dvalue = tsum/npartial_t;
  }
}

/* ---------------------------------------------------------------------- */

void Thermo::compute_pave()
{
  compute_press();

  if (firststep == 0) {
    if (npartial_p == 0) {
      ncurrent_p = 0;
      npartial_p = 1;
      psum = ppast[0] = dvalue;
    }
    dvalue = psum/npartial_p;
  } else {
    ncurrent_p++;
    if (ncurrent_p == nwindow) ncurrent_p = 0;
    if (npartial_p == nwindow) psum -= ppast[ncurrent_p];
    else npartial_p++;
    ppast[ncurrent_p] = dvalue;
    psum += ppast[ncurrent_p];
    dvalue = psum/npartial_p;
  }
}

/* ---------------------------------------------------------------------- */

void Thermo::compute_eave()
{
  compute_etotal();

  if (firststep == 0) {
    if (npartial_e == 0) {
      ncurrent_e = 0;
      npartial_e = 1;
      esum = epast[0] = dvalue;
    }
    dvalue = esum/npartial_e;
  } else {
    ncurrent_e++;
    if (ncurrent_e == nwindow) ncurrent_e = 0;
    if (npartial_e == nwindow) esum -= epast[ncurrent_e];
    else npartial_e++;
    epast[ncurrent_e] = dvalue;
    esum += epast[ncurrent_e];
    dvalue = esum/npartial_e;
  }
}

/* ---------------------------------------------------------------------- */

void Thermo::compute_peave()
{
  compute_pe();

  if (firststep == 0) {
    if (npartial_pe == 0) {
      ncurrent_pe = 0;
      npartial_pe = 1;
      pesum = pepast[0] = dvalue;
    }
    dvalue = pesum/npartial_pe;
  } else {
    ncurrent_pe++;
    if (ncurrent_pe == nwindow) ncurrent_pe = 0;
    if (npartial_pe == nwindow) pesum -= pepast[ncurrent_pe];
    else npartial_pe++;
    pepast[ncurrent_pe] = dvalue;
    pesum += pepast[ncurrent_pe];
    dvalue = pesum/npartial_pe;
  }
}

/* ---------------------------------------------------------------------- */

//***** Indenter Properties *****//
void Thermo::compute_indX() { 
   int IndexLocal = FindIndexFixIndent(); 
#warning FIXME: Only and only if indD is the first of all indenter thermal outputs, it will use the right index!
   if (IndexFixIndent == IndexFixIndentBad) dvalue=0.0;
   else dvalue = ((FixIndent *) modify->fix[IndexFixIndent])->IndX; 
}
void Thermo::compute_indY() { 
   if (IndexFixIndent == IndexFixIndentBad) dvalue=0.0;
   else dvalue = ((FixIndent *) modify->fix[IndexFixIndent])->IndY; 
}
void Thermo::compute_indZ() { 
   if (IndexFixIndent == IndexFixIndentBad) dvalue=0.0;
   else dvalue = ((FixIndent *) modify->fix[IndexFixIndent])->IndZ; 
}

/* ---------------------------------------------------------------------- */

// Indentation Force
void Thermo::compute_indFX() { 
   if (IndexFixIndent == IndexFixIndentBad) dvalue=0.0;
   else dvalue = ((FixIndent *) modify->fix[IndexFixIndent])->FIndX; 
}
void Thermo::compute_indFY() { 
   if (IndexFixIndent == IndexFixIndentBad) dvalue=0.0;
   else dvalue = ((FixIndent *) modify->fix[IndexFixIndent])->FIndY; 
}
void Thermo::compute_indFZ() { 
   if (IndexFixIndent == IndexFixIndentBad) dvalue=0.0;
   else dvalue = ((FixIndent *) modify->fix[IndexFixIndent])->FIndZ; 
}

/* ---------------------------------------------------------------------- */

// Indenter Contact Properties
void Thermo::compute_indNc() { 
   if (IndexFixIndent == IndexFixIndentBad) dvalue=0;
   else dvalue = ((FixIndent *) modify->fix[IndexFixIndent])->Nc; 
}
void Thermo::compute_indAcAtomarMax() { 
   if (IndexFixIndent == IndexFixIndentBad) dvalue=0.0;
   else dvalue = ((FixIndent *) modify->fix[IndexFixIndent])->AcAtomarMax; 
}
void Thermo::compute_indAcAtomarProjected() { 
   if (IndexFixIndent == IndexFixIndentBad) dvalue=0.0;
   else dvalue = ((FixIndent *) modify->fix[IndexFixIndent])->AcAtomarProjected; 
}
void Thermo::compute_indAcBrinell() {
   if (IndexFixIndent == IndexFixIndentBad) dvalue=0.0;
   else dvalue = ((FixIndent *) modify->fix[IndexFixIndent])->AcBrinell;
}
void Thermo::compute_indAcMeyer() {
   if (IndexFixIndent == IndexFixIndentBad) dvalue=0.0;
   else dvalue = ((FixIndent *) modify->fix[IndexFixIndent])->AcMeyer;
}
void Thermo::compute_indAcElliptic() {
   if (IndexFixIndent == IndexFixIndentBad) dvalue=0.0;
   else dvalue = ((FixIndent *) modify->fix[IndexFixIndent])->AcElliptic;
}
void Thermo::compute_indPlasticDepth() {
   if (IndexFixIndent == IndexFixIndentBad) dvalue=0.0;
   else dvalue = ((FixIndent *) modify->fix[IndexFixIndent])->PlasticDepth;
}

/* ---------------------------------------------------------------------- */

// Defect Concentrations
void Thermo::compute_cAckUNKNOWN() {
   int IndexLocal = FindIndexComputeAckland ();
#warning FIXME: Only and only if cAckUNKNOWN is the first of all ackland concentrations, it will use the right index!
   if (IndexComputeAckland == IndexComputeAcklandBad) dvalue=0.;
   else dvalue = ((ComputeAcklandAtom *) modify->compute[IndexComputeAckland])->ConcentrationUNKNOWN;
}

void Thermo::compute_cAckICO() {
   if (IndexComputeAckland == IndexComputeAcklandBad) dvalue=0.;
   else dvalue = ((ComputeAcklandAtom *) modify->compute[IndexComputeAckland])->ConcentrationICO;
}

void Thermo::compute_cAckBCC() {
   if (IndexComputeAckland == IndexComputeAcklandBad) dvalue=0.;
   else dvalue = ((ComputeAcklandAtom *) modify->compute[IndexComputeAckland])->ConcentrationBCC;
}
void Thermo::compute_cAckFCC() {
   if (IndexComputeAckland == IndexComputeAcklandBad) dvalue=0.;
   else dvalue = ((ComputeAcklandAtom *) modify->compute[IndexComputeAckland])->ConcentrationFCC;
}
void Thermo::compute_cAckHCP() {
   if (IndexComputeAckland == IndexComputeAcklandBad) dvalue=0.;
   else dvalue = ((ComputeAcklandAtom *) modify->compute[IndexComputeAckland])->ConcentrationHCP;
}


//***** *****//
