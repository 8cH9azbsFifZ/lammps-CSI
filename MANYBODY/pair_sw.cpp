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
   Contributing author: Aidan Thompson (SNL)
------------------------------------------------------------------------- */

#include "math.h"
#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "pair_sw.h"
#include "atom.h"
#include "force.h"
#include "comm.h"
#include "update.h"
#include "memory.h"
#include "neighbor.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;

#define MAXLINE 1024
#define DELTA 4

/* ---------------------------------------------------------------------- */

PairSW::PairSW(LAMMPS *lmp) : Pair(lmp)
{
  neigh_half_every = 0;
  neigh_full_every = 1;
  single_enable = 0;
  one_coeff = 1;

  nelements = 0;
  elements = NULL;
  nparams = 0;
  maxparam = 0;
  params = NULL;
  elem2param = NULL;
}

/* ----------------------------------------------------------------------
   check if allocated, since class can be destructed when incomplete
------------------------------------------------------------------------- */

PairSW::~PairSW()
{
  if (elements)
    for (int i = 0; i < nelements; i++) delete [] elements[i];
  delete [] elements;
  memory->sfree(params);
  memory->destroy_3d_int_array(elem2param);

  if (allocated) {
    memory->destroy_2d_int_array(setflag);
    memory->destroy_2d_double_array(cutsq);
    delete [] map;
  }
}

/* ---------------------------------------------------------------------- */

void PairSW::compute(int eflag, int vflag)
{
  int i,j,k,m,n,itag,jtag,itype,jtype,ktype,iparam,numneigh;
  double xtmp,ytmp,ztmp,delx,dely,delz;
  double rsq,rsq1,rsq2,eng,fforce;
  double delr1[3],delr2[3],fj[3],fk[3];
  int *neighs;
  double **f;

  eng_vdwl = 0.0;
  if (vflag) for (i = 0; i < 6; i++) virial[i] = 0.0;

  if (vflag == 2) f = update->f_pair;
  else f = atom->f;
  double **x = atom->x;
  int *tag = atom->tag;
  int *type = atom->type;
  int nlocal = atom->nlocal;

  // loop over full neighbor list of my atoms

  for (i = 0; i < nlocal; i++) {
    itag = tag[i];
    itype = map[type[i]];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];

    // two-body interactions, skip half of them

    neighs = neighbor->firstneigh_full[i];
    numneigh = neighbor->numneigh_full[i];

    for (m = 0; m < numneigh; m++) {
      j = neighs[m];
      jtag = tag[j];

      if (itag > jtag) {
	if ((itag+jtag) % 2 == 0) continue;
      } else if (itag < jtag) {
	if ((itag+jtag) % 2 == 1) continue;
      } else {
	if (x[j][2] < ztmp) continue;
	else if (x[j][2] == ztmp && x[j][1] < ytmp) continue;
	else if (x[j][2] == ztmp && x[j][1] == ytmp && x[j][0] < xtmp)
	    continue;
      }

      jtype = map[type[j]];

      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx*delx + dely*dely + delz*delz;

      iparam = elem2param[itype][jtype][jtype];
      if (rsq > params[iparam].cutsq) continue;

      twobody(&params[iparam],rsq,fforce,eflag,eng);

      if (eflag) eng_vdwl += eng;
      f[i][0] += fforce*delx;
      f[i][1] += fforce*dely;
      f[i][2] += fforce*delz;
      f[j][0] -= fforce*delx;
      f[j][1] -= fforce*dely;
      f[j][2] -= fforce*delz;
    }

    // three-body interactions
    // cannot test I-J distance against cutoff outside of 2nd loop
    // b/c must use I-J-K cutoff for both rij and rik

    for (m = 0; m < numneigh-1; m++) {
      j = neighs[m];
      jtype = map[type[j]];

      for (n = m+1; n < numneigh; n++) {
	k = neighs[n];
	ktype = map[type[k]];
	iparam = elem2param[itype][jtype][ktype];

	delr1[0] = x[j][0] - xtmp;
	delr1[1] = x[j][1] - ytmp;
	delr1[2] = x[j][2] - ztmp;
	rsq1 = delr1[0]*delr1[0] + delr1[1]*delr1[1] + delr1[2]*delr1[2];
	if (rsq1 > params[iparam].cutsq) continue;

	delr2[0] = x[k][0] - xtmp;
	delr2[1] = x[k][1] - ytmp;
	delr2[2] = x[k][2] - ztmp;
	rsq2 = delr2[0]*delr2[0] + delr2[1]*delr2[1] + delr2[2]*delr2[2];
	if (rsq2 > params[iparam].cutsq) continue;

	threebody(&params[iparam],rsq1,rsq2,delr1,delr2,fj,fk,eflag,eng);

	if (eflag) eng_vdwl += eng;
	f[i][0] -= fj[0] + fk[0];
	f[i][1] -= fj[1] + fk[1];
	f[i][2] -= fj[2] + fk[2];
	f[j][0] += fj[0];
	f[j][1] += fj[1];
	f[j][2] += fj[2];
	f[k][0] += fk[0];
	f[k][1] += fk[1];
	f[k][2] += fk[2];
      }
    }
  }
  if (vflag == 2) virial_compute();
}

/* ---------------------------------------------------------------------- */

void PairSW::allocate()
{
  allocated = 1;
  int n = atom->ntypes;

  setflag = memory->create_2d_int_array(n+1,n+1,"pair:setflag");
  cutsq = memory->create_2d_double_array(n+1,n+1,"pair:cutsq");

  map = new int[n+1];
}

/* ----------------------------------------------------------------------
   global settings 
------------------------------------------------------------------------- */

void PairSW::settings(int narg, char **arg)
{
  if (narg != 0) error->all("Illegal pair_style command");
}

/* ----------------------------------------------------------------------
   set coeffs for one or more type pairs
------------------------------------------------------------------------- */

void PairSW::coeff(int narg, char **arg)
{
  int i,j,n;

  if (!allocated) allocate();

  if (narg != 3 + atom->ntypes)
    error->all("Incorrect args for pair coefficients");

  // insure I,J args are * *

  if (strcmp(arg[0],"*") != 0 || strcmp(arg[1],"*") != 0)
    error->all("Incorrect args for pair coefficients");

  // read args that map atom types to elements in potential file
  // map[i] = which element the Ith atom type is, -1 if NULL
  // nelements = # of unique elements
  // elements = list of element names

  if (elements) {
    for (i = 0; i < nelements; i++) delete [] elements[i];
    delete [] elements;
  }
  elements = new char*[atom->ntypes];
  for (i = 0; i < atom->ntypes; i++) elements[i] = NULL;

  nelements = 0;
  for (i = 3; i < narg; i++) {
    if (strcmp(arg[i],"NULL") == 0) {
      map[i-2] = -1;
      continue;
    }
    for (j = 0; j < nelements; j++)
      if (strcmp(arg[i],elements[j]) == 0) break;
    map[i-2] = j;
    if (j == nelements) {
      n = strlen(arg[i]) + 1;
      elements[j] = new char[n];
      strcpy(elements[j],arg[i]);
      nelements++;
    }
  }

  // read potential file and initialize potential parameters
  
  read_file(arg[2]);
  setup();

  // clear setflag since coeff() called once with I,J = * *

  n = atom->ntypes;
  for (int i = 1; i <= n; i++)
    for (int j = i; j <= n; j++)
      setflag[i][j] = 0;

  // set setflag i,j for type pairs where both are mapped to elements

  int count = 0;
  for (int i = 1; i <= n; i++)
    for (int j = i; j <= n; j++)
      if (map[i] >= 0 && map[j] >= 0) {
	setflag[i][j] = 1;
	count++;
      }

  if (count == 0) error->all("Incorrect args for pair coefficients");
}

/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

double PairSW::init_one(int i, int j)
{
  if (setflag[i][j] == 0) error->all("All pair coeffs are not set");

  return cutmax;
}

/* ---------------------------------------------------------------------- */

void PairSW::init_style()
{
  if (atom->tag_enable == 0)
    error->all("Pair style Stillinger-Weber requires atom IDs");
  if (force->newton_pair == 0)
    error->all("Pair style Stillinger-Weber requires newton pair on");
}

/* ---------------------------------------------------------------------- */

void PairSW::read_file(char *file)
{
  int params_per_line = 13;
  char **words = new char*[params_per_line+1];

  if (params) delete [] params;
  params = NULL;
  nparams = 0;

  // open file on proc 0

  FILE *fp;
  if (comm->me == 0) {
    fp = fopen(file,"r");
    if (fp == NULL) {
      char str[128];
      sprintf(str,"Cannot open Stillinger-Weber potential file %s",file);
      error->one(str);
    }
  }

  // read each set of params from potential file
  // one set of params can span multiple lines
  // store params if all 3 element tags are in element list

  int n,nwords,ielement,jelement,kelement;
  char line[MAXLINE],*ptr;
  int eof = 0;

  while (1) {
    if (comm->me == 0) {
      ptr = fgets(line,MAXLINE,fp);
      if (ptr == NULL) {
	eof = 1;
	fclose(fp);
      } else n = strlen(line) + 1;
    }
    MPI_Bcast(&eof,1,MPI_INT,0,world);
    if (eof) break;
    MPI_Bcast(&n,1,MPI_INT,0,world);
    MPI_Bcast(line,n,MPI_CHAR,0,world);

    // strip comment, skip line if blank

    if (ptr = strchr(line,'#')) *ptr = '\0';
    nwords = atom->count_words(line);
    if (nwords == 0) continue;

    // concatenate additional lines until have params_per_line words

    while (nwords < params_per_line) {
      n = strlen(line);
      if (comm->me == 0) {
        ptr = fgets(&line[n],MAXLINE-n,fp);
        if (ptr == NULL) {
	  eof = 1;
	  fclose(fp);
        } else n = strlen(line) + 1;
      }
      MPI_Bcast(&eof,1,MPI_INT,0,world);
      if (eof) break;
      MPI_Bcast(&n,1,MPI_INT,0,world);
      MPI_Bcast(line,n,MPI_CHAR,0,world);
      if (ptr = strchr(line,'#')) *ptr = '\0';
      nwords = atom->count_words(line);
    }

    if (nwords != params_per_line)
      error->all("Incorrect format in Stillinger-Weber potential file");

    // words = ptrs to all words in line

    nwords = 0;
    words[nwords++] = strtok(line," \t\n\r\f");
    while (words[nwords++] = strtok(NULL," \t\n\r\f")) continue;

    // ielement,jelement,kelement = 1st args
    // if all 3 args are in element list, then parse this line
    // else skip to next entry in file

    for (ielement = 0; ielement < nelements; ielement++)
      if (strcmp(words[0],elements[ielement]) == 0) break;
    if (ielement == nelements) continue;
    for (jelement = 0; jelement < nelements; jelement++)
      if (strcmp(words[1],elements[jelement]) == 0) break;
    if (jelement == nelements) continue;
    for (kelement = 0; kelement < nelements; kelement++)
      if (strcmp(words[2],elements[kelement]) == 0) break;
    if (kelement == nelements) continue;

    // load up parameter settings and error check their values

    if (nparams == maxparam) {
      maxparam += DELTA;
      params = (Param *) memory->srealloc(params,maxparam*sizeof(Param),
					  "pair:params");
    }

    params[nparams].ielement = ielement;
    params[nparams].jelement = jelement;
    params[nparams].kelement = kelement;
    params[nparams].epsilon = atof(words[3]);
    params[nparams].sigma = atof(words[4]);
    params[nparams].littlea = atof(words[5]);
    params[nparams].lambda = atof(words[6]);
    params[nparams].gamma = atof(words[7]);
    params[nparams].costheta = atof(words[8]);
    params[nparams].biga = atof(words[9]);
    params[nparams].bigb = atof(words[10]);
    params[nparams].powerp = atof(words[11]);
    params[nparams].powerq = atof(words[12]);

    if (params[nparams].epsilon < 0.0 || params[nparams].sigma < 0.0 || 
	params[nparams].littlea < 0.0 || params[nparams].lambda < 0.0 ||
	params[nparams].gamma < 0.0 || params[nparams].biga < 0.0 || 
	params[nparams].bigb < 0.0 || params[nparams].powerp < 0.0 ||
	params[nparams].powerq < 0.0)
      error->all("Illegal Stillinger-Weber parameter");

    nparams++;
  }

  delete [] words;
}

/* ---------------------------------------------------------------------- */

void PairSW::setup()
{
  int i,j,k,m,n;

  // set elem2param for all triplet combinations
  // must be a single exact match to lines read from file
  // do not allow for ACB in place of ABC

  if (elem2param) memory->destroy_3d_int_array(elem2param);
  elem2param = memory->create_3d_int_array(nelements,nelements,nelements,
					   "pair:elem2param");

  for (i = 0; i < nelements; i++)
    for (j = 0; j < nelements; j++)
      for (k = 0; k < nelements; k++) {
	n = -1;
	for (m = 0; m < nparams; m++) {
	  if (i == params[m].ielement && j == params[m].jelement && 
	      k == params[m].kelement) {
	    if (n >= 0) error->all("Potential file has duplicate entry");
	    n = m;
	  }
	}
	if (n < 0) error->all("Potential file is missing an entry");
	elem2param[i][j][k] = n;
      }


  // compute parameter values derived from inputs

  for (m = 0; m < nparams; m++) {
    params[m].cut = params[m].sigma*params[m].littlea;
    params[m].cutsq = params[m].cut*params[m].cut;
    params[m].sigma_gamma = params[m].sigma*params[m].gamma;
    params[m].lambda_epsilon = params[m].lambda*params[m].epsilon;
    params[m].lambda_epsilon2 = 2.0*params[m].lambda*params[m].epsilon;
    params[m].c1 = params[m].biga*params[m].epsilon * 
      params[m].powerp*params[m].bigb * 
      pow(params[m].sigma,params[m].powerp);
    params[m].c2 = params[m].biga*params[m].epsilon*params[m].powerq * 
      pow(params[m].sigma,params[m].powerq);
    params[m].c3 = params[m].biga*params[m].epsilon*params[m].bigb * 
      pow(params[m].sigma,params[m].powerp+1.0);
    params[m].c4 = params[m].biga*params[m].epsilon * 
      pow(params[m].sigma,params[m].powerq+1.0);
    params[m].c5 = params[m].biga*params[m].epsilon*params[m].bigb * 
      pow(params[m].sigma,params[m].powerp);
    params[m].c6 = params[m].biga*params[m].epsilon *
      pow(params[m].sigma,params[m].powerq);
  }

  // set cutmax to max of all params

  cutmax = 0.0;
  for (m = 0; m < nparams; m++)
    if (params[m].cut > cutmax) cutmax = params[m].cut;
}  

/* ---------------------------------------------------------------------- */

void PairSW::twobody(Param *param, double rsq, double &fforce,
		     int eflag, double &eng)
{
  double r,rinvsq,rp,rq,rainv,rainvsq,expsrainv;

  r = sqrt(rsq);
  rinvsq = 1.0/rsq;
  rp = pow(r,-param->powerp);
  rq = pow(r,-param->powerq);
  rainv = 1.0 / (r - param->cut);
  rainvsq = rainv*rainv*r; 
  expsrainv = exp(param->sigma * rainv);
  fforce = (param->c1*rp - param->c2*rq +
	    (param->c3*rp -param->c4*rq) * rainvsq) * expsrainv * rinvsq;
  if (eflag) eng = (param->c5*rp - param->c6*rq) * expsrainv;
}

/* ---------------------------------------------------------------------- */

void PairSW::threebody(Param *param, double rsq1, double rsq2,
		       double *delr1, double *delr2,
		       double *fj, double *fk, int eflag, double &eng)
{
  double r1,rinvsq1,rainv1,gsrainv1,gsrainvsq1,expgsrainv1;
  double r2,rinvsq2,rainv2,gsrainv2,gsrainvsq2,expgsrainv2;
  double rinv12,cs,delcs,delcssq,facexp,facrad,frad1,frad2;
  double facang,facang12,csfacang,csfac1,csfac2;

  r1 = sqrt(rsq1);
  rinvsq1 = 1.0/rsq1;
  rainv1 = 1.0/(r1 - param->cut);
  gsrainv1 = param->sigma_gamma * rainv1;
  gsrainvsq1 = gsrainv1*rainv1/r1; 
  expgsrainv1 = exp(gsrainv1);

  r2 = sqrt(rsq2);
  rinvsq2 = 1.0/rsq2;
  rainv2 = 1.0/(r2 - param->cut);
  gsrainv2 = param->sigma_gamma * rainv2;
  gsrainvsq2 = gsrainv2*rainv2/r2; 
  expgsrainv2 = exp(gsrainv2);

  rinv12 = 1.0/(r1*r2);
  cs = (delr1[0]*delr2[0] + delr1[1]*delr2[1] + delr1[2]*delr2[2]) * rinv12;
  delcs = cs - param->costheta;
  delcssq = delcs*delcs;

  facexp = expgsrainv1*expgsrainv2;
  facrad = param->lambda_epsilon * facexp*delcssq;
  frad1 = facrad*gsrainvsq1;
  frad2 = facrad*gsrainvsq2;
  facang = param->lambda_epsilon2 * facexp*delcs;
  facang12 = rinv12*facang;
  csfacang = cs*facang;
  csfac1 = rinvsq1*csfacang;
  
  fj[0] = delr1[0]*(frad1+csfac1)-delr2[0]*facang12;
  fj[1] = delr1[1]*(frad1+csfac1)-delr2[1]*facang12;
  fj[2] = delr1[2]*(frad1+csfac1)-delr2[2]*facang12;
  
  csfac2 = rinvsq2*csfacang;

  fk[0] = delr2[0]*(frad2+csfac2)-delr1[0]*facang12;
  fk[1] = delr2[1]*(frad2+csfac2)-delr1[1]*facang12;
  fk[2] = delr2[2]*(frad2+csfac2)-delr1[2]*facang12;

  if (eflag) eng = facrad;
}
