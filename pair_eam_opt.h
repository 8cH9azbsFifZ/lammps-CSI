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
   Contributing authors:  
     James Fischer, High Performance Technologies, Inc.
     Charles Cornwell, High Performance Technologies, Inc.
     David Richie, Stone Ridge Technology
     Vincent Natol, Stone Ridge Technology
------------------------------------------------------------------------- */

#ifndef PAIR_EAM_OPT_H
#define PAIR_EAM_OPT_H

#include "math.h"
#include "stdlib.h"
#include "pair_eam.h"
#include "atom.h"
#include "update.h"
#include "force.h"
#include "neighbor.h"
#include "comm.h"
#include "memory.h"

#define MIN(a,b) ((a) < (b) ? (a) : (b))
#define MAX(a,b) ((a) > (b) ? (a) : (b))

namespace LAMMPS_NS {

// use virtual public since this class is parent in multiple inheritance

class PairEAMOpt : virtual public PairEAM {
 public:
  PairEAMOpt(class LAMMPS *);
  virtual ~PairEAMOpt() {}
  void compute(int, int);

 private:
  template < int EFLAG, int VFLAG, int NEWTON_PAIR > void eval();
};

template < int EFLAG, int VFLAG, int NEWTON_PAIR >
void PairEAMOpt::eval()
{
  typedef struct { double x,y,z; } vec3_t;
  
  typedef struct {
    double rhor0i,rhor1i,rhor2i,rhor3i;
    double rhor0j,rhor1j,rhor2j,rhor3j;
  } fast_alpha_t;
  
  typedef struct {
    double frho0,frho1,frho2,frho3,frho4,frho5,frho6;
    double _pad[1];
  } fast_beta_t;
  
  typedef struct {
    double rhor4i,rhor5i,rhor6i;
    double rhor4j,rhor5j,rhor6j;
    double z2r0,z2r1,z2r2,z2r3,z2r4,z2r5,z2r6;
    double _pad[3];
  } fast_gamma_t;
  
  double** __restrict__ f;
  double* __restrict__ coeff;
  
  // grow energy array if necessary
  
  if (atom->nmax > nmax) {
    memory->sfree(rho);
    memory->sfree(fp);
    nmax = atom->nmax;
    rho = (double *) memory->smalloc(nmax*sizeof(double),"pair:rho");
    fp = (double *) memory->smalloc(nmax*sizeof(double),"pair:fp");
  }
  
  eng_vdwl = 0.0;
  if (VFLAG) for (int i = 0; i < 6; i++) virial[i] = 0.0;
  
  if (VFLAG == 2) f = update->f_pair;
  else f = atom->f;
  
  double** __restrict__ x = atom->x;
  int* __restrict__ type = atom->type;
  int nlocal = atom->nlocal;
  
  vec3_t* __restrict__ xx = (vec3_t*)x[0];
  vec3_t* __restrict__ ff = (vec3_t*)f[0];
  
  double tmp_cutforcesq = cutforcesq;
  double tmp_rdr = rdr;
  int nr2 = nr-2;
  int nr1 = nr-1;
  
  int** __restrict__ firstneigh = neighbor->firstneigh;
  int* __restrict__ num = neighbor->numneigh;
  
  int ntypes = atom->ntypes;
  int ntypes2 = ntypes*ntypes;
  
  fast_alpha_t* __restrict__ fast_alpha = 
    (fast_alpha_t*) malloc(ntypes2*(nr+1)*sizeof(fast_alpha_t));
  for( int i = 0; i < ntypes; i++) for( int j = 0; j < ntypes; j++) {
    fast_alpha_t* __restrict__ tab = &fast_alpha[i*ntypes*nr+j*nr];
    for(int m = 1; m <= nr; m++) {
      tab[m].rhor0i =  rhor_spline[type2rhor[i+1][j+1]][m][6];
      tab[m].rhor1i =  rhor_spline[type2rhor[i+1][j+1]][m][5];
      tab[m].rhor2i =  rhor_spline[type2rhor[i+1][j+1]][m][4];
      tab[m].rhor3i =  rhor_spline[type2rhor[i+1][j+1]][m][3];
      tab[m].rhor0j =  rhor_spline[type2rhor[j+1][i+1]][m][6];
      tab[m].rhor1j =  rhor_spline[type2rhor[j+1][i+1]][m][5];
      tab[m].rhor2j =  rhor_spline[type2rhor[j+1][i+1]][m][4];
      tab[m].rhor3j =  rhor_spline[type2rhor[j+1][i+1]][m][3];
    }
  }
  fast_alpha_t* __restrict__ tabeight = fast_alpha;
  
  fast_gamma_t* __restrict__ fast_gamma = 
    (fast_gamma_t*) malloc(ntypes2*(nr+1)*sizeof(fast_gamma_t));
  for( int i = 0; i < ntypes; i++) for( int j = 0; j < ntypes; j++) {
    fast_gamma_t* __restrict__ tab = &fast_gamma[i*ntypes*nr+j*nr];
    for(int m = 1; m <= nr; m++) {
      tab[m].rhor4i =  rhor_spline[type2rhor[i+1][j+1]][m][2];
      tab[m].rhor5i =  rhor_spline[type2rhor[i+1][j+1]][m][1];
      tab[m].rhor6i =  rhor_spline[type2rhor[i+1][j+1]][m][0];
      tab[m].rhor4j =  rhor_spline[type2rhor[j+1][i+1]][m][2];
      tab[m].rhor5j =  rhor_spline[type2rhor[j+1][i+1]][m][1];
      tab[m].rhor6j =  rhor_spline[type2rhor[j+1][i+1]][m][0];
      tab[m].z2r0 =  z2r_spline[type2z2r[i+1][j+1]][m][6];
      tab[m].z2r1 =  z2r_spline[type2z2r[i+1][j+1]][m][5];
      tab[m].z2r2 =  z2r_spline[type2z2r[i+1][j+1]][m][4];
      tab[m].z2r3 =  z2r_spline[type2z2r[i+1][j+1]][m][3];
      tab[m].z2r4 =  z2r_spline[type2z2r[i+1][j+1]][m][2];
      tab[m].z2r5 =  z2r_spline[type2z2r[i+1][j+1]][m][1];
      tab[m].z2r6 =  z2r_spline[type2z2r[i+1][j+1]][m][0];
    }
  }
  fast_gamma_t* __restrict__ tabss = fast_gamma;
  
  // zero out density
  
  if (NEWTON_PAIR) {
    int m = nlocal + atom->nghost;
    for (int i = 0; i < m; i++) rho[i] = 0.0;
  } else for (int i = 0; i < nlocal; i++) rho[i] = 0.0;
  
  // rho = density at each atom
  // loop over neighbors of my atoms
  
  for (int i = 0; i < nlocal; i++) {
    double xtmp = xx[i].x;
    double ytmp = xx[i].y;
    double ztmp = xx[i].z;
    int itype = type[i] - 1;
    int* __restrict__ neighs = firstneigh[i];
    int numneigh = num[i];
    
    double tmprho = rho[i];
    
    fast_alpha_t* __restrict__ tabeighti = &tabeight[itype*ntypes*nr];
    for (int k = 0; k < numneigh; k++) {
      int j = neighs[k];
      
      double delx = xtmp - xx[j].x;
      double dely = ytmp - xx[j].y;
      double delz = ztmp - xx[j].z;
      double rsq = delx*delx + dely*dely + delz*delz;
      
      if (rsq < tmp_cutforcesq) {
	
	int jtype = type[j] - 1;
	
	double p = sqrt(rsq)*tmp_rdr;
	if ( (int)p <= nr2 ) {
	  int m = (int)p + 1;
	  p -= (double)((int)p);
	  fast_alpha_t& a = tabeighti[jtype*nr+m];
	  tmprho += ((a.rhor3j*p+a.rhor2j)*p+a.rhor1j)*p+a.rhor0j;
	  if (NEWTON_PAIR || j < nlocal) {
	    rho[j] += ((a.rhor3i*p+a.rhor2i)*p+a.rhor1i)*p+a.rhor0i;
	  }
	} else {
	  fast_alpha_t& a = tabeighti[jtype*nr+nr1];
	  tmprho += a.rhor3j+a.rhor2j+a.rhor1j+a.rhor0j;
	  if (NEWTON_PAIR || j < nlocal) {
	    rho[j] += a.rhor3i+a.rhor2i+a.rhor1i+a.rhor0i;
	  }
	  
	}
	
      }
    }
    rho[i] = tmprho;
  }
  
  // communicate and sum densities
  
  if (NEWTON_PAIR) comm->reverse_comm_pair(this);
  
  // fp = derivative of embedding energy at each atom
  // phi = embedding energy at each atom
  
  for (int i = 0; i < nlocal; i++) {
    double p = rho[i]*rdrho;
    int m = MIN((int)p,nrho-2);
    p -= (double)m;
    ++m;
    coeff = frho_spline[type2frho[type[i]]][m];
    fp[i] = (coeff[0]*p + coeff[1])*p + coeff[2];
    if (EFLAG) eng_vdwl += ((coeff[3]*p + coeff[4])*p + coeff[5])*p + coeff[6];
  }
  
  // communicate derivative of embedding function
  
  comm->comm_pair(this);
  
  // compute forces on each atom
  // loop over neighbors of my atoms
  
  for (int i = 0; i < nlocal; i++) {
    double xtmp = xx[i].x;
    double ytmp = xx[i].y;
    double ztmp = xx[i].z;
    int itype1 = type[i] - 1;
    int* __restrict__ neighs = firstneigh[i];
    int numneigh = num[i];
    
    double tmpfx = 0.0;
    double tmpfy = 0.0;
    double tmpfz = 0.0;
    
    fast_gamma_t* __restrict__ tabssi = &tabss[itype1*ntypes*nr];
    for (int k = 0; k < numneigh; k++) {
      int j = neighs[k];
      
      double delx = xtmp - xx[j].x;
      double dely = ytmp - xx[j].y;
      double delz = ztmp - xx[j].z;
      double rsq = delx*delx + dely*dely + delz*delz;
      
      if (rsq < tmp_cutforcesq) {
	int jtype = type[j] - 1;
	double r = sqrt(rsq);
	double rhoip,rhojp,z2,z2p;
	double p = r*tmp_rdr;
	if ( (int)p <= nr2 ) {
	  int m = (int) p + 1;
	  p -= (double)((int) p);
	  
	  fast_gamma_t& a = tabssi[jtype*nr+m];
	  rhoip = (a.rhor6i*p + a.rhor5i)*p + a.rhor4i;
	  rhojp = (a.rhor6j*p + a.rhor5j)*p + a.rhor4j;
	  z2 = ((a.z2r3*p + a.z2r2)*p + a.z2r1)*p + a.z2r0;
	  z2p = (a.z2r6*p + a.z2r5)*p + a.z2r4;
	  
	} else {
	  
	  fast_gamma_t& a = tabssi[jtype*nr+nr1];
	  rhoip = a.rhor6i + a.rhor5i + a.rhor4i;
	  rhojp = a.rhor6j + a.rhor5j + a.rhor4j;
	  z2 = a.z2r3 + a.z2r2 + a.z2r1 + a.z2r0;
	  z2p = a.z2r6 + a.z2r5 + a.z2r4;
	}
	
	// rhoip = derivative of (density at atom j due to atom i)
	// rhojp = derivative of (density at atom i due to atom j)
	// phi = pair potential energy
	// phip = phi'
	// z2 = phi * r
	// z2p = (phi * r)' = (phi' r) + phi
	// psip needs both fp[i] and fp[j] terms since r_ij appears in two
	//   terms of embed eng: Fi(sum rho_ij) and Fj(sum rho_ji)
	//   hence embed' = Fi(sum rho_ij) rhojp + Fj(sum rho_ji) rhoip
	
	double recip = 1.0/r;
	double phi = z2*recip;
	double phip = z2p*recip - phi*recip;
	double psip = fp[i]*rhojp + fp[j]*rhoip + phip;
	double fforce = -psip*recip;
	
	tmpfx += delx*fforce;
	tmpfy += dely*fforce;
	tmpfz += delz*fforce;
	if (NEWTON_PAIR || j < nlocal) {
	  ff[j].x -= delx*fforce;
	  ff[j].y -= dely*fforce;
	  ff[j].z -= delz*fforce;
	}
	
	if (EFLAG) {
	  if (NEWTON_PAIR || j < nlocal) eng_vdwl += phi;
	  else eng_vdwl += 0.5*phi;
	}
	
	if (VFLAG == 1) {
	  if (NEWTON_PAIR == 0 && j >= nlocal) fforce *= 0.5;
	  virial[0] += delx*delx*fforce;
	  virial[1] += dely*dely*fforce;
	  virial[2] += delz*delz*fforce;
	  virial[3] += delx*dely*fforce;
	  virial[4] += delx*delz*fforce;
	  virial[5] += dely*delz*fforce;
	}
      }
    }
    ff[i].x += tmpfx;
    ff[i].y += tmpfy;
    ff[i].z += tmpfz;
  }
  
  free(fast_alpha); fast_alpha = 0;
  free(fast_gamma); fast_gamma = 0;
  
  if (VFLAG == 2) virial_compute();
}

}

#endif
