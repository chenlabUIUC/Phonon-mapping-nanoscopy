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
   Changed from fix_nve.cpp May 14th 2016 version. Guang Shi, July 2016
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Brownian Dynamics integrator. Euler Algorithm.
------------------------------------------------------------------------- */

#ifdef FIX_CLASS

FixStyle(bd/body/nve,FixBDBodyNVE)

#else

#ifndef LMP_FIX_BD_BODY_NVE_H
#define LMP_FIX_BD_BODY_NVE_H

#include "fix.h"

namespace LAMMPS_NS {

class FixBDBodyNVE : public Fix {
 public:
  FixBDBodyNVE(class LAMMPS *, int, char **);
  virtual ~FixBDBodyNVE() {}
  int setmask();
  virtual void init();
  virtual void initial_integrate(int);
  virtual void reset_dt();

 protected:
  double t_target,t_period;
  double dtv,dtf;
  double gfactor;
  int mass_require;
  double dtq; //Shu:210825

  class RanMars *random;
  int seed;

  class AtomVecBody *avec;

  //Shu:210906
  char paramFileName[1024];
  char paramFileName_head[1024];
  int atom_tgt;
  int DT;
  int dump_time;
  //.
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

*/
