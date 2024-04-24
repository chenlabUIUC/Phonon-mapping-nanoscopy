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

#include <math.h>
#include "math_extra.h"
#include <stdio.h>
#include <string.h>
#include "fix_bd_body_nve.h"
#include "atom.h"
#include "comm.h"
#include "force.h"
#include "update.h"
#include "error.h"
#include "random_mars.h"

#include "atom_vec_body.h"//Shu:210818
#include "body_rounded_polygon.h"//Shu:210818

//Shu:210906
#include <iostream>
#include <fstream>
#include <math.h>
#include "update.h"
#define _DUMP_F

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixBDBodyNVE::FixBDBodyNVE(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  // if (strcmp(style,"nve/sphere") != 0 && narg <= 5)
  //   error->all(FLERR,"Illegal fix nve command");

  // t_target = force->numeric(FLERR,arg[3]); // set temperature
  // t_period = force->numeric(FLERR,arg[4]); // same as t_period in fix_langevin_overdamp.cpp
  // seed = force->inumeric(FLERR,arg[5]); //seed for random number generator. integer
  t_target = utils::numeric(FLERR,arg[3],false,lmp);
  t_period = utils::numeric(FLERR,arg[4],false,lmp);
  seed = utils::inumeric(FLERR,arg[5],false,lmp);
  //Shu:210906
  atom_tgt =  utils::inumeric(FLERR,arg[6],false,lmp);
  sprintf(paramFileName_head, arg[7]);
  DT = utils::inumeric(FLERR,arg[8],false,lmp);
  dump_time = utils::inumeric(FLERR,arg[9],false,lmp);
  //.

  if (t_target <= 0.0) error->all(FLERR,"Fix bd temperature must be > 0.0");
  if (t_period <= 0.0) error->all(FLERR,"Fix bd period must be > 0.0");
  if (seed <= 0) error->all(FLERR,"Illegal fix bd command");

  // initialize Marsaglia RNG with processor-unique seed

  random = new RanMars(lmp,seed + comm->me);

  dynamic_group_allow = 1;
  time_integrate = 1;
}

/* ---------------------------------------------------------------------- */

int FixBDBodyNVE::setmask()
{
  int mask = 0;
  mask |= INITIAL_INTEGRATE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixBDBodyNVE::init()
{
  dtv = update->dt;  // timestep
  dtf = t_period * force->ftm2v;
  gfactor = sqrt(2*force->boltz*t_target/t_period/dtv/force->mvv2e) / force->ftm2v;

  avec = (AtomVecBody *) atom->style_match("body");
}

/* ----------------------------------------------------------------------
   allow for both per-type and per-atom mass
------------------------------------------------------------------------- */

void FixBDBodyNVE::initial_integrate(int vflag)
{
  double dtfm;
  double randf;

  AtomVecBody::Bonus *bonus;
  double **torque = atom->torque;
  double **angmom = atom->angmom;
  double omega[3],tran[3];
  double *inertia, *quat; //Shu:210729
  int ibonus;
  double ex_space[3],ey_space[3],ez_space[3];

  // update v and x of atoms in group

  double **x = atom->x;
  double **v = atom->v;
  double **f = atom->f;
  double *rmass = atom->rmass;
  double *mass = atom->mass;
  int *type = atom->type;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  double fx_rand, fy_rand, fz_rand, T_rand, rand_factor;
  double x0_old, x1_old, x2_old;
  // std::cout << "~~~~~~~~~~~~~~~~~~~~~~~~~\n";
  // std::cout << f[0][0];
  // std::cout << f[0][1];
  // std::cout << f[0][2] << "\n";
  // std::cout << f[1][0];
  // std::cout << f[1][1];
  // std::cout << f[1][2] << "\n";
  // std::cout << f[2][0];
  // std::cout << f[2][1];
  // std::cout << f[2][2] << "\n";
  // std::cout << f[3][0];
  // std::cout << f[3][1];
  // std::cout << f[3][2] << "\n";
  // std::cout << "~~~~~~~~~~~~~~~~~~~~~~~~~\n";


  if (igroup == atom->firstgroup) nlocal = atom->nfirst;
  // printf("Omega = [%f, %f, %f]\n", omega[0], omega[1], omega[2] );
  // printf("torque is %f, %f \n",torque[0][2], torque[1][2] );
  // printf("rmass = %d\n",rmass );
  if (rmass) {
    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) {
        dtfm = dtf / rmass[i];
        randf = sqrt(rmass[i]) * gfactor;
        // std::cout << "dtfm: " << dtfm << std::endl;
        // std::cout << "randf: " << randf << std::endl;
        // std::cout << "~~~~~~~~~~~~~~~~~~~~~~~~~\n";
        // std::cout << f[i][0] << "   ";
        // std::cout << f[i][1] << "   ";
        // std::cout << f[i][2] << "\n";
        // // std::cout << f[1][0];
        // // std::cout << f[1][1];
        // // std::cout << f[1][2] << "\n";
        // // std::cout << f[2][0];
        // // std::cout << f[2][1];
        // // std::cout << f[2][2] << "\n";
        // // std::cout << f[3][0];
        // // std::cout << f[3][1];
        // // std::cout << f[3][2] << "\n";
        // std::cout << "~~~~~~~~~~~~~~~~~~~~~~~~~\n";
        v[i][0] = dtfm * (f[i][0]+57.8/16/sqrt(dtv)*random->gaussian()); // 57.8/16, 2.8
        v[i][1] = dtfm * (f[i][1]+57.8/16/sqrt(dtv)*random->gaussian()); // 57.8/16
        v[i][2] = dtfm * (f[i][2]+57.8/16/sqrt(dtv)*random->gaussian()); // 57.8/16

        // v[i][0] = dtfm * f[i][0];
        // v[i][1] = dtfm * f[i][1];
        // v[i][2] = dtfm * f[i][2];



        x0_old = x[i][0];
        x1_old = x[i][1];
        x2_old = x[i][2];

        x[i][0] += dtv *v[i][0];
        x[i][1] += dtv *v[i][1];
        x[i][2] += dtv *v[i][2];

        ibonus = atom->body[i];
        bonus = &avec->bonus[ibonus];
        inertia = bonus->inertia;
        // printf("inertia = [%f, %f, %f]\n", inertia[0], inertia[1], inertia[2] );
        angmom[i][2] = dtf * (torque[i][2]+ 0.0/sqrt(dtv)*random->gaussian()); // 0.7341
        // angmom[i][2] = dtf * torque[i][2];
        dtq = 0.5 * dtv;
        quat = bonus->quat;
        // printf("1. Omega = [%f, %f, %f]\n", omega[0], omega[1], omega[2] );
        MathExtra::mq_to_omega(angmom[i],quat,inertia,omega);
        // printf("1. v = [%f, %f, %f]\n", v[i][0], v[i][1], v[i][2] );
        // printf("2. Omega = [%f, %f, %f]\n", omega[0], omega[1], omega[2] );
        MathExtra::richardson(quat,angmom[i],omega,inertia,dtq);

        #ifdef _DUMP_F
        if(update->ntimestep > dump_time && update->ntimestep % DT == 0){
          // if(atom->tag[i] == atom_tgt  ){
              //open file and write fx, fy, fz;
              // double **x = atom->x;
              std::ofstream outfile;
              sprintf(paramFileName, "%s_bd_body_%d.txt",paramFileName_head,update->ntimestep);
              outfile.open(paramFileName, std::ios_base::app); // append instead of overwrite
              if ( !outfile.is_open() ) {
              std::cout << "Couldn't open file"<<paramFileName<<"\n";
              }
              outfile << atom->tag[i] << " "<<x0_old<<" "<<x1_old<< " "<<x[i][0]<<" "<<x[i][1]<<" "<<v[i][0]<<" "<<v[i][1]<<" "<<angmom[i][2] <<" "
              <<dtfm * f[i][0]<<" "<<dtfm * f[i][1]<<" "<<t_period*force->mvv2e/inertia[0] * torque[i][2]<<" "<<fx_rand<<" "<<fy_rand<<" "<<T_rand<<" "<<dtfm<<" "<<dtv<<" "<<rand_factor<<"\n";

              outfile.close();
          // }
        }
        #endif

      }

  } else {
    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) {
        dtfm = dtf / mass[type[i]];
        // randf = sqrt(mass[type[i]]) * gfactor;
        // x[i][0] += dtv * dtfm * (f[i][0]+randf*random->gaussian());
        // x[i][1] += dtv * dtfm * (f[i][1]+randf*random->gaussian());
        // x[i][2] += dtv * dtfm * (f[i][2]+randf*random->gaussian());
        // v[i][0] = dtfm * (f[i][0]+randf*random->gaussian());
        // v[i][1] = dtfm * (f[i][1]+randf*random->gaussian());
        // v[i][2] = dtfm * (f[i][2]+randf*random->gaussian());

        v[i][0] = dtfm * f[i][0];
        v[i][1] = dtfm * f[i][1];
        v[i][2] = dtfm * f[i][2];
        std::cout << "~~~~~~~~~~~~~~~~~~~~~~~~~";
        std::cout << v;

        x[i][0] += dtv *v[i][0];
        x[i][1] += dtv *v[i][1];
        x[i][2] += dtv *v[i][2];

        ibonus = atom->body[i];
        bonus = &avec->bonus[ibonus];
        // quat = bonus->quat;
        inertia = bonus->inertia;
        // angmom[i][2] = dtf * (torque[i][2]+ sqrt(inertia[0])*gfactor*random->gaussian());
        angmom[i][2] = dtf * torque[i][2];
        dtq = 0.5 * dtv;
        quat = bonus->quat;
        MathExtra::mq_to_omega(angmom[i],quat,inertia,omega);
        MathExtra::richardson(quat,angmom[i],omega,inertia,dtq);

        #ifdef _DUMP_F
        if(update->ntimestep > dump_time && update->ntimestep % DT == 0){
          // if(atom->tag[i] == atom_tgt  ){
              //open file and write fx, fy, fz;
              // double **x = atom->x;
              std::ofstream outfile;
              sprintf(paramFileName, "%s_bd_body_%d.txt",paramFileName_head,update->ntimestep);
              outfile.open(paramFileName, std::ios_base::app); // append instead of overwrite
              if ( !outfile.is_open() ) {
              std::cout << "Couldn't open file"<<paramFileName<<"\n";
              }
              outfile << atom->tag[i]<< " "<<x0_old<<" "<<x1_old << " "<<x[i][0]<<" "<<x[i][1]<<" "<<v[i][0]<<" "<<v[i][1]<<" "<<angmom[i][2] <<" "
              <<dtfm*f[i][0]<<" "<<dtfm*f[i][1]<<" "<<t_period*force->mvv2e/inertia[0] * torque[i][2]<<" "<<fx_rand<<" "<<fy_rand<<" "<<T_rand<<" "<<dtfm<<" "<<dtv<<" "<<rand_factor<<"\n";
              outfile.close();
          // }
        }
        #endif
      }
  }
}

/* ---------------------------------------------------------------------- */

void FixBDBodyNVE::reset_dt()
{
  dtv = update->dt;
  dtf = t_period * force->ftm2v;
  gfactor = sqrt(2*force->boltz*t_target/t_period/dtv/force->mvv2e) / force->ftm2v;
}
