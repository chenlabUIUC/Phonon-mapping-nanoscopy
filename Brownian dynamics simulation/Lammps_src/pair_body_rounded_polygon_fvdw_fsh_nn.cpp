// clang-format off
/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
  author: Shu Wang
------------------------------------------------------------------------- */

#include "pair_body_rounded_polygon_fvdw_fsh_nn.h"

#include <cmath>
#include <cstring>
#include "math_extra.h"
#include "atom.h"
#include "atom_vec_body.h"
#include "body_rounded_polygon.h"
#include "comm.h"
#include "force.h"
#include "fix.h"
#include "modify.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "memory.h"
#include "error.h"

//Shu:180911
#include <iostream>
#include <fstream>
#include <math.h>
#include "update.h"
#include <Eigen/Dense>

// #include <iterator>
// #include <map>
// #include <tuple>

# define M_PI           3.14159265358979323846  /* pi */
// #define _DEBUG_SHU
#define MINI_THETA 1e-3
#define MINI_THETA_x_diag 1e-3
#define MINI_THETA_x1 0.0
#define MINI_THETA_x2 1.0
#define MINI_THETA_x3 4.0
#define MINI_THETA_x4 5.0
#define MINI_THETA_t 0.2
#define MINI_THETA_y 5.0
#define MINI_THETA_t_diag 1e-2

//.


using namespace LAMMPS_NS;

#define DELTA 10000
#define EPSILON 1e-3
#define MAX_CONTACTS 4  // maximum number of contacts for 2D models
#define EFF_CONTACTS 2  // effective contacts for 2D models

//#define _CONVEX_POLYGON
// #define _POLYGON_DEBUG
#define _DUMP_F
// #define _DUMP_F_more
#define DUMP_F_allpair 1

enum {INVALID=0,NONE=1,VERTEXI=2,VERTEXJ=3,EDGE=4};

/* ---------------------------------------------------------------------- */

PairBodyRoundedPolygonFvdwFshNN::PairBodyRoundedPolygonFvdwFshNN(LAMMPS *lmp) : Pair(lmp)
{
  dmax = nmax = 0;
  discrete = nullptr;
  dnum = dfirst = nullptr;

  edmax = ednummax = 0;
  edge = nullptr;
  ednum = edfirst = nullptr;

  enclosing_radius = nullptr;
  rounded_radius = nullptr;
  maxerad = nullptr;

  single_enable = 0;
  restartinfo = 0;

  c_n = 0.1;
  c_t = 0.2;
  mu = 0.0;
  delta_ua = 1.0;
}

/* ---------------------------------------------------------------------- */

PairBodyRoundedPolygonFvdwFshNN::~PairBodyRoundedPolygonFvdwFshNN()
{
  memory->destroy(discrete);
  memory->destroy(dnum);
  memory->destroy(dfirst);

  memory->destroy(edge);
  memory->destroy(ednum);
  memory->destroy(edfirst);

  memory->destroy(enclosing_radius);
  memory->destroy(rounded_radius);

  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(cutsq);

    memory->destroy(k_n);
    memory->destroy(k_na);
    memory->destroy(maxerad);
  }
}

/* ---------------------------------------------------------------------- */

void PairBodyRoundedPolygonFvdwFshNN::compute(int eflag, int vflag)
{
  int i,j,ii,jj,inum,jnum,itype,jtype;
  int ni,nj,npi,npj,ifirst,jfirst;
  int nei,nej,iefirst,jefirst;
  double xtmp,ytmp,ztmp,delx,dely,delz,evdwl;
  double rsq,r,radi,radj,k_nij,k_naij;
  double facc[3];
  int *ilist,*jlist,*numneigh,**firstneigh;

  evdwl = 0.0;
  ev_init(eflag,vflag);

  double **x = atom->x;
  double **v = atom->v;
  double **f = atom->f;
  double **torque = atom->torque;
  double **angmom = atom->angmom;
  double *radius = atom->radius;
  tagint* tag = atom->tag;
  int *body = atom->body;
  int *type = atom->type;
  int nlocal = atom->nlocal;
  int nall = nlocal + atom->nghost;
  int newton_pair = force->newton_pair;
  double fvdw_i[4],fvdw_j[4],fsh_i[4],fsh_j[4],ibody_coor[8],jbody_coor[8];

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  // grow the per-atom lists if necessary and initialize

  if (atom->nmax > nmax) {
    memory->destroy(dnum);
    memory->destroy(dfirst);
    memory->destroy(ednum);
    memory->destroy(edfirst);
    memory->destroy(enclosing_radius);
    memory->destroy(rounded_radius);
    nmax = atom->nmax;
    memory->create(dnum,nmax,"pair:dnum");
    memory->create(dfirst,nmax,"pair:dfirst");
    memory->create(ednum,nmax,"pair:ednum");
    memory->create(edfirst,nmax,"pair:edfirst");
    memory->create(enclosing_radius,nmax,"pair:enclosing_radius");
    memory->create(rounded_radius,nmax,"pair:rounded_radius");
  }

  ndiscrete = nedge = 0;
  for (i = 0; i < nall; i++)
    dnum[i] = ednum[i] = 0;

  // loop over neighbors of my atoms

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    itype = type[i];
    radi = radius[i];
    jlist = firstneigh[i];
    jnum = numneigh[i];

    if (body[i] >= 0) {
      if (dnum[i] == 0) body2space(i);
      npi = dnum[i];
      ifirst = dfirst[i];
      nei = ednum[i];
      iefirst = edfirst[i];
    }

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      j &= NEIGHMASK;

      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = 0; //ztmp - x[j][2]; //Shu:210922.
      rsq = delx*delx + dely*dely + delz*delz;

      jtype = type[j];
      radj = radius[j];

      // body/body interactions

      evdwl = 0.0;
      facc[0] = facc[1] = facc[2] = 0;

      if (body[i] < 0 || body[j] < 0) continue;

      if (dnum[j] == 0) body2space(j);
      npj = dnum[j];
      jfirst = dfirst[j];
      nej = ednum[j];
      jefirst = edfirst[j];

      k_nij = k_n[itype][jtype];
      k_naij = k_na[itype][jtype];

      // no interaction

      r = sqrt(rsq);

      //Shu:210802 ************************************************
      if(flag_fvdw > 0 || flag_fsh > 0){
        if(r < cut_global ){

          int ib=i,jb=j;
          // generate coor for ibody
          for(int id_vx = 0; id_vx < dnum[ib]; id_vx++){ //loop for each vertex2edge
            for(int id_xy = 0; id_xy < 2; id_xy++){  //loop for each dim
              ibody_coor[2*id_vx+id_xy] = discrete[ dfirst[ib]+id_vx ][ id_xy ]+ x[ib][id_xy];
              jbody_coor[2*id_vx+id_xy] = discrete[ dfirst[jb]+id_vx ][ id_xy ]+ x[jb][id_xy];
            }
          } //end of coor generate

          // printf("i=%d,j=%d,x[%d]=[%f,%f], x[%d]=[%f,%f],\n",i,j,ib,x[ib][0],x[ib][1],jb,x[jb][0],x[jb][1] );
          // compute_fvdw_fsh(fvdw_i, x[jb], x[ib], jbody_coor, ibody_coor); //force applid on i
          // f[ib][0] += fvdw_i[0];
          // f[ib][1] += fvdw_i[1];
          // torque[ib][2] += fvdw_i[2];

          // printf("finished f%d\n",i );

          compute_fvdw_fsh(fvdw_i, fsh_i, x[ib], x[jb], ibody_coor, jbody_coor); //force applid on i
          compute_fvdw_fsh(fvdw_j, fsh_j, x[jb], x[ib], jbody_coor, ibody_coor); //force applid on j
          // std::cout << "x[ib]: " << x[ib][0] << "  " << x[ib][1] << "  " << x[ib][2] << '\n';
          // std::cout << "x[jb]: " << x[jb][0] << "  " << x[jb][1] << "  " << x[jb][2] << '\n';
          // printf("before: fvdw_i=[%f,%f,%f],fvdw_j=[%f,%f,%f]\n",fvdw_i[0],fvdw_i[1],fvdw_i[2],fvdw_j[0],fvdw_j[1],fvdw_j[2]);
          // printf("before: fsh_i=[%f,%f,%f],fsh_j=[%f,%f,%f]\n",fsh_i[0],fsh_i[1],fsh_i[2],fsh_j[0],fsh_j[1],fsh_j[2]);
          fvdw_j[0] = (-fvdw_i[0]+ fvdw_j[0])*0.5;
          fvdw_j[1] = (-fvdw_i[1]+ fvdw_j[1])*0.5;
          fsh_j[0]  = (-fsh_i[0]+ fsh_j[0])*0.5;
          fsh_j[1]  = (-fsh_i[1]+ fsh_j[1])*0.5;


          fvdw_i[0] = -fvdw_j[0];
          fvdw_i[1] = -fvdw_j[1];
          fsh_i[0] = -fsh_j[0];
          fsh_i[1] = -fsh_j[1];

          // printf("after: fvdw_i=[%f,%f,%f],fvdw_j=[%f,%f,%f]\n",fvdw_i[0],fvdw_i[1],fvdw_i[2],fvdw_j[0],fvdw_j[1],fvdw_j[2]);
          // printf("after: fsh_i=[%f,%f,%f],fsh_j=[%f,%f,%f]\n",fsh_i[0],fsh_i[1],fsh_i[2],fsh_j[0],fsh_j[1],fsh_j[2]);

          f[jb][0] += fvdw_j[0]*flag_fvdw +fsh_j[0]*flag_fsh;
          f[jb][1] += fvdw_j[1]*flag_fvdw +fsh_j[1]*flag_fsh;
          torque[jb][2] += fvdw_j[2]*flag_fvdw + fsh_j[2]*flag_fsh;

          f[ib][0] += fvdw_i[0]*flag_fvdw +fsh_i[0]*flag_fsh;
          f[ib][1] += fvdw_i[1]*flag_fvdw +fsh_i[1]*flag_fsh;
          torque[ib][2] += fvdw_i[2]*flag_fvdw + fsh_i[2]*flag_fsh;

          #ifdef _DUMP_F
          if(update->ntimestep > dump_time && update->ntimestep % DT == 0){
            // if(DUMP_F_allpair || tag[i] == atom_tgt ||tag[j] == atom_tgt){
            // if(tag[i] == 38 ||tag[j] == 38 ||tag[i] == 42 ||tag[j] == 42){
            // if(tag[i] == 26 ||tag[j] == 26 ||tag[i] == 12 ||tag[j] == 12){
            //open file and write fx, fy, fz;
              std::ofstream outfile;
              sprintf(file_output, "%s_Fvdw_Fsh_%d.txt",file_output_head,update->ntimestep);
              outfile.open(file_output, std::ios_base::app); // append instead of overwrite
              if ( !outfile.is_open() ) {
              std::cout << "Couldn't open file"<<file_output<<"\n";
              }
              outfile <<tag[i]<<" "<<tag[j]<< " "
              <<fvdw_i[0] << " " <<fvdw_i[1] << " " <<fvdw_i[2] << " "
              <<fsh_i[0] << " " <<fsh_i[1] << " " <<fsh_i[2] <<" "
              <<f[ib][0] << " " <<f[ib][1] << " " <<torque[ib][2]<<"\n";

              outfile <<tag[j]<<" "<<tag[i]<< " "
              <<fvdw_j[0] << " " <<fvdw_j[1] << " " <<fvdw_j[2] << " "
              <<fsh_j[0] << " " <<fsh_j[1] << " " <<fsh_j[2]  <<" "
              <<f[jb][0] << " " <<f[jb][1] << " " <<torque[jb][2]<<"\n";
              outfile.close();
           // }
         }

          #endif
          // printf("1. torque is %f, %f \n",torque[0][2], torque[1][2] );

          // if(fvdw_i[0]>test_flarge || fvdw_i[1]>test_flarge ||fvdw_i[2]>test_flarge
          //   ||fvdw_j[0]>test_flarge || fvdw_j[1]>test_flarge ||fvdw_j[2]>test_flarge ){
          //   printf("time=%d,i=%d,j=%d,xi=[%f,%f], xj=[%f,%f],fi=[%f,%f,%f],fj=[%f,%f,%f],\n",
          //   update->ntimestep,i,j,x[ib][0],x[ib][1],x[jb][0],x[jb][1],fvdw_i[0],fvdw_i[1],fvdw_i[2],fvdw_j[0],fvdw_j[1],fvdw_j[2] );
          // }

        }

      }//end Fvdw

      //Shu.

      if (r > radi + radj + cut_inner) continue;

      if (npi == 1 && npj == 1) {
        sphere_against_sphere(i, j, delx, dely, delz, rsq,
                            k_nij, k_naij, x, v, f, evflag);
        continue;
      }

      // reset vertex and edge forces

      for (ni = 0; ni < npi; ni++) {
        discrete[ifirst+ni][3] = 0;
        discrete[ifirst+ni][4] = 0;
        discrete[ifirst+ni][5] = 0;
      }

      for (nj = 0; nj < npj; nj++) {
        discrete[jfirst+nj][3] = 0;
        discrete[jfirst+nj][4] = 0;
        discrete[jfirst+nj][5] = 0;
      }

      for (ni = 0; ni < nei; ni++) {
        edge[iefirst+ni][2] = 0;
        edge[iefirst+ni][3] = 0;
        edge[iefirst+ni][4] = 0;
      }

      for (nj = 0; nj < nej; nj++) {
        edge[jefirst+nj][2] = 0;
        edge[jefirst+nj][3] = 0;
        edge[jefirst+nj][4] = 0;
      }

      int num_contacts, done;
      double delta_a, j_a;
      Contact contact_list[MAX_CONTACTS];

      num_contacts = 0;

      // check interaction between i's vertices and j' edges

      vertex_against_edge(i, j, k_nij, k_naij, x, f, torque, tag,
                          contact_list, num_contacts, evdwl, facc);



      // check interaction between j's vertices and i' edges

      vertex_against_edge(j, i, k_nij, k_naij, x, f, torque, tag,
                          contact_list, num_contacts, evdwl, facc);
      // printf("2. torque is %f, %f \n",torque[0][2], torque[1][2] );
      // printf("num_contacts =  %d \n", num_contacts);

      if (num_contacts >= 2) {

        #ifdef _DEBUG_SHU
        if (num_contacts > 2) printf("num_contacts=%d\n",num_contacts );
        #endif
        // find the first two distinct contacts

        done = 0;
        for (int m = 0; m < num_contacts-1; m++) {
          for (int n = m+1; n < num_contacts; n++) {
            delta_a = contact_separation(contact_list[m], contact_list[n]);
            if (delta_a > 0) {
              j_a = delta_a / (EFF_CONTACTS * delta_ua);
              if (j_a < 1.0) j_a = 1.0;
              // j_a = 0;

              // scale the force at both contacts

              contact_forces(contact_list[m], j_a, x, v, angmom, f, torque,
                             evdwl, facc);
              contact_forces(contact_list[n], j_a, x, v, angmom, f, torque,
                             evdwl, facc);
              done = 1;

              #ifdef _POLYGON_DEBUG
              printf("  Two separate contacts %d and %d: delta_a = %f; j_a = %f\n",
                m, n, delta_a, j_a);
              printf("    %d: vertex %d of body %d and edge %d of body %d; "
                     "xv = %f %f %f; xe = %f %f %f\n",
                     m, contact_list[m].vertex, contact_list[m].ibody,
                     contact_list[m].edge, contact_list[m].jbody,
                     contact_list[m].xv[0], contact_list[m].xv[1],
                     contact_list[m].xv[2], contact_list[m].xe[0],
                     contact_list[m].xe[1], contact_list[m].xe[2]);
              printf("    %d: vertex %d of body %d and edge %d of body %d; "
                     "xv = %f %f %f; xe = %f %f %f\n",
                     n, contact_list[n].vertex, contact_list[n].ibody,
                     contact_list[n].edge, contact_list[n].jbody,
                     contact_list[n].xv[0], contact_list[n].xv[1],
                     contact_list[n].xv[2], contact_list[n].xe[0],
                     contact_list[n].xe[1], contact_list[n].xe[2]);
              #endif

              break;
            }
          }
          if (done == 1) break;
        }


      } else if (num_contacts == 1) {

        // if there's only one contact, it should be handled here
        // since forces/torques have not been accumulated from vertex2edge()

        contact_forces(contact_list[0], 1.0, x, v, angmom, f, torque, evdwl, facc);

        #ifdef _POLYGON_DEBUG
        printf("One contact between vertex %d of body %d and edge %d of body %d:\n",
                contact_list[0].vertex, tag[contact_list[0].ibody],
                contact_list[0].edge, tag[contact_list[0].jbody]);
        printf("xv = %f %f %f; xe = %f %f %f\n",
               contact_list[0].xv[0], contact_list[0].xv[1], contact_list[0].xv[2],
               contact_list[0].xe[0], contact_list[0].xe[1], contact_list[0].xe[2]);
        #endif
      }

      #ifdef _POLYGON_DEBUG
      int num_overlapping_contacts = 0;
      for (int m = 0; m < num_contacts-1; m++) {
        for (int n = m+1; n < num_contacts; n++) {
          double l = contact_separation(contact_list[m], contact_list[n]);
          if (l < EPSILON) num_overlapping_contacts++;
        }
      }
      printf("There are %d contacts detected, %d of which overlap.\n",
             num_contacts, num_overlapping_contacts);
      #endif

      if (evflag) ev_tally_xyz(i,j,nlocal,newton_pair,evdwl,0.0,
                               facc[0],facc[1],facc[2],delx,dely,delz);

    } // end for jj
  }

  if (vflag_fdotr) virial_fdotr_compute();
  // printf("3. torque is %f, %f \n",torque[0][2], torque[1][2] );
}

/* ----------------------------------------------------------------------
   allocate all arrays
------------------------------------------------------------------------- */

void PairBodyRoundedPolygonFvdwFshNN::allocate()
{
  allocated = 1;
  int n = atom->ntypes;

  memory->create(setflag,n+1,n+1,"pair:setflag");
  for (int i = 1; i <= n; i++)
    for (int j = i; j <= n; j++)
      setflag[i][j] = 0;

  memory->create(cutsq,n+1,n+1,"pair:cutsq");

  memory->create(k_n,n+1,n+1,"pair:k_n");
  memory->create(k_na,n+1,n+1,"pair:k_na");
  memory->create(maxerad,n+1,"pair:maxerad");
}

/* ----------------------------------------------------------------------
   global settings
------------------------------------------------------------------------- */

void PairBodyRoundedPolygonFvdwFshNN::settings(int narg, char **arg)
{
  if (narg < 8) error->all(FLERR,"Illegal pair_style command");

  c_n = utils::numeric(FLERR,arg[0],false,lmp);
  c_t = utils::numeric(FLERR,arg[1],false,lmp);
  mu = utils::numeric(FLERR,arg[2],false,lmp);
  delta_ua = utils::numeric(FLERR,arg[3],false,lmp);
  cut_inner = utils::numeric(FLERR,arg[4],false,lmp);

  if (delta_ua < 0) delta_ua = 1;

  //Shu:210802
  cut_global = utils::numeric(FLERR,arg[5],false,lmp); //cutoff for Fvdw
  flag_fvdw = utils::inumeric(FLERR,arg[8],false,lmp);
  flag_fsh = utils::inumeric(FLERR,arg[9],false,lmp);

  double x,y,Fx,Fy,T,E;
  int temp;
  // printf("-----before resize matrix\n" );
  //
  param_fvdw_m_0.resize(3,256);
  param_fvdw_v_0.resize(256);
  param_fvdw_m_2.resize(256,128);
  param_fvdw_v_2.resize(128);
  param_fvdw_m_4.resize(128,64);
  param_fvdw_v_4.resize(64);
  param_fvdw_m_6.resize(64,3);
  param_fvdw_v_6.resize(3);
  // printf("-----half resize matrix\n" );

  param_fsh_m_0.resize(3,256);
  param_fsh_v_0.resize(256);
  param_fsh_m_2.resize(256,128);
  param_fsh_v_2.resize(128);
  param_fsh_m_4.resize(128,64);
  param_fsh_v_4.resize(64);
  param_fsh_m_6.resize(64,3);
  param_fsh_v_6.resize(3);
  // printf("-----after resize matrix\n" );

  sprintf(paramFileName, "%s%d.txt",arg[6],0);
  read_param(param_fvdw_m_0, param_fvdw_v_0,paramFileName,4,256);
  printf("param_fvdw_m_0=[%f,%f,%f]\n",param_fvdw_m_0(0,0),param_fvdw_m_0(0,1),param_fvdw_m_0(1,1) );

  sprintf(paramFileName, "%s%d.txt",arg[6],2);
  read_param(param_fvdw_m_2, param_fvdw_v_2,paramFileName,257,128);
  printf("param_fvdw_m_2=[%f,%f,%f]\n",param_fvdw_m_2(0,0),param_fvdw_m_2(0,1),param_fvdw_m_2(1,1) );

  sprintf(paramFileName, "%s%d.txt",arg[6],4);
  read_param(param_fvdw_m_4, param_fvdw_v_4,paramFileName,129,64);
  printf("param_fvdw_m_4=[%f,%f,%f]\n",param_fvdw_m_4(0,0),param_fvdw_m_4(0,1),param_fvdw_m_4(1,1) );

  sprintf(paramFileName, "%s%d.txt",arg[6],6);
  read_param(param_fvdw_m_6, param_fvdw_v_6,paramFileName,65,3);

  sprintf(paramFileName, "%s%d.txt",arg[7],0);
  read_param(param_fsh_m_0, param_fsh_v_0, paramFileName, 4, 256);

  sprintf(paramFileName, "%s%d.txt",arg[7],2);
  read_param(param_fsh_m_2, param_fsh_v_2, paramFileName, 257, 128);

  sprintf(paramFileName, "%s%d.txt",arg[7],4);
  read_param(param_fsh_m_4, param_fsh_v_4, paramFileName, 129, 64);

  sprintf(paramFileName, "%s%d.txt",arg[7],6);
  read_param(param_fsh_m_6, param_fsh_v_6, paramFileName, 65, 3);
}

/* ----------------------------------------------------------------------
   set coeffs for one or more type pairs
------------------------------------------------------------------------- */

void PairBodyRoundedPolygonFvdwFshNN::coeff(int narg, char **arg)
{
  if (narg < 7 || narg > 8)
    error->all(FLERR,"Incorrect args for pair coefficients");
  if (!allocated) allocate();

  int ilo,ihi,jlo,jhi;
  utils::bounds(FLERR,arg[0],1,atom->ntypes,ilo,ihi,error);
  utils::bounds(FLERR,arg[1],1,atom->ntypes,jlo,jhi,error);

  double k_n_one = utils::numeric(FLERR,arg[2],false,lmp);
  double k_na_one = utils::numeric(FLERR,arg[3],false,lmp);

  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
    for (int j = MAX(jlo,i); j <= jhi; j++) {
      k_n[i][j] = k_n_one;
      k_na[i][j] = k_na_one;
      setflag[i][j] = 1;
      count++;
    }
  }
  //Shu:210806
  atom_tgt =  utils::inumeric(FLERR,arg[4],false,lmp);
  sprintf(file_output_head, arg[5]);
  DT = utils::inumeric(FLERR,arg[6],false,lmp);
  dump_time = utils::inumeric(FLERR,arg[7],false,lmp);

  if (count == 0) error->all(FLERR,"Incorrect args for pair coefficients");
}

/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */

void PairBodyRoundedPolygonFvdwFshNN::init_style()
{
  avec = (AtomVecBody *) atom->style_match("body");
  if (!avec)
    error->all(FLERR,"Pair body/rounded/polygon requires atom style body");
  if (strcmp(avec->bptr->style,"rounded/polygon") != 0)
    error->all(FLERR,"Pair body/rounded/polygon requires "
               "body style rounded/polygon");
  bptr = (BodyRoundedPolygon *) avec->bptr;

  if (force->newton_pair == 0)
    error->all(FLERR,"Pair style body/rounded/polygon requires "
               "newton pair on");

  if (comm->ghost_velocity == 0)
    error->all(FLERR,"Pair body/rounded/polygon requires "
               "ghost atoms store velocity");

  neighbor->request(this);

  // find the maximum enclosing radius for each atom type

  int i, itype;
  double eradi;
  int* body = atom->body;
  int* type = atom->type;
  int ntypes = atom->ntypes;
  int nlocal = atom->nlocal;

  if (atom->nmax > nmax) {
    memory->destroy(dnum);
    memory->destroy(dfirst);
    memory->destroy(ednum);
    memory->destroy(edfirst);
    memory->destroy(enclosing_radius);
    memory->destroy(rounded_radius);
    nmax = atom->nmax;
    memory->create(dnum,nmax,"pair:dnum");
    memory->create(dfirst,nmax,"pair:dfirst");
    memory->create(ednum,nmax,"pair:ednum");
    memory->create(edfirst,nmax,"pair:edfirst");
    memory->create(enclosing_radius,nmax,"pair:enclosing_radius");
    memory->create(rounded_radius,nmax,"pair:rounded_radius");
  }

  ndiscrete = nedge = 0;
  for (i = 0; i < nlocal; i++)
    dnum[i] = ednum[i] = 0;

  double *merad = nullptr;
  memory->create(merad,ntypes+1,"pair:merad");
  for (i = 1; i <= ntypes; i++)
    maxerad[i] = merad[i] = 0;

  int ipour;
  for (ipour = 0; ipour < modify->nfix; ipour++)
    if (strcmp(modify->fix[ipour]->style,"pour") == 0) break;
  if (ipour == modify->nfix) ipour = -1;

  int idep;
  for (idep = 0; idep < modify->nfix; idep++)
    if (strcmp(modify->fix[idep]->style,"deposit") == 0) break;
  if (idep == modify->nfix) idep = -1;

  for (i = 1; i <= ntypes; i++) {
    merad[i] = 0.0;
    if (ipour >= 0) {
      itype = i;
      merad[i] =
        *((double *) modify->fix[ipour]->extract("radius",itype));
    }
    if (idep >= 0) {
      itype = i;
      merad[i] =
        *((double *) modify->fix[idep]->extract("radius",itype));
    }
  }

  for (i = 0; i < nlocal; i++) {
    itype = type[i];
    if (body[i] >= 0) {
      if (dnum[i] == 0) body2space(i);
      eradi = enclosing_radius[i];
      if (eradi > merad[itype]) merad[itype] = eradi;
    } else
      merad[itype] = 0;
  }

  MPI_Allreduce(&merad[1],&maxerad[1],ntypes,MPI_DOUBLE,MPI_MAX,world);

  memory->destroy(merad);
}

/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

double PairBodyRoundedPolygonFvdwFshNN::init_one(int i, int j)
{
  k_n[j][i] = k_n[i][j];
  k_na[j][i] = k_na[i][j];

  return (maxerad[i]+maxerad[j]);
}

/* ----------------------------------------------------------------------
   convert N sub-particles in body I to space frame using current quaternion
   store sub-particle space-frame displacements from COM in discrete list
------------------------------------------------------------------------- */

void PairBodyRoundedPolygonFvdwFshNN::body2space(int i)
{
  int ibonus = atom->body[i];
  AtomVecBody::Bonus *bonus = &avec->bonus[ibonus];
  int nsub = bptr->nsub(bonus);
  double *coords = bptr->coords(bonus);
  int body_num_edges = bptr->nedges(bonus);
  double* edge_ends = bptr->edges(bonus);
  double eradius = bptr->enclosing_radius(bonus);
  double rradius = bptr->rounded_radius(bonus);

  // get the number of sub-particles (vertices)
  // and the index of the first vertex of my body in the list

  dnum[i] = nsub;
  dfirst[i] = ndiscrete;

  // grow the vertex list if necessary
  // the first 3 columns are for coords, the last 3 for forces

  if (ndiscrete + nsub > dmax) {
    dmax += DELTA;
    memory->grow(discrete,dmax,6,"pair:discrete");
  }

  double p[3][3];
  MathExtra::quat_to_mat(bonus->quat,p);

  for (int m = 0; m < nsub; m++) {
    MathExtra::matvec(p,&coords[3*m],discrete[ndiscrete]);
    discrete[ndiscrete][3] = 0;
    discrete[ndiscrete][4] = 0;
    discrete[ndiscrete][5] = 0;
    ndiscrete++;
  }

  // get the number of edges (vertices)
  // and the index of the first edge of my body in the list

  ednum[i] = body_num_edges;
  edfirst[i] = nedge;

  // grow the edge list if necessary
  // the first 2 columns are for vertex indices within body, the last 3 for forces

  if (nedge + body_num_edges > edmax) {
    edmax += DELTA;
    memory->grow(edge,edmax,5,"pair:edge");
  }

  if ((body_num_edges > 0) && (edge_ends == nullptr))
    error->one(FLERR,"Inconsistent edge data for body of atom {}",
                                 atom->tag[i]);

  for (int m = 0; m < body_num_edges; m++) {
    edge[nedge][0] = static_cast<int>(edge_ends[2*m+0]);
    edge[nedge][1] = static_cast<int>(edge_ends[2*m+1]);
    edge[nedge][2] = 0;
    edge[nedge][3] = 0;
    edge[nedge][4] = 0;
    nedge++;
  }

  enclosing_radius[i] = eradius;
  rounded_radius[i] = rradius;
}

/* ----------------------------------------------------------------------
   Interaction between two spheres with different radii
   according to the 2D model from Fraige et al.
---------------------------------------------------------------------- */

void PairBodyRoundedPolygonFvdwFshNN::sphere_against_sphere(int i, int j,
                       double delx, double dely, double delz, double rsq,
                       double k_n, double k_na, double** /*x*/, double** v,
                       double** f, int evflag)
{
  double rradi,rradj;
  double vr1,vr2,vr3,vnnr,vn1,vn2,vn3,vt1,vt2,vt3;
  double rij,rsqinv,R,fx,fy,fz,fn[3],ft[3],fpair,shift,energy;
  int nlocal = atom->nlocal;
  int newton_pair = force->newton_pair;

  rradi = rounded_radius[i];
  rradj = rounded_radius[j];

  rsqinv = 1.0/rsq;
  rij = sqrt(rsq);
  R = rij - (rradi + rradj);
  shift = k_na * cut_inner;

  energy = 0;

  if (R <= 0) {           // deformation occurs
    fpair = -k_n * R - shift;
    energy = (0.5 * k_n * R + shift) * R;
  } else if (R <= cut_inner) {   // not deforming but cohesive ranges overlap
    fpair = k_na * R - shift;
    energy = (-0.5 * k_na * R + shift) * R;
  } else fpair = 0.0;

  fx = delx*fpair/rij;
  fy = dely*fpair/rij;
  fz = delz*fpair/rij;

  if (R <= EPSILON) { // in contact

    // relative translational velocity

    vr1 = v[i][0] - v[j][0];
    vr2 = v[i][1] - v[j][1];
    vr3 = v[i][2] - v[j][2];

    // normal component

    vnnr = vr1*delx + vr2*dely + vr3*delz;
    vn1 = delx*vnnr * rsqinv;
    vn2 = dely*vnnr * rsqinv;
    vn3 = delz*vnnr * rsqinv;

    // tangential component

    vt1 = vr1 - vn1;
    vt2 = vr2 - vn2;
    vt3 = vr3 - vn3;

    // normal friction term at contact

    fn[0] = -c_n * vn1;
    fn[1] = -c_n * vn2;
    fn[2] = -c_n * vn3;

    // tangential friction term at contact,
    // excluding the tangential deformation term for now

    ft[0] = -c_t * vt1;
    ft[1] = -c_t * vt2;
    ft[2] = -c_t * vt3;

    fx += fn[0] + ft[0];
    fy += fn[1] + ft[1];
    fz += fn[2] + ft[2];
  }

  f[i][0] += fx;
  f[i][1] += fy;
  f[i][2] += fz;

  if (newton_pair || j < nlocal) {
    f[j][0] -= fx;
    f[j][1] -= fy;
    f[j][2] -= fz;
  }

  if (evflag) ev_tally_xyz(i,j,nlocal,newton_pair,
                           energy,0.0,fx,fy,fz,delx,dely,delz);
}

/* ----------------------------------------------------------------------
   Determine the interaction mode between i's vertices against j's edges

   i = atom i (body i)
   j = atom j (body j)
   x      = atoms' coordinates
   f      = atoms' forces
   torque = atoms' torques
   tag    = atoms' tags
   contact_list = list of contacts
   num_contacts = number of contacts between i's vertices and j's edges
   Return:
     interact = 0 no interaction at all
                1 there's at least one case where i's vertices interacts
                  with j's edges
---------------------------------------------------------------------- */

int PairBodyRoundedPolygonFvdwFshNN::vertex_against_edge(int i, int j,
                                                double k_n, double k_na,
                                                double** x, double** f,
                                                double** torque, tagint* tag,
                                                Contact* contact_list,
                                                int &num_contacts,
                                                double &evdwl, double* facc)
{
  int ni, npi, ifirst;
  int nj, jfirst, nej, jefirst;
  double xpi[3], xpj[3], dist, eradj, rradi, rradj;
  double fx, fy, fz, energy;
  int interact;

  npi = dnum[i];
  ifirst = dfirst[i];
  rradi = rounded_radius[i];

  jfirst = dfirst[j];
  nej = ednum[j];
  jefirst = edfirst[j];
  eradj = enclosing_radius[j];
  rradj = rounded_radius[j];

  energy = 0;
  interact = 0;

  // loop through body i's vertices

  for (ni = 0; ni < npi; ni++) {

    // convert body-fixed coordinates to space-fixed, xi

    xpi[0] = x[i][0] + discrete[ifirst+ni][0];
    xpi[1] = x[i][1] + discrete[ifirst+ni][1];
    xpi[2] = x[i][2] + discrete[ifirst+ni][2];

    // compute the distance from the vertex to the COM of body j

    distance(xpi, x[j], dist);

    #ifdef _POLYGON_DEBUG
    printf("Distance between vertex %d of body %d (%0.1f %0.1f %0.1f) "
           "to body %d's COM: %f (cut = %0.1f)\n",
           ni, xpi[0], xpi[1], xpi[2], atom->tag[i], atom->tag[j], dist,
           eradj + rradi + rradj + cut_inner);
    #endif

    // the vertex is within the enclosing circle (sphere) of body j,
    // possibly interacting

    if (dist > eradj + rradj + rradi + cut_inner) continue;

    int mode, contact, p2vertex;
    double d, R, hi[3], t, delx, dely, delz, fpair, shift;
    double rij;

    // loop through body j's edges

    for (nj = 0; nj < nej; nj++) {

      // compute the distance between the edge nj to the vertex xpi

      mode = compute_distance_to_vertex(j, nj, x[j], rradj,
                                        xpi, rradi, cut_inner,
                                        d, hi, t, contact);

      if (mode == INVALID || mode == NONE) continue;

      if (mode == VERTEXI || mode == VERTEXJ) {

        interact = 1;

        // vertex i interacts with a vertex of the edge, but does not contact

        if (mode == VERTEXI) p2vertex = edge[jefirst+nj][0];
        else if (mode == VERTEXJ) p2vertex = edge[jefirst+nj][1];

        // double xj[3];
        // p2.body2space(p2vertex, xj);
        xpj[0] = x[j][0] + discrete[jfirst+p2vertex][0];
        xpj[1] = x[j][1] + discrete[jfirst+p2vertex][1];
        xpj[2] = x[j][2] + discrete[jfirst+p2vertex][2];

        delx = xpi[0] - xpj[0];
        dely = xpi[1] - xpj[1];
        delz = 0; //xpi[2] - xpj[2]; //Shu:210922.

        // R = surface separation = rij shifted by the rounded radii
        // R = rij - (p1.rounded_radius + p2.rounded_radius);
        // note: the force is defined for R, not for rij
        // R > rc:     no interaction between vertex ni and p2vertex
        // 0 < R < rc: cohesion between vertex ni and p2vertex
        // R < 0:      deformation between vertex ni and p2vertex

        rij = sqrt(delx*delx + dely*dely + delz*delz);

        R = rij - (rradi + rradj);
        shift = k_na * cut_inner;

        // the normal frictional term -c_n * vn will be added later
        F_rep = 0; //Shu:210806
        F_att = 0; //Shu:210806

        if (R <= 0) {           // deformation occurs
          fpair = -k_n * R - shift;
          F_rep = fpair; //Shu:210806

          energy += (0.5 * k_n * R + shift) * R;
        } else if (R <= cut_inner) {   // not deforming but cohesive ranges overlap
          fpair = k_na * R - shift;
          F_att = fpair; //Shu:210806
          energy += (-0.5 * k_na * R + shift) * R;
        } else fpair = 0.0;

        fx = delx*fpair/rij;
        fy = dely*fpair/rij;
        fz = delz*fpair/rij;

        //Shu:210806
        F_rep_x = delx*F_rep/rij;
        F_rep_y = dely*F_rep/rij;

        F_att_x = delx*F_att/rij;
        F_att_y = dely*F_att/rij;

        #ifdef _POLYGON_DEBUG
        printf("  Interaction between vertex %d of %d and vertex %d of %d:",
               ni, tag[i], p2vertex, tag[j]);
        printf("    mode = %d; contact = %d; d = %f; rij = %f, t = %f\n",
               mode, contact, d, rij, t);
        printf("    R = %f; cut_inner = %f\n", R, cut_inner);
        printf("    fpair = %f\n", fpair);
        #endif

        // add forces to body i and body j directly
        // avoid double counts this pair of vertices
        // i and j can be either local or ghost atoms (bodies)
        // probably need more work here when the vertices' interaction
        // are not symmetric, e.g. j interacts with the edge
        // consisting of i but in mode = EDGE instead of VERTEX*.
        // OR, for the time being assume that the edge length is
        // sufficiently greater than the rounded radius to distinguish
        // vertex-vertex from vertex-edge contact modes.
        // Special case: when i is a sphere, also accumulate

        if (tag[i] < tag[j] || npi == 1) {

          f[i][0] += fx;
          f[i][1] += fy;
          f[i][2] += fz;
          double tx_i;
          tx_i = sum_torque(x[i], xpi, fx, fy, fz, torque[i]);

          #ifdef _DUMP_F

          if(update->ntimestep > dump_time && update->ntimestep % DT == 0){

            // printf("----time=%d,dump_time=%d,DT=%d,atom_tgt=%d\n",
            // update->ntimestep, dump_time, DT,atom_tgt);

            // if(DUMP_F_allpair || tag[i] == atom_tgt ||tag[j] == atom_tgt){
            if(tag[i] == 38 ||tag[j] == 38 ||tag[i] == 42 ||tag[j] == 42){
            // if(tag[i] == 26 ||tag[j] == 26 ||tag[i] == 12 ||tag[j] == 12){
                //open file and write fx, fy, fz;
                std::ofstream outfile;
                sprintf(file_output, "%s_nonlap_%d.txt",file_output_head,update->ntimestep);
                outfile.open(file_output, std::ios_base::app); // append instead of overwrite
                if ( !outfile.is_open() ) {
                std::cout << "Couldn't open file"<<file_output<<"\n";
              }
              // else{
              //   std::cout << "Opened:"<<file_output<<"\n";
              // }
                // outfile << tag[i] << " " << tag[j]<< " "<<x[i][0]<<" "<<x[i][1]<<" "<<x[j][0]<<" "
                // <<x[j][1]<<" "<< fx<<" "<<fy<<" "<<tx_i <<" "<<f[i][0]<<" "<<f[i][1]<<" "<<torque[i][2]<<" "
                // << F_att_x<<" "<<F_att_y<<" "<< F_rep_x<<" "<<F_rep_y<<"\n";
                outfile << tag[i] << " " << tag[j]<< " "<< fx<<" "<<fy<<" "<<tx_i <<" "
                <<f[i][0]<<" "<<f[i][1]<<" "<<torque[i][2]<<" "
                <<x[i][0]<<" "<<x[i][1]<<" "
                // << F_att_x<<" "<<F_att_y<<" "<< F_rep_x<<" "<<F_rep_y<<"\n";
                << F_att_x<<" "<<F_att_y<<" "<< F_rep_x<<" "<<F_rep_y<<" "
                << F_rep<<" "<<F_att<<" "<< delx<<" "<<dely<<" "<<delz<<" "
                << R<<" "<<rij<<" "<< rradi<<" "<<rradj<<" "
                <<"\n";
                outfile.close();
            }
          }

          #endif

          f[j][0] -= fx;
          f[j][1] -= fy;
          f[j][2] -= fz;
          double tx_j;
          tx_j = sum_torque(x[j], xpj, -fx, -fy, -fz, torque[j]);

          #ifdef _DUMP_F

          if(update->ntimestep > dump_time && update->ntimestep % DT == 0){
            // if(DUMP_F_allpair || tag[i] == atom_tgt ||tag[j] == atom_tgt  ){
            if(tag[i] == 38 ||tag[j] == 38 ||tag[i] == 42 ||tag[j] == 42){
            // if(tag[i] == 26 ||tag[j] == 26 ||tag[i] == 12 ||tag[j] == 12){
                //open file and write fx, fy, fz;
                std::ofstream outfile;
                sprintf(file_output, "%s_nonlap_%d.txt",file_output_head,update->ntimestep);
                outfile.open(file_output, std::ios_base::app); // append instead of overwrite
                if ( !outfile.is_open() ) {
                std::cout << "Couldn't open file"<<file_output<<"\n";
                }
                // outfile << tag[j] << " " << tag[i]<< " "<<x[j][0]<<" "<<x[j][1]<<" "<<x[i][0]<<" "
                // <<x[i][1]<<" "<< -fx<<" "<<-fy<<" "<<tx_j <<" "<<f[j][0]<<" "<<f[j][1]<<" "<<torque[j][2]<<" "
                // << -F_att_x<<" "<<-F_att_y<<" "<< -F_rep_x<<" "<<-F_rep_y<<"\n";
                outfile << tag[j] << " " << tag[i]<< " "<<-fx<<" "<<-fy<<" "<<tx_j <<" "
                <<f[j][0]<<" "<<f[j][1]<<" "<<torque[j][2]<<" "
                <<x[j][0]<<" "<<x[j][1]<<" "
                // << -F_att_x<<" "<<-F_att_y<<" "<< -F_rep_x<<" "<<-F_rep_y<<"\n";
                << -F_att_x<<" "<<-F_att_y<<" "<< -F_rep_x<<" "<<-F_rep_y<<" "
                << F_rep<<" "<<F_att<<" "<< delx<<" "<<dely<<" "<<delz<<" "
                << R<<" "<<rij<<" "<< rradi<<" "<<rradj<<" "
                <<"\n";

                outfile.close();
            }
          }
          #endif

          facc[0] += fx; facc[1] += fy; facc[2] += fz;

          #ifdef _POLYGON_DEBUG
          printf("    from vertex-vertex: "
                 "force on vertex %d of body %d: fx %f fy %f fz %f\n"
                 "      torque body %d: %f %f %f\n"
                 "      torque body %d: %f %f %f\n", ni, tag[i], fx, fy, fz,
            tag[i],torque[i][0],torque[i][1],torque[i][2],
            tag[j],torque[j][0],torque[j][1],torque[j][2]);
        #endif
        }

        #ifdef _CONVEX_POLYGON
        // done with the edges from body j,
        // given that vertex ni interacts with only one vertex
        //   from one edge of body j
        break;
        #endif

      } else if (mode == EDGE) {

        interact = 1;

        // vertex i interacts with the edge

        delx = xpi[0] - hi[0];
        dely = xpi[1] - hi[1];
        delz = 0; //xpi[2] - hi[2]; //Shu:210922.

        // R = surface separation = d shifted by the rounded radii
        // R = d - (p1.rounded_radius + p2.rounded_radius);
        // Note: the force is defined for R, not for d
        // R > rc:     no interaction between vertex i and edge j
        // 0 < R < rc: cohesion between vertex i and edge j
        // R < 0:      deformation between vertex i and edge j
        // rij = sqrt(delx*delx + dely*dely + delz*delz);

        R = d - (rradi + rradj);
        shift = k_na * cut_inner;

        // the normal frictional term -c_n * vn will be added later
        F_rep = 0; //Shu:210806
        F_att = 0; //Shu:210806

        if (R <= 0) {           // deformation occurs
          fpair = -k_n * R - shift;
          F_rep = fpair;
          energy += (0.5 * k_n * R + shift) * R;
        } else if (R <= cut_inner) {   // not deforming but cohesive ranges overlap
          fpair = k_na * R - shift;
          F_att = fpair;
          energy += (-0.5 * k_na * R + shift) * R;
        } else fpair = 0.0;

        fx = delx*fpair/d;
        fy = dely*fpair/d;
        fz = delz*fpair/d;

        F_rep_x = delx*F_rep/d;
        F_rep_y = dely*F_rep/d;

        F_att_x = delx*F_att/d;
        F_att_y = dely*F_att/d;

        #ifdef _POLYGON_DEBUG
        printf("  Interaction between vertex %d of %d and edge %d of %d:",
               ni, tag[i], nj, tag[j]);
        printf("    mode = %d; contact = %d; d = %f; t = %f\n",
               mode, contact, d, t);
        printf("    R = %f; cut_inner = %f\n", R, cut_inner);
        printf("    fpair = %f\n", fpair);
        #endif

        if (contact == 1) {

          // vertex ni of body i contacts with edge nj of body j

          contact_list[num_contacts].ibody = i;
          contact_list[num_contacts].jbody = j;
          contact_list[num_contacts].vertex = ni;
          contact_list[num_contacts].edge = nj;
          contact_list[num_contacts].xv[0] = xpi[0];
          contact_list[num_contacts].xv[1] = xpi[1];
          contact_list[num_contacts].xv[2] = xpi[2];
          contact_list[num_contacts].xe[0] = hi[0];
          contact_list[num_contacts].xe[1] = hi[1];
          contact_list[num_contacts].xe[2] = hi[2];
          contact_list[num_contacts].separation = R;
          num_contacts++;

          // store forces to vertex ni and the edge nj
          // to be rescaled later

          discrete[ifirst+ni][3] = fx;
          discrete[ifirst+ni][4] = fy;
          discrete[ifirst+ni][5] = fz;

          edge[jefirst+nj][2] = -fx;
          edge[jefirst+nj][3] = -fy;
          edge[jefirst+nj][4] = -fz;

          #ifdef _POLYGON_DEBUG
          printf("  Stored forces at vertex and edge for accumulating later.\n");
          #endif

        } else { // no contact

          // accumulate force and torque to both bodies directly

          f[i][0] += fx;
          f[i][1] += fy;
          f[i][2] += fz;
          double tx_i;
          tx_i = sum_torque(x[i], xpi, fx, fy, fz, torque[i]);
          #ifdef _DUMP_F
          if(update->ntimestep > dump_time && update->ntimestep % DT == 0){
            // if(DUMP_F_allpair || atom->tag[i] == atom_tgt ||atom->tag[j] == atom_tgt  ){
            if(atom->tag[i] == 38 ||atom->tag[j] == 38 ||atom->tag[i] == 42 ||atom->tag[j] == 42){
            // if(atom->tag[i] == 26 ||atom->tag[j] == 26 ||atom->tag[i] == 12 ||atom->tag[j] == 12){
                //open file and write fx, fy, fz;
                std::ofstream outfile;
                sprintf(file_output, "%s_nonlap_%d.txt",file_output_head,update->ntimestep);
                outfile.open(file_output, std::ios_base::app); // append instead of overwrite
                if ( !outfile.is_open() ) {
                std::cout << "Couldn't open file"<<file_output<<"\n";
              }
                outfile << atom->tag[i] << " " << atom->tag[j]<< " "<< fx<<" "<<fy<<" "<<tx_i <<" "
                <<f[i][0]<<" "<<f[i][1]<<" "<<torque[i][2]<<" "
                <<x[i][0]<<" "<<x[i][1]<<" "
                << F_att_x<<" "<<F_att_y<<" "<< F_rep_x<<" "<<F_rep_y<<" "
                << F_rep<<" "<<F_att<<" "<< delx<<" "<<dely<<" "<<delz<<" "
                << R<<" "<<rij<<" "<< rradi<<" "<<rradj<<" "
                <<"\n";
                outfile.close();
            }
          }
          #endif

          f[j][0] -= fx;
          f[j][1] -= fy;
          f[j][2] -= fz;
          double tx_j;
          tx_j = sum_torque(x[j], hi, -fx, -fy, -fz, torque[j]);
          #ifdef _DUMP_F
          if(update->ntimestep > dump_time && update->ntimestep % DT == 0){
            // if(DUMP_F_allpair || atom->tag[i] == atom_tgt ||atom->tag[j] == atom_tgt  ){
            if(atom->tag[i] == 38 ||atom->tag[j] == 38 ||atom->tag[i] == 42 ||atom->tag[j] == 42){
            // if(atom->tag[i] == 26 ||atom->tag[j] == 26 ||atom->tag[i] == 12 ||atom->tag[j] == 12){
                //open file and write fx, fy, fz;
                std::ofstream outfile;
                sprintf(file_output, "%s_nonlap_%d.txt",file_output_head,update->ntimestep);
                outfile.open(file_output, std::ios_base::app); // append instead of overwrite
                if ( !outfile.is_open() ) {
                std::cout << "Couldn't open file"<<file_output<<"\n";
                }
                outfile << atom->tag[j] << " " << atom->tag[i]<< " "<< -fx <<" "<<-fy <<" "<<tx_j <<" "
                <<f[j][0]<<" "<<f[j][1]<<" "<<torque[j][2]<<" "
                <<x[j][0]<<" "<<x[j][1]<<" "
                // << -F_att_x<<" "<<-F_att_y<<" "<< -F_rep_x<<" "<<-F_rep_y<<"\n";
                << -F_att_x<<" "<<-F_att_y<<" "<< -F_rep_x<<" "<< -F_rep_y<<" "
                << F_rep<<" "<<F_att<<" "<< delx<<" "<<dely<<" "<<delz<<" "
                << R<<" "<<rij<<" "<< rradi<<" "<<rradj<<" "
                <<"\n";
                outfile.close();
            }
          }
          #endif

          facc[0] += fx; facc[1] += fy; facc[2] += fz;

          #ifdef _POLYGON_DEBUG
          printf("    from vertex-edge, no contact: "
                 "force on vertex %d of body %d: fx %f fy %f fz %f\n"
                 "      torque body %d: %f %f %f\n"
                 "      torque body %d: %f %f %f\n", ni, tag[i], fx, fy, fz,
                 tag[i],torque[i][0],torque[i][1],torque[i][2],
                 tag[j],torque[j][0],torque[j][1],torque[j][2]);
          #endif
        } // end if contact

        #ifdef _CONVEX_POLYGON
        // done with the edges from body j,
        // given that vertex ni interacts with only one edge from body j
        break;
        #endif
      } // end if mode

    } // end for looping through the edges of body j

  } // end for looping through the vertices of body i

  evdwl += energy;

  return interact;
}

/* -------------------------------------------------------------------------
  Compute the distance between an edge of body i and a vertex from
  another body
  Input:
    ibody      = body i (i.e. atom i)
    edge_index = edge index of body i
    xmi        = atom i's coordinates (body i's center of mass)
    x0         = coordinate of the tested vertex from another body
    x0_rounded_radius = rounded radius of the tested vertex
    cut_inner  = cutoff for vertex-vertex and vertex-edge interaction
  Output:
    d          = Distance from a point x0 to an edge
    hi         = coordinates of the projection of x0 on the edge
    t          = ratio to determine the relative position of hi
                 wrt xi and xj on the segment
  contact      = 0 no contact between the queried vertex and the edge
                 1 contact detected
  return
    INVALID if the edge index is invalid
    NONE    if there is no interaction
    VERTEXI if the tested vertex interacts with the first vertex of the edge
    VERTEXJ if the tested vertex interacts with the second vertex of the edge
    EDGE    if the tested vertex interacts with the edge
------------------------------------------------------------------------- */

int PairBodyRoundedPolygonFvdwFshNN::compute_distance_to_vertex(int ibody,
                                                int edge_index,
                                                double *xmi,
                                                double rounded_radius,
                                                double* x0,
                                                double x0_rounded_radius,
                                                double cut_inner,
                                                double &d,
                                                double hi[3],
                                                double &t,
                                                int &contact)
{
  if (edge_index >= ednum[ibody]) return INVALID;

  int mode,ifirst,iefirst,npi1,npi2;
  double xi1[3],xi2[3],u[3],v[3],uij[3];
  double udotv, magv, magucostheta;
  double delx,dely,delz;

  ifirst = dfirst[ibody];
  iefirst = edfirst[ibody];
  npi1 = static_cast<int>(edge[iefirst+edge_index][0]);
  npi2 = static_cast<int>(edge[iefirst+edge_index][1]);

  // compute the space-fixed coordinates for the vertices of the edge

  xi1[0] = xmi[0] + discrete[ifirst+npi1][0];
  xi1[1] = xmi[1] + discrete[ifirst+npi1][1];
  xi1[2] = xmi[2] + discrete[ifirst+npi1][2];

  xi2[0] = xmi[0] + discrete[ifirst+npi2][0];
  xi2[1] = xmi[1] + discrete[ifirst+npi2][1];
  xi2[2] = xmi[2] + discrete[ifirst+npi2][2];

  // u = x0 - xi1

  u[0] = x0[0] - xi1[0];
  u[1] = x0[1] - xi1[1];
  u[2] = x0[2] - xi1[2];

  // v = xi2 - xi1

  v[0] = xi2[0] - xi1[0];
  v[1] = xi2[1] - xi1[1];
  v[2] = xi2[2] - xi1[2];

  // dot product between u and v = magu * magv * costheta

  udotv = u[0] * v[0] + u[1] * v[1] + u[2] * v[2];
  magv = sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
  magucostheta = udotv / magv;

  // uij is the unit vector pointing from xi to xj

  uij[0] = v[0] / magv;
  uij[1] = v[1] / magv;
  uij[2] = v[2] / magv;

  // position of the projection of x0 on the line (xi, xj)

  hi[0] = xi1[0] + magucostheta * uij[0];
  hi[1] = xi1[1] + magucostheta * uij[1];
  hi[2] = xi1[2] + magucostheta * uij[2];

  // distance from x0 to the line (xi, xj) = distance from x0 to hi

  distance(hi, x0, d);

  // determine the interaction mode
  // for 2D: a vertex can interact with one edge at most
  // for 3D: a vertex can interact with one face at most

  mode = NONE;
  contact = 0;

  if (d > rounded_radius + x0_rounded_radius + cut_inner) {

    // if the vertex is far away from the edge

    mode = NONE;

  } else {

    // check if x0 (the queried vertex) and xmi (the body's center of mass)
    // are on the different sides of the edge

    #ifdef _CONVEX_POLYGON
    int m = opposite_sides(xi1, xi2, x0, xmi);
    #else
    int m = 1;
    #endif

    if (m == 0) {

      // x0 and xmi are on not the opposite sides of the edge
      // leave xpi for another edge to detect

      mode = NONE;

    } else {

      // x0 and xmi are on the different sides
      // t is the ratio to detect if x0 is closer to the vertices xi or xj

      if (fabs(xi2[0] - xi1[0]) > EPSILON)
        t = (hi[0] - xi1[0]) / (xi2[0] - xi1[0]);
      else if (fabs(xi2[1] - xi1[1]) > EPSILON)
        t = (hi[1] - xi1[1]) / (xi2[1] - xi1[1]);
      else if (fabs(xi2[2] - xi1[2]) > EPSILON)
        t = (hi[2] - xi1[2]) / (xi2[2] - xi1[2]);

      double contact_dist = rounded_radius + x0_rounded_radius;
      if (t >= 0 && t <= 1) {
        mode = EDGE;
        if (d < contact_dist + EPSILON)
          contact = 1;

      } else { // t < 0 || t > 1: closer to either vertices of the edge

        if (t < 0) {
          // measure the distance from x0 to xi1
          delx = x0[0] - xi1[0];
          dely = x0[1] - xi1[1];
          delz = 0; //x0[2] - xi1[2]; //Shu:210922.
          double dx0xi1 = sqrt(delx*delx + dely*dely + delz*delz);

          if (dx0xi1 > contact_dist + cut_inner)
            mode = NONE;
          else
            mode = VERTEXI;
        } else {
          // measure the distance from x0 to xi2
          delx = x0[0] - xi2[0];
          dely = x0[1] - xi2[1];
          delz = 0; //x0[2] - xi2[2]; //Shu:210922.
          double dx0xi2 = sqrt(delx*delx + dely*dely + delz*delz);
          if (dx0xi2 > contact_dist + cut_inner)
            mode = NONE;
          else
            mode = VERTEXJ;
        }
      } // end if t >= 0 && t <= 1
    } // end if x0 and xmi are on the same side of the edge
  }

  return mode;
}

/* ----------------------------------------------------------------------
  Compute contact forces between two bodies
  modify the force stored at the vertex and edge in contact by j_a
  sum forces and torque to the corresponding bodies
  fn = normal friction component
  ft = tangential friction component (-c_t * v_t)
------------------------------------------------------------------------- */

void PairBodyRoundedPolygonFvdwFshNN::contact_forces(Contact& contact, double j_a,
                       double** x, double** v, double** angmom, double** f,
                       double** torque, double &/*evdwl*/, double* facc)
{
  int ibody,jbody,ibonus,jbonus,ifirst,jefirst,ni,nj;
  double fx,fy,fz,delx,dely,delz,rsq,rsqinv;
  double vr1,vr2,vr3,vnnr,vn1,vn2,vn3,vt1,vt2,vt3;
  double fn[3],ft[3],vi[3],vj[3];
  double *quat, *inertia;
  AtomVecBody::Bonus *bonus;

  ibody = contact.ibody;
  jbody = contact.jbody;

  // compute the velocity of the vertex in the space-fixed frame

  ibonus = atom->body[ibody];
  bonus = &avec->bonus[ibonus];
  quat = bonus->quat;
  inertia = bonus->inertia;
  total_velocity(contact.xv, x[ibody], v[ibody], angmom[ibody],
                 inertia, quat, vi);

  // compute the velocity of the point on the edge in the space-fixed frame

  jbonus = atom->body[jbody];
  bonus = &avec->bonus[jbonus];
  quat = bonus->quat;
  inertia = bonus->inertia;
  total_velocity(contact.xe, x[jbody], v[jbody], angmom[jbody],
                 inertia, quat, vj);

  // vector pointing from the vertex to the point on the edge

  delx = contact.xv[0] - contact.xe[0];
  dely = contact.xv[1] - contact.xe[1];
  delz = 0; //contact.xv[2] - contact.xe[2]; //Shu:210922.
  rsq = delx*delx + dely*dely + delz*delz;

  rsqinv = 1.0/rsq;

  // relative translational velocity

  vr1 = vi[0] - vj[0];
  vr2 = vi[1] - vj[1];
  vr3 = vi[2] - vj[2];

  // normal component

  vnnr = vr1*delx + vr2*dely + vr3*delz;
  vn1 = delx*vnnr * rsqinv;
  vn2 = dely*vnnr * rsqinv;
  vn3 = delz*vnnr * rsqinv;

  // tangential component

  vt1 = vr1 - vn1;
  vt2 = vr2 - vn2;
  vt3 = vr3 - vn3;

  // normal friction term at contact

  fn[0] = -c_n * vn1;
  fn[1] = -c_n * vn2;
  fn[2] = -c_n * vn3;

  // tangential friction term at contact
  // excluding the tangential deformation term for now

  ft[0] = -c_t * vt1;
  ft[1] = -c_t * vt2;
  ft[2] = -c_t * vt3;

  // only the cohesive force is scaled by j_a
  // mu * fne = tangential friction deformation during gross sliding
  // see Eq. 4, Fraige et al.

  ifirst = dfirst[ibody];
  ni = contact.vertex;

  fx = discrete[ifirst+ni][3] * j_a + fn[0] + ft[0] +
    mu * discrete[ifirst+ni][3];
  fy = discrete[ifirst+ni][4] * j_a + fn[1] + ft[1] +
    mu * discrete[ifirst+ni][4];
  fz = discrete[ifirst+ni][5] * j_a + fn[2] + ft[2] +
    mu * discrete[ifirst+ni][5];
  f[ibody][0] += fx;
  f[ibody][1] += fy;
  f[ibody][2] += fz;
  double tx_i;
  tx_i = sum_torque(x[ibody], contact.xv, fx, fy, fz, torque[ibody]);
  #ifdef _DUMP_F
  if(update->ntimestep > dump_time && update->ntimestep % DT == 0){
    // if(DUMP_F_allpair || atom->tag[ibody] == atom_tgt ||atom->tag[jbody] == atom_tgt  ){
    if(atom->tag[ibody] == 38 ||atom->tag[jbody] == 38 ||atom->tag[ibody] == 42 ||atom->tag[jbody] == 42){
    // if(atom->tag[ibody] == 26 ||atom->tag[jbody] == 26 ||atom->tag[ibody] == 12 ||atom->tag[jbody] == 12){
        //open file and write fx, fy, fz;
        std::ofstream outfile;
        sprintf(file_output, "%s_Fmu_%d.txt",file_output_head,update->ntimestep);
        outfile.open(file_output, std::ios_base::app); // append instead of overwrite
        if ( !outfile.is_open() ) {
        std::cout << "Couldn't open file"<<file_output<<"\n";
        }
        outfile << atom->tag[ibody] << " " << atom->tag[jbody]<< " "
        << fx<<" "<<fy<<" "<<tx_i <<" "
        <<f[ibody][0]<<" "<<f[ibody][1]<<" "<<torque[ibody][2]<<" "
        <<x[ibody][0]<<" "<<x[ibody][1]<<" "
        <<j_a <<" "<< mu<<" "
        <<discrete[ifirst+ni][3]<<" "<<fn[0]<<" "<<ft[0]<<" "
        <<discrete[ifirst+ni][4]<<" "<<fn[1]<<" "<<ft[1]<<" "
        <<discrete[ifirst+ni][5]<<" "<<fn[2]<<" "<<ft[2] <<"\n";
        outfile.close();
    }
  }
  #endif

  // accumulate forces to the vertex only

  facc[0] += fx; facc[1] += fy; facc[2] += fz;

  // only the cohesive force is scaled by j_a
  // mu * fne = tangential friction deformation during gross sliding
  // Eq. 4, Fraige et al.

  jefirst = edfirst[jbody];
  nj = contact.edge;

  fx = edge[jefirst+nj][2] * j_a - fn[0] - ft[0] +
    mu * edge[jefirst+nj][2];
  fy = edge[jefirst+nj][3] * j_a - fn[1] - ft[1] +
    mu * edge[jefirst+nj][3];
  fz = edge[jefirst+nj][4] * j_a - fn[2] - ft[2] +
    mu * edge[jefirst+nj][4];
  f[jbody][0] += fx;
  f[jbody][1] += fy;
  f[jbody][2] += fz;
  double tx_j;
  tx_j = sum_torque(x[jbody], contact.xe, fx, fy, fz, torque[jbody]);

  #ifdef _DUMP_F
  if(update->ntimestep > dump_time && update->ntimestep % DT == 0){
    // if(DUMP_F_allpair || atom->tag[ibody] == atom_tgt ||atom->tag[jbody] == atom_tgt  ){
    if(atom->tag[ibody] == 38 ||atom->tag[jbody] == 38 ||atom->tag[ibody] == 42 ||atom->tag[jbody] == 42){
    // if(atom->tag[ibody] == 26 ||atom->tag[jbody] == 26 ||atom->tag[ibody] == 12 ||atom->tag[jbody] == 12){
        //open file and write fx, fy, fz;
        std::ofstream outfile;
        sprintf(file_output, "%s_Fmu_%d.txt",file_output_head,update->ntimestep);
        outfile.open(file_output, std::ios_base::app); // append instead of overwrite
        if ( !outfile.is_open() ) {
        std::cout << "Couldn't open file"<<file_output<<"\n";
        }
        outfile << atom->tag[jbody] << " " << atom->tag[ibody]<< " "
        << fx<<" "<<fy<<" "<<tx_j <<" "
        <<f[jbody][0]<<" "<<f[jbody][1]<<" "<<torque[jbody][2]<<" "
        <<x[jbody][0]<<" "<<x[jbody][1]<<" "
        <<j_a <<" "<< mu<<" "
        <<edge[jefirst+nj][2]<<" "<<fn[0]<<" "<<ft[0]<<" "
        <<edge[jefirst+nj][3]<<" "<<fn[1]<<" "<<ft[1]<<" "
        <<edge[jefirst+nj][4]<<" "<<fn[2]<<" "<<ft[2] <<"\n";
        outfile.close();
    }
  }
  #endif

  #ifdef _POLYGON_DEBUG
  printf("From contact forces: vertex fx %f fy %f fz %f\n"
         "      torque body %d: %f %f %f\n"
         "      torque body %d: %f %f %f\n",
         discrete[ifirst+ni][3], discrete[ifirst+ni][4], discrete[ifirst+ni][5],
         atom->tag[ibody],torque[ibody][0],torque[ibody][1],torque[ibody][2],
         atom->tag[jbody],torque[jbody][0],torque[jbody][1],torque[jbody][2]);
  #endif
}

/* ----------------------------------------------------------------------
  Determine the length of the contact segment, i.e. the separation between
  2 contacts, should be extended for 3D models.
------------------------------------------------------------------------- */

double PairBodyRoundedPolygonFvdwFshNN::contact_separation(const Contact& c1,
                                                  const Contact& c2)
{
  double x1 = c1.xv[0];
  double y1 = c1.xv[1];
  double x2 = c1.xe[0];
  double y2 = c1.xe[1];
  double x3 = c2.xv[0];
  double y3 = c2.xv[1];

  double delta_a = 0.0;
  if (fabs(x2 - x1) > EPSILON) {
    double A = (y2 - y1) / (x2 - x1);
    delta_a = fabs(y1 - A * x1 - y3 + A * x3) / sqrt(1 + A * A);
  } else {
    delta_a = fabs(x1 - x3);
  }

  return delta_a;
}

/* ----------------------------------------------------------------------
  Accumulate torque to body from the force f=(fx,fy,fz) acting at point x
------------------------------------------------------------------------- */

double PairBodyRoundedPolygonFvdwFshNN::sum_torque(double* xm, double *x, double fx,
                                        double fy, double fz, double* torque)
{
  double rx = x[0] - xm[0];
  double ry = x[1] - xm[1];
  double rz = x[2] - xm[2];
  double tx = ry * fz - rz * fy;
  double ty = rz * fx - rx * fz;
  double tz = rx * fy - ry * fx;
  torque[0] += tx;
  torque[1] += ty;
  torque[2] += tz;
  return tz; //Shu:210806
}

/* ----------------------------------------------------------------------
  Test if two points a and b are in opposite sides of the line that
  connects two points x1 and x2
------------------------------------------------------------------------- */

int PairBodyRoundedPolygonFvdwFshNN::opposite_sides(double* x1, double* x2,
                                           double* a, double* b)
{
  double m_a = (x1[1] - x2[1])*(a[0] - x1[0]) + (x2[0] - x1[0])*(a[1] - x1[1]);
  double m_b = (x1[1] - x2[1])*(b[0] - x1[0]) + (x2[0] - x1[0])*(b[1] - x1[1]);
  // equal to zero when either a or b is inline with the line x1-x2
  if (m_a * m_b <= 0)
    return 1;
  else
    return 0;
}

/* ----------------------------------------------------------------------
  Calculate the total velocity of a point (vertex, a point on an edge):
    vi = vcm + omega ^ (p - xcm)
------------------------------------------------------------------------- */

void PairBodyRoundedPolygonFvdwFshNN::total_velocity(double* p, double *xcm,
                              double* vcm, double *angmom, double *inertia,
                              double *quat, double* vi)
{
  double r[3],omega[3],ex_space[3],ey_space[3],ez_space[3];
  r[0] = p[0] - xcm[0];
  r[1] = p[1] - xcm[1];
  r[2] = p[2] - xcm[2];
  MathExtra::q_to_exyz(quat,ex_space,ey_space,ez_space);
  MathExtra::angmom_to_omega(angmom,ex_space,ey_space,ez_space,
                             inertia,omega);
  vi[0] = omega[1]*r[2] - omega[2]*r[1] + vcm[0];
  vi[1] = omega[2]*r[0] - omega[0]*r[2] + vcm[1];
  vi[2] = omega[0]*r[1] - omega[1]*r[0] + vcm[2];
}

/* ---------------------------------------------------------------------- */

void PairBodyRoundedPolygonFvdwFshNN::distance(const double* x2, const double* x1,
                                      double& r)
{
  r = sqrt((x2[0] - x1[0]) * (x2[0] - x1[0])
    + (x2[1] - x1[1]) * (x2[1] - x1[1])
    + (x2[2] - x1[2]) * (x2[2] - x1[2]));
}

/* ---------------------------------------------------------------------- */
void PairBodyRoundedPolygonFvdwFshNN::compute_fvdw_fsh(double * fvdw_j, double * fsh_j, double *xi, double *xj,
      double * ibody_coor, double * jbody_coor) //Shu:210802
{
  //not distinguish itype and jtype now
  // ibody coor -> X1+theta
  // printf("in compute_fvdw_fsh_fsh:1451\n");
  double dx, dy, theta_i, theta_j;
  int i_low,j_low;
  double yi_low, yj_low;
  yi_low = ibody_coor[1]; i_low = 0;
  yj_low = jbody_coor[1]; j_low = 0;
  for(int t = 1; t < 4; t++){
    if(ibody_coor[1+2*t] < yi_low) { i_low = t; yi_low = ibody_coor[1+2*t];}
    if(jbody_coor[1+2*t] < yj_low) { j_low = t; yj_low = jbody_coor[1+2*t];}
  }
  if(i_low == 0){
    dy = ibody_coor[i_low*2+1] - ibody_coor[i_low*2+7];
    dx = ibody_coor[i_low*2]   - ibody_coor[i_low*2+6];
  }else{
    dy = ibody_coor[i_low*2+1] - ibody_coor[i_low*2-1];
    dx = ibody_coor[i_low*2]   - ibody_coor[i_low*2-2];
  }
  // printf("i: dx = %f, dy=%f\n", dx, dy);
  if(dx < 0.001) theta_i = atan2(dy,dx);
  else theta_i = M_PI*0.5;
  // theta_i = atan2(dy, dx);

  // printf("compute_fvdw_fsh:1471 theta_i=%f\n", theta_i );
  if(j_low == 0){
    dy = jbody_coor[j_low*2+1] - jbody_coor[j_low*2+7];
    dx = jbody_coor[j_low*2] - jbody_coor[j_low*2+6];
  }else{
    dy = jbody_coor[j_low*2+1] - jbody_coor[j_low*2-1];
    dx = jbody_coor[j_low*2  ] - jbody_coor[j_low*2-2];
  }
  // printf("j: dx = %f, dy=%f\n", dx, dy);
  if(dx < 0.001) theta_j = atan2(dy,dx);
  else theta_j = M_PI*0.5;
  // theta_j = atan2(dy, dx);

  // printf("compute_fvdw_fsh:1481 theta_j=%f\n", theta_j );
  // compute Fvdw
  double table_input[3], rotation;
  rotation = coordinate_transfer(table_input, xi, xj, theta_i, theta_j);
  // nn training theta in [0, pi/2)
  // map exp data to our force field dx = dx_true(4.4/6);
  // table_input[0] = table_input[0]*4.4/6;
  // table_input[1] = table_input[1]*4.4/6;

  // if(abs(table_input[0])< MINI_THETA) table_input[0] = 0;
  // printf("compute_fvdw_fsh: xi=[%f,%f];xj=[%f,%f];theta_i=%f;theta_j=%f;\ntable_input=[%f,%f,%f]\n",
  // xi[0],xi[1],xj[0],xj[1],theta_i,theta_j,table_input[0],table_input[1],table_input[2] );
  #ifdef _DEBUG_SHU
  printf("compute_fvdw_fsh: xi=[%f,%f];xj=[%f,%f];theta_i=%f;theta_j=%f;\ntable_input=[%f,%f,%f]\n",
  xi[0],xi[1],xj[0],xj[1],theta_i,theta_j,table_input[0],table_input[1],table_input[2] );
  #endif

  if(table_input[1] < 4.01){
    fvdw_j[0] = 0;
    fvdw_j[1] = 0;
    fvdw_j[2] = 0;

    fsh_j[0] = 0;
    fsh_j[1] = 0;
    fsh_j[2] = 0;
    // cal_fsh_nn(fsh_j, table_input);
  }
  else{
    cal_fvdw_nn(fvdw_j, table_input);
    // printf("after cal_fvdw_nn\n");
    if((table_input[1] > 5.0 )||((table_input[0]*table_input[0] + table_input[1] * table_input[1]) > 7.1*7.1 )) {
      fsh_j[0] = 0.0;
      fsh_j[1] = 0.0;
      fsh_j[2] = 0.0;
    }
    else{
      cal_fsh_nn(fsh_j, table_input);
    }
  }



  // printf("after cal_fsh_nn\n");
  // #ifdef _DUMP_F
  // if(update->ntimestep > dump_time && update->ntimestep % DT == 0){
  //   //open file and write fx, fy, fz;
  //   std::ofstream outfile;
  //   sprintf(file_output, "%s_Fvdw_Fsh_%d.txt",file_output_head,update->ntimestep);
  //   outfile.open(file_output, std::ios_base::app); // append instead of overwrite
  //   if ( !outfile.is_open() ) {
  //   std::cout << "Couldn't open file"<<file_output<<"\n";
  //   }
  //   outfile <<"0 "<<table_input[0]<< " "<<table_input[1]<< " "<<table_input[2]<< " "
  //   <<fvdw_j[0] << " " <<fvdw_j[1] << " " <<fvdw_j[2] << " "
  //   <<fsh_j[0] << " " <<fsh_j[1] << " " <<fsh_j[2]  <<"\n";
  //   outfile.close();
  // }
  // #endif

  #ifdef _DEBUG_SHU
  printf("-----------: fvdw_j=[%f,%f,%f],fsh_j=[%f,%f,%f]\n",fvdw_j[0],fvdw_j[1],fvdw_j[2],fsh_j[0],fsh_j[1],fsh_j[2]);
  #endif
  // printf("-----------: fvdw_j=[%f,%f,%f],fsh_j=[%f,%f,%f]\n",fvdw_j[0],fvdw_j[1],fvdw_j[2],fsh_j[0],fsh_j[1],fsh_j[2]);
  double fvdw_j_new[2],fsh_j_new[2];
  global_to_local(fvdw_j_new, fvdw_j, -1*rotation);
  global_to_local(fsh_j_new, fsh_j, -1*rotation);
  fvdw_j[0] = fvdw_j_new[0];
  fvdw_j[1] = fvdw_j_new[1];
  fsh_j[0] = fsh_j_new[0];
  fsh_j[1] = fsh_j_new[1];
  // printf("------Global: fvdw_j=[%f,%f,%f];rotation=%f\n",fvdw_j[0],fvdw_j[1],fvdw_j[2],rotation);
  #ifdef _DEBUG_SHU
  printf("------Global: fvdw_j=[%f,%f,%f];rotation=%f\n",fvdw_j[0],fvdw_j[1],fvdw_j[2],rotation);
  #endif

  // #ifdef _DUMP_F
  // if(update->ntimestep > dump_time && update->ntimestep % DT == 0){
  //   //open file and write fx, fy, fz;
  //   std::ofstream outfile;
  //   sprintf(file_output, "%s_Fvdw_Fsh_%d.txt",file_output_head,update->ntimestep);
  //   outfile.open(file_output, std::ios_base::app); // append instead of overwrite
  //   if ( !outfile.is_open() ) {
  //   std::cout << "Couldn't open file"<<file_output<<"\n";
  //   }
  //   outfile <<"1 "<<table_input[0]<< " "<<table_input[1]<< " "<<table_input[2]<< " "
  //   <<fvdw_j[0] << " " <<fvdw_j[1] << " " <<fvdw_j[2] << " "
  //   <<fsh_j[0] << " " <<fsh_j[1] << " " <<fsh_j[2]  <<"\n";
  //   outfile.close();
  // }
  // #endif
  // printf("compute_fvdw_fsh:1491 fvdw_j=[%f,%f,%f]\n", fvdw_j[0],fvdw_j[1],fvdw_j[2] );

  // if(update->ntimestep == 100){
    // printf("------time=%d,input=[%f,%f,%f],output=[%f,%f,%f]\n",
    // update->ntimestep,table_input[0],table_input[1],table_input[2],fvdw_j[0],fvdw_j[1],fvdw_j[2] );
  // }

}

/* ---------------------------------------------------------------------- */
void PairBodyRoundedPolygonFvdwFshNN::global_to_local(double *xnew, double *x, double theta)
{
  xnew[0] = x[0]*cos(theta) + x[1]*sin(theta);
  xnew[1] =-x[0]*sin(theta) + x[1]*cos(theta);
}

/* ---------------------------------------------------------------------- */
double PairBodyRoundedPolygonFvdwFshNN::coordinate_transfer(double *table_input,
  double *X1, double *X2, double theta1,double theta2)
{
  double x2[2][2], x_new1[2], x_new2[2], xc[2], theta, rotation;
  int k = 0;

  x2[0][0] = X2[0] - X1[0];
  x2[0][1] = X2[1] - X1[1];
  x2[1][0] = X2[0] + 1*cos(theta2) - X1[0];
  x2[1][1] = X2[1] + 1*sin(theta2) - X1[1];
  // printf("coordinate_transfer:1518, X1[0]=%f,X1[1]=%f, X2[0]=%f,X2[1]=%f\n", X1[0],X1[1], X2[0],X2[1] );
  x_new1[0] = 0; x_new1[1] = 0;
  x_new2[0] = 0; x_new2[1] = 0;
  while ( (abs(x_new1[1])< abs(x_new1[0])) || (x_new1[1]<=0.0001)){
    rotation = theta1+k*M_PI/2;
    global_to_local(x_new1, x2[0], rotation);
    global_to_local(x_new2, x2[1], rotation);
    xc[0] = x_new2[0] - x_new1[0];
    xc[1] = x_new2[1] - x_new1[1];

    if(abs(xc[0]) < MINI_THETA) theta = 0.5*M_PI;
    else theta = atan(xc[1]/xc[0]);

    if(abs(theta) < MINI_THETA) theta = 0;

    if(theta<0)
        theta += M_PI/2;

    k = k+1;
    // printf("coordinate_transfer:1518, X1[0]=%f,X1[1]=%f, X2[0]=%f,X2[1]=%f\n", X1[0],X1[1], X2[0],X2[1] );
    // printf("coordinate_transfer:1535 k=%d,x_new1[0]=%f,x_new1[1]=%f \n",k,x_new1[0],x_new1[1] );
  }

  table_input[0] = x_new1[0];
  table_input[1] = x_new1[1];
  table_input[2] = theta;
  // printf("coordinate_transfer:1549, X1[0]=%f,X1[1]=%f, X2[0]=%f,X2[1]=%f,theta1=%f,theta2=%f,table_input=[%f,%f,%f]\n", X1[0],X1[1], X2[0],X2[1],theta1,theta2,table_input[0],table_input[1],table_input[2]);
  return rotation;

}
/* ---------------------------------------------------------------------- */
void PairBodyRoundedPolygonFvdwFshNN::read_param(MatrixXd& param_m, VectorXd& param_v, char* paramFileName, int nrow, int ncol){
  int i,j;
  std::ifstream MyReadFile;
  MyReadFile.open(paramFileName);
  if ( MyReadFile.is_open() ) {
    std::cout << "File is opened "<< paramFileName<<"\n";
  }
  else {
    std::cout << "Couldn't open file"<<paramFileName<<"\n";
    exit(0);
  }
  // double tmp;
  // std::cout << "The matrix param_m is of size "
  //           << param_m.rows() << "x" << param_m.cols() << std::endl;

  for(i=0; i<nrow-1; i++){
    for(j=0; j<ncol; j++){
      MyReadFile >> param_m(i,j);
    }
  }
  i = nrow - 1;
  for(j=0; j<ncol; j++){
    MyReadFile >> param_v(j);
    // printf("i=%d,j=%d,value=%f\n",i,j,param_v(j) );
  }


  MyReadFile.close();
}

/* ---------------------------------------------------------------------- */
void PairBodyRoundedPolygonFvdwFshNN::cal_fvdw_nn(double * output, double * input){
  VectorXd FxyT0, FxyT2, FxyT4, FxyT6, FxyT0_tanh, FxyT2_tanh, FxyT4_tanh, FxyT6_tanh, Input_v;
  Input_v.resize(3);
  FxyT0.resize(param_fsh_v_0.size());
  FxyT2.resize(param_fsh_v_2.size());
  FxyT4.resize(param_fsh_v_4.size());
  FxyT6.resize(param_fsh_v_6.size());
  FxyT0_tanh.resize(param_fsh_v_0.size());
  FxyT2_tanh.resize(param_fsh_v_2.size());
  FxyT4_tanh.resize(param_fsh_v_4.size());
  FxyT6_tanh.resize(param_fsh_v_6.size());
  #ifdef _DEBUG_SHU
  printf("time=%d, Fvdw input=[%f,%f,%f]\n",update->ntimestep,input[0],input[1],input[2]);
  #endif

  int flag_sign = 1;
  if(input[0] < 0){
    flag_sign = -1;
    Input_v(0) = abs(input[0]);
    Input_v(1) = input[1];
    Input_v(2) = M_PI*0.5 - input[2];
    if( abs(Input_v(2)-M_PI*0.5) < MINI_THETA) Input_v(2) = 0;
  }else{
    Input_v(0) = input[0];
    Input_v(1) = input[1];
    Input_v(2) = input[2];
    if( abs(Input_v(2)-M_PI*0.5) < MINI_THETA) Input_v(2) = 0;
  }

  FxyT0 = param_fvdw_m_0.transpose()*Input_v+ param_fvdw_v_0;
  FxyT0_tanh = FxyT0.array().tanh();
  FxyT2 = param_fvdw_m_2.transpose()*FxyT0_tanh + param_fvdw_v_2;
  FxyT2_tanh = FxyT2.array().tanh();
  FxyT4 = param_fvdw_m_4.transpose()*FxyT2_tanh + param_fvdw_v_4;
  FxyT4_tanh = FxyT4.array().tanh();
  FxyT6 = param_fvdw_m_6.transpose()*FxyT4_tanh + param_fvdw_v_6;

  // printf("input_v = [%f, %f, %f]\n", Input_v(0), Input_v(1), Input_v(2) );

  if(Input_v(0) > MINI_THETA_x1 && Input_v(0)< MINI_THETA_x2 && Input_v(2)< MINI_THETA_t  && (Input_v(1)< MINI_THETA_y) )  {
    // printf("lalala1\n");
    // FxyT6(0) = Input_v(0);
    FxyT6(2) = 0.01*Input_v(2);
  }

  if(Input_v(0) > MINI_THETA_x3 && Input_v(0)< MINI_THETA_x4 && Input_v(2)< MINI_THETA_t  && (Input_v(1)< MINI_THETA_y) )  {
    // printf("lalala1\n");
    // FxyT6(0) = Input_v(0);
    FxyT6(2) = 0.01*Input_v(2);
  }

  if(Input_v(0) > MINI_THETA_x1 && Input_v(0)< MINI_THETA_x2 && (M_PI*0.5 - Input_v(2) )< MINI_THETA_t && (Input_v(1)< MINI_THETA_y) ){
    // printf("lalala2\n");
    // FxyT6(0) = Input_v(0);
    FxyT6(2) = -0.01*(M_PI*0.5 - Input_v(2) );
  }

  if(Input_v(0) > MINI_THETA_x3 && Input_v(0)< MINI_THETA_x4 && (M_PI*0.5 - Input_v(2) )< MINI_THETA_t && (Input_v(1)< MINI_THETA_y) ){
    // printf("lalala2\n");
    // FxyT6(0) = Input_v(0);
    FxyT6(2) = -0.01*(M_PI*0.5 - Input_v(2) );
  }

  // When |Dx-Dy|<epsilon and |theta|<epsilon, we have Fx = Fy = (Fx+Fy)/2, T = 0
  if( abs(Input_v(0) - Input_v(1)) < MINI_THETA_x_diag){
    if(abs(Input_v(2)) < MINI_THETA_t_diag){
    FxyT6(0) = 0.5*( abs(FxyT6(0)) + abs(FxyT6(1)) );
    FxyT6(1) = 0.5*( abs(FxyT6(0)) + abs(FxyT6(1)) );
    FxyT6(2) = Input_v(2);
    }
    if((M_PI*0.5 - Input_v(2)) < MINI_THETA_t_diag){
    FxyT6(0) = 0.5*( abs(FxyT6(0)) + abs(FxyT6(1)) );
    FxyT6(1) = 0.5*( abs(FxyT6(0)) + abs(FxyT6(1)) );
    FxyT6(2) = -(M_PI*0.5 - Input_v(2) );
    }
  }

  if (( abs(Input_v(0)) < MINI_THETA_x_diag )&(abs(Input_v(2)) < MINI_THETA_t_diag)) {
    FxyT6(0) = 0.0;
    FxyT6(2) = 0.0;
  }


  output[0]= FxyT6(0)*flag_sign;
  output[1]= abs(FxyT6(1)); // make sure Fy is negtive in longrange
  output[2]= FxyT6(2)*flag_sign;

  // #ifdef _DEBUG_SHU
  // printf("Fvdw input=[%f,%f,%f]\n",Input_v(0),Input_v(1),Input_v(2) );
  // printf("Fvdw output=[%f,%f,%f]\n",output[0],output[1],output[2] );
  // #endif

  // #endif
}

/* ---------------------------------------------------------------------- */
void PairBodyRoundedPolygonFvdwFshNN::cal_fsh_nn(double * output, double * input){
  VectorXd FxyT0, FxyT2, FxyT4, FxyT6, FxyT0_tanh, FxyT2_tanh, FxyT4_tanh, FxyT6_tanh, Input_v;
  Input_v.resize(3);
  FxyT0.resize(param_fsh_v_0.size());
  FxyT2.resize(param_fsh_v_2.size());
  FxyT4.resize(param_fsh_v_4.size());
  FxyT6.resize(param_fsh_v_6.size());
  FxyT0_tanh.resize(param_fsh_v_0.size());
  FxyT2_tanh.resize(param_fsh_v_2.size());
  FxyT4_tanh.resize(param_fsh_v_4.size());
  FxyT6_tanh.resize(param_fsh_v_6.size());
  #ifdef _DEBUG_SHU
  printf("time=%d, Fsh input=[%f,%f,%f]\n",update->ntimestep,input[0],input[1],input[2]);
  #endif

  int flag_sign = 1;
  if(input[0] < 0){
    flag_sign = -1;
    Input_v(0) = -input[0];
    Input_v(1) = input[1];
    Input_v(2) = M_PI*0.5 - input[2];
    if( abs(Input_v(2)-M_PI*0.5) < MINI_THETA) Input_v(2) = 0;
  }
  else{
    Input_v(0) = input[0];
    Input_v(1) = input[1];
    Input_v(2) = input[2];
    if( abs(Input_v(2)-M_PI*0.5) < MINI_THETA) Input_v(2) = 0;
  }

  FxyT0 = param_fsh_m_0.transpose()*Input_v + param_fsh_v_0;
  FxyT0_tanh = FxyT0.array().tanh();

  FxyT2 = param_fsh_m_2.transpose()*FxyT0_tanh   + param_fsh_v_2;
  FxyT2_tanh = FxyT2.array().tanh();

  FxyT4 = param_fsh_m_4.transpose()*FxyT2_tanh   + param_fsh_v_4;
  FxyT4_tanh = FxyT4.array().tanh();

  FxyT6 = param_fsh_m_6.transpose()*FxyT4_tanh   + param_fsh_v_6;
  // FxyT6_tanh = FxyT6.array().tanh();



  // When |Dx-Dy|<epsilon and |theta|<epsilon, we have Fx = Fy = (Fx+Fy)/2, T = 0
  if( abs(Input_v(0) - Input_v(1)) < MINI_THETA_x_diag){
    if(abs(Input_v(2)) < MINI_THETA_t_diag){
    FxyT6(0) = - 0.5*( abs(FxyT6(0)) + abs(FxyT6(1)) );
    FxyT6(1) = - 0.5*( abs(FxyT6(0)) + abs(FxyT6(1)) );
    FxyT6(2) = Input_v(2);
    }
    if((M_PI*0.5 - Input_v(2)) < MINI_THETA_t_diag){
    FxyT6(0) = - 0.5*( abs(FxyT6(0)) + abs(FxyT6(1)) );
    FxyT6(1) = - 0.5*( abs(FxyT6(0)) + abs(FxyT6(1)) );
    FxyT6(2) = -(M_PI*0.5 - Input_v(2) );
    }
  }


  // // printf("NN input=[%f,%f,%f]\n",FxyT6(0), FxyT6(1), FxyT6(2) );
  // if((abs(Input_v(0) - 0.5) < 0.5) && ( abs(Input_v(1) - 5.6) < 0.05 ) && ( abs(Input_v(2) - 1.35) < 0.15 ) ){
  //   // if(FxyT6(2) > 0.00001){
  //   //   FxyT6(2) =  - 0.1;
  //   // }
  //   // if(FxyT6(2) < -0.00001){
  //     FxyT6(2) = - 0.1;
  // }
    // FxyT6(0) = FxyT6(0) + 0.1* Input_v(0);
    // FxyT6(1) = FxyT6(1) + 0.01* Input_v(1);


  if(Input_v(0) > MINI_THETA_x1 && Input_v(0)< MINI_THETA_x2 && Input_v(2)< MINI_THETA_t  && (Input_v(1)< MINI_THETA_y) ){
    // FxyT6(0) = -Input_v(0);
    FxyT6(2) = 0.2*Input_v(2);
  }

  if(Input_v(0) > MINI_THETA_x3 && Input_v(0)< MINI_THETA_x4 && Input_v(2)< MINI_THETA_t  && (Input_v(1)< MINI_THETA_y) ){
    // FxyT6(0) = -Input_v(0);
    FxyT6(2) = 0.2*Input_v(2);
  }

  if(Input_v(0) > MINI_THETA_x1 && Input_v(0)< MINI_THETA_x2 && (M_PI*0.5 - Input_v(2) )< MINI_THETA_t  && (Input_v(1)< MINI_THETA_y) ){
    // FxyT6(0) = Input_v(0);
    FxyT6(2) = -0.2*(M_PI*0.5 - Input_v(2) );
  }

  if(Input_v(0) > MINI_THETA_x3 && Input_v(0)< MINI_THETA_x4 && (M_PI*0.5 - Input_v(2) )< MINI_THETA_t  && (Input_v(1)< MINI_THETA_y) ){
    // FxyT6(0) = Input_v(0);
    FxyT6(2) = -0.2*(M_PI*0.5 - Input_v(2) );
  }

  output[0]= FxyT6(0)*flag_sign;
  output[1]= -abs(FxyT6(1));
  output[2]= FxyT6(2)*flag_sign;
  // printf("flag sign = %d\n",flag_sign );

  if(sqrt(output[0]*output[0] + output[1]*output[1] ) < 0.0001){
    output[0]= 0;
    output[1]= 0;
    output[2]= 0;
  }
  // #ifdef _DEBUG_SHU
  // printf("Fsh input=[%f,%f,%f]\n",Input_v(0),Input_v(1),Input_v(2) );
  // printf("Fsh output=[%f,%f,%f]\n",output[0],output[1],output[2] );
  // #endif

  #ifdef _DUMP_F_more
  if(update->ntimestep > dump_time && update->ntimestep % DT == 0){
    //open file and write fx, fy, fz;
    std::ofstream outfile;
    sprintf(file_output, "%s_Fvdw_Fsh_%d.txt",file_output_head,update->ntimestep);
    outfile.open(file_output, std::ios_base::app); // append instead of overwrite
    if ( !outfile.is_open() ) {
    std::cout << "Couldn't open file"<<file_output<<"\n";
    }
    outfile <<"Input: "<<Input_v(0)<<" "<<Input_v(1)<<" "<<Input_v(2)<< "\n";

    outfile <<"param_fsh_m_0="<<param_fsh_m_0(0,0)<< " "<<param_fsh_m_0(0,1)<< " "<<param_fsh_m_0(1,1)<<" "<<param_fsh_v_0(0)<<"\n"
    <<"param_fsh_m_2="<<param_fsh_m_2(0,0)<< " "<<param_fsh_m_2(0,1)<< " "<<param_fsh_m_2(1,1)<<" "<<param_fsh_v_2(0)<<"\n"
    <<"param_fsh_m_4="<<param_fsh_m_4(0,0)<< " "<<param_fsh_m_4(0,1)<< " "<<param_fsh_m_4(1,1)<<" "<<param_fsh_v_4(0)<<"\n"
    <<"param_fsh_m_6="<<param_fsh_m_6(0,0)<< " "<<param_fsh_m_6(0,1)<< " "<<param_fsh_m_6(1,1)<<" "<<param_fsh_v_6(0) <<"\n";

    outfile <<"fsh FxyT0="<<FxyT0(0)<< " "<<FxyT0(1)<< " "<<FxyT0(2)<<"\n"
    <<" FxyT0_tanh="<<FxyT0_tanh(0)<< " "<<FxyT0_tanh(1)<< " "<<FxyT0_tanh(2)<<"\n"
    <<" FxyT2="<<FxyT2(0)<< " "<<FxyT2(1)<< " "<<FxyT2(2)<<"\n"
    <<" FxyT2_tanh="<<FxyT2_tanh(0)<< " "<<FxyT2_tanh(1)<< " "<<FxyT2_tanh(2)<<"\n"
    <<" FxyT4="<<FxyT4(0)<< " "<<FxyT4(1)<< " "<<FxyT4(2)<<"\n"
    <<" FxyT4_tanh="<<FxyT4_tanh(0)<< " "<<FxyT4_tanh(1)<< " "<<FxyT4_tanh(2)<<"\n"
    <<" FxyT6="<<FxyT6(0)<< " "<<FxyT6(1)<< " "<<FxyT6(2)<<"\n"
    <<" FxyT6_tanh="<<FxyT6_tanh(0)<< " "<<FxyT6_tanh(1)<< " "<<FxyT6_tanh(2)<<"\n";
    outfile.close();
  }
  #endif
}
