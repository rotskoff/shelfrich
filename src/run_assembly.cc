/* run_assembly.cc Runs assembly simulations.
 * Copyright (C) 2018 Grant M. Rotskoff
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 **/

#include <stdio.h>
#include <math.h>
#include <gsl/gsl_rng.h>
#include "cargo.h"
#include "graph.h"
#include "montecarlo.h"

gsl_rng *r;
int main(int argc, char **argv) {

  // parse inputs
  if (argc!=5) {
    fprintf(stderr,"usage: ./assemble seed lnK lnz Rc\n");
    exit(-1);
  }

  // initialize the rng
  const gsl_rng_type *t;
  t = gsl_rng_taus2;
  r = gsl_rng_alloc(t);
  long unsigned int seed = atoi(argv[1]);
  srand((unsigned)seed);
  gsl_rng_set(r, seed);

  const double pi = 4*atan(1);


  // set up the cargo
  double eps = 1.0;
  double mu = -5.0*eps; // 6 is coexistence
  double gamma = 5.;
  double gammao = 100.;
  double T = 1.;
  int Nmaxc = 200;
  cargo cargo(Nmaxc, eps, mu, gamma, gammao, T);
  cargo.empty();

  int Rc = atoi(argv[4]);

  // create a cluster
  for (int i=-2*Rc; i<2*Rc; i++) {
    for (int j=-2*Rc; j<2*Rc; j++) {
      for (int k=-2*Rc; k<2*Rc; k++) {
        double x = i/2.0; double y = j/2.0; double z = k/2.0;
        if ((x*x+y*y+z*z) < Rc*Rc) {
          int is = i + Nmaxc/2;
          int js = j + Nmaxc/2;
          int ks = k + Nmaxc/2;
          if (((is+js+ks)%2)==0) {
            int site = is + js*Nmaxc + ks*Nmaxc*Nmaxc;
            if (cargo.lat[site]!=1) {
              cargo.lat[site] = 1;
              cargo.sites.insert(site);
              cargo.N += 1;
            }
          }
        }
      }
    }
  }

  cargo.update_surf();
  // initialize the graph
  graph g;
  g.epsilon = 500.;
  g.l0 = 1;
  g.lmax = 1.5;
  g.lmin = 0.5;
  g.lfuse = 0.5;
  // defines the insertion radius for monomer insertion
  g.xi = g.lfuse;
  int Nmax = 10000;
  int Ntype = 1;
  g.initialize(Nmax, 6*Nmax, Ntype);
  g.kappa[0] = 12.5;
  g.min_bond[0] = 0.;
  g.T = 1.;
  double lnK = atof(argv[2]);
  double lnz = atof(argv[3]);
  g.K = exp(lnK);
  g.z = exp(lnz);
  g.vfuse = 4./3*pi*pow(g.lfuse,3);

  fprintf(stdout, "Graph initialized.\n");
  dump_log(g);

  // build the initial configuration
  make_hexamer(g);
  g.center_at(0, 0., Rc);

  g.update_t();
  g.orient();
  g.update_normals();
  cargo.compute_occupied(g);

  double etot = g.compute_energy();
  fprintf(stdout, "Graph has %d vertices and %d edges\n",g.Nv,g.Ne);
  fprintf(stdout, "Initial energy is %f\n", etot);
  // do some mc

  int ind, e;
  // set up an output file
  FILE *efile, *dfile;
  efile = fopen("energy.dat","w");
  dfile = fopen("defects.dat","w");

  double kc0 = 0.000075;
  double ks0 = 0.0075;

  int sweep = 0;
  while (g.Nperim>0)
  {
    // equilibrate the structure with vertex moves
    for (int step=0; step<g.Nv; step++)
    {
      vertex_move_cargo(g,cargo,r);
    }


    for (int ei=0; ei<g.Ne; ei++)
    {
      int valid = g.update_edges(ei);
      if (valid!=1) {
        fprintf(stderr, "Invalid configuration at sweep %d\n",sweep);
        exit(EXIT_FAILURE);
      }
    }

    // attempt fusion / fission
    // technically, atttempted with kfuse0*Nperim, but kfuse0=1
    // so if there's at least one move, we attempt it
    attempt_fusion(g, cargo, r);
    attempt_fission(g, cargo, r);


    // attempt a cargo move
    double pc_attempt = kc0*cargo.Nsurf;
    if (gsl_rng_uniform(r)<pc_attempt)
    {
      cargo.update_surf();
      cargo.mc_step_controller(r);
    }

    // attempt to delete a monomer
    double ps_attempt = ks0*g.Nperim;
    if (gsl_rng_uniform(r)<ps_attempt)
    {
      if (sweep>0 && g.Nperim>0) {
        ind = gsl_rng_uniform_int(r, g.Nperim);
        e = g.esurf[ind];
        delete_monomer(g, cargo, e, r);
      }
    }

    // attempt to insert a monomer
    ps_attempt = ks0*g.Nperim;
    if (gsl_rng_uniform(r)<ps_attempt && g.Nperim>0) {
      ind = gsl_rng_uniform_int(r, g.Nperim);
      e = g.esurf[ind];
      // delete monomers that violate overlap considerations
      if (g.check_overlap(e)==1) {
        delete_monomer(g, cargo, e, r);
      }
      else {
        if (g.check_wedge(e)==1) {
          insert_monomer_wedge(g, cargo, e, r);
        } else {
          insert_monomer_free(g, cargo, e, r);
        }
      }
    }

    // compute the energy
    if (sweep % 100 == 0)
    {
      fprintf(efile,"%d %d %d %f\n", sweep, cargo.N, g.Nm, g.compute_energy());
      fflush(efile);
      dump_defects(g, dfile, sweep);
    }

    // dump output files
    if (sweep%1000 == 0)
    {
      fprintf(stdout,"%d %d %d %f\n", sweep, cargo.N, g.Nm, g.compute_energy());
      cargo.dump_xml(g,sweep);
    }

    for (int vi=0; vi<g.Nv; vi++) {
      if (g.surface_vertex(vi)==0 && g.nneigh[vi]<5) {
        fprintf(stderr, "Catastrophic error: %d\n", vi);
        dump_norm_xml(g,sweep);
        exit(-1);
      }
    }

    sweep++;
  }
  int relax_sweep = 0;
  while (relax_sweep<100000)
  {

    // equilibrate the structure with vertex moves
    for (int step=0; step<g.Nv; step++)
    {
      vertex_move_cargo(g,cargo,r);
    }
    sweep++;
    relax_sweep++;
  }
  // close everything out
  cargo.dump_xml(g,sweep);
  fprintf(efile,"%d %d %d %f\n", sweep, cargo.N, g.Nm, g.compute_energy());
  dump_defects(g, dfile, sweep);
  fclose(efile);
  fclose(dfile);
  return 0;
}
