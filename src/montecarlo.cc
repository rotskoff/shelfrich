/* montecarlo.cc Implements the Metropolis Monte Carlo moves for assembly.
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


#include "montecarlo.h"
#include "math.h"
#include "graph.h"
#include <algorithm>

void vertex_move(graph &g, bool override, gsl_rng *r)
{
  double d = (0.5-gsl_rng_uniform(r))*2.0/sqrt(g.epsilon);
  int vi = gsl_rng_uniform_int(r,g.Nv);
  int dim = gsl_rng_uniform_int(r,3);
  double acc = gsl_rng_uniform(r);
  double etrial = 0;
  // make a perturbation in a random direction to v
  double eref = g.vertex_energy(vi);
  g.v[vi][dim] += d;
  // update the edge lengths
  int valid = g.update_edges(vi);
  if (valid == 1 || override == 1) {
    etrial = g.vertex_energy(vi);
    double de = etrial - eref;
    if (de>0) {
      if (acc < exp(-de/g.T)) {
        g.energy += de;
      } else {
        g.v[vi][dim] -= d;
        g.update_edges(vi);
      }
    } else {
      g.energy += de;
    }
  }
  else {
    g.v[vi][dim] -= d;
    g.update_edges(vi);
  }
  return;
}

void vertex_move_cargo(graph &g, cargo &c, gsl_rng *r)
{

  double d = (0.5-gsl_rng_uniform(r))*2.0/sqrt(g.epsilon);
  int vi = gsl_rng_uniform_int(r,g.Nv);
  int dim = gsl_rng_uniform_int(r,3);
  double acc = gsl_rng_uniform(r);
  double etrial = 0, de = 0;
  // make a perturbation in a random direction to v
  double eref = g.vertex_energy(vi) + c.cargo_graph_energy_vertex(vi,g);
  g.v[vi][dim] += d;
  // update the edge lengths
  int valid = g.update_edges(vi);
  if (valid == 1) {
    etrial = g.vertex_energy(vi) + c.cargo_graph_energy_vertex(vi,g);
    de = etrial - eref;
    if (de>0) {
      if (acc < exp(-de/g.T)) {
        g.energy += de;
      } else {
        g.v[vi][dim] -= d;
        g.update_edges(vi);
      }
    } else {
      g.energy += de;
    }
  }
  else {
    g.v[vi][dim] -= d;
    g.update_edges(vi);
  }
  return;
}

int attempt_fusion(graph &g, cargo &c, gsl_rng *r)
{
  // choose the vertices
  int vi, vj;
  g.update_fusion_set();
  if (g.Nfusion==0) {
    return 0;
  }
  int iset = gsl_rng_uniform_int(r, g.Nfusion);
  set< pair <int,int> >::iterator it;
  it = g.vfusion.begin();
  advance(it,iset);
  pair<int, int> pij = *it;
  vi = pij.first;
  vj = pij.second;
  // choose a random fusable vertex
  if (vi==vj) {
    return 0;
  }

  // store the initial state
  vector<double> vi_init, vj_init;
  vector<int> vi_neigh_init, vj_neigh_init;
  int vi_nneigh_init, vj_nneigh_init;
  g.save_local_state(vi, vi_init, vi_neigh_init, vi_nneigh_init);
  g.save_local_state(vj, vj_init, vj_neigh_init, vj_nneigh_init);
  g.adjust_indices(vi, vj, vi_nneigh_init, vi_neigh_init);
  g.adjust_indices(vj, vi, vj_nneigh_init, vj_neigh_init);

  // check the total number of resulting neighbors
  int vk;
  int nshared = 0, nshared_edge = 0;
  for (int nj=0; nj<g.nneigh[vj]; nj++) {
    vk = g.neigh[vj][nj];
    nshared_edge += g.check_next_neigh(vi, vk);
    for (int ni=0; ni<g.nneigh[vi]; ni++) {
      if (g.neigh[vi][ni]==vk) {
        nshared += 1;
        if (g.surface_edge(g.eind[vi][vk])==0 || g.surface_edge(g.eind[vj][vk])==0) {
          return 0;
        }
      }
    }
  }

  // check if the edge fusion criterion is met
  for (int ni=0; ni<g.nneigh[vi]; ni++) {
    int ei = g.eind[vi][g.neigh[vi][ni]];
    if (g.surface_edge(ei)==1) {
      for (int nj=0; nj<g.nneigh[vj]; nj++) {
        int ej = g.eind[vj][g.neigh[vj][nj]];
        if (g.surface_edge(ej)==1) {
          if (g.check_edge_dist(ei, ej)==1 && nshared==0) {
            return attempt_edge_fusion(g, c, vi, vj, r);
          }
        }
      }
    }
  }

  if (nshared==0 && nshared_edge==0) {
    return 0;
  }


  g.add_midpoint_vertex(vi,vj);

  // check distance constraints
  if (g.check_distance_constraints(vi)==0 || g.check_distance_constraints(vj)==0)
  {
    g.delete_vertex(g.Nv-1);
    return 0;
  }

  double e_init = g.compute_energy();

  int n_added = 0;
  g.add_new_edges(vi, g.Nv-1, n_added);
  g.add_new_edges(vj, g.Nv-1, n_added);

  int valid = g.update_edges(g.Nv-1);
  if (valid == -1 || g.check_defect(g.Nv-1)==0) {
    g.remove_added(n_added);
    return 0;
  }

  g.remove_edges(vi);
  g.remove_edges(vj);
  g.delete_vertex(vi);
  if (vi<vj) {
    g.delete_vertex(vj-1);
  } else g.delete_vertex(vj);

  // check overlap
  for (int ne=0; ne<n_added; ne++) {
    // prevent sevenfold defects
    int nvi, nvj, svi, svj;
    nvi = g.nneigh[g.e[g.Ne-1-ne][0]];
    svi = g.surface_vertex(g.e[g.Ne-1-ne][0]);
    nvj = g.nneigh[g.e[g.Ne-1-ne][1]];
    svj = g.surface_vertex(g.e[g.Ne-1-ne][1]);
    if (g.check_overlap(g.Ne-1-ne)==1 || (nvi>6 && svi==0) || (nvj>6 && svj==0)) {
      // delete the new vertex / edges
      g.remove_added(n_added);
      // restore the saved state
      g.restore_local_state(vi_init, vi_neigh_init, vi_nneigh_init);
      g.restore_local_state(vj_init, vj_neigh_init, vj_nneigh_init);
      g.update_topology();
      return 0;
    }
  }

  // compute the energy difference
  g.update_topology();
  double e_final = g.compute_energy();
  double de = e_final - e_init;

  if (gsl_rng_uniform(r) < g.K/g.vfuse*exp(-de/g.T)) {
    fprintf(stderr, "Fusion %d-%d\n", vi, vj);
    return 1;
  }
  else {
    // delete the added edges
    g.remove_added(n_added);
    g.restore_local_state(vi_init, vi_neigh_init, vi_nneigh_init);
    g.restore_local_state(vj_init, vj_neigh_init, vj_nneigh_init);
    g.update_topology();
    return 0;
  }
}

int attempt_edge_fusion(graph &g, cargo &c, int vi, int vj, gsl_rng *r)
{
  if (vi==vj) {
    return 1;
  }
  int vi1=-1, vj1=-1; // define the other fusion pair
  for (int ni=0; ni<g.nneigh[vi]; ni++) {
    int ei = g.eind[vi][g.neigh[vi][ni]];
    if (g.surface_edge(ei)==1) {
      for (int nj=0; nj<g.nneigh[vj]; nj++) {
        int ej = g.eind[vj][g.neigh[vj][nj]];
        if (g.surface_edge(ej)==1) {
          if (g.check_edge_dist(ei, ej)==1) {
            vi1 = g.e[ei][0];
            if (vi1==vi) { vi1 = g.e[ei][1]; }
            vj1 = g.e[ej][0];
            if (vj1==vj) { vj1 = g.e[ej][1]; }
          }
        }
      }
    }
  }

  if (vi1==-1 || vj1==-1) {
    fprintf(stderr, "Edge fusion not possible on surface edge\n");
    exit(EXIT_FAILURE);
  }

  double e_init = g.compute_energy();

  g.delete_edge(g.eind[vi][vi1]);
  g.delete_edge(g.eind[vj][vj1]);

  // store the initial state
  vector<double> vi_init, vj_init, vi1_init, vj1_init;
  vector<int> vi_neigh_init, vj_neigh_init, vi1_neigh_init, vj1_neigh_init;
  int vi_nneigh_init, vj_nneigh_init, vi1_nneigh_init, vj1_nneigh_init;
  g.save_local_state(vi, vi_init, vi_neigh_init, vi_nneigh_init);
  g.save_local_state(vj, vj_init, vj_neigh_init, vj_nneigh_init);
  g.save_local_state(vi1, vi1_init, vi1_neigh_init, vi1_nneigh_init);
  g.save_local_state(vj1, vj1_init, vj1_neigh_init, vj1_nneigh_init);
  g.adjust_indices_edge_fusion(vi, vi1, vj, vj1, vi_neigh_init);
  g.adjust_indices_edge_fusion(vj, vj1, vi, vi1, vj_neigh_init);
  g.adjust_indices_edge_fusion(vi1, vi, vj, vj1, vi1_neigh_init);
  g.adjust_indices_edge_fusion(vj1, vj, vi, vi1, vj1_neigh_init);


  g.add_midpoint_vertex(vi, vj);
  if (g.check_distance_constraints(vi)==0 ||
      g.check_distance_constraints(vj)==0 )
  {
    g.delete_vertex(g.Nv-1);
    g.add_edge(vi,vi1);
    g.add_edge(vj,vj1);
    return 0;
  }

  g.add_midpoint_vertex(vi1, vj1);
  if (g.check_distance_constraints(vi1)==0 ||
      g.check_distance_constraints(vj1)==0 )
  {
    g.delete_vertex(g.Nv-1);
    g.delete_vertex(g.Nv-1);
    g.add_edge(vi,vi1);
    g.add_edge(vj,vj1);
    return 0;
  }


  int n_added = 0, success=1;
  success*=g.add_new_edges(vi, g.Nv-2, n_added);
  success*=g.add_new_edges(vj, g.Nv-2, n_added);
  int valid = g.update_edges(g.Nv-2);
  if (valid == -1 || g.check_defect(g.Nv-2)==0 || success==0) {
    g.remove_added(n_added);
    g.delete_vertex(g.Nv-1);
    g.add_edge(vi,vi1);
    g.add_edge(vj,vj1);
    g.update_topology();
    return 0;
  }

  success*=g.add_new_edges(vi1, g.Nv-1, n_added);
  success*=g.add_new_edges(vj1, g.Nv-1, n_added);
  // add the mutual edge
  int last_add = g.add_edge(g.Nv-1, g.Nv-2); // no guarantee that these are in range
  success *= last_add;
  n_added+=last_add;
  valid = g.update_edges(g.Nv-1);
  if (valid == -1 || g.check_defect(g.Nv-1)==0 || success==0) {
    g.remove_added(n_added);
    g.delete_vertex(g.Nv-1);
    g.add_edge(vi,vi1);
    g.add_edge(vj,vj1);
    g.update_topology();
    return 0;
  }

  g.remove_edges(vi);
  g.remove_edges(vj);
  g.remove_edges(vi1);
  g.remove_edges(vj1);

  // sort to delete largest first
  int varr[] = {vi, vj, vi1, vj1};
  vector<int> vvec (varr, varr+sizeof(varr)/sizeof(int));
  sort(vvec.begin(), vvec.end());
  for (int i=3; i>-1; i--) {
    int vdel = vvec[i];
    g.delete_vertex(vdel);
  }

  // check overlap
  for (int ne=0; ne<n_added; ne++) {
    // prevent sevenfold defects
    int nvi, nvj, svi, svj;
    nvi = g.nneigh[g.e[g.Ne-1-ne][0]];
    svi = g.surface_vertex(g.e[g.Ne-1-ne][0]);
    nvj = g.nneigh[g.e[g.Ne-1-ne][1]];
    svj = g.surface_vertex(g.e[g.Ne-1-ne][1]);
    if (g.check_overlap(g.Ne-1-ne)==1 || (nvi>6 && svi==0) || (nvj>6 && svj==0)) {
      g.remove_added(n_added);
      g.delete_vertex(g.Nv-1);
      g.restore_local_state(vi_init, vi_neigh_init, vi_nneigh_init);
      g.restore_local_state(vj_init, vj_neigh_init, vj_nneigh_init);
      g.restore_local_state(vi1_init, vi1_neigh_init, vi1_nneigh_init);
      g.restore_local_state(vj1_init, vj1_neigh_init, vj1_nneigh_init);
      g.add_edge(g.Nv-2,g.Nv-4);
      g.add_edge(g.Nv-1,g.Nv-3);
      g.update_topology();
      return 0;
    }
  }
  // compute the energy difference
  g.update_topology();
  double e_final = g.compute_energy();
  double de = e_final - e_init;

  // factor of 0.5, can choose either pair vi,vj or vi1, vj1
  if (gsl_rng_uniform(r) < 0.5*(g.K*g.K)/(g.vfuse*g.vfuse)*exp(-de/g.T)) {
    return 1;
  }
  else {
    g.remove_added(n_added);
    g.delete_vertex(g.Nv-1);
    g.restore_local_state(vi_init, vi_neigh_init, vi_nneigh_init);
    g.restore_local_state(vj_init, vj_neigh_init, vj_nneigh_init);
    g.restore_local_state(vi1_init, vi1_neigh_init, vi1_nneigh_init);
    g.restore_local_state(vj1_init, vj1_neigh_init, vj1_nneigh_init);
    g.add_edge(g.Nv-2,g.Nv-4); // vi1, vi edge
    g.add_edge(g.Nv-1,g.Nv-3); // vj1, vj edge
    g.update_topology();
    return 0;
  }
}

int attempt_fission(graph &g, cargo &c, gsl_rng *r)
{
  int vi, vj;
  bool edge_fission=false;
  // choose a random fissionable vertex and partner
  g.update_fission_set();
  if (g.Nfission==0) {
    return 0;
  }
  int iset = gsl_rng_uniform_int(r, g.Nfission);
  set< pair <int, int> >::iterator it;
  it = g.vfission.begin();
  advance(it,iset);
  pair<int, int> pij = *it;
  vi = pij.first;
  vj = pij.second;
  //vj = g.select_random_fission_partner(vi, r);
  // no fission on the initial hexamer
  if (vi<7 || vj<7)
  {
    return 0;
  }

  // find the edge to attempt fission
  if (g.surface_vertex(vi)==1 && g.surface_vertex(vj)==1)
  {
    edge_fission = true;
  }

  if (edge_fission==true) {
    return attempt_edge_fission(g, c, vi, vj, r);
  }

  if (g.surface_vertex(vi)==0 && g.surface_vertex(vj)==1) {
    int vit = vi;
    vi = vj;
    vj = vit;
  } else if (g.surface_vertex(vi)==0 && g.surface_vertex(vj)==0) {
    return 0;
  }

  double e_init = g.compute_energy();

  vector<double> vi_init;
  vector<int> vi_neigh_init;
  int vi_nneigh_init;
  g.save_local_state(vi, vi_init, vi_neigh_init, vi_nneigh_init);
  g.adjust_indices(vi, g.Nvmax, vi_nneigh_init, vi_neigh_init); // second arg prevents double adjustment

  g.add_fission_vertices(vi, r);
  int split_edge = g.eind[vi][vj];
  int success = 1;
  double dx1, dx2;
  success *= g.add_edge(g.Nv-2, vj);
  success *= g.add_edge(g.Nv-1, vj);
  int vc1, vc2;
  vc1 = g.t[split_edge][0];
  vc2 = g.t[split_edge][1];
  // distribute the neighbors of vi
  dx1 = dist(g.v[g.Nv-2], g.v[vc1]);
  dx2 = dist(g.v[g.Nv-1], g.v[vc1]);
  if (dx1<dx2) {
    success *= g.connect_neighs(vc1, vi, vj, g.Nv-2);
    success *= g.connect_neighs(vc2, vi, vj, g.Nv-1);
  } else {
    success *= g.connect_neighs(vc2, vi, vj, g.Nv-2);
    success *= g.connect_neighs(vc1, vi, vj, g.Nv-1);
  }

  // check constraints
  int valid1 = g.update_edges(g.Nv-2);
  int valid2 = g.update_edges(g.Nv-1);
  if (valid1!=1 || valid2!=1 || success==0) {
    g.remove_newest_vertex(2);
    return 0;
  }

  // delete vi
  g.remove_edges(vi);
  g.delete_vertex(vi);

  // compute the energy difference
  g.update_topology();
  double e_final = g.compute_energy();
  double de = e_final - e_init;

  // accept or reject
  double acc = gsl_rng_uniform(r);
  if (acc < g.vfuse/g.K * exp(-de/g.T))
  {
    // delete vi
    fprintf(stderr, "Fission of %d\n", vi);
    return 1;
  }
  else {
    g.remove_newest_vertex(2);
    g.restore_local_state(vi_init, vi_neigh_init, vi_nneigh_init);
    g.update_topology();
    return 0;
  }
}

int attempt_edge_fission(graph &g, cargo &c, int vi, int vj, gsl_rng *r)
{
  vector<double> vi_init, vj_init;
  vector<int> vi_neigh_init, vj_neigh_init;
  int vi_nneigh_init, vj_nneigh_init;
  g.save_local_state_fission(vj, vi, vi_init, vi_neigh_init, vi_nneigh_init);
  g.save_local_state_fission(vi, vj, vj_init, vj_neigh_init, vj_nneigh_init);
  g.adjust_indices(vi, vj, vi_nneigh_init, vi_neigh_init);
  g.adjust_indices(vj, vi, vj_nneigh_init, vj_neigh_init);


  int success = 1;
  int split_edge = g.eind[vi][vj];

  if (g.surface_edge(split_edge)==1)
  {
    return 0;
  }

  double e_init = g.compute_energy();

  int vc1 = g.t[split_edge][0], vc2 = g.t[split_edge][1];

  g.add_fission_vertices(vi,r);
  success *= g.connect_neighs(vc1, vi, vj, g.Nv-2);
  success *= g.connect_neighs(vc2, vi, vj, g.Nv-1);

  g.add_fission_vertices(vj,r);
  success *= g.connect_neighs(vc1, vj, vi, g.Nv-2);
  success *= g.connect_neighs(vc2, vj, vi, g.Nv-1);

  success *= g.add_edge(g.Nv-1, g.Nv-3);
  success *= g.add_edge(g.Nv-2, g.Nv-4);

  // check constraints
  int valid1 = g.update_edges(g.Nv-1);
  int valid2 = g.update_edges(g.Nv-2);
  int valid3 = g.update_edges(g.Nv-3);
  int valid4 = g.update_edges(g.Nv-4);
  if (valid1!=1 || valid2!=1 || valid3!=1 || valid4!=1 || success==0) {
    g.remove_newest_vertex(4);
    return 0;
  }

  // delete vi
  g.remove_edges(vi);
  g.remove_edges(vj);
  g.delete_vertex(vi);
  if (vi<vj) {
    g.delete_vertex(vj-1);
  } else g.delete_vertex(vj);

  if (g.check_topologically_connected(g.Nv-1, g.Nv-2)==0)
  {
    g.remove_newest_vertex(4);
    g.restore_local_state(vi_init, vi_neigh_init, vi_nneigh_init);
    g.restore_local_state(vj_init, vj_neigh_init, vj_nneigh_init);
    g.add_edge(g.Nv-1, g.Nv-2);
    g.update_topology();
    return 0;
  }

  // compute the energy difference
  g.update_topology();
  double e_final = g.compute_energy();
  double de = e_final - e_init;

  // accept or reject
  if (gsl_rng_uniform(r) < g.vfuse*g.vfuse/(g.K*g.K)*exp(-de/g.T)) {
    fprintf(stderr, "Fission of %d-%d\n", vi, vj);
    return 1;
  }
  else {
    g.restore_local_state(vi_init, vi_neigh_init, vi_nneigh_init);
    g.restore_local_state(vj_init, vj_neigh_init, vj_nneigh_init);
    g.add_edge(g.Nv-1, g.Nv-2);
    g.update_topology();
    return 0;
  }
}

int insert_monomer_wedge(graph &g, cargo &c, int edge, gsl_rng *r) {
  //g.update_surface_counts();
  int vi = g.e[edge][0];
  int vj = g.e[edge][1];
  // attempt to add edge vk-vl
  int vk = vi;
  int vl = g.get_wedge_pair(vi,vj);
  int vm = vj;
  if (vl==-1) {
    vk = vj;
    vl = g.get_wedge_pair(vj,vi);
    vm = vi;
  }
  if (vl==-1) {
    fprintf(stderr, "Wedge move attempted on non-wedge edge!\n");
    exit(EXIT_FAILURE);
  }
  // check corner case diamond shapes an occupied edge
  for (int ni=0; ni<g.nneigh[vk]; ni++) {
    int vn = g.neigh[vk][ni];
    if (g.check_neigh(vl, vn)==1 && vn!=vm) {
      // ensure that both edges of vl-vn and vk-vn are unocc
      if (g.surface_edge(g.eind[vl][vn])==0
          || g.surface_edge(g.eind[vk][vn])==0) {
        return 0;
      }
    }
  }
  // no sevenfold defects (see also, check wedge)
  if (g.nneigh[vm]>6) {
    return 0;
  }
  c.compute_occupied(g);
  int added = g.add_edge(vk, vl);
  if (added!=1) {
    fprintf(stderr, "Wedge move exceeds length on %d-%d!\n", vk, vl);
    fprintf(stderr, "vi and vj: %d-%d!\n", vi, vj);
    exit(EXIT_FAILURE);
  }
  g.update_topology();
  double de = g.stretch_energy(g.Ne-1);
  de += g.bend_energy(g.eind[vk][vm]);
  de += g.bend_energy(g.eind[vm][vl]);
  if (g.umbrella_sampling==true) {
    de += 0.5*g.kumb*((g.Nm+1-g.Ncenter)*(g.Nm+1-g.Ncenter)-(g.Nm-g.Ncenter)*(g.Nm-g.Ncenter));
  }
  c.compute_occupied(g);
  double dec;
  if (g.t[edge][0]==vl) {
    dec = c.cargo_graph_energy(vi, vj, 0, g);
  } else dec =  c.cargo_graph_energy(vi, vj, 1, g);
  double crit = 2.*g.z*g.K*g.K*g.K*exp((-de-dec)/g.T);
  if (gsl_rng_uniform(r) < crit) {
    g.Nm+=1;
    return 1;
  }
  else {
    g.delete_edge(g.Ne-1);
    g.update_topology();
    c.compute_occupied(g);
    return 0;
  }
}

int insert_monomer_free(graph &g, cargo &c, int edge, gsl_rng *r) {
  //g.update_surface_counts();
  double *vnew = new double[3];
  new_center(g, vnew, edge);
  vnew[0] += (0.5-gsl_rng_uniform(r))*g.xi;
  vnew[1] += (0.5-gsl_rng_uniform(r))*g.xi;
  vnew[2] += (0.5-gsl_rng_uniform(r))*g.xi;
  g.add_vertex(vnew[0],vnew[1],vnew[2]);
  delete[] vnew;

  int vi, vj;
  vi = g.e[edge][0];
  vj = g.e[edge][1];

  if (dist(g.v[vi],g.v[g.Nv-1]) > g.lmax)
  {
    g.v[g.Nv-1][0] = 0;
    g.v[g.Nv-1][1] = 0;
    g.v[g.Nv-1][2] = 0;
    g.Nv -= 1;
    return 0;
  }

  if (dist(g.v[vj],g.v[g.Nv-1]) > g.lmax)
  {
    g.v[g.Nv-1][0] = 0;
    g.v[g.Nv-1][1] = 0;
    g.v[g.Nv-1][2] = 0;
    g.Nv -= 1;
    return 0;
  }

  // attempt to add edges
  int n_added = 0;
  n_added+=g.add_edge(vi,g.Nv-1);
  n_added+=g.add_edge(vj,g.Nv-1);
  if (n_added!=2) {
    g.remove_added(n_added);
    return 0;
  }

  g.update_topology();


  int e1 = g.Ne-1;
  int e2 = g.Ne-2;
  int valid = g.update_edges(g.Nv-1);
  if (valid == -1) {
    g.delete_edge(g.Ne-1);
    g.delete_edge(g.Ne-1);
    g.Nv -= 1;
    g.update_topology();
    return 0;
  }
  if (g.check_overlap(e1)==1 || g.check_overlap(e2)==1) {
    g.delete_edge(g.Ne-1);
    g.delete_edge(g.Ne-1);
    g.Nv -= 1;
    g.update_topology();
    return 0;
  }

  double vp = 4/3.*4*atan(1.)*g.xi*g.xi*g.xi;

  double crit = g.z*g.K*g.K*vp;
  double de = g.vertex_energy(g.Nv-1);
  if (g.umbrella_sampling==true) {
    de += 0.5*g.kumb*((g.Nm+1-g.Ncenter)*(g.Nm+1-g.Ncenter)-(g.Nm-g.Ncenter)*(g.Nm-g.Ncenter));
  }
  c.compute_occupied(g);
  int ti;
  if (g.t[edge][0]==g.Nv-1) {
    ti=0;
  } else ti = 1;
  double dec = c.cargo_graph_energy(vi,vj,ti,g);
  crit *= exp((-de-dec)/g.T);
  if (gsl_rng_uniform(r) < crit) {
    g.Nm+=1;
    return 1;
  }
  else {
    g.delete_edge(e1);
    g.delete_edge(e2);
    g.Nv -= 1;
    g.update_topology();
    c.compute_occupied(g);
    return 0;
  }
}



int delete_monomer(graph &g, cargo &c, int edge, gsl_rng *r)
{
  if (edge<12) {
    return 0;
  }
  int e0 = g.e[edge][0];
  int e1 = g.e[edge][1];
  if (!(g.nneigh[e0]==2 || g.nneigh[e1]==2))
  {
    // two edge monomer
    int vm = g.t[edge][0];
    if (vm==-1) {
      vm = g.t[edge][1];
    }
    if (g.surface_vertex(vm)==1) {
      // vertex connected
      return 0;
    }
    if (g.check_overlap(edge)==1) {
      g.delete_edge(edge);
      g.Nm-=1;
      g.update_topology();
      c.compute_occupied(g);
      return 1;
    }
    double eext = g.bend_energy(g.eind[e0][vm]) + g.bend_energy(g.eind[e1][vm]);
    double eint = g.stretch_energy(edge);
    if (g.umbrella_sampling==true) {
      eint -= 0.5*g.kumb*((g.Nm-1-g.Ncenter)*(g.Nm-1-g.Ncenter)-(g.Nm-g.Ncenter)*(g.Nm-g.Ncenter));
    }
    vm = g.t[edge][0];
    if (vm==-1) {
      vm = g.t[edge][1];
    }
    int ti;
    if (g.t[g.eind[e0][vm]][0]==e1) {
      ti = 0;
    } else ti = 1;
    double dec = c.cargo_graph_energy(e0,vm,ti,g);
    double crit = exp((eint+eext+dec)/g.T)/(2.0*g.z*g.K*g.K*g.K);
    if (gsl_rng_uniform(r) < crit) {
      g.delete_edge(edge);
      g.Nm-=1;
      g.update_topology();
      c.compute_occupied(g);
      return 1;
    }
    else return 0;
  }
  int edel;
  if (g.nneigh[e0]==2) {
    edel = e0;
  } else edel = e1;

  if (g.check_overlap(edge)==1) {
    while (g.nneigh[edel]>0) {
      g.delete_edge(g.eind[edel][g.neigh[edel][0]]);
    }
    g.delete_vertex(edel);
    if (g.nneigh[e1]==1) {
      g.delete_edge(g.eind[e1][g.neigh[e1][0]]);
      g.delete_vertex(e1);
    }
    g.Nm-=1;
    g.update_topology();
    c.compute_occupied(g);
    return 1;
  }


  double de = g.vertex_energy(edel);
  if (g.umbrella_sampling==true) {
    de -= 0.5*g.kumb*((g.Nm-1-g.Ncenter)*(g.Nm-1-g.Ncenter)-(g.Nm-g.Ncenter)*(g.Nm-g.Ncenter));
  }
  double vp = 4/3.*4*atan(1.)*g.xi*g.xi*g.xi;
  int vm = g.t[edge][0];
  if (vm==-1) {
    vm = g.t[edge][1];
  }
  int ti;
  if (g.t[g.eind[e0][vm]][0]==e1) {
    ti = 0;
  } else ti = 1;
  double dec = c.cargo_graph_energy(e0, vm, ti, g);
  double crit = exp((de+dec)/g.T)/(2*vp*g.z*g.K*g.K);
  if (gsl_rng_uniform(r)<crit) {
    while (g.nneigh[edel]>0) {
      g.delete_edge(g.eind[edel][g.neigh[edel][0]]);
    }
    g.delete_vertex(edel);
    if (g.nneigh[e1]==1) {
      g.delete_edge(g.eind[e1][g.neigh[e1][0]]);
      g.delete_vertex(e1);
    }
    g.Nm-=1;
    g.update_topology();
    c.compute_occupied(g);
    return 1;
  }
  else {
    return 0;
  }
}

void new_center(graph &g, double *vnew, int edge) {
  int vi, vj, vk;
  double h=1.0, w = 0.5, theta = 0.;
  vi = g.e[edge][0];
  vj = g.e[edge][1];
  vk = g.t[edge][0];
  if (vk == -1) vk = g.t[edge][1];
  double x,y,z,vx,vy,vz,ex,ey,ez,nm,vxr,vyr,vzr;
  // random direction in the plane
  x = g.v[vi][0] + w *(g.v[vj][0]-g.v[vi][0]);
  y = g.v[vi][1] + w *(g.v[vj][1]-g.v[vi][1]);
  z = g.v[vi][2] + w *(g.v[vj][2]-g.v[vi][2]);
  // random length for the displacement
  vx = -h*(g.v[vk][0] - x);
  vy = -h*(g.v[vk][1] - y);
  vz = -h*(g.v[vk][2] - z);
  // random bending angle for the new triangle
  // rotate the vector around the edge by theta
  ex = g.v[vj][0] - g.v[vi][0];
  ey = g.v[vj][1] - g.v[vi][1];
  ez = g.v[vj][2] - g.v[vi][2];
  nm = sqrt(ex*ex + ey*ey + ez*ez);
  ex /= nm;
  ey /= nm;
  ez /= nm;
  double ct = cos(theta);
  double st = sin(theta);
  double mct = 1-cos(theta);
  vxr = vx*(ct+ex*ex*mct) + vy*(ex*ey*mct-ez*st) + vz*(ex*ez*mct+ey*st);
  vyr = vx*(ey*ex*mct+ez*st) + vy*(ct+ey*ey*mct) + vz*(ey*ez*mct-ex*st);
  vzr = vx*(ez*ex*mct-ey*st) + vy*(ez*ey*mct+ex*st) + vz*(ct+ez*ez*mct);
  vx = x+vxr;
  vy = y+vyr;
  vz = z+vzr;
  vnew[0] = vx;
  vnew[1] = vy;
  vnew[2] = vz;
}
