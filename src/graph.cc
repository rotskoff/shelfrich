/* graph.cc Implements elastic shell structure.
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
#include <stdlib.h>
#include <math.h>
#include <set>
#include <queue>
#include <vector>
#include <algorithm>
#include "graph.h"

const double pi = 4.0*atan(1.0);

graph::graph()
{
  Nv = 0;
  Ne = 0;
  Nm = 0;
  Nperim = 0;
  Nfusion = 0;
  Nfission = 0;
  Nvmax = 0;
  Nemax = 0;
  v = NULL;
  vcs = NULL;
  e = NULL;
  esurf = NULL;
  etype = NULL;
  l = NULL;
  eind = NULL;
  energy = 0;
  z = 0;
  K = 0;
  xi = 0;
  epsilon = 0;
  l0 = 0;
  kappa = NULL;
  lmax = 0;
  lmin = 0;
  lfuse = 0;
  vfuse = 0;
  umbrella_sampling = false;
  Ncenter = 0;
  kumb= 0;
}

graph::~graph()
{
  if (Nv>0) {
    for (int i=0; i<Nvmax; i++)
    {
      delete[] v[i];
      delete[] eind[i];
      delete[] neigh[i];
      delete[] fusion_neigh[i];
    }
    delete[] v;
    delete[] vcs;
    delete[] neigh;
    delete[] fusion_neigh;
    delete[] nneigh;
    delete[] fusion_nneigh;
    delete[] eind;
  }
  if (Ne>0) {
    for (int i=0; i<Nemax; i++)
    {
      delete[] e[i];
      delete[] n[i][0];
      delete[] n[i][1];
      delete[] tri[i][0];
      delete[] tri[i][1];
      delete[] n[i];
      delete[] tri[i];
      delete[] t[i];
      delete[] evec[i];
    }
    delete[] e;
    delete[] n;
    delete[] t;
    delete[] tri;
    delete[] evec;
    delete[] esurf;
    delete[] etype;
    delete[] l;
  }
  if (Ntype>0) {
    delete[] kappa;
    delete[] min_bond;
  }
}

void graph::initialize(int Nvmax0, int Nemax0, int Ntype0)
{
  Nvmax = Nvmax0;
  Nemax = Nemax0;
  Ntype = Ntype0;
  v = new double*[Nvmax];
  vcs = new int[Nvmax];
  e = new int*[Nemax];
  eind = new int*[Nvmax];
  t = new int*[Nemax];
  tri = new int**[Nemax];
  n = new double**[Nemax];
  neigh = new int*[Nvmax];
  fusion_neigh = new int*[Nvmax];
  evec = new double*[Nemax];
  nneigh = new int[Nvmax];
  fusion_nneigh = new int[Nvmax];
  for (int i=0; i<Nvmax; i++)
  {
    v[i] = new double[3];
    eind[i] = new int[Nvmax];
    neigh[i] = new int[10];
    fusion_neigh[i] = new int[10];
    nneigh[i] = 0;
    fusion_nneigh[i] = 0;
  }
  l = new double[Nemax];
  for (int i=0; i<Nemax; i++)
  {
    e[i] = new int[2];
    evec[i] = new double[3];
    t[i] = new int[2];
    t[i][0] = -1;
    t[i][1] = -1;
    tri[i] = new int*[2];
    tri[i][0] = new int[3];
    tri[i][1] = new int[3];
    n[i] = new double*[2];
    n[i][0] = new double[3];
    n[i][1] = new double[3];
    l[i] = 0;
  }
  esurf = new int[Nemax];
  etype = new int[Nemax];
  kappa = new double[Ntype];
  min_bond = new double[Ntype];
}

void graph::add_vertex(double x, double y, double z)
{
  if (Nv<Nvmax) {
    v[Nv][0] = x;
    v[Nv][1] = y;
    v[Nv][2] = z;
    nneigh[Nv] = 0;
    Nv+=1;
  } else {
    fprintf(stderr,"Data structure cannot accommodate another vertex.\n");
    exit(EXIT_FAILURE);
  }
}

void graph::delete_vertex(int vi)
{
  if (nneigh[vi]>0) {
    fprintf(stderr, "Error: cannot delete vertex with neighbors.\n");
    exit(-1);
  }

  // shift all vertices above vi
  for(int vj=vi; vj<Nv; vj++)
  {
    v[vj][0] = v[vj+1][0];
    v[vj][1] = v[vj+1][1];
    v[vj][2] = v[vj+1][2];
    vcs[vj] = vcs[vj+1];
    // update the index
    nneigh[vj] = nneigh[vj+1];
    for (int ni=0; ni<nneigh[vj]; ni++) {
      if (neigh[vj+1][ni]>vi) {
        neigh[vj][ni] = neigh[vj+1][ni]-1;
      } else neigh[vj][ni] = neigh[vj+1][ni];
    }
  }

  // update edges involving the other vertexes
  for(int vj=0; vj<vi; vj++)
  {
    for (int ni=0; ni<nneigh[vj]; ni++) {
      if (neigh[vj][ni]>vi) {
        neigh[vj][ni] = neigh[vj][ni]-1;
      } else neigh[vj][ni] = neigh[vj][ni];
    }
  }

  Nv-=1;

  // adjust edges
  for (int edge=0; edge<Ne; edge++) {
    if (e[edge][0] > vi) e[edge][0]-=1;
    if (e[edge][1] > vi) e[edge][1]-=1;
    if (t[edge][0] > vi) t[edge][0]-=1;
    if (t[edge][1] > vi) t[edge][1]-=1;
    int v1,v2;
    v1 = e[edge][0];
    v2 = e[edge][1];
    eind[v1][v2] = edge;
    eind[v2][v1] = edge;
  }
  update_t();
}

int graph::add_edge(int v1, int v2)
{
  // only add the edge if it is new
  int connected = 0;
  for (int ni=0; ni<nneigh[v1]; ni++) {
    if (neigh[v1][ni]==v2) {
      connected=1;
      return 0;
    }
  }
  if (dist(v[v1],v[v2])>lmax || dist(v[v1],v[v2])<lmin) {
    return 0;
  }

  if (connected == 0) {
    if (Ne<Nemax) {
      e[Ne][0] = v1;
      e[Ne][1] = v2;
      eind[v1][v2] = Ne;
      eind[v2][v1] = Ne;
      etype[Ne] = 0; // default type
      evec[Ne][0] = v[v2][0] - v[v1][0];
      evec[Ne][1] = v[v2][1] - v[v1][1];
      evec[Ne][2] = v[v2][2] - v[v1][2];
      l[Ne] = sqrt ( evec[Ne][0]*evec[Ne][0] +
                     evec[Ne][1]*evec[Ne][1] +
                     evec[Ne][2]*evec[Ne][2] );
      neigh[v1][nneigh[v1]] = v2;
      nneigh[v1]++;
      neigh[v2][nneigh[v2]] = v1;
      nneigh[v2]++;
      Ne+=1;
      return 1;
    } else {
      fprintf(stderr,"Ne is %d\n",Ne);
      fprintf(stderr,"Data structure cannot accommodate another edge.\n");
      exit(EXIT_FAILURE);
    }
  }
  return 0;
}

int graph::force_add_edge(int v1, int v2)
{
  // only add the edge if it is new
  int connected = 0;
  for (int ni=0; ni<nneigh[v1]; ni++) {
    if (neigh[v1][ni]==v2) {
      connected=1;
      return 0;
    }
  }

  if (connected == 0) {
    if (Ne<Nemax) {
      e[Ne][0] = v1;
      e[Ne][1] = v2;
      eind[v1][v2] = Ne;
      eind[v2][v1] = Ne;
      etype[Ne] = 0; // default type
      evec[Ne][0] = v[v2][0] - v[v1][0];
      evec[Ne][1] = v[v2][1] - v[v1][1];
      evec[Ne][2] = v[v2][2] - v[v1][2];
      l[Ne] = sqrt ( evec[Ne][0]*evec[Ne][0] +
                     evec[Ne][1]*evec[Ne][1] +
                     evec[Ne][2]*evec[Ne][2] );
      neigh[v1][nneigh[v1]] = v2;
      nneigh[v1]++;
      neigh[v2][nneigh[v2]] = v1;
      nneigh[v2]++;
      Ne+=1;
      return 1;
    } else {
      fprintf(stderr,"Ne is %d\n",Ne);
      fprintf(stderr,"Data structure cannot accommodate another edge.\n");
      exit(EXIT_FAILURE);
    }
  }
  return 0;
}

void graph::delete_edge(int edge)
{
  int vi, vj;
  vi = e[edge][0];
  vj = e[edge][1];
  etype[edge]=0;
  // shift all the arrays
  for (int i=edge; i<Ne; i++) {
    e[i][0] = e[i+1][0];
    e[i][1] = e[i+1][1];
    etype[i] = etype[i+1];
    eind[e[i][0]][e[i][1]] = i;
    eind[e[i][1]][e[i][0]] = i;
    t[i][0] = t[i+1][0];
    t[i][1] = t[i+1][1];
    tri[i][0][0] = tri[i+1][0][0];
    tri[i][0][1] = tri[i+1][0][1];
    tri[i][0][2] = tri[i+1][0][2];
    tri[i][1][0] = tri[i+1][1][0];
    tri[i][1][1] = tri[i+1][1][1];
    tri[i][1][2] = tri[i+1][1][2];
    n[i][0][0] = n[i+1][0][0];
    n[i][0][1] = n[i+1][0][1];
    n[i][0][2] = n[i+1][0][2];
    n[i][1][0] = n[i+1][1][0];
    n[i][1][1] = n[i+1][1][1];
    n[i][1][2] = n[i+1][1][2];
    evec[i][0] = evec[i+1][0];
    evec[i][1] = evec[i+1][1];
    evec[i][2] = evec[i+1][2];
    l[i] = l[i+1];
  }

  // remove vj from the neighbors of vi
  int i=0;
  int nindex = 0;
  while (i<nneigh[vi]) {
    if (neigh[vi][i] == vj) {
      nindex = i;
      break;
    } else i++;
  }
  for (int i=nindex; i<nneigh[vi]-1; i++) {
    neigh[vi][i] = neigh[vi][i+1];
  }
  nneigh[vi] -= 1;

  // remove vi from the neighbors of vj
  i=0;
  nindex = 0;
  while (i<nneigh[vj]) {
    if (neigh[vj][i] == vi) {
      nindex = i;
      break;
    } else i++;
  }
  for (int i=nindex; i<nneigh[vj]-1; i++) {
    neigh[vj][i] = neigh[vj][i+1];
  }
  nneigh[vj] -= 1;

  // decrease the number of edges
  Ne-=1;
  update_t();
}

int graph::update_edges(int vi)
{
  int edge, oedge, vj;
  for (int j=0; j<nneigh[vi]; j++)
  {
    vj = neigh[vi][j];
    edge = eind[vi][vj];
    evec[edge][0] = v[e[edge][1]][0] - v[e[edge][0]][0];
    evec[edge][1] = v[e[edge][1]][1] - v[e[edge][0]][1];
    evec[edge][2] = v[e[edge][1]][2] - v[e[edge][0]][2];
    l[edge] = norm(evec[edge]);
    compute_normals(edge);
    // check the edge length constraint
    if (l[edge]<lmin || l[edge]>lmax) {
      return -1;
    }
  }
  for (int j=0; j<nneigh[vi]; j++) {
    vj = neigh[vi][j];
    edge = eind[vi][vj];
    // update normals for the triangles specified by opposite edges
    if (t[edge][0] != -1) {
      oedge = eind[vj][t[edge][0]];
      compute_normals(oedge);
    }
    if (t[edge][1] != -1) {
      oedge = eind[vj][t[edge][1]];
      compute_normals(oedge);
    }
  }
  return 1;
}

void graph::recenter()
{
  double ox,oy,oz;
  ox = v[0][0];
  oy = v[0][1];
  oz = v[0][2];
  for (int vi=0; vi<Nv; vi++)
  {
    v[vi][0] -= ox;
    v[vi][1] -= oy;
    v[vi][2] -= oz;
  }
}

void graph::center_at(double x, double y, double z)
{
  double ox,oy,oz;
  ox = v[0][0];
  oy = v[0][1];
  oz = v[0][2];
  for (int vi=0; vi<Nv; vi++)
  {
    v[vi][0] -= ox-x;
    v[vi][1] -= oy-y;
    v[vi][2] -= oz-z;
  }
}

int graph::surface_vertex(int vi)
{
  // loop over edges containing the vertex and check if they're on
  // the surface by checking the triangles
  for (int ni = 0; ni<nneigh[vi]; ni++) {
    int vj = neigh[vi][ni];
    int edge = eind[vi][vj];
    if ((t[edge][0] == -1) || (t[edge][1] == -1)) {
      return 1;
    }
  }
  return 0;
}

int graph::surface_edge(int edge)
{
  if ((t[edge][0] == -1) || (t[edge][1] == -1)) {
    return 1;
  } else return 0;
}

int graph::check_neigh(int vi, int vj)
{
  for (int ni=0; ni<nneigh[vi]; ni++)
  {
    if (neigh[vi][ni]==vj) {
      return 1;
    }
  }
  return 0;
}

int graph::check_next_neigh(int vi, int vj)
{
  for (int ni=0; ni<nneigh[vi]; ni++)
  {
    int vn = neigh[vi][ni];
    for (int nj=0; nj<nneigh[vj]; nj++)
    {
      if (neigh[vj][nj]==vn) {
        return 1;
      }
    }
  }
  return 0;
}

int graph::check_adjacent(int vi, int vj)
{
  int vk, vl, edge, oedge;
  for (int ni=0; ni<nneigh[vi]; ni++)
  {
    vk = neigh[vi][ni];
    edge = eind[vi][vk];
    vl = t[edge][0];
    oedge = eind[vk][vl];
    if (t[oedge][0]==vj || t[oedge][1]==vj)
    {
      return 1;
    }
    if ( check_neigh(t[oedge][0],vj)==1 ||
         check_neigh(t[oedge][1],vj)==1  )
    {
      return 1;
    }
  }
  return 0;
}

int graph::check_defect(int vi) {
  if (surface_vertex(vi)==0 && (nneigh[vi]<5 || nneigh[vi]>7)) {
    return 0;
  } else return 1;
}

void graph::update_t()
{
  int tcount;
  int v1, v2;
  Nperim = 0;
  for (int edge=0; edge<Ne; edge++)
  {
    v1 = e[edge][0];
    v2 = e[edge][1];
    tcount = 0;
    for (int i=0; i<nneigh[v1]; i++) {
      for (int j=0; j<nneigh[v2]; j++) {
        if (neigh[v1][i] == neigh[v2][j])
        {
          t[edge][tcount] = neigh[v1][i];
          tcount++;
        }
      }
    }
    // if on the surface, other triangle doesn't exist
    if (tcount == 1) {
      t[edge][tcount] = -1;
      esurf[Nperim] = edge;
      Nperim++;
    }
  }
}

int graph::check_wedge(int ei)
{
  int vi = e[ei][0];
  int vj = e[ei][1];
  int vk;
  for (int ni=0; ni<nneigh[vi]; ni++) {
    vk = neigh[vi][ni];
    if (check_neigh(vj,vk)==0 && surface_edge(eind[vi][vk])==1) {
      if (check_adjacent(vj,vk)==0 && dist(v[vj],v[vk])<lmax
          && dist(v[vj],v[vk])>lmin && nneigh[vi]<7) {
        return 1;
      }
    }
  }
  for (int nj=0; nj<nneigh[vj]; nj++) {
    vk = neigh[vj][nj];
    if (check_neigh(vi,vk)==0 && surface_edge(eind[vj][vk])==1) {
      if (check_adjacent(vi,vk)==0 && dist(v[vi], v[vk])<lmax
          && dist(v[vi],v[vk])>lmin && nneigh[vj]<7) {
        return 1;
      }
    }
  }
  return 0;
}

int graph::get_wedge_pair(int vi, int vm)
{
  int vk;
  for (int nm=0; nm<nneigh[vm]; nm++) {
    vk = neigh[vm][nm];
    if (check_neigh(vi,vk)==0 && surface_edge(eind[vm][vk])==1) {
      if (check_adjacent(vi,vk)==0 && dist(v[vi],v[vk])<lmax
          && dist(v[vi],v[vk])>lmin) {
        return vk;
      }
    }
  }
  return -1;
}

int graph::check_fusion_vertex_pair(int vi, int vj)
{
  if (surface_vertex(vi)==0 || surface_vertex(vj)==0)
  {
    return 0;
  }
  if (vi!=vj && check_neigh(vi,vj)==0 && check_adjacent(vi,vj)==0)
  {
    if (dist(v[vi],v[vj]) < lfuse) {
      return 1;
    } else return 0;
  } else return 0;
}

void graph::add_fusion_pair(int vi, int vj) {
  pair< set< pair<int,int> >::iterator, bool > success;
  if (check_fusion_vertex_pair(vi, vj)==1) {
    pair<int, int> pij(vi,vj);
    success = vfusion.insert(pij);
    if (success.second) {
      Nfusion+=1;
      fusion_neigh[vi][fusion_nneigh[vi]] = vj;
      fusion_nneigh[vi] += 1;
      fusion_neigh[vj][fusion_nneigh[vj]] = vi;
      fusion_nneigh[vj] += 1;
    }
    /*
    success = vfusion.insert(vj);
    if (success.second) {
      Nfusion+=1;
      fusion_neigh[vj][fusion_nneigh[vj]] = vi;
      fusion_nneigh[vj] += 1;
    }
    */
  }
}

void graph::update_fusion_set()
{
  // reset the nneigh array; fusion_neigh doesn't need to be reset
  // because the nneigh array will be determine which values are accessed
  Nfusion = 0;
  update_topology();
  fill(fusion_nneigh, fusion_nneigh+Nvmax, 0);
  vfusion.clear();
  for (int i=0; i<Nperim; i++) {
    int ei = esurf[i];
    int vi0 = e[ei][0];
    int vi1 = e[ei][1];
    for (int j=i+1; j<Nperim; j++) {
      int ej = esurf[j];
      int vj0 = e[ej][0];
      int vj1 = e[ej][1];
      add_fusion_pair(vi0, vj0);
      add_fusion_pair(vi0, vj1);
      add_fusion_pair(vi1, vj0);
      add_fusion_pair(vi1, vj1);
    }
  }
}

void graph::save_local_state(int vi, vector<double> &vi_init, vector<int> &vi_neigh_init, int &vi_nneigh_init)
{
  vi_init.insert(vi_init.end(), &v[vi][0], &v[vi][3]);
  vi_neigh_init.insert(vi_neigh_init.end(), &neigh[vi][0], &neigh[vi][10]);
  vi_nneigh_init = nneigh[vi];
}

void graph::save_local_state_fission(int vj, int vi, vector<double> &vi_init, vector<int> &vi_neigh_init, int &vi_nneigh_init)
{
  vi_init.insert(vi_init.end(), &v[vi][0], &v[vi][3]);
  for (int ni=0; ni<nneigh[vi]; ni++)
  {
    if (neigh[vi][ni]!=vj)
    {
      vi_neigh_init.push_back(neigh[vi][ni]);
    }
  }
  vi_nneigh_init = nneigh[vi]-1;
}

void graph::adjust_indices(int vi, int vj, int vi_nneigh_init, vector<int> &vi_neigh_init)
{
  int vn;
  for (int ni=0; ni<vi_nneigh_init; ni++) {
    vn = vi_neigh_init[ni];
    if (vn>vi) {
      vi_neigh_init[ni]-=1;
    }
    if (vn>vj) {
      vi_neigh_init[ni]-=1;
    }
  }
}



void graph::adjust_indices_edge_fusion(int vi, int vi1, int vj, int vj1,
                                       vector<int> &vi_neigh_init)
{
  int ic=0;
  for (int ni=0; ni<nneigh[vi]; ni++) {
    if (neigh[vi][ni]!=vi1) {
      vi_neigh_init[ic] = neigh[vi][ni];
      if (neigh[vi][ni]>vi) {
        vi_neigh_init[ic]-=1;
      }
      if (neigh[vi][ni]>vj) {
        vi_neigh_init[ic]-=1;
      }
      if (neigh[vi][ni]>vi1) {
        vi_neigh_init[ic]-=1;
      }
      if (neigh[vi][ni]>vj1) {
        vi_neigh_init[ic]-=1;
      }
      ic+=1;
    }
  }
}

void graph::restore_local_state(vector<double> &vi_init, vector<int> &vi_neigh_init, int &vi_nneigh_init)
{
  int success = 0;
  add_vertex(vi_init[0],vi_init[1],vi_init[2]);
  for (int ni=0; ni<vi_nneigh_init; ni++) {
    success += add_edge(Nv-1, vi_neigh_init[ni]);
  }
  if (success != vi_nneigh_init) {
    dump_xml(0);
    fprintf(stderr, "Cannot restore\n");
    exit(-1);
  }
}

void graph::add_midpoint_vertex(int vi, int vj)
{
  double x = v[vi][0] + 0.5*(v[vj][0]-v[vi][0]);
  double y = v[vi][1] + 0.5*(v[vj][1]-v[vi][1]);
  double z = v[vi][2] + 0.5*(v[vj][2]-v[vi][2]);
  add_vertex(x,y,z);
}

int graph::check_distance_constraints(int vi)
{
  double d;
  int vn;
  for (int ni=0; ni<nneigh[vi]; ni++)
  {
    vn = neigh[vi][ni];
    d = dist(v[Nv-1],v[vn]);
    if (d<lmin || d>lmax) {
      return 0;
    }
  }
  return 1;
}

int graph::add_new_edges(int vi, int vtarget, int &n_added)
{
  int vk;
  for (int ni=0; ni<nneigh[vi]; ni++)
  {
    vk = neigh[vi][ni];
    n_added += add_edge(vtarget,vk);
  }
  if (n_added == nneigh[vi])
  {
    return 1;
  } else return 0;
}

void graph::update_fission_set()
{
  Nfission = 0;
  update_topology();
  vfission.clear();
  pair< set< pair <int,int> >::iterator, bool> success;
  for (int i=0; i<Nperim; i++) {
    int ei = esurf[i];
    int vi = e[ei][0];
    int vj = e[ei][1];
    {
      for (int ni=0; ni<nneigh[vi]; ni++)
      {
        int vk = neigh[vi][ni];
        if (check_fission_edge(vi, vk)==1) {
          pair <int, int> pik(vi,vk);
          success = vfission.insert(pik);
          Nfission += success.second;
        }
      }
      for (int nj=0; nj<nneigh[vj]; nj++)
      {
        int vk = neigh[vj][nj];
        if (check_fission_edge(vj, vk)==1) {
          if (vk!=vi) { // don't overcount edge fission
            pair <int, int> pjk(vj,vk);
            success = vfission.insert(pjk);
            Nfission += success.second;
          }
        }
      }
    }
  }
}

int graph::check_fission_vertex(int vi)
{
  for (int ni=0; ni<nneigh[vi]; ni++)
  {
    int vj = neigh[vi][ni];
    if (check_fission_edge(vi, vj)==1) {
      return 1;
    }
  }
  return 0;
}

int graph::check_fission_edge(int vi, int vj)
{
  if (surface_vertex(vi)==1 && surface_vertex(vj)==0) {
    return 1;
  }
  else if (surface_vertex(vi)==1 && surface_vertex(vj)==1) {
    int eij = eind[vi][vj];
    if (surface_edge(eij)==0)
    {
      int vk = t[eij][0], vl = t[eij][1];
      if (check_topologically_connected_blocking(vi, vj, vk, vl)==0) {
        return 0;
      }
      else return 1;
    }
    else return 0;
  }
  else return 0;
}

int graph::select_random_fission_partner(int vi, gsl_rng *r)
{
  if (check_fission_vertex(vi)==0) {
    fprintf(stderr, "Attempting to select a fission partner for invalid vi.");
    exit(-1);
  }
  int vj = neigh[vi][gsl_rng_uniform_int(r,nneigh[vi])];
  while (check_fission_edge(vi,vj)==0) {
    vj = neigh[vi][gsl_rng_uniform_int(r,nneigh[vi])];
  }
  return vj;
}

void graph::add_fission_vertices(int vi, gsl_rng *r)
{
  double rad, theta, phi, dx, dy, dz;
  rad = lfuse * gsl_rng_uniform(r);
  theta = pi * gsl_rng_uniform(r);
  phi = 2 * pi * gsl_rng_uniform(r);
  dx = rad * sin(theta) * cos(phi);
  dy = rad * sin(theta) * sin(phi);
  dz = rad * cos(theta);
  // choose a point in the sphere of radius lfuse
  // find the other new vertex v2
  double v1x, v1y, v1z, v2x, v2y, v2z;
  v1x = v[vi][0]+dx;
  v1y = v[vi][1]+dy;
  v1z = v[vi][2]+dz;
  v2x = v[vi][0]-dx;
  v2y = v[vi][1]-dy;
  v2z = v[vi][2]-dz;
  add_vertex(v1x, v1y, v1z);
  add_vertex(v2x, v2y, v2z);
}

int graph::connect_neighs(int vc, int vi, int vj, int vnew)
{
  int success = 1;
  int vt = vj;
  int vk = vc;
  int count = 0, vkn = 0, vl;
  while (vk!=-1) {
    success *= add_edge(vk, vnew);
    vl = vk;
    vkn = t[eind[vi][vk]][0];
    if (vkn==vt) {
      vk = t[eind[vi][vk]][1];
    } else vk = vkn;
    vt = vl;
    count+=1;
    if (success==0){
      return 0;
    }
    if (count > 8) {
      // catches malformed structures and exits
      fprintf(stderr, "%d %d %d %d\n", vc, vi, vj, vnew);
      fprintf(stderr, "Fission edge addition failed.\n");
      dump_xml(1);
      exit(EXIT_FAILURE);
    }
  }
  return success;
}

int graph::check_topologically_connected(int vi, int vj)
{
  queue<int> pathq;
  set<int> visited;
  pathq.push(vi);
  int vn, vnn;
  while (pathq.empty()==false) {
    vn = pathq.front();
    pathq.pop();
    if (vn==vj) {
      return 1;
    } else {
      visited.insert(vn);
      for (int ni=0; ni<nneigh[vn]; ni++)
      {
        vnn = neigh[vn][ni];
        if (visited.count(vnn)==0) {
          pathq.push(neigh[vn][ni]);
        }
      }
    }
  }
  return 0;
}

int graph::check_topologically_connected_blocking(int vi, int vj, int vk, int vl)
{
  queue<int> pathq;
  set<int> visited;
  pathq.push(vi);
  visited.insert(vk);
  visited.insert(vl);
  int vn, vnn;
  while (pathq.empty()==false) {
    vn = pathq.front();
    pathq.pop();
    if (vn==vj) {
      return 1;
    } else {
      visited.insert(vn);
      for (int ni=0; ni<nneigh[vn]; ni++)
      {
        vnn = neigh[vn][ni];
        if (visited.count(vnn)==0) {
          pathq.push(neigh[vn][ni]);
        }
      }
    }
  }
  return 0;
}

void graph::restore_local_state_fission(int vj, vector<double> &vi_init, vector<int> &vi_neigh_init, int &vi_nneigh_init)
{
  add_vertex(vi_init[0],vi_init[1],vi_init[2]);
  for (int ni=0; ni<vi_nneigh_init; ni++) {
    if (vi_neigh_init[ni]!=vj) {
      fprintf(stderr, "restoring neigh %d\n", vi_neigh_init[ni]);
      force_add_edge(Nv-1, vi_neigh_init[ni]);
    }
  }
}

void graph::remove_edges(int vi)
{
  while (nneigh[vi]>0) {
    delete_edge(eind[vi][neigh[vi][0]]);
  }
}

void graph::remove_newest_vertex(int n_times)
{
  for (int i=0; i<n_times; i++)
  {
    remove_edges(Nv-1);
    delete_vertex(Nv-1);
  }
}

void graph::remove_added(int n_added)
{
  for (int nj=0; nj<n_added; nj++)
  {
    delete_edge(Ne-1);
  }
  delete_vertex(Nv-1);
}

void graph::get_normal(int v1, int v2, int v3, double *n) {
  double *ev1 = new double[3];
  double *ev2 = new double[3];
  vsub(v[v1], v[v2], ev1);
  vsub(v[v1], v[v3], ev2);
  cross(ev1, ev2, n);
  double nm = norm(n);
  n[0] /= nm;
  n[1] /= nm;
  n[2] /= nm;
  delete[] ev1;
  delete[] ev2;
}

void graph::compute_normals(int i)
{
  /* v3 and v4 need not both be defined */
  if (t[i][0]==-1 || t[i][1]==-1)
  {
    n[i][0][0] = 0;
    n[i][0][1] = 0;
    n[i][0][2] = 0;
    n[i][1][0] = 0;
    n[i][1][1] = 0;
    n[i][1][2] = 0;
    return;
  }
  else {
    // set the first normal
    int v1 = tri[i][0][0];
    int v2 = tri[i][0][1];
    int v3 = tri[i][0][2];
    if (v1==0 && v2==0 && v3==0) {
      get_normal(e[i][0], e[i][1], t[i][0], n[i][0]);
    } else {
      get_normal(v1, v2, v3, n[i][0]);
    }
    // second normal
    v1 = tri[i][1][0];
    v2 = tri[i][1][1];
    v3 = tri[i][1][2];
    if (v1==0 && v2==0 && v3==0) {
      get_normal(e[i][0], e[i][1], t[i][1], n[i][1]);
    } else {
      get_normal(v1, v2, v3, n[i][1]);
    }
  }
}

void graph::set_orientation(int edge, int ti, int v1, int v2, int v3)
{
  tri[edge][ti][0] = v1;
  tri[edge][ti][1] = v2;
  tri[edge][ti][2] = v3;
  int e0 = e[edge][0];
  int e1 = e[edge][1];
  int t0 = t[edge][ti];
  int et1 = eind[e0][t0];
  int i=0;
  if (t[et1][i]!=e1) { i=1; }
  tri[et1][i][0] = v1;
  tri[et1][i][1] = v2;
  tri[et1][i][2] = v3;
  int et2 = eind[e1][t0];
  i=0;
  if (t[et2][i]!=e0) { i=1; }
  tri[et2][i][0] = v1;
  tri[et2][i][1] = v2;
  tri[et2][i][2] = v3;
}

int graph::check_oriented(int ei) {
  if (surface_edge(ei)==0) {
    if ((tri[ei][0][0]==0 && tri[ei][0][1]==0 && tri[ei][0][2]==0) ||
        (tri[ei][1][0]==0 && tri[ei][1][1]==0 && tri[ei][1][2]==0) ) {
      return 0;
    }
  }
  else {
    if ((tri[ei][0][0]==0 && tri[ei][0][1]==0 && tri[ei][0][2]==0) &&
        (tri[ei][1][0]==0 && tri[ei][1][1]==0 && tri[ei][1][2]==0) ) {
      return 0;
    }
  }
  return 1;
}

void graph::orient()
{
  // propagate the orientation
  int v1, v2, v3;
  for (int i=0; i<Ne; i++) {
    tri[i][0][0] = 0;
    tri[i][0][1] = 0;
    tri[i][0][2] = 0;
    tri[i][1][0] = 0;
    tri[i][1][1] = 0;
    tri[i][1][2] = 0;
  }
  set_orientation(0,0,e[0][0], e[0][1], t[0][0]);
  queue<int> orientq;
  for (int ei=0; ei<Ne; ei++) {
    if (check_oriented(ei)==0) {
      orientq.push(ei);
    }
  }

  int nattempt=0;
  int i;
  while (!orientq.empty()) {
    i = orientq.front();
    if (!(tri[i][0][0]==0 && tri[i][0][1]==0 && tri[i][0][2]==0)
        && (tri[i][1][0]==0 && tri[i][1][1]==0 && tri[i][1][2]==0)) {
      v1 = e[i][0];
      v2 = e[i][1];
      v3 = t[i][1];
      if (v3==-1) {
        orientq.pop();
      }
      else if ((tri[i][0][0]==v1 && tri[i][0][1]==v2) ||
               (tri[i][0][1]==v1 && tri[i][0][2]==v2) ||
               (tri[i][0][2]==v1 && tri[i][0][0]==v2) )
      {
        set_orientation(i, 1, v2, v1, v3);
        orientq.pop();
      }
      else if ((tri[i][0][0]==v2 && tri[i][0][1]==v1) ||
               (tri[i][0][1]==v2 && tri[i][0][2]==v1) ||
               (tri[i][0][2]==v2 && tri[i][0][0]==v1) )
      {
        set_orientation(i, 1, v1, v2, v3);
        orientq.pop();
      }
    }
    else if (!(tri[i][1][0]==0 && tri[i][1][1]==0 && tri[i][1][2]==0)
             && (tri[i][0][0]==0 && tri[i][0][1]==0 && tri[i][0][2]==0)) {
      v1 = e[i][0];
      v2 = e[i][1];
      v3 = t[i][0];
      if (v3==-1) {
        orientq.pop();
      }
      else if ((tri[i][1][0]==v1 && tri[i][1][1]==v2) ||
               (tri[i][1][1]==v1 && tri[i][1][2]==v2) ||
               (tri[i][1][2]==v1 && tri[i][1][0]==v2) ) {
        set_orientation(i, 0, v2, v1, v3);
        orientq.pop();
      }
      else if ((tri[i][1][0]==v2 && tri[i][1][1]==v1) ||
               (tri[i][1][1]==v2 && tri[i][1][2]==v1) ||
               (tri[i][1][2]==v2 && tri[i][1][0]==v1) )
      {
        set_orientation(i, 0, v1, v2, v3);
        orientq.pop();
      }
    }
    else {
      if (check_oriented(i)==0) {
        orientq.push(i);
        orientq.pop();
      }
      else {
        orientq.pop();
      }
    }
    nattempt++;
    if (nattempt==100000) {
      fprintf(stderr,"Failed to orient, Queue has %ld items\n",orientq.size());
      dump_xml(1);
      exit(EXIT_FAILURE);
    }
  }
}

void graph::update_normals()
{
  for (int edge=0; edge<Ne; edge++) {
    compute_normals(edge);
  }
}

void graph::update_topology()
{
  update_t();
  orient();
  update_normals();
}

/*********************************************************************
   energy helper functions
*********************************************************************/
double graph::stretch_energy(int edge)
{
  return 0.5 * epsilon * (l[edge] - l0) * (l[edge] - l0);
}

double graph::bend_energy(int edge)
{
  double bendE=0., ndot;
  int tp = etype[edge]; // type 0 is always no curvature
  if (t[edge][1] != -1 && t[edge][0] != -1) {
    if (tp==0) {
      ndot = dot(n[edge][0],n[edge][1]);
      bendE = kappa[tp] * (1 - ndot);
    }
    else {
      ndot = dot(n[edge][0],n[edge][1]);
      if (ndot>1) { ndot = 1.0; } // prevents numerical precision errors
      double theta = acos(ndot);
      bendE = kappa[tp] * (1 - cos(theta-min_bond[tp]));
    }
  }
  return bendE;
}

double graph::compute_energy()
{
  double tot_eng = 0;
  for (int i=0; i<Ne; i++)
  {
    // add the contribution from stretching the bonds
    tot_eng += stretch_energy(i);
    // compute the normal vectors for the triangles along edge i
    tot_eng += bend_energy(i);
  }
  return tot_eng;
}

double graph::vertex_energy(int v1)
{
  double v_eng = 0;
  int edge, oedge;
  for (int i=0; i<nneigh[v1]; i++)
  {
    edge = eind[v1][neigh[v1][i]];
    v_eng += stretch_energy(edge);
    v_eng += bend_energy(edge);
    // also need bending energy from neighboring triangles
    if (t[edge][0] != -1) {
      oedge = eind[neigh[v1][i]][t[edge][0]];
      v_eng += 0.5 * bend_energy(oedge);
    }
    if (t[edge][1] != -1) {
      oedge = eind[neigh[v1][i]][t[edge][1]];
      v_eng += 0.5 * bend_energy(oedge);
    }
  }
  return v_eng;
}

void make_sheet(graph &g, int Nside)
{
  int vcount = 0;
  for (int i=0; i<Nside; i++)
  {
    for (int j=0; j<Nside; j++)
    {
      g.add_vertex((j-0.5*(i%2)), sqrt(3)/2.*i, 0);
      g.nneigh[vcount]=0;
      vcount++;
    }
  }
  vcount=0;
  for (int i=0; i<Nside; i++)
  {
    for (int j=0; j<Nside; j++)
    {
      if (j<Nside-1) {
        g.add_edge(vcount,vcount+1);
      }
      if (i<Nside-1) {
        if (i%2==0)
        {
          if (j<Nside-1)
          {
            g.add_edge(vcount,vcount+Nside);
            g.add_edge(vcount,vcount+Nside+1);
          }
          else {
            g.add_edge(vcount,vcount+Nside);
          }
        }
        else
        {
          if (j>0) {
            g.add_edge(vcount,vcount+Nside-1);
            g.add_edge(vcount,vcount+Nside);
          }
          else {
            g.add_edge(vcount,vcount+Nside);
          }
        }
      }
      vcount++;
    }
  }
}

void make_hexamer(graph &g)
{
  g.add_vertex(0,0,0);
  int vorigin = g.Nv-1;
  double pi = atan(1)*4.;
  int vcount = g.Nv;
  for (int i=0; i<6; i++)
  {
    g.add_vertex(cos(2*pi/6*i),sin(2*pi/6*i),0);
    g.add_edge(vcount,vcount-1);
    g.add_edge(vcount,vorigin);
    vcount+=1;
  }
  g.add_edge(g.Nv-6,g.Nv-1);
  g.Nm += 6;
}

void make_pentamer(graph &g)
{
  g.add_vertex(0,0,0);
  double pi = atan(1)*4.;
  int vcount = 1;
  for (int i=0; i<5; i++)
  {
    g.add_vertex(cos(2*pi/5*i),sin(2*pi/5*i),-0.5);
    g.add_edge(vcount,vcount-1);
    g.add_edge(vcount,0);
    vcount+=1;
  }
  g.add_edge(1,5);
}

void make_cap(graph &g, int Nside)
{
  double s = Nside/2.0;
  double t = s*(1.0 + sqrt(5))/2.0;
  // build the vertices of the bottom cap
  g.add_vertex(t,0,s);
  g.add_vertex(0,-s,t);
  g.add_vertex(s,-t,0);
  g.add_vertex(t,0,-s);
  g.add_vertex(s,t,0);
  g.add_vertex(0,s,t);

  // populate the faces of the bottom cap
  double* dedge = new double[3];
  double* tedge = new double[3];
  double* vj = new double[3];
  double* vk = new double[3];
  for (int i=1; i<6; i++) {
    vsub(g.v[0],g.v[i],dedge);
    if (i==5) {
      vsub(g.v[5],g.v[1],tedge);
    } else {
      vsub(g.v[i],g.v[i+1],tedge);
    }
    for (int j=1; j<Nside+1; j++) {
      vj[0] = j*dedge[0]/Nside + g.v[0][0];
      vj[1] = j*dedge[1]/Nside + g.v[0][1];
      vj[2] = j*dedge[2]/Nside + g.v[0][2];
      if (j<Nside) {
        g.add_vertex(vj[0],vj[1],vj[2]);
      }
      for (int k=1; k<j; k++) {
        vk[0] = vj[0] + k*tedge[0]/Nside;
        vk[1] = vj[1] + k*tedge[1]/Nside;
        vk[2] = vj[2] + k*tedge[2]/Nside;
        g.add_vertex(vk[0],vk[1],vk[2]);
      }
    }
  }

  //delete the vertices
  for (int i=1; i<6; i++) {
    g.delete_vertex(1);
  }

  // add the edges
  for (int i=0; i<g.Nv; i++) {
    for (int j=0; j<g.Nv; j++) {
      if (dist(g.v[i],g.v[j])==1.0) {
        g.add_edge(i,j);
      }
    }
  }
  fprintf(stderr, "vcount %d and ecount %d\n", g.Nv, g.Ne);
}

void make_icosahedron(graph &g, int Nside)
{
  double s = Nside/2.0;
  double t = s*(1.0 + sqrt(5))/2.0;

  // build the vertices of the bottom cap
  g.add_vertex(t,0,s);
  g.add_vertex(0,-s,t);
  g.add_vertex(s,-t,0);
  g.add_vertex(t,0,-s);
  g.add_vertex(s,t,0);
  g.add_vertex(0,s,t);

  // build the vertices of the top cap
  g.add_vertex(-t,0,-s);
  g.add_vertex(-s,-t,0);
  g.add_vertex(0,-s,-t);
  g.add_vertex(0,s,-t);
  g.add_vertex(-s,t,0);
  g.add_vertex(-t,0,s);

  // populate the faces of the bottom cap
  double* dedge = new double[3];
  double* tedge = new double[3];
  double* vj = new double[3];
  double* vk = new double[3];
  for (int i=1; i<6; i++) {
    vsub(g.v[0],g.v[i],dedge);
    if (i==5) {
      vsub(g.v[5],g.v[1],tedge);
    } else {
      vsub(g.v[i],g.v[i+1],tedge);
    }
    for (int j=1; j<Nside+1; j++) {
      vj[0] = j*dedge[0]/Nside + g.v[0][0];
      vj[1] = j*dedge[1]/Nside + g.v[0][1];
      vj[2] = j*dedge[2]/Nside + g.v[0][2];
      g.add_vertex(vj[0],vj[1],vj[2]);
      for (int k=1; k<j; k++) {
        vk[0] = vj[0] + k*tedge[0]/Nside;
        vk[1] = vj[1] + k*tedge[1]/Nside;
        vk[2] = vj[2] + k*tedge[2]/Nside;
        g.add_vertex(vk[0],vk[1],vk[2]);
      }
    }
  }

  // bottom cap
  for (int i=7; i<12; i++) {
    vsub(g.v[6],g.v[i],dedge);
    if (i==11) {
      vsub(g.v[11],g.v[7],tedge);
    } else {
      vsub(g.v[i],g.v[i+1],tedge);
    }
    for (int j=1; j<Nside+1; j++) {
      vj[0] = j*dedge[0]/Nside + g.v[6][0];
      vj[1] = j*dedge[1]/Nside + g.v[6][1];
      vj[2] = j*dedge[2]/Nside + g.v[6][2];
      g.add_vertex(vj[0],vj[1],vj[2]);
      for (int k=1; k<j; k++) {
        vk[0] = vj[0] + k*tedge[0]/Nside;
        vk[1] = vj[1] + k*tedge[1]/Nside;
        vk[2] = vj[2] + k*tedge[2]/Nside;
        g.add_vertex(vk[0],vk[1],vk[2]);
      }
    }
  }

  // side edges
  int va, vb, vc;
  for (int i=1; i<5; i++) {
    va = i;
    vb = i+6;
    vc = i+1;
    vsub(g.v[va],g.v[vb],dedge);
    vsub(g.v[vb],g.v[vc],tedge);
    for (int j=1; j<Nside; j++) {
      vj[0] = j*dedge[0]/Nside + g.v[va][0];
      vj[1] = j*dedge[1]/Nside + g.v[va][1];
      vj[2] = j*dedge[2]/Nside + g.v[va][2];
      g.add_vertex(vj[0],vj[1],vj[2]);
      for (int k=1; k<j; k++) {
        vk[0] = vj[0] + k*tedge[0]/Nside;
        vk[1] = vj[1] + k*tedge[1]/Nside;
        vk[2] = vj[2] + k*tedge[2]/Nside;
        g.add_vertex(vk[0],vk[1],vk[2]);
      }
    }
  }

  for (int i=7; i<11; i++) {
    va = i;
    vb = i-6+1;
    vc = i+1;
    vsub(g.v[va],g.v[vb],dedge);
    vsub(g.v[vb],g.v[vc],tedge);
    for (int j=1; j<Nside; j++) {
      vj[0] = j*dedge[0]/Nside + g.v[va][0];
      vj[1] = j*dedge[1]/Nside + g.v[va][1];
      vj[2] = j*dedge[2]/Nside + g.v[va][2];
      g.add_vertex(vj[0],vj[1],vj[2]);
      for (int k=1; k<j; k++) {
        vk[0] = vj[0] + k*tedge[0]/Nside;
        vk[1] = vj[1] + k*tedge[1]/Nside;
        vk[2] = vj[2] + k*tedge[2]/Nside;
        g.add_vertex(vk[0],vk[1],vk[2]);
      }
    }
  }

  va = 5;
  vb = 11;
  vc = 1;
  vsub(g.v[va],g.v[vb],dedge);
  vsub(g.v[vb],g.v[vc],tedge);
  for (int j=1; j<Nside; j++) {
    vj[0] = j*dedge[0]/Nside + g.v[va][0];
    vj[1] = j*dedge[1]/Nside + g.v[va][1];
    vj[2] = j*dedge[2]/Nside + g.v[va][2];
    g.add_vertex(vj[0],vj[1],vj[2]);
    for (int k=1; k<j; k++) {
      vk[0] = vj[0] + k*tedge[0]/Nside;
      vk[1] = vj[1] + k*tedge[1]/Nside;
      vk[2] = vj[2] + k*tedge[2]/Nside;
      g.add_vertex(vk[0],vk[1],vk[2]);
    }
  }
  va = 11;
  vb = 1;
  vc = 7;
  vsub(g.v[va],g.v[vb],dedge);
  vsub(g.v[vb],g.v[vc],tedge);
  for (int j=1; j<Nside; j++) {
    vj[0] = j*dedge[0]/Nside + g.v[va][0];
    vj[1] = j*dedge[1]/Nside + g.v[va][1];
    vj[2] = j*dedge[2]/Nside + g.v[va][2];
    g.add_vertex(vj[0],vj[1],vj[2]);
    for (int k=1; k<j; k++) {
      vk[0] = vj[0] + k*tedge[0]/Nside;
      vk[1] = vj[1] + k*tedge[1]/Nside;
      vk[2] = vj[2] + k*tedge[2]/Nside;
      g.add_vertex(vk[0],vk[1],vk[2]);
    }
  }


  // add the edges
  for (int i=0; i<g.Nv; i++) {
    for (int j=0; j<g.Nv; j++) {
      if (dist(g.v[i],g.v[j])==1.0) {
        g.add_edge(i,j);
      }
    }
  }
}

void dump_defects(graph &g, FILE *f, int time)
{
  int N5=0, N7=0;
  for (int vi=0; vi<g.Nv; vi++) {
    if (g.surface_vertex(vi)==0 && g.nneigh[vi]!=6) {
      if (g.nneigh[vi]==5) {
        N5 += 1;
      } else if (g.nneigh[vi]==7) {
        N7 += 1;
      }
    }
  }
  fprintf(f, "%d %d %d\n", time, N5, N7);
  fflush(f);
}

void dump_xyz(graph &g, int time)
{
  char filename[80];
  sprintf(filename, "lat_%07d.xyz", time);
  FILE *f;
  f = fopen(filename, "w");
  fprintf(f,"%d\n",g.Nv);
  fprintf(f,"frame %d\n",time);
  for (int vi = 0; vi < g.Nv; vi++)
  {
    fprintf(f,"O %8.3f %8.3f %8.3f\n", g.v[vi][0], g.v[vi][1], g.v[vi][2]);
  }
  fclose(f);
}

void dump_xml_frame(graph &g, FILE* f, int time)
{
  fprintf(f,"<configuration time_step=\"%d\">\n",time);
  fprintf(f,"<position num=\"%d\">\n",g.Nvmax);
  for (int vi=0; vi<g.Nv; vi++)
  {
    fprintf(f,"%8.3f %8.3f %8.3f\n", g.v[vi][0], g.v[vi][1], g.v[vi][2]);
  }
  for (int vi=g.Nv; vi<g.Nvmax; vi++) {
    fprintf(f,"%8.3f %8.3f %8.3f\n", 0., 0., 0.);
  }
  fprintf(f,"</position>\n");
  fprintf(f,"<bond>\n");
  for (int edge=0; edge<g.Ne; edge++)
  {
    fprintf(f, "bond %d %d\n", g.e[edge][0], g.e[edge][1]);
  }
  fprintf(f,"</bond>\n");
  fprintf(f,"</configuration>\n");
}

void graph::dump_xml(int time)
{
  char filename[80];
  sprintf(filename, "lat_%07d.xml", time);
  FILE *f;
  f = fopen(filename, "w");
  fprintf(f,"<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n");
  fprintf(f,"<hoomd_xml version=\"1.4\">\n");
  fprintf(f,"<configuration time_step=\"%d\">\n",time);
  fprintf(f,"<position num=\"%d\">\n",Nvmax);
  for (int vi=0; vi<Nv; vi++)
  {
    fprintf(f,"%8.3f %8.3f %8.3f\n", v[vi][0], v[vi][1], v[vi][2]);
  }
  for (int vi=Nv; vi<Nvmax; vi++) {
    fprintf(f,"%8.3f %8.3f %8.3f\n", 0., 0., 0.);
  }
  fprintf(f,"</position>\n");
  fprintf(f,"<type num=\"%d\">\n",Nvmax);
  for (int vi=0; vi<Nvmax; vi++)
  {
    fprintf(f,"%d\n", vcs[vi]);
  }
  fprintf(f,"</type>\n");

  fprintf(f,"<bond num=\"%d\">\n",Nemax);
  for (int edge=0; edge<Ne; edge++)
  {
    fprintf(f, "bond %d %d\n", e[edge][0], e[edge][1]);
  }
  fprintf(f,"</bond>\n");
  fprintf(f,"</configuration>\n");
  fprintf(f,"</hoomd_xml>\n");
  fclose(f);
}

void dump_surf_xml(graph &g, int time)
{
  char filename[80];
  sprintf(filename, "latsurf_%07d.xml", time);
  FILE *f;
  f = fopen(filename, "w");
  fprintf(f,"<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n");
  fprintf(f,"<hoomd_xml version=\"1.4\">\n");
  fprintf(f,"<configuration time_step=\"%d\">\n",time);
  fprintf(f,"<position num=\"%d\">\n",g.Nvmax);
  for (int vi=0; vi<g.Nv; vi++)
  {
    fprintf(f,"%8.3f %8.3f %8.3f\n", g.v[vi][0], g.v[vi][1], g.v[vi][2]);
  }
  for (int vi=g.Nv; vi<g.Nvmax; vi++) {
    fprintf(f,"%8.3f %8.3f %8.3f\n", 0., 0., 0.);
  }
  fprintf(f,"</position>\n");
  fprintf(f,"<bond>\n");
  for (int edge=0; edge<g.Ne; edge++)
  {
    if (g.etype[edge]==1) {
      fprintf(f, "surfbond %d %d\n", g.e[edge][0], g.e[edge][1]);
    }
  }
  for (int edge=0; edge<g.Nperim; edge++)
  {
    int eind = g.esurf[edge];
    if (g.check_wedge(eind)==1) {
      fprintf(f, "surfbond %d %d\n", g.e[eind][0], g.e[eind][1]);
    }
  }
  fprintf(f,"</bond>\n");
  fprintf(f,"</configuration>\n");
  fprintf(f,"</hoomd_xml>\n");
  fclose(f);
}

void dump_norm_xml(graph &g, int time)
{
  char filename[80];
  sprintf(filename, "latnorm_%07d.xml", time);
  FILE *f;
  f = fopen(filename, "w");
  fprintf(f,"<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n");
  fprintf(f,"<hoomd_xml version=\"1.4\">\n");
  fprintf(f,"<configuration time_step=\"%d\">\n",time);
  fprintf(f,"<position num=\"%d\">\n",g.Nvmax);
  for (int vi=0; vi<g.Nv; vi++)
  {
    fprintf(f,"%8.3f %8.3f %8.3f\n", g.v[vi][0], g.v[vi][1], g.v[vi][2]);
  }
  int vc;
  double x,y,z;
  for (int vi=0; vi<g.Ne; vi++) {
    vc = g.e[vi][0];
    x = g.v[vc][0] + 0.5*g.evec[vi][0];
    y = g.v[vc][1] + 0.5*g.evec[vi][1];
    z = g.v[vc][2] + 0.5*g.evec[vi][2];
    fprintf(f,"%8.3f %8.3f %8.3f\n", x+g.n[vi][0][0],
            y+g.n[vi][0][1],
            z+g.n[vi][0][2]);
    fprintf(f,"%8.3f %8.3f %8.3f\n", x+g.n[vi][1][0],
            y+g.n[vi][1][1],
            z+g.n[vi][1][2]);

  }
  for (int vi=g.Nv+2*g.Ne; vi<g.Nvmax; vi++) {
    fprintf(f,"%8.3f %8.3f %8.3f\n", 0., 0., 0.);
  }
  fprintf(f,"</position>\n");
  fprintf(f,"<bond>\n");
  for (int edge=0; edge<g.Ne; edge++)
  {
    fprintf(f, "surfbond %d %d\n", g.e[edge][0], g.e[edge][1]);
    //}
  }
  fprintf(f,"</bond>\n");
  fprintf(f,"</configuration>\n");
  fprintf(f,"</hoomd_xml>\n");
  fclose(f);
}

void dump_log(graph &g) {
  FILE *f;
  f = fopen("input_parameters.out", "w");
  fprintf(f,"epsilon:\t %f\n",g.epsilon);
  fprintf(f,"kappa:\t %f\n",g.kappa[0]);
  fprintf(f,"z:\t %f\n",g.z);
  fprintf(f,"K:\t %f\n",g.K);
  fprintf(f,"xi:\t %f\n",g.xi);
  fprintf(f,"l0:\t %f\n",g.l0);
  fprintf(f,"lmax:\t %f\n",g.lmax);
  fprintf(f,"lmin:\t %f\n",g.lmin);
  fprintf(f,"lfuse:\t %f\n",g.lfuse);
  fclose(f);
}


double dist(double *v1, double *v2)
{
  return sqrt( (v2[0]-v1[0])*(v2[0]-v1[0]) +
               (v2[1]-v1[1])*(v2[1]-v1[1]) +
               (v2[2]-v1[2])*(v2[2]-v1[2]) );
}

double norm(double *v)
{
  return sqrt( v[0]*v[0] + v[1]*v[1] + v[2]*v[2] );
}

double dot(double *v1, double *v2)
{
  return v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2];
}

void cross(double *v1, double *v2, double *res)
{
  res[0] = v2[1]*v1[2] - v2[2]*v1[1];
  res[1] = v2[2]*v1[0] - v2[0]*v1[2];
  res[2] = v2[0]*v1[1] - v2[1]*v1[0];
}

void vsub(double *v1, double *v2, double *res)
{
  res[0] = v2[0] - v1[0];
  res[1] = v2[1] - v1[1];
  res[2] = v2[2] - v1[2];
}

void graph::com(int vi, int vj, int vk, double *vcom) {
  vcom[0] = (v[vi][0]+v[vj][0]+v[vk][0])/3.0;
  vcom[1] = (v[vi][1]+v[vj][1]+v[vk][1])/3.0;
  vcom[2] = (v[vi][2]+v[vj][2]+v[vk][2])/3.0;
}

void graph::nvec(int vi, int vj, int ti, double *nveci, double s) {
  double *comi = new double[3];
  int edge = eind[vi][vj];
  com(vi, vj, t[edge][ti], comi);
  nveci[0] = comi[0] - s*n[edge][ti][0];
  nveci[1] = comi[1] - s*n[edge][ti][1];
  nveci[2] = comi[2] - s*n[edge][ti][2];
  delete[] comi;
}

void graph::outernvec(int vi, int vj, int ti, double *nveci, double s) {
  double *comi = new double[3];
  int edge = eind[vi][vj];
  com(vi, vj, t[edge][ti], comi);
  nveci[0] = comi[0] + s*n[edge][ti][0];
  nveci[1] = comi[1] + s*n[edge][ti][1];
  nveci[2] = comi[2] + s*n[edge][ti][2];
  delete[] comi;
}

int graph::check_edge_dist(int ei, int ej) {
  double de0, de1;
  int ei0, ei1;
  int ej0, ej1;
  ei0 = e[ei][0];
  ei1 = e[ei][1];
  ej0 = e[ej][0];
  ej1 = e[ej][1];
  double da = dist(v[ei0],v[ej0]);
  double db = dist(v[ei0],v[ej1]);
  if (da > db) {
    de0 = db;
    de1 = dist(v[ei1],v[ej0]);
    if (de0<lfuse && de1<lfuse) {
      return 1;
    } else return 0;
  } else {
    de0 = da;
    de1 = dist(v[ei1],v[ej1]);
    if (de0<lfuse && de1<lfuse) {
      return 1;
    } else return 0;
  }
}

int graph::check_overlap(int ei) {
  int i0 = t[ei][0];
  int i1 = t[ei][1];
  double ln = 0.5*l0; // overlap distance cutoff
  double dn = 0;
  // centers of mass of the triangles
  double *comi0 = new double[3];
  double *comi1 = new double[3];


  if (i0!=-1) {
    com(e[ei][0], e[ei][1], i0, comi0);
  }
  if (i1!=-1) {
    com(e[ei][0], e[ei][1], i1, comi1);
  }
  double *comj = new double[3];
  for (int ej=0; ej<Ne; ej++) {
    if (e[ej][0]!=e[ei][0] && e[ej][1]!=e[ei][0]
        && e[ej][0]!=e[ei][1] && e[ej][1]!=e[ei][1]) {
      if (t[ej][0]!=-1) {
        com(e[ej][0], e[ej][1], t[ej][0], comj);
        if (i0!=-1) {
          dn = dist(comi0,comj);
          if (dn<ln && dn>0.) {
            return 1;
          }
        }
        if (i1!=-1) {
          dn = dist(comi1,comj);
          if (dn<ln && dn>0.) {
            return 1;
          }
        }
      }
      if (t[ej][1]!=-1) {
        com(e[ej][0], e[ej][1], t[ej][1], comj);
        if (i0!=-1) {
          dn = dist(comi0,comj);
          if (dn<ln && dn>0.) {
            return 1;
          }
        }
        if (i1!=-1) {
          dn = dist(comi1,comj);
          if (dn<ln && dn>0.) {
            return 1;
          }
        }
      }
    }
  }
  delete[] comi0; delete[] comi1; delete[] comj;
  return 0;
}
