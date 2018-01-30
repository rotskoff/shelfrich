/* cargo.cc Implements the fcc lattice gas cargo.
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

#include <math.h>
#include "cargo.h"
#include "graph.h"

cargo::cargo()
{
  Nsurf = 0;
  N = 0;
  Nmax = 0;
  epsilon = 0.0;
  mu = 0.0;
  gamma = 0.0;
  gammao = 0.0;
  T = 0.0;
  energy = 0.0;
  kumb = 0;
  Ncenter = 0;
  umbrella_sampling=false;
}

cargo::cargo(int Nmax1, double epsilon1,
             double mu1, double gamma1, double gammao1, double T1)
{
  Nmax = Nmax1;
  epsilon = epsilon1;
  mu = mu1;
  gamma = gamma1;
  gammao = gammao1;
  T = T1;
}

cargo::~cargo(){
}

void cargo::empty()
{
  lat.resize(Nmax*Nmax*Nmax,0);
  occ.resize(Nmax*Nmax*Nmax,0);
  sites.clear();
  surf.clear();
  occsurf.clear();
  N = 0;
  Nsurf = 0;
}

int cargo::site_index(int i, int j, int k)
{
  return i + j * Nmax + k * Nmax * Nmax;
}


int cargo::site_from_coord(double x, double y, double z)
{
  int i = round(2*x) + Nmax/2;
  int j = round(2*y) + Nmax/2;
  int k = round(2*z) + Nmax/2;
  if (((i+j+k)%2)==1) {
    double di = i/2.0 - Nmax/4. - x;
    double dj = j/2.0 - Nmax/4. - y;
    double dk = k/2.0 - Nmax/4. - z;
    if ((fabs(di) > fabs(dj)) && (fabs(di) > fabs(dk))) {
      if (di<0) {
        i-=1;
      } else i+=1;
      return i + j * Nmax + k * Nmax * Nmax;
    } else if ((fabs(dj) > fabs(di)) && (fabs(dj) > fabs(dk))) {
      if (dj<0) {
        j-=1;
      } else j+=1;
      return i + j * Nmax + k * Nmax * Nmax;
    } else if ((fabs(dk) > fabs(di)) && (fabs(dk) > fabs(dj))) {
      if (dk<0) {
        k-=1;
      } else k+=1;
      return i + j * Nmax + k * Nmax * Nmax;
    }
  }
  return i + j * Nmax + k * Nmax * Nmax;
}

int cargo::nneigh(int site)
{
  int i,j,k;
  k = site / (Nmax*Nmax);
  j = (site - k*Nmax*Nmax) / Nmax;
  i = site - j*Nmax - k*Nmax*Nmax;
  int nn=0;
  /* FCC */
  nn += lat[site_index(i-1,j-1,k)];
  nn += lat[site_index(i-1,j+1,k)];
  nn += lat[site_index(i+1,j-1,k)];
  nn += lat[site_index(i+1,j+1,k)];
  // j+1 neighs
  nn += lat[site_index(i,j-1,k+1)];
  nn += lat[site_index(i-1,j,k+1)];
  nn += lat[site_index(i+1,j,k+1)];
  nn += lat[site_index(i,j+1,k+1)];
  // j-1 neighs
  nn += lat[site_index(i,j-1,k-1)];
  nn += lat[site_index(i-1,j,k-1)];
  nn += lat[site_index(i+1,j,k-1)];
  nn += lat[site_index(i,j+1,k-1)];
  return nn;
}

int cargo::check_isolated_neigh(int site)
{
  int i,j,k;
  k = site / (Nmax*Nmax);
  j = (site - k*Nmax*Nmax) / Nmax;
  i = site - j*Nmax - k*Nmax*Nmax;
  if (nneigh(lat[site_index(i-1,j-1,k)])==0) {
    return 1;
  }
  if (nneigh(lat[site_index(i-1,j+1,k)])==0) {
    return 1;
  }
  if (nneigh(lat[site_index(i+1,j-1,k)])==0) {
    return 1;
  }
  if (nneigh(lat[site_index(i+1,j+1,k)])==0) {
    return 1;
  }
  if (nneigh(lat[site_index(i,j-1,k+1)])==0) {
    return 1;
  }
  if (nneigh(lat[site_index(i-1,j,k+1)])==0) {
    return 1;
  }
  if (nneigh(lat[site_index(i+1,j,k+1)])==0) {
    return 1;
  }
  if (nneigh(lat[site_index(i,j+1,k+1)])==0) {
    return 1;
  }
  if (nneigh(lat[site_index(i,j-1,k-1)])==0) {
    return 1;
  }
  if (nneigh(lat[site_index(i-1,j,k-1)])==0) {
    return 1;
  }
  if (nneigh(lat[site_index(i+1,j,k-1)])==0) {
    return 1;
  }
  if (nneigh(lat[site_index(i,j+1,k-1)])==0) {
    return 1;
  }
  return 0;
}

double cargo::compute_energy()
{
  int isite;
  double E=0; // the energy
  // Loop over sites
  set<int>::iterator it;
  for (it=sites.begin(); it!=sites.end(); it++)
  {
    isite = *it;
    E -= mu; // chemical potential contribution
    E += occ[isite];
    E -= 0.5 * epsilon * nneigh(isite);
  }
  energy = E;
  return E;
}

double cargo::compute_de(int site)
{
  int n;
  if (lat[site]==1) {
    n = -1;
  } else n = 1;
  double de = -mu * n;
  de +=  n * occ[site]; // occ is positive is outnvec inside site
  de += -n * epsilon * nneigh(site);
  if (umbrella_sampling==true) {
    de += (0.5*kumb*(N+n-Ncenter)*(N+n-Ncenter))-(0.5*kumb*(N-Ncenter)*(N-Ncenter));
  }
  return de;
}

bool cargo::mc_step(int site, gsl_rng *r)
{
  double de=0.;
  // do not delete the initial seed
  if (lat[site]==1 && N==1) {
    update_surf();
    return 0;
  }
  de = compute_de(site);
  double rand = gsl_rng_uniform(r);
  double Nsurf_init = (double)Nsurf;
  lat[site] = (lat[site]+1) % 2;
  if (lat[site]==1) {
    sites.insert(site);
  } else {
    sites.erase(site);
  }
  update_surf();
  double gen = ((double)Nsurf)/Nsurf_init;
  if (rand < gen*exp(-de/T)) {
    N += 2*lat[site] - 1;
    energy += de;
    return 1;
  }
  else {
    lat[site] = (lat[site]+1) % 2;
    if (lat[site]==1) {
      sites.insert(site);
    } else {
      sites.erase(site);
    }
    update_surf();
    return 0;
  }
}

void cargo::mc_step_controller(gsl_rng *r)
{
  int isite, site;
  isite = gsl_rng_uniform_int(r, Nsurf);
  set<int>::iterator it;
  it = surf.begin();
  advance(it,isite);
  site = *it;
  mc_step(site,r);
}

void cargo::update_surf()
{
  N = 0;
  Nsurf = 0;
  surf.clear();
  Nsurfocc = 0;
  occsurf.clear();
  int isite;
  set<int>::iterator it;
  for (it=sites.begin(); it!=sites.end(); it++)
  {
    isite = *it;
    N += 1;
    int nn = nneigh(isite);
    if (nn==0) {
      lat[isite]=0;
      sites.erase(isite);
      N-=1;
    }
  }
  for (it=sites.begin(); it!=sites.end(); it++) {
    // loop over neighbors and add them all
    isite = *it;
    int k = isite / (Nmax*Nmax);
    int j = (isite - k*Nmax*Nmax) / Nmax;
    int i = isite - j*Nmax - k*Nmax*Nmax;
    if ((i+j+k)%2==1) {
      fprintf(stderr,"vacancy occ\n");
      exit(-1);
    }
    int nsite;
    int nn = nneigh(isite);
    pair<set<int>::iterator, bool> success;
    /* FCC */
    if (nn<12) {
      success = surf.insert(isite);
      if (success.second) {
        Nsurf+=1;
      }
      // occupied sites
      success = occsurf.insert(isite);
      if (success.second) {
        Nsurfocc += 1;
      }
    }
    nsite = site_index(i-1,j-1,k);
    if (nneigh(nsite)<12) {
      success = surf.insert(nsite);
      if (success.second) {
        Nsurf+=1;
      }
    }
    nsite = site_index(i-1,j+1,k);
    if (nneigh(nsite)<12) {
      success = surf.insert(nsite);
      if (success.second) {
        Nsurf+=1;
      }
    }
    nsite = site_index(i+1,j-1,k);
    if (nneigh(nsite)<12) {
      success = surf.insert(nsite);
      if (success.second) {
        Nsurf+=1;
      }
    }
    nsite = site_index(i+1,j+1,k);
    if (nneigh(nsite)<12) {
      success = surf.insert(nsite);
      if (success.second) {
        Nsurf+=1;
      }
    }
    nsite = site_index(i,j-1,k+1);
    if (nneigh(nsite)<12) {
      success = surf.insert(nsite);
      if (success.second) {
        Nsurf+=1;
      }
    }
    nsite = site_index(i-1,j,k+1);
    if (nneigh(nsite)<12) {
      success = surf.insert(nsite);
      if (success.second) {
        Nsurf+=1;
      }
    }
    nsite = site_index(i+1,j,k+1);
    if (nneigh(nsite)<12) {
      success = surf.insert(nsite);
      if (success.second) {
        Nsurf+=1;
      }
    }
    nsite = site_index(i,j+1,k+1);
    if (nneigh(nsite)<12) {
      success = surf.insert(nsite);
      if (success.second) {
        Nsurf+=1;
      }
    }

    nsite = site_index(i,j-1,k-1);
    if (nneigh(nsite)<12) {
      success = surf.insert(nsite);
      if (success.second) {
        Nsurf+=1;
      }
    }
    nsite = site_index(i,j+1,k-1);
    if (nneigh(nsite)<12) {
      success = surf.insert(nsite);
      if (success.second) {
        Nsurf+=1;
      }
    }
    nsite = site_index(i-1,j,k-1);
    if (nneigh(nsite)<12) {
      success = surf.insert(nsite);
      if (success.second) {
        Nsurf+=1;
      }
    }
    nsite = site_index(i+1,j,k-1);
    if (nneigh(nsite)<12) {
      success = surf.insert(nsite);
      if (success.second) {
        Nsurf+=1;
      }
    }
  }
}

double cargo::cargo_graph_energy(int vi, int vj, int ti, graph &g)
{
  if (g.surface_edge(g.eind[vi][vj])==1 || g.check_neigh(vi,vj)==0) {
    return 0;
  }
  double *vx = new double[3];
  double *nx = new double[3];
  double evi = 0;
  g.outernvec(vi,vj,ti,nx,0.5);
  int osite = site_from_coord(nx[0],nx[1],nx[2]);
  evi += gammao*(lat[osite]);
  g.nvec(vi,vj,ti,nx,0.25);
  int isite = site_from_coord(nx[0],nx[1],nx[2]);
  evi += -gamma*(lat[isite]);
  delete[] nx;
  delete[] vx;
  return evi;
}

double cargo::cargo_graph_energy_vertex(int vi, graph &g)
{
  double evi = 0;
  for (int ni=0; ni<g.nneigh[ni]; ni++) {
    int vj = g.neigh[vi][ni];
    int ei = g.eind[vi][vj];
    evi += cargo_graph_energy(vi,vj,0,g);
    evi += cargo_graph_energy(vi,vj,1,g);
    if (g.t[ei][0]!=-1) {
      evi += 0.5*cargo_graph_energy(vj,g.t[ei][0],0,g);
    }
    if (g.t[ei][1]!=-1) {
      evi += 0.5*cargo_graph_energy(vj,g.t[ei][1],1,g);
    }
  }
  return evi;
}

void cargo::compute_occupied(graph &g)
{
  fill(occ.begin(), occ.end(), 0);
  double *nx = new double[3];
  int vi, vj;
  for (int ei=0; ei<g.Ne; ei++) {
    vi = g.e[ei][0];
    vj = g.e[ei][1];
    if (g.surface_edge(ei)==0) {
      // get the overcounting scalar
      double sfac = 1.;
      int vk = g.t[ei][0];
      if (g.surface_edge(g.eind[vi][vk])==0 &&
          g.surface_edge(g.eind[vj][vk])==0) {
        sfac = 1./3.;
      }
      else if (g.surface_edge(g.eind[vi][vk])==0 ||
               g.surface_edge(g.eind[vj][vk])==0) {
        sfac = 1./2.;
      }

      g.outernvec(vi,vj,0,nx,0.5);
      set_occ(site_from_coord(nx[0],nx[1],nx[2]), sfac*gammao);
      g.nvec(vi,vj,0,nx,0.25);
      set_occ(site_from_coord(nx[0],nx[1],nx[2]), -sfac*gamma);
      // next one
      sfac = 1.;
      vk = g.t[ei][1];
      if (g.surface_edge(g.eind[vi][vk])==0 &&
          g.surface_edge(g.eind[vj][vk])==0) {
        sfac = 1./3.;
      }
      else if (g.surface_edge(g.eind[vi][vk])==0 ||
               g.surface_edge(g.eind[vj][vk])==0) {
        sfac = 1./2.;
      }
      g.outernvec(vi,vj,1,nx,0.5);
      set_occ(site_from_coord(nx[0],nx[1],nx[2]), sfac*gammao);
      g.nvec(vi,vj,1,nx,0.25);
      set_occ(site_from_coord(nx[0],nx[1],nx[2]), -sfac*gamma);
    }
  }
  delete[] nx;
}

void cargo::set_occ(int site, double occ_val)
{
  occ[site] += occ_val;
}

void cargo::dump_xml(graph &g, int time)
{
  char filename[80];
  sprintf(filename, "lat_%07d.xml", time);
  FILE *f;
  f = fopen(filename, "w");
  fprintf(f,"<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n");
  fprintf(f,"<hoomd_xml version=\"1.4\">\n");
  fprintf(f,"<configuration time_step=\"%d\">\n",time);
  fprintf(f,"<position num=\"%ld\">\n",g.Nvmax+sites.size());
  for (int vi=0; vi<g.Nv; vi++)
  {
    fprintf(f,"%8.3f %8.3f %8.3f\n", g.v[vi][0], g.v[vi][1], g.v[vi][2]);
  }
  for (int vi=g.Nv; vi<g.Nvmax; vi++) {
    fprintf(f,"%8.3f %8.3f %8.3f\n", 0., 0., 0.);
  }
  set<int>::iterator it;
  int site;
  for (it=sites.begin(); it!=sites.end(); it++)
  {
    site = *it;
    int i,j,k; double x,y,z;
    k = site / (Nmax*Nmax);
    j = (site - k*Nmax*Nmax) / Nmax;
    i = site - j*Nmax - k*Nmax*Nmax;
    /* FCC */
    x = 0.5*i-Nmax/4.0;
    y = 0.5*j-Nmax/4.0;
    z = 0.5*k-Nmax/4.0;
    fprintf(f,"%8.3f %8.3f %8.3f\n",x,y,z);
  }
  fprintf(f,"</position>\n");
  fprintf(f,"<bond num=\"%d\">\n",g.Nemax);
  for (int edge=0; edge<g.Ne; edge++)
  {
    fprintf(f, "bond %d %d\n", g.e[edge][0], g.e[edge][1]);
  }
  fprintf(f,"</bond>\n");
  fprintf(f,"</configuration>\n");
  fprintf(f,"</hoomd_xml>\n");
  fclose(f);
}

void cargo::dump_occ(graph &g, int time)
{
  char filename[80];
  sprintf(filename, "latocc_%07d.xml", time);
  FILE *f;
  f = fopen(filename, "w");
  fprintf(f,"<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n");
  fprintf(f,"<hoomd_xml version=\"1.4\">\n");
  fprintf(f,"<configuration time_step=\"%d\">\n",time);
  fprintf(f,"<position num=\"%ld\">\n",g.Nvmax+sites.size());
  //fprintf(f,"<position num=\"%d\">\n",g.Nvmax+N+Nmax*Nmax*Nmax);
  for (int vi=0; vi<g.Nv; vi++)
  {
    fprintf(f,"%8.3f %8.3f %8.3f\n", g.v[vi][0], g.v[vi][1], g.v[vi][2]);
  }
  int ecount = 0;
  for (int vi=g.Nv+ecount; vi<g.Nvmax; vi++) {
    fprintf(f,"%8.3f %8.3f %8.3f\n", 0., 0., 0.);
  }
  set<int>::iterator it;
  int site;
  for (it=sites.begin(); it!=sites.end(); it++)
  {
    site = *it;
    int i,j,k; double x,y,z;
    k = site / (Nmax*Nmax);
    j = (site - k*Nmax*Nmax) / Nmax;
    i = site - j*Nmax - k*Nmax*Nmax;
    z = k-Nmax/2.0 + 0.5*((j+1)%2);
    y = 0.5*j-Nmax/4.0;
    x = i-Nmax/2.0 + 0.5*((j+1)%2);
    if (occ[site]!=0) {
      fprintf(f,"%8.3f %8.3f %8.3f\n",x,y,z);
    }
    else {
      fprintf(f,"%8.3f %8.3f %8.3f\n",0.,0.,0.);
    }
  }
  fprintf(f,"</position>\n");
  fprintf(f,"<bond num=\"%d\">\n",g.Nemax);
  for (int edge=0; edge<g.Ne; edge++)
  {
    fprintf(f, "bond %d %d\n", g.e[edge][0], g.e[edge][1]);
  }
  fprintf(f,"</bond>\n");
  fprintf(f,"</configuration>\n");
  fprintf(f,"</hoomd_xml>\n");
  fclose(f);
}

void cargo::dump_surf(graph &g, int time)
{
  char filename[80];
  sprintf(filename, "latsurf_%07d.xml", time);
  FILE *f;
  f = fopen(filename, "w");
  fprintf(f,"<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n");
  fprintf(f,"<hoomd_xml version=\"1.4\">\n");
  fprintf(f,"<configuration time_step=\"%d\">\n",time);
  int count = 0;
  set<int>::iterator it;
  int site;
  for (it=surf.begin(); it!=surf.end(); it++) {
    site = *it;
    if (lat[site]==1) {
      count++;
    }
  }
  fprintf(f,"<position num=\"%d\">\n",g.Nvmax+count);
  for (int vi=0; vi<g.Nv; vi++)
  {
    fprintf(f,"%8.3f %8.3f %8.3f\n", g.v[vi][0], g.v[vi][1], g.v[vi][2]);
  }
  for (int vi=g.Nv; vi<g.Nvmax; vi++) {
    fprintf(f,"%8.3f %8.3f %8.3f\n", 0., 0., 0.);
  }
  for (it=surf.begin(); it!=surf.end(); it++)
  {
    site = *it;
    int i,j,k; double x,y,z;
    if (lat[site]==1) {
      k = site / (Nmax*Nmax);
      j = (site - k*Nmax*Nmax) / Nmax;
      i = site - j*Nmax - k*Nmax*Nmax;
      /*
         z = k-Nmax/2.0 + 0.5*((j+1)%2);
         y = 0.5*j-Nmax/4.0;
         x = i-Nmax/2.0 + 0.5*((j+1)%2);
       */
      x = i-Nmax/2+0.5;
      y = j-Nmax/2+0.5;
      z = k-Nmax/2+0.5;
      fprintf(f,"%8.3f %8.3f %8.3f\n",x,y,z);
    }
  }
  fprintf(f,"</position>\n");
  fprintf(f,"<bond num=\"%d\">\n",g.Nemax);
  for (int edge=0; edge<g.Ne; edge++)
  {
    fprintf(f, "bond %d %d\n", g.e[edge][0], g.e[edge][1]);
  }
  fprintf(f,"</bond>\n");
  fprintf(f,"</configuration>\n");
  fprintf(f,"</hoomd_xml>\n");
  fclose(f);
}
