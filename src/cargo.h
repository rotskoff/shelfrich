/* cargo.h Implements the fcc lattice gas cargo.
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

#ifndef CARGO_H
#define CARGO_H

#include <stdio.h>
#include <cstdlib>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <vector>
#include <set>
#include "graph.h"

using namespace std;

/** Ising lattice gas class for cargo. */
class cargo {
public:
/** Lattice parameters.
 *  Set the size of the system. Declare the arrays.
 */
int Nsurf;   /**< Current number of surface lattice sites. */
int Nsurfocc;   /**< Current number of surface lattice sites. */
int N;   /**< Current number of filled lattice sites. */
int Nmax;   /**< Maximum extent in any direction Nmax/2. */
vector <int> lat;   /*** full state of the lattice  */
vector <double> occ;   /*** graph occupancy state of the lattice */
set <int> sites;   /*** array of occupied sites */
set <int> surf;   /*** array of site indices of the occupied surface sites.*/
set <int> occsurf;   /*** array of site indices of the occupied surface sites.*/

/** Cargo lattice gas parameters.
 *  Set the energetic parameters of the cargo and cargo - graph interaction.
 */
double epsilon;   /**< Cargo - cargo coupling strength. \f$\epsilon = 4J\f$ */
double mu;   /**< Cargo chemical potential. \f$\mu = 2h - 4Jz\f$ */
double gamma;   /**< Cargo - graph coupling strength. */
double gammao;
double T;   /**< The temperature of the cargo lattice gas. */
double energy;   /**< Current value of the total energy. */

bool umbrella_sampling; /**< Enable umbrella sampling. */
int Ncenter; /**< The center of the harmonic potential for umbrella sampling */
double kumb; /**<  The force constant for the umbrella potential */

/** Default constructor. */
cargo();
/** Default destructor. */
~cargo();
/** Parameter filling constructor. */
cargo(int Nmax, double epsilon1, double mu1, double gamma1, double gammao1, double T1);

/** Set all lattice sites to 0. */
void empty();
/** Get the index of site i,j,k */
int site_index(int i, int j, int k);
/** Get the index from coordinates*/
int site_from_coord(double x, double y, double z);

/** Compute the number of occupied nearest neighbors.
 * @param i indices of the lattice site.
 * @return The number of occupied neighbors. */
int nneigh(int site);

/** Check if the site has no occupied neighbors */
int check_isolated_neigh(int site);

/** Compute the total energy of the system.
    @return A double with the current value of the energy. */
double compute_energy();

/** Compute the energy difference resulting from a virtual spin flip.
 * @param i indices of the lattice site.
 * @return The energy difference resulting from the flip. */
double compute_de(int site);

/** Perform a single site Monte Carlo spin flip using Metropolis dynamics.
 * @param site An int giving the site index.
 * @param r A pointer to a random number generator. */
bool mc_step(int site, gsl_rng *r);

/** Attempt Nmax single site Monte Carlo spin flips using Metropolis dynamics.
 * @param r A pointer to a random number generator. */
void mc_step_controller(gsl_rng *r);

/** Update the surface sites of the cargo globule. */
void update_surf();

/** Compute the interaction between the triangle defined by vi, vj, and
    the vertex t[eij][ti].
 * @param vi Edge vertex.
 * @param vj Edge vertex.
 * @param ti Index of the triangle.
 * @param g  A pointer to a graph object. */
double cargo_graph_energy(int vi, int vj, int ti, graph &g);

/** Compute the interaction energy between all triangles associated with
    vertex vi and the cargo. */
double cargo_graph_energy_vertex(int vi, graph &g);

/** Compute the set of occupied lattice sites. */
void compute_occupied(graph &g);

/** Set the occ_val array. */
void set_occ(int site, double occ_val);

/** Dump a HOOMD style xml file with the graph and cargo. */
void dump_xml(graph &g, int time);

/** Dump the occupied cargo sites. */
void dump_occ(graph &g, int time);

/** Dump the surface cargo sites. */
void dump_surf(graph &g, int time);
};

#endif
