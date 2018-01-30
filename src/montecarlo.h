/* montecarlo.h Implements the Metropolis Monte Carlo moves for assembly.
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


#ifndef MONTECARLO_H
#define MONTECARLO_H

#include "cargo.h"
#include "graph.h"
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

/** Vertex move for use without the cargo structure.
 * @param g The graph structure.
 * @param override Ignore the length constraints, e.g., for equilibration.
 * @param r Point to a random number generator. */
void vertex_move(graph &g, bool override, gsl_rng *r);

/** Vertex move incorporating cargo-graph interaction and fusion attempts.
 * @param g The graph structure.
 * @param c The cargo structure.
 * @param d The magnitude of the displacement.
 * @param r Point to a random number generator. */
void vertex_move_cargo(graph &g, cargo &c, gsl_rng *r);

/** Attempt a fusion move.
 * @param g, Pointer to the graph structure.
 * @param c, Pointer to the cargo structure.
 * @param vi, Surface vertex.
 * @param vj, Surface vertex.
 * @param r,  RNG
 * @return 1 if move is accepted, 0 otherwise */
int attempt_fusion(graph &g, cargo &c, gsl_rng *r);

/** Attempt a fusion move.
 * @param g, Pointer to the graph structure.
 * @param c, Pointer to the cargo structure.
 * @param vi, Surface vertex.
 * @param vj, Surface vertex.
 * @param r,  RNG
 * @return 1 if move is accepted, 0 otherwise */
int attempt_edge_fusion(graph &g, cargo &c, int vi, int vj, gsl_rng *r);

/** Attempt a fission move.
 * @param g, Pointer to the graph structure.
 * @param c, Pointer to the cargo structure.
 * @param vi, Vertex.
 * @param r,  RNG
 * @return 1 if move is accepted, 0 otherwise */
int attempt_fission(graph &g, cargo &c, gsl_rng *r);

/** Attempt a fusion move.
 * @param g, Pointer to the graph structure.
 * @param c, Pointer to the cargo structure.
 * @param vi, Surface vertex.
 * @param vj, Surface vertex.
 * @param r,  RNG
 * @return 1 if move is accepted, 0 otherwise */
int attempt_edge_fission(graph &g, cargo &c, int vi, int vj, gsl_rng *r);

/** Attempt insertion on a wedge edge.
 * @param g The graph structure.
 * @param c The cargo structure.
 * @param edge The attempt edge.
 * @param r Pointer to a random number generator. */
int insert_monomer_wedge(graph &g, cargo &c, int edge, gsl_rng *r);

/** Attempt insertion on a free edge.
 * @param g The graph structure.
 * @param c The cargo structure.
 * @param edge The attempt edge.
 * @param r Pointer to a random number generator. */
int insert_monomer_free(graph &g, cargo &c, int edge, gsl_rng *r);

/** Attempt monomer deletion on a surface edge.
 * @param g The graph structure.
 * @param c The cargo structure.
 * @param edge The attempt edge.
 * @param r Pointer to a random number generator. */
int delete_monomer(graph &g, cargo &c, int edge, gsl_rng *r);

/** Generate idealized coordinates for a new monomer insertion.
 * @param g The graph structure.
 * @param vnew The new vertex.
 * @param edge The attempt edge. */
void new_center(graph &g, double *vnew, int edge);

#endif
