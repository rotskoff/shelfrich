/* graph.h Implements elastic shell structure.
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



#ifndef GRAPH_H
#define GRAPH_H
#include <stdio.h>
#include <set>
#include <gsl/gsl_rng.h>


using namespace std;
/** The growing triangulated network with a discretized elastic Hamiltonian */
struct graph
{
  int Nvmax; /**< Maximum number of vertices. */
  int Nemax; /**< Maximum number of edges. */
  int Nv; /**< Current number of vertices. */
  int Ne; /**< Current number of edges. */
  int Nm; /**< Current number of monomers. */
  int Nperim; /**< Current number of surface edges */
  int Nfusion; /**< Number of fusion eligible vertices */
  int Nfission; /**< Number of fission eligible vertices */
  double vb; /**< volume of random perturbations added to the monomers */
  int Ntype; /**< Number of bond types. */
  double *min_bond; /**< the minimum of the quadratic bending potential. */
  double **v; /**< Current positions of the vertices. */
  int *vcs; /**< Vertex-indexed array of the center vertices. 1 if center. */
  int **e; /**< Pairs i,j of vertex indices that are connected. */
  int *etype; /**< Edge-indexed list of edge types. */
  int **eind; /**< Returns the edge index of for the vertices i,j. */
  int **t; /**< The triangles that share an edge e. Index by 0 and 1. */
  int ***tri; /**< The oriented triangles for the computation of normals. */
  double ***n; /**< The normal vectors to triangles sharing edge e. */
  double **evec; /**< The edge vector indexed by edges. */
  int **neigh; /**< Vertex-indexed array of neighbors of vertex i. */
  int *nneigh; /**< Vertex-indexed array of number of neighbors of vertex i. */
  int *esurf; /**< Array of edge indices of surface edges. */
  double *l; /**< Edge-indexed array of bond lengths of the edges in e. */

  set< pair <int,int> >vfusion; /** Set of fusable vertices */
  set< pair <int,int> >vfission; /** Set of fissionable vertices */
  int **fusion_neigh; /**< Vertex-indexed array of neighbors of vertex i. */
  int *fusion_nneigh; /**< Vertex-indexed array of number of neighbors of vertex i. */

  double energy; /**< Current total energy of the graph. */
  double z; /**< dimensional generalized insertion affinity. */
  double K; /**< Volume scaled chemical potential contribution for monomers. */
  double xi;/**< Typical bond length fluctuation. Natural length scale. */
  double epsilon; /**< Elastic spring constant for stretching energy. */
  double l0; /**< equilibrium bond length. */
  double *kappa; /**< Spring constant for bending fluctuations. */
  double lmax; /**< Hard constraint on max edge length. */
  double lmin; /**< Hard constraint on min edge length. */
  double lfuse; /**< The distance at which fusion moves are proposed. */
  double vfuse; /**< The distance at which fusion moves are proposed. */
  double T; /**< the temperature of the lattice. */

  bool umbrella_sampling; /**< turn on umbrella sampling */
  int Ncenter; /**< number of monomers at the center of the harmonic well */
  double kumb; /**< umbrella sampling spring constant*/

  /** Default constructor */
  graph();

  /** Default destructor */
  ~graph();

  /* functions of the graph */
  /** Initialize the arrays for the system.
   * @param Nvmax The maximum number of vertices.
   * @param Nemax The maximum number of edges.
   * @param Ntype The number of edge types. */
  void initialize(int Nvmax, int Nemax, int Ntype);

  /** Center the 0-index vertex at the origin. */
  void recenter();

  /** Center the 0-index vertex at the specified point. */
  void center_at(double x, double y, double z);

  /** Create a new vertex at the specified point. */
  void add_vertex(double x, double y, double z);

  /** Delete the vertex indexed by vi.
   * @param vi The vertex to be deleted. */
  void delete_vertex(int vi);

  /** Add an edge between vertices vi and vj */
  int add_edge(int vi, int vj);

  /** Add an edge between vertices vi and vj ignoring constraints */
  int force_add_edge(int vi, int vj);

  /** Add the edge indexed by edge */
  void delete_edge(int edge);

  /** Update all edges connected to the vertex vi. */
  int update_edges(int vi);

  /** Check if the vertex is on the surface. */
  int surface_vertex(int vi);

  /** Check if the edge is on the surface. */
  int surface_edge(int edge);

  /** Check if two vertices are connected by an edge.
   * @param vi, Vertex.
   * @param vj, Vertex.
   * @return 1 if neighbors, 0 otherwise */
  int check_neigh(int vi, int vj);

  /** Check if vj is a neighbor of a neighbor of vi */
  int check_next_neigh(int vi, int vj);

  /** Check if two vertices are part of a single triangle.
   * @param vi, Vertex.
   * @param vj, Vertex.
   * @return 1 if part of a triangle, 0 otherwise */
  int check_adjacent(int vi, int vj);

  /** Check if a vertex is 5, 6, or 7 coordinated */
  int check_defect(int vi);

  /** Check if a surface vertex pair is fusable. */
  int check_fusion_vertex_pair(int vi, int vj);

  /** Update the set of fusable vertex pairs. */
  void add_fusion_pair(int vi, int vj);

  /** Update the set of fusable vertices and fusion_neigh arrays */
  void update_fusion_set();

  /** Fusion helper function, copy the coordinates of vi, its neighs
   * @param vj, Vertex.
   * @param vi_init, Vector pointer for coordinates
   * @param vi_init_neigh, Vector pointer for neighbors
   * @param vi_init_nneigh, Int pointer for number of neighbors. */
  void save_local_state(int vi, vector<double> &vi_init,
                        vector<int> &vi_neigh_init, int &vi_nneigh_init);

  /** Fusion helper function, adjust neighbor counts for vertex deletion
   * @param vi, Base vertex (neighs of vi)
   * @param vj, Additional vertex to be deleted
   * @param vi_init_neigh, Vector pointer for neighbors. */
  void adjust_indices(int vi, int vj, int vi_nneigh_init, vector<int> &vi_neigh_init);

  /** Fusion helper function, adjust neighbor counts for vertex deletion
   * @param vi, Base vertex (neighs of vi)
   * @param vj, Additional vertex to be deleted
   * @param vi1, Additional vertex to be deleted
   * @param vj1, Additional vertex to be deleted
   * @param vi_init_neigh, Vector pointer for neighbors. */
  void adjust_indices_edge_fusion(int vi, int vj, int vi1, int vj1,
                                  vector<int> &vi_neigh_init);

  /** Fusion helper function, restore vi to the graph
   * @param vi_init, Vector pointer for coordinates
   * @param vi_init_neigh, Vector pointer for neighbors
   * @param vi_init_nneigh, Int pointer for number of neighbors. */
  void restore_local_state(vector<double> &vi_init, vector<int> &vi_neigh_init,
                           int &vi_nneigh_init);

  /** Fusion helper function, add point at midpoint */
  void add_midpoint_vertex(int vi, int vj);

  /** Fusion helper function, check new vertex is allowed */
  int check_distance_constraints(int vi);

  /** Fusion helper function, add edges of vi to Nv-1
   * @param vi, Vertex
   * @param n_added, Pointer to the number of additions */
  int add_new_edges(int vi, int vtarget, int &n_added);

  /** Check if a vertex has any edges eligible for fission.
   * @param vj, Vertex.
   * @return 1 if fissionable, 0 otherwise */
  int check_fission_vertex(int vi);

  /** Check if a vertex has any edges eligible for fission.
   * @param vi, Surface vertex.
   * @param vj, Vertex.
   * @return 1 if fissionable pair, 0 otherwise */
  int check_fission_edge(int vi, int vj);

  /** Update the set of fissionable vertices */
  void update_fission_set();

  /** Select a vertex as a fission partner
   * @param vi, Fission vertex.
   * @param r, RNG  pointer.
   * @return vj vertex */
  int select_random_fission_partner(int vi, gsl_rng *r);

  /** Fission helper function, connect neighbors after splitting
   * @param vc, shared neighbor of vi, vj
   * @param vi, Fission vertex.
   * @param vj, Fission partner.
   * @param vnew, Added vertex to connect edges to.
   * @return 1 if successful edge additions */
  int connect_neighs(int vc, int vi, int vj, int vnew);

  /** Edge fission helper function, check if vertices are topologically
      connected after an edge fission move. This is basically just a
      breadth first search across the vertices.
   * @param vi, Vertex
   * @param vj, Vertex
   * @return 1 if topologically connected, 0 otherwise */
  int check_topologically_connected(int vi, int vj);

  /** Edge fission helper function, check if vertices are topologically
      connected after an edge fission move. This is basically just a
      breadth first search across the vertices. Block vertices vk, vl.
   * @param vi, Vertex
   * @param vj, Vertex
   * @param vk, Vertex excluded from search.
   * @param vl, Vertex excluded from search.
   * @return 1 if topologically connected, 0 otherwise */
  int check_topologically_connected_blocking(int vi, int vj, int vk, int vl);

  /** Fission helper function, adds random pair of vertices
   * @param vi, vertex. */
  void add_fission_vertices(int vi, gsl_rng *r);

  void save_local_state_fission(int vj, int vi, vector<double> &vi_init,
                                vector<int> &vi_neigh_init, int &vi_nneigh_init);

  void restore_local_state_fission(int vj, vector<double> &vi_init,
                                   vector<int> &vi_neigh_init, int &vi_nneigh_init);

  /** Fission helper function, remove edges without topology update
   * @param vi, vertex. */
  void remove_edges(int vi);

  /** Fission helper function, remove n latest vertices
   * @param n_times, number of vertices to remove. */
  void remove_newest_vertex(int n_times);

  /** Fusion helper function, remove the most recently added edges and last vi
   * @param n_added, number of edges to remove. */
  void remove_added(int n_added);

  /** Update the triangle and surface arrays. */
  void update_t();

  /** Compute the normal vector defined by a 3 vertex triangle.
   * @param vi Vertex.
   * @param vj Vertex.
   * @param vk Vertex.
   * @param n Pointer to the normal vector. Modified in place. */
  void get_normal(int v1, int v2, int v3, double *n);

  /** Compute the normal vectors associated with vertex vi. */
  void compute_normals(int vi);

  /**  Set the orientation as a cycle (v1, v2, v3)
  * @param edge The edge to specify.
  * @param ti   The t index 0 or 1
  * @param v1   The first element of the cycle.
  * @param v2   The first element of the cycle.
  * @param v3   The first element of the cycle. */
  void set_orientation(int edge, int ti, int v1, int v2, int v3);

  /** Check if edge ei is oriented. */
  int check_oriented(int ei);

  /** Orient all the triangles in the graph consistently. */
  void orient();

  /** Compute the normal vectors for each edge in the graph. */
  void update_normals();

  /** Update the triangles, orient, and compute the normals. */
  void update_topology();

  /** Check if two edges are within fusing distance on both vertices. */
  int check_edge_dist(int ei, int ej);

  /** Check if edge is eligible for a wedge monomer insertion. */
  int check_wedge(int ei);

  /** Get the vertex to which vi should be paired in insertion attempt. */
  int get_wedge_pair(int vi, int vj);

  /** Set Nmono1, Nmono2, Nwedge */
  void update_surface_counts();

  /** Compute the harmonic energy for the specified edge.
   * \f$E_{\textrm{stretch}} = \frac{\epsilon}{2} (l_{\rm edge} - l_0)^2\f$
   * @param edge The index of the edge.
   * @return \f$E_{\textrm{stretch}}\f$ */
  double stretch_energy(int edge);

  /** Compute the harmonic energy for the specified edge.
   * \f$E_{\textrm{bend}} = \frac{\kappa}{2} 1-\cos(\theta_{\rm edge})\f$
   * @param edge The index of the edge.
   * @return \f$E_{\textrm{bend}}\f$ */
  double bend_energy(int edge);

  /** Compute the total energy. */
  double compute_energy();

  /** Compute the energy involving the vertex vi. */
  double vertex_energy(int vi);

  /** Check if edge ei causes a steric overlap with the graph */
  int check_overlap(int ei);

  /** Output a HOOMD-xml format file at time t */
  void dump_xml(int time);

  /** Get the center of mass of triangle defined by vi, vj, vk */
  void com(int vi, int vj, int vk, double *com);

  /** Get the inward facing normal vector associated with triangle defined by vi, vj, vk */
  void nvec(int vi, int vj, int vk, double *nveci, double s);

  /** Get the outward normal vector associated with triangle defined by vi, vj, vk */
  void outernvec(int vi, int vj, int vk, double *nveci, double s);
};

// other functions
/** Build a sheet of Nside x Nside monomers */
void make_sheet(graph &g, int Nside);

/** Make a single hexamer centered at the origin */
void make_hexamer(graph &g);

/** Make a single pentamer centered at the origin */
void make_pentamer(graph &g);

/** Make an icosahedral cap centered at the origin */
void make_cap(graph &g, int Nside);

/** Make a complete icosahedron for testing energies, order params */
void make_icosahedron(graph &g, int Nside);

/** Dump the number of fivefold and sevenfold defects to a file */
void dump_defects(graph &g, FILE *f, int time);

/** Dump the surface bonds to a file. */
void dump_surf_xml(graph &g, int time);

/** Dump the coordinates of the normal vectors in addition to the graph coordinates */
void dump_norm_xml(graph &g, int time);

void dump_xml_frame(graph &g, FILE *f, int time);

/** Dump a log file with all simulation parameters */
void dump_log(graph &g);

/** norm of vector v */
double norm(double *v);

/** dist between vectors v1 and v2 */
double dist(double *v1, double *v2);

/** dot product between vectors v1 and v2 */
double dot(double *v1, double *v2);

/** cross product between vectors v1 and v2 */
void cross(double *v1, double *v2, double *res);

/** subtract vector v2 from v1 */
void vsub(double *v1, double *v2, double *res);

#endif
