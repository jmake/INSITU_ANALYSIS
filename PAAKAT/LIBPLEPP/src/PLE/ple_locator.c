/*============================================================================
 * Locate points in a representation associated with a mesh
 *============================================================================*/

/*
  This file is part of the "Parallel Location and Exchange" library,
  intended to provide mesh or particle-based code coupling services.

  Copyright (C) 2005-2011  EDF

  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License as published by the Free Software Foundation; either
  version 2.1 of the License, or (at your option) any later version.

  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public
  License along with this library; if not, write to the Free Software
  Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
*/

/*!
 * \file ple_locator.c
 *
 * \brief Locate points in a representation associated with a mesh.
 */

/*----------------------------------------------------------------------------*/

#include "ple_config.h"

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#if defined(PLE_HAVE_MPI)
#include <mpi.h>
#if !defined(MPI_VERSION) /* Defined in up-to-date MPI versions */
#  define MPI_VERSION 1
#endif
#endif

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "ple_config_defs.h"
#include "ple_defs.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "ple_locator.h"

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force back Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */

/*============================================================================
 * Local macro definitions
 *============================================================================*/

/* In case <math.h> does not define HUGE_VAL, use a "safe" value */
#if !defined(HUGE_VAL)
#define HUGE_VAL 1.0e+30
#endif

/* Geometric operation macros*/

#define _MODULE(vect) \
  sqrt(vect[0] * vect[0] + vect[1] * vect[1] + vect[2] * vect[2])

/*============================================================================
 * Type definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Structure defining a locator
 *----------------------------------------------------------------------------*/

struct _ple_locator_t {

  /* Basic information */
  /*-------------------*/

  double    tolerance;         /* Associated tolerance */

  int       dim;               /* Spatial dimension */

  int       locate_closest;    /* 0: no special handling of unlocated points
                                  1: local (but not distant) points not located
                                     within the tolerance will be located on
                                     the closest element
                                  2: distant (and possibly local) points not
                                     located within the tolerance will be
                                     located on the closest element */

#if defined(PLE_HAVE_MPI)
  int       async_exchange;    /* flag for asynchronous variables exchange */
  MPI_Comm  comm;              /* Associated MPI communicator */
#endif

  int       n_ranks;           /* Number of MPI ranks of distant location */
  int       start_rank;        /* First MPI rank of distant location */

  int       n_intersects;      /* Number of intersecting distant ranks */
  int      *intersect_rank;    /* List of intersecting distant ranks */
  double   *intersect_extents; /* List of intersecting distant extents */

  ple_lnum_t   *local_points_idx;   /* Start index of local points per rank
                                       (size: n_intersects + 1)*/
  ple_lnum_t   *distant_points_idx; /* Start index of distant points per rank
                                       (size: n_intersects + 1)*/

  ple_lnum_t   *local_point_ids;        /* Local point index for data received
                                           (with blocs starting at
                                           local_points_idx[] indexes,
                                           0 to n-1 numbering) */

  ple_lnum_t   *distant_point_location; /* Location of distant points by parent
                                           element number (with blocs starting
                                           at distant_points_idx[] indexes) */
  ple_coord_t  *distant_point_coords;   /* Coordinates of distant points
                                           (with blocs starting at
                                           distant_points_idx[]*dim indexes) */

ple_coord_t  *distant_point_props;    // 2016Jun13  

  ple_lnum_t    n_interior;         /* Number of local points located */
  ple_lnum_t   *interior_list;      /* List (1 to n numbering) of points
                                       located */
  ple_lnum_t    n_exterior;         /* Number of local points not located */
  ple_lnum_t   *exterior_list;      /* List (1 to n numbering) of points
                                       not located */

  /* Timing information (2 or 4 fields/time; 0: total; 1: communication;
     2/3: closest point location total/communication) */

  double  location_wtime[4];       /* Location Wall-clock time */
  double  location_cpu_time[4];    /* Location CPU time */
  double  exchange_wtime[2];       /* Variable exchange Wall-clock time */
  double  exchange_cpu_time[2];    /* Variable exchange CPU time */
};

/*============================================================================
 * Local function pointer type documentation
 *============================================================================*/

#ifdef DOXYGEN_ONLY

/*!
 * \brief Query number of extents and compute extents of a mesh representation.
 *
 * For future optimizations, computation of extents should not be limited
 * to mesh extents, but to 1 to n extents, allowing different extent
 * refinements, from global mesh to individual element extents.
 *
 * The minimum required functionality for this function is to compute
 * whole mesh extents, but it could also return extents of individual
 * elements, or intermediate extents of mesh subdivisions or coarsened
 * elements. As such, it takes an argument indicating the maximum
 * local number of extents it should compute (based on the size of
 * the extents array argument), but returns the number of extents
 * really computed, which may be lower (usually 1 for mesh extents,
 * possibly even 0 if the local mesh is empty). If n_max_extents = 1,
 * the whole mesh extents should be computed.
 *
 * If n_max_extents is set to a negative value (-1), no extents are computed,
 * but the function returns the maximum number of extents it may compute.
 * This query mode allows for the caller to allocate the correct amount
 * of memory for a subsequent call.
 *
 * \param[in] mesh          pointer to mesh representation structure
 * \param[in] n_max_extents maximum number of sub-extents (such as element
 *                          extents) to compute, or -1 to query
 * \param[in] tolerance     addition to local extents of each element:
 *                          extent = base_extent * (1 + tolerance)
 * \param[in] extents       extents associated with the mesh or elements (or
 *                          even aggregated elements in case of coarser
 *                          representation):
 *                          x_min_0, y_min_0, ..., x_max_i, y_max_i, ...
 *                          (size: 2*dim*n_max_extents), ignored in query mode
 *
 * \return the number of extents computed
 */

typedef ple_lnum_t
(ple_mesh_extents_t) (const void  *mesh,
                      ple_lnum_t   n_max_extents,
                      double       tolerance,
                      double       extents[]);

/*!
 * \brief Find elements in a given local mesh containing points: updates the
 * location[] and distance[] arrays associated with a set of points
 * for points that are in an element of this mesh, or closer to one
 * than to previously encountered elements.
 *
 * param[in]      this_nodal   pointer to nodal mesh representation structure
 * param[in]      tolerance    associated tolerance
 * param[in]      n_points     number of points to locate
 * param[in]      point_coords point coordinates
 * param[in, out] location     number of element containing or closest to each
 *                             point (size: n_points)
 * param[in, out] distance     distance from point to element indicated by
 *                             location[]: < 0 if unlocated, 0 - 1 if inside,
 *                             and > 1 if outside a volume element, or
 *                             absolute distance to a surface element
 *                             (size: n_points)
 */

typedef void
(ple_mesh_elements_contain_t) (const void         *mesh,
                               double              tolerance,
                               ple_lnum_t          n_points,
                               const ple_coord_t   point_coords[],
                               ple_lnum_t          location[],
                               float               distance[]);

/*!
 * \brief Find elements in a given local mesh closest to points: updates the
 * location[] and distance[] arrays associated with a set of points for
 * points that are closer to an element of this mesh than to previously
 * encountered elements.
 *
 * parameters:
 *   mesh         <-- pointer to mesh representation structure
 *   n_points     <-- number of points to locate
 *   point_coords <-- point coordinates
 *   location     <-> number of element containing or closest to each point
 *                    (size: n_points)
 *   distance     <-> distance from point to element indicated by location[]:
 *                    < 0 if unlocated, or absolute distance to a surface
 *                    element (size: n_points)
 */

typedef void
(ple_mesh_elements_closest_t) (const void         *mesh,
                               ple_lnum_t          n_points,
                               const ple_coord_t   point_coords[],
                               ple_lnum_t          location[],
                               float               distance[]);

/*!
 * \brief Function pointer type for user definable logging/profiling
 * type functions
 */

typedef int
(ple_locator_log_t) (int         event,
                     int         data,
                     const char *string);

#endif /* DOXYGEN_ONLY */

/*============================================================================
 * Static global variables
 *============================================================================*/

#if defined(PLE_HAVE_MPI)

/* maximum number of exchanging ranks for which we use asynchronous
   MPI sends and receives instead of MPI_SendRecv */

static int _ple_locator_async_threshold = 0;

/* global logging function */

static ple_locator_log_t   *_ple_locator_log_func = NULL;

/* global variables associated with communication logging */

static int  _ple_locator_log_start_p_comm = 0;
static int  _ple_locator_log_end_p_comm = 0;
static int  _ple_locator_log_start_g_comm = 0;
static int  _ple_locator_log_end_g_comm = 0;

#endif

/*============================================================================
 * Private function definitions
 *============================================================================*/

#if defined(PLE_HAVE_MPI)

/*----------------------------------------------------------------------------
 * Log communication start.
 *
 * parameters:
 *   start_p_comm <-- event number for the start of the "send/recv"
 *                    or "global send/recv" state
 *   timing       <-> 0: wall-clock total; 1: CPU total;
 *                    2: wall-clock timer start; 3: CPU timer start
 *----------------------------------------------------------------------------*/

inline static void
_locator_trace_start_comm(int      start_p_comm,
                          double   timing[4])
{
  timing[2] = ple_timer_wtime();
  timing[3] = ple_timer_cpu_time();

  if(_ple_locator_log_func != NULL)
    _ple_locator_log_func(start_p_comm, 0, NULL);
}

/*----------------------------------------------------------------------------
 * Log communication end.
 *
 * parameters:
 *   end_p_comm  <-- event number for the end of the "send/recv" or
 *                   "global send/recv" state
 *   timing      <-> 0: wall-clock total; 1 CPU total;
 *                   2: wall-clock timer start; 3: CPU timer start
 *----------------------------------------------------------------------------*/

inline static void
_locator_trace_end_comm(int      end_p_comm,
                        double   timing[4])
{
  if(_ple_locator_log_func != NULL)
    _ple_locator_log_func(end_p_comm, 0, NULL);

  timing[0] += ple_timer_wtime() - timing[2];
  timing[1] += ple_timer_cpu_time() - timing[3];
}

#endif /* defined(PLE_HAVE_MPI) */

/*----------------------------------------------------------------------------
 * Test if two extents intersect
 *
 * parameters:
 *   dim             <-- spatial (coordinates) dimension
 *   extents_1       <-- extents: x_min, y_min, ..., x_max, y_max, ...
 *                       size: dim*2
 *   extents_2       <-- extents: x_min, y_min, ..., x_max, y_max, ...
 *                       size: dim*2
 *
 * returns:
 *   true if extents intersect, false otherwise
 *----------------------------------------------------------------------------*/

inline static _Bool
_intersect_extents(int           dim,
                   const double  extents_1[],
                   const double  extents_2[])
{
  int i;
  _Bool retval = true;

  for (i = 0; i < dim; i++) {
    if (   (extents_1[i] > extents_2[i + dim])
        || (extents_2[i] > extents_1[i + dim])) {
      retval = false;
      break;
    }
  }

  return retval;
}

/*----------------------------------------------------------------------------
 * Test if a point is within given extents
 *
 * parameters:
 *   dim             <-- spatial (coordinates) dimension
 *   coords          <-- coordinates: x, y, ...
 *                       size: dim
 *   extents         <-- extents: x_min, y_min, ..., x_max, y_max, ...
 *                       size: dim*2
 *
 * returns:
 *   true if point lies within extents, false otherwise
 *----------------------------------------------------------------------------*/

inline static _Bool
_within_extents(int                dim,
                const ple_coord_t  coords[],
                const double       extents[])
{
  int i;
  _Bool retval = true;

  for (i = 0; i < dim; i++) {
    if (   (coords[i] < extents[i])
        || (coords[i] > extents[i + dim])) {
      retval = false;
      break;
    }
  }

  return retval;
}

/*----------------------------------------------------------------------------
 * Compute minimum distances between points of different extents
 *
 * parameters:
 *   dim             <-- spatial (coordinates) dimension
 *   extents_1       <-- extents: x_min, y_min, ..., x_max, y_max, ...
 *                       size: dim*2
 *   extents_2       <-- extents: x_min, y_min, ..., x_max, y_max, ...
 *                       size: dim*2
 *   distances       <-- minimum and maximum estimated distance
 *                       between 2 points of extents_1 and extents_2
 *
 * returns:
 *   minimum possible distance between points of different extents
 *----------------------------------------------------------------------------*/

inline static double
_extents_min_distance(int           dim,
                      const double  extents_1[],
                      const double  extents_2[])
{
  int i;
  double _c_dist_min[3] = {HUGE_VAL, HUGE_VAL, HUGE_VAL};

  for (i = 0; i < dim; i++) {
    if (extents_1[i] > extents_2[i + dim])
      _c_dist_min[i] = extents_1[i] - extents_2[i + dim];
    else if (extents_2[i] > extents_1[i + dim])
      _c_dist_min[i] = extents_2[i] - extents_1[i + dim];
    else
      _c_dist_min[i] = 0.; /* Intersect */
  }

  return _MODULE(_c_dist_min);
}

/*----------------------------------------------------------------------------
 * Compute maximum distances between points of different extents
 *
 * parameters:
 *   dim             <-- spatial (coordinates) dimension
 *   extents_1       <-- extents: x_min, y_min, ..., x_max, y_max, ...
 *                       size: dim*2
 *   extents_2       <-- extents: x_min, y_min, ..., x_max, y_max, ...
 *                       size: dim*2
 *   distances       <-- minimum and maximum estimated distance
 *                       between 2 points of extents_1 and extents_2
 *
 * returns:
 *   maximum possible distance between points of different extents
 *----------------------------------------------------------------------------*/

inline static double
_extents_max_distance(int           dim,
                      const double  extents_1[],
                      const double  extents_2[])
{
  int i;
  double _c_dist_max[3] = {0., 0., 0.};

  for (i = 0; i < dim; i++) {

    if (extents_1[i] > extents_2[i + dim])
      _c_dist_max[i] = extents_1[i + dim] - extents_2[i];
    else if (extents_2[i] > extents_1[i + dim])
      _c_dist_max[i] = extents_2[i + dim] - extents_1[i];
    else {
      double d1 = fabs(extents_1[i + dim] - extents_2[i]);
      double d2 = fabs(extents_1[i] - extents_2[i + dim]);
      _c_dist_max[i] = PLE_MAX(d1, d2);
    }

  }

  return _MODULE(_c_dist_max);
}

/*----------------------------------------------------------------------------
 * Compute extents of a point set
 *
 * parameters:
 *   dim          <-- space dimension of points to locate
 *   n_points     <-- number of points to locate
 *   point_list   <-- optional indirection array to point_coords
 *                    (1 to n_points numbering)
 *   point_coords <-- coordinates of points to locate
 *                    (dimension: dim * n_points)
 *   extents      --> extents associated with mesh:
 *                    x_min, y_min, ..., x_max, y_max, ... (size: 2*dim)
 *----------------------------------------------------------------------------*/

static void
_point_extents(int                  dim,
               ple_lnum_t           n_points,
               const ple_lnum_t     point_list[],
               const ple_coord_t    point_coords[],
               double               extents[])
{
  int i;
  ple_lnum_t j, coord_idx;

  /* initialize extents in case mesh is empty or dim < 3 */
  for (i = 0; i < dim; i++) {
    extents[i]       =  HUGE_VAL;
    extents[i + dim] = -HUGE_VAL;
  }

  /* Compute extents */

  if (point_list != NULL) {

    for (j = 0; j < n_points; j++) {
      coord_idx = point_list[j] - 1;
      for (i = 0; i < dim; i++) {
        if (extents[i]       > point_coords[(coord_idx * dim) + i])
          extents[i]       = point_coords[(coord_idx * dim) + i];
        if (extents[i + dim] < point_coords[(coord_idx * dim) + i])
          extents[i + dim] = point_coords[(coord_idx * dim) + i];
      }
    }
  }

  else {

    for (coord_idx = 0; coord_idx < n_points; coord_idx++) {
      for (i = 0; i < dim; i++) {
        if (extents[i]       > point_coords[(coord_idx * dim) + i])
          extents[i]       = point_coords[(coord_idx * dim) + i];
        if (extents[i + dim] < point_coords[(coord_idx * dim) + i])
          extents[i + dim] = point_coords[(coord_idx * dim) + i];
      }
    }
  }

}

#if defined(PLE_HAVE_MPI)

/*----------------------------------------------------------------------------
 * Update intersection rank information once location is done
 * point set.
 *
 * parameters:
 *   this_locator      <-> pointer to locator structure
 *   n_points          <-- number of points to locate
 *   location_rank_id  <-> rank on which a point is located in, id
 *                         in updated locator intersection rank info out
 *                         (1 to n_points numbering)
 *----------------------------------------------------------------------------*/

static void
_update_location_ranks(ple_locator_t  *this_locator,
                       ple_lnum_t      n_points,
                       ple_lnum_t      location_rank_id[])
{
  int i, k;
  ple_lnum_t j;

  int loc_vals[2], max_vals[2];
  int comm_size;
  int *send_flag = NULL, *recv_flag = NULL, *intersect_rank_id = NULL;

  double comm_timing[4] = {0., 0., 0., 0.};

  /* Initialization */

  MPI_Comm_size(this_locator->comm, &comm_size);

  PLE_MALLOC(send_flag, comm_size, int);
  PLE_MALLOC(recv_flag, comm_size, int);

  for (i = 0; i < comm_size; i++)
    send_flag[i] = 0;

  for (j = 0; j < n_points; j++) {
    if (location_rank_id[j] > -1)
      send_flag[location_rank_id[j]] = 1;
  }

  /* As exchange detection is asymetric, synchronize it */

  _locator_trace_start_comm(_ple_locator_log_start_g_comm, comm_timing);

  MPI_Alltoall(send_flag, 1, MPI_INT, recv_flag, 1, MPI_INT,
               this_locator->comm);

  _locator_trace_end_comm(_ple_locator_log_end_g_comm, comm_timing);

  /* Update number of "intersects" and matching rank info */

  this_locator->n_intersects = 0;

  for (i = 0; i < this_locator->n_ranks; i++) {
    j = i + this_locator->start_rank;
    if (send_flag[j] == 1 || recv_flag[j] == 1)
      this_locator->n_intersects++;
  }

  PLE_REALLOC(this_locator->intersect_rank,
              this_locator->n_intersects,
              int);

  for (i = 0, k = 0; i < this_locator->n_ranks; i++) {
    j = i + this_locator->start_rank;
    if (send_flag[j] == 1 || recv_flag[j] == 1)
      this_locator->intersect_rank[k++] = j;
  }

  PLE_FREE(send_flag);
  PLE_FREE(recv_flag);

  /* Now convert location rank id to intersect rank
     (communication ordering not yet optimized) */

  PLE_MALLOC(intersect_rank_id, this_locator->n_ranks, int);

  for (i = 0; i < this_locator->n_ranks; i++)
    intersect_rank_id[i] = -1;

  for (i = 0; i < this_locator->n_intersects; i++) {
    intersect_rank_id[  this_locator->intersect_rank[i]
                      - this_locator->start_rank] = i;
  }

  for (j = 0; j < n_points; j++) {
    k = location_rank_id[j] - this_locator->start_rank;
    if (k > -1)
      location_rank_id[j] = intersect_rank_id[k];
  }

  PLE_FREE(intersect_rank_id);

  _locator_trace_start_comm(_ple_locator_log_start_g_comm, comm_timing);

  loc_vals[0] = this_locator->n_intersects;
  loc_vals[1] = _ple_locator_async_threshold;

  MPI_Allreduce(loc_vals, max_vals, 2, MPI_INT, MPI_MAX,
                this_locator->comm);

  if (max_vals[0] <= max_vals[1])
    this_locator->async_exchange = 1;
  else
    this_locator->async_exchange = 0;

  _locator_trace_end_comm(_ple_locator_log_end_g_comm, comm_timing);

  this_locator->location_wtime[1] += comm_timing[0];
  this_locator->location_cpu_time[1] += comm_timing[1];
}

/*----------------------------------------------------------------------------
 * Location of points not yet located on the closest elements.
 *
 * parameters:
 *   this_locator      <-> pointer to locator structure
 *   mesh              <-- pointer to mesh representation structure
 *   dim               <-- spatial dimension
 *   n_points          <-- number of points to locate
 *   point_list        <-- optional indirection array to point_coords
 *                         (1 to n_points numbering)
 *   point_coords      <-- coordinates of points to locate
 *                         (dimension: dim * n_points)
 *   location          <-> number of distant element containing or closest
 *                         to each point, or -1 (size: n_points)
 *   location_rank_id  <-> rank id for distant element containing or closest
 *                         to each point, or -1
 *   distance          <-> distance from point to element indicated by
 *                         location[]: < 0 if unlocated, 0 - 1 if inside,
 *                         > 1 if outside (size: n_points)
 *   mesh_extents_f    <-- function computing mesh or mesh subset extents
 *   locate_closest_f  <-- function locating the closest local elements
 *----------------------------------------------------------------------------*/

static void
_locate_closest_distant(ple_locator_t                *this_locator,
                        const void                   *mesh,
                        int                           dim,
                        ple_lnum_t                    n_points,
                        const ple_lnum_t              point_list[],
                        const ple_coord_t             point_coords[],
                        ple_lnum_t                    location[],
                        ple_lnum_t                    location_rank_id[],
                        float                         distance[],
                        ple_mesh_extents_t           *mesh_extents_f,
                        ple_mesh_elements_closest_t  *locate_closest_f)
{
  int i, k, stride2;
  int dist_rank, dist_index;
  ple_lnum_t j;
  ple_lnum_t n_coords_loc, n_coords_dist;
  ple_lnum_t *location_loc, *location_dist;
  ple_lnum_t *send_index;
  ple_coord_t *coords_dist, *send_coords;
  float *min_dist = NULL;
  float *distance_dist, *distance_loc;
  double w_start, w_end, cpu_start, cpu_end;
  double e_extents[6] = {1., 0., 0., 0., 0., 0.}; /* initially "impossible" */
  double p_extents[6];

  MPI_Status status;

  int comm_size = 0;
  ple_lnum_t _n_points = 0, n_exchanges = 0;
  double min_max_dist = HUGE_VAL;
  double comm_timing[4] = {0., 0., 0., 0.};

  int *send_flag = NULL, *recv_flag = NULL;
  ple_lnum_t *_point_ids = NULL, *_point_list = NULL, *exchange_rank = NULL;
  double *recvbuf = NULL;

  /* Initialize timing */

  w_start = ple_timer_wtime();
  cpu_start = ple_timer_cpu_time();

  /* Count points not yet located */

  if (this_locator->locate_closest > 0) {

    for (j = 0; j < n_points; j++) {
      if (location[j] == -1)
        _n_points += 1;
    }

    PLE_MALLOC(_point_ids, _n_points, ple_lnum_t);
    PLE_MALLOC(_point_list, _n_points, ple_lnum_t);

    if (point_list == NULL) {
      for (j = 0, _n_points = 0; j < n_points; j++) {
        if (location[j] == -1) {
          _point_ids[_n_points] = j;
          _point_list[_n_points] = j+1;
          _n_points++;
        }
      }
    }
    else {
      for (j = 0, _n_points = 0; j < n_points; j++) {
        if (location[j] == -1) {
          _point_ids[_n_points] = j;
          _point_list[_n_points] = point_list[j];
          _n_points++;
        }
      }
    }

  }

  /* Recompute mesh extents only if local location function is available */

  if (locate_closest_f != NULL)
    mesh_extents_f(mesh, 1, 0.1, e_extents);

  /* Recompute point extents using only points not yet located */

  _point_extents(dim,
                 _n_points,
                 _point_ids,
                 point_coords,
                 p_extents);

  /* Exchange extent information */

  MPI_Comm_size(this_locator->comm, &comm_size);

  stride2 = dim * 2; /* Stride for one type of extent */

  PLE_MALLOC(min_dist, this_locator->n_ranks, float);

  PLE_MALLOC(recvbuf, stride2*comm_size, double);

  _locator_trace_start_comm(_ple_locator_log_start_g_comm, comm_timing);

  MPI_Allgather(e_extents, stride2, MPI_DOUBLE, recvbuf, stride2, MPI_DOUBLE,
                this_locator->comm);

  _locator_trace_end_comm(_ple_locator_log_end_g_comm, comm_timing);

  /* Count and mark possible exchanges */

  n_exchanges = 0;

  PLE_MALLOC(exchange_rank, this_locator->n_ranks, int);

  PLE_MALLOC(send_flag, comm_size, int);
  PLE_MALLOC(recv_flag, comm_size, int);

  for (i = 0; i < comm_size; i++)
    send_flag[i] = 0;

  if (this_locator->locate_closest > 0) {

    for (i = 0; i < this_locator->n_ranks; i++) {
      j = this_locator->start_rank + i;
      double _max_dist = _extents_max_distance(dim,
                                               p_extents,
                                               recvbuf + (j*stride2));
      if (_max_dist < min_max_dist)
        min_max_dist = _max_dist;
    }

    for (i = 0; i < this_locator->n_ranks; i++) {
      j = this_locator->start_rank + i;
      double _min_dist = _extents_min_distance(dim,
                                               p_extents,
                                               recvbuf + (j*stride2));
      if (_min_dist > min_max_dist)
        send_flag[j] = 0;
      else if (_min_dist <= min_max_dist)
        send_flag[j] = 1;
      min_dist[i] = _min_dist;
    }

  }

  PLE_FREE(recvbuf);

  /* As exchange detection is asymetric, synchronize it */

  _locator_trace_start_comm(_ple_locator_log_start_g_comm, comm_timing);

  MPI_Alltoall(send_flag, 1, MPI_INT, recv_flag, 1, MPI_INT,
               this_locator->comm);

  _locator_trace_end_comm(_ple_locator_log_end_g_comm, comm_timing);

  for (i = 0; i < this_locator->n_ranks; i++) {
    j = i + this_locator->start_rank;
    if (send_flag[j] == 1 || recv_flag[j] == 1)
      exchange_rank[n_exchanges++] = j;
  }

  PLE_FREE(send_flag);
  PLE_FREE(recv_flag);

  /* Prepare communication */

  PLE_MALLOC(send_coords, _n_points * dim, ple_coord_t);
  PLE_MALLOC(send_index, _n_points, ple_lnum_t);

  /* First loop on distant ranks */
  /*-----------------------------*/

  for (i = 0; i < n_exchanges; i++) {

    double _min_dist;

    dist_index = i; /* Ordering (communication schema) not yet optimized */
    dist_rank  = exchange_rank[dist_index];

    /* Prepare and send coords that should fit in each send buffer */
    /* Reset buffers for current exchange rank */

    n_coords_loc = 0;
    _min_dist = min_dist[dist_rank - this_locator->start_rank];

    /* Build partial buffer */

    for (j = 0; j < _n_points; j++) {

      ple_lnum_t _point_id = _point_ids[j];
      ple_lnum_t coord_idx = _point_list[j] - 1;

      if (distance[_point_id] < -0.1 || distance[_point_id] > _min_dist) {

        send_index[n_coords_loc] = _point_id;
        for (k = 0; k < dim; k++)
          send_coords[n_coords_loc*dim + k]
            = point_coords[dim*coord_idx + k];

        n_coords_loc += 1;
      }

    }

    /* Send then receive partial buffer */

    _locator_trace_start_comm(_ple_locator_log_start_p_comm, comm_timing);

    MPI_Sendrecv(&n_coords_loc, 1, PLE_MPI_LNUM, dist_rank, PLE_MPI_TAG,
                 &n_coords_dist, 1, PLE_MPI_LNUM, dist_rank,
                 PLE_MPI_TAG, this_locator->comm, &status);

    _locator_trace_end_comm(_ple_locator_log_end_p_comm, comm_timing);

    PLE_MALLOC(coords_dist, n_coords_dist*dim, ple_coord_t);

    _locator_trace_start_comm(_ple_locator_log_start_p_comm, comm_timing);

    MPI_Sendrecv(send_coords, (int)(n_coords_loc*dim),
                 PLE_MPI_COORD, dist_rank, PLE_MPI_TAG,
                 coords_dist, (int)(n_coords_dist*dim),
                 PLE_MPI_COORD, dist_rank, PLE_MPI_TAG,
                 this_locator->comm, &status);

    _locator_trace_end_comm(_ple_locator_log_end_p_comm, comm_timing);

    /* Now locate received coords on local rank */

    PLE_MALLOC(location_dist, n_coords_dist, ple_lnum_t);
    PLE_MALLOC(distance_dist, n_coords_dist, float);

    for (j = 0; j < n_coords_dist; j++) {
      location_dist[j] = -1;
      distance_dist[j] = -1.0;
    }

    if (locate_closest_f != NULL)
      locate_closest_f(mesh,
                       n_coords_dist,
                       coords_dist,
                       location_dist,
                       distance_dist);

    PLE_FREE(coords_dist);

    /* Exchange location return information with distant rank */

    PLE_MALLOC(location_loc, n_coords_loc, ple_lnum_t);
    PLE_MALLOC(distance_loc, n_coords_loc, float);

    _locator_trace_start_comm(_ple_locator_log_start_p_comm, comm_timing);

    MPI_Sendrecv(location_dist, (int)n_coords_dist,
                 PLE_MPI_LNUM, dist_rank, PLE_MPI_TAG,
                 location_loc, (int)n_coords_loc,
                 PLE_MPI_LNUM, dist_rank, PLE_MPI_TAG,
                 this_locator->comm, &status);

    MPI_Sendrecv(distance_dist, (int)n_coords_dist,
                 MPI_FLOAT, dist_rank, PLE_MPI_TAG,
                 distance_loc, (int)n_coords_loc,
                 MPI_FLOAT, dist_rank, PLE_MPI_TAG,
                 this_locator->comm, &status);

    _locator_trace_end_comm(_ple_locator_log_end_p_comm, comm_timing);

    PLE_FREE(location_dist);
    PLE_FREE(distance_dist);

    /* Now update location information */

    for (j = 0; j < n_coords_loc; j++) {

      ple_lnum_t l = send_index[j];

      if (   (distance_loc[j] > -0.1)
          && (distance_loc[j] < distance[l] || distance[l] < -0.1)) {
        location_rank_id[l] = dist_rank;
        location[l] = location_loc[j];
        distance[l] = distance_loc[j];
      }

    }

    PLE_FREE(location_loc);
    PLE_FREE(distance_loc);
  }

  /* Free temporary arrays */

  PLE_FREE(send_coords);
  PLE_FREE(send_index);

  PLE_FREE(min_dist);
  PLE_FREE(exchange_rank);

  PLE_FREE(_point_ids);
  PLE_FREE(_point_list);

  /* Finalize timing */

  w_end = ple_timer_wtime();
  cpu_end = ple_timer_cpu_time();

  this_locator->location_wtime[2] += (w_end - w_start);
  this_locator->location_cpu_time[2] += (cpu_end - cpu_start);

  this_locator->location_wtime[1] += comm_timing[0];
  this_locator->location_cpu_time[1] += comm_timing[1];
  this_locator->location_wtime[3] += comm_timing[0];
  this_locator->location_cpu_time[3] += comm_timing[1];
}

/*----------------------------------------------------------------------------
 * Prepare locator for use with a given mesh representation and point set.
 *
 * parameters:
 *   this_locator      <-> pointer to locator structure
 *   mesh              <-- pointer to mesh representation structure
 *   dim               <-- spatial dimension
 *   n_points          <-- number of points to locate
 *   point_list        <-- optional indirection array to point_coords
 *                         (1 to n_points numbering)
 *   point_coords      <-- coordinates of points to locate
 *                         (dimension: dim * n_points)
 *   distance          --> optional distance from point to matching element:
 *                         < 0 if unlocated; 0 - 1 if inside and > 1 if
 *                         outside a volume element, or absolute distance
 *                         to a surface element (size: n_points)
 *   mesh_extents_f    <-- pointer to function computing mesh or mesh
 *                         subset or element extents
 *   locate_inside_f   <-- function locating points in or on elements
 *   locate_closest_f  <-- function locating the closest local elements
 *----------------------------------------------------------------------------*/

static void
_locate_all_distant(ple_locator_t                *this_locator,
                    const void                   *mesh,
                    int                           dim,
                    int                           dim_props,
                    ple_lnum_t                    n_points,
                    const ple_lnum_t              point_list[],
                    const ple_coord_t             point_coords[],
                    float                         distance[],
                    ple_mesh_extents_t           *mesh_extents_f,
                    ple_mesh_elements_contain_t  *locate_inside_f,
                    ple_mesh_elements_closest_t  *locate_closest_f,  
                    const ple_coord_t             point_props[]
                   )
{
  int i, k, stride;
  int dist_rank, dist_index;
  ple_lnum_t j;
  ple_lnum_t n_coords_loc, n_coords_dist, n_interior, n_exterior;
  ple_lnum_t coord_idx, start_idx;
  ple_lnum_t *location_loc, *location_dist;
  ple_lnum_t *location, *location_rank_id;
  ple_lnum_t *send_id, *send_location;
  ple_lnum_t *location_count, *location_shift;
  ple_coord_t *coords_dist, *send_coords;
ple_coord_t *send_props; // 2016Jun13  
  float *_distance = distance;
  float *distance_dist, *distance_loc;
  double extents[6];

  MPI_Status status;

  double comm_timing[4] = {0., 0., 0., 0.};

  /* Initialization */

  stride = dim * 2;

  PLE_MALLOC(send_props, n_points * dim_props, ple_coord_t);  // 2016Jun13 
 
  PLE_MALLOC(send_coords, n_points * dim, ple_coord_t);
  PLE_MALLOC(send_id, n_points, ple_lnum_t);

  PLE_MALLOC(location, n_points, ple_lnum_t);
  PLE_MALLOC(location_rank_id, n_points, ple_lnum_t);

  if (distance == NULL)
    PLE_MALLOC(_distance, n_points, float);

  for (j = 0; j < n_points; j++) {
    location[j] = -1;
    location_rank_id[j] = -1;
    _distance[j] = -1.0;
  }

  /* First loop on possibly intersecting distant ranks */
  /*---------------------------------------------------*/

  for (i = 0; i < this_locator->n_intersects; i++) {

    dist_index = i; /* Ordering (communication schema) not yet optimized */
    dist_rank  = this_locator->intersect_rank[dist_index];

    /* Prepare and send coords that should fit in each send buffer */
    /* Reset buffers for current intersect rank */

    n_coords_loc = 0;

    for (k = 0; k < stride; k++)
      extents[k] = this_locator->intersect_extents[dist_index*stride + k];

    /* Build partial buffer */

    for (j = 0; j < n_points; j++) {

      if (point_list != NULL)
        coord_idx = point_list[j] - 1;
      else
        coord_idx = j;

      if (_within_extents(dim,
                          &(point_coords[dim*coord_idx]),
                          extents) == true) {

        send_id[n_coords_loc] = j;
        for (k = 0; k < dim; k++)
          send_coords[n_coords_loc*dim + k]
            = point_coords[dim*coord_idx + k];

        n_coords_loc += 1;
      }

    }

    /* Send then receive partial buffer */

    dist_rank = this_locator->intersect_rank[dist_index];

    _locator_trace_start_comm(_ple_locator_log_start_p_comm, comm_timing);

    MPI_Sendrecv(&n_coords_loc, 1, PLE_MPI_LNUM, dist_rank, PLE_MPI_TAG,
                 &n_coords_dist, 1, PLE_MPI_LNUM, dist_rank,
                 PLE_MPI_TAG, this_locator->comm, &status);

    _locator_trace_end_comm(_ple_locator_log_end_p_comm, comm_timing);

    PLE_MALLOC(coords_dist, n_coords_dist*dim, ple_coord_t);

    _locator_trace_start_comm(_ple_locator_log_start_p_comm, comm_timing);

    MPI_Sendrecv(send_coords, (int)(n_coords_loc*dim),
                 PLE_MPI_COORD, dist_rank, PLE_MPI_TAG,
                 coords_dist, (int)(n_coords_dist*dim),
                 PLE_MPI_COORD, dist_rank, PLE_MPI_TAG,
                 this_locator->comm, &status);

    _locator_trace_end_comm(_ple_locator_log_end_p_comm, comm_timing);

    /* Now locate received coords on local rank */

    PLE_MALLOC(location_dist, n_coords_dist, ple_lnum_t);
    PLE_MALLOC(distance_dist, n_coords_dist, float);

    for (j = 0; j < n_coords_dist; j++) {
      location_dist[j] = -1;
      distance_dist[j] = -1.0;
    }

    locate_inside_f(mesh,
                    this_locator->tolerance,
                    n_coords_dist,
                    coords_dist,
                    location_dist,
                    distance_dist);

    PLE_FREE(coords_dist);

    /* Exchange location return information with distant rank */

    PLE_MALLOC(location_loc, n_coords_loc, ple_lnum_t);
    PLE_MALLOC(distance_loc, n_coords_loc, float);

    _locator_trace_start_comm(_ple_locator_log_start_p_comm, comm_timing);

    MPI_Sendrecv(location_dist, (int)n_coords_dist,
                 PLE_MPI_LNUM, dist_rank, PLE_MPI_TAG,
                 location_loc, (int)n_coords_loc,
                 PLE_MPI_LNUM, dist_rank, PLE_MPI_TAG,
                 this_locator->comm, &status);

    MPI_Sendrecv(distance_dist, (int)n_coords_dist,
                 MPI_FLOAT, dist_rank, PLE_MPI_TAG,
                 distance_loc, (int)n_coords_loc,
                 MPI_FLOAT, dist_rank, PLE_MPI_TAG,
                 this_locator->comm, &status);

    _locator_trace_end_comm(_ple_locator_log_end_p_comm, comm_timing);

    PLE_FREE(location_dist);
    PLE_FREE(distance_dist);

    /* Now update location information */

    for (j = 0; j < n_coords_loc; j++) {

      ple_lnum_t l = send_id[j];

      if (   (distance_loc[j] > -0.1)
          && (distance_loc[j] < _distance[l] || _distance[l] < -0.1)) {
        location_rank_id[l] = dist_rank;
        location[l] = location_loc[j];
        _distance[l] = distance_loc[j];
      }

    }

    PLE_FREE(location_loc);
    PLE_FREE(distance_loc);

  }

  /* If option activated, second search to associate points not yet
     located to the closest elements */

  if (this_locator->locate_closest > 0)
    _locate_closest_distant(this_locator,
                            mesh,
                            dim,
                            n_points,
                            point_list,
                            point_coords,
                            location,
                            location_rank_id,
                            _distance,
                            mesh_extents_f,
                            locate_closest_f);

  /* Update info on communicating ranks and matching ids */
  /*----------------------------------------------------*/

  _update_location_ranks(this_locator,
                         n_points,
                         location_rank_id);

  /* Reorganize location information */
  /*---------------------------------*/

  /* Now that location is done, the location[] array contains
     either -1 if a point was not located, or a local index
     (associated with the corresponding rank); the distance[] array
     is not needed anymore now that all comparisons have been done */

  if (_distance != distance)
    PLE_FREE(_distance);

  PLE_MALLOC(location_shift, this_locator->n_intersects, ple_lnum_t);
  PLE_MALLOC(location_count, this_locator->n_intersects, ple_lnum_t);

  for (i = 0; i < this_locator->n_intersects; i++)
    location_count[i] = 0;

  n_exterior = 0;
  for (j = 0; j < n_points; j++) {
    if (location_rank_id[j] > -1)
      location_count[location_rank_id[j]] += 1;
    else
      n_exterior += 1;
  }

  this_locator->n_interior = n_points - n_exterior;
  PLE_MALLOC(this_locator->interior_list, this_locator->n_interior, ple_lnum_t);

  this_locator->n_exterior = n_exterior;
  PLE_MALLOC(this_locator->exterior_list, this_locator->n_exterior, ple_lnum_t);

  if (this_locator->n_intersects > 0)
    location_shift[0] = 0;
  for (i = 1; i < this_locator->n_intersects; i++)
    location_shift[i] = location_shift[i-1] + location_count[i-1];

  for (i = 0; i < this_locator->n_intersects; i++)
    location_count[i] = 0;

  /* send_id[] will now contain information for all blocks */
  for (j = 0; j < n_points; j++)
    send_id[j] = -1;

  PLE_MALLOC(send_location, n_points, ple_lnum_t);

  n_interior = 0;
  n_exterior = 0;
  for (j = 0; j < n_points; j++) {
    const int l_rank = location_rank_id[j];
    if (l_rank > -1) {
      send_id[location_shift[l_rank] + location_count[l_rank]] = j;
      location_count[l_rank] += 1;
      this_locator->interior_list[n_interior] = j + 1;
      n_interior += 1;
    }
    else {
      this_locator->exterior_list[n_exterior] = j + 1;
      n_exterior += 1;
    }
  }

  /* Second loop on possibly intersecting distant ranks */
  /*----------------------------------------------------*/

  /* Count and organize total number of local and distant points */

  PLE_MALLOC(this_locator->local_points_idx,
             this_locator->n_intersects + 1,
             ple_lnum_t);

  PLE_MALLOC(this_locator->distant_points_idx,
             this_locator->n_intersects + 1,
             ple_lnum_t);

  this_locator->local_points_idx[0] = 0;
  this_locator->distant_points_idx[0] = 0;

  for (i = 0; i < this_locator->n_intersects; i++) {

    dist_index = i; /* Ordering (communication schema) not yet optimized */
    dist_rank = this_locator->intersect_rank[dist_index];

    n_coords_loc = location_count[i];

    this_locator->local_points_idx[i+1]
      = this_locator->local_points_idx[i] + n_coords_loc;

    _locator_trace_start_comm(_ple_locator_log_start_p_comm, comm_timing);

    MPI_Sendrecv(&n_coords_loc, 1, PLE_MPI_LNUM, dist_rank, PLE_MPI_TAG,
                 &n_coords_dist, 1, PLE_MPI_LNUM, dist_rank,
                 PLE_MPI_TAG, this_locator->comm, &status);

    _locator_trace_end_comm(_ple_locator_log_end_p_comm, comm_timing);

    this_locator->distant_points_idx[i+1]
      = this_locator->distant_points_idx[i] + n_coords_dist;

  }

  /* Third loop on possibly intersecting distant ranks */
  /*----------------------------------------------------*/

  PLE_MALLOC(this_locator->local_point_ids,
             this_locator->local_points_idx[this_locator->n_intersects],
             ple_lnum_t);

  PLE_MALLOC(this_locator->distant_point_location,
             this_locator->distant_points_idx[this_locator->n_intersects],
             ple_lnum_t);

  PLE_MALLOC(this_locator->distant_point_coords,
             this_locator->distant_points_idx[this_locator->n_intersects] * dim,
             ple_coord_t);

PLE_MALLOC(this_locator->distant_point_props,                                 // 2016Jun13, ALLOCATE  
           this_locator->distant_points_idx[this_locator->n_intersects] * dim_props,
           ple_coord_t);

  for (i = 0; i < this_locator->n_intersects; i++) {

    dist_index = i; /* Ordering (communication schema) not yet optimized */
    dist_rank = this_locator->intersect_rank[dist_index];

    n_coords_loc =    this_locator->local_points_idx[i+1]
                    - this_locator->local_points_idx[i];

    n_coords_dist =   this_locator->distant_points_idx[i+1]
                    - this_locator->distant_points_idx[i];

    start_idx = this_locator->local_points_idx[i];

    for (j = 0; j < n_coords_loc; j++) {

      coord_idx = send_id[location_shift[i] + j];
      this_locator->local_point_ids[start_idx + j] = coord_idx;
      send_location[j] = location[coord_idx];
      if (point_list != NULL) {
        for (k = 0; k < dim; k++)
          send_coords[j*dim + k]
            = point_coords[dim*(point_list[coord_idx] - 1) + k];
      }
      else {
        for (k = 0; k < dim; k++)
        {
          send_coords[j*dim + k] = point_coords[dim*coord_idx + k];
        }  
        for (k = 0; k < dim_props; k++)
        {
          send_props[ j*dim_props + k] = point_props[ dim_props*coord_idx + k]; // 2016Jun13, copy 
        }  
      }
    }

    _locator_trace_start_comm(_ple_locator_log_start_p_comm, comm_timing);

    MPI_Sendrecv(send_location, (int)n_coords_loc,
                 PLE_MPI_LNUM, dist_rank, PLE_MPI_TAG,
                 (this_locator->distant_point_location
                  + this_locator->distant_points_idx[i]), (int)n_coords_dist,
                 PLE_MPI_LNUM, dist_rank, PLE_MPI_TAG,
                 this_locator->comm, &status);

    MPI_Sendrecv(send_coords, (int)(n_coords_loc*dim),
                 PLE_MPI_COORD, dist_rank, PLE_MPI_TAG,
                 (this_locator->distant_point_coords
                  + (this_locator->distant_points_idx[i]*dim)),
                 (int)(n_coords_dist*dim),
                 PLE_MPI_COORD, dist_rank, PLE_MPI_TAG,
                 this_locator->comm, &status);

MPI_Sendrecv(send_props, (int)(n_coords_loc*dim_props),  // 2016Jun13, MPI_Sendrecv  
             PLE_MPI_COORD, dist_rank, PLE_MPI_TAG,
             (this_locator->distant_point_props + (this_locator->distant_points_idx[i]*dim_props)),
             (int)(n_coords_dist*dim_props),
             PLE_MPI_COORD, dist_rank, PLE_MPI_TAG,
             this_locator->comm, &status);

    _locator_trace_end_comm(_ple_locator_log_end_p_comm, comm_timing);

  }

  PLE_FREE(location_count);
  PLE_FREE(location_shift);

  PLE_FREE(send_id);
  PLE_FREE(send_location);
  PLE_FREE(send_coords);
PLE_FREE(send_props); // 2018Mar29 - fix memory leak - matias

  PLE_FREE(location_rank_id);

  PLE_FREE(location);

  this_locator->location_wtime[1] += comm_timing[0];
  this_locator->location_cpu_time[1] += comm_timing[1];
}

#endif /* defined(PLE_HAVE_MPI) */

/*----------------------------------------------------------------------------
 * Location of points not yet located on the closest elements.
 *
 * parameters:
 *   this_locator      <-> pointer to locator structure
 *   mesh              <-- pointer to mesh representation structure
 *   dim               <-- spatial dimension
 *   n_points          <-- number of points to locate
 *   point_list        <-- optional indirection array to point_coords
 *                         (1 to n_points numbering)
 *   point_coords      <-- coordinates of points to locate
 *                         (dimension: dim * n_points)
 *   location          <-> number of distant element containing or closest
 *                         to each point, or -1 (size: n_points)
 *   distance          <-> distance from point to element indicated by
 *                         location[]: < 0 if unlocated, 0 - 1 if inside,
 *                         > 1 if outside (size: n_points)
 *   locate_inside_f   <-- function locating points in or on elements
 *   locate_closest_f  <-- function locating the closest local elements
 *----------------------------------------------------------------------------*/

static void
_locate_closest_local(ple_locator_t                *this_locator,
                      const void                   *mesh,
                      int                           dim,
                      ple_lnum_t                    n_points,
                      const ple_lnum_t              point_list[],
                      const ple_coord_t             point_coords[],
                      ple_lnum_t                    location[],
                      float                         distance[],
                      ple_mesh_elements_closest_t  *locate_closest_f)
{
  ple_lnum_t i, j;
  double w_start, w_end, cpu_start, cpu_end;

  ple_lnum_t _n_points = 0;
  ple_lnum_t  *_location = NULL;
  ple_coord_t *_point_coords = NULL;
  float       *_distance = NULL;

  /* Initialize timing */

  w_start = ple_timer_wtime();
  cpu_start = ple_timer_cpu_time();

  /* Count points not yet located */

  for (i = 0; i < n_points; i++) {
    if (location[i] == -1)
      _n_points += 1;
  }

  PLE_MALLOC(_location, _n_points, ple_lnum_t);
  PLE_MALLOC(_distance, _n_points, float);
  PLE_MALLOC(_point_coords, _n_points*dim, ple_coord_t);

  for (i = 0; i < n_points; i++) {
    _location[i] = -1;
    _distance[i] = -1;
  }

  if (point_list == NULL) {
    for (i = 0, _n_points = 0; i < n_points; i++) {
      if (location[i] == -1) {
        for (j = 0; j < dim; j++)
          _point_coords[_n_points*dim + j] = point_coords[i*dim + j];
        _n_points++;
      }
    }
  }
  else {
    for (j = 0, _n_points = 0; j < n_points; j++) {
      if (location[i] == -1) {
        ple_lnum_t k = point_list[j] - 1;
        for (j = 0; j < dim; j++)
          _point_coords[_n_points*dim + j] = point_coords[k*dim + j];
        _n_points++;
      }
    }
  }

  locate_closest_f(mesh,
                   _n_points,
                   _point_coords,
                   _location,
                   _distance);

  /* Copy location information to input/output arrays */

  for (i = 0, _n_points = 0; i < n_points; i++) {
    if (location[i] == -1) {
      location[i] = _location[_n_points];
      distance[i] = _distance[_n_points];
      _n_points++;
    }
  }

  PLE_FREE(_location);
  PLE_FREE(_distance);
  PLE_FREE(_point_coords);

  /* Finalize timing */

  w_end = ple_timer_wtime();
  cpu_end = ple_timer_cpu_time();

  this_locator->location_wtime[2] += (w_end - w_start);
  this_locator->location_cpu_time[2] += (cpu_end - cpu_start);
}

/*----------------------------------------------------------------------------
 * Prepare locator for use with a given mesh representation and point set.
 *
 * parameters:
 *   this_locator      <-> pointer to locator structure
 *   mesh              <-- pointer to mesh representation structure
 *   dim               <-- spatial dimension
 *   n_points          <-- number of points to locate
 *   point_list        <-- optional indirection array to point_coords
 *                         (1 to n_points numbering)
 *   point_coords      <-- coordinates of points to locate
 *                         (dimension: dim * n_points)
 *   distance          --> optional distance from point to matching element:
 *                         < 0 if unlocated; 0 - 1 if inside and > 1 if
 *                         outside a volume element, or absolute distance
 *                         to a surface element (size: n_points)
 *   locate_closest_f  <-- function locating the closest local elements
 *----------------------------------------------------------------------------*/

static void
_locate_all_local(ple_locator_t                *this_locator,
                  const void                   *mesh,
                  int                           dim,
                  int                           dim_props,
                  ple_lnum_t                    n_points,
                  const ple_lnum_t              point_list[],
                  const ple_coord_t             point_coords[],
                  float                         distance[],
                  ple_mesh_elements_contain_t  *locate_inside_f,
                  ple_mesh_elements_closest_t  *locate_closest_f, 
                  const ple_coord_t             point_props[]
                 )
{
  int l, stride;
  ple_lnum_t j, k;
  ple_lnum_t n_coords, n_interior, n_exterior, coord_idx;
  ple_lnum_t *location;
  ple_lnum_t location_count;
  ple_coord_t *coords;
  float *_distance = distance;
  double extents[6];

  /* Initialization */

  stride = dim * 2;

  PLE_MALLOC(coords, n_points * dim, ple_coord_t);

  /* Initialize location information */
  /*---------------------------------*/

  n_coords = 0;

  for (k = 0; k < stride; k++)
    extents[k] = this_locator->intersect_extents[k];

  /* Build partial buffer */

  for (j = 0; j < n_points; j++) {

    if (point_list != NULL)
      coord_idx = point_list[j] - 1;
    else
      coord_idx = j;

    if (_within_extents(dim,
                        &(point_coords[dim*coord_idx]),
                        extents) == true) {

      for (k = 0; k < dim; k++)
        coords[n_coords*dim + k]
          = point_coords[dim*coord_idx + k];

      n_coords += 1;
    }

  }

  PLE_REALLOC(coords, n_coords * dim, ple_coord_t);

 /*  Now locate coords */

  PLE_MALLOC(location, n_coords, ple_lnum_t);

  if (distance == NULL)
    PLE_MALLOC(_distance, n_coords, float);

  for (j = 0; j < n_coords; j++) {
    location[j] = -1;
    _distance[j] = -1.0;
  }

  locate_inside_f(mesh,
                  this_locator->tolerance,
                  n_coords,
                  coords,
                  location,
                  _distance);

  if (this_locator->locate_closest > 0)
    _locate_closest_local(this_locator,
                          mesh,
                          dim,
                          n_points,
                          point_list,
                          point_coords,
                          location,
                          _distance,
                          locate_closest_f);

  /* Reorganize location information */
  /*---------------------------------*/

  /* Now that location is done, the location[] array contains
     either -1 if a point was not located, or a local index;
     the distance[] array is not needed anymore now that all comparisons have
     been done */

  if (_distance != distance)
    PLE_FREE(_distance);

  location_count = 0;

  n_exterior = 0;
  for (j = 0; j < n_coords; j++) {
    if (location[j] > -1)
      location_count += 1;
    else
      n_exterior += 1;
  }

  this_locator->n_interior = n_coords - n_exterior;
  PLE_MALLOC(this_locator->interior_list, this_locator->n_interior, ple_lnum_t);

  this_locator->n_exterior = (n_points - n_coords) + n_exterior;
  PLE_MALLOC(this_locator->exterior_list, this_locator->n_exterior, ple_lnum_t);

  /* Organize total number of "local" and "distant" points */

  PLE_MALLOC(this_locator->local_points_idx, 2, ple_lnum_t);
  PLE_MALLOC(this_locator->distant_points_idx, 2, ple_lnum_t);

  this_locator->local_points_idx[0] = 0;
  this_locator->local_points_idx[1] = location_count;

  this_locator->distant_points_idx[0] = 0;
  this_locator->distant_points_idx[1] = location_count;

  this_locator->local_point_ids = NULL; /* Not needed for single-process */

  PLE_MALLOC(this_locator->distant_point_location, location_count, ple_lnum_t);
  PLE_MALLOC(this_locator->distant_point_coords, n_coords * dim, ple_coord_t);

PLE_MALLOC(this_locator->distant_point_props, n_coords * dim_props, ple_coord_t); // 2016Jun13, _locate_all_local  


  location_count = 0;
  n_interior = 0;
  n_exterior = 0;

  for (j = 0, k = 0; j < n_points; j++) {

    if (point_list != NULL)
      coord_idx = point_list[j] - 1;
    else
      coord_idx = j;

    if (_within_extents(dim,
                        &(point_coords[dim*coord_idx]),
                        extents) == true) {

      if (location[k] > -1) {
        this_locator->distant_point_location[location_count] = location[k];
        for (l = 0; l < dim; l++) {
          this_locator->distant_point_coords[location_count*dim + l] = point_coords[coord_idx*dim + l];
          this_locator->distant_point_props[ location_count*dim_props + l] = point_props[ coord_idx*dim_props + l]; // 2016Jun13, _locate_all_local Reorganize  
        }
        location_count += 1;
        this_locator->interior_list[n_interior] = j + 1;
        n_interior += 1;

      }
      else {
        this_locator->exterior_list[n_exterior] = j + 1;
        n_exterior += 1;
      }

      k += 1;

    }
    else {
      this_locator->exterior_list[n_exterior] = j + 1;
      n_exterior += 1;
    }

  }

  PLE_FREE(location);
  PLE_FREE(coords);
}

#if defined(PLE_HAVE_MPI)

/*----------------------------------------------------------------------------
 * Distribute variable defined on distant points to processes owning
 * the original points (i.e. distant processes).
 *
 * The exchange is symmetric if both variables are defined, receive
 * only if distant_var is NULL, or send only if local_var is NULL.
 *
 * parameters:
 *   this_locator  <-- pointer to locator structure
 *   distant_var   <-> variable defined on distant points (ready to send)
 *   local_var     <-> variable defined on local points (received)
 *   local_list    <-- optional indirection list (1 to n) for local_var
 *   datatype      <-- variable type
 *   stride        <-- dimension (1 for scalar, 3 for interlaced vector)
 *   reverse       <-- if true, exchange is reversed
 *                     (receive values associated with distant points
 *                     from the processes owning the original points)
 *----------------------------------------------------------------------------*/

static void
_exchange_point_var_distant(ple_locator_t     *this_locator,
                            void              *distant_var,
                            void              *local_var,
                            const ple_lnum_t  *local_list,
                            MPI_Datatype       datatype,
                            size_t             stride,
                            _Bool              reverse)
{
  int i, dist_v_count, loc_v_count, size;
  int dist_rank, dist_index;
  int dist_v_flag, loc_v_flag;
  ple_lnum_t n_points_loc, n_points_loc_max, n_points_dist;
  size_t dist_v_idx;
  void *dist_v_ptr;
  void *loc_v_buf;

  double comm_timing[4] = {0., 0., 0., 0.};

  MPI_Aint lb, extent;
  MPI_Status status;

  /* Check extent of datatype */

#if (MPI_VERSION >= 2)
  MPI_Type_get_extent(datatype, &lb, &extent);
#else
  MPI_Type_extent(datatype, &extent);
#endif
  MPI_Type_size(datatype, &size);

  if (extent != size)
    ple_error(__FILE__, __LINE__, 0,
              _("_exchange_point_var() is not implemented for use with\n"
                "MPI datatypes associated with structures using padding\n"
                "(for which size != extent)."));

  /* Initialization */

  n_points_loc_max = 0;

  for (i = 0; i < this_locator->n_intersects; i++) {
    n_points_loc =    this_locator->local_points_idx[i+1]
                    - this_locator->local_points_idx[i];
    if (n_points_loc > n_points_loc_max)
      n_points_loc_max = n_points_loc;
  }

  PLE_MALLOC(loc_v_buf, n_points_loc_max*size*stride, char);

  /* Loop on MPI ranks */
  /*-------------------*/

  for (i = 0; i < this_locator->n_intersects; i++) {

    const ple_lnum_t *_local_point_ids
      = this_locator->local_point_ids + this_locator->local_points_idx[i];

    dist_index = i; /* Ordering (communication schema) not yet optimized */
    dist_rank = this_locator->intersect_rank[dist_index];

    n_points_loc =    this_locator->local_points_idx[i+1]
                    - this_locator->local_points_idx[i];

    n_points_dist =   this_locator->distant_points_idx[i+1]
                    - this_locator->distant_points_idx[i];

    if (distant_var != NULL && n_points_dist > 0)
      dist_v_flag = 1;
    else
      dist_v_flag = 0;

    _locator_trace_start_comm(_ple_locator_log_start_p_comm, comm_timing);

    MPI_Sendrecv(&dist_v_flag, 1, MPI_INT, dist_rank, PLE_MPI_TAG,
                 &loc_v_flag, 1, MPI_INT, dist_rank, PLE_MPI_TAG,
                 this_locator->comm, &status);

    _locator_trace_end_comm(_ple_locator_log_end_p_comm, comm_timing);

    if (loc_v_flag == 1 && (local_var == NULL || n_points_loc == 0))
      ple_error(__FILE__, __LINE__, 0,
                _("Incoherent arguments to different instances in "
                  "_exchange_point_var().\n"
                  "Send and receive operations do not match "
                  "(dist_rank = %d\n)\n"), dist_rank);

    dist_v_idx = this_locator->distant_points_idx[i] * stride*size;
    dist_v_count = n_points_dist * stride * dist_v_flag;

    if (loc_v_flag > 0)
      loc_v_count = n_points_loc*stride;
    else
      loc_v_count = 0;

    /* Exchange information */

    if (distant_var != NULL)
      dist_v_ptr = (void *)(((char *)distant_var) + dist_v_idx);
    else
      dist_v_ptr = NULL;

    if (reverse == false) {

      _locator_trace_start_comm(_ple_locator_log_start_p_comm, comm_timing);

      MPI_Sendrecv(dist_v_ptr, dist_v_count, datatype, dist_rank, PLE_MPI_TAG,
                   loc_v_buf, loc_v_count, datatype, dist_rank, PLE_MPI_TAG,
                   this_locator->comm, &status);

      _locator_trace_end_comm(_ple_locator_log_end_p_comm, comm_timing);

      if (loc_v_flag > 0) {
        if (local_list == NULL) {
          int k;
          size_t l;
          const size_t nbytes = stride*size;
          for (k = 0; k < n_points_loc; k++) {
            char *local_v_p = (char *)local_var + _local_point_ids[k]*nbytes;
            const char *loc_v_buf_p = (const char *)loc_v_buf + k*nbytes;
            for (l = 0; l < nbytes; l++)
              local_v_p[l] = loc_v_buf_p[l];
          }
        }
        else {
          int k;
          size_t l;
          const size_t nbytes = stride*size;
          for (k = 0; k < n_points_loc; k++) {
            char *local_v_p =   (char *)local_var
                              + (local_list[_local_point_ids[k]] - 1)*nbytes;
            const char *loc_v_buf_p = (const char *)loc_v_buf + k*nbytes;
            for (l = 0; l < nbytes; l++)
              local_v_p[l] = loc_v_buf_p[l];
          }
        }
      }

    }
    else { /* if (reverse == true) */

      if (loc_v_flag > 0) {
        if (local_list == NULL) {
          int k;
          size_t l;
          const size_t nbytes = stride*size;
          for (k = 0; k < n_points_loc; k++) {
            const char *local_v_p
              = (const char *)local_var + _local_point_ids[k]*nbytes;
            char *loc_v_buf_p = (char *)loc_v_buf + k*nbytes;
            for (l = 0; l < nbytes; l++)
              loc_v_buf_p[l] = local_v_p[l];
          }
        }
        else {
          int k;
          size_t l;
          const size_t nbytes = stride*size;
          for (k = 0; k < n_points_loc; k++) {
            const char *local_v_p
              = (const char *)local_var
                + (local_list[_local_point_ids[k]] - 1)*nbytes;
            char *loc_v_buf_p = (char *)loc_v_buf + k*nbytes;
            for (l = 0; l < nbytes; l++)
              loc_v_buf_p[l] = local_v_p[l];
          }
        }
      }

      _locator_trace_start_comm(_ple_locator_log_start_p_comm, comm_timing);

      MPI_Sendrecv(loc_v_buf, loc_v_count, datatype, dist_rank, PLE_MPI_TAG,
                   dist_v_ptr, dist_v_count, datatype, dist_rank, PLE_MPI_TAG,
                   this_locator->comm, &status);

      _locator_trace_end_comm(_ple_locator_log_end_p_comm, comm_timing);

    }

  } /* End of loop on MPI ranks */

  PLE_FREE(loc_v_buf);

  this_locator->exchange_wtime[1] += comm_timing[0];
  this_locator->exchange_cpu_time[1] += comm_timing[1];
}

/*----------------------------------------------------------------------------
 * Distribute variable defined on distant points to processes owning
 * the original points (i.e. distant processes).
 *
 * The exchange is symmetric if both variables are defined, receive
 * only if distant_var is NULL, or send only if local_var is NULL.
 *
 * This variant of the function uses asynchronous MPI calls.
 *
 * parameters:
 *   this_locator  <-- pointer to locator structure
 *   distant_var   <-> variable defined on distant points (ready to send)
 *   local_var     <-> variable defined on local points (received)
 *   local_list    <-- optional indirection list (1 to n) for local_var
 *   datatype      <-- variable type
 *   stride        <-- dimension (1 for scalar, 3 for interlaced vector)
 *   reverse       <-- if true, exchange is reversed
 *                     (receive values associated with distant points
 *                     from the processes owning the original points)
 *----------------------------------------------------------------------------*/

static void
_exchange_point_var_distant_asyn(ple_locator_t     *this_locator,
                                 void              *distant_var,
                                 void              *local_var,
                                 const ple_lnum_t  *local_list,
                                 MPI_Datatype       datatype,
                                 size_t             stride,
                                 _Bool              reverse)
{
  int i, dist_v_count, loc_v_count, size;
  int dist_rank;
  ple_lnum_t n_points_loc, n_points_loc_tot, n_points_dist;
  size_t dist_v_idx;
  unsigned char *dist_v_ptr, *loc_v_ptr;

  MPI_Aint lb, extent;
  void *loc_v_buf = NULL;
  int *dist_v_flag = NULL, *loc_v_flag = NULL;
  MPI_Status *status = NULL;
  MPI_Request *request = NULL;

  double comm_timing[4] = {0., 0., 0., 0.};

  /* Check extent of datatype */

#if (MPI_VERSION >= 2)
  MPI_Type_get_extent(datatype, &lb, &extent);
#else
  MPI_Type_extent(datatype, &extent);
#endif
  MPI_Type_size(datatype, &size);

  if (extent != size)
    ple_error(__FILE__, __LINE__, 0,
              _("_exchange_point_var() is not implemented for use with\n"
                "MPI datatypes associated with structures using padding\n"
                "(for which size != extent)."));

  /* Initialization */

  n_points_loc_tot
    = this_locator->local_points_idx[this_locator->n_intersects];

  PLE_MALLOC(loc_v_flag, this_locator->n_intersects, int);
  PLE_MALLOC(dist_v_flag, this_locator->n_intersects, int);
  PLE_MALLOC(request, this_locator->n_intersects*2, MPI_Request);
  PLE_MALLOC(status, this_locator->n_intersects*2, MPI_Status);

  PLE_MALLOC(loc_v_buf, n_points_loc_tot*size*stride, char);

  /* First loop on distant ranks for argument checks */
  /*-------------------------------------------------*/

  _locator_trace_start_comm(_ple_locator_log_start_p_comm, comm_timing);

  for (i = 0; i < this_locator->n_intersects; i++) {

    dist_rank = this_locator->intersect_rank[i];

    n_points_dist =   this_locator->distant_points_idx[i+1]
                    - this_locator->distant_points_idx[i];

    if (distant_var != NULL && n_points_dist > 0)
      dist_v_flag[i] = 1;
    else
      dist_v_flag[i] = 0;

    _locator_trace_start_comm(_ple_locator_log_start_p_comm, comm_timing);

    MPI_Irecv(loc_v_flag + i, 1, MPI_INT, dist_rank, PLE_MPI_TAG,
              this_locator->comm, &request[i*2]);
    MPI_Isend(dist_v_flag + i, 1, MPI_INT, dist_rank, PLE_MPI_TAG,
              this_locator->comm, &request[i*2+1]);
  }

  MPI_Waitall(this_locator->n_intersects*2, request, status);

  _locator_trace_end_comm(_ple_locator_log_end_p_comm, comm_timing);

  PLE_FREE(dist_v_flag);

  for (i = 0; i < this_locator->n_intersects; i++) {

    dist_rank = this_locator->intersect_rank[i];

    n_points_loc =   this_locator->local_points_idx[i+1]
                   - this_locator->local_points_idx[i];

    if (loc_v_flag[i] == 1 && (local_var == NULL || n_points_loc == 0))
      ple_error(__FILE__, __LINE__, 0,
                _("Incoherent arguments to different instances in "
                  "_exchange_point_var().\n"
                  "Send and receive operations do not match "
                  "(dist_rank = %d\n)\n"), dist_rank);
  }

  /* Loop on distant ranks for exchange of data in standard mode */
  /*-------------------------------------------------------------*/

  if (reverse == false) {

    _locator_trace_start_comm(_ple_locator_log_start_p_comm, comm_timing);

    loc_v_ptr = loc_v_buf;

    for (i = 0; i < this_locator->n_intersects; i++) {

      dist_rank = this_locator->intersect_rank[i];

      n_points_loc =    this_locator->local_points_idx[i+1]
                      - this_locator->local_points_idx[i];

      n_points_dist =   this_locator->distant_points_idx[i+1]
                      - this_locator->distant_points_idx[i];

      dist_v_idx = this_locator->distant_points_idx[i] * stride*size;
      dist_v_count = n_points_dist * stride * dist_v_flag[i];

      if (loc_v_flag[i] > 0)
        loc_v_count = n_points_loc*stride;
      else
        loc_v_count = 0;

      /* Exchange information */

      if (distant_var != NULL)
        dist_v_ptr = ((unsigned char *)distant_var) + dist_v_idx;
      else
        dist_v_ptr = NULL;

      MPI_Irecv(loc_v_ptr, loc_v_count, datatype, dist_rank, PLE_MPI_TAG,
                this_locator->comm, &request[i*2]);
      MPI_Isend(dist_v_ptr, dist_v_count, datatype, dist_rank, PLE_MPI_TAG,
                this_locator->comm, &request[i*2+1]);

      loc_v_ptr += loc_v_count;
    }

    MPI_Waitall(this_locator->n_intersects*2, request, status);

    _locator_trace_end_comm(_ple_locator_log_end_p_comm, comm_timing);
  }

  /* Loop on distant ranks for preparation or retrieval of buffer data */
  /*-------------------------------------------------------------------*/

  loc_v_ptr = loc_v_buf;

  for (i = 0; i < this_locator->n_intersects; i++) {

    const ple_lnum_t *_local_point_ids
      = this_locator->local_point_ids + this_locator->local_points_idx[i];

    dist_rank = this_locator->intersect_rank[i];

    n_points_loc =    this_locator->local_points_idx[i+1]
      - this_locator->local_points_idx[i];

    n_points_dist =   this_locator->distant_points_idx[i+1]
                    - this_locator->distant_points_idx[i];

    dist_v_idx = this_locator->distant_points_idx[i] * stride*size;
    dist_v_count = n_points_dist * stride * dist_v_flag[i];

    if (loc_v_flag[i] > 0)
      loc_v_count = n_points_loc*stride;
    else
      loc_v_count = 0;

    /* Exchange information */

    if (distant_var != NULL)
      dist_v_ptr = ((unsigned char *)distant_var) + dist_v_idx;
    else
      dist_v_ptr = NULL;

    dist_rank = this_locator->intersect_rank[i];

    n_points_loc =    this_locator->local_points_idx[i+1]
                    - this_locator->local_points_idx[i];

    n_points_dist =   this_locator->distant_points_idx[i+1]
                    - this_locator->distant_points_idx[i];

    if (reverse == false) {

      if (loc_v_flag[i] > 0) {
        if (local_list == NULL) {
          int k;
          size_t l;
          const size_t nbytes = stride*size;
          for (k = 0; k < n_points_loc; k++) {
            char *local_v_p = (char *)local_var + _local_point_ids[k]*nbytes;
            const char *loc_v_buf_p = (const char *)loc_v_ptr + k*nbytes;
            for (l = 0; l < nbytes; l++)
              local_v_p[l] = loc_v_buf_p[l];
          }
        }
        else {
          int k;
          size_t l;
          const size_t nbytes = stride*size;
          for (k = 0; k < n_points_loc; k++) {
            char *local_v_p =   (char *)local_var
                              + (local_list[_local_point_ids[k]] - 1)*nbytes;
            const char *loc_v_buf_p = (const char *)loc_v_ptr + k*nbytes;
            for (l = 0; l < nbytes; l++)
              local_v_p[l] = loc_v_buf_p[l];
          }
        }
      }

    }
    else { /* if (reverse == true) */

      if (loc_v_flag[i] > 0) {
        if (local_list == NULL) {
          int k;
          size_t l;
          const size_t nbytes = stride*size;
          for (k = 0; k < n_points_loc; k++) {
            const char *local_v_p
              = (const char *)local_var + _local_point_ids[k]*nbytes;
            char *loc_v_buf_p = (char *)loc_v_ptr + k*nbytes;
            for (l = 0; l < nbytes; l++)
              loc_v_buf_p[l] = local_v_p[l];
          }
        }
        else {
          int k;
          size_t l;
          const size_t nbytes = stride*size;
          for (k = 0; k < n_points_loc; k++) {
            const char *local_v_p
              = (const char *)local_var
                + (local_list[_local_point_ids[k]] - 1)*nbytes;
            char *loc_v_buf_p = (char *)loc_v_ptr + k*nbytes;
            for (l = 0; l < nbytes; l++)
              loc_v_buf_p[l] = local_v_p[l];
          }
        }
      }

    }

    loc_v_ptr += loc_v_count;

  }

  /* Loop on distant ranks for exchange of data in reverse mode */
  /*------------------------------------------------------------*/

  if (reverse == true) {

    _locator_trace_start_comm(_ple_locator_log_start_p_comm, comm_timing);

    loc_v_ptr = loc_v_buf;

    for (i = 0; i < this_locator->n_intersects; i++) {

      dist_rank = this_locator->intersect_rank[i];

      n_points_loc =    this_locator->local_points_idx[i+1]
                      - this_locator->local_points_idx[i];

      n_points_dist =   this_locator->distant_points_idx[i+1]
                      - this_locator->distant_points_idx[i];

      dist_v_idx = this_locator->distant_points_idx[i] * stride*size;
      dist_v_count = n_points_dist * stride * dist_v_flag[i];

      if (loc_v_flag[i] > 0)
        loc_v_count = n_points_loc*stride;
      else
        loc_v_count = 0;

      /* Exchange information */

      if (distant_var != NULL)
        dist_v_ptr = ((unsigned char *)distant_var) + dist_v_idx;
      else
        dist_v_ptr = NULL;

      MPI_Irecv(dist_v_ptr, dist_v_count, datatype, dist_rank, PLE_MPI_TAG,
                this_locator->comm, &request[i*2]);
      MPI_Isend(loc_v_buf, loc_v_count, datatype, dist_rank, PLE_MPI_TAG,
                this_locator->comm, &request[i*2+1]);

      loc_v_ptr += loc_v_count;
    }

    MPI_Waitall(this_locator->n_intersects*2, request, status);

    _locator_trace_end_comm(_ple_locator_log_end_p_comm, comm_timing);

  }

  /* Free temporary arrays */

  PLE_FREE(loc_v_buf);

  PLE_FREE(loc_v_flag);
  PLE_FREE(request);
  PLE_FREE(status);

  this_locator->exchange_wtime[1] += comm_timing[0];
  this_locator->exchange_cpu_time[1] += comm_timing[1];
}

#endif /* defined(PLE_HAVE_MPI) */

/*----------------------------------------------------------------------------
 * Distribute variable defined on "distant points" to the original ("local")
 * points.
 *
 * The exchange is symmetric if both variables are defined, receive
 * only if distant_var is NULL, or send only if local_var is NULL.
 *
 * parameters:
 *   this_locator  <-- pointer to locator structure
 *   distant_var   <-> variable defined on distant points (ready to send)
 *   local_var     <-> variable defined on local points (received)
 *   local_list    <-- optional indirection list (1 to n) for local_var
 *   type_size     <-- sizeof (float or double) variable type
 *   stride        <-- dimension (1 for scalar, 3 for interlaced vector)
 *   reverse       <-- if true, exchange is reversed
 *                     (receive values associated with distant points
 *                     from the processes owning the original points)
 *----------------------------------------------------------------------------*/

static void
_exchange_point_var_local(ple_locator_t     *this_locator,
                          void              *distant_var,
                          void              *local_var,
                          const ple_lnum_t  *local_list,
                          size_t             type_size,
                          size_t             stride,
                          _Bool              reverse)
{
  ple_lnum_t i;
  size_t j;
  ple_lnum_t n_points_loc;

  const size_t nbytes = stride*type_size;

  /* Initialization */

  if (this_locator->n_interior == 0)
    return;

  n_points_loc =   this_locator->local_points_idx[1]
                 - this_locator->local_points_idx[0];

  assert(n_points_loc == (  this_locator->distant_points_idx[1]
                          - this_locator->distant_points_idx[0]));

  /* Exchange information */

  if (reverse == false) {

    if (local_list == NULL)
      memcpy(local_var, distant_var, n_points_loc*nbytes);

    else {
      for (i = 0; i < n_points_loc; i++) {
        char *local_var_p = (char *)local_var + (local_list[i] - 1)*nbytes;
        const char *distant_var_p = (const char *)distant_var + i*nbytes;
        for (j = 0; j < nbytes; j++)
          local_var_p[j] = distant_var_p[j];
      }
    }

  }
  else { /* if (reverse == true) */

    if (local_list == NULL)
      memcpy(distant_var, local_var, n_points_loc*nbytes);

    else {
      for (i = 0; i < n_points_loc; i++) {
        const char *local_var_p
          = (const char *)local_var + (local_list[i] - 1)*nbytes;
        char *distant_var_p = (char *)distant_var + i*nbytes;
        for (j = 0; j < nbytes; j++)
          distant_var_p[j] = local_var_p[j];
      }
    }

  }
}

/*----------------------------------------------------------------------------
 * Return timing information.
 *
 * When location on closest elements to force location of all points is
 * active, location times include a total value, followed by the value
 * associated with the location of closest elements stage.
 *
 * parameters:
 *   this_locator      <-- pointer to locator structure
 *   time_type         <-- 0 for total times, 1 for communication times
 *   location_wtime    --> Location Wall-clock time (or NULL)
 *   location_cpu_time --> Location CPU time (or NULL)
 *   exchange_wtime    --> Variable exchange Wall-clock time (or NULL)
 *   exchange_cpu_time --> Variable exchange CPU time (or NULL)
 *----------------------------------------------------------------------------*/

static void
_get_times(const ple_locator_t  *this_locator,
           int                   time_type,
           double               *location_wtime,
           double               *location_cpu_time,
           double               *exchange_wtime,
           double               *exchange_cpu_time)
{
  const ple_locator_t  *_locator = this_locator;

  if (this_locator != NULL) {

    if (location_wtime != NULL) {
      location_wtime[0] = _locator->location_wtime[time_type];
      location_wtime[1] = _locator->location_wtime[time_type + 2];
    }
    if (location_cpu_time != NULL) {
      location_cpu_time[0] = _locator->location_cpu_time[time_type];
      location_cpu_time[1] = _locator->location_cpu_time[time_type + 2];
    }
    if (exchange_wtime != NULL)
      *exchange_wtime = _locator->exchange_wtime[time_type];
    if (exchange_cpu_time != NULL)
      *exchange_cpu_time = _locator->exchange_cpu_time[time_type];

  }
  else {

    if (location_wtime != NULL) {
      location_wtime[0] = 0.;
      location_wtime[1] = 0.;
    }
    if (location_cpu_time != NULL) {
      location_cpu_time[0] = 0.;
      location_cpu_time[1] = 0.;
    }
    if (exchange_wtime != NULL)
      *exchange_wtime = 0.;
    if (exchange_cpu_time != NULL)
      *exchange_cpu_time = 0.;
  }
}

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Creation of a locator structure.
 *
 * Note that depending on the choice of ranks of the associated communicator,
 * distant ranks may in fact be truly distant or not. If n_ranks = 1 and
 * start_rank is equal to the current rank in the communicator, the locator
 * will work only locally.
 *
 * \param[in] tolerance  addition to local extents of each element:
 *                       extent = base_extent * (1 + tolerance)
 * \param[in] comm       associated MPI communicator
 * \param[in] n_ranks    number of MPI ranks associated with distant location
 * \param[in] start_rank first MPI rank associated with distant location
 *
 * \return pointer to locator
 */
/*----------------------------------------------------------------------------*/

#if defined(PLE_HAVE_MPI)
ple_locator_t *
ple_locator_create(double    tolerance,
                   MPI_Comm  comm,
                   int       n_ranks,
                   int       start_rank)
#else
ple_locator_t *
ple_locator_create(double  tolerance)
#endif
{
  int  i;
  ple_locator_t  *this_locator;

  PLE_MALLOC(this_locator, 1, ple_locator_t);

  this_locator->tolerance = tolerance;
  this_locator->dim = 0;

  this_locator->locate_closest = 0;

#if defined(PLE_HAVE_MPI)
  this_locator->comm = comm;
  this_locator->n_ranks = n_ranks;
  this_locator->start_rank = start_rank;
#else
  this_locator->n_ranks = 1;
  this_locator->start_rank = 0;
#endif

  this_locator->n_intersects = 0;
  this_locator->intersect_rank = NULL;
  this_locator->intersect_extents = NULL;

  this_locator->local_points_idx = NULL;
  this_locator->distant_points_idx = NULL;

  this_locator->local_point_ids = NULL;

  this_locator->distant_point_location = NULL;
  this_locator->distant_point_coords = NULL;

this_locator->distant_point_props = NULL; // 2016Jun13, create 

  this_locator->n_interior = 0;
  this_locator->interior_list = NULL;

  this_locator->n_exterior = 0;
  this_locator->exterior_list = NULL;

  for (i = 0; i < 4; i++) {
    this_locator->location_wtime[i] = 0.;
    this_locator->location_cpu_time[i] = 0.;
  }

  for (i = 0; i < 2; i++) {
    this_locator->exchange_wtime[i] = 0.;
    this_locator->exchange_cpu_time[i] = 0.;
  }

  return this_locator;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Destruction of a locator structure.
 *
 * \param[in, out] this_locator locator to destroy
 *
 * \return NULL pointer
 */
/*----------------------------------------------------------------------------*/

ple_locator_t *
ple_locator_destroy(ple_locator_t  *this_locator)
{
  if (this_locator != NULL) {

    PLE_FREE(this_locator->local_points_idx);
    PLE_FREE(this_locator->distant_points_idx);

    if (this_locator->local_point_ids != NULL)
      PLE_FREE(this_locator->local_point_ids);

    PLE_FREE(this_locator->distant_point_location);
    PLE_FREE(this_locator->distant_point_coords);

PLE_FREE(this_locator->distant_point_props); // 2016Jun13, destroy  

    PLE_FREE(this_locator->intersect_rank);
    PLE_FREE(this_locator->intersect_extents);

    PLE_FREE(this_locator->interior_list);
    PLE_FREE(this_locator->exterior_list);

    PLE_FREE(this_locator);
  }

  return NULL;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Prepare locator for use with a given mesh representation.
 *
 * \param[in, out] this_locator     pointer to locator structure
 * \param[in]      mesh             pointer to mesh representation structure
 * \param[in]      dim              spatial dimension of mesh and points to
 *                                  locate
 * \param[in]      n_points         number of points to locate
 * \param[in]      point_list       optional indirection array to point_coords
 *                                  (1 to n_points numbering)
 * \param[in]      point_coords     coordinates of points to locate
 *                                  (dimension: dim * n_points)
 * \param[out]     distance         optional distance from point to matching
 *                                  element: < 0 if unlocated; 0 - 1 if inside
 *                                  and > 1 if outside a volume element, or
 *                                  absolute distance to a surface element
 *                                  (size: n_points)
 * \param[in]      mesh_extents_f   pointer to function computing mesh or mesh
 *                                  subset or element extents
 * \param[in]      locate_inside_f  pointer to function wich updates the
 *                                  location[] and distance[] arrays associated
 *                                  with a set of points for points that are in
 *                                  an element of this mesh, or closer to one
 *                                  than to previously encountered elements.
 * \param[in]      locate_closest_f pointer to function locating the closest
 *                                  local elements if points not located on an
 *                                  element within the tolerance should be
 *                                  located on the closest element,
 *                                  NULL otherwise
 */
/*----------------------------------------------------------------------------*/

void
ple_locator_set_mesh(ple_locator_t                *this_locator,
                     const void                   *mesh,
                     int                           dim,
                     int                           dim_props,
                     ple_lnum_t                    n_points,
                     const ple_lnum_t              point_list[],
                     const ple_coord_t             point_coords[],
                     float                         distance[],
                     ple_mesh_extents_t           *mesh_extents_f,
                     ple_mesh_elements_contain_t  *locate_inside_f,
                     ple_mesh_elements_closest_t  *locate_closest_f, 
                     const ple_coord_t             point_props[]      // 2016Jun13 
                    )
{
  int i;
  int stride2;
  double tolerance;
  double w_start, w_end, cpu_start, cpu_end;
  double extents[12];

#if defined(PLE_HAVE_MPI)
  int j;
  int stride4;
  int comm_rank, comm_size;
  int n_intersects;
  int  *intersect_rank;
  double *recvbuf;
#endif

  double comm_timing[4] = {0., 0., 0., 0.};
  int mpi_flag = 0;

  /* Initialize timing */

  w_start = ple_timer_wtime();
  cpu_start = ple_timer_cpu_time();

  /* Other initializations */

  this_locator->locate_closest = 0;
  if (locate_closest_f != NULL)
    this_locator->locate_closest = 2;

  this_locator->dim = dim;

  this_locator->n_intersects = 0;

  tolerance = PLE_MAX(this_locator->tolerance, 1.e-3);

  mesh_extents_f(mesh,
                 1,
                 tolerance,
                 extents);

  _point_extents(dim,
                 n_points,
                 point_list,
                 point_coords,
                 extents + 2*dim);

  /* Release information if previously present */

  if (this_locator->intersect_rank != NULL)
    PLE_FREE(this_locator->intersect_rank);
  if (this_locator->intersect_extents != NULL)
    PLE_FREE(this_locator->intersect_extents);

  if (this_locator->local_points_idx != NULL)
    PLE_FREE(this_locator->local_points_idx);
  if (this_locator->distant_points_idx != NULL)
    PLE_FREE(this_locator->distant_points_idx);

  if (this_locator->local_point_ids != NULL)
    PLE_FREE(this_locator->local_point_ids);

  if (this_locator->distant_point_location != NULL)
    PLE_FREE(this_locator->distant_point_location);
  if (this_locator->distant_point_coords != NULL)
    PLE_FREE(this_locator->distant_point_coords);

if (this_locator->distant_point_props != NULL) PLE_FREE(this_locator->distant_point_props);  // 2016Jun13, DEALLOCATE  

  if (this_locator->interior_list != NULL)
    PLE_FREE(this_locator->interior_list);
  if (this_locator->exterior_list != NULL)
    PLE_FREE(this_locator->exterior_list);

  /* Prepare locator (MPI version) */
  /*-------------------------------*/

#if defined(PLE_HAVE_MPI)

  this_locator->async_exchange = 0;

  MPI_Initialized(&mpi_flag);

  if (mpi_flag && this_locator->comm == MPI_COMM_NULL)
    mpi_flag = 0;

  if (mpi_flag) {

    int globflag[3];
    int locflag[3] = {-1, -1, -1};

    /* Check that at least one of the local or distant nodal meshes
       is non-NULL, and at least one of the local or distant
       point sets is non null */

    if (mesh != NULL)
      locflag[0] = dim;

    if (n_points > 0) {
      locflag[1] = dim;
    }

    if (locate_closest_f != NULL)
      locflag[2] = 1;

    _locator_trace_start_comm(_ple_locator_log_start_g_comm, comm_timing);

    MPI_Allreduce(locflag, globflag, 3, MPI_INT, MPI_MAX,
                  this_locator->comm);

    _locator_trace_end_comm(_ple_locator_log_end_g_comm, comm_timing);

    if (globflag[0] < 0 || globflag[1] < 0)
      return;
    else if (mesh != NULL && globflag[1] != dim)
      ple_error(__FILE__, __LINE__, 0,
                _("Locator trying to use distant space dimension %d\n"
                  "with local space dimension %d\n"),
                globflag[1], dim);
    else if (mesh == NULL && globflag[0] != dim)
      ple_error(__FILE__, __LINE__, 0,
                _("Locator trying to use local space dimension %d\n"
                  "with distant space dimension %d\n"),
                dim, globflag[0]);

    if (this_locator->locate_closest == 0 && globflag[2] == 1)
      this_locator->locate_closest = 1;

    /* Exchange extent information */

    MPI_Comm_rank(this_locator->comm, &comm_rank);
    MPI_Comm_size(this_locator->comm, &comm_size);

    stride2 = dim * 2; /* Stride for one type of extent */
    stride4 = dim * 4; /* Stride for element and vertex
                          extents, end-to-end */

    PLE_MALLOC(recvbuf, stride4*comm_size, double);

    _locator_trace_start_comm(_ple_locator_log_start_g_comm, comm_timing);

    MPI_Allgather(extents, stride4, MPI_DOUBLE, recvbuf, stride4, MPI_DOUBLE,
                  this_locator->comm);

    _locator_trace_end_comm(_ple_locator_log_end_g_comm, comm_timing);

    /* Count and mark possible overlaps */

    n_intersects = 0;
    PLE_MALLOC(intersect_rank, this_locator->n_ranks, int);

    for (i = 0; i < this_locator->n_ranks; i++) {
      j = this_locator->start_rank + i;
      if (  (_intersect_extents(dim,
                                extents + (2*dim),
                                recvbuf + (j*stride4)) == true)
          || (_intersect_extents(dim,
                                 extents,
                                 recvbuf + (j*stride4) + (2*dim)) == true)) {
        intersect_rank[n_intersects] = j;
        n_intersects += 1;
      }

    }

    this_locator->n_intersects = n_intersects;
    PLE_MALLOC(this_locator->intersect_rank,
               this_locator->n_intersects,
               int);
    PLE_MALLOC(this_locator->intersect_extents,
               this_locator->n_intersects * stride2,
               double);

    for (i = 0; i < this_locator->n_intersects; i++) {

      this_locator->intersect_rank[i] = intersect_rank[i];

      /* Copy only distant element (and not point) extents */

      for (j = 0; j < stride2; j++)
        this_locator->intersect_extents[i*stride2 + j]
          = recvbuf[intersect_rank[i]*stride4 + j];

    }

    /* Free temporary memory */

    PLE_FREE(intersect_rank);
    PLE_FREE(recvbuf);

    _locate_all_distant(this_locator,
                        mesh,
                        dim,
                        dim_props,
                        n_points,
                        point_list,
                        point_coords,
                        distance,
                        mesh_extents_f,
                        locate_inside_f,
                        locate_closest_f, point_props // 2016Jun13  
                       );

  }

#endif

  /* Prepare locator (local version) */
  /*---------------------------------*/

  if (!mpi_flag) {

    if (mesh == NULL || n_points == 0)
      return;

    stride2 = dim * 2;

    /* Count and mark possible overlaps */

    if (_intersect_extents(dim,
                           extents,
                           extents + (2*dim)) == true) {

      this_locator->n_intersects = 1;

      PLE_MALLOC(this_locator->intersect_rank, 1, int);
      PLE_MALLOC(this_locator->intersect_extents, stride2, double);

      this_locator->intersect_rank[0] = 0;

      for (i = 0; i < stride2; i++)
        this_locator->intersect_extents[i] = extents[i];

      _locate_all_local(this_locator,
                        mesh,
                        dim,
                        dim_props,
                        n_points,
                        point_list,
                        point_coords,
                        distance,
                        locate_inside_f,
                        locate_closest_f, point_props
                       );
    }

    else {

      this_locator->n_exterior = n_points;
      PLE_MALLOC(this_locator->exterior_list,
                 this_locator->n_exterior,
                 ple_lnum_t);
      for (i = 0; i < this_locator->n_exterior; i++)
        this_locator->exterior_list[i] = i + 1;

    }

  }

  /* Update local_point_ids values */
  /*-------------------------------*/

  if (   this_locator->n_interior > 0
      && this_locator->local_point_ids != NULL) {

    ple_lnum_t  *reduced_index;

    PLE_MALLOC(reduced_index, n_points, ple_lnum_t);

    for (i = 0; i < n_points; i++)
      reduced_index[i] = -1;

    assert(  this_locator->local_points_idx[this_locator->n_intersects]
           == this_locator->n_interior);

    for (i = 0; i < this_locator->n_interior; i++)
      reduced_index[this_locator->interior_list[i] - 1] = i;

    /* Update this_locator->local_point_ids[] so that it refers
       to an index in a dense [0, this_locator->n_interior] subset
       of the local points */

    for (i = 0; i < this_locator->n_interior; i++)
      this_locator->local_point_ids[i]
        = reduced_index[this_locator->local_point_ids[i]];

    PLE_FREE(reduced_index);

  }

  /* If an initial point list was given, update
     this_locator->interior_list and this_locator->exterior_list
     so that they refer to the same point set as that initial
     list (and not to an index within the selected point set) */

  if (point_list != NULL) {

    for (i = 0; i < this_locator->n_interior; i++)
      this_locator->interior_list[i]
        = point_list[this_locator->interior_list[i] - 1];

    for (i = 0; i < this_locator->n_exterior; i++)
      this_locator->exterior_list[i]
        = point_list[this_locator->exterior_list[i] - 1];

  }

  /* Finalize timing */

  w_end = ple_timer_wtime();
  cpu_end = ple_timer_cpu_time();

  this_locator->location_wtime[0] += (w_end - w_start);
  this_locator->location_cpu_time[0] += (cpu_end - cpu_start);

  this_locator->location_wtime[1] += comm_timing[0];
  this_locator->location_cpu_time[1] += comm_timing[1];
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Return number of distant points after locator initialization.
 *
 * \param[in] this_locator pointer to locator structure
 *
 * \return number of distant points.
 */
/*----------------------------------------------------------------------------*/

ple_lnum_t
ple_locator_get_n_dist_points(const ple_locator_t  *this_locator)
{
  ple_lnum_t retval = 0;

  if (this_locator != NULL) {
    if (this_locator->n_intersects != 0)
      retval = this_locator->distant_points_idx[this_locator->n_intersects];
  }

  return retval;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Return an array of local element numbers containing (or nearest to)
 *  each distant point after locator initialization.
 *
 * \param[in] this_locator pointer to locator structure
 *
 * \return local element numbers associated with distant points
 *        (1 to n numbering).
 */
/*----------------------------------------------------------------------------*/

const ple_lnum_t *
ple_locator_get_dist_locations(const ple_locator_t  *this_locator)
{
  const ple_lnum_t * retval = NULL;

  if (this_locator != NULL) {
    if (this_locator->n_ranks != 0)
      retval = this_locator->distant_point_location;
  }

  return retval;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Return an array of coordinates of each distant point after
 * locator initialization.
 *
 * \param[in] this_locator pointer to locator structure
 *
 * \return coordinate array associated with distant points (interlaced).
 */
/*----------------------------------------------------------------------------*/

const ple_coord_t *
ple_locator_get_dist_coords(const ple_locator_t  *this_locator)
{
  const ple_coord_t * retval = NULL;

  if (this_locator != NULL) {
    if (this_locator->n_intersects != 0)
      retval = this_locator->distant_point_coords;
  }

  return retval;
}

const ple_coord_t *
ple_locator_get_dist_props(const ple_locator_t  *this_locator)
{
  const ple_coord_t * retval = NULL;

  if (this_locator != NULL) {
    if (this_locator->n_intersects != 0)
      retval = this_locator->distant_point_props;                // 2016Jun13  
  }

  return retval;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Return number of points located after locator initialization.
 *
 * \param[in] this_locator pointer to locator structure
 *
 * \return number of points located.
 */
/*----------------------------------------------------------------------------*/

ple_lnum_t
ple_locator_get_n_interior(const ple_locator_t  *this_locator)
{
  ple_lnum_t retval = 0;

  if (this_locator != NULL)
    retval = this_locator->n_interior;

  return retval;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Return list of points located after locator initialization.
 * This list defines a subset of the point set used at initialization.
 *
 * \param[in] this_locator pointer to locator structure
 *
 * \return list of points located (1 to n numbering).
 */
/*----------------------------------------------------------------------------*/

const ple_lnum_t *
ple_locator_get_interior_list(const ple_locator_t  *this_locator)
{
  return this_locator->interior_list;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Return number of points not located after locator initialization.
 *
 * \param[in] this_locator pointer to locator structure
 *
 * \return number of points not located.
 */
/*----------------------------------------------------------------------------*/

ple_lnum_t
ple_locator_get_n_exterior(const ple_locator_t  *this_locator)
{
  return this_locator->n_exterior;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Return list of points not located after locator initialization.
 * This list defines a subset of the point set used at initialization.
 *
 * \param[in] this_locator pointer to locator structure
 *
 * \return list of points not located (1 to n numbering).
 */
/*----------------------------------------------------------------------------*/

const ple_lnum_t *
ple_locator_get_exterior_list(const ple_locator_t  *this_locator)
{
  return this_locator->exterior_list;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Discard list of points not located after locator initialization.
 * This list defines a subset of the point set used at initialization.
 *
 * \param[in] this_locator pointer to locator structure
 */
/*----------------------------------------------------------------------------*/

void
ple_locator_discard_exterior(ple_locator_t  *this_locator)
{
  this_locator->n_exterior = 0;
  PLE_FREE(this_locator->exterior_list);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Distribute variable defined on distant points to processes owning
 * the original points (i.e. distant processes).
 *
 * The exchange is symmetric if both variables are defined, receive
 * only if distant_var is NULL, or send only if local_var is NULL.
 *
 * The caller should have defined the values of distant_var[] for the
 * distant points, whose coordinates are given by
 * ple_locator_get_dist_coords(), and which are located in the elements
 * whose numbers are given by ple_locator_get_dist_locations().
 *
 * The local_var[] is defined at the located points (those whose
 * numbers are returned by ple_locator_get_interior_list().
 *
 * \param[in]      this_locator pointer to locator structure
 * \param[in, out] distant_var  variable defined on distant points
 *                              (ready to send); size: n_dist_points*stride
 * \param[in, out] local_var    variable defined on located local points
 *                              (received); size: n_interior*stride
 * \param[in]      local_list   optional indirection list (1 to n) for local_var
 * \param[in]      type_size    sizeof (float or double) variable type
 * \param[in]      stride       dimension (1 for scalar,
 *                              3 for interleaved vector)
 * \param[in]      reverse      if nonzero, exchange is reversed
 *                              (receive values associated with distant points
 *                              from the processes owning the original points)
 */
/*----------------------------------------------------------------------------*/

void
ple_locator_exchange_point_var(ple_locator_t     *this_locator,
                               void              *distant_var,
                               void              *local_var,
                               const ple_lnum_t  *local_list,
                               size_t             type_size,
                               size_t             stride,
                               int                reverse)
{
  double w_start, w_end, cpu_start, cpu_end;

  int mpi_flag = 0;
  _Bool _reverse = reverse;

  /* Initialize timing */

  w_start = ple_timer_wtime();
  cpu_start = ple_timer_cpu_time();

#if defined(PLE_HAVE_MPI)

  MPI_Initialized(&mpi_flag);

  if (mpi_flag && this_locator->comm == MPI_COMM_NULL)
    mpi_flag = 0;

  if (mpi_flag) {

    MPI_Datatype datatype = MPI_DATATYPE_NULL;

    if (type_size == sizeof(double))
      datatype = MPI_DOUBLE;
    else if (type_size == sizeof(float))
      datatype = MPI_FLOAT;
    else
      ple_error(__FILE__, __LINE__, 0,
                _("type_size passed to ple_locator_exchange_point_var() does\n"
                  "not correspond to double or float."));

    assert (datatype != MPI_DATATYPE_NULL);

    if (this_locator->async_exchange == 0)
      _exchange_point_var_distant(this_locator,
                                  distant_var,
                                  local_var,
                                  local_list,
                                  datatype,
                                  stride,
                                  _reverse);

    else
      _exchange_point_var_distant_asyn(this_locator,
                                       distant_var,
                                       local_var,
                                       local_list,
                                       datatype,
                                       stride,
                                       _reverse);

  }

#endif /* defined(PLE_HAVE_MPI) */

  if (!mpi_flag)
    _exchange_point_var_local(this_locator,
                              distant_var,
                              local_var,
                              local_list,
                              type_size,
                              stride,
                              _reverse);

  /* Finalize timing */

  w_end = ple_timer_wtime();
  cpu_end = ple_timer_cpu_time();

  this_locator->exchange_wtime[0] += (w_end - w_start);
  this_locator->exchange_cpu_time[0] += (cpu_end - cpu_start);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Return timing information.
 *
 * In parallel mode, this includes communication time.
 *
 * When location on closest elements to force location of all points is
 * active, location times include a total value, followed by the value
 * associated with the location of closest elements stage.
 *
 * \param[in]  this_locator      pointer to locator structure
 * \param[out] location_wtime    Location Wall-clock time (size: 2 or NULL)
 * \param[out] location_cpu_time Location CPU time (size: 1 or NULL)
 * \param[out] exchange_wtime    Variable exchange Wall-clock time
 *                               (size: 1 or NULL)
 * \param[out] exchange_cpu_time Variable exchange CPU time (size: 2 or NULL)
 */
/*----------------------------------------------------------------------------*/

void
ple_locator_get_times(const ple_locator_t  *this_locator,
                      double               *location_wtime,
                      double               *location_cpu_time,
                      double               *exchange_wtime,
                      double               *exchange_cpu_time)
{
  _get_times(this_locator,
             0,
             location_wtime, location_cpu_time,
             exchange_wtime, exchange_cpu_time);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Return communication timing information.
 *
 * In serial mode, return times are always zero.
 *
 * When location on closest elements to force location of all points is
 * active, location times include a total value, followed by the value
 * associated with the location of closest elements stage.
 *
 * parameters:
 * \param[in]  this_locator      pointer to locator structure
 * \param[out] location_wtime    Location Wall-clock time (size: 2 or NULL)
 * \param[out] location_cpu_time Location CPU time (size: 1 or NULL)
 * \param[out] exchange_wtime    Variable exchange Wall-clock time
 *                               (size: 1 or NULL)
 * \param[out] exchange_cpu_time Variable exchange CPU time (size: 2 or NULL)
 */
/*----------------------------------------------------------------------------*/

void
ple_locator_get_comm_times(const ple_locator_t  *this_locator,
                           double               *location_wtime,
                           double               *location_cpu_time,
                           double               *exchange_wtime,
                           double               *exchange_cpu_time)
{
  _get_times(this_locator,
             1,
             location_wtime, location_cpu_time,
             exchange_wtime, exchange_cpu_time);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Dump printout of a locator structure.
 *
 * \param this_locator pointer to structure that should be dumped
 */
/*----------------------------------------------------------------------------*/

void
ple_locator_dump(const ple_locator_t  *this_locator)
{
  int  i;
  ple_lnum_t  j, k;
  const ple_lnum_t  *idx, *index, *loc;
  const ple_coord_t  *coords;

  const ple_locator_t  *_locator = this_locator;

  if (this_locator == NULL)
    return;

  /* Basic information */
  /*-------------------*/

  ple_printf("\n"
             "Locator:\n\n"
             "Tolerance:                             %f\n"
             "Spatial dimension:                     %d\n"
             "Locate on closest:                     %d\n"
             "Number of ranks of distant location:   %d\n"
             "First rank of distant location:        %d\n"
             "Number of intersecting distant ranks:  %d\n",
             _locator->tolerance, _locator->dim,
             (int)_locator->locate_closest,
             _locator->n_ranks, _locator->start_rank,
             _locator->n_intersects);

#if defined(PLE_HAVE_MPI)
  if (_locator->comm != MPI_COMM_NULL)
    ple_printf("Asynchronous exchange:                 %d\n"
               "\n"
               "Associated MPI communicator:           %ld\n",
               _locator->async_exchange,
               (long)(_locator->comm));
#endif

  /* Arrays indexed by rank */
  /*------------------------*/

  for (i = 0; i < _locator->n_intersects; i++) {

    ple_printf("\n"
               "  Intersection %d with distant rank %d\n\n",
               i+1, _locator->intersect_rank[i]);

    ple_printf("    Distant rank extents:\n");

    k = i * (_locator->dim) * 2;
    for (j = 0; j < _locator->dim; j++)
      ple_printf("    [%12.5e, %12.5e]\n",
                 _locator->intersect_extents[k + j],
                 _locator->intersect_extents[k + _locator->dim + j]);

  }

  if (_locator->n_interior > 0) {

    if (_locator->local_point_ids != NULL) {

      ple_printf("\n  Local point ids (for receiving):\n\n");
      idx = _locator->local_points_idx;
      index = _locator->local_point_ids;
      for (i = 0; i < _locator->n_intersects; i++) {
        if (idx[i+1] > idx[i]) {
          ple_printf("%6d (idx = %10d) %10d\n",
                     i + 1, idx[i], index[idx[i]]);
          for (j = idx[i] + 1; j < idx[i + 1]; j++)
            ple_printf("                          %10d\n", index[j]);
        }
        else {
          ple_printf("%6d (idx = %10d)\n", i + 1, idx[i]);
        }
        ple_printf("   end (idx = %10d)\n", idx[_locator->n_intersects]);
      }

    }

  }

  if (_locator->distant_points_idx != NULL) {

    idx = _locator->distant_points_idx;
    loc = _locator->distant_point_location;
    coords = _locator->distant_point_coords;

    if (idx[_locator->n_intersects] > 0)
      ple_printf("\n  Distant point location:\n\n");

    for (i = 0; i < _locator->n_intersects; i++) {

      if (idx[i+1] > idx[i]) {

        if (_locator->dim == 1) {
          ple_printf("%6d (idx = %10d) %10d [%12.5e]\n",
                     i + 1, _locator->intersect_rank[i], idx[i],
                     loc[idx[i]], coords[idx[i]]);
          for (j = idx[i] + 1; j < idx[i + 1]; j++)
            ple_printf("                          %10d [%12.5e]\n",
                       loc[j], coords[j]);
        }
        else if (_locator->dim == 2) {
          ple_printf("%6d (idx = %10d) %10d [%12.5e, %12.5e]\n",
                     i + 1, idx[i], loc[idx[i]],
                     coords[2*idx[i]], coords[2*idx[i]+1]);
          for (j = idx[i] + 1; j < idx[i + 1]; j++)
            ple_printf("                          %10d [%12.5e, %12.5e]\n",
                       loc[j], coords[2*j], coords[2*j+1]);
        }
        else if (_locator->dim == 3) {
          ple_printf("%6d (idx = %10d) %10d [%12.5e, %12.5e, %12.5e]\n",
                     i + 1, idx[i], loc[idx[i]],
                     coords[3*idx[i]], coords[3*idx[i]+1], coords[3*idx[i]+2]);
          for (j = idx[i] + 1; j < idx[i + 1]; j++)
            ple_printf("                          "
                       "%10d [%12.5e, %12.5e, %12.5e]\n",
                       loc[j], coords[3*j], coords[3*j+1], coords[3*j+2]);
        }

      } /* if (idx[i+1] > idx[i]) */

    }

    if (idx[_locator->n_intersects] > 0)
      ple_printf("   end (idx = %10d)\n", idx[_locator->n_intersects]);
  }

  /* Local arrays */
  /*--------------*/

  ple_printf("\n"
             "  Number of local points successfully located:  %d\n\n",
             _locator->n_interior);

  for (j = 0; j < _locator->n_interior; j++)
    ple_printf("    %10d\n", _locator->interior_list[j]);

  if  (_locator->n_interior > 0)
    ple_printf("\n");

  ple_printf("  Number of local points not located:  %d\n",
             _locator->n_exterior);

  for (j = 0; j < _locator->n_exterior; j++)
    ple_printf("    %10d\n", _locator->exterior_list[j]);

  if  (_locator->n_exterior > 0)
    ple_printf("\n");

  /* Timing information */
  /*--------------------*/

  ple_printf("  Location Wall-clock time: %12.5f (comm: %12.5f)\n",
             _locator->location_wtime[0], _locator->location_wtime[0]);

  ple_printf("  Location CPU time:        %12.5f (comm: %12.5f)\n",
             _locator->location_cpu_time[0], _locator->location_cpu_time[0]);

  ple_printf("  Exchange Wall-clock time: %12.5f (comm: %12.5f)\n",
             _locator->exchange_wtime[0], _locator->exchange_wtime[0]);

  ple_printf("  Exchange CPU time:        %12.5f (comm: %12.5f)\n",
             _locator->exchange_cpu_time[0], _locator->exchange_cpu_time[0]);

}

#if defined(PLE_HAVE_MPI)

/*----------------------------------------------------------------------------*/
/*!
 * \brief Get the maximum number of exchanging ranks for which we use
 * asynchronous MPI sends and receives instead of MPI_SendRecv.
 *
 * \return the maximum number of ranks allowing asynchronous exchanges
 */
/*----------------------------------------------------------------------------*/

int
ple_locator_get_async_threshold(void)
{
  return _ple_locator_async_threshold;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Set the maximum number of exchanging ranks for which we use
 * asynchronous MPI sends and receives instead of MPI_SendRecv.
 *
 * \param threshold maximum number of ranks allowing asynchronous exchanges
 */
/*----------------------------------------------------------------------------*/

void
ple_locator_set_async_threshold(int threshold)
{
  _ple_locator_async_threshold = threshold;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Register communication logging functions for locator instrumentation.
 *
 * By default, locators are not instrumented.
 *
 * Functions using MPE may be defined and used, but other similar systems
 * may be used.
 *
 * \param[in] fct          pointer to logging function
 * \param[in] start_p_comm point to point communication start event number
 * \param[in] end_p_comm   point to point communication end event number
 * \param[in] start_g_comm global communication start event number
 * \param[in] end_g_comm   global communication end event number
 */
/*----------------------------------------------------------------------------*/

void
ple_locator_set_comm_log(ple_locator_log_t  *log_function,
                         int                 start_p_comm,
                         int                 end_p_comm,
                         int                 start_g_comm,
                         int                 end_g_comm)
{
  _ple_locator_log_func = log_function;

  _ple_locator_log_start_p_comm = start_p_comm;
  _ple_locator_log_end_p_comm = end_p_comm;
  _ple_locator_log_start_g_comm = start_g_comm;
  _ple_locator_log_end_g_comm = end_g_comm;
}

#endif /* defined(PLE_HAVE_MPI) */

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
}
#endif /* __cplusplus */
