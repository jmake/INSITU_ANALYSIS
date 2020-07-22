#include "outside_ple.h"
//=======================================================================||===//
//=======================================================================||===//
#ifndef BUILD_OCTREE_INL
#define BUILD_OCTREE_INL


typedef struct 
{
  ple_lnum_t  octant_id[8];   // Ids of sub-octants in octree array 
  ple_lnum_t  idx[9];         // Start index of point list for each octant 
  ple_lnum_t  n_points;       // Number of points in octree 

} _octant_t;


typedef struct 
{
  size_t       n_points;      // Number of points in octree 
  size_t       n_nodes;       // Current number of nodes in octree 
  size_t       n_nodes_max;   // Maximum number of nodes in octree 

  double       extents[6];    // Associated extents 

  ple_lnum_t  *point_ids;     // Id's of points sorted by octree (size: n_points + 1) 
                                 
  _octant_t   *nodes;         // Array of octree nodes
                              // (size: n_nodes_max) 

} _octree_t;


_octree_t
build_octree(ple_lnum_t n_points, const ple_coord_t  point_coords[]); 


#if !defined(HUGE_VAL)
  #define HUGE_VAL 1.0e+30
#endif
/*
#ifndef __cplusplus
 typedef unsigned char bool;
 static const bool false = 0;
 static const bool true = 1;
 //#include <stdbool.h>
#endif

#ifdef __cplusplus
extern "C" {  // c/fortran visibility!!
#endif

void
syr_cfd_point_location_contain(const void         *mesh,
                               double              tolerance,
                               ple_lnum_t          n_points,
                               const ple_coord_t   point_coords[],
                               ple_lnum_t          location[],
                               float               distance[]); 

ple_lnum_t
syr_cfd_point_location_extents(const void  *mesh,
                               ple_lnum_t   n_max_extents,
                               double       tolerance,
                               double       extents[]); 

void
syr_cfd_point_location_closest(const void         *mesh,
                               ple_lnum_t          n_points,
                               const ple_coord_t   point_coords[],
                               ple_lnum_t          location[],
                               float               distance[]);

#ifdef __cplusplus
}
#endif
*/
#endif 
//=======================================================================||===//
//=======================================================================||===//
