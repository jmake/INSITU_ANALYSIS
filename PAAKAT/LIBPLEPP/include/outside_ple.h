/*
* From: syrthes-kernel/src/syr_cfd_point_location.c
*/

#include "assert.h"
#include "ple_defs.h"
#include "ple_locator.h"
#include "ple_coupling.h"
//=======================================================================||===//
//=======================================================================||===//
#ifndef OUTSIDE_PLE_H
#define OUTSIDE_PLE_H


//=================================================================| MESH |===//
/*
  SATURNE: /src/fvm/fvm_nodal_priv._fvm_nodal_t 
  SYRTHES: /src/syrthes-kernel/include/syr_cfd_mesh_priv._syr_cfd_mesh_t 
  Structure defining a mesh
*/

const int  syr_cfd_mesh_n_vertices_element[] = {2,   // Edge 
                                                3,   // Triangle 
                                                4,   // Tetrahedron  
                                                8,   // Hexahedron. 2019FEB21  
                                               };
typedef enum
{
  SYR_CFD_EDGE,            /* Edge */
  SYR_CFD_TRIA,            /* Triangle */
  SYR_CFD_TETRA,           /* Tetrahedron */
  FVM_CELL_HEXA,           // Hexahedron (brick). 2019FEB21  
  SYR_CFD_N_ELEMENT_TYPES  /* Number of element types */
} syr_cfd_element_t;



typedef struct 
{
  int                 dim;           // Spatial dimension 
  syr_cfd_element_t   element_type;  // Element type 
  ple_lnum_t          n_vertices;    // Number of vertices 
  ple_lnum_t          n_elements;    // Number of elements 
  ple_coord_t        *vertex_coords; // (x1, y1, z1, x2, y2, z2, ...)  
  ple_lnum_t         *vertex_num;    // element vertex numbers (1 to n); size: n_elements*(element_dim+1) 
} syr_cfd_mesh_t; //, fvm_nodal_t; 

//=======================================================================||===//
//=======================================================================||===//
/*
  SATURNE: ./src/fvm/fvm_defs.h

typedef enum 
{
  FVM_EDGE,               // Edge 
  FVM_FACE_TRIA,          // Triangle 
  FVM_FACE_QUAD,          // Quadrangle 
  FVM_FACE_POLY,          // Simple Polygon 
  FVM_CELL_TETRA,         // Tetrahedron 
  FVM_CELL_PYRAM,         // Pyramid 
  FVM_CELL_PRISM,         // Prism (pentahedron) 
  FVM_CELL_HEXA,          // Hexahedron (brick) 
  FVM_CELL_POLY,          // Simple Polyhedron (convex or quasi-convex)
  FVM_N_ELEMENT_TYPES     // Number of element types
} fvm_element_t;
*/

/*
  SATURNE: ./src/fvm/fvm_nodal.c

const int  fvm_nodal_n_vertices_element[] = {2,   // Edge
                                             3,   // Triangle
                                             4,   // Quadrangle 
                                             0,   // Simple polygon 
                                             4,   // Tetrahedron 
                                             5,   // Pyramid 
                                             6,   // Prism 
                                             8,   // Hexahedron 
                                             0};  // Simple polyhedron 

*/
//=======================================================================||===//
//=======================================================================||===//


//=======================================================================||===//
//=======================================================================||===//

/* 
PLE_COUPLING_NO_SYNC; PLE_COUPLING_TS_MIN; PLE_COUPLING_TS_LEADER; PLE_COUPLING_UNSTEADY 
*/
static int _cs_coupling_sync_flag = PLE_COUPLING_TS_MIN; 

/*
const char *sync_name[2] = {"point-to-point or not synchronized",
                            "group synchronized"};
*/

#define CS_MIN(a,b)   ((a) < (b) ?  (a) : (b))  /* Minimum of a et b */

//=======================================================================||===//
//=======================================================================||===//


//=======================================================================||===//
//=======================================================================||===//
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

void
_locate_in_cell_3d(                                                // 2019Mar06  
                   syr_cfd_element_t   elt_type,
                   const ple_lnum_t    element_vertex_num[],
                   const ple_coord_t   vertex_coords[],
                   const ple_coord_t   point[],
                   double              tolerance,
                   double*             shapef
                   );

#ifdef __cplusplus
}
#endif
//=======================================================================||===//
//=======================================================================||===//


#endif
