#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "outside_ple.h"

//#include "locate_in_hexa.c" // 2019FEB21   

/*
  /home/bsc21/bsc21704/EDF/SYRTHES/syrthes4.0.1_impi/src/syrthes-kernel/src
*/


//
//   elt_type    <-- type of element
//   uvw[]       <-- parametric coordinates
//   shapef[]    --> barycenter's coordinates
//   deriv [][]  --> derivative of shape function
//
static void
_compute_shapef_3d(syr_cfd_element_t elt_type,
                   const double   uvw[3],
                   double         shapef[8],
                   double         deriv[8][3])

{
  assert(elt_type == FVM_CELL_HEXA);

  memset(shapef, 0.0, 8*sizeof(double));

  shapef[0] = (1.0 - uvw[0]) * (1.0 - uvw[1]) * (1.0 - uvw[2]);
  shapef[1] = uvw[0] * (1.0 - uvw[1]) * (1.0 - uvw[2]);
  shapef[2] = uvw[0] * uvw[1] * (1.0 - uvw[2]);
  shapef[3] = (1.0 - uvw[0]) * uvw[1] * (1.0 - uvw[2]);
  shapef[4] = (1.0 - uvw[0]) * (1.0 - uvw[1]) * uvw[2];
  shapef[5] = uvw[0] * (1.0 - uvw[1]) * uvw[2];
  shapef[6] = uvw[0] * uvw[1] * uvw[2];
  shapef[7] = (1.0 - uvw[0]) * uvw[1] * uvw[2];

  if (deriv != NULL)
  {
      deriv[0][0] = -(1.0 - uvw[1]) * (1.0 - uvw[2]);
      deriv[0][1] = -(1.0 - uvw[0]) * (1.0 - uvw[2]);
      deriv[0][2] = -(1.0 - uvw[0]) * (1.0 - uvw[1]);
      deriv[1][0] =  (1.0 - uvw[1]) * (1.0 - uvw[2]);
      deriv[1][1] = -uvw[0] * (1.0 - uvw[2]);
      deriv[1][2] = -uvw[0] * (1.0 - uvw[1]);
      deriv[2][0] =  uvw[1] * (1.0 - uvw[2]);
      deriv[2][1] =  uvw[0] * (1.0 - uvw[2]);
      deriv[2][2] = -uvw[0] * uvw[1];
      deriv[3][0] = -uvw[1] * (1.0 - uvw[2]);
      deriv[3][1] =  (1.0 - uvw[0]) * (1.0 - uvw[2]);
      deriv[3][2] = -(1.0 - uvw[0]) * uvw[1];
      deriv[4][0] = -(1.0 - uvw[1]) * uvw[2];
      deriv[4][1] = -(1.0 - uvw[0]) * uvw[2];
      deriv[4][2] =  (1.0 - uvw[0]) * (1.0 - uvw[1]);
      deriv[5][0] =  (1.0 - uvw[1]) * uvw[2];
      deriv[5][1] = -uvw[0] * uvw[2];
      deriv[5][2] =  uvw[0] * (1.0 - uvw[1]);
      deriv[6][0] =  uvw[1] * uvw[2];
      deriv[6][1] =  uvw[0] * uvw[2];
      deriv[6][2] =  uvw[0] * uvw[1];
      deriv[7][0] = -uvw[1] * uvw[2];
      deriv[7][1] =  (1.0 - uvw[0]) * uvw[2];
      deriv[7][2] =  (1.0 - uvw[0]) * uvw[1];
  }

} // _compute_shapef_3d

//-----------------------------------------------------------------------||---//
static int
_inverse_3x3(double  m[3][3],
             double  b[3],
             double  x[3])
{
  double det, det_inv, x0, x1, x2;
  double _epsilon_denom = 1.e-28;

  det =   m[0][0]*(m[1][1]*m[2][2] - m[2][1]*m[1][2])
        - m[1][0]*(m[0][1]*m[2][2] - m[2][1]*m[0][2])
        + m[2][0]*(m[0][1]*m[1][2] - m[1][1]*m[0][2]);

  if(fabs(det) < _epsilon_denom)
  {
    return 1;
  }
  else
    det_inv = 1./det;

  // Use local variables to ensure no aliasing
  x0 = (  b[0]*(m[1][1]*m[2][2] - m[2][1]*m[1][2])
        - b[1]*(m[0][1]*m[2][2] - m[2][1]*m[0][2])
        + b[2]*(m[0][1]*m[1][2] - m[1][1]*m[0][2])) * det_inv;

  x1 = (  m[0][0]*(b[1]*m[2][2] - b[2]*m[1][2])
        - m[1][0]*(b[0]*m[2][2] - b[2]*m[0][2])
        + m[2][0]*(b[0]*m[1][2] - b[1]*m[0][2])) * det_inv;

  x2 = (  m[0][0]*(m[1][1]*b[2] - m[2][1]*b[1])
        - m[1][0]*(m[0][1]*b[2] - m[2][1]*b[0])
        + m[2][0]*(m[0][1]*b[1] - m[1][1]*b[0])) * det_inv;

  // Copy local variables to output

  x[0] = x0; x[1] = x1; x[2] = x2;

  return 0;
}

/*----------------------------------------------------------------------------
 * Calculation of natural coordinates of the element
 *
 *    elt_type            <-- type of element
 *    point_coords        <-- point coordinates
 *    vertex_coords[]     <-- pointer to element vertex coordinates
 *    tolerance           <-- location tolerance factor
 *    uvw[]               --> parametric coordinates of point in element
 *
 * returns.
 *    uvw[]
 *----------------------------------------------------------------------------*/
static int
_compute_uvw(syr_cfd_element_t  elt_type,
             const ple_coord_t  point_coords[],
             double             vertex_coords[8][3],
             double             tolerance,
             double             uvw[3])

{
  int i, j, n_elt_vertices, iter;
  int max_iter = 20;
  double dist;
  double a[3][3], b[3], x[3], shapef[8], dw[8][3];

  n_elt_vertices = syr_cfd_mesh_n_vertices_element[elt_type];

  assert(elt_type == FVM_CELL_HEXA); // || elt_type == FVM_CELL_PRISM || elt_type == FVM_CELL_PYRAM);

  // Use Newton-method to determine parametric coordinates and shape function
  for (i = 0; i < 3; i++) uvw[i] = 0.5;

  for (iter = 0; iter < max_iter; iter++)
  {
    _compute_shapef_3d(elt_type, uvw, shapef, dw);

    b[0] = - point_coords[0];
    b[1] = - point_coords[1];
    b[2] = - point_coords[2];

    for (i = 0; i < 3; i++) for (j = 0; j < 3; j++) a[i][j] = 0.0;

    for (i = 0; i < n_elt_vertices; i++)
    {
      b[0] += (shapef[i] * vertex_coords[i][0]);
      b[1] += (shapef[i] * vertex_coords[i][1]);
      b[2] += (shapef[i] * vertex_coords[i][2]);

      for (j = 0; j < 3; j++)
      {
        a[0][j] -= (dw[i][j] * vertex_coords[i][0]);
        a[1][j] -= (dw[i][j] * vertex_coords[i][1]);
        a[2][j] -= (dw[i][j] * vertex_coords[i][2]);
      }

    }

    if (_inverse_3x3(a, b, x)) return 0;

    dist = 0.0;

    for (i = 0; i < 3; i++)
    {
      dist   += x[i] * x[i];
      uvw[i] += x[i];
    }

    if (dist <= (tolerance * tolerance)) return 1;
  } // for

  return 0;
}

/*----------------------------------------------------------------------------
 * Calculation of shape functions iteratively
 *
 * This subroutine is executed for all the elements and points which are
 * potentially associated with the coupling.
 * Shape functions are used for localization and for the interpoaltion
 *
 *    elt_type             <-- type of element
 *    element_vertex_num[] <-- id (index) element vertices
 *    vertex_coords[]      <-- array of coordinates
 *    point[]              <-- coordinates point
 *    tolerance            <-- tolerance
 *    shapef               --> shape functions (weigth for the calculation)
 * returns.
 *    shape functions
 *----------------------------------------------------------------------------*/
void
_locate_in_cell_3d(
                   syr_cfd_element_t   elt_type,
                   const ple_lnum_t    element_vertex_num[],
                   const ple_coord_t   vertex_coords[],
                   const ple_coord_t   point[],
                   double              tolerance,
                   double*             shapef
                   )
{
  int j, k, n_vertices;
  ple_lnum_t coord_idx, vertex_id;

  double uvw[3], dist, max_dist;
  double  _vertex_coords[8][3];

  n_vertices = syr_cfd_mesh_n_vertices_element[elt_type];

  // Initialize local element coordinates copy
  for (vertex_id = 0; vertex_id < n_vertices; vertex_id++)
  {
    coord_idx = element_vertex_num[vertex_id] - 1; //<-NOTE: fortran numeration
    for (j = 0; j < 3; j++)  _vertex_coords[vertex_id][j] = vertex_coords[(coord_idx * 3) + j];
  }

  memset(shapef, 0.0, 8*sizeof(double));
  if (elt_type == SYR_CFD_TETRA)
  {
    // Shape functions may be computed directly with tetrahedra
    printf("LOCATE_IN_CELL_3D: ERROR WE HAVE A TETRA!!!");
    exit(0);
  }
  else
  {
    // For cell shapes other than tetrahedra, find shape functions iteratively
    if(_compute_uvw(elt_type, point, _vertex_coords, tolerance, uvw) )
    {
      max_dist = -1.0;

      // For hexahedra, no need to compute shape functions, as the 3 parametric coordinates are simpler to use
      if(elt_type == FVM_CELL_HEXA)
      {
        _compute_shapef_3d(elt_type, uvw, shapef, NULL);

        for (j = 0; j < 3; j++)
        {
          dist = 2.0 * fabs(uvw[j] - 0.5);
          if(max_dist < dist) max_dist = dist;
        }
      }
      else
      {
         // For pyramids ands prisms, we need to compute shape functions
         printf("LOCATE_IN_CELL_3D: PYRAMIDS AND PRISMS NOT IMPLEMENTED");
         exit(0);
      }

    }
    else
    {
      printf("LOCATE_IN_CELL_3D: ERROR calculation natural coordinates HEXA:");
      printf("       <u,v,w>: %f, %f, %f \n", uvw[0], uvw[1], uvw[2]);
      exit(0);
    }

  }
} // _locate_in_cell_3d

//-----------------------------------------------------------------------||---//
static void
__print_cell__(
                   ple_lnum_t         __elt_num,
                   syr_cfd_element_t   elt_type,
                   const ple_lnum_t    element_vertex_num[],
                   const ple_coord_t   vertex_coords[]
                   )
{
  int j, k, n_vertices;
  ple_lnum_t vertex_id;

  n_vertices = syr_cfd_mesh_n_vertices_element[elt_type];

  printf(" %d) ", __elt_num );
  for (vertex_id = 0; vertex_id < n_vertices; vertex_id++) printf(" %d ", element_vertex_num[vertex_id] );
  printf(" \n");

}
/*----------------------------------------------------------------------------
 * Main subroutine for hexahedron
 *
 * This subroutine locates a point inside an hexaedron
 *
 * returns.
 *    __location
 *    __distance
 *----------------------------------------------------------------------------*/
static void
_locate_in_hexa(ple_lnum_t         __elt_num,               //  input: for (i = 0; i < mesh->n_elements; i++)
                const ple_lnum_t   __element_vertex_num[],  //  input: vertex_num + i*stride
                const ple_coord_t  __vertex_coords[],       //  input: vertex_coords
                const ple_coord_t  __point_coords[],        //  input: vertex_coords_j ??
                ple_lnum_t         __n_points_in_extents,   //  input: _query_octree -> n_points_in_extents
                const ple_lnum_t   __points_in_extents[],   //  input: _query_octree -> points_in_extents
                double             __tolerance,             //  input: 1e-8
                ple_lnum_t         __location[],            // output:
                float              __distance[])            // output:
{
  double vol6;
  double dist, max_dist;
  int i, j, k;

/*
  double coord[3];
  for (i = 0; i < __n_points_in_extents; i++)
  {
    j = __points_in_extents[i];
    for(k=0; k<3; k++) coord[k] = __point_coords[i*3 + k];

    printf("%d) ", i);
    for(k=0; k<3; k++) printf("%f ", coord[k] );
    printf("\n ");
  }
// */

//  if(__n_points_in_extents) printf("%d) n_points_in_extents:%d \n ", __elt_num, __n_points_in_extents);

/*
  __print_cell__(
                     __elt_num,
                     FVM_CELL_HEXA,
                     __element_vertex_num,
                     __vertex_coords
                 );
*/

  for (k = 0; k < __n_points_in_extents; k++)
  {
    i = __points_in_extents[k];

    double coord[3];
    for(j=0; j<3; j++) coord[j] = __point_coords[i*3 + j];

    double shapef[8];
    //printf("<element>: %d\n", __elt_num);
    _locate_in_cell_3d(
                        FVM_CELL_HEXA,
                        __element_vertex_num,
                        __vertex_coords,
                        coord,
                        __tolerance,
                        shapef
                      );
    //
    // Criterion to find if the point is inside or outside the element
    //
    max_dist = -1.0;
    for (j = 0; j < 4; j++)
    {
      dist = 2.0 * PLE_ABS(shapef[j] - 0.5);
      if (max_dist < dist) max_dist = dist;
    }
    //
    // Localization (We only localize if the point is inside the element)
    //
    if( (max_dist > -0.5 && max_dist < (1. + 2.*__tolerance)) && (max_dist < __distance[i] || __distance[i] < 0))
    {
      __location[i] = __elt_num;
      __distance[i] = max_dist;
    }

  } // for

}  // _locate_in_hexa

/*
#include "ple_defs.h"
#include "syr_cfd_mesh.h"
#include "syr_cfd_mesh_priv.h"

#include "syr_cfd_point_location.h"
*/

/*============================================================================
 * Local macro definitions
 *============================================================================*/

/* In case <math.h> does not define HUGE_VAL, use a "safe" value */
#if !defined(HUGE_VAL)
#define HUGE_VAL 1.0e+30
#endif

/* Geometric operation macros*/

enum {X, Y, Z};

#define _DOT_PRODUCT(vect1, vect2) \
  (vect1[X] * vect2[X] + vect1[Y] * vect2[Y] + vect1[Z] * vect2[Z])

#define _MODULE(vect) \
  sqrt(vect[X] * vect[X] + vect[Y] * vect[Y] + vect[Z] * vect[Z])

#define _CROSS_PRODUCT(prod_vect, vect1, vect2)  \
  (prod_vect[X] = vect1[Y] * vect2[Z] - vect2[Y] * vect1[Z], \
   prod_vect[Y] = vect2[X] * vect1[Z] - vect1[X] * vect2[Z], \
   prod_vect[Z] = vect1[X] * vect2[Y] - vect2[X] * vect1[Y])

#define _DOT_PRODUCT_2D(vect1, vect2) \
  (vect1[X] * vect2[X] + vect1[Y] * vect2[Y])


#if defined(ENABLE_NLS)

#include <libintl.h>
#define _(String) gettext(String)
#define gettext_noop(String) String
#define N_(String) gettext_noop(String)

#else

#define _(String) String
#define N_(String) String
#define textdomain(Domain)
#define bindtextdomain(Package, Directory)

#endif

/*============================================================================
 * Type definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Structure defining a local octree (3d)
 *----------------------------------------------------------------------------*/

typedef struct {

  ple_lnum_t  octant_id[8];   /* Ids of sub-octants in octree array */
  ple_lnum_t  idx[9];         /* Start index of point list for each octant */
  ple_lnum_t  n_points;       /* Number of points in octree */

} _octant_t;

typedef struct {

  size_t       n_points;      /* Number of points in octree */
  size_t       n_nodes;       /* Current number of nodes in octree */
  size_t       n_nodes_max;   /* Maximum number of nodes in octree */

  double       extents[6];    /* Associated extents */

  ple_lnum_t  *point_ids;     /* Id's of points sorted by octree
                                 (size: n_points + 1) */
  _octant_t   *nodes;         /* Array of octree nodes
                                 (size: n_nodes_max) */

} _octree_t;

/*----------------------------------------------------------------------------
 * Structure defining a local quadtree (2d)
 *----------------------------------------------------------------------------*/

typedef struct {

  ple_lnum_t  quadrant_id[4]; /* Id of sub-quadrants in quadtree array */
  ple_lnum_t  idx[5];         /* Start index of point list for each quadrant */
  ple_lnum_t  n_points;       /* Number of points in quadtree */

} _quadrant_t;

typedef struct {

  size_t        n_points;     /* Number of points in quadtree */
  size_t        n_nodes;      /* Current number of nodes in quadtree */
  size_t        n_nodes_max;  /* Maximum number of nodes in quadtree */

  double        extents[4];   /* Associated extents */

  ple_lnum_t   *point_ids;    /* Id's of points sorted by quadtree
                                 (size: n_points + 1) */
  _quadrant_t  *nodes;        /* Array of quadtree nodes
                                 (size: n_nodes_max) */

} _quadtree_t;

/*============================================================================
 * Static global variables
 *============================================================================*/

static double      _epsilon_denom = 1.e-28; /* Minimum denominator */

static ple_lnum_t  _octree_threshold = 4; /* Number of points in octree node
                                             under which the node is final */

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Adjust extents with sub-extents
 *
 * parameters:
 *   dim         <-- spatial (coordinates) dimension
 *   sub_extents <-> extents associated with element or section:
 *                   x_min, y_min, ..., x_max, y_max, ... (size: 2*dim)
 *   extents     <-> optional section or mesh extents, to be updated:
 *                   x_min, y_min, ..., x_max, y_max, ... (size: 2*dim);
 *                   NULL if unused
 *----------------------------------------------------------------------------*/

inline static void
_update_extents(int               dim,
                double           *sub_extents,
                double           *extents)
{
  int i;

  for (i = 0; i < dim; i++) {
    if (sub_extents[i] < extents[i])
      extents[i] = sub_extents[i];
    if (sub_extents[i+dim] > extents[i+dim])
      extents[i+dim] = sub_extents[i+dim];
  }
}

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
 *   1 if extents intersect, 0 otherwise
 *----------------------------------------------------------------------------*/

inline static int
_intersect_extents(int           dim,
                   const double  extents_1[],
                   const double  extents_2[])
{
  int i;
  int retval = 1;

  for (i = 0; i < dim; i++) {
    if (   (extents_1[i] > extents_2[i + dim])
        || (extents_2[i] > extents_1[i + dim])) {
      retval = 0;
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
 *   1 if point lies within extents, 0 otherwise
 *----------------------------------------------------------------------------*/

inline static int
_within_extents(int                dim,
                const ple_coord_t  coords[],
                const double       extents[])
{
  int i;
  int retval = 1;

  for (i = 0; i < dim; i++) {
    if (   (coords[i] < extents[i])
        || (coords[i] > extents[i + dim])) {
      retval = 0;
      break;
    }
  }

  return retval;
}

/*----------------------------------------------------------------------------
 * Update element extents with a given vertex
 *
 * parameters:
 *   dim               <-- spatial (coordinates) dimension
 *   vertex_id         <-- vertex index (0 to n-1)
 *   vertex_coords     <-- pointer to vertex coordinates
 *   elt_extents       <-> extents associated with element:
 *                         x_min, y_min, ..., x_max, y_max, ... (size: 2*dim)
 *   elt_initialized   <-> are extents already initialized for this vertex
 *                         (for all element vertices except the first) ?
 *----------------------------------------------------------------------------*/

inline static void
_update_elt_extents(int                 dim,
                    ple_lnum_t          vertex_id,
                    const ple_coord_t   vertex_coords[],
                    double              elt_extents[],
                    int                *elt_initialized)
{
  ple_lnum_t  i;
  ple_lnum_t coord_idx = vertex_id;

  if (*elt_initialized == 0) {
    for (i = 0; i < dim; i++) {
      elt_extents[i]       = vertex_coords[(coord_idx * dim) + i];
      elt_extents[i + dim] = vertex_coords[(coord_idx * dim) + i];
    }
    *elt_initialized = 1;
  }
  else {
    for (i = 0; i < dim; i++) {
      if (elt_extents[i]       > vertex_coords[(coord_idx * dim) + i])
        elt_extents[i]       = vertex_coords[(coord_idx * dim) + i];
      if (elt_extents[i + dim] < vertex_coords[(coord_idx * dim) + i])
        elt_extents[i + dim] = vertex_coords[(coord_idx * dim) + i];
    }

  }

}

/*----------------------------------------------------------------------------
 * Adjust element extents with search tolerance and update global extents
 *
 * parameters:
 *   dim         <-- spatial (coordinates) dimension
 *   elt_dim     <-- element dimension
 *   tolerance   <-- addition to local extents of each element:
 *                   extent = base_extent * (1 + tolerance)
 *   elt_extents <-> extents associated with element:
 *                   x_min, y_min, ..., x_max, y_max, ... (size: 2*dim)
 *----------------------------------------------------------------------------*/

inline static void
_elt_extents_finalize(int               dim,
                      int               elt_dim,
                      double            tolerance,
                      double           *elt_extents)
{
  int i;
  double delta[3];

  for (i = 0; i < dim; i++)
    delta[i] = (elt_extents[i+dim] - elt_extents[i]) * tolerance;

  if (elt_dim < dim) {
    double delta_max = delta[0];  /* for 1d or 2d elements, ensure */
    for (i = 0; i < dim; i++) {   /* search extent "thickness" */
      if (delta[i] > delta_max)
        delta_max = delta[i];
    }
    for (i = 0; i < dim; i++)
      delta[i] = delta_max;
  }

  for (i = 0; i < dim; i++) {
    elt_extents[i]     = elt_extents[i]     - delta[i];
    elt_extents[i+dim] = elt_extents[i+dim] + delta[i];
  }
}

/*----------------------------------------------------------------------------
 * Compute extents of a point set
 *
 * parameters:
 *   dim          <-- space dimension of points to locate
 *   n_points     <-- number of points to locate
 *   point_index  <-- optional indirection array to point_coords
 *                    (1 to n_points numbering)
 *   point_coords <-- coordinates of points to locate
 *                    (dimension: dim * n_points)
 *   extents      --> extents associated with mesh:
 *                    x_min, y_min, ..., x_max, y_max, ... (size: 2*dim)
 *----------------------------------------------------------------------------*/

static void
_point_extents(const int            dim,
               const ple_lnum_t     n_points,
               const ple_lnum_t     point_index[],
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

  if (point_index != NULL) {

    for (j = 0; j < n_points; j++) {
      coord_idx = point_index[j] - 1;
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

/*----------------------------------------------------------------------------
 * Build a local octree's leaves.
 *
 * parameters:
 *   extents            <-> extents associated with node:
 *                          x_min, y_min, z_min, x_max, y_max, z_max (size: 6)
 *   point_coords       <-- point coordinates
 *   point_ids_tmp      <-- temporary point indexes
 *   pos_tmp            <-- temporary point position in octree
 *   octree             <-> current octree structure
 *   point_range        <-> start and past-the end index in point_idx
 *                          for current node (size: 2)
 *----------------------------------------------------------------------------*/

static void
_build_octree_leaves(const double        extents[],
                     const ple_coord_t   point_coords[],
                     ple_lnum_t         *point_ids_tmp,
                     _octree_t          *octree,
                     ple_lnum_t          point_range[2])
{
  ple_lnum_t i, j, k, _n_nodes, _n_points, tmp_size;

  ple_lnum_t count[8], idx[9], octant_id[8];
  double mid[3], sub_extents[6];
  _octant_t  *_node;

  int octant_mask[3] = {4, 2, 1}; /* pow(2, 2), pow(2, 1), pow(2,0) */

  _n_nodes = octree->n_nodes;
  tmp_size = octree->n_nodes;

  /* Resize octree if necesary */

  if (octree->n_nodes >= octree->n_nodes_max) {
    if (octree->n_nodes == 0) {
      octree->n_nodes = 1;
      octree->n_nodes_max = 8;
    }
    octree->n_nodes_max *= 2;
    PLE_REALLOC(octree->nodes, octree->n_nodes_max, _octant_t);
  }

  /* Number of points */

  _n_points = point_range[1] - point_range[0];

  /* Extents center */

  for (j = 0; j < 3; j++)
    mid[j]= (extents[j] + extents[j + 3]) * 0.5;

  for (j = 0; j < 8; j++) {
    count[j] = 0;
    octant_id[j] = -1;
  }

  /* Count points in each octant */

  for (i = point_range[0]; i < point_range[1]; i++) {

    for (j = 0, k = 0; j < 3; j++) {
      if (point_coords[octree->point_ids[i]*3 + j] > mid[j])
        k += octant_mask[j];
    }

    count[k] += 1;
  }

  /* Build index */

  idx[0] = 0;
  for (j = 0; j < 8; j++)
    idx[j+1] = idx[j] + count[j];

  for (j = 0; j < 8; j++)
    count[j] = 0;

  for (i = point_range[0], j = 0; i < point_range[1]; i++) {

    for (j = 0, k = 0; j < 3; j++) {
      if (point_coords[octree->point_ids[i]*3 + j] > mid[j])
        k += octant_mask[j];
    }

    point_ids_tmp[idx[k] + count[k]] = octree->point_ids[i];
    count[k] += 1;
  }

  for (i = point_range[0], j = 0; i < point_range[1]; i++, j++)
    octree->point_ids[i] = point_ids_tmp[j];

  for (i = 0; i < 9; i++)
    idx[i] = point_range[0] + idx[i];

  /* Build leaves recursively */

  for (i = 0; i < 8; i++) {

    if (count[i] > _octree_threshold) {

      tmp_size++;

      octant_id[i] = tmp_size;

      if (i < 4) {
        sub_extents[0] = extents[0];
        sub_extents[3] = mid[0];
      }
      else {
        sub_extents[0] = mid[0];
        sub_extents[3] = extents[3];
      }
      /* 1.0e-12 term in assert() used to allow for
         truncation error in for xmin = xmax case */
      assert(sub_extents[0] < sub_extents[3] + 1.0e-12);

      if (i%4 < 2) {
        sub_extents[1] = extents[1];
        sub_extents[4] = mid[1];
      }
      else {
        sub_extents[1] = mid[1];
        sub_extents[4] = extents[4];
      }
      assert(sub_extents[1] < sub_extents[4] + 1.0e-12);

      if (i%2 < 1) {
        sub_extents[2] = extents[2];
        sub_extents[5] = mid[2];
      }
      else {
        sub_extents[2] = mid[2];
        sub_extents[5] = extents[5];
      }
      assert(sub_extents[2] < sub_extents[5] + 1.0e-12);

      octree->n_nodes = tmp_size;

      _build_octree_leaves(sub_extents,
                           point_coords,
                           point_ids_tmp,
                           octree,
                           idx + i);

      tmp_size = octree->n_nodes;
    }

  }

  /* Finalize node */

  _node = octree->nodes + _n_nodes;

  for (i = 0; i < 9; i++)
    _node->idx[i] = idx[i];

  for (i = 0; i < 8; i++)
    _node->octant_id[i] = octant_id[i];

  _node->n_points = _n_points;
}

/*----------------------------------------------------------------------------
 * Build an octree structure to locate 3d points in mesh.
 *
 * parameters:
 *   n_points        <-- number of points to locate
 *   point_coords    <-- point coordinates
 *
 * returns:
 *   pointer to local octree structure
 *----------------------------------------------------------------------------*/

static _octree_t
_build_octree(ple_lnum_t         n_points,
              const ple_coord_t  point_coords[])
{
  size_t i;
  ple_lnum_t point_range[2];
  _octree_t _octree;

  int *point_ids_tmp = NULL;

  /* Initialization */

  point_range[0] = 0;
  point_range[1] = n_points;

  _octree.n_points = n_points;
  _octree.n_nodes = 0;
  _octree.n_nodes_max = 0;
  _octree.nodes = NULL;
  _octree.point_ids = NULL;

  if (n_points > 0) {

    _point_extents(3,
                   n_points,
                   NULL,
                   point_coords,
                   _octree.extents);

    PLE_MALLOC(_octree.point_ids, _octree.n_points, ple_lnum_t);

    for (i = 0; i < _octree.n_points; i++)
      _octree.point_ids[i] = i;

    PLE_MALLOC(point_ids_tmp, n_points, int);

    _build_octree_leaves(_octree.extents,
                         point_coords,
                         point_ids_tmp,
                         &_octree,
                         point_range);

    PLE_FREE(point_ids_tmp);

  }

  return _octree;
}

/*----------------------------------------------------------------------------
 * Free an octree structure.
 *
 * parameters:
 *   octree <-> octree structure whose elements are to be freed
 *
 * returns:
 *   pointer to local octree structure
 *----------------------------------------------------------------------------*/

static void
_free_octree(_octree_t *octree)
{

  octree->n_points = 0;
  octree->n_nodes = 0;
  octree->n_nodes_max = 0;

  PLE_FREE(octree->nodes);
  PLE_FREE(octree->point_ids);
}

/*----------------------------------------------------------------------------
 * locate points in box defined by extents in an octant.
 *
 * parameters:
 *   search_extents  <-- extents associated with element:
 *   point_coords    <-- point coordinates
 *                       x_min, y_min, z_min, x_max, y_max, z_max (size: 6)
 *   octree          <-- point octree
 *   node_extents    <-- extents of octant
 *   node_id         <-- id of node in octree
 *   loc_point_ids   --> ids of located points (max size: octree->n_points)
 *   n_loc_points    --> number of points located
 *----------------------------------------------------------------------------*/

static void
_query_octree_node(const double        extents[],
                   const ple_coord_t   point_coords[],
                   const _octree_t    *octree,
                   const double        node_extents[],
                   int                 node_id,
                   ple_lnum_t         *loc_point_ids,
                   ple_lnum_t         *n_loc_points)
{
  _octant_t *node;
  int i, j, k;
  int dim = 3;
  double sub_extents[6], mid[3];

  node = octree->nodes + node_id;

  if (_intersect_extents(dim, node_extents, extents)) {

    for (j = 0; j < dim; j++)
      mid[j]= (node_extents[j] + node_extents[j + dim]) * 0.5;

    /* Loop on node leaves */

    for (i = 0; i < 8; i++) {

      /* Compute octant extents */

      if (i < 4) {
        sub_extents[0] = node_extents[0];
        sub_extents[3] = mid[0];
      }
      else {
        sub_extents[0] = mid[0];
        sub_extents[3] = node_extents[3];
      }

      if (i%4 < 2) {
        sub_extents[1] = node_extents[1];
        sub_extents[4] = mid[1];
      }
      else {
        sub_extents[1] = mid[1];
        sub_extents[4] = node_extents[4];
      }
      if (i%2 < 1) {
        sub_extents[2] = node_extents[2];
        sub_extents[5] = mid[2];
      }
      else {
        sub_extents[2] = mid[2];
        sub_extents[5] = node_extents[5];
      }

      /* Search recursively if octant is not final */

      if (node->octant_id[i] > -1)
        _query_octree_node(extents,
                           point_coords,
                           octree,
                           sub_extents,
                           node->octant_id[i],
                           loc_point_ids,
                           n_loc_points);

      else {

        /* Update list of points located */

        if (_intersect_extents(dim, sub_extents, extents)) {

          for (k = node->idx[i]; k < node->idx[i+1]; k++) {

            ple_lnum_t point_id = octree->point_ids[k];

            if (_within_extents(dim, point_coords + point_id*dim, extents)) {
              loc_point_ids[*n_loc_points] = point_id;
              (*n_loc_points)++;
            }

          }
        }

      }

    } /* End of loop on node leaves */

  }

}

/*----------------------------------------------------------------------------
 * locate points in box defined by extents using an octree.
 *
 * parameters:
 *   search_extents  <-- extents associated with element:
 *                       x_min, y_min, ..., x_max, y_max, ... (size: 2*dim)
 *   point_coords    <-- point coordinates
 *   octree          <-- point octree
 *   n_loc_points    --> number of points located
 *   loc_point_ids   --> ids of located points (max size: octree->n_points)
 *----------------------------------------------------------------------------*/

static void
_query_octree(const double        extents[],
              const ple_coord_t   point_coords[],
              const _octree_t    *octree,
              ple_lnum_t         *n_loc_points,
              ple_lnum_t          loc_point_ids[])
{
  *n_loc_points = 0;

  if (octree->n_points > 0)
    _query_octree_node(extents,
                       point_coords,
                       octree,
                       octree->extents,
                       0,
                       loc_point_ids,
                       n_loc_points);
}

/*----------------------------------------------------------------------------
 * Build a local quadtree's leaves.
 *
 * parameters:
 *   extents            <-> extents associated with node:
 *                          x_min, y_min, x_max, y_max, (size: 4)
 *   point_coords       <-- point coordinates
 *   point_ids_tmp      <-- temporary point indexes
 *   pos_tmp            <-- temporary point position in quadtree
 *   octree             <-> current quadtree structure
 *   point_range        <-> start and past-the end index in point_idx
 *                          for current node (size: 2)
 *----------------------------------------------------------------------------*/

static void
_build_quadtree_leaves(const double        extents[],
                       const ple_coord_t   point_coords[],
                       ple_lnum_t         *point_ids_tmp,
                       _quadtree_t        *quadtree,
                       ple_lnum_t          point_range[2])
{
  ple_lnum_t i, j, k, _n_nodes, _n_points, tmp_size;

  ple_lnum_t count[4], idx[5], quadrant_id[4];
  double mid[2], sub_extents[4];
  _quadrant_t  *_node;

  int quadrant_mask[2] = {2, 1}; /* pow(2, 1), pow(2,0) */

  _n_nodes = quadtree->n_nodes;
  tmp_size = quadtree->n_nodes;

  /* Resize quadtree if necesary */

  if (quadtree->n_nodes >= quadtree->n_nodes_max) {
    if (quadtree->n_nodes == 0) {
      quadtree->n_nodes = 1;
      quadtree->n_nodes_max = 4;
    }
    quadtree->n_nodes_max *= 2;
    PLE_REALLOC(quadtree->nodes, quadtree->n_nodes_max, _quadrant_t);
  }

  /* Number of points */

  _n_points = point_range[1] - point_range[0];

  /* Extents center */

  for (j = 0; j < 2; j++)
    mid[j]= (extents[j] + extents[j + 2]) * 0.5;

  for (j = 0; j < 4; j++) {
    count [j] = 0;
    quadrant_id[j] = -1;
  }

  /* Count points in each quadrant */

  for (i = point_range[0]; i < point_range[1]; i++) {

    for (j = 0, k = 0; j < 2; j++) {
      if (point_coords[quadtree->point_ids[i]*2 + j] > mid[j])
        k += quadrant_mask[j];
    }

    count[k] += 1;
  }

  /* Build index */

  idx[0] = 0;
  for (j = 0; j < 4; j++)
    idx[j+1] = idx[j] + count[j];

  for (j = 0; j < 4; j++)
    count[j] = 0;

  for (i = point_range[0], j = 0; i < point_range[1]; i++) {

    for (j = 0, k = 0; j < 2; j++) {
      if (point_coords[quadtree->point_ids[i]*2 + j] > mid[j])
        k += quadrant_mask[j];
    }

    point_ids_tmp[idx[k] + count[k]] = quadtree->point_ids[i];
    count[k] += 1;
  }

  for (i = point_range[0], j = 0; i < point_range[1]; i++, j++)
    quadtree->point_ids[i] = point_ids_tmp[j];

  for (i = 0; i < 5; i++)
    idx[i] = point_range[0] + idx[i];

  /* Build leaves recursively */

  for (i = 0; i < 4; i++) {

    if (count[i] > _octree_threshold) {

      tmp_size++;

      quadrant_id[i] = tmp_size;

      if (i < 2) {
        sub_extents[0] = extents[0];
        sub_extents[2] = mid[0];
      }
      else {
        sub_extents[0] = mid[0];
        sub_extents[2] = extents[2];
      }
      assert(sub_extents[0] < sub_extents[2] + 1.0e-12);

      if (i%2 < 1) {
        sub_extents[1] = extents[1];
        sub_extents[3] = mid[1];
      }
      else {
        sub_extents[1] = mid[1];
        sub_extents[3] = extents[3];
      }
      assert(sub_extents[1] < sub_extents[3] + 1.0e-12);

      quadtree->n_nodes = tmp_size;

      _build_quadtree_leaves(sub_extents,
                             point_coords,
                             point_ids_tmp,
                             quadtree,
                             idx + i);

      tmp_size = quadtree->n_nodes;
    }

  }

  /* Finalize node */

  _node = quadtree->nodes + _n_nodes;

  for (i = 0; i < 5; i++)
    _node->idx[i] = idx[i];

  for (i = 0; i < 4; i++)
    _node->quadrant_id[i] = quadrant_id[i];

  _node->n_points = _n_points;
}

/*----------------------------------------------------------------------------
 * Build a quadtree structure to locate 2d points in mesh.
 *
 * parameters:
 *   n_points        <-- number of points to locate
 *   point_coords    <-- point coordinates
 *
 * returns:
 *   pointer to local quadtree structure
 *----------------------------------------------------------------------------*/

static _quadtree_t
_build_quadtree(ple_lnum_t         n_points,
                const ple_coord_t  point_coords[])
{
  size_t i;
  ple_lnum_t point_range[2];
  _quadtree_t _quadtree;

  int *point_ids_tmp = NULL;

  /* Initialization */

  point_range[0] = 0;
  point_range[1] = n_points;

  _quadtree.n_points = n_points;
  _quadtree.n_nodes = 0;
  _quadtree.n_nodes_max = 0;
  _quadtree.nodes = NULL;
  _quadtree.point_ids = NULL;

  if (n_points > 0) {

    _point_extents(2,
                   n_points,
                   NULL,
                   point_coords,
                   _quadtree.extents);

    PLE_MALLOC(_quadtree.point_ids, _quadtree.n_points, ple_lnum_t);

    for (i = 0; i < _quadtree.n_points; i++)
      _quadtree.point_ids[i] = i;

    PLE_MALLOC(point_ids_tmp, n_points, int);

    _build_quadtree_leaves(_quadtree.extents,
                           point_coords,
                           point_ids_tmp,
                           &_quadtree,
                           point_range);

    PLE_FREE(point_ids_tmp);

  }

  return _quadtree;
}

/*----------------------------------------------------------------------------
 * Free a quadtree structure.
 *
 * parameters:
 *   quadtree <-> quadtree structure whose elements are to be freed
 *
 * returns:
 *   pointer to local quadtree structure
 *----------------------------------------------------------------------------*/

static void
_free_quadtree(_quadtree_t *quadtree)
{

  quadtree->n_points = 0;
  quadtree->n_nodes = 0;
  quadtree->n_nodes_max = 0;

  PLE_FREE(quadtree->nodes);
  PLE_FREE(quadtree->point_ids);
}

/*----------------------------------------------------------------------------
 * locate points in box defined by extents in a quadrant.
 *
 * parameters:
 *   search_extents  <-- extents associated with element:
 *   point_coords    <-- point coordinates
 *                       x_min, y_min, x_max, y_max (size: 4)
 *   quadtree        <-- point quadtree
 *   node_extents    <-- extents of quadrant
 *   node_id         <-- if of node in octree
 *   loc_point_ids   --> ids of located points (max size: quadtree->n_points)
 *   n_loc_points    --> number of points located
 *----------------------------------------------------------------------------*/

static void
_query_quadtree_node(const double        extents[],
                     const ple_coord_t   point_coords[],
                     const _quadtree_t  *quadtree,
                     const double        node_extents[],
                     int                 node_id,
                     ple_lnum_t         *loc_point_ids,
                     ple_lnum_t         *n_loc_points)
{
  _quadrant_t *node;
  int i, j, k;
  double sub_extents[4], mid[2];

  node = quadtree->nodes + node_id;

  if (_intersect_extents(2, node_extents, extents)) {

    for (j = 0; j < 2; j++)
      mid[j]= (node_extents[j] + node_extents[j + 2]) * 0.5;

    /* Loop on node leaves */

    for (i = 0; i < 4; i++) {

      /* Compute quadrant extents */

      if (i < 2) {
        sub_extents[0] = node_extents[0];
        sub_extents[2] = mid[0];
      }
      else {
        sub_extents[0] = mid[0];
        sub_extents[2] = node_extents[2];
      }

      if (i%2 < 1) {
        sub_extents[1] = node_extents[1];
        sub_extents[3] = mid[1];
      }
      else {
        sub_extents[1] = mid[1];
        sub_extents[3] = node_extents[3];
      }

      /* Search recursively if quadrant is not final */

      if (node->quadrant_id[i] > -1)
        _query_quadtree_node(extents,
                             point_coords,
                             quadtree,
                             sub_extents,
                             node->quadrant_id[i],
                             loc_point_ids,
                             n_loc_points);

      else {

        /* Update list of points located */

        if (_intersect_extents(2, sub_extents, extents)) {

          for (k = node->idx[i]; k < node->idx[i+1]; k++) {

            ple_lnum_t point_id = quadtree->point_ids[k];

            if (_within_extents(2, point_coords + point_id*2, extents)) {
              loc_point_ids[*n_loc_points] = point_id;
              (*n_loc_points)++;
            }

          }
        }

      }

    } /* End of loop on node leaves */

  }

}

/*----------------------------------------------------------------------------
 * locate points in box defined by extents using an octree.
 *
 * parameters:
 *   search_extents  <-- extents associated with element:
 *   point_coords    <-- point coordinates
 *                       x_min, y_min, x_max, y_max (size: 4)
 *   quadtree        <-- point quadtree
 *   n_loc_points    --> number of points located
 *   loc_point_ids   --> ids of located points (max size: octree->n_points)
 *----------------------------------------------------------------------------*/

static void
_query_quadtree(const double        extents[],
                const ple_coord_t   point_coords[],
                const _quadtree_t  *quadtree,
                ple_lnum_t         *n_loc_points,
                ple_lnum_t          loc_point_ids[])
{
  *n_loc_points = 0;

  if (quadtree->n_points > 0)
    _query_quadtree_node(extents,
                         point_coords,
                         quadtree,
                         quadtree->extents,
                         0,
                         loc_point_ids,
                         n_loc_points);
}

/*----------------------------------------------------------------------------
 * Locate points on a 2d edge, and update the location[] and distance[]
 * arrays associated with the point set.
 *
 * parameters:
 *   elt_num             <-- number of element corresponding to extents
 *   element_vertex_num  <-- element vertex numbers
 *   vertex_coords       <-- pointer to vertex coordinates
 *   point_coords        <-- point coordinates
 *   n_points_in_extent  <-- number of points in extents
 *   points_in_extent    <-- ids of points in extents
 *   n_points_in_extent  <-- number of points in extents
 *   points_in_extent    <-- ids of points in extents
 *   tolerance           <-- associated tolerance (considered infinite if < 0)
 *   location            <-> number of element containing or closest to each
 *                           point (size: n_points)
 *   distance            <-> distance from point to element indicated by
 *                           location[]: < 0 if unlocated, absolute distance
 *                           to element if located (size: n_points)
 *----------------------------------------------------------------------------*/

static void
_locate_on_edge_2d(ple_lnum_t           elt_num,
                   const ple_lnum_t     element_vertex_num[],
                   const ple_coord_t    vertex_coords[],
                   const ple_coord_t    point_coords[],
                   ple_lnum_t           n_points_in_extents,
                   const ple_lnum_t     points_in_extents[],
                   double               tolerance,
                   ple_lnum_t           location[],
                   float                distance[])
{
  ple_lnum_t  i, j, k, coord_idx_0, coord_idx_1;

  double u[2], v[2];
  double uv, len2, isop_0;
  double dist2, epsilon2, vertex_dist2;

  /* vertex index of the edge studied */

  coord_idx_0 = element_vertex_num[0] - 1;
  coord_idx_1 = element_vertex_num[1] - 1;

  /* Calculate edge vector and length */

  for (j = 0; j < 2; j++)
    u[j] =   vertex_coords[(coord_idx_1*2) + j]
           - vertex_coords[(coord_idx_0*2) + j];

  len2 = _DOT_PRODUCT_2D(u, u);

  if (tolerance < 0.0)
    epsilon2 = HUGE_VAL;

  else
    epsilon2 = len2*tolerance*tolerance;

  /* Loop on points resulting from extent query */

  for (k = 0; k < n_points_in_extents; k++) {

    i =  points_in_extents[k];

    vertex_dist2 = distance[i]*distance[i];

    /* Calculate linear coordinates of projection of point on edge axis */

    for (j = 0; j < 2; j++)
      v[j] = point_coords[i*2 + j] - vertex_coords[(coord_idx_0*2) + j];

    uv = u[0]*v[0] + u[1]*v[1];

    if (len2 >= _epsilon_denom)
      isop_0 = uv / len2;
    else
      isop_0 = 0.5; /* for degenerate edges, use center */

    /* Set v to be the vector from the point to the closest point on
       the segment (if isop_0 < 0, v is already that vector) */

    if (isop_0 >= 1) {
      for (j = 0; j < 2; j++)
        v[j] = point_coords[i*2 + j] - vertex_coords[(coord_idx_1*2) + j];
    }
    else if (isop_0 > 0) {
      for (j = 0; j < 2; j++)
        v[j] -= isop_0*u[j];
    }

    /* Distance between point to locate and its projection */

    dist2 = _DOT_PRODUCT_2D(v, v);

    if (dist2 < epsilon2 && (dist2 < vertex_dist2 || distance[i] < 0.0)) {
      location[i] = elt_num;
      distance[i] = sqrt(dist2);
    }

  } /* End of loop on points resulting from extent query */

}

/*----------------------------------------------------------------------------
 * Locate points in a 3d triangle, and updates the location[]
 * and distance[] arrays associated with the point set.
 *
 * Barycentric coordinates are used to locate the projection of points.
 *
 * parameters:
 *   elt_num             <-- number of element corresponding to triangle
 *   element_vertex_num  <-- element vertex numbers
 *   vertex_coords       <-- pointer to vertex coordinates
 *   point_coords        <-- point coordinates
 *   n_points_in_extents <-- number of points in extents
 *   points_in_extents   <-- ids of points in extents
 *   tolerance           <-- associated tolerance (considered infinite if < 0)
 *   location            <-> number of element containing or closest to each
 *                           point (size: n_points)
 *   distance            <-> distance from point to element indicated by
 *                           location[]: < 0 if unlocated, absolute distance
 *                           to element if located (size: n_points)
 *----------------------------------------------------------------------------*/

static void
_locate_on_triangle_3d(ple_lnum_t           elt_num,
                       const ple_lnum_t     element_vertex_num[],
                       const ple_coord_t    vertex_coords[],
                       const ple_coord_t    point_coords[],
                       ple_lnum_t           n_points_in_extents,
                       const ple_lnum_t     points_in_extents[],
                       const double         tolerance,
                       ple_lnum_t           location[],
                       float                distance[])
{
  ple_lnum_t  i, j, k, coord_idx_0, coord_idx_1, coord_idx_2;

  double t[3], u[3], v[3], w[3], vect_tmp[3];
  double uu, vv, uv, ut, vt, ww, det, tmp_max;
  double epsilon2, dist2, vertex_dist2, isop_0, isop_1;

  double tolerance2 = tolerance*tolerance;

  /* vertex index of the triangle studied */

  coord_idx_0 = element_vertex_num[0] - 1;
  coord_idx_1 = element_vertex_num[1] - 1;
  coord_idx_2 = element_vertex_num[2] - 1;

    /* Calculate triangle-constant values for barycentric coordinates */

  for (j = 0; j < 3; j++) {
    u[j] = - vertex_coords[(coord_idx_0*3) + j]
           + vertex_coords[(coord_idx_1*3) + j];
    v[j] = - vertex_coords[(coord_idx_0*3) + j]
           + vertex_coords[(coord_idx_2*3) + j];
    w[j] =   vertex_coords[(coord_idx_1*3) + j]
           - vertex_coords[(coord_idx_2*3) + j];
  }

  uu = _DOT_PRODUCT(u, u);
  vv = _DOT_PRODUCT(v, v);
  ww = _DOT_PRODUCT(w, w);
  uv = _DOT_PRODUCT(u, v);

  det = (uu*vv - uv*uv);

  /* epsilon2 is based on maximum edge length (squared) */

  tmp_max = PLE_MAX(vv, ww);

  if (tolerance < 0.)
    epsilon2 = HUGE_VAL;
  else
    epsilon2 = PLE_MAX(uu, tmp_max) * tolerance2;

  /* Loop on points resulting from extent query */

  for (k = 0; k < n_points_in_extents; k++) {

    i =  points_in_extents[k];

    vertex_dist2 = distance[i]*distance[i];

    /* Calculation of the barycenter coordinates for the projected node */

    for (j = 0; j < 3; j++)
      t[j] = - vertex_coords[(coord_idx_0*3) + j]
             + point_coords[i*3 + j];

    ut = _DOT_PRODUCT(u, t);
    vt = _DOT_PRODUCT(v, t);

    if (det >= _epsilon_denom) {
      isop_0 = (ut*vv - vt*uv) / det;
      isop_1 = (uu*vt - uv*ut) / det;
    }
    else { /* for degenerate triangles, use triangle center */
      isop_0 = 0.5;
      isop_1 = 0.5;
    }

    _CROSS_PRODUCT(vect_tmp, u, v);

    /* if the projected point is not on triangle, we project it
       on the nearest edge or node */

    if (isop_0 < 0.)
      isop_0 = 0.;

    if (isop_1 < 0.)
      isop_1 = 0.;

    if ((1.0 - isop_0 - isop_1) < 0.) {
      isop_0 = isop_0 / (isop_0 + isop_1);
      isop_1 = isop_1 / (isop_0 + isop_1);
    }

    /* re-use vect_tmp */

    for (j = 0; j < 3; j++)
      vect_tmp[j] =   vertex_coords[coord_idx_0*3 + j]
                    + u[j]*isop_0
                    + v[j]*isop_1
                    - point_coords[i*3 + j];

    /* distance between point to locate and its projection */

    dist2 = _DOT_PRODUCT(vect_tmp, vect_tmp);

    if (dist2 < epsilon2 && (dist2 < vertex_dist2 || distance[i] < 0.0)) {
      location[i] = elt_num;
      distance[i] = sqrt(dist2);
    }

  } /* End of loop on points resulting from extent query */
}

/*----------------------------------------------------------------------------
 * Locate points in a 2d triangle, and updates the location[]
 * and distance[] arrays associated with the point set.
 *
 * Barycentric coordinates are used to locate the projection of points.
 *
 * parameters:
 *   elt_num             <-- number of element corresponding to extents
 *   element_vertex_num  <-- element vertex numbers
 *   vertex_coords       <-- pointer to vertex coordinates
 *   point_coords        <-- point coordinates
 *   n_points_in_extent  <-- number of points in extents
 *   points_in_extent    <-- ids of points in extents
 *   n_points_in_extent  <-- number of points in extents
 *   points_in_extent    <-- ids of points in extents
 *   tolerance           <-- associated tolerance
 *   location            <-> number of element containing or closest to each
 *                           point (size: n_points)
 *   distance            <-> distance from point to element indicated by
 *                           location[]: < 0 if unlocated, 0 - 1 if inside,
 *                           > 1 if outside (size: n_points)
 *----------------------------------------------------------------------------*/

static void
_locate_in_triangle_2d(ple_lnum_t           elt_num,
                       const ple_lnum_t     element_vertex_num[],
                       const ple_coord_t    vertex_coords[],
                       const ple_coord_t    point_coords[],
                       ple_lnum_t           n_points_in_extents,
                       const ple_lnum_t     points_in_extents[],
                       double               tolerance,
                       ple_lnum_t           location[],
                       float                distance[])
{
  ple_lnum_t  i, j, k, coord_idx_0, coord_idx_1, coord_idx_2;

  double t[2], u[2], v[2], shapef[3];
  double uu, vv, uv, ut, vt, det;
  double dist, max_dist, isop_0, isop_1;

  /* vertex index of the triangle studied */

  coord_idx_0 = element_vertex_num[0] - 1;
  coord_idx_1 = element_vertex_num[1] - 1;
  coord_idx_2 = element_vertex_num[2] - 1;

  /* Calculate triangle-constant values for barycentric coordinates */

  for (j = 0; j < 2; j++) {
    u[j] = - vertex_coords[(coord_idx_0*2) + j]
           + vertex_coords[(coord_idx_1*2) + j];
    v[j] = - vertex_coords[(coord_idx_0*2) + j]
           + vertex_coords[(coord_idx_2*2) + j];
  }

  uu = _DOT_PRODUCT_2D(u, u);
  vv = _DOT_PRODUCT_2D(v, v);
  uv = _DOT_PRODUCT_2D(u, v);

  det = (uu*vv - uv*uv);

  if (det < _epsilon_denom)
    return;

  /* Loop on points resulting from extent query */

  for (k = 0; k < n_points_in_extents; k++) {

    i =  points_in_extents[k];

    /* Calculation of the barycenter coordinates for the projected node */

    for (j = 0; j < 2; j++)
      t[j] = - vertex_coords[(coord_idx_0*2) + j]
             + point_coords[i*2 + j];

    ut = u[0]*t[0] + u[1]*t[1];
    vt = v[0]*t[0] + v[1]*t[1];

    isop_0 = (ut*vv - vt*uv) / det;
    isop_1 = (uu*vt - uv*ut) / det;

    shapef[0] = 1. - isop_0 - isop_1;
    shapef[1] =      isop_0;
    shapef[2] =               isop_1;

    max_dist = -1.0;

    for (j = 0; j < 3; j++){

      dist = 2.*PLE_ABS(shapef[j] - 0.5);

      if (max_dist < dist)
        max_dist = dist;
    }

    if (   (max_dist > -0.5 && max_dist < (1. + 2.*tolerance))
        && (max_dist < distance[i] || distance[i] < 0)) {
      location[i] = elt_num;
      distance[i] = max_dist;
    }

  } /* End of loop on points resulting from extent query */
}

/*----------------------------------------------------------------------------
 * Locate points in a tetrahedron whose coordinates are pre-computed,
 * updating the location[] and distance[] arrays associated with a set
 * of points.
 *
 * parameters:
 *   elt_num             <-- element number
 *   element_vertex_num  <-- element vertex numbers
 *   vertex_coords       <-- pointer to vertex coordinates
 *   point_coords        <-- point coordinates
 *   n_points_in_extents <-- number of points in element extents
 *   points_in_extents   <-- ids of points in extents
 *   tolerance           <-- associated tolerance
 *   location            <-> number of element containing or closest to each
 *                           point (size: n_points)
 *   distance            <-> distance from point to element indicated by
 *                           location[]: < 0 if unlocated, 0 - 1 if inside,
 *                           > 1 if outside (size: n_points)
 *----------------------------------------------------------------------------*/

static void
_locate_in_tetra(ple_lnum_t         elt_num,
                 const ple_lnum_t   element_vertex_num[],
                 const ple_coord_t  vertex_coords[],
                 const ple_coord_t  point_coords[],
                 ple_lnum_t         n_points_in_extents,
                 const ple_lnum_t   points_in_extents[],
                 double             tolerance,
                 ple_lnum_t         location[],
                 float              distance[])
{
  double vol6;
  double dist, max_dist;
  int i, j, k;

  double isop_0, isop_1, isop_2;
  double t00, t10, t20, t01, t02, t03, t11, t12, t13, t21, t22, t23;
  double v01[3], v02[3], v03[3], shapef[4];

  ple_coord_t  tetra_coords[4][3];

  for (i = 0; i < 4; i++) {
    ple_lnum_t coord_idx = element_vertex_num[i] -1;
    for (j = 0; j < 3; j++)
      tetra_coords[i][j] = vertex_coords[(coord_idx * 3) + j];
  }

  for (i = 0; i < 3; i++) {
    v01[i] = tetra_coords[1][i] - tetra_coords[0][i];
    v02[i] = tetra_coords[2][i] - tetra_coords[0][i];
    v03[i] = tetra_coords[3][i] - tetra_coords[0][i];
  }

  vol6 = fabs(  v01[0] * (v02[1]*v03[2] - v02[2]*v03[1])
              - v02[0] * (v01[1]*v03[2] - v01[2]*v03[1])
              + v03[0] * (v01[1]*v02[2] - v01[2]*v02[1]));

  if (vol6 < _epsilon_denom)
    return;

  for (k = 0; k < n_points_in_extents; k++) {

    i = points_in_extents[k];

    t00  =   point_coords[i*3]     - tetra_coords[0][0];
    t10  =   point_coords[i*3 + 1] - tetra_coords[0][1];
    t20  =   point_coords[i*3 + 2] - tetra_coords[0][2];

    t01  = - tetra_coords[0][0] + tetra_coords[1][0];
    t02  = - tetra_coords[0][0] + tetra_coords[2][0];
    t03  = - tetra_coords[0][0] + tetra_coords[3][0];

    t11  = - tetra_coords[0][1] + tetra_coords[1][1];
    t12  = - tetra_coords[0][1] + tetra_coords[2][1];
    t13  = - tetra_coords[0][1] + tetra_coords[3][1];

    t21  = - tetra_coords[0][2] + tetra_coords[1][2];
    t22  = - tetra_coords[0][2] + tetra_coords[2][2];
    t23  = - tetra_coords[0][2] + tetra_coords[3][2];

    isop_0 = (  t00 * (t12*t23 - t13*t22)
              - t10 * (t02*t23 - t22*t03)
              + t20 * (t02*t13 - t12*t03)) / vol6;
    isop_1 = (- t00 * (t11*t23 - t13*t21)
              + t10 * (t01*t23 - t21*t03)
              - t20 * (t01*t13 - t03*t11)) / vol6;
    isop_2 = (  t00 * (t11*t22 - t21*t12)
              - t10 * (t01*t22 - t21*t02)
              + t20 * (t01*t12 - t11*t02)) / vol6;

    shapef[0] = 1. - isop_0 - isop_1 - isop_2;
    shapef[1] =      isop_0;
    shapef[2] =               isop_1;
    shapef[3] =                        isop_2;

    max_dist = -1.0;

    for (j = 0; j < 4; j++){

      dist = 2.*PLE_ABS(shapef[j] - 0.5);

      if (max_dist < dist)
        max_dist = dist;
    }

    if (   (max_dist > -0.5 && max_dist < (1. + 2.*tolerance))
        && (max_dist < distance[i] || distance[i] < 0)) {
      location[i] = elt_num;
      distance[i] = max_dist;
    }

  }
}

/*----------------------------------------------------------------------------
 * Find elements in a given section containing 3d points: updates the
 * location[] and distance[] arrays associated with a set of points
 * for points that are in an element of this section, or closer to one
 * than to previously encountered elements.
 *
 * parameters:
 *   mesh              <-- pointer to mesh representation structure
 *   tolerance         <-- associated tolerance
 *   point_coords      <-- point coordinates
 *   octree            <-- point octree
 *   points_in_extents <-- array for query of ids of points in extents
 *                         (size: octree->n_points, less usually needed)
 *   location          <-> number of element containing or closest to each
 *                         point (size: n_points)
 *   distance          <-> distance from point to element indicated by
 *                         location[]: < 0 if unlocated; 0 - 1 if inside,
 *                         and > 1 if outside a volume element, or absolute
 *                         distance to a surface element (size: n_points)
 *----------------------------------------------------------------------------*/

static void
_mesh_locate_3d(const syr_cfd_mesh_t  *mesh,
                double                 tolerance,
                const ple_coord_t      point_coords[],
                _octree_t             *octree,
                ple_lnum_t             points_in_extents[],
                ple_lnum_t             location[],
                float                  distance[])
{
  ple_lnum_t  i, j;
  double elt_extents[6];

  ple_lnum_t n_points_in_extents = 0;

  const int stride = syr_cfd_mesh_n_vertices_element[mesh->element_type];
  const int elt_dim = stride - 1; /* true for simplexes */

  if (mesh->n_elements == 0)
    return;

  for (i = 0; i < mesh->n_elements; i++) {

    int elt_initialized = 0;

    for (j = 0; j < stride; j++) {
      ple_lnum_t vertex_id = mesh->vertex_num[i*stride + j] - 1;
      _update_elt_extents(3,
                          vertex_id,
                          mesh->vertex_coords,
                          elt_extents,
                          &elt_initialized);
    }
    _elt_extents_finalize(3,
                          elt_dim,
                          tolerance,
                          elt_extents);

    _query_octree(elt_extents,
                  point_coords,
                  octree,
                  &n_points_in_extents,
                  points_in_extents);

    if (mesh->element_type == SYR_CFD_TETRA)
      _locate_in_tetra(i + 1,
                       mesh->vertex_num + i*stride,
                       mesh->vertex_coords,
                       point_coords,
                       n_points_in_extents,
                       points_in_extents,
                       tolerance,
                       location,
                       distance);

    else if (mesh->element_type == SYR_CFD_TRIA)
      _locate_on_triangle_3d(i + 1,
                             mesh->vertex_num + i*stride,
                             mesh->vertex_coords,
                             point_coords,
                             n_points_in_extents,
                             points_in_extents,
                             tolerance,
                             location,
                             distance);

    else if (mesh->element_type == FVM_CELL_HEXA)    // 2019FEB21    
    {
      _locate_in_hexa(i + 1,
                      mesh->vertex_num + i*stride,
                      mesh->vertex_coords,
                      point_coords,
                      n_points_in_extents,
                      points_in_extents,
                      tolerance,
                      location,
                      distance);
    }
    else                                            // 2019FEB21    
      assert(0); // mesh->element_type ??  

  }
}

/*----------------------------------------------------------------------------
 * Find elements in a given section closest to 3d points: updates the
 * location[] and distance[] arrays associated with a set of points
 * for points that are in an element of this section, or closer to one
 * than to previously encountered elements.
 *
 * parameters:
 *   mesh              <-- pointer to mesh representation structure
 *   point_coords      <-- point coordinates
 *   n_point_ids       <-- number of points to locate
 *   point_id          <-- ids of points to locate
 *   location          <-> number of element containing or closest to each
 *                         point (size: n_points)
 *   distance          <-> distance from point to element indicated by
 *                         location[]: < 0 if unlocated, or absolute
 *                         distance to a surface element (size: n_points)
 *----------------------------------------------------------------------------*/

static void
_mesh_closest_3d(const syr_cfd_mesh_t  *mesh,
                 const ple_coord_t      point_coords[],
                 ple_lnum_t             n_point_ids,
                 const ple_lnum_t       point_ids[],
                 ple_lnum_t             location[],
                 float                  distance[])
{
  ple_lnum_t  i;

  const int stride = syr_cfd_mesh_n_vertices_element[mesh->element_type];

  if (mesh->n_elements == 0)
    return;

  if (mesh->element_type == SYR_CFD_TETRA)
    ple_error(__FILE__, __LINE__, 0,
              _("Locating volume elements closest to points not handled yet"));

  else if (mesh->element_type == SYR_CFD_TRIA) {

    for (i = 0; i < mesh->n_elements; i++)
      _locate_on_triangle_3d(i + 1,
                             mesh->vertex_num + i*stride,
                             mesh->vertex_coords,
                             point_coords,
                             n_point_ids,
                             point_ids,
                             (HUGE_VAL / 4.),
                             location,
                             distance);

  }
}

/*----------------------------------------------------------------------------
 * Find elements in a given section containing 2d points: updates the
 * location[] and distance[] arrays associated with a set of points
 * for points that are in an element of this section, or closer to one
 * than to previously encountered elements.
 *
 * parameters:
 *   mesh              <-- pointer to mesh representation structure
 *   tolerance         <-- associated tolerance
 *   point_coords      <-- point coordinates
 *   quadtree          <-- point quadtree
 *   points_in_extents <-- array for query of ids of points in extents
 *                         (size: quadtree->n_points, less usually needed)
 *   location          <-> number of element containing or closest to each
 *                         point (size: n_points)
 *   distance          <-> distance from point to element indicated by
 *                         location[]: < 0 if unlocated, 0 - 1 if inside,
 *                         > 1 if outside (size: n_points)
 *----------------------------------------------------------------------------*/

static void
_mesh_locate_2d(const syr_cfd_mesh_t  *mesh,
                double                 tolerance,
                const ple_coord_t      point_coords[],
                _quadtree_t           *quadtree,
                ple_lnum_t             points_in_extents[],
                ple_lnum_t             location[],
                float                  distance[])
{
  ple_lnum_t  i, j;
  double elt_extents[4];

  ple_lnum_t n_points_in_extents = 0;

  const int stride = syr_cfd_mesh_n_vertices_element[mesh->element_type];
  const int elt_dim = stride - 1; /* true for simplexes */

  if (mesh->n_elements == 0)
    return;

  /* Main loop on elements */

  for (i = 0; i < mesh->n_elements; i++) {

    int elt_initialized = 0;

    for (j = 0; j < stride; j++) {
      ple_lnum_t vertex_id = mesh->vertex_num[i*stride + j] - 1;
      _update_elt_extents(2,
                          vertex_id,
                          mesh->vertex_coords,
                          elt_extents,
                          &elt_initialized);
    }
    _elt_extents_finalize(2,
                          elt_dim,
                          tolerance,
                          elt_extents);

    _query_quadtree(elt_extents,
                    point_coords,
                    quadtree,
                    &n_points_in_extents,
                    points_in_extents);

    if (mesh->element_type == SYR_CFD_TRIA)
      _locate_in_triangle_2d(i + 1,
                             mesh->vertex_num + i*stride,
                             mesh->vertex_coords,
                             point_coords,
                             n_points_in_extents,
                             points_in_extents,
                             tolerance,
                             location,
                             distance);

    else if (mesh->element_type == SYR_CFD_EDGE)
      _locate_on_edge_2d(i + 1,
                         mesh->vertex_num + i*stride,
                         mesh->vertex_coords,
                         point_coords,
                         n_points_in_extents,
                         points_in_extents,
                         tolerance,
                         location,
                         distance);

  } /* End of loop on elements */
}

/*----------------------------------------------------------------------------
 * Find elements in a given section closest to 3d points: updates the
 * location[] and distance[] arrays associated with a set of points
 * for points that are in an element of this section, or closer to one
 * than to previously encountered elements.
 *
 * parameters:
 *   mesh              <-- pointer to mesh representation structure
 *   point_coords      <-- point coordinates
 *   n_point_ids       <-- number of points to locate
 *   point_id          <-- ids of points to locate
 *   location          <-> number of element containing or closest to each
 *                         point (size: n_points)
 *   distance          <-> distance from point to element indicated by
 *                         location[]: < 0 if unlocated, or absolute
 *                         distance to a line element (size: n_points)
 *----------------------------------------------------------------------------*/

static void
_mesh_closest_2d(const syr_cfd_mesh_t  *mesh,
                 const ple_coord_t      point_coords[],
                 ple_lnum_t             n_point_ids,
                 const ple_lnum_t       point_id[],
                 ple_lnum_t             location[],
                 float                  distance[])
{
  ple_lnum_t  i;

  const int stride = syr_cfd_mesh_n_vertices_element[mesh->element_type];

  /* Return immediately if nothing to do for this rank */

  if (mesh->n_elements == 0 || mesh->element_type != SYR_CFD_EDGE)
    return;

  /* Main loop on elements */

  for (i = 0; i < mesh->n_elements; i++) {

    /* Locate on edge */

    _locate_on_edge_2d(i + 1,
                       mesh->vertex_num + i*stride,
                       mesh->vertex_coords,
                       point_coords,
                       n_point_ids,
                       point_id,
                       -1.0,
                       location,
                       distance);

  } /* End of loop on elements */
}

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Compute extents of a mesh representation
 *
 * parameters:
 *   mesh          <-- pointer to mesh representation structure
 *   n_max_extents <-- maximum number of sub-extents (such as element extents)
 *                     to compute, or -1 to query
 *   tolerance     <-- addition to local extents of each element:
 *                     extent = base_extent * (1 + tolerance)
 *   extents       <-> extents associated with mesh:
 *                     x_min, y_min, ..., x_max, y_max, ... (size: 2*dim)
 *
 * returns:
 *   the number of extents computed
 *----------------------------------------------------------------------------*/

ple_lnum_t
syr_cfd_point_location_extents(const void  *mesh,
                               ple_lnum_t   n_max_extents,
                               double       tolerance,
                               double       extents[])
{
//const syr_cfd_mesh_t  *m = mesh;
  const syr_cfd_mesh_t  *m = (syr_cfd_mesh_t*)mesh; 

  ple_lnum_t retval = 0;

  if (m == NULL)
    return 0;

  /* In query mode, return maximum extents available */

  if (n_max_extents < 0)
    retval = m->n_elements;

  /* If n_max_extents < n_elements, we return global mesh extents */

  else if (n_max_extents > 0 && n_max_extents < m->n_elements) {

    int i;
    ple_lnum_t j, k, vertex_id;
    double elt_extents[6];

    int dim = m->dim;
    int stride = syr_cfd_mesh_n_vertices_element[m->element_type];
    int elt_dim = stride - 1; /* True for simplexes */

    /* initialize extents in case mesh is empty or dim < 3 */
    for (i = 0; i < dim; i++) {
      extents[i]       =  HUGE_VAL;
      extents[i + dim] = -HUGE_VAL;
    }

    /* Compute extents */

    for (j = 0; j < m->n_elements; j++) {

      int elt_initialized = 0;

      for (k = 0; k < stride; k++) {

        vertex_id = m->vertex_num[j*stride + k] - 1;

        _update_elt_extents(dim,
                            vertex_id,
                            m->vertex_coords,
                            elt_extents,
                            &elt_initialized);

      }

      _elt_extents_finalize(dim,
                            elt_dim,
                            tolerance,
                            elt_extents);
      _update_extents(dim, elt_extents, extents);
    }

    retval = 1;
  }

  /* If n_max_extents = n_elements, we return element extents */

  else if (n_max_extents == m->n_elements) {

    ple_lnum_t j, k, vertex_id;

    int dim = m->dim;
    int stride = syr_cfd_mesh_n_vertices_element[m->element_type];
    int elt_dim = stride - 1; /* True for simplexes */

    /* Compute extents */

    for (j = 0; j < m->n_elements; j++) {

      double *elt_extents = extents + dim*2*j;
      int elt_initialized = 0;

      for (k = 0; k < stride; k++) {

        vertex_id = m->vertex_num[j*stride + k] - 1;

        _update_elt_extents(dim,
                            vertex_id,
                            m->vertex_coords,
                            elt_extents,
                            &elt_initialized);

      }

      _elt_extents_finalize(dim,
                            elt_dim,
                            tolerance,
                            elt_extents);
    }

    retval = n_max_extents;
  }

  return retval;
}

/*----------------------------------------------------------------------------
 * Find elements in a given mesh containing points: updates the
 * location[] and distance[] arrays associated with a set of points
 * for points that are in an element of this mesh, or closer to one
 * than to previously encountered elements.
 *
 * parameters:
 *   mesh         <-- pointer to mesh representation structure
 *   tolerance    <-- associated tolerance
 *   n_points     <-- number of points to locate
 *   point_coords <-- point coordinates
 *   location     <-> number of element containing or closest to each
 *                    point (size: n_points)
 *   distance     <-> distance from point to element indicated by
 *                    location[]: < 0 if unlocated, 0 - 1 if inside,
 *                    and > 1 if outside a volume element, or absolute
 *                    distance to a surface element (size: n_points)
 *----------------------------------------------------------------------------*/

void
syr_cfd_point_location_contain(const void         *mesh,
                               double              tolerance,
                               ple_lnum_t          n_points,
                               const ple_coord_t   point_coords[],
                               ple_lnum_t          location[],
                               float               distance[])
{
  if (mesh == NULL)
    return;

  else {

    ple_lnum_t  *points_in_extents = NULL;

  //const syr_cfd_mesh_t  *m = mesh;
    const syr_cfd_mesh_t  *m = (syr_cfd_mesh_t*)mesh;


    /* Build point query list
       (max size: n_points, usually much less) */

    PLE_MALLOC(points_in_extents, n_points, ple_lnum_t);

    /* Use octree for 3d point location */

    if (m->dim == 3) {

      _octree_t  octree = _build_octree(n_points, point_coords);

      /* Locate for all sections */

      _mesh_locate_3d(m,
                      tolerance,
                      point_coords,
                      &octree,
                      points_in_extents,
                      location,
                      distance);

      _free_octree(&octree);
    }

    /* Use quadtree for 2d point location */

    else if (m->dim == 2) {

      _quadtree_t  quadtree = _build_quadtree(n_points, point_coords);

      _mesh_locate_2d(m,
                      tolerance,
                      point_coords,
                      &quadtree,
                      points_in_extents,
                      location,
                      distance);

      _free_quadtree(&quadtree);
    }

    PLE_FREE(points_in_extents);
  }
}

/*----------------------------------------------------------------------------
 * Find elements in a given mesh closest to points: updates the
 * location[] and distance[] arrays associated with a set of points
 * for points that are closer to an element of this mesh than to previously
 * encountered elements.
 *
 * This function currently only handles elements of lower dimension than
 * the spatial dimension.
 *
 * parameters:
 *   mesh         <-- pointer to mesh representation structure
 *   n_points     <-- number of points to locate
 *   point_coords <-- point coordinates
 *   location     <-> number of element containing or closest to each
 *                    point (size: n_points)
 *   distance     <-> distance from point to element indicated by
 *                    location[]: < 0 if unlocated, or absolute
 *                    distance to a surface element (size: n_points)
 *----------------------------------------------------------------------------*/

void
syr_cfd_point_location_closest(const void         *mesh,
                               ple_lnum_t          n_points,
                               const ple_coord_t   point_coords[],
                               ple_lnum_t          location[],
                               float               distance[])
{
  if (mesh == NULL)
    return;

  else {

    ple_lnum_t j;
    ple_lnum_t  *point_ids = NULL;

  //const syr_cfd_mesh_t  *m = mesh;
    const syr_cfd_mesh_t  *m = (syr_cfd_mesh_t*)mesh;

    PLE_MALLOC(point_ids, n_points, ple_lnum_t);
    for (j = 0; j < n_points; j++)
      point_ids[j] = j;

    /* Use brute force for closest 3d point location */

    if (m->dim == 3)
      _mesh_closest_3d(m,
                       point_coords,
                       n_points,
                       point_ids,
                       location,
                       distance);

    /* Use brute force for closest 2d point location */

    else if (m->dim == 2)
      _mesh_closest_2d(m,
                       point_coords,
                       n_points,
                       point_ids,
                       location,
                       distance);

    if (point_ids != NULL)
      PLE_FREE(point_ids);
  }
}

/*----------------------------------------------------------------------------*/

#undef _DOT_PRODUCT
#undef _MODULE
#undef _CROSS_PRODUCT


/*===========================================================================
 * Functions for coupling between Syrthes and Code_Saturne
 * AUTHORS  : J. Bonelle, Y Fournier
 *
 * Library: Syrthes                                   Copyright EDF 2008-2010
 *===========================================================================*/
