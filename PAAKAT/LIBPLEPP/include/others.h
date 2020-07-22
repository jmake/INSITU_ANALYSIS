//#include "build_octree.h"
  #include "cs_alya.h"
//=======================================================================||===//
//=======================================================================||===//

/*----------------------------------------------------------------------------
 *   elt_type    <-- type of element
 *   uvw[]       <-- parametric coordinates
 *   shapef[]    --> barycenter's coordinates
 *   deriv [][]  --> derivative of shape function
 *----------------------------------------------------------------------------*/

static void
_compute_shapef_3d(fvm_element_t  elt_type,
                   const double   uvw[3],
                   double         shapef[8],
                   double         deriv[8][3])

{
  switch (elt_type) 
  {
  case FVM_CELL_HEXA:
    shapef[0] = (1.0 - uvw[0]) * (1.0 - uvw[1]) * (1.0 - uvw[2]);
    shapef[1] = uvw[0] * (1.0 - uvw[1]) * (1.0 - uvw[2]);
    shapef[2] = uvw[0] * uvw[1] * (1.0 - uvw[2]);
    shapef[3] = (1.0 - uvw[0]) * uvw[1] * (1.0 - uvw[2]);
    shapef[4] = (1.0 - uvw[0]) * (1.0 - uvw[1]) * uvw[2];
    shapef[5] = uvw[0] * (1.0 - uvw[1]) * uvw[2];
    shapef[6] = uvw[0] * uvw[1] * uvw[2];
    shapef[7] = (1.0 - uvw[0]) * uvw[1] * uvw[2];

    if (deriv != NULL) {
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
    break;

  case FVM_CELL_PRISM:
    shapef[0] = (1.0 - uvw[0] - uvw[1]) * (1.0 - uvw[2]);
    shapef[1] = uvw[0] * (1.0 - uvw[2]);
    shapef[2] = uvw[1] * (1.0 - uvw[2]);
    shapef[3] = (1.0 - uvw[0] - uvw[1]) * uvw[2];
    shapef[4] = uvw[0] * uvw[2];
    shapef[5] = uvw[1] * uvw[2];

    if (deriv != NULL) {
      deriv[0][0] = -(1.0 - uvw[2]);
      deriv[0][1] = -(1.0 - uvw[2]);
      deriv[0][2] = -(1.0 - uvw[0] - uvw[1]);
      deriv[1][0] =  (1.0 - uvw[2]);
      deriv[1][1] =  0.0;
      deriv[1][2] = -uvw[0];
      deriv[2][0] =  0.0;
      deriv[2][1] =  (1.0 - uvw[2]);
      deriv[2][2] = -uvw[1];
      deriv[3][0] = -uvw[2];
      deriv[3][1] = -uvw[2];
      deriv[3][2] =  (1.0 - uvw[0] - uvw[1]);
      deriv[4][0] =  uvw[2];
      deriv[4][1] =  0.0;
      deriv[4][2] =  uvw[0];
      deriv[5][0] =  0.0;
      deriv[5][1] =  uvw[2];
      deriv[5][2] =  uvw[1];
    }
    break;

  case FVM_CELL_PYRAM:
    shapef[0] = (1.0 - uvw[0]) * (1.0 - uvw[1]) * (1.0 - uvw[2]);
    shapef[1] = uvw[0] * (1.0 - uvw[1]) * (1.0 - uvw[2]);
    shapef[2] = uvw[0] * uvw[1] * (1.0 - uvw[2]);
    shapef[3] = (1.0 - uvw[0]) * uvw[1] * (1.0 - uvw[2]);
    shapef[4] = uvw[2];

    if (deriv != NULL) {
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
      deriv[4][0] =  0.0;
      deriv[4][1] =  0.0;
      deriv[4][2] =  1.0;
    }
    break;

  /*
  default:

    bft_error(__FILE__, __LINE__, 0,
              _("_compute_shapef: unhandled element type %s\n"),
              fvm_element_type_name[elt_type]);
  */
  }

}


static int
_inverse_3x3(double  m[3][3],
             double  b[3],
             double  x[3])
{
  double det, det_inv, x0, x1, x2;

  det =   m[0][0]*(m[1][1]*m[2][2] - m[2][1]*m[1][2])
        - m[1][0]*(m[0][1]*m[2][2] - m[2][1]*m[0][2])
        + m[2][0]*(m[0][1]*m[1][2] - m[1][1]*m[0][2]);

  if(fabs(det) < _epsilon_denom) 
  {
    return 1;
/*
  if (CS_ABS(det) < _epsilon_denom)
    return 1;
*/
  }
  else
    det_inv = 1./det;

  /* Use local variables to ensure no aliasing */

  x0 = (  b[0]*(m[1][1]*m[2][2] - m[2][1]*m[1][2])
        - b[1]*(m[0][1]*m[2][2] - m[2][1]*m[0][2])
        + b[2]*(m[0][1]*m[1][2] - m[1][1]*m[0][2])) * det_inv;

  x1 = (  m[0][0]*(b[1]*m[2][2] - b[2]*m[1][2])
        - m[1][0]*(b[0]*m[2][2] - b[2]*m[0][2])
        + m[2][0]*(b[0]*m[1][2] - b[1]*m[0][2])) * det_inv;

  x2 = (  m[0][0]*(m[1][1]*b[2] - m[2][1]*b[1])
        - m[1][0]*(m[0][1]*b[2] - m[2][1]*b[0])
        + m[2][0]*(m[0][1]*b[1] - m[1][1]*b[0])) * det_inv;

  /* Copy local variables to output */

  x[0] = x0; x[1] = x1; x[2] = x2;

  return 0;
}


static int
_compute_uvw(fvm_element_t      elt_type,
             const ple_coord_t   point_coords[],
             double             vertex_coords[8][3],
             double             tolerance,
             double             uvw[3])

{
  int i, j, n_elt_vertices, iter;
  int max_iter = 20;
  double dist;
  double a[3][3], b[3], x[3], shapef[8], dw[8][3];

  n_elt_vertices = fvm_nodal_n_vertices_element[elt_type];

  assert(   elt_type == FVM_CELL_HEXA
         || elt_type == FVM_CELL_PRISM
         || elt_type == FVM_CELL_PYRAM);

  // Use Newton-method to determine parametric coordinates and shape function

  for (i = 0; i < 3; i++)
    uvw[i] = 0.5;

  for (iter = 0; iter < max_iter; iter++) {

    _compute_shapef_3d(elt_type, uvw, shapef, dw);

    b[0] = - point_coords[0];
    b[1] = - point_coords[1];
    b[2] = - point_coords[2];

    for (i = 0; i < 3; i++) {
      for (j = 0; j < 3; j++)
        a[i][j] = 0.0;
    }

    for (i = 0; i < n_elt_vertices; i++) {

      b[0] += (shapef[i] * vertex_coords[i][0]);
      b[1] += (shapef[i] * vertex_coords[i][1]);
      b[2] += (shapef[i] * vertex_coords[i][2]);

      for (j = 0; j < 3; j++) {
        a[0][j] -= (dw[i][j] * vertex_coords[i][0]);
        a[1][j] -= (dw[i][j] * vertex_coords[i][1]);
        a[2][j] -= (dw[i][j] * vertex_coords[i][2]);
      }

    }

    if (_inverse_3x3(a, b, x))
      return 0;

    dist = 0.0;

    for (i = 0; i < 3; i++) {
      dist += x[i] * x[i];
      uvw[i] += x[i];
    }

    if (dist <= (tolerance * tolerance))
      return 1;

  }

  return 0;
}


static void
_locate_in_cell_3d(ple_lnum_t          elt_num,
                   fvm_element_t       elt_type,
                   const ple_lnum_t    element_vertex_num[],
                   const ple_lnum_t   *parent_vertex_num,
                   const ple_coord_t   vertex_coords[],
                   const ple_coord_t   point_coords[],
                   ple_lnum_t          n_points_in_extents,
                   const ple_lnum_t    points_in_extents[],
                   double              tolerance,
                   ple_lnum_t          location[],
                   float               distance[])
{
  int i, j, k, n_vertices;
  ple_lnum_t coord_idx, vertex_id;

  double uvw[3], dist, shapef[8],max_dist;
  double  _vertex_coords[8][3];

  n_vertices = fvm_nodal_n_vertices_element[elt_type];

  /* Initialize local element coordinates copy */

  for (vertex_id = 0; vertex_id < n_vertices; vertex_id++) {

    if (parent_vertex_num == NULL)
      coord_idx = element_vertex_num[vertex_id] -1;
    else
      coord_idx = parent_vertex_num[element_vertex_num[vertex_id] - 1] - 1;

    for (j = 0; j < 3; j++)
      _vertex_coords[vertex_id][j] = vertex_coords[(coord_idx * 3) + j];

  }

  /* Shape functions may be computed directly with tetrahedra */

  if (elt_type == FVM_CELL_TETRA)
  {
/*
    _locate_in_tetra(elt_num,
                     _vertex_coords,
                     point_coords,
                     n_points_in_extents,
                     points_in_extents,
                     tolerance,
                     location,
                     distance);
*/

      _locate_in_tetra(elt_num,
                       element_vertex_num,
                       vertex_coords,
                       point_coords,
                       n_points_in_extents,
                       points_in_extents,
                       tolerance,
                       location,
                       distance);


  /* For cell shapes other than tetrahedra, find shape functions iteratively */
  }
  else {

    for (k = 0; k < n_points_in_extents; k++) {

      i = points_in_extents[k];

      if (_compute_uvw(elt_type,
                       point_coords + 3*i,
                       _vertex_coords,
                       tolerance,
                       uvw)) {

        max_dist = -1.0;

        /* For hexahedra, no need to compute shape functions, as
           the 3 parametric coordinates are simpler to use */

        if (elt_type == FVM_CELL_HEXA) {

          for (j = 0; j < 3; j++){

            dist = 2.*fabs(uvw[j] - 0.5);
            //dist = 2.*CS_ABS(uvw[j] - 0.5);

            if (max_dist < dist)
              max_dist = dist;
          }

        }

        /* For pyramids ands prisms, we need to compute shape functions */

        else {

          _compute_shapef_3d(elt_type, uvw, shapef, NULL);

          for (j = 0; j < n_vertices; j++){

            dist = 2.*fabs(shapef[j] - 0.5);
            //dist = 2.*CS_ABS(shapef[j] - 0.5);

          if (max_dist < dist)
            max_dist = dist;
          }

        }

        /* For all element types, update location and distance arrays */

        if (   (max_dist > -0.5 && max_dist < (1. + 2.*tolerance))
            && (max_dist < distance[i] || distance[i] < 0)) {
          location[i] = elt_num;
          distance[i] = max_dist;
        }

      }

    } /* End of loop on points in extents */

  }

}


/*----------------------------------------------------------------------------
 * Structure defining a mesh section
 *----------------------------------------------------------------------------*/

typedef int      cs_lnum_t;     /* Local integer index or number */
typedef double   cs_coord_t;    /* Real number (coordinate value) */


typedef struct _fvm_nodal_section_t 
{
  /* Basic information */
  /*-------------------*/

  int         entity_dim;          /* Entity dimension */

  cs_lnum_t   n_elements;          /* Number of elements */

  fvm_element_t  type;             /* Element types */

  /* Connectivity */
  /*--------------*/

  size_t      connectivity_size;   /* Size of vertex_num array;
                                      for strided elements:
                                       (n_elements * stride)
                                      for polygons:
                                       (vertex_index[n_elements])
                                      for polyhedra:
                                       (vertex_index[n_faces]) */

  int         stride;              /* Element size for regular elements
                                      (0 for polygons and polyhedra) */

  cs_lnum_t   n_faces;             /* Number of faces defining polyhedra */

  /* Pointers to connectivity arrays, which may be shared */

  const cs_lnum_t   *face_index;   /* polyhedron -> faces index (O to n-1);
                                      size: n_elements + 1 */
  const cs_lnum_t   *face_num;     /* polyhedron -> face numbers (1 to n, signed,
                                      > 0 for outwards pointing face normal
                                      < 0 for inwards pointing face normal);
                                      size: face_index[n_elements] */

  const cs_lnum_t   *vertex_index; /* polygon face -> vertices index (O to n-1);
                                      size: n_faces + 1 */

  const cs_lnum_t   *vertex_num;   /* vertex numbers (1 to n);
                                      size: connectivity_size */

  /* Pointers to local connectivity arrays, if owner */

  cs_lnum_t   *_face_index;        /* face_index if owner, NULL if shared */
  cs_lnum_t   *_face_num;          /* face_num if owner, NULL if shared */
  cs_lnum_t   *_vertex_index;      /* vertex_index if owner, NULL if shared */
  cs_lnum_t   *_vertex_num;        /* vertex numbers if owner, NULL if shared */

  /* Pointers to group class ids, if present */

  int         *gc_id;              /* Group class id, NULL if implicit 0 */

  /* Auxiliary structure used to define subdivision of elements into
     simpler element types (usually polygons to triangles and
     polyhedra to tetrahedra and pyramids) */

//  fvm_tesselation_t  *tesselation;

  /* Numbering */
  /*-----------*/

  const cs_lnum_t   *parent_element_num; /* Local numbers (1 to n) of local
                                            elements in the parent mesh,
                                            associated with the section's
                                            elements.

                                            This array is necessary to redis-
                                            tribute output fields when the
                                            section has been either associated
                                            with an unsorted mixed mesh,
                                            renumbered, or is associated with a
                                            subset of a more complete mesh,
                                            such as a clip plane. When used for
                                            a subset, it also defines the lists
                                            of elements of the parent mesh
                                            belonging to that subset.

                                            This array is present only when non
                                            "trivial" (i.e. not 1, 2, ..., n). */

  cs_lnum_t     *_parent_element_num;    /* pointer to parent_element_num if
                                            owner, NULL otherwise */

//  fvm_io_num_t  *global_element_num;     /* Global element numbers */

} fvm_nodal_section_t;



/*----------------------------------------------------------------------------
 * Find elements in a given section containing 3d points: updates the
 * location[] and distance[] arrays associated with a set of points
 * for points that are in an element of this section, or closer to one
 * than to previously encountered elements.
 *
 * parameters:
 *   this_section      <-- pointer to mesh section representation structure
 *   parent_vertex_num <-- pointer to parent vertex numbers (or NULL)
 *   vertex_coords     <-- pointer to vertex coordinates
 *   tolerance         <-- associated tolerance
 *   base_element_num  <-- < 0 for location relative to parent element numbers,
 *                         number of elements in preceding sections of same
 *                         element dimension + 1 otherwise
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
_nodal_section_locate_3d(const fvm_nodal_section_t  *this_section,
                         const cs_lnum_t            *parent_vertex_num,
                         const cs_coord_t            vertex_coords[],
                         double                      tolerance,
                         cs_lnum_t                   base_element_num,
                         const cs_coord_t            point_coords[],
                         _octree_t                  *octree,
                         cs_lnum_t                   points_in_extents[],
                         cs_lnum_t                   location[],
                         float                       distance[])
{
  cs_lnum_t   i, j, vertex_id, elt_num, triangle_vertices[6];
  int n_triangles;
  double elt_extents[6];

  cs_lnum_t n_points_in_extents = 0;

  // If section contains polyhedra 
/*
  if (this_section->type == FVM_CELL_POLY)

    _polyhedra_section_locate(this_section,
                              parent_vertex_num,
                              vertex_coords,
                              tolerance,
                              base_element_num,
                              point_coords,
                              octree,
                              points_in_extents,
                              location,
                              distance);

  // If section contains polygons 

  else if (this_section->type == FVM_FACE_POLY)

    _polygons_section_locate_3d(this_section,
                                parent_vertex_num,
                                vertex_coords,
                                tolerance,
                                base_element_num,
                                point_coords,
                                octree,
                                points_in_extents,
                                location,
                                distance);

  // If section contains regular elements 

  else {
*/
{
    for (i = 0; i < this_section->n_elements; i++) 
    {

      bool elt_initialized = false;

      if (base_element_num < 0) 
      {
        if (this_section->parent_element_num != NULL)
          elt_num = this_section->parent_element_num[i];
        else
          elt_num = i + 1;
      }
      else
        elt_num = base_element_num + i;

      for (j = 0; j < this_section->stride; j++) 
      {
        vertex_id = this_section->vertex_num[i*this_section->stride + j] - 1;
/*
        _update_elt_extents(3,
                            vertex_id,
                            parent_vertex_num,
                            vertex_coords,
                            elt_extents,
                            &elt_initialized);
*/
      }

      _elt_extents_finalize(3,
                            this_section->entity_dim,
                            tolerance,
                            elt_extents);

      _query_octree(elt_extents,
                    point_coords,
                    octree,
                    &n_points_in_extents,
                    points_in_extents);

      if (this_section->entity_dim == 3)

        _locate_in_cell_3d(elt_num,
                           this_section->type,
                           this_section->vertex_num + i*this_section->stride,
                           parent_vertex_num,
                           vertex_coords,
                           point_coords,
                           n_points_in_extents,
                           points_in_extents,
                           tolerance,
                           location,
                           distance);

      else if (this_section->entity_dim == 2) 
      {
/*
        if (this_section->type == FVM_FACE_QUAD)

          n_triangles = fvm_triangulate_quadrangle(3,
                                                   vertex_coords,
                                                   parent_vertex_num,
                                                   (  this_section->vertex_num
                                                    + i*this_section->stride),
                                                   triangle_vertices);

        else {

          assert(this_section->type == FVM_FACE_TRIA);

          n_triangles = 1;
          for (j = 0; j < 3; j++)
            triangle_vertices[j]
              = this_section->vertex_num[i*this_section->stride + j];


        }

        _locate_on_triangles_3d(elt_num,
                                n_triangles,
                                triangle_vertices,
                                parent_vertex_num,
                                vertex_coords,
                                point_coords,
                                n_points_in_extents,
                                points_in_extents,
                                tolerance,
                                location,
                                distance);
      }
*/
      }
      else if (this_section->entity_dim == 1) 
      {
/*
        assert(this_section->type == FVM_EDGE);

        _locate_on_edge_3d(elt_num,
                           this_section->vertex_num + i*this_section->stride,
                           parent_vertex_num,
                           vertex_coords,
                           point_coords,
                           n_points_in_extents,
                           points_in_extents,
                           tolerance,
                           location,
                           distance);
*/
      }

    }

  }
}

