/**
  @class CommDom

  @brief  communication domains [safe fluids interchange]  
  @author J. M. Zavala Aké
  @see      commdom.hpp  

  2014Jan29, __mpi_reduce__, __mpi_gather_double__
*/

#include <iostream> 
#include <cstring>
#include <vector>
#include <map>
using namespace std;


#ifndef COMMDOM_H
/**  @cond */
#define COMMDOM_H 


//=======================================================================||===//
//=======================================================================||===//
#include "outside_ple.h"

//#define __DIM__ 3
typedef ple_coord_t   REAL;
typedef ple_lnum_t    INTEGER;
/*
  ple_coord_t        *vertex_coords
  ple_lnum_t         *vertex_num
*/
/**  @endcond */


//=======================================================================||===//
//==============================================================| CommDom |===//
class CommDom
{
  public: 
    CommDom(); 
   ~CommDom();
   
    void init(); 
    void end(); 

   /*!
    *  \brief A useful method.
    *  \param level an integer setting how useful to be
    *  \return Output that is extra useful
    * 
    *   This method does unbelievably useful things.  
    *   And returns exceptionally useful results.
    *   Use it everyday with good health.
    */
    void set_app_type(string, int=-1); 
    void set_app_name(string); 

    void __genesis__(); 

    void set_world_comm(MPI_Comm); 

/*!
 * \brief
 *
 * \param[in]  cosa01
 * \param[out] cosa02
 * \return NULL
 * \sa
 * \note
 * \warning
 */
    void name_to_id(int*, int*, MPI_Comm*); 
    
    void set_applications(string); 
    void set_iteration(string); 
    void create_commij(MPI_Comm); 
    void create_interaction(); 
    void __create_interaction__(); 


    // Python 
    MPI_Comm set_mpi_comms(); 
    MPI_Comm get_mpi_commij(string);
    int      get_mpi_commij_size(); 
    void     sendrecv_int(int*, int, int*, int, MPI_Comm, MPI_Comm); 


    // COMMS 
    void get_commij_size(int*); 
    void print_commij_names(); 
    void get_commij(string, MPI_Comm*); 


    // MPI
    void __mpi_sendrecv_int__(int*, int, int*, int, MPI_Comm, MPI_Comm); 
    void __mpi_sendrecv_real__(double*, int, double*, int, MPI_Comm, MPI_Comm);
    void __mpi_sendrecv_char__(char*, int, char*, int, MPI_Comm, MPI_Comm);
    void __mpi_bcast_real__( double*, int, MPI_Comm, MPI_Comm);

    void __mpi_reduce__(double*, double*, int, MPI_Comm, MPI_Comm, MPI_Op, MPI_Datatype);
    void __mpi_reduce_min_int__(int*, int*, MPI_Comm, MPI_Comm); 
    void __mpi_reduce_min_double__(double*, double*, MPI_Comm, MPI_Comm); 
    void __mpi_reduce_sum_double__(double*, double*, MPI_Comm, MPI_Comm); 
    int  __mpi_get_roots__(MPI_Comm, MPI_Comm, int*, int*);
    void __mpi_gather_double__(    double*, double*, int*, MPI_Comm, MPI_Comm);


    // LOCATOR    
    void __locator_init_meshes__();
    void __locator_end_meshes__();

    void __locator_set_meshi__(int, int, double*, int*); 
    void __locator_set_meshj__(int, double*); 
    void __locator_set_meshij__(ple_locator_t*, syr_cfd_mesh_t*, int, int, int, const int*, double*, float*);
    void __locator_set_mesh__(ple_locator_t*); 

    void locator_create(MPI_Comm, int, int, double=1e-3);
    void locator_create2(MPI_Comm, MPI_Comm, double, int); // 2019FEB21  
    void locator_destroy(); 

    ple_locator_t*  __locator_get__(); 
    void __locator_get_interior_list__(int*); 
    void __locator_get_dist_locations__(int*);
    void __locator_get_dist_coords__(double*);
    void __locator_get_dist_props__(double*); //MATIAS
    void __locator_exchange_double_scalar__(double*, double*, int=1); 
    void __locator_send_double_scalar__(double *, int=1); 
    void __locator_recv_double_scalar__(double *, int=1); 


    void locator_set_mesh(   int, int, double*, int*,       int, double*); 
    void locator_set_cs_mesh(int, int, double*, int*, int*, int, double*, int, double*, int);  // 2016Mar22 

//    void __alya2cs_type__(int*, int*, int); 
//    void __get_alya2cs_type__(int, int*); 


    //void exchange_point_var_dsend(double*, int, int);
    //void exchange_point_var_drecv(double*, int, int);
    int get_n_dist_points(); 
    int get_n_interior(); 
    void save_dist_coords(int, MPI_Comm=MPI_COMM_NULL); 
    void exchange(int=-1); 


    // SYRTHES
    void __syr_cfd_coupling_face_vars__();
    void _send_nodal_var(int=1); 
    void _send_nodal_var_vec();
    void _send_elt_var(); 


    // AUX 
    void __set_create__(string, MPI_Comm, int*); 
    void __set_get_info__(); 
    void __print_comm__(MPI_Comm);
    void __error__(bool, string);
    void __strncmp__(const char*, const char*, int*, int*); 
    const char* string_f2c(const char*, int);
    void __to_vtk__(string, int, const double*);

    int         __hexa_interpolation__(double*, int*, double*, double*); // 2019Mar06 
    void       __tetra_interpolation__(double*, int*, double*, double*);
    void __triangle_2d_interpolation__(double*, int*, double*, double*); // 2016Mar22  
    void      __simple_interpolation__(double*, int*, double*, double*); // 2016Mar22  

    void __set_argvs__(string);
    void __get_argvs__(string&);
    void __analyse_argvs__(string);

void __coupling_sync_apps__(int, int, int*, double*); 
//void __coupling_sync_apps__(int, int, int*, double*, ple_coupling_mpi_set_t); 
//int _cs_coupling_sync_flag;
void __set_synchronization_flags();
const char* get_synch_flags(int);
//ple_coupling_mpi_set_t *_cs_glob_coupling_mpi_app_world;

    // used by fortran wrapper 
    vector<string> __split_string__(const string&, const string&, const bool=true);
    string         __get_app_name__(); 
    string         __get_app_type__(); 
    int            __get_n_apps__(); 
    int            __get_n_types__(); 
    int            __get_friends__(string); 
    MPI_Comm       __get_world_comm__(); 

#if SATURNE  // 2016MAY04 
    void __alya2cs_type__(int*, int*, int);
    void __get_alya2cs_type__(int, int*);
  //void  __create_nodal_t__(REAL* = NULL, int=0, INTEGER* =NULL, INTEGER* =NULL, int=0);
    void  __create_nodal_t__(void* =NULL, int = 0, int* =NULL, int* =NULL, int=0);
    void __interpolation__(int, double, double*,  int*, double*, double*);
    void __get_parent__(int*, int);
#endif 

  private:
    string                           IAM; 
    vector< vector<int> >            CplngIJ;

    map<string, int>                 WhoIam; 
    map<string, int>::iterator       WhoIamIT, Whoi, Whoj;

    map<string, MPI_Comm>            CommIJ;
    map<string, MPI_Comm>::iterator  CommIT;   

    map<string, vector<int> >            RootIJ; 
    map<string, vector<int> >::iterator  RootIT;

    MPI_Comm                        world_comm; 
    MPI_Comm                        quitze_comm;
    MPI_Comm                        commij; 

    
    string                          app_namej;
    string                          app_namei; 
    string                          app_type; 
    string                          app_name; 

    int                             world_size,  world_rank;
    int                             local_size,  local_rank;

    int                             n_apps, app_id; 
    int                             n_types; 


    // COUPLING 
    ple_coupling_mpi_set_t               *cplng; 
    vector<ple_coupling_mpi_set_info_t>   Info;


    // LOCATOR 
    ple_locator_t              *locatorj; 
    //vector<ple_locator_t>     LocatorIJ;
    
    vector<syr_cfd_mesh_t>      LocalMesh;
    int                         n_meshes; 

    int                         n_coords_dist; 
    REAL                       *coords_dist; 
    float                      *distance_dist;

    REAL                       *props_dist;

    int                         n_locs; 

    vector<string>              argvs;
    string                      argv_found; 

    map<int, string>            Synch;
    map<int, string>::iterator  SynchIT;   

    int *parent_element_num; // saturne

    map<int,int>            Alya2Cs;
    map<int,int>::iterator  Alya2CsIT;

   int dimension;    // 2016Mar22
   int dimension_props;

   int element_type ; //2019FEB21  
}; 


#endif // COMMDOM_H
