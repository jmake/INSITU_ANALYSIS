#include "commdom_wrapper.h"
#include "commdom.hpp"
#include <iterator>

//=======================================================================||===//
//=======================================================================||===//
/**
  @page pagename05 Cwrapper
  @brief  communication domains [safe fluids interchange]  
  @author J. M. Zavala Aké
  @see      commdom_wrapper.cpp  
  @include  commdom_wrapper.cpp  
*/


static CommDom* ptr_class = NULL; 

extern "C" 
{
  void 
  commdom_create_()
  {
    ptr_class = new CommDom(); //[num]; 
    ptr_class->init();
  }

  void 
  commdom_delete_()
  {
    delete ptr_class; 
    //ptr_class = NULL; 
  } 


  void 
  commdom_set_argvs_(const char* ftype, int* ntype)
  {
    std::string  ctype( ptr_class->string_f2c(ftype, ntype[0]) );
    ptr_class->__set_argvs__(ctype);
  } 


  void 
  commdom_get_argvs_(char* ftype)
  {
    std::string  ctype; 

    ctype = "";
    ptr_class->__get_argvs__(ctype);
    if(ctype.size() > 0) strcpy(ftype, ctype.c_str()); 
    //cout<<"--> \'"<< ctype <<"\' \n\n"; 
  } 


  void 
  commdom_analyse_argvs_(const char* ftype, int* ntype)
  {
    std::string  ctype( ptr_class->string_f2c(ftype, ntype[0]) );
    ptr_class->__analyse_argvs__(ctype); 
  } 


  void 
  commdom_set_names_(const char* ftype, int* ntype, const char* fname, int* nname)
  {
    std::string  ctype( ptr_class->string_f2c(ftype, ntype[0]) );
    std::string  cname( ptr_class->string_f2c(fname, nname[0]) );

    const std::string     to_split( cname ); 
    const vector<string>  words = ptr_class->__split_string__(to_split, "/"); 
    if(words.size()>1)    cname = words[words.size()-1]; 

    ptr_class->set_app_name(cname);
    ptr_class->set_app_type(ctype);
/*
    ptr_class->__genesis__();
*/
  }


  void 
  commdom_create_commij_(int* fcomm, int* lcomm)
  {
    int app_id = -1;
    int n_apps = -1; 
    MPI_Comm  world_comm = MPI_Comm_f2c( fcomm[0] );
    ptr_class->set_world_comm(world_comm);
    // NOTE: 
    // there exit the posibility that 
    // will be changed 'world_comm'

    MPI_Comm  new_world_comm = MPI_COMM_NULL; 
    new_world_comm = ptr_class->__get_world_comm__();
    fcomm[0] = MPI_Comm_c2f( new_world_comm );


    MPI_Comm  local_comm = MPI_COMM_NULL; 
    ptr_class->name_to_id(&app_id, &n_apps, &local_comm);

    //ptr_class->__print_comm__(local_comm); 
    //ptr_class->__print_comm__(world_comm); 
    //ptr_class->print_commij_names(); 

    ptr_class->__create_interaction__(); 
    ptr_class->create_commij(local_comm);

    lcomm[0] = MPI_Comm_c2f(local_comm);

  } 


  void 
  commdom_get_n_types_(int* n_types)
  {
    n_types[0] = ptr_class->__get_n_types__(); 
  }


  void 
  commdom_get_n_apps_(int* n_apps)
  {
    n_apps[0] = ptr_class->__get_n_apps__(); 
  }


  void
  commdom_my_surname_(const char* fname, int* nname, int* ok)
  {
    std::string    cname( ptr_class->string_f2c(fname, nname[0]) );
    std::string  surname( ptr_class->__get_app_type__() );

    int size = surname.size(); 

    ptr_class->__strncmp__(cname.c_str(), surname.c_str(), &size, ok); 
  }


  void
  commdom_who_iam_(const char* fname, int* nname, int* ok)
  {
    std::string  cname( ptr_class->string_f2c(fname, nname[0]) );
    std::string    iam( ptr_class->__get_app_name__() );

    int size = iam.size(); 
    ptr_class->__strncmp__(cname.c_str(), iam.c_str(), &size, ok); 

  }


  void
  commdom_who_areyou_(const char* fname, int* nname, int* ok)
  {
    std::string  cname( ptr_class->string_f2c(fname, nname[0]) );
    ok[0] = ptr_class->__get_friends__(cname); 
  }


  void 
  commdom_get_commij_size_(int* commij_size)
  {
    ptr_class->get_commij_size(commij_size); 
  }


  void 
  commdom_get_commij_(char* fnamej, int* nnamej, int* fcommij)
  {
    std::string cnamej( ptr_class->string_f2c(fnamej, nnamej[0]) );
    
    MPI_Comm commij = MPI_COMM_NULL; 
    ptr_class->get_commij(cnamej, &commij); 
    fcommij[0] = MPI_Comm_c2f(commij);
  }


  void commdom_sendrecv_int_(int* send, int* n_send, int* recv, int* n_recv, int* fcommi, int* fcommij)
  {
    MPI_Comm commi  = MPI_Comm_f2c( fcommi [0]);
    MPI_Comm commij = MPI_Comm_f2c( fcommij[0]);

    ptr_class->__mpi_sendrecv_int__(send,  n_send[0], 
                                    recv,  n_recv[0], 
                                    commi, commij);
  }

  void commdom_sendrecv_real_(double* send, int* n_send, double* recv, int* n_recv, int* fcommi, int* fcommij)
  {
    MPI_Comm commi  = MPI_Comm_f2c( fcommi [0]);
    MPI_Comm commij = MPI_Comm_f2c( fcommij[0]);

    ptr_class->__mpi_sendrecv_real__(send,  n_send[0],
                                     recv,  n_recv[0],
                                     commi, commij);
  }


  void commdom_sendrecv_char_(char* send, int* n_send, char* recv, int* n_recv, int* fcommi, int* fcommij)
  {
    MPI_Comm commi  = MPI_Comm_f2c( fcommi [0]);
    MPI_Comm commij = MPI_Comm_f2c( fcommij[0]);

    std::string csend( ptr_class->string_f2c(send, n_send[0]) );
  //std::string crecv( ptr_class->string_f2c(recv, n_recv[0]) );

    char *aux_char = new char[ csend.size() ];  
    memcpy( aux_char, csend.c_str(), csend.size() ); 

    ptr_class->__mpi_sendrecv_char__(aux_char, csend.size(),
                                         recv, n_recv[0],
                                        commi, commij);
    delete[] aux_char; 
  }


  void commdom_bcast_real_(double* send, int* n_send, int* fcommi, int* fcommij)
  {
    MPI_Comm commi  = MPI_Comm_f2c( fcommi [0]);
    MPI_Comm commij = MPI_Comm_f2c( fcommij[0]);

    ptr_class->__mpi_bcast_real__(send, n_send[0], commi, commij);

  }


  void commdom_reduce_min_int_(int* send, int* recv, int* fcommi, int* fcommij)
  {
    MPI_Comm commi  = MPI_Comm_f2c( fcommi [0]);
    MPI_Comm commij = MPI_Comm_f2c( fcommij[0]);

    ptr_class->__mpi_reduce_min_int__(send, recv, commi, commij); 
  }


  void commdom_reduce_min_real_(double* send, double* recv, int* fcommi, int* fcommij)
  {
    MPI_Comm commi  = MPI_Comm_f2c( fcommi [0]);
    MPI_Comm commij = MPI_Comm_f2c( fcommij[0]);

    ptr_class->__mpi_reduce_min_double__(send, recv, commi, commij); 
  }


  void commdom_reduce_sum_real_(double* send, double* recv, int* fcommi, int* fcommij)
  {
    MPI_Comm commi  = MPI_Comm_f2c( fcommi [0]);
    MPI_Comm commij = MPI_Comm_f2c( fcommij[0]);

    ptr_class->__mpi_reduce_sum_double__(send, recv, commi, commij); 
  }


  void commdom_gather_real_(double* send, double* recv, int* size, int* fcommi, int* fcommij)
  {
    MPI_Comm commi  = MPI_Comm_f2c( fcommi [0]);
    MPI_Comm commij = MPI_Comm_f2c( fcommij[0]);

    ptr_class->__mpi_gather_double__(send, recv, size, commi, commij);
  }


  void commdom_sync_apps_(int* i_time, int* n_time, double* time) 
  {
    int flags = 0;     
    ptr_class->__coupling_sync_apps__(flags, i_time[0], n_time, time); 
  }

}
//=======================================================================||===//


//=======================================================================||===//
extern "C" 
{
  void commdom_locator_destroy_()
  {
    ptr_class->locator_destroy();
  }


  void commdom_locator_create_(int* fcommi, int* fcommij) 
  {
  }


  void commdom_locator_create2_(int* fcommi, int* fcommij, double* tol, int* element_type) 
  {
    MPI_Comm commi  = MPI_Comm_f2c( fcommi [0]);
    MPI_Comm commij = MPI_Comm_f2c( fcommij[0]);

    ptr_class->locator_create2(commi, commij, tol[0], element_type[0]);
  }


  void commdom_locator_set_mesh_(
                                int*        n_vertices_i, 
                                int*        n_elements_i, 
                                double*     vertex_coords_i, 
                                int*        vertex_num_i, 
                                int*        n_vertices_j, 
                                double*     vertex_coords_j
                                )
  {

    ptr_class->locator_set_mesh(n_vertices_i[0], n_elements_i[0], 
                                vertex_coords_i, vertex_num_i, 
                                n_vertices_j[0], vertex_coords_j);
  }


  void commdom_locator_set_cs_mesh_(
                                   int*        n_vertices_i,
                                   int*        n_elements_i,
                                   double*     vertex_coords_i,
                                   int*        vertex_num_i,
                                   int*        vertex_type_i,
                                   int*        n_vertices_j,
                                   double*     vertex_coords_j, 
                                   int*        dim,               // < 2016Mar22 
                                   double*     vertex_props_j,
                                   int*        dim_props 
                                   )
  {

    ptr_class->locator_set_cs_mesh(n_vertices_i[0], n_elements_i[0],
                                   vertex_coords_i, vertex_num_i, vertex_type_i,
                                   n_vertices_j[0], vertex_coords_j, dim[0], vertex_props_j, dim_props[0]);
  }



  void commdom_locator_get_n_dist_points_(int* size)
  { 
    size[0] = ptr_class->get_n_dist_points();
  }


  void commdom_locator_get_n_interior_(int* size)
  { 
    size[0] = ptr_class->get_n_interior();
  }


  void commdom_locator_get_get_n_intersects_(int* size)
  {
  // size[0] = ptr_class->get_n_intersects(); SEE: Thirdparties/libple/PLEPP/Wrappers/commdom_wrapper.cpp 
  }


  void commdom_locator_get_interior_list_(int* ids)
  { 

    ptr_class->__locator_get_interior_list__(ids); 

  }


  void commdom_locator_get_dist_locations_(int* ids)
  {
    ptr_class->__locator_get_dist_locations__(ids);
  }


  void commdom_locator_get_dist_coords_(double* coords)
  {
    ptr_class->__locator_get_dist_coords__(coords);
  }


  void commdom_locator_get_dist_props_(double* props)
  {
    ptr_class->__locator_get_dist_props__(props);
  }


  void commdom_locator_exchange_double_scalar_(double *var_ij, double *var_ji)
  {
    ptr_class->__locator_exchange_double_scalar__(var_ij, var_ji); 
  }


  void commdom_locator_exchange_double_stride_(double *var_ij, double *var_ji, int* stride)
  {
    ptr_class->__locator_exchange_double_scalar__(var_ij, var_ji, stride[0]); 
  }


  void commdom_locator_send_double_stride_(double *var_ij, int* stride)
  {
    ptr_class->__locator_send_double_scalar__(var_ij, stride[0]); 
  }


  void commdom_locator_recv_double_stride_(double *var_ji, int* stride)
  {
    ptr_class->__locator_recv_double_scalar__(var_ji, stride[0]); 
  }


  void commdom_locator_save_dist_coords_(int* id, int* fcommij)
  {

    if( (id[0]>=0)&&(fcommij[0]>0) )
    {
      MPI_Comm commij = MPI_Comm_f2c( fcommij[0] );
      ptr_class->save_dist_coords(id[0], commij); 
    }
    else 
    if(id[0]>=0)
    {
      ptr_class->save_dist_coords(id[0]); 
    }

    ptr_class->exchange();

  }


  void commdom_locator_tetra_interpolation_(double* coords_i, int* vertices_i, 
                                            double* point_j, double* tetrahedral_coords_j)
  {
    ptr_class->__tetra_interpolation__(coords_i, vertices_i, point_j, tetrahedral_coords_j);
  }


  void commdom_locator_interpolation_(
                                      int*   alya_type,  
                                      double*       tol,  
                                      double* coords_i, int* vertices_i,
                                      double* point_j,  double* shapef_j)
  {
#if SATURNE
    int cs_type[1];

    ptr_class->__get_alya2cs_type__(alya_type[0], cs_type); 
/*
cout<<"cs_type:"<< cs_type[0] <<" ";  
cout<<"alya_type:"<< alya_type[0] <<" ";
cout<<"\n";
*/
    ptr_class->__interpolation__(   cs_type[0],
                                        tol[0],
                                      coords_i,
                                    vertices_i,
                                       point_j,
                                      shapef_j); 

#else
    ptr_class->__simple_interpolation__(  coords_i, 
                                        vertices_i, 
                                           point_j, 
                                          shapef_j);
#endif
  } 


  double* commdom_locator_get_pointer2()
  {
/*  
    CommDom CD = ptr_class[0];
    int size = CD.get_mpi_commij_size();
    std::cout<<" [commdom_locator_get_pointer2] size:" << size <<" \n";

    return ptr_class;
*/
std::vector<double> bar(2,3.14159); 
return &bar[0]; 

  }

  void commdom_locator_get_pointer4_(double*  ptr)
  {
 
 std::cout<<"[commdom_locator_get_pointer4] f:"<< ptr[0] <<" ->";

std::vector<double> bar(2,1.0);

ptr[0] = bar[0];  
 
 std::cout<<" "<< ptr[0] <<" \n";
  }


  void commdom_locator_get_pointer1_(void*  ptr) // recv :) . send :( 
  {

// std::cout<<"[commdom_locator_get_pointer1] f:"<< ptr[0] <<" ->"; 

std::vector<double> bar(1,2.71828);

ptr = &bar[0];   

// std::cout<<" "<< ptr[0] <<" \n";

//double bar = 3;
//ptr = &bar;  

/* 
    CommDom* CD = (CommDom*) commdom_locator_get_pointer2();

    int size = CD[0].get_mpi_commij_size();
    std::cout<<" [commdom_locator_get_pointer1] size:" << size <<" \n";
*/
  }
 

  void commdom_locator_get_pointer3_(void* ptr)
  {
/* 
    CommDom* CD = (CommDom*) ptr; //commdom_locator_get_pointer2();
    
    int size = CD[0].get_mpi_commij_size();
    std::cout<<" [commdom_locator_get_pointer3] size:" << size <<" \n";
*/

double received = *( (double*)ptr ); 

std::cout<<" [commdom_locator_get_pointer3] received:" << received <<" \n";
  }




}

/* 
//=======================================================================||===//
extern "C"
{
// static CommDom* ptr_class = NULL;

  void plepp_locator_get_pointer_(void* ptr)
  {
    ptr = ptr_class; 
  }


  void plepp_locator_destroy_(void* ptr)
  {
    assert(ptr);

    CommDom* ptr_plepp = (CommDom*) ptr;     
    ptr_plepp->locator_destroy(); 

  }


  void plepp_locator_create2_(void* ptr, int* fcommi, int* fcommij, double* tol, int* element_type)
  { 
    MPI_Comm commi  = MPI_Comm_f2c( fcommi [0]);
    MPI_Comm commij = MPI_Comm_f2c( fcommij[0]);

    assert(ptr);
    CommDom* ptr_plepp = (CommDom*) ptr;
    ptr_plepp->locator_create2(commi, commij, tol[0], element_type[0]);

  }




}  
*/

//=======================================================================||===//
//=======================================================================||===//
