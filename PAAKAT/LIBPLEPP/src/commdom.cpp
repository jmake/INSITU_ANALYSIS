#include "commdom.hpp"
#include <algorithm>  // std::search 
#include <iomanip>    // setw 
#include <cmath>      // fabs
#include <fstream>
#include <sstream>     
#include <numeric>
#ifdef SATURNE
  #include "cs_alya.h"
#include "cs_alya.c"
#endif 
//-----------------------------------------------------------------------||---//
// NOTE: syr only tetras !!
/*
ple_mesh_extents_t           *mesh_extents_f   = syr_cfd_point_location_extents;
ple_mesh_elements_contain_t  *locate_inside_f  = syr_cfd_point_location_contain;
ple_mesh_elements_closest_t  *locate_closest_f = syr_cfd_point_location_closest;
*/
ple_mesh_extents_t           *mesh_extents_f   = NULL;
ple_mesh_elements_contain_t  *locate_inside_f  = NULL;
ple_mesh_elements_closest_t  *locate_closest_f = NULL;
//-----------------------------------------------------------------------||---//


#ifndef __PRINT__
#define __PRINT__ 0 
#endif 
//=======================================================================||===//
//=======================================================================||===//


//=======================================================================||===//
//================================================================| COMMs |===//

//-----------------------------------------------------------------------||---//
CommDom::CommDom()
{
  if(__PRINT__) cout<<"+CommDom Create \n"; 
}


CommDom::~CommDom()
{
  if(__PRINT__) cout<<"+CommDom \'"<< IAM <<"\' End \n"; 
}


void 
CommDom::init()
{  
       IAM = "";
  app_name = "";

  quitze_comm = MPI_COMM_NULL;
 //local_comm = MPI_COMM_NULL;
   world_comm = MPI_COMM_NULL;
       commij = MPI_COMM_NULL;

   local_size = -1; 
   local_rank = -1;
   world_size = -1; 
   world_rank = -1;

  app_namei = ""; 
  app_namej = ""; 

  cplng    = NULL; 

  n_meshes =  0; 
  n_locs   =  0; 
  n_apps   = -1; 
  n_types  = -1; 
  app_id   = -1; 

//  _cs_coupling_sync_flag = PLE_COUPLING_NO_SYNC;
//  _cs_coupling_sync_flag = 0; /* Synchronized */

  dimension = -1;

  if(__PRINT__) cout<<"|_CommDom.init() \n"; 
}


void 
CommDom::end()
{  
  if(__PRINT__) cout<<"|_CommDom.end() \n"; 
}
//-----------------------------------------------------------------------||---//

//-----------------------------------------------------------------------||---//
void
CommDom::set_app_type(string last_name, int length) 
{
  //last_name.erase(std::remove_if(last_name.begin(), last_name.end(), ::isspace), last_name.end());
  if(length==-1) 
  {
    app_type = last_name; 
  }
  else
  {
  /*
    //  cout<< last_name.size() <<" \n\n";
    cout<<"(1)->\'"<< last_name <<"\' \n\n";
    app_type = string_f2c(last_name.c_str(), last_name.size());
    cout<<"(2)->\'"<< app_type <<"\' \n\n";
  */
  }
  
  if(__PRINT__) cout<<"|_CommDom.set_app_type \n"; 
}


void
CommDom::set_app_name(string first_name)
{
  app_name = first_name; 
  IAM      = first_name; 

  if(__PRINT__) cout<<"|_CommDom.set_app_name \n"; 
}


void 
CommDom::__genesis__() 
{
//
//  |_ ./src/syrthes-kernel/src/syrthes.c
//    |_ syrthes_initmpi
//
//  if exit "--name"
//  cfd.do_coupling = 1; // mainsyrthes.c 
//

  MPI_Comm  Gucumatz = MPI_COMM_WORLD; 

  int w_rank, w_size;
  MPI_Comm_rank(Gucumatz, &w_rank);
  MPI_Comm_size(Gucumatz, &w_size);

  int id_type = ple_coupling_mpi_name_to_id(Gucumatz, app_type.c_str());
/*
int id_name = ple_coupling_mpi_name_to_id(Gucumatz, app_name.c_str());
*/
  quitze_comm = MPI_COMM_NULL;
  if(id_type > -1) 
  {

    MPI_Comm_split(Gucumatz, id_type, w_rank, &quitze_comm);

    int sync_flag = _cs_coupling_sync_flag;
    ple_coupling_mpi_set_t *cplng_aux = 
                           ple_coupling_mpi_set_create(sync_flag, app_type.c_str(), app_name.c_str(), Gucumatz, quitze_comm);
    
    
    n_types = ple_coupling_mpi_set_n_apps(cplng_aux);

    cout<<"|_["; 
    //cout<<"n_types:"<< n_types; 
    cout<< w_rank+1 <<"."<< app_type; 
    //cout<<"."<< app_name;  
    cout<< "] ";
    //cout<< id_type+1 <<"/"<< id_name+1 <<" ";  
    __print_comm__(quitze_comm);
    cout<<"\n";
  }
  else
  {
    n_types = 1; 
  }

  //cout<< "<--------";
  //cout<<"\n";
}


void 
CommDom::set_world_comm(MPI_Comm  comm)
{
  world_comm = comm; 
  if(quitze_comm != MPI_COMM_NULL) world_comm = quitze_comm; // communicator given by app_type  

  MPI_Comm_size(world_comm, &world_size);
  MPI_Comm_rank(world_comm, &world_rank); 

  if(__PRINT__) cout<<"|_CommDom.set_world_comm \n"; 
}


void 
CommDom::name_to_id(int* appi, int* napps, MPI_Comm* lcomm)
{
  MPI_Comm local_comm = world_comm;
 
  app_id = ple_coupling_mpi_name_to_id(world_comm, app_name.c_str());
  if(app_id > -1) MPI_Comm_split(world_comm, app_id, world_rank, &local_comm); 

  int gsize = -1;
  MPI_Comm_rank(world_comm, &gsize);

  if(gsize <= 4) 
  {
    cout<<"|_["<< world_rank+1 <<"."<< app_type <<"."<< app_name << "] ";
    __print_comm__(local_comm); cout<<"\n";
  } 

  __set_create__(app_name, local_comm, &n_apps); 
  if(app_id < 0){ app_id = 0; n_apps= 1;}

   appi[0] = app_id; 
  napps[0] = n_apps; 
  lcomm[0] = local_comm; 

  __set_get_info__(); 

  if(__PRINT__) cout<<"|_CommDom.name_to_id \n"; 
}


void
CommDom::__set_create__(string namei, MPI_Comm commi, int* napps) 
{
  int sync_flag = _cs_coupling_sync_flag; // & PLE_COUPLING_UNSTEADY 

  cplng = ple_coupling_mpi_set_create(sync_flag, app_type.c_str(), namei.c_str(), world_comm, commi);
  napps[0] = ple_coupling_mpi_set_n_apps(cplng);
}


void
CommDom::__set_get_info__()
{
  for(int i = 0; i < n_apps; i++) // n_apps <-- name_to_id() 
  { 
    ple_coupling_mpi_set_info_t  infoi = ple_coupling_mpi_set_get_info(cplng, i);
    Info.push_back( infoi ); 

    WhoIam[infoi.app_name] = i; 
    //CplngIJ.push_back( vector<int>() ); 
  }
}


void
CommDom::__create_interaction__()
{
//  cout<<"|_["<< world_rank+1 <<".\'"<< app_type <<"\'.\'"<< app_name << "\'] "; 

  string namei = app_name; 

//  cout<<" | "; 

  for(int i = 0; i < n_apps; i++) CplngIJ.push_back( vector<int>() ); 
  
  for(Whoj = WhoIam.begin(); Whoj!=WhoIam.end(); ++Whoj) 
  {
    string  namej = Whoj->first; 
    if( (namei==namej) )
    { 
    }
    else
    {
      bool you_are_friendly = true; //find( friendsi.begin(), friendsi.end(), namej) != friendsi.end(); 
      if(you_are_friendly==true)
      {
//        cout<<"["<< namej <<"]";
        CplngIJ[ WhoIam[namei] ].push_back( WhoIam[namej] ); 
      }
      else
      {
      }
    }
  }
//  cout<<" | "; 


  int appi = WhoIam[app_name]; 
//  cout<<"\n";

  if(__PRINT__) cout<<"|_CommDom.__create_interaction__ \n"; 
}


void
CommDom::set_applications(string last_name)
{
}


void
CommDom::set_iteration(string namej)
{
}

void
CommDom::create_interaction()
{
}
//-----------------------------------------------------------------------||---//


//=======================================================================||===//
//================================================================| COMMs |===//
//-----------------------------------------------------------------------||---//
void
CommDom::create_commij(MPI_Comm lcomm) 
{
  int appi = -1, appj = -1; 
  
  appi = WhoIam[app_name];  // <-- Sooo... important!!
  //if(CplngIJ[appi].size()<1)  __error__(true, "interaction doesnt exist");

  string namei = app_name; //cout<<"|_["<< namei;
  for(int j = 0; j < CplngIJ[appi].size(); j++)
  {
    appj = CplngIJ[appi][j]; 

    string namej = Info[appj].app_name; 
//    if(world_size<=4) cout<<"|_["<< namei <<"."<< Info[appi].root_rank+1 <<"]->["<< namej <<"."<< Info[appj].root_rank+1 <<"]: "; 

    vector<int> rangei(2); 
    rangei[0] = Info[appi].root_rank; 
    rangei[1] = Info[appi].root_rank + Info[appi].n_ranks; 
    
    vector<int> rangej(2); 
    rangej[0] = Info[appj].root_rank; 
    rangej[1] = Info[appj].root_rank + Info[appj].n_ranks; 
/*
    if(world_size<=4) 
    {
      cout<<"|";
      cout<<"["<< appi+1 <<"]("<< rangei[0]+1 <<","<< rangei[1] <<")->"; 
      cout<<"["<< appj+1 <<"]("<< rangej[0]+1 <<","<< rangej[1] <<")| "; 
      //cout<<"|\n";
    }
*/ 
    int commij_size = -1; 
    int commij_rank = -1;

    commij = MPI_COMM_NULL;
    ple_coupling_mpi_intracomm_create(world_comm, lcomm, rangej[0], &commij, rangei.data(), rangej.data());  
    CommIJ[namej] = commij; 

    MPI_Comm_size(commij, &commij_size);
    MPI_Comm_rank(commij, &commij_rank); 
/* 
    if(world_size<=4)
    {
      cout<< commij_rank+1 <<"/"<< world_rank+1 <<","<< commij_size <<"/"<< world_size; 
      cout<<"\n";
    }
*/
    int rooti = -1; 
    int rootj = -1; 
    __mpi_get_roots__(lcomm, commij, &rooti, &rootj); 
    RootIJ[namej].push_back(rooti);
    RootIJ[namej].push_back(rootj);

  }
  MPI_Barrier(lcomm); 

  if(__PRINT__)
  {
    cout<<"|_CommDom.create_commij: |"; 
    //__print_comm__(lcomm); 
    //__print_comm__(world_comm); 
    print_commij_names(); 
    cout<<"|\n";
  }
}
//-----------------------------------------------------------------------||---//


//-----------------------------------------------------------------------||---//
void
CommDom::get_commij_size(int* size) 
{
  size[0] = CommIJ.size(); 
}


void
CommDom::print_commij_names() 
{
  cout<<"["<< app_type <<"."<< app_name <<"."<< world_rank+1 <<"] -> "; 
  for(CommIT = CommIJ.begin(); CommIT != CommIJ.end(); ++CommIT) 
  cout<<"["<< CommIT->first <<"."<< app_type <<"."<< MPI_Comm_c2f(CommIT->second) <<"]"; 
}


void
CommDom::get_commij(string  namej, MPI_Comm*  commij) 
{
  commij[0] = MPI_COMM_NULL; 

  CommIT = CommIJ.find(namej); 
  if(CommIT == CommIJ.end()) 
  {
    cout<<"\n\nERROR: ["; 
    cout<< app_type <<".";  
    cout<< namej;
    cout<<"] not found. INTERNSHIPS->"; 
    print_commij_names(); 
    cout<<"!!\n\n";
    exit(1);
  }
  else
  {
//    if(local_size<=4) cout<<"|_["<< IAM <<"."<< world_rank+1 <<"]<->["<< namej <<"] "; 
    commij[0] = CommIJ[namej]; 
  }
}


MPI_Comm
CommDom::set_mpi_comms() // used by python 
{
  MPI_Comm lcomm = MPI_COMM_NULL;
  int appi  =-1; 
  int napps =-1; 
  
  name_to_id(&appi, &napps, &lcomm); 
  __create_interaction__(); 
  create_commij(lcomm);

  return lcomm; 
}


int
CommDom::get_mpi_commij_size() // used by python 
{
  int size;
  get_commij_size(&size);
  return size; 
}


MPI_Comm
CommDom::get_mpi_commij(string namej) // used by python 
{
  MPI_Comm  commij = MPI_COMM_NULL; 
  get_commij(namej, &commij);
  return commij; 
}
//-----------------------------------------------------------------------||---//


//=======================================================================||===//
//===================================================================| MPI|===//
//-----------------------------------------------------------------------||---//
int
CommDom::__mpi_get_roots__(MPI_Comm commi, MPI_Comm commij, 
                           int* rooti, int* rootj)
{
  int ranki  = -1;
  int rankij = -1;
  int sizei  = -1;
  int sizeij = -1;

  MPI_Comm_size(commi,  &sizei );
  MPI_Comm_size(commij, &sizeij);
  MPI_Comm_rank(commi,  &ranki );
  MPI_Comm_rank(commij, &rankij);

  int rootij = -666; 
  if(ranki==rankij)
  {
    rooti[0] = 0;
    rootj[0] = sizei;
    rootij = 1;
  }
  else
  {
    rooti[0] = sizeij-sizei;
    rootj[0] = 0;
  }

  //if(rankij==0) rootij = 1; 

  return rootij; 
}
//-----------------------------------------------------------------------||---//

//-----------------------------------------------------------------------||---//
void
CommDom::sendrecv_int(int* send,  int n_send,
                      int* recv,  int n_recv,
                      MPI_Comm  commi,
                      MPI_Comm  commij)

{
  __mpi_sendrecv_int__(send,  n_send,
                       recv,  n_recv,
                       commi, commij);
}


void
CommDom::__mpi_sendrecv_int__(int* varij, int n_varij, 
                              int* varji, int n_varji, 
                              MPI_Comm commi, 
                              MPI_Comm commij)
{
  int rooti = -1; 
  int rootj = -1; 
  __mpi_get_roots__(commi, commij, &rooti, &rootj); 

  int ranki  = -1;
  MPI_Comm_rank(commi, &ranki); 

  int rut = -1; 
  if(ranki==0)
  { 
    rut = rootj; 
    MPI_Status  status; 
    const int  syr_cfd_coupling_tag = 'C'+'S'+'_'+'C'+'O'+'U'+'P'+'L'+'A'+'G'+'E';

    MPI_Sendrecv(varij, n_varij, MPI_INTEGER, rut, syr_cfd_coupling_tag, 
                 varji, n_varji, MPI_INTEGER, rut, syr_cfd_coupling_tag,
                 commij, &status); 

  }

  //MPI_Barrier(commi);
}


void
CommDom::__mpi_sendrecv_real__(double* varij, int n_varij,
                               double* varji, int n_varji,
                               MPI_Comm commi,
                               MPI_Comm commij)
{
  int rooti = -1; 
  int rootj = -1; 
  __mpi_get_roots__(commi, commij, &rooti, &rootj); 

  int ranki  = -1;
  MPI_Comm_rank(commi, &ranki); 

  int rut = -1;
  if(ranki==0)
  {
    rut = rootj;
    MPI_Status  status;
    const int  syr_cfd_coupling_tag = 'C'+'S'+'_'+'C'+'O'+'U'+'P'+'L'+'A'+'G'+'E';

    MPI_Sendrecv(varij, n_varij, MPI_DOUBLE, rut, syr_cfd_coupling_tag,
                 varji, n_varji, MPI_DOUBLE, rut, syr_cfd_coupling_tag,
                 commij, &status);
  }

}



void
CommDom::__mpi_sendrecv_char__(char* varij, int n_varij,
                               char* varji, int n_varji,
                               MPI_Comm commi,
                               MPI_Comm commij)
{
  int rooti = -1; 
  int rootj = -1; 
  __mpi_get_roots__(commi, commij, &rooti, &rootj); 

  int ranki  = -1;
  MPI_Comm_rank(commi, &ranki); 

  int rut = -1;
  if(ranki==0)
  {
    rut = rootj;
    MPI_Status  status;
    const int  syr_cfd_coupling_tag = 'C'+'S'+'_'+'C'+'O'+'U'+'P'+'L'+'A'+'G'+'E';

    MPI_Sendrecv( varij, n_varij, MPI_CHAR, rut, syr_cfd_coupling_tag,
                  varji, n_varji, MPI_CHAR, rut, syr_cfd_coupling_tag,
                  commij, &status);
  }

}


void
CommDom::__mpi_bcast_real__(double*     input,
                            int       n_input,
                            MPI_Comm  commi,
                            MPI_Comm  commij)
{
  int rooti = -1;
  int rootj = -1;
  __mpi_get_roots__(commi, commij, &rooti, &rootj);

  int ranki  = -1;
  MPI_Comm_rank(commi, &ranki);

  double   varij = -666;
  int    n_varij = -666;
  if(ranki == 0)
  {
      varij =   input[0];
    n_varij = n_input;
  }

  MPI_Datatype  datatype = MPI_DOUBLE;

  MPI_Bcast(&n_varij,       1, MPI_INTEGER, 0, commi);
  MPI_Bcast(  &varij, n_varij,    datatype, 0, commi);

  input[0] = varij;

}
//-----------------------------------------------------------------------||---//

//=======================================================================||===//
//==============================================================| LOCATOR |===//

//-----------------------------------------------------------------------||---//
void
CommDom::locator_create(
                        MPI_Comm   commij,
                        int        n_rank,
                        int     root_rank, 
                        double        tol
                        )
{
  __error__(true, "USE locator_create2(commi, commij); !!"); 
}


void
CommDom::locator_create2(
                         MPI_Comm   commi,
                         MPI_Comm   commij, 
                         double     tol,
                         int        etype
                        )
{
  n_locs += 1;
  if(n_locs>1) __error__(true, "n_locs>1. Only ONE locator can be create. Multi-locators must be implemented !!"); 

  int ranki  = -1;
  int rankij = -1;
  int sizei  = -1;
  int sizeij = -1;

  MPI_Comm_size(commi,  &sizei );
  MPI_Comm_size(commij, &sizeij);
  MPI_Comm_rank(commi,  &ranki );
  MPI_Comm_rank(commij, &rankij);

  int rooti=-1, rootj=-1; 
  if(ranki==rankij)
  {
    rooti = 0;
    rootj = sizei;
  }
  else
  {
    rooti = sizeij-sizei;
    rootj = 0;
  }


  int root_rank = rootj;        //= Info[appj].n_ranks;
  int n_rank    = sizeij-sizei; //= Info[appj].root_rank;

  locatorj = ple_locator_create(tol, commij, n_rank, root_rank);

  __locator_end_meshes__();
  __locator_init_meshes__(); 

  assert( etype==4 || etype==6 );  // 2019FEB21  
  element_type = etype;            // 2019FEB21  
}


ple_locator_t*
CommDom::__locator_get__()
{
  return locatorj;
}


void
CommDom::locator_destroy()
{
  //cout<<"locator_destroy\n"; 

  n_locs = 0; 
  ple_locator_destroy(locatorj); 

  n_meshes = 0; 
  __locator_end_meshes__(); 
  //__locator_init_meshes__(); 
}


void
CommDom::__locator_init_meshes__() 
{
  LocalMesh.push_back( syr_cfd_mesh_t() ); 
  n_meshes += 1; 

  //if(n_meshes>1) __error__(false, "n_meshes>1. Only ONE locator can be create. Multi-meshes must be implemented !!"); 
  if(n_meshes>1) cout<<"\n  n_meshes="<< n_meshes << ". Only ONE locator can be create. Multi-meshes must be implemented !!\n\n"; 

  // MESHi 
  LocalMesh[0].dim           = dimension; // __DIM__; 
//LocalMesh[0].element_type  = SYR_CFD_TETRA or SYR_CFD_TRIA;
  LocalMesh[0].n_vertices    = 0; 
  LocalMesh[0].n_elements    = 0; 
  LocalMesh[0].vertex_coords = NULL;
  LocalMesh[0].vertex_num    = NULL; 

  // MESHj
  n_coords_dist = 0; 
  coords_dist   = NULL; 
  distance_dist = NULL;
}


void
CommDom::__locator_end_meshes__() 
{
  LocalMesh.clear(); 
}
//-----------------------------------------------------------------------||---//

//-----------------------------------------------------------------------||---//
void
CommDom::locator_set_mesh(
                          int         n_vertices_i, 
                          int         n_elements_i, 
                          REAL*       vertex_coords_i, //double*     vertex_coords_i, 
                          INTEGER*    vertex_num_i,    //int*        vertex_num_i, 
                          int         n_vertices_j, 
                          REAL*       vertex_coords_j  //double*     vertex_coords_j
                          ) 
{
  if(dimension<=0)
  {
    cout<<"dimension:"<< dimension <<"\n";
    __error__(true, " ERROR: 'locator_set_mesh' dimension<=0 !!");
  }

  __locator_set_meshi__(n_vertices_i, n_elements_i, vertex_coords_i, vertex_num_i);
  __locator_set_meshj__(n_vertices_j, vertex_coords_j);
  __locator_set_mesh__( locatorj ); 
}


void 
CommDom::locator_set_cs_mesh(
                             int         n_vertices_i,    // <-
                             int         n_elements_i,    // <-
                             REAL*       vertex_coords_i, // <-
                             INTEGER*    vertex_num_i,    // <- 
                             INTEGER*    vertex_type_i,   // <-
                             int         n_vertices_j,    // <- 
                             REAL*       vertex_coords_j, // <-
                             int         dim,
                             REAL*       vertex_props_j,
                             int         dim_props 
                            )
{
  dimension = dim;
  dimension_props = dim_props;

#if SATURNE
  __locator_set_meshj__(n_vertices_j, vertex_coords_j);  

  Alya2Cs[10] = 1; //TRI03 = 10, FVM_FACE_TRIA
  Alya2Cs[12] = 2; //QUA04 = 12, FVM_FACE_QUAD
  Alya2Cs[30] = 4; //TET04 = 30, FVM_CELL_TETRA 
  Alya2Cs[32] = 5; //PYR05 = 32, FVM_CELL_PYRAM
  Alya2Cs[34] = 6; //PEN06 = 34, FVM_CELL_PRISM
  Alya2Cs[37] = 7; //HEX08 = 37, FVM_CELL_HEXA

  INTEGER *cs_type_i = new INTEGER[n_elements_i];
  __alya2cs_type__(vertex_type_i, cs_type_i, n_elements_i); 

  __create_nodal_t__(vertex_coords_i, // <- 
                     n_vertices_i,    // <-
                     vertex_num_i,    // <- 
                     cs_type_i,       //  
                     n_elements_i);   // <- 

  delete[] cs_type_i;  
#else 
   props_dist = vertex_props_j;  
   locator_set_mesh(   n_vertices_i, 
                       n_elements_i, 
                       vertex_coords_i, 
                       vertex_num_i,  
                       n_vertices_j, 
                       vertex_coords_j); 
#endif 
}


void
CommDom::__locator_set_meshi__(int         n_vertices_i, 
                               int         n_elements_i, 
                               REAL*       vertex_coords_i, 
                               INTEGER*    vertex_num_i 
                              ) 
{
  LocalMesh[0].dim           = dimension;       //  __DIM__; 

  if(element_type==4)
  {
    if(dimension==3) LocalMesh[0].element_type  = SYR_CFD_TETRA;
    if(dimension==2) LocalMesh[0].element_type  = SYR_CFD_TRIA;
  }
  if(element_type==6)
  {
    if(dimension==3) LocalMesh[0].element_type = FVM_CELL_HEXA;
    if(dimension==2) assert(0);
  }

  LocalMesh[0].n_vertices    = n_vertices_i;
  LocalMesh[0].n_elements    = n_elements_i;

  if(n_vertices_i>=0) LocalMesh[0].vertex_coords = vertex_coords_i;
  if(n_elements_i>=0) LocalMesh[0].vertex_num    = vertex_num_i;
}


void
CommDom::__locator_set_meshj__(int       n_vertices_j, 
                               REAL*  vertex_coords_j) 
{

  n_coords_dist = 0; 
  coords_dist   = NULL; 
  distance_dist = NULL;

  if(n_vertices_j > 0) 
  {
    n_coords_dist = n_vertices_j;
      coords_dist = vertex_coords_j;
    distance_dist = new float[n_coords_dist];
  }

}


void
CommDom::__locator_set_mesh__(ple_locator_t  *locator) 
{
  int               dim = dimension; // __DIM__; 
  int               dim_props = dimension_props;
  const int *point_list = NULL; 
  
  __locator_set_meshij__(locator, 
                        &LocalMesh[0], 
                         dim,
                         dim_props, 
                         n_coords_dist, 
                         point_list, 
                         coords_dist, 
                         distance_dist);
  delete[] distance_dist; 

}


void
CommDom::__locator_set_meshij__(ple_locator_t      *locator, 
                                syr_cfd_mesh_t     *Mesh, 
                                int                 dim,
                                int                 dim_props, 
                                int                 n_coords_dist, 
                                const int          *point_list, 
                                double             *coords_dist, 
                                float              *distance_dist 
                              ) 
{
      mesh_extents_f   = syr_cfd_point_location_extents;
      locate_inside_f  = syr_cfd_point_location_contain;
    //locate_closest_f = syr_cfd_point_location_closest;
      locate_closest_f = NULL;               // Why dosent it work? R. 3d not implemented 

      ple_locator_set_mesh(locator,          //*ple_locator_t
                           Mesh,             // void*
                           dim,              // coord_dim!!
                           dim_props, 
                           n_coords_dist,    // ple_lnum_t 
                           point_list,       // const ple_lnum_t[] 
                           coords_dist,      // const ple_coord_t[]
                           distance_dist,    // float[]
                           mesh_extents_f,   // ple_mesh_extents_t.          cs_coupling_mesh_extents
                           locate_inside_f,  // ple_mesh_elements_contain_t. cs_coupling_point_in_mesh. 
                           locate_closest_f, props_dist);// ple_mesh_elements_closest_t. 
                                             // locate_on_closest = cs_coupling_point_closest_mesh. optional 

    //ple_locator_dump(locator); 

    //syr_cfd_coupling._interpolation_init ?? 
    ple_locator_exchange_point_var(locator,
                                   NULL,
                                   distance_dist,
                                   NULL,
                                   sizeof(float),
                                   1,
                                   1);


    //int appj = WhoIam[app_name];  // <-- Sooo... important!!
    //int rut = Info[appj].root_rank;
    if(world_size<=8)
    {
        cout<<"|_["<< app_type <<"."<< app_name <<"] ";
        cout<<"n_elemets: "; 
        cout<< Mesh[0].n_elements <<", "; 
        cout<<"n_exterior: ";
        cout<< ple_locator_get_n_exterior(locator)<<", ";
        cout<<"n_interior: ";
        cout<< ple_locator_get_n_interior(locator)<<", ";
        cout<<"n_dist_points: "; 
        cout<< ple_locator_get_n_dist_points(locator); 
        cout<<"\n";
    }

    __set_synchronization_flags();
}
//-----------------------------------------------------------------------||---//

//-----------------------------------------------------------------------||---//
void
CommDom::__locator_get_interior_list__(int* ids) 
{
  int          n_interior   = ple_locator_get_n_interior(locatorj);
  const int  *interior_list = ple_locator_get_interior_list(locatorj); 

  for(int i=0; i<n_interior; i++) ids[i] = interior_list[i]; 
}


void
CommDom::__locator_get_dist_locations__(int* ids)      // NOTE-> to_send: list cell_i/point_j 
{
  int        n_dist_points   = ple_locator_get_n_dist_points(locatorj);
  const int *dist_locations  = ple_locator_get_dist_locations(locatorj);

  for(int i=0; i<n_dist_points; i++) ids[i] = dist_locations[i];
}


void
CommDom::__locator_get_dist_coords__(double* coords)   // NOTE-> to_send: list  point_j 
{
  int           n_dist_points = ple_locator_get_n_dist_points(locatorj);
  const double   *dist_coord  = ple_locator_get_dist_coords(locatorj);

//for(int i=0; i<n_dist_points*__DIM__; i++) coords[i] = dist_coord[i];
  for(int i=0; i<n_dist_points*dimension; i++) coords[i] = dist_coord[i];
}


void //MATIAS
CommDom::__locator_get_dist_props__(double* coords)   // NOTE-> to_send: list  point_j 
{
  int           n_dist_points = ple_locator_get_n_dist_points(locatorj);
  const double   *dist_coord  = ple_locator_get_dist_props(locatorj);

//for(int i=0; i<n_dist_points*__DIM__; i++) coords[i] = dist_coord[i];
  for(int i=0; i<n_dist_points*dimension_props; i++) coords[i] = dist_coord[i];
}


int
CommDom::get_n_dist_points() 
{
  return ple_locator_get_n_dist_points(locatorj); 
}


int
CommDom::get_n_interior() 
{
  return ple_locator_get_n_interior(locatorj); 
}
//-----------------------------------------------------------------------||---//


//-----------------------------------------------------------------------||---//
void 
CommDom::__locator_exchange_double_scalar__(double *var_ij, double *var_ji, int stride) 
{
  //var_ij = new double[stride*n_dist_points];
  //var_ji = new double[stride*n_interior];
  //int stride  = 1;
  //cout<<"stride:"<< stride <<"\n";

  if(stride<1) __error__(true, "exchange_double_scalar. stride<1!!"); 

  int reverse   = 0;
  int type_size = sizeof(double);

  ple_locator_exchange_point_var(locatorj,
                                 var_ij, // to be sent
                                 var_ji, // to be received
                                 NULL,   // ?? 
                                 type_size, stride, reverse);
}


void 
CommDom::__locator_send_double_scalar__(double *var_ij, int stride) 
{
/*
  if(stride<1) __error__(true, "exchange_double_scalar. stride<1!!"); 

  int       reverse = 0;
  int     type_size = sizeof(double);
  double    *var_ji = NULL;

  ple_locator_exchange_point_var(locatorj,
                                 var_ij, // to be sent
                                 var_ji, // to be received
                                 NULL,   // ?? 
                                 type_size, stride, reverse);
*/

  double *var_ji = NULL;
  __locator_exchange_double_scalar__(var_ij, var_ji, stride); 
}


void 
CommDom::__locator_recv_double_scalar__(double *var_ji, int stride) 
{
/*
  if(stride<1) __error__(true, "exchange_double_scalar. stride<1!!"); 

  int      reverse = 0;
  int    type_size = sizeof(double);
  double   *var_ij = NULL;

  ple_locator_exchange_point_var(locatorj,
                                 var_ij, // to be sent
                                 var_ji, // to be received
                                 NULL,   // ?? 
                                 type_size, stride, reverse);
*/

  double *var_ij = NULL;
  __locator_exchange_double_scalar__(var_ij, var_ji, stride); 
}
//-----------------------------------------------------------------------||---//


//=======================================================================||===//
//==============================================================| SYRTHES |===//
void
CommDom::__syr_cfd_coupling_face_vars__() 
{
/*
  syr_cfd_coupling.syr_cfd_coupling_face_vars(syr_cfd_coupling_t  *coupling,
                           const double        *t_node,
                           double              *t_face,
                           double              *h_face)

  PLE_MALLOC(recv_buf, ip->n_coupled_elts*2, double);

  ple_locator_exchange_point_var(ip->locator,
                                 NULL,
                                 recv_buf,
                                 NULL,
                                 sizeof(double),
                                 2,
                                 0);

*/

      // exchange 
      int type_size = sizeof(double);
      int stride    = 2;
      int reverse   = 0;

      double *var_ij = NULL; // to be sent
      double *var_ji = NULL; // to be received

      int n_interior = ple_locator_get_n_interior(locatorj);
      var_ji = new double[stride*n_interior];
      for(int i=0; i<stride*n_interior;    i++) var_ji[i] = 269.0; 

      ple_locator_exchange_point_var(locatorj,
                                     var_ij,
                                     var_ji,
                                     NULL, // ?? 
                                     type_size, stride, reverse);

}


void
CommDom::_send_nodal_var(int stride)
{
/*
  ple_lnum_t n_dist = ple_locator_get_n_dist_points(ip->locator);
  PLE_MALLOC(var_send, n_dist, double);

  ple_locator_exchange_point_var(ip->locator,
                                 var_send,
                                 NULL,
                                 NULL,
                                 sizeof(double),
                                 1,
                                 0);
*/

      // exchange 
      double *var_ij = NULL; // to be sent
      double *var_ji = NULL; // to be received

/*
      int n_dist_points = ple_locator_get_n_dist_points(locatorj);
      var_ij = new double[stride*n_dist_points];
      for(int i=0; i<stride*n_dist_points; i++) var_ij[i] = 0.0; 
*/

      int n_interior = ple_locator_get_n_interior(locatorj);
      var_ji = new double[stride*n_interior];
      for(int i=0; i<stride*n_interior; i++) var_ji[i] = 0.0;

/*
      int type_size = sizeof(double);
      int reverse   = 0;

      ple_locator_exchange_point_var(locatorj,
                                     var_ij,
                                     var_ji,
                                     NULL, // ?? 
                                     type_size, stride, reverse);
*/


  __locator_exchange_double_scalar__(var_ij, var_ji, stride); 

  double tf_avrg = 0.0;
  double hf_avrg = 0.0;
  for (int i = 0; i < n_interior; i++) 
  {
    tf_avrg += var_ji[i*2];
    hf_avrg += var_ji[i*2 + 1];
  }

tf_avrg = (n_interior==0)?(0.0):(tf_avrg/n_interior);
hf_avrg = (n_interior==0)?(0.0):(hf_avrg/n_interior); 

printf("\n"); 
printf("-----> "); 
printf("tf_avrg: %f, ", tf_avrg); 
printf("hf_avrg: %f  ", hf_avrg); 
printf("<-----"); 
printf("\n"); 

}


void
CommDom::_send_nodal_var_vec()
{
   _send_nodal_var(2); 
}


void
CommDom::_send_elt_var()
{
      // exchange 
      int type_size = sizeof(double);
      int stride    = 1;
      int reverse   = 0;

      double *var_ij = NULL; // to be sent
      double *var_ji = NULL; // to be received

      int n_dist_points = ple_locator_get_n_dist_points(locatorj);
      var_ij = new double[stride*n_dist_points];
      for(int i=0; i<stride*n_dist_points; i++) var_ij[i] = 0.0; 

      ple_locator_exchange_point_var(locatorj,
                                     var_ij,
                                     var_ji,
                                     NULL, // ?? 
                                     type_size, stride, reverse);
}

/*
void
CommDom::__syr_cfd_coupling_ensure_conservativity__()
{

}
*/

//=======================================================================||===//
//==================================================================| AUXs|===//

//-----------------------------------------------------------------------||---//
void
CommDom::__print_comm__(MPI_Comm  comm)
{
  int size=-1, rank=-1;
  MPI_Comm_size(comm, &size);
  MPI_Comm_rank(comm, &rank); 

  MPI_Fint  fcomm = MPI_Comm_c2f(comm); 

  cout<<"["<< rank+1 <<"/"<< size <<"]"; 
  //cout<<"["<< rank+1 <<"/"<< size <<"."<< fcomm <<"]"; 
}
//-----------------------------------------------------------------------||---//


//-----------------------------------------------------------------------||---//
void
CommDom::__error__(bool ok, string message)
{
  if(ok==true) 
  {
    cout<<"\n\nERROR: "<< message  <<"!! \n\n"; 
    exit(1);
  }
}


void
CommDom::__strncmp__(const char* x, const char* y, int* size, int* ok)
{
  ok[0] = -1;
  if(strncmp(x, y, size[0]) == 0) 
  {
    ok[0] = 1;
  }
  else
  {
   //printf(" \'%s\', \'%s\' \n", x, y); 
   /*
    cout<<"\n\nWARNING: ";
    cout<<"\'"<< x <<"\' != \'"<< y <<"\'";   
    cout<<"!! \n\n";
    //exit(1);
   */
  } 
}


const char* 
CommDom::string_f2c(const char* licName, int nameLen)
{
/*
  int i = fchar_length-1; 
  while( ((fchar[i]==' ')||(fchar[i]==0)) && (i>=0) ) i--; 
  i+=1; 
    
  std::string  cstr; 
  for(int j=0; j<i; j++) cstr.push_back(fchar[j]); 

  return cstr.c_str(); 
*/ 


  char* name;
  name = (char*)malloc(sizeof(char)*(nameLen));
  for(int i = 0; i<nameLen; i++) name[i] = licName[i];
  name[nameLen] = '\0';

  return name ; 
}


vector<string> 
CommDom::__split_string__(const string& s, const string& delim, const bool keep_empty) 
{
    vector<string> result;
    if(delim.empty()) 
    {
        result.push_back(s);
        return result;
    }
    string::const_iterator substart = s.begin(), subend;

    while(true) 
    {
        subend = search(substart, s.end(), delim.begin(), delim.end());
        string temp(substart, subend);
        if(keep_empty || !temp.empty()) 
        {
            result.push_back(temp);
        }
        
        if(subend == s.end()) 
        {
            break;
        }
        
        substart = subend + delim.size();
    }
    return result;
}
//-----------------------------------------------------------------------||---//


//-----------------------------------------------------------------------||---//
string
CommDom::__get_app_type__()
{
  return app_type; 
}


string
CommDom::__get_app_name__()
{
  return IAM; 
}


int 
CommDom::__get_n_apps__()
{
  return n_apps; 
}


int 
CommDom::__get_n_types__()
{
  return n_types; 
}


MPI_Comm 
CommDom::__get_world_comm__()
{
  return world_comm; 
}


int 
CommDom::__get_friends__(string  namej) 
{
  int ok[] = {1}; 

  CommIT = CommIJ.find(namej); 
  if(CommIT == CommIJ.end()) 
  {
    cout<<"WARNINIG: ["<< namej <<"] not found!! \n"; 
    ok[0] = -1; 
  }

  return ok[0]; 
}
//-----------------------------------------------------------------------||---//


//-----------------------------------------------------------------------||---//
void
CommDom::__to_vtk__(string filename, int npts, const double* pts) 
{
  std::ofstream  fout( filename.c_str() ); 
  fout<<"# vtk DataFile Version 3.0\n";
  fout<<"chingon!!"<<"\n"; 
  fout<<"ASCII \nDATASET POLYDATA \nPOINTS "<< npts <<" float \n";
  for(int i=0, k=0; i<npts; i++) 
  {
    for(int j=0; j<dimension; j++, k++) fout<< pts[k] <<" ";                 // 2016Mar22  
    if(dimension==2) fout<< 0.0 <<" ";
    fout<<"\n";
  }
  fout.close();   
}
//-----------------------------------------------------------------------||---//


//-----------------------------------------------------------------------||---//
void
CommDom::__simple_interpolation__(double* coords, int* vertices, double* point, double* shapef)  // 2016Mar22   
{
  if(element_type==4)
  {
    if(dimension==3)       __tetra_interpolation__(coords, vertices, point, shapef);
    if(dimension==2) __triangle_2d_interpolation__(coords, vertices, point, shapef);
  }
  if(element_type==6)
  {
    if(dimension==3)       __hexa_interpolation__(coords, vertices, point, shapef);
    if(dimension==2) assert(0); // quad 
  }

}

/*
void
CommDom::__simple_interpolation__(double* coords, int* vertices, double* point, double* shapef)  // 2016Mar22   
{
  if(dimension==3) 
          __tetra_interpolation__(coords, vertices, point, shapef); 
  if(dimension==2)
    __triangle_2d_interpolation__(coords, vertices, point, shapef);   
} 
*/

void
CommDom::__tetra_interpolation__(double* coords, int* vertices, double* point, double* shapef)
{
  double _epsilon_denom = 1.e-28; 
  double isop_0, isop_1, isop_2;
  double t00, t10, t20, t01, t02, t03, t11, t12, t13, t21, t22, t23;
  double v01[3], v02[3], v03[3]; //, shapef[4];
  double vol6;
  int i, j, k;


  ple_coord_t  tetra_coords[4][3];
  for(i = 0; i < 4; i++)
  {
    ple_lnum_t   coord_idx = vertices[i] - 1; //<-NOTE: fortran numeration 
    for(j = 0; j < 3; j++) tetra_coords[i][j] = coords[(coord_idx * 3) + j];
  }

  // vol[i] = (<x_j,y_j,z_j> - <x_0,y_0,z_0>)_i = <x_j-x_0, y_j-y_0, z_j-z_0>_i = R_j - R_0 = R_0j 
  for(i = 0; i < 3; i++)
  {
    v01[i] = tetra_coords[1][i] - tetra_coords[0][i];
    v02[i] = tetra_coords[2][i] - tetra_coords[0][i];
    v03[i] = tetra_coords[3][i] - tetra_coords[0][i];
  }

  // 6*Vol = J = A.BxC = R_01 . R_02 x R_03
  vol6 = fabs(  v01[0] * (v02[1]*v03[2] - v02[2]*v03[1])
              - v02[0] * (v01[1]*v03[2] - v01[2]*v03[1])
              + v03[0] * (v01[1]*v02[2] - v01[2]*v02[1]));
  if(vol6 < _epsilon_denom) return;


  // Tetrahedral coordinates 
  t00  = - tetra_coords[0][0] + point[0];
  t10  = - tetra_coords[0][1] + point[1];
  t20  = - tetra_coords[0][2] + point[2];

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

  shapef[0] = 1.0 - isop_0 - isop_1 - isop_2;
  shapef[1] =       isop_0;
  shapef[2] =                isop_1;
  shapef[3] =                         isop_2;

}

void
CommDom::__triangle_2d_interpolation__(double* vertex_coords, int* element_vertex_num, double* point_coords, double* shapef)
{
#define _DOT_PRODUCT_2D(vect1, vect2) (vect1[0] * vect2[0] + vect1[1] * vect2[1])

  // FROM: '_locate_in_triangle_2d'
  double _epsilon_denom = 1.e-28;
  ple_lnum_t  i, j, k, coord_idx_0, coord_idx_1, coord_idx_2;

  double t[2], u[2], v[2]; //, shapef[3];
  double uu, vv, uv, ut, vt, det;
  double dist, max_dist, isop_0, isop_1;

  // vertex index of the triangle studied 
  coord_idx_0 = element_vertex_num[0] - 1;
  coord_idx_1 = element_vertex_num[1] - 1;
  coord_idx_2 = element_vertex_num[2] - 1;

  // Calculate triangle-constant values for barycentric coordinates 
  for (j = 0; j < 2; j++) 
  {
    u[j] = - vertex_coords[(coord_idx_0*2) + j] + vertex_coords[(coord_idx_1*2) + j];
    v[j] = - vertex_coords[(coord_idx_0*2) + j] + vertex_coords[(coord_idx_2*2) + j];
  }

  uu = _DOT_PRODUCT_2D(u, u);
  vv = _DOT_PRODUCT_2D(v, v);
  uv = _DOT_PRODUCT_2D(u, v);

  det = (uu*vv - uv*uv);

  if (det < _epsilon_denom) return;

/*
  // Loop on points resulting from extent query 

  for (k = 0; k < n_points_in_extents; k++) 
  { 
    i =  points_in_extents[k];

    // Calculation of the barycenter coordinates for the projected node 
    for (j = 0; j < 2; j++) t[j] = - vertex_coords[(coord_idx_0*2) + j] + point_coords[i*2 + j];
*/
    for (j = 0; j < 2; j++) t[j] = - vertex_coords[(coord_idx_0*2) + j] + point_coords[j];

    ut = u[0]*t[0] + u[1]*t[1];
    vt = v[0]*t[0] + v[1]*t[1];

    isop_0 = (ut*vv - vt*uv) / det;
    isop_1 = (uu*vt - uv*ut) / det;

    shapef[0] = 1. - isop_0 - isop_1;
    shapef[1] =      isop_0;
    shapef[2] =               isop_1;
/*
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
 
  } // k  
*/
}


int
CommDom::__hexa_interpolation__(double* vertex_coords, int* element_vertex_num, double* point, double* shapef)
{

  int n_elt_vertices = syr_cfd_mesh_n_vertices_element[FVM_CELL_HEXA];
  double         tol = 1e-8;

  vector<double> __shapef(n_elt_vertices,HUGE_VAL);

  _locate_in_cell_3d(FVM_CELL_HEXA,
                     element_vertex_num,
                     vertex_coords,
                     point,
                     tol,
                     __shapef.data()
                     );

  memcpy(shapef, __shapef.data(), __shapef.size()*sizeof(double));

  return 0; 
}

//-----------------------------------------------------------------------||---//

//----------------------------------------------------------------| argvs |---//
void
CommDom::__set_argvs__(string dummy)
{
  argvs.push_back(dummy); 
}


void
CommDom::__get_argvs__(string& dummy)
{
  dummy = argv_found; 
}


void
CommDom::__analyse_argvs__(string token)
{  
  argv_found = ""; 

  string s;
  
  int argc = argvs.size(); 
  int numarg = 0;
  while(++numarg < argc) 
  {
    s = argvs[numarg];
    if( strncmp(s.c_str(), token.c_str(), token.size()) == 0 ) 
    {
      //cout<<"["<< token <<"] = "; 
      
      if(++numarg < argc) 
      {
        s   = argvs[numarg];
        int len = s.size();
        char* log_name = (char *)malloc(sizeof(char)*(len+1));
        strncpy(log_name, s.c_str(), len);
        log_name[len] = '\0';

        //cout<<"\'"<< log_name <<"\'"; 

        argv_found = s; 
      }

      //cout<<"\n"; 
    }
  }

}
//-----------------------------------------------------------------------||---//


//=======================================================================||===//
//================================================================| SYNCH |===//

/*
   +./src/base/cs_coupling.c
   |_'cs_coupling_discover_mpi_apps'
     |_ple_coupling_mpi_set_create(_cs_coupling_sync_flag | PLE_COUPLING_TS_MIN, 
   
   +./src/base/cs_coupling/cs_coupling_sync_apps

   +./src/base/cs_syr4_coupling.c
   |_'cs_syr4_coupling_init_comm'
     if (strncmp(op_name_recv, op_name_send, 16))

   |_'cs_syr4_coupling_init_mesh'
      if (strncmp(op_name_recv, "coupling:start", 14))

*/
void
CommDom::__coupling_sync_apps__(int      flags,
                                int      current_ts_id,
                                int     *max_ts_id,
                                double  *ts) 
{

  double _cs_coupling_ts_multiplier = 1.0;
  ple_coupling_mpi_set_t *_cs_glob_coupling_mpi_app_world = cplng; 

  if (_cs_glob_coupling_mpi_app_world != NULL) 
  {
    int i;

    int sync_flags = 0;
    int stop_mask = _cs_coupling_sync_flag & PLE_COUPLING_STOP;
    int leader_id = -1;
    double ts_min = -1.;

    double _ts = *ts * _cs_coupling_ts_multiplier;

    int n_apps
      = ple_coupling_mpi_set_n_apps(_cs_glob_coupling_mpi_app_world);
    int app_id
      = ple_coupling_mpi_set_get_app_id(_cs_glob_coupling_mpi_app_world);

    const int *app_status = NULL;
    const double *app_ts = NULL;

    ple_coupling_mpi_set_info_t ai;

    // Set synchronization flag 

    app_status
      = ple_coupling_mpi_set_get_status(_cs_glob_coupling_mpi_app_world);

    sync_flags = app_status[app_id] | flags | stop_mask;

    if (current_ts_id >= *max_ts_id)
      sync_flags = sync_flags | PLE_COUPLING_STOP;
    else {
      sync_flags = sync_flags | PLE_COUPLING_NEW_ITERATION;
      if (current_ts_id == *max_ts_id - 1)
        sync_flags = sync_flags | PLE_COUPLING_LAST;
    }

    if (flags & PLE_COUPLING_REDO_ITERATION) {
      if (sync_flags & PLE_COUPLING_NEW_ITERATION)
        sync_flags -= PLE_COUPLING_NEW_ITERATION;
      if (sync_flags & PLE_COUPLING_STOP)
        sync_flags -= PLE_COUPLING_STOP;
    }

    // Synchronize applications 

    ple_coupling_mpi_set_synchronize(_cs_glob_coupling_mpi_app_world,
                                     sync_flags,
                                     _ts);

    app_status
      = ple_coupling_mpi_set_get_status(_cs_glob_coupling_mpi_app_world);
    app_ts
      = ple_coupling_mpi_set_get_timestep(_cs_glob_coupling_mpi_app_world);

    // Check if we should use the smallest time step 

    if (app_status[app_id] & PLE_COUPLING_TS_MIN)
      ts_min = _ts;

    for (i = 0; i < n_apps; i++)
    {
      
      if(app_status[i] & PLE_COUPLING_TS_MIN) 
      {
        if(ts_min > 0) ts_min = CS_MIN(ts_min, app_ts[i]);
      }
      else
      {
        __error__(true, " into cs_coupling_discover_mpi_apps use: \'_cs_coupling_sync_flag | PLE_COUPLING_TS_MIN\' ");
      }

    }


    /* Loop on applications */
/*
    for (i = 0; i < n_apps; i++) {

      if (app_status[i] & PLE_COUPLING_NO_SYNC)
        continue;

      // Handle leader or minimum time step update 

      if (app_status[i] & PLE_COUPLING_TS_LEADER) {
        if (leader_id > -1) {
          ple_coupling_mpi_set_info_t ai_prev
            = ple_coupling_mpi_set_get_info(_cs_glob_coupling_mpi_app_world,
                                            ts_min);
          ai = ple_coupling_mpi_set_get_info(_cs_glob_coupling_mpi_app_world, i);
          bft_error
            (__FILE__, __LINE__, 0,
             _("\nApplication \"%s\" (%s) tried to set the group time step, but\n"
               "application \"%s\" (%s) has already done so."),
             ai.app_name, ai.app_type, ai_prev.app_name, ai_prev.app_type);
        }
        else {
          leader_id = i;
          *ts = app_ts[i] / _cs_coupling_ts_multiplier;
        }
      }
      else if (app_status[i] & PLE_COUPLING_TS_MIN) {
        if (ts_min > 0)
          ts_min = CS_MIN(ts_min, app_ts[i]);
      }

      // Handle time stepping behavior 

      if (app_status[i] & PLE_COUPLING_STOP) {
        if (*max_ts_id > current_ts_id) {
          ai = ple_coupling_mpi_set_get_info(_cs_glob_coupling_mpi_app_world, i);
          bft_printf
            (_("\nApplication \"%s\" (%s) requested calculation stop.\n"),
             ai.app_name, ai.app_type);
          *max_ts_id = current_ts_id;
        }
      }
      else if (app_status[i] & PLE_COUPLING_REDO_ITERATION) {
        ai = ple_coupling_mpi_set_get_info(_cs_glob_coupling_mpi_app_world, i);
        bft_error
          (__FILE__, __LINE__, 0,
           _("\nApplication \"%s\" (%s) requested restarting iteration,\n"
             "but this is not currently handled."),
           ai.app_name, ai.app_type);
      }
      else if (! (app_status[i] & PLE_COUPLING_NEW_ITERATION)) {
        ai = ple_coupling_mpi_set_get_info(_cs_glob_coupling_mpi_app_world, i);
        bft_error
          (__FILE__, __LINE__, 0,
           _("\nApplication \"%s\" (%s) synchronized with status flag %d,\n"
             "which does not specify a known behavior."),
           ai.app_name, ai.app_type, app_status[i]);
      }

      if (app_status[i] & PLE_COUPLING_LAST) {
        if (*max_ts_id > current_ts_id + 1) {
          ai = ple_coupling_mpi_set_get_info(_cs_glob_coupling_mpi_app_world, i);
          bft_printf
            (_("\nApplication \"%s\" (%s) requested last iteration.\n"),
             ai.app_name, ai.app_type);
          *max_ts_id = current_ts_id + 1;
        }
      }

    } // end of loop on applications 
*/
    if (ts_min > 0)
      *ts = ts_min / _cs_coupling_ts_multiplier;
  }


}


void
CommDom::__set_synchronization_flags() 
{

  Synch[PLE_COUPLING_INIT]           = "INIT";
  Synch[PLE_COUPLING_NO_SYNC]        = "NO_SYNC";
  Synch[PLE_COUPLING_STOP]           = "STOP";  
  Synch[PLE_COUPLING_LAST]           = "LAST";
  Synch[PLE_COUPLING_NEW_ITERATION]  = "NEW_ITERATION";
  Synch[PLE_COUPLING_REDO_ITERATION] = "REDO_ITERATION";  
  Synch[PLE_COUPLING_TS_LEADER]      = "PLE_COUPLING_TS_LEADER"; 
}


const char*
CommDom::get_synch_flags(int key) 
{
  string type = "";

  SynchIT = Synch.find(key); 
  if(SynchIT == Synch.end()) 
  {
    cout<<"\n\nERROR: ["; 
    cout<< app_type <<"] ";  
    cout<< key;
    cout<<" not found !!\n\n";
    exit(1);
  }
  else
  {
    type = Synch[key]; 
  }

  return type.c_str(); 
/*
printf("%d %d %d %d %d %d %d %d   \n", PLE_COUPLING_INIT, PLE_COUPLING_NO_SYNC, PLE_COUPLING_STOP, PLE_COUPLING_LAST,
                       PLE_COUPLING_NEW_ITERATION, PLE_COUPLING_REDO_ITERATION, 
                       PLE_COUPLING_TS_MIN , PLE_COUPLING_TS_LEADER ); 
*/
}


//=======================================================================||===//
//=====================================================| to be eliminated |===//
void
CommDom::exchange(int id)
{
/*
    // exchange 
      double *var_ij = NULL; // to be sent
      double *var_ji = NULL; // to be received
      int stride    = 1;
      int reverse   = 0;

      int n_dist_points = ple_locator_get_n_dist_points(locatorj);
      var_ij = new double[stride*n_dist_points];
      for(int i=0; i<stride*n_dist_points; i++) var_ij[i] = id+1; //(i+1);

      int n_interior    = ple_locator_get_n_interior(locatorj);
      var_ji = new double[stride*n_interior];
      for(int i=0; i<stride*n_interior;    i++) var_ji[i] = -666.66; 

      //__locator_exchange_point_var_send_double__(locator, var_ij, stride, reverse); 

      int type_size = sizeof(double);

      ple_locator_exchange_point_var(locatorj,
                                     var_ij,
                                     var_ji,
                                     NULL, // ?? 
                                     type_size, stride, reverse);


  cout<<"|_["<< app_type <<"."<< app_name <<"."<< id+1 <<"] ";
  cout<<"(to_send:"<< n_dist_points <<") "; 
  cout<<"(to_recv:"<< n_interior <<") "; 
  cout<<"\n";
*/
}


void
CommDom::save_dist_coords(int id, MPI_Comm commij)
{

  int rankij=-1;
  if( (commij != MPI_COMM_NULL) || !(commij) )
  {
    MPI_Comm_rank(commij, &rankij);
    id = rankij+1;
  }


  int appj = WhoIam[app_name];  // <-- Sooo... important!!

  std::string        filename(app_name);
  std::stringstream  srut;
  srut<< std::setfill('0') << std::setw(3) << id;
  filename += "_"+srut.str();
  filename += ".vtk";


  //int          n_interior   = ple_locator_get_n_interior(locatorj);
  //const int  *interior_list = ple_locator_get_interior_list(locatorj); 
  //cout<<"#n_interior: "<< n_interior <<" [ "; 
  //for(int i=0; i<n_interior; i++) fout<< interior_list[i] <<" ";
  //cout<<"]";

  int           n_dist_points    = ple_locator_get_n_dist_points(locatorj);
  const double   *dist_coord     = ple_locator_get_dist_coords(locatorj);
  //const int      *dist_locations = ple_locator_get_dist_locations(locatorj);

  __to_vtk__(filename, n_dist_points, dist_coord);

}


void
CommDom::__mpi_reduce__(double*         inputi,
                        double*        outputi,
                        int             n_data,
                        MPI_Comm         commi,
                        MPI_Comm        commij,
                        MPI_Op              op,
                        MPI_Datatype  datatype)

{
  /*
  MPI_MAX          maximum, max
  MPI_MIN          minimum, min
  MPI_SUM          sum
  MPI_PROD         product
  MPI_Datatype     MPI_DOUBLE, MPI_INTEGER;
  */
  int rooti = -666;
  int rootj = -666;
  __mpi_get_roots__(commi, commij, &rooti, &rootj);

  int ranki  = -666;
  MPI_Comm_rank(commi, &ranki);

  int rut = -666;
  if(ranki == 0)   // <- local   roots!!
  {
    rut   = rooti; // <- general roots!!
  }
  MPI_Bcast(&rut,   1, MPI_INTEGER, 0, commi); // all processor knows which are the general roots

  MPI_Allreduce(inputi, outputi, n_data, datatype, op, commij); // OK!!
}


void
CommDom::__mpi_reduce_min_int__(int*     inputi,
                                int*    outputi,
                                MPI_Comm  commi,
                                MPI_Comm  commij)
{

  //MPI_MAX          maximum, max
  //MPI_MIN          minimum, min
  //MPI_SUM          sum
  //MPI_PROD         product
  MPI_Op              op = MPI_MIN;
  MPI_Datatype  datatype = MPI_INTEGER;

  int rooti  = -666;
  int rootj  = -666;
  int rootij = __mpi_get_roots__(commi, commij, &rooti, &rootj);

  int ranki  = -666;
  MPI_Comm_rank(commi, &ranki);

  int varij = -666;
  int   rut = -666;
  if(ranki == 0)
  {
    rut   = rootj;
    varij = inputi[0];
  }
  MPI_Bcast(&rut,   1, MPI_INTEGER, 0, commi);
  MPI_Bcast(&varij, 1,    datatype, 0, commi);

  int rankij  = -666;
  MPI_Comm_rank(commij, &rankij);

  int result = -666;
  MPI_Reduce(&varij, &result, 1, datatype, op, 0, commij); // Why dosent work?? 
/*
if(rootij==1) 
{
  if(ranki==0)
  {
    MPI_Send(&result, 1, datatype, rootj, 0, commij);
  }
}
else
{
  if(ranki==0)
  {  
    MPI_Recv( &result, 1, datatype,     0, 0, commij, MPI_STATUS_IGNORE);  
  }
}
*/

  if(ranki==0)
  {
    if(rankij==0) MPI_Send(&result, 1, datatype, rootj, 0, commij);
    else          MPI_Recv( &result, 1, datatype,     0, 0, commij, MPI_STATUS_IGNORE);
  }
  MPI_Bcast(&result, 1, datatype,      0, commi);

cout<<" ---------> ";
cout<< rankij <<" ";
//cout<< ranki  <<" ";
//cout<< varij  <<" ";
cout<< result <<" ";
//cout<< rooti <<" ";
//cout<< rootj <<" ";
//cout<< rootij <<" ";
//cout<< sizeij <<" ";
cout<<" \n";

  outputi[0] = result;

/* GOOD 
  MPI_Bcast(&rut,   1, MPI_INTEGER, 0, commi); 
  MPI_Bcast(&varij, 1,    datatype, 0, commi);

  int result = -666; 
//MPI_Reduce(&varij, &result, 1, datatype, op, rut, commij); // Why dosent work?? 
  MPI_Allreduce(&varij, &result, 1, datatype, op, commij); 

  outputi[0] = result; 
*/
}


void
CommDom::__mpi_reduce_sum_double__(double*    inputi,
                                   double*   outputi,
                                   MPI_Comm    commi,
                                   MPI_Comm   commij)
{
  int             n_data = 1;
  MPI_Op              op = MPI_SUM;
  MPI_Datatype  datatype = MPI_DOUBLE;

  __mpi_reduce__( inputi,
                  outputi,
                  n_data,
                  commi,
                  commij,
                  op,
                  datatype);

}


void
CommDom::__mpi_reduce_min_double__(double*    inputi,
                                   double*   outputi,
                                   MPI_Comm    commi,
                                   MPI_Comm   commij)
{
  int             n_data = 1;
  MPI_Op              op = MPI_MIN;
  MPI_Datatype  datatype = MPI_DOUBLE;

  __mpi_reduce__( inputi,
                  outputi,
                  n_data,
                  commi,
                  commij,
                  op,
                  datatype);

}


void
CommDom::__mpi_gather_double__(double*    inputi,
                               double*   rbuf,
                                  int*       num,
                               MPI_Comm    commi,
                               MPI_Comm   commij)
{
  /*
  MPI_MAX          maximum, max
  MPI_MIN          minimum, min
  MPI_SUM          sum
  MPI_PROD         product
  */
  MPI_Op              op = MPI_MIN;
  MPI_Datatype  datatype = MPI_DOUBLE;

  double  varij = -666;
  //double result = -666;
  //double*  rbuf =  NULL; 

  int rooti = -666;
  int rootj = -666;
  __mpi_get_roots__(commi, commij, &rooti, &rootj);

  int ranki   = -666;
  int rankij  = -666;
  MPI_Comm_rank(commi, &ranki);

  int sizei   = -666;
  int sizeij  = -666;
  MPI_Comm_size(commi, &sizei);


  int rut = -666;
  if(ranki == 0) // <- local root!!
  {
    rut   = rooti;
    varij = inputi[0];
    for(int i=0; i<sizei; i++) rbuf[i] = 0.0;
  }
  MPI_Bcast( &rut,    1, MPI_INTEGER, 0, commi);

  MPI_Gather( inputi, 1, datatype,
                rbuf, 1, datatype,
                                0,
                            commi);

}

#if SATURNE
void
CommDom::__alya2cs_type__(
                            int*   alya_types,
                            int*     cs_types,
                            int      n_elems)
{
  int alya_type = -1;
  int   cs_type = -1;
  for(int i=0; i<n_elems; i++) 
  {
   cs_types[i] = -1;  
   __get_alya2cs_type__(alya_types[i], cs_types+i); 
/*
cout<< i+1 <<") ";
cout<<"cs_type:"<< cs_types[i] <<" ";  
cout<<"alya_type:"<< alya_types[i] <<" ";
cout<<"\n";
*/
  } 
}


void 
CommDom::__get_alya2cs_type__(int alya_type, int* cs_type) 
{
  Alya2CsIT = Alya2Cs.find(alya_type);
  if(Alya2CsIT == Alya2Cs.end())
  {
    cout<<"WARNINIG: [alya_type:"<< alya_type <<"] not found!! \n";
    __error__(true, "interaction doesnt exist");
  }

  cs_type[0] = Alya2Cs[ alya_type ];
}


//=======================================================================||===//
//=======================================================================||===//
/*
|
|_ mesh/cs_mesh.h
 |_ cs_mesh_t 
   |_n_cells;                     // Number of cells 
   |_n_i_faces;                   // Number of interior faces 
   |_n_b_faces;                   // Number of boundary faces 
   |_n_vertices;                  // Number of vertices 
   |_*vtx_coord;                  // Vertex coordinates
   |_*i_face_cells;               // Interior faces -> cells connectivity
   |_*b_face_cells;               // Boundary faces -> cells connectivity
   |_*i_face_vtx_idx;             // Interior faces -> vertices index 
   |_*i_face_vtx_lst;             // Interior faces -> vertices connectivity 
   |_*b_face_vtx_idx;             // Boundary faces -> vertices index 
   |_*b_face_vtx_lst;             // Boundary faces -> vertices connectivity 

|
|_ base/cs_syr4_coupling
  |_ typedef cs_syr4_coupling_t
       cs_syr4_coupling_ent_t  *faces
       cs_syr4_coupling_ent_t  *cells
         fvm_nodal_t             *elts
         fvm_nodal_section_t     **sections;
|
|_ base/cs_syr4_coupling.cs_syr4_coupling_init_mesh
  |_ base/cs_syr4_coupling._create_coupled_ent
    |
    |_ /mesh/cs_mesh_connect.cs_mesh_connect_[cells|faces]_to_nodal -> fvm_nodal_t
      |_ cs_mesh_connect_get_cell_faces(mesh,                       -> cell_face_idx, cell_face_num
      |_ fvm/fvm_nodal.fvm_nodal_create                             -> extr_mesh [fvm_nodal_t]  
      |_ fvm/fvm_nodal_from_desc.fvm_nodal_from_desc_add_cells      
        |_ fvm_nodal_t                                              <- extr_mesh 
        |_ extr_cell_count
        |_ NULL                                  \
        |_ 2                                     |
        |_ [0, n_b_faces, n_b_faces + n_i_faces] |--__related to cell_type  
        |_ [b_face_vtx_idx, i_face_vtx_idx]      |--
        |_ [b_face_vtx_lst, i_face_vtx_lst]      |
        |_ cell_face_idx                         | 
        |_ cell_face_num                         /
        |_ cell_family = NULL??                  ??
        |_ cell_list   = NULL??                  ??
        |_ &polyhedra_faces = NULL??             |--FVM_CELL_POLY
|_ face_families = NULL?? 
|_ NULL
      |_ fvm/fvm_nodal_set_shared_vertices
      |_ fvm/fvm_nodal_set_group_class_set
      |_ fvm/fvm_nodal_order_cells
      |_ fvm/fvm_nodal_init_io_num
      |_ fvm/fvm_nodal_order_vertices
      |_ fvm/fvm_nodal_init_io_num
    |
    |_ location_elts = coupling_ent->elts; <- fvm_nodal_t
    |_ ple_locator_set_mesh(..., location_elts, ...)
|
|_
  |_fvm/fvm_nodal_from_desc.fvm_nodal_from_desc_add_cells
|_src/base/cs_coupling.c
  |
  |_cs_coupling_mesh_extents
    |_fvm_nodal_extents(const void  *mesh, <-- (const fvm_nodal_t *)mesh
                        double      tolerance,
                        double      extents[])
  |
  |_cs_coupling_point_in_mesh
    |_fvm_point_location_nodal(const void         *mesh,
                               double              tolerance,
                               ple_lnum_t          n_points,
                               const ple_coord_t   point_coords[],
                               ple_lnum_t          location[],
                               float               distance[])
  |
  |_cs_coupling_point_closest_mesh
    |_fvm_point_location_closest_nodal(const void        *mesh, 
                                       ple_lnum_t         n_points 
                                       const ple_coord_t  point_coords[], 
                                       ple_lnum_t         location[], 
                                       float              distance[])

      |_ fvm/fvm_nodal_set_shared_vertices
      |_ fvm/fvm_nodal_set_group_class_set
      |_ fvm/fvm_nodal_order_cells
      |_ fvm/fvm_nodal_init_io_num
      |_ fvm/fvm_nodal_order_vertices
      |_ fvm/fvm_nodal_init_io_num

*/
void
CommDom::__create_nodal_t__(void*  coords,
                            int    n_coords,
                            int*   cs_elements,
                            int*   cs_types,
                            int    n_elems)
{
  //----------------------------------------------------------------------||--//
  //----------------------------------------------------------------------||--//
  fvm_nodal_t  *CsMesh = NULL;
//CsMesh = fvm_nodal_create("COMMDOM", __DIM__);
  CsMesh = fvm_nodal_create("COMMDOM", dimension );
  CsMesh->n_cells    = n_elems;
  CsMesh->n_sections = 0;

  CsMesh->_vertex_coords    = NULL;
  CsMesh->parent_vertex_num = NULL;
  CsMesh->global_vertex_num = NULL;
  //----------------------------------------------------------------------||--//

  //----------------------------------------------------------------------||--//
  map<int,fvm_element_t> CsType;
  //CsType[0] = FVM_EDGE;
  CsType[1] = FVM_FACE_TRIA;
  CsType[2] = FVM_FACE_QUAD;
  //CsType[3] = FVM_FACE_POLY; 
  CsType[4] = FVM_CELL_TETRA;
  CsType[5] = FVM_CELL_PYRAM;
  CsType[6] = FVM_CELL_PRISM;
  CsType[7] = FVM_CELL_HEXA;
  //CsType[8] = FVM_CELL_POLY;
  //CsType[9] = FVM_N_ELEMENT_TYPES;

  map<int,string> Cs2Alya;
  Cs2Alya[1] = "TRI03"; // 10; //TRI03 = 10, FVM_FACE_TRIA
  Cs2Alya[2] = "QUA04"; // 12; //QUA04 = 12, FVM_FACE_QUAD
  Cs2Alya[4] = "TET04"; // 30; //TET04 = 30, FVM_CELL_TETRA 
  Cs2Alya[5] = "PYR05"; // 32; //PYR05 = 32, FVM_CELL_PYRAM
  Cs2Alya[6] = "PEN06"; // 34; //PEN06 = 34, FVM_CELL_PRISM
  Cs2Alya[7] = "HEX08"; // 37; //HEX08 = 37, FVM_CELL_HEXA
  //----------------------------------------------------------------------||--//

  //----------------------------------------------------------------------||--//
  vector<int> cs_stride(n_elems+1, 0);

  vector<int> n_elements_type(FVM_N_ELEMENT_TYPES, 0);
  for(int i = 0; i < CsMesh->n_cells; i++)
  {
    int cell_idx = cs_types[i];
    n_elements_type[cell_idx] += 1;

    int stride = fvm_nodal_n_vertices_element[cell_idx];
    cs_stride[i+1] = stride + cs_stride[i];
  }
  //----------------------------------------------------------------------||--//

  //----------------------------------------------------------------------||--//  
  CsMesh->n_vertices = n_coords;
  //CsMesh->_vertex_coords = (cs_coord_t*) malloc( CsMesh->n_vertices*CsMesh->dim );
  CsMesh->_vertex_coords = new cs_coord_t[CsMesh->n_vertices*CsMesh->dim];
  memset(CsMesh->_vertex_coords, 0, CsMesh->n_vertices*CsMesh->dim);
  CsMesh->vertex_coords = CsMesh->_vertex_coords;

  cs_coord_t* l_coords = (cs_coord_t*) coords;
  for(int i=0; i<CsMesh->n_vertices*CsMesh->dim; i++) CsMesh->_vertex_coords[i] = l_coords[i];
  //----------------------------------------------------------------------||--//

  //----------------------------------------------------------------------||--//
  map<fvm_element_t, int>            TypeId;
  map<fvm_element_t, int>::iterator  TypeIdIT;
  for(int type_id = 0, idx=0; type_id < FVM_N_ELEMENT_TYPES; type_id++)
  {
    if(n_elements_type[type_id]>0)
    {
      fvm_element_t cell = CsType[type_id];
      TypeId[ cell ] = idx++;
    }
  }

  vector<int> cs_stride_sections(TypeId.size()+1, 0);

  CsMesh->n_sections = TypeId.size();
  CsMesh->sections = (fvm_nodal_section_t**) malloc( CsMesh->n_sections );
  for(TypeIdIT = TypeId.begin(); TypeIdIT != TypeId.end(); TypeIdIT++)
  {
    fvm_element_t cell = TypeIdIT->first;
    int           idx  = TypeIdIT->second;
    int     n_elements = n_elements_type[ cell ];

    CsMesh->sections[idx]                    = fvm_nodal_section_create( cell );
    CsMesh->sections[idx]->n_elements        = n_elements;
    cs_stride_sections[idx+1]                = cs_stride_sections[idx] + n_elements;

    int stride            = CsMesh->sections[idx]->stride;
    int connectivity_size = stride * n_elements;
    CsMesh->sections[idx]->connectivity_size = connectivity_size;

    
    ////CsMesh->sections[idx]->_parent_element_num = (cs_lnum_t*) malloc( n_elements ); 
    //CsMesh->sections[idx]->_parent_element_num = new cs_lnum_t[ n_elements ]; 
    //CsMesh->sections[idx]->parent_element_num = CsMesh->sections[idx]->_parent_element_num;
    //for(int i = 0; i < n_elements; i++)
    //  CsMesh->sections[idx]->_parent_element_num[i] = -(i+1);
    

    //CsMesh->sections[idx]->_vertex_num = (cs_lnum_t*) malloc( CsMesh->sections[idx]->connectivity_size ); 
    CsMesh->sections[idx]->_vertex_num = new cs_lnum_t[ CsMesh->sections[idx]->connectivity_size ];
    CsMesh->sections[idx]->vertex_num  = CsMesh->sections[idx]->_vertex_num;
    for(int i = 0; i < CsMesh->sections[idx]->connectivity_size; i++)
      CsMesh->sections[idx]->_vertex_num[i] = -(i+1);

    n_elements_type[cell] = 0;
  }
  //----------------------------------------------------------------------||--//

  //----------------------------------------------------------------------||--//
  fvm_element_t         cell;
  fvm_nodal_section_t*  section = NULL;

  //vector<int> 
  parent_element_num = new int[CsMesh->n_cells];
  for(int i=0, idx=0; i < CsMesh->n_cells; i++)
  {
    cell = (fvm_element_t) cs_types[i];
    idx  = TypeId[cell];

    section = CsMesh->sections[idx];

    int s_stride = cs_stride_sections[idx];
    parent_element_num[ s_stride+n_elements_type[cell] ] = i + 1;
    //section->_parent_element_num[ n_elements_type[cell] ] = i + 1;

    int l_stride = section->stride*n_elements_type[cell];
    int g_stride = cs_stride[i];
    for(int j=0; j<section->stride; j++) section->_vertex_num[l_stride+j] = cs_elements[g_stride + j];

    n_elements_type[cell] += 1;
  }

  //----------------------------------------------------------------------||--//
  if(local_size<=4)
  { 
    cout<<"|_["<< app_type <<"."<< app_name <<"] ";
    for(TypeIdIT = TypeId.begin(); TypeIdIT != TypeId.end(); TypeIdIT++)
    {
      fvm_element_t cell = TypeIdIT->first;
      cout<< Cs2Alya[cell] <<":"<< n_elements_type[cell] <<" ";  
    }
    cout<<"\n";
  }
  //----------------------------------------------------------------------||--//

  //----------------------------------------------------------------------||--//  
  mesh_extents_f   = cs_coupling_mesh_extents;
  locate_inside_f  = cs_coupling_point_in_mesh;
  locate_closest_f = cs_coupling_point_closest_mesh;

//double *extents = new double[2*__DIM__];
//memset(extents, 0.0, 2*__DIM__);
  double *extents = new double[2*dimension];
  memset(extents, 0.0, 2*dimension);
  mesh_extents_f(CsMesh, 1, 1e-3, extents);

  if(local_size<=4)
  {
    cout<<"|_["<< app_type <<"."<< app_name <<"] ";
    cout<<"Range: ";
  //for(int i=0; i<__DIM__; i++) cout<<"x"<< i <<":["<< extents[i] <<","<< extents[i+__DIM__] <<"] ";
    for(int i=0; i<dimension; i++) cout<<"x"<< i <<":["<< extents[i] <<","<< extents[i+dimension] <<"] ";
    cout<<"\n";
  }

  int *location = new int[n_coords];
  memset(location, 0.0, n_coords);

  float *distance = new float[n_coords];
  memset(distance, 0.0, n_coords);

  locate_inside_f(CsMesh, 1e-3,
                  n_coords, (cs_coord_t*) coords,
                  location, distance);

  
  //fvm_point_location.c:3648: Fatal error.
  //Locating volume elements closest to points not handled yet
  //
  //locate_closest_f(CsMesh,
  //                n_coords, (cs_coord_t*) coords,
  //                 location, distance);
  
  //----------------------------------------------------------------------||--//
  //----------------------------------------------------------------------||--//

      locate_closest_f = NULL;
      const int *point_list = NULL;
      ple_locator_set_mesh(locatorj,         //*ple_locator_t
                           CsMesh,           // void*
                           dimension;        // __DIM__,          // coord_dim!! 
                           n_coords_dist,    // ple_lnum_t 
                           point_list,       // const ple_lnum_t[] 
                           coords_dist,      // const ple_coord_t[]
                           distance_dist,    // float[]
                           mesh_extents_f,   // ple_mesh_extents_t.          cs_coupling_mesh_extents
                           locate_inside_f,  // ple_mesh_elements_contain_t. cs_coupling_point_in_mesh. 
                           locate_closest_f);// ple_mesh_elements_closest_t. 
                                             // locate_on_closest = cs_coupling_point_closest_mesh. optional 
//ple_locator_dump(locatorj); 
  //----------------------------------------------------------------------||--//

  //----------------------------------------------------------------------||--//
  if(world_size<=8)
  {
        cout<<"|_["<< app_type <<"."<< app_name <<"] ";
        cout<<"n_cells: ";
        cout<< CsMesh->n_cells <<", ";;
        cout<<"n_pts: ";
        cout<< n_coords_dist <<", ";;
        cout<<"n_exterior: ";
        cout<< ple_locator_get_n_exterior(locatorj)<<", ";
        cout<<"n_interior: ";
        cout<< ple_locator_get_n_interior(locatorj)<<", ";
        cout<<"n_dist_points: ";
        cout<< ple_locator_get_n_dist_points(locatorj);
        cout<<"\n";
  } 
  //----------------------------------------------------------------------||--//
  //fvm_nodal_dump(CsMesh); 
  fvm_nodal_destroy(CsMesh);
  //----------------------------------------------------------------------||--//
}

void
CommDom::__get_parent__(int* _parent_element_num, int n_cells)
{
  for(int i=0; i < n_cells; i++) _parent_element_num[i] = parent_element_num[i];
  delete[] parent_element_num;
}


void
CommDom::__interpolation__(
                           int      cs_type,
                           double       tol,
                           double*   coords,
                           int*    vertices,
                           double*    point,
                           double*    shapef
                           )
{

  //CommDom::__tetra_interpolation__(double* coords, int* vertices, double* point, double* shapef)
  //vm_element_t  cell = (fvm_element_t) cs_type; 
/*
                                      int*   alya_type,
                                      double*       tol,
                                      double* coords_i, int* vertices_i,
                                      double* point_j,  double* shapef_j)

*/

  get_shape_func((fvm_element_t) cs_type,
                  vertices,
                  (cs_coord_t*) coords,
                  (cs_coord_t*) point,
                  tol,
                  shapef);

}
#endif // SATURNE 


//=======================================================================||===//
//=======================================================================||===//
