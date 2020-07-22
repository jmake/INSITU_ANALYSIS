/*
  2019DIC02. KTH, STOCKHOLM, SWEDEN. Migue Zavala 
  FROM   : 
           /home/bsc21/bsc21704/z2019_2/REPOSITORY/HEMELB02_01/ToUCL02_01/HEMELB/EXAMPLES/CPLNG03_1
           ../../Plepp04_1/plepp.hpp
*/

#ifndef INSITU_COMMDOM_H
#define INSITU_COMMDOM_H

#include "read_file.hpp"
#include "commdom.hpp"

#include <map>
#include <cmath> 

/* 
//--------------------------------------------------------------------------||--//
const char* String_f2c(const char* licName, int nameLen)
{
  char* name;
  name = (char*)malloc(sizeof(char)*(nameLen));
  for(int i = 0; i<nameLen; i++) name[i] = licName[i];
  name[nameLen] = '\0';

  return name ;
}


//--------------------------------------------------------------------------||--//
int n_node = 4; // tetras:4, tria:3 

void
propi2propj(
            std::vector<double>  props_i, 
            std::vector<int>     vertices, 
            int                  n_dist_j, 
            std::vector<int>     elemts_i, 
            std::vector<double>  tetracoords_j, 
            std::vector<double> &props_j 
           )
{
  int    ii, ielem;
  int      vertices_i[n_node];
  double vol_coords_j[n_node];
  double       prop_j[n_node];

  if(n_dist_j > 0)
  {
    for(int ii=0; ii<n_dist_j; ii++ )
    {
      ielem = elemts_i[ii] - 1; //  <- Fortran style ?? 
      for(int jj=0; jj<n_node; jj++)   vertices_i[jj]  =      vertices[n_node*ielem + jj] - 1; // <- Fortran style -> C style !!  
      for(int jj=0; jj<n_node; jj++) vol_coords_j[jj]  = tetracoords_j[n_node*ii    + jj];
      for(int jj=0; jj<n_node; jj++)       prop_j[jj]  =       props_i[ vertices_i[jj]  ];
//for(auto x:prop_j) std::cout<< x <<" ";  
      props_j[ii] = 0.0;
      for(int jj=0; jj<n_node; jj++)      props_j[ii] += prop_j[jj] * vol_coords_j[jj];
//std::cout<< props_j[ii] <<" ";  
    }
  }
}
*/

//----------------------------------------------------------------------||---//
class ReadFile
{ 
  public :
 ~ReadFile(){}
  ReadFile(){}
  
  template<class T>
  std::vector<T> LoadFile(std::string fname, int dim, int print)
  { 
    read_log_file DATA;
    DATA.set_name( fname );
    DATA.run();
    
    vector< vector<double> > vdata( DATA.get_vdata() );
    return PrintVector<T>( vdata, dim, 1, print);
  }

  
  template<class T>
  std::vector<T> PrintVector( vector< vector<double> > data, int dim, int init_col=0, int print=0)
  { 
    std::vector<T> array;  
    for(int i=0,k=0; i<data.size(); i++)
    { 
      vector<double> row( data[i] );
      if(print) std::cout<<" "<< i ;  
      for(int j=init_col; j<init_col+dim; j++, k++)
      { 
        if(print) std::cout<<" "<< row[j] ;
        array.push_back( row[j] );
      }
      if(print) std::cout<<" \n";
    }
    
    return array;
  }
};


//--------------------------------------------------------------------------||--//
class AnalyseArgvs 
{
  public : 
 ~AnalyseArgvs()
  {
    CommByLocalValue.clear();
    Dic.clear(); 
    argvs.clear();  
  };


  AnalyseArgvs(){}; 


  void __set_argvs__(string dummy)
  {
    argvs.push_back(dummy);
  }

  
  void __print_argvs__(std::string sms="\t[__print_argvs__]") 
  {
    for(const auto& item : Dic) 
    {
      std::vector<std::string> V( item.second ); 
      cout<< sms <<" Arg('"<< item.first <<"'):[";
      for (auto v : V) std::cout<<"'"<<  v <<"' ";     
      cout<<"] ";
      cout<<"\n";
    }
  }


  void __analyse_argvs__(string _token, int print)
  {

    for(int i=0; i<argvs.size(); i++) 
    {
      string argv = argvs[i]; 

      if(argv.find(_token) != std::string::npos) 
      {
        if(std::strncmp(argv.c_str(), _token.c_str(), _token.size()) == 0)
        {
          //std::cout<<"'"<< argv <<"' Found! \n";
          string next( (i+1>=argvs.size())?("NoFound!"):(argvs[i+1]) );
	  Dic[_token].push_back( next ); 
        } 
      }

      token = _token;
      if(print) __print_argvs__() ;
    } 
  }


  string __get_argvs__()
  {
    std::string found; 

    if(Dic.find(token) != Dic.end())
    {
      std::vector<std::string> V( Dic[token] );
      if(V.size() > 1) {__print_argvs__("\tERROR (nArgs>1)"); exit(0);} 
      found = V[0];   
    }

    return found;  
  }


  MPI_Comm SimpleSplit(MPI_Comm gcomm, int print) 
  {
    // X.0.
    assert(gcomm != MPI_COMM_NULL);  
    MPI_Barrier(gcomm);

    // X.0. In each rank, values 'LocalValue' are collected (Gathered) in 'AllValues' 
    std::vector<std::string> LocalValue = Dic[token]; 
    std::vector<std::string> AllValues( AllGatherVString(LocalValue,gcomm) ) ;
    //for(auto value:AllValues) std::cout<<" '"<< value <<"' ";

    // X.0. Ranks with same 'LocalValue's are identified in 'RanksByLocalValue' 
    std::map<std::string, std::vector<int> > RanksByLocalValue = GetDicFromVector3( AllValues );

    // X.0. Creating ALL local comms by using 'RanksByLocalValue'  
    for(const auto& item : RanksByLocalValue)
    {
      MPI_Comm lcomm = GroupComm(gcomm, item.second); 
      CommByLocalValue[item.first] = lcomm; 
      /*  
      int nranks = -1;
      if(lcomm != MPI_COMM_NULL) MPI_Comm_size(lcomm, &nranks); 

      cout<<" '"<< item.first <<"' ("<< nranks << ") :[";
      for (auto v : item.second) std::cout<<"'"<<  v <<"' ";
      cout<<"] ";
      cout<<"\n";
      */ 
    }

    // X.0. Extracting local comm 'lcomm'  
    assert(Dic.find(token) != Dic.end());  
	    
    std::string   key = Dic[token][0];  
    MPI_Comm    lcomm = CommByLocalValue[key]; 

    if(print && (lcomm != MPI_COMM_NULL) )
    { 
      int nranks = -1, irank=-1; 
      MPI_Comm_size(lcomm, &nranks);
      MPI_Comm_rank(lcomm, &irank );

      cout<<" '"<< key <<"' ("<< nranks <<"."<< irank << ") Ranks:[ ";
      for(auto v : RanksByLocalValue[key]) std::cout<< v <<" ";
      cout<<"] \n";
    }
    AllValues.clear();
    RanksByLocalValue.clear();

    return lcomm;  
  }


  std::vector<std::string> AllGatherVString(std::vector<std::string> sendBufferAux, MPI_Comm comm)
  {
    // X.0.
    std::vector<char> sendBuffer = string_to_char(sendBufferAux);  

    // X.0.
    MPI_Datatype Tmpi = MPI_CHAR; 

    int nranks = -1;  
    MPI_Comm_size(comm, &nranks);

    int sendLength = sendBuffer.size();
    std::vector<int> recvLengths(nranks, -1);

    MPI_Allgather(&sendLength,     1, MPI_INTEGER,
                  &recvLengths[0], 1, MPI_INTEGER, comm);

    std::vector<int> offsets(nranks + 1, 0);
    for (int i = 0; i < nranks; ++i) offsets[i + 1] = recvLengths[i] + offsets[i] ;

    std::vector<char> recvBuffer( offsets[nranks] );
    if( recvBuffer.size() )
    {
      MPI_Allgatherv(&sendBuffer[0],  sendLength,                  Tmpi,
                     &recvBuffer[0], &recvLengths[0], &offsets[0], Tmpi, comm);
    }

    // X.0. char_to_string   
    std::vector<std::string> recvBufferAux;   
    for (int i = 0; i < nranks; ++i) 
    {
      int b = offsets[i + 0]; 
      int t = offsets[i + 1];

      std::string s; 
      for(int j=b; j<t; j++) s+= recvBuffer[j];
      recvBufferAux.push_back(s);  
    }

    return recvBufferAux;
  }


  MPI_Comm GroupComm(MPI_Comm orig_comm, std::vector<int> subranks) 
  {
    MPI_Group orig_group = MPI_GROUP_EMPTY;  
    MPI_Comm_group(orig_comm, &orig_group);  

    MPI_Group subgroup   = MPI_GROUP_EMPTY;  
    MPI_Group_incl(orig_group, subranks.size(), &(subranks[0]), &subgroup);   
 
    MPI_Comm subcommunicator = MPI_COMM_NULL;  
    MPI_Comm_create(orig_comm, subgroup, &subcommunicator);  
    return subcommunicator;  
  }



  MPI_Comm SpliComm(int color, int rank, MPI_Comm worl_comm=MPI_COMM_WORLD)
  {
    MPI_Comm row_comm = MPI_COMM_NULL;
    MPI_Comm_split(worl_comm, color, rank, &row_comm);
    return row_comm;
  }



  std::vector<char> string_to_char(const std::vector<std::string>& strings)
  {
    std::vector<char> cstrings;
    cstrings.reserve(strings.size());
    for(std::string s: strings)
      for(size_t i = 0; i < strlen(s.c_str()); ++i)  cstrings.push_back(s.c_str()[i]);

    return cstrings;
  }


  template <class T>
  std::map<T, std::vector<int> >
  GetDicFromVector3( std::vector<T> vector )
  {
    std::map<T, std::vector<int> > dic;
    for(int i=0; i<vector.size(); i++) dic[ vector[i] ].push_back(i);
    return dic;
  }



  private :
  vector<string>  argvs;
  string          argv_found;
  string          token;

  std::map<std::string, std::vector<std::string> > Dic;
  std::map<std::string, MPI_Comm> CommByLocalValue;

};  


//--------------------------------------------------------------------------||--//
class Coupling 
{
  public:
 ~Coupling()
  {
    CD.locator_destroy(); 
  } 


  Coupling(std::string type, std::string namei, std::string namej) 
  {
    IAM     =  type;  
    NAMEi   =  namei; 
    NAMEj   =  namej; 

    global_comm = MPI_COMM_NULL; //MPI_COMM_WORLD; 
     local_comm = MPI_COMM_NULL;
         commij = MPI_COMM_NULL;
  }


  MPI_Comm Init(MPI_Comm gcomm) 
  {
    assert(gcomm); 
    global_comm = gcomm;

    CD = CommDom();
    CD.init(); 
    CD.set_app_type(IAM); 
    CD.set_world_comm( global_comm ); 
    CD.set_app_name(NAMEi); 

    local_comm = CD.set_mpi_comms();
    MPI_Comm_rank(local_comm, &local_rank);

    commij     = CD.get_mpi_commij(NAMEj);
    return local_comm ;  
  }


  MPI_Comm GetCommij() 
  {
    assert( commij );
    return commij; 
  }


  int LRoot()
  {
    return !local_rank;   
  }


  void GBarrier()
  {
    MPI_Barrier(global_comm);
  }


  void sendrecv_int(int* send, int n_send, int* recv, int n_recv, int print=0)   
  {
    assert( commij ); 
    assert( local_comm ); 

    CD.__mpi_sendrecv_int__(send, n_send, recv, n_recv, local_comm, commij); 
    MPI_Bcast(recv, n_recv, MPI_INT, 0, local_comm);

    if(!local_rank && print)
    std::cout<<" '"<< NAMEi <<"."<< local_rank+1 <<"'"<<
               " Send:'"    << send[0]           <<"'"<<
               " Recv:'"    << recv[0]           <<"' \n";
  }


  void sendrecv_double(double* send, int n_send, double* recv, int n_recv, int print=0)
  {
    assert( commij );
    assert( local_comm );

    CD.__mpi_sendrecv_real__(send, n_send, recv, n_recv, local_comm, commij);
    MPI_Bcast(recv, n_recv, MPI_DOUBLE, 0, local_comm);

    if(!local_rank && print)
    std::cout<<" '"<< NAMEi <<"."<< local_rank+1 <<"'"<<
               " Send:'"    << send[0]           <<"'"<<
               " Recv:'"    << recv[0]           <<"' \n";
  }


  std::string sendrecv_string(std::string csend, int print=0) 
  {
    // Send|Recv 'int'  
  //int send[1] = { csend.size() };  
    int recv[1] = { -1 };
    int send[1]; // = { csend.size() };
    send[0] = csend.size(); 

    sendrecv_int(send, 1, recv, 1);  

    // Send|Recv 'string' 
    std::vector<char> aux_send( csend.begin(), csend.end() );
    std::vector<char> aux_recv( recv[0], '-' );

    CD.__mpi_sendrecv_char__(aux_send.data(), aux_send.size(),
                             aux_recv.data(), aux_recv.size(),
                             local_comm, commij);
    MPI_Bcast(aux_recv.data(), aux_recv.size(), MPI_CHAR, 0, local_comm);

    if(!local_rank && print)
    std::cout<<" '"<< NAMEi <<"."<< local_rank+1 <<"'"<<
               " Send:'"    << aux_send.data()   <<"'"<<
               " Recv:'"    << aux_recv.data()   <<"' \n";

    return std::string(aux_recv.data());  
  }


  protected :
  std::map<std::string,MPI_Comm>            Commij;
  std::map<std::string,MPI_Comm>::iterator  It;

  MPI_Comm       commij; 
  MPI_Comm   local_comm; 
  MPI_Comm  global_comm; 

  string IAM   ; 
  string NAMEi ;
  string NAMEj ;  

  CommDom  CD;  

  int local_rank ; 
  int n_time_steps[1] ;

}; // Coupling     



//--------------------------------------------------------------------------||--//
class Locatization : public Coupling 
{
  public : 
 ~Locatization(){};


  Locatization(std::string type, std::string namei, std::string namej, int dim) : Coupling(type, namei, namej)
  {
    CELL_DIM = dim; 
  }


  void CreateMesh(vector<double> ple_pts, vector<int> ple_cells, int CELL_SIZE, int CELL_TYPE)
  {
    assert(commij);
    assert(commij != MPI_COMM_NULL);  
    assert(local_comm);

    int            n_vertices_i = 0;
    double* vertex_coords_i_ptr = NULL;

    int            n_vertices_j = 0;
    double* vertex_coords_j_ptr = NULL; 
    double*  vertex_props_j_ptr = NULL; //MATMATMAT 

    int       n_elements_i = 0;
    int*  vertex_num_i_ptr = NULL;
    int* vertex_type_i_ptr = NULL;

           n_vertices_i = ple_pts.size() / CELL_DIM; 
    vertex_coords_i_ptr = ple_pts.data(); 

           n_vertices_j = n_vertices_i;
    vertex_coords_j_ptr = vertex_coords_i_ptr; 

           n_elements_i = ple_cells.size() / CELL_SIZE; 
    vertex_num_i_ptr    = ple_cells.data();    
      vertex_type_i_ptr = NULL; 

  //std::vector<double> dummy( n_vertices_i, 0.0); 
  //vertex_props_j_ptr = &dummy[0];

    CD.locator_create2(local_comm, commij, 1e-3, CELL_TYPE);
    CD.locator_set_cs_mesh(n_vertices_i,
                           n_elements_i,
                           vertex_coords_i_ptr,
                           vertex_num_i_ptr,
                           vertex_type_i_ptr, 
                           n_vertices_j,
                           vertex_coords_j_ptr, CELL_DIM, //NULL, 0); 
                           vertex_props_j_ptr, 0); 

    CD.save_dist_coords(0, local_comm);

    // ToSend     
    n_send = CD.get_n_dist_points();
    dist_locations_i = vector<int   >(n_send           ,       -1); // Touched Local Elements 
    dist_coords_j    = vector<double>(n_send * CELL_DIM, HUGE_VAL);
  //dist_props_j     = vector<double>(n_send * 1, HUGE_VAL);  
    CD.__locator_get_dist_locations__( dist_locations_i.data() );
    CD.__locator_get_dist_coords__(       dist_coords_j.data() );
  //CD.__locator_get_dist_props__(         dist_props_j.data() ); 

    // ToRecv   
    n_recv = CD.get_n_interior();
    interior_list_j = vector<int>(n_recv, -1);
    CD.__locator_get_interior_list__( interior_list_j.data() ); // Touched Local Points 

    vertex_coords_i = ple_pts; 
    //vertex_num_i = ple_cells; // ??  

  }


  void CheckInput(int      n_vertices_i,             //  points.size() / DIM; 
                  int      n_elements_i,             //  types.size() 
                  double* _vertex_coords_i_ptr,      //  points.data() 
                  int*    _vertex_num_i_ptr,         //  connectivity.data() 
                  int*     vertex_type_i_ptr,        //  types.data()  
                  int     _n_vertices_j,             //  n_vertices_i 
                  double*  vertex_coords_j_ptr,      // _vertex_coords_i_ptr 
		  int      DIM,  
                  double*  vertex_props_j_ptr,       //  dummy.data()  
		  int      DIM_PROPS,                // 
		  int print = 1 
		 )
  {
    if(print) 
    {
      std::cout<<" \n";

      std::cout<<"\t[CreateMesh]"
               <<" nPts:"<< n_vertices_i
               <<" nCells:"<< n_elements_i
               <<" Dim:" << DIM
               <<" \n";

      int cell_size = 3;  
      for(int i=0,k=0; i<n_elements_i; i++)  
      {
        int t = vertex_type_i_ptr[i]; 
        std::cout<<" ielem:"<< i <<" type:"<< t ; 

        std::cout<<" elem:[";
        for(int j=0; j<cell_size; j++, k++)
        {
          int icell = _vertex_num_i_ptr[k] - 1; 
          std::cout<< icell <<":";

          double *pt = _vertex_coords_i_ptr + icell * DIM;    
          std::cout<<"(";
          for(int l=0; l<DIM; l++) std::cout<<" "<< pt[l]; 
          std::cout<<" ) ";

        }
        std::cout<<"] \n"; 
      } 

      std::cout<<" \n";
      MPI_Barrier(global_comm);
    } // if 
  }


  int GetNumberOfFarPointsFound() //GetNumberOfSending() //GetNumberOfFarPointsFound() 
  {
    return n_send; // CD.get_n_dist_points() 
  }


  int GetNumberOfLocalPointsFound() //GetNumberOfReceiving() 
  {
    return n_recv; // CD.get_n_interior()  
  }


  vector<int> GetTouchedElements()
  {
    return dist_locations_i;     
  }


  vector<double> GetFarPointsFound(std::string format, int print=1) 
  {
    int dim = CELL_DIM; 
 
    std::vector<double> result; 
    if(format=="ple") result = dist_coords_j;  // n_send * CELL_DIM   
    else 
    if(format=="vtk") 
    {
      int VTK_DIM = 3; 
      result = std::vector<double>( n_send * VTK_DIM);  
      for(int i=0, k=0; i<n_send; i++)
        for(int j=0; j<CELL_DIM; j++, k++) result[i*VTK_DIM + j] = dist_coords_j[k];
      dim = VTK_DIM;  
    }  

    if(print)
    {
      for(int i=0, k=0; i<n_send; i++)
      {
        std::cout<<" "<< i <<") [";
        for(int j=0; j<dim; j++, k++) std::cout<<" "<< result[k];
        std::cout<<" ] \n";
      }
    }

    std::cout<<"\t[GetFarPointsFound] nPts:"<< result.size() / dim <<" \n";
    return result; 
  } 


  void Exchange(std::vector<double> var_ij, std::vector<double>& array)
  {
    assert(vertex_coords_i.size()/CELL_DIM == array.size()); 
 
    std::vector<double> var_ji = Exchange(var_ij); 
  //for(auto x:var_ji) std::cout<< x <<" ";

    for(int i=0; i<array.size();  i++) array[i] = 0; 
    for(int i=0; i<var_ji.size(); i++) array[interior_list_j[i] - 1] = var_ji[i]; 
  //for(auto x:array) std::cout<< x <<" ";
  }


  std::vector<double> Exchange(std::vector<double> var_ij)  
  {
    std::vector<double> var_ji(n_recv, HUGE_VAL);
/* 
    std::cout<<" var_ij:"<< var_ij.size()
             <<" n_send:"<< n_send  
             <<" var_ji:"<< var_ji.size()  
             <<" n_recv:"<< n_recv 
             <<"\n";
*/
    assert(n_send == var_ij.size());   
  //for(auto x:var_ij) std::cout<< x <<" ";

    std::cout<<"\t[Exchange] ..."; 
    CD.__locator_exchange_double_scalar__(&var_ij[0], &var_ji[0], 1);
    std::cout<<" ok!\n"; 

  //for(auto x:var_ji) std::cout<< x <<" ";
    return var_ji; 
  }


  protected :
  int    CELL_DIM = -1; 

  vector<int>    vertex_num_i;  
  vector<double> vertex_coords_i;   
  vector<double> shapef_j; 

  int n_send = -1;
  vector<int>    dist_locations_i; 
  vector<double>    dist_coords_j; 
  vector<double>     dist_props_j; 
 
  int n_recv = -1;
  vector<int> interior_list_j; 

}; // Locatization


//--------------------------------------------------------------------------||--//


/* 
//--------------------------------------------------------------------------||--//
class Nek5KCoupling : public Locatization 
{
  public :
 ~Nek5KCoupling(){};

  Nek5KCoupling(std::string type, std::string namei, std::string namej, int dim) : Locatization(type, namei, namej, dim)
  {

  }

}; // Nek5KCoupling
*/

//--------------------------------------------------------------------------||--//
class HemelbCoupling : public Coupling  
{
  public :
 ~HemelbCoupling(){};


  HemelbCoupling(std::string type, std::string namei, std::string namej, int dim) : Coupling(type, namei, namej)
  {

  }


  //---------------------------------------------------------------| main |---//
  int Test(int print) 
  {
    int lrank = -1, grank = -1, gsize = -1, froot = -1, lroot = -1;
    MPI_Comm_rank( local_comm, &lrank);
    MPI_Comm_rank(global_comm, &grank);
    MPI_Comm_size(global_comm, &gsize);

    lroot = (lrank)?(-69):(grank);  // if(lrank==0) then lroot==grank   
    sendrecv_int(&lroot, 1, &froot, 1, 0);

    MPI_Bcast(&lroot, 1, MPI_INT, 0, local_comm);  
    if(print) std::cout<<"\t["<< NAMEi <<"] local_root:"<< lroot <<"/"<< gsize <<" far_root:"<< froot<<"/"<< gsize << " \n"; 
    return froot;  
  }


  double Exchange(double send, int print)  
  {
    int nCoupledLetsHere  = 1;
    int nCoupledLetsThere = 0;
    sendrecv_int(&nCoupledLetsHere, 1, &nCoupledLetsThere, 1, 0);

    std::vector<double> coupledField(nCoupledLetsHere);
    std::vector<double> coupledFieldThere(nCoupledLetsThere);

    coupledField[0] = send;
    sendrecv_double(&coupledField[0], coupledField.size(), &coupledFieldThere[0], coupledFieldThere.size(), print);
    //for(auto r: coupledFieldThere) std::cout<<" "<< r;

    return coupledFieldThere[0];  
  }


  //----------------------------------| SimulationMaster::RunSimulation() |---//  
  int SimulationSetTotalTimeSteps(int nsteps) 
  {
    int recv = -1;  
    sendrecv_int(&nsteps, 1, &recv, 1, 0);
    if( LRoot() ) std::cout<<"\t["<< IAM <<"] nSteps:"<< recv <<" \n"; 
    return recv; 
  }


  //-----------------------------------------| LocalPropertyOutput::Write |---//


};





//--------------------------------------------------------------------------||--//
#endif
/*
  SEE : 
    /home/bsc21/bsc21704/z2019_2/REPOSITORY/ALYA/ALYA_2019JUL12/Thirdparties/libple/PLEPP/Wrappers/Cpp

*/ 
