#include "plepp.hpp"

//--------------------------------------------------------------------------||--//
const char* String_f2c(const char* licName, int nameLen)
{
  char* name;
  name = (char*)malloc(sizeof(char)*(nameLen));
  for(int i = 0; i<nameLen; i++) name[i] = licName[i];
  name[nameLen] = '\0';

  return name ;
}


static CommDom*             CD = NULL;
static HemelbCoupling *coupling = NULL;

//=======================================================================||===//
extern "C"
{

  void
  plepp_create_(
                const char* ftype,  int* ntype,
                const char* fnamei, int* nnamei,
                const char* fnamej, int* nnamej,
                int*        fdim, 
                int*        fcomm  
	       )  
  {
    //std::cout<<"\t[plepp_create] \n"; 

    // X.0 
    std::string namej( String_f2c(fnamej, nnamej[0]) ); // El order de namej, namei da segmentation!!
    std::string namei( String_f2c(fnamei, nnamei[0]) );
    std::string  type( String_f2c( ftype,  ntype[0]) );
    int dim = fdim[0];
 
    // X.0 
    coupling = new HemelbCoupling(type, namei, namej, dim);
    assert(coupling);

    MPI_Comm gcomm = MPI_Comm_f2c( fcomm[0] );
    MPI_Comm lcomm = coupling->Init( gcomm );

    fcomm[0] = MPI_Comm_c2f(lcomm);
  }    


  void plepp_test_() //int* print)
  {
//    coupling->Test( print[0] );
    int far_root = coupling->Test( coupling->LRoot() );
//    assert( (far_root==0)||(far_root==2) );
  }


  void plepp_exchange_(double* send, double* recv, int* print)   
  {
    recv[0] = coupling->Exchange(send[0],0);  
  }


  void plepp_sendrecv_real_(double* send, int* n_send, double* recv, int* n_recv, int* print)
  {
    coupling->sendrecv_double(send, n_send[0], recv, n_recv[0], print[0]);  
  }


  void plepp_sendrecv_int_(int* send, int* n_send, int* recv, int* n_recv, int* print)
  {
    coupling->sendrecv_int(send, n_send[0], recv, n_recv[0], print[0]);
  }


  void plepp_sendrecv_char_(char* send, int* n_send, char* recv, int* n_recv, int* print)  
  {
    std::string csend( String_f2c(send, n_send[0])); 
    std::string crecv( coupling->sendrecv_string(csend, print[0]) );
  }



}


//=======================================================================||===//
static AnalyseArgvs*  ptr_class1 = NULL;

extern "C"
{

  void
  plepp_init_argvs_() 
  {
    ptr_class1 = new AnalyseArgvs();
  }


  void
  plepp_set_argvs_(const char* ftype, int* ntype)
  {
    std::string  ctype( String_f2c(ftype, ntype[0]) );
    ptr_class1->__set_argvs__(ctype);
  }


  void
  plepp_get_argvs_(char* ftype, int* print)
  {
    std::string  ctype = ptr_class1->__get_argvs__();

    int len = ctype.size(); 
    char* log_name = (char *)malloc(sizeof(char)*(len+1));
    strncpy(log_name, ctype.c_str(), len);
    log_name[len] = '\0';

    strcpy(ftype, log_name);   
    if(print[0]==1) cout<<"\t[plepp_get_argvs] \'"<< ftype <<"\' \n\n";
  }


  void
  plepp_analyse_argvs_(const char* ftype, int* ntype, int* print)
  {
    //ptr_class1->__print_argvs__(); 
    std::string  ctype( String_f2c(ftype, ntype[0]) );
    if(ctype.size() > 0) ptr_class1->__analyse_argvs__(ctype,print[0]);
  }


  void
  plepp_split_(int* gcommf, int* lcommf, int* print)   
  {
    MPI_Comm gcomm = MPI_Comm_f2c( gcommf[0] );
    MPI_Comm lcomm = MPI_COMM_NULL;

    lcomm = ptr_class1->SimpleSplit(gcomm,print[0]); 
    assert(lcomm != MPI_COMM_NULL);  

    int nranks = -1;
    MPI_Comm_size(lcomm, &nranks);

    lcommf[0] = MPI_Comm_c2f(lcomm);
  }



}


 /* 
//=======================================================================||===//
static InputParser *  ptr_class2 = NULL; 

extern "C"
{
  void
  plepp_init_parser_(int argc, char **argv)
  {
    ptr_class2 = new InputParser();
  }



} 
*/



//=======================================================================||===//
