/*  
  2019DIC02. KTH, STOCKHOLM, SWEDEN. Migue Zavala 
  FROM : 
  /afs/pdc.kth.se/home/m/miguelza/z2020_1/REPOSITORY/NEK5K/MPMD01_1/FromMAC/NEK5_PLEPP/NEK5_PLEPP/Example02
*/
#include <iostream>
#include <vector>
#include <map>
#include <mpi.h>
#include <cmath>
#include <assert.h>
#include <algorithm>
#include <sstream>
#include <iomanip>
#include <fstream>
#include <limits>

#include "plepp.hpp"

std::string IAM;
std::string JAM;

int main(int argc, char** argv)
{
  //---------------------------------------------------------------------||---//
  assert(argc==2);
  IAM = argv[1];

  //--------------------------------------------------------------| PLEPP |---//
  MPI_Init(NULL, NULL);
  MPI_Comm  lcomm = MPI_COMM_NULL;
  MPI_Comm  gcomm = MPI_COMM_WORLD;
// /*
  Nek5KCoupling *coupling = NULL;
  if(IAM=="EXPLORER") coupling = new Nek5KCoupling("INSITU", IAM, "NEK5K"   , 2);
  if(IAM=="NEK5K")    coupling = new Nek5KCoupling("INSITU", IAM, "EXPLORER", 2);
  assert(coupling); 
  lcomm = coupling->Init(gcomm);

  // X.0 
  //coupling->sendrecv_string(IAM,1); 
    int lrank = -1, grank = -1, groot = -1, lroot = -1;
    MPI_Comm_rank(lcomm, &lrank);
    MPI_Comm_rank(gcomm, &grank);

    lroot = (grank-lrank)?(grank):(-69);
    coupling->sendrecv_int(&lroot, 1, &groot, 1, 1); 
//  */

//MPI_Barrier(gcomm);


//  coupling->sendrecv_double(double* send, int n_send, double* recv, int n_recv, int print=0)
/* 
  std::string csend(IAM); 
  coupling->sendrecv_string(csend,1); 
*/

  //---------------------------------------------------------------------||---//
  //if(!local_rank) std::cout<<" \n";
  MPI_Finalize();
  return 0;
}
