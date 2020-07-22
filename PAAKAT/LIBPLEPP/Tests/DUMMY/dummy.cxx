/*
  2019NOV14. BSC, BARCELONA, SPAIN 
  Migue Zavala 
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
#include <cctype> // tolower 

#include "plepp.hpp"  
#include "insitu_py.hpp"

std::string GetLowers(std::string in)
{
  std::string out(in); 
  std::transform(out.begin(), out.end(), out.begin(), [](char c){ return std::tolower(c); });
  return out; 
}


int main(int argc, char** argv)
{
  //---------------------------------------------------------------------||---//
  assert(argc==2);
  std::string IAM = argv[1];


  //--------------------------------------------------------------| PLEPP |---//
  MPI_Init(NULL, NULL);

  HemelbCoupling *hemelbCoupling = NULL; 
  if(IAM.find("ALYA"  ) != std::string::npos) hemelbCoupling = new HemelbCoupling("PLEPP","ALYA",  "HEMELB", 0); 
  if(IAM.find("HEMELB") != std::string::npos) hemelbCoupling = new HemelbCoupling("PLEPP","HEMELB","ALYA"  , 0);
  assert(hemelbCoupling);

  MPI_Comm global_comm = MPI_COMM_WORLD; 
  MPI_Comm  local_comm = hemelbCoupling->Init( global_comm );  

  int far_root = hemelbCoupling->Test( hemelbCoupling->LRoot() );
  assert( (far_root==0)||(far_root==2) ); 

#ifdef WITHALYA  
#else     
  //---------------------------------------------------------------------||---//
  ofstream myfile;
  if( hemelbCoupling->LRoot() )  
  {
    myfile.open("iam_"+GetLowers(IAM) + ".dat");
    myfile<<"## itime recv";
    myfile<< scientific << showpoint;
    myfile<< setprecision( 10 );
    myfile<<" \n";
  }
  //hemelbCoupling->GBarrier();


  //----------------------------------| SimulationMaster::RunSimulation() |---//  
  int dummy     =  3;    
  int timeSteps = hemelbCoupling->SimulationSetTotalTimeSteps( dummy ); 


  //-----------------------------------------| LocalPropertyOutput::Write |---//
  for(int iSteps=0; iSteps<timeSteps; iSteps++)
  {
    std::vector<double> recv;  
    std::vector<double> send = {-1.0,-2.0}; 
    for(auto s: send) recv.push_back( hemelbCoupling->Exchange(s,0) ); 

    if( hemelbCoupling->LRoot() ) 
    { 
      myfile<<" "<< iSteps;
      for(auto r: recv) myfile<<" "<< r; 
      myfile<<" \n"; 
    }
  }

 
  //---------------------------------------------------------------------||---//
  if( hemelbCoupling->LRoot() )
  {
    std::cout<<"\t["<< IAM <<"] OK!! \n"; 
    myfile.close();
  }
#endif 

  MPI_Finalize();
  return 0;
}
//=======================================================================||===//
//=======================================================================||===//
