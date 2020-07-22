#include <iostream> // cin, cout, endl, cerr
#include <vector>   // vector
#include <map>
#include <mpi.h>
#include <unistd.h> // getopt 

#include <fstream>
#include <numeric>      // std::iota
#include <algorithm>    // std::min_element, std::max_element, std::sort 
#include <math.h>       /* atan2 */
#include <functional>

#include "simplest.hpp"  // VtkInterpolation  


int main(int argc, char* argv[])
{
  //--------------------------------------------------------------------||---//
  assert(argc==2);  
  std::string fname( argv[1] );   

  MPI_Init(NULL, NULL);
  MPI_Comm  local_comm = MPI_COMM_WORLD; 

  //vtkDataObject *vtu = Reader(fname, "vtu"); 
  vtkDataObject *vtu = VtmReader(fname);
  Simplest(vtu, 0); 

  vtkDataObject *obj = ExtractSurface01_1(vtu); 

  std::string name   = "test";  
  std::string type   = obj->GetClassName();
//  PVtkObjectWriter(obj, name, type); 

  //--------------------------------------------------------------------||---//
  MPI_Finalize();
  return 0;
}
/*
*/
