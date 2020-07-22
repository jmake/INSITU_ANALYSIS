#include <vtkIdTypeArray.h>
#include <vtkUnsignedCharArray.h>
#include <vtkPoints.h>
#include <vtkPointData.h>
#include <vtkUnstructuredGrid.h>
#include <vtkCellArray.h>
#include <vtkDoubleArray.h>
#include <vtkCPDataDescription.h>
#include <vtkCPInputDataDescription.h>
#include <vtkCPProcessor.h>
#include <vtkMPIController.h>
#include <vtkMPI.h>

#include <mpi.h>
#include <vector>
#include <numeric>      // std::iota
#include <vtkMPIController.h> 
#include "insitu.hpp"  

typedef uint64_t UINT;

void process(vtkMultiProcessController*, void*);

InSituCpp* inSituCpp = NULL; 

extern "C"
{ 
  void buildvtkgrid_(int     *numberOfPoints,
                     double  *pointsData,
                     int     *numberOfCells,
                     int     *cellsTypes,
                     int     *cellsData,
                     int     *cellsDataMax,  
                     int     *dim, 
		     int     *debug   
		     )
  {
    assert(inSituCpp); // <- It is assumed that 'nekcatalystinitialize' has already been called!! 

    int print = debug[0];

    if( !inSituCpp->GetMesh() )  
    {
      // Points   
      vtkPoints *points = vtkPoints::New();
      points->SetNumberOfPoints(numberOfPoints[0]);
      if(dim[0]==2)
      {
        for(int i=0, k=0; i<numberOfPoints[0]; i++, k+=2) points->SetPoint(i, pointsData[k+0], pointsData[k+1], 0.0            );
      }
      if(dim[0]==3)
      {
        for(int i=0, k=0; i<numberOfPoints[0]; i++, k+=3) points->SetPoint(i, pointsData[k+0], pointsData[k+1], pointsData[k+2]);
      }

      vtkUnstructuredGrid *VTKGrid = vtkUnstructuredGrid::New();
      VTKGrid->SetPoints(points);
      points->Delete();

      // Cells 
      std::map<int,VTKCellType>  Alya2Cs;
      Alya2Cs[10] = VTK_TRIANGLE;   //TRI03 = 10
      Alya2Cs[12] = VTK_QUAD;       //QUA04 = 12
      Alya2Cs[30] = VTK_TETRA;      //TET04 = 30
      Alya2Cs[32] = VTK_PYRAMID;    //PYR05 = 32
      Alya2Cs[34] = VTK_WEDGE;      //PEN06 = 34
      Alya2Cs[37] = VTK_HEXAHEDRON; //HEX08 = 37

      int vtk_size_max = cellsDataMax[0];
      int nCells       = numberOfCells[0]; 
      for(int i=0, offset=0; i<nCells; i++)
      {
        int alya_type  = cellsTypes[i]; 
        int  vtk_type  = Alya2Cs[alya_type]; 
        int  vtk_size  = -1; 
	switch(vtk_type) 
        {
          case VTK_TETRA : 
            vtk_size = 4;  
            break;
          case VTK_PYRAMID :
            vtk_size = 5;
            break;
          case VTK_WEDGE :
            vtk_size = 6;
            break;
	  case VTK_HEXAHEDRON : 
            vtk_size = 8;
            break;
          default : 
	    exit(0); 
	}
	assert(vtk_size <= vtk_size_max); 

        std::vector<int> idx(vtk_size_max);  
        std::iota(idx.begin(), idx.end(), offset);

        std::vector<vtkIdType> Cells;  
        for(auto j : idx) Cells.push_back( cellsData[j] );

        if( print )
        {
	  std::cout<< i <<") "
		   <<" alya:"<< alya_type 
	           <<" vtk:" << vtk_type  
		   <<" size:"<< vtk_size
                   <<" max_size:"<< vtk_size_max
		   <<" offset:"<< offset 
		   <<" V:[ ";   
          for(auto c :  Cells) std::cout<< c <<" ";
          std::cout<<"] \n"; 
	}

        for(int j=0; j<Cells.size(); j++) Cells[j]-=1;  
        VTKGrid->InsertNextCell(vtk_type, vtk_size, Cells.data());   
        offset += vtk_size_max;   
      }

      // Mesh 
      inSituCpp->SetMesh(VTKGrid, dim[0]);
    } 

  } 


  void add_scalar_field_(double *data, char *name, int *ndata) 
  {
    if(inSituCpp) inSituCpp->UseArray( std::string(name), ndata[0], 1, data);
  }


  void add_vector_field_(double *xdata, double *ydata, double *zdata, int *dim, char *name, int *ndata) 
  {
    exit(1); 
  }


  void nekcatalystinitialize_(int* fcomm, int* frequency)
  {
    MPI_Comm ccomm = MPI_Comm_f2c(fcomm[0]); //{MPI_COMM_WORLD}; 

    std::string fileName = "output";
    int outputFrequency  = frequency[0]; 

    inSituCpp = new InSituCpp(); 
    inSituCpp->Initialize(ccomm, outputFrequency, process);
  }


  void nekcatalystfinalize_( ) 
  {
    if(inSituCpp) inSituCpp->Finalize();  
  }


  void nekcatalystcoprocess_(double* time, unsigned int *timeStep, int *save) 
  {
    if(inSituCpp) inSituCpp->Coprocess(time[0], timeStep[0], save[0]); 
  }

} // extern 


///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
void process(vtkMultiProcessController* controller, void* vtkNotUsed(arg))
{
  // SEE : ParallelProcessing/Generic/Cxx/TaskParallelism.cxx, classvtkMultiProcessController  
  int rank = controller->GetLocalProcessId();
  int size = controller->GetNumberOfProcesses();

  int sendBuffer = 1;
  int recvBuffer = 0;
  int length     = 1;
  controller->AllReduce(&sendBuffer, &recvBuffer, length, vtkCommunicator::SUM_OP);
  assert(size==recvBuffer);

  controller->Barrier();
  //std::cout<<"[process] rank/size:"<< rank <<"."<< recvBuffer <<"\n";
  //assert( !rank ); 
}

 
