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
#include <vtkMPIController.h> 
#include "insitu_py.hpp"  

typedef uint64_t UINT;

//void TestProcessController(vtkMultiProcessController*, void*);

InSituCpp* inSituCpp = NULL; 

extern "C"
{ 
  void paakat_buildvtkgrid_(const double *x, const double *y, const double *z,
                     const int *lx1, const int *ly1, const int *lz1, 
                     const int *lelt, const int *dim)
  {
    assert(inSituCpp); // <- It is assumed that 'paakat_initialize' has already been called!! 

if( !inSituCpp->GetMesh() )  
{ 
    const vtkIdType points_per_el = (vtkIdType) *lx1 * (vtkIdType) * ly1 * (vtkIdType) *lz1;

    vtkUnstructuredGrid* VTKGrid = vtkUnstructuredGrid::New();
    vtkPoints *points = vtkPoints::New();
  
    for (UINT i = 0; i < points_per_el * (*lelt); i++) points->InsertNextPoint((float) x[i], (float) y[i], (float) z[i]);
    VTKGrid->SetPoints(points);  
    points->Delete();

    const UINT h_p_lelt = (((*lx1) - 1) * ((*ly1) - 1) * ((*dim) == 3 ? ((*lz1 - 1)) : 1));
    const vtkIdType h_lelt = h_p_lelt * (*lelt);  
    const vtkIdType total = ((*dim) == 3 ? (9 * h_lelt) : (5 * h_lelt));
  
    vtkIdTypeArray *id_list = vtkIdTypeArray::New();
    id_list->SetNumberOfValues(total);
    vtkIdType *idp  = id_list->GetPointer(0);

    vtkUnsignedCharArray *cell_type = vtkUnsignedCharArray::New();
    cell_type->SetNumberOfValues(h_lelt);
    unsigned char *ctp = cell_type->GetPointer(0);

    vtkIdTypeArray *cell_location = vtkIdTypeArray::New();
    cell_location->SetNumberOfValues(h_lelt);
    vtkIdType *clp = cell_location->GetPointer(0);

    vtkIdType hex_idx = 0;
    vtkIdType el_idx = 0;
    for (UINT i = 0; i < (*lelt); i++, el_idx++) 
    {
      vtkIdType points_offset = points_per_el * el_idx;
      for (UINT ii = 0; ii < ((*lx1) - 1); ii++) 
      {
        for (UINT jj = 0; jj < ((*ly1) - 1); jj++) 
	{
	  if ((*dim) == 2) 
          {
	    *idp++ = 4;
	    *idp++ = jj * (*lx1) + ii + points_offset;
	    *idp++ = jj * (*lx1) + (ii + 1) + points_offset;
	    *idp++ = (jj + 1) * (*lx1) + (ii + 1) + points_offset;
	    *idp++ = (jj + 1) * (*lx1) + (ii) + points_offset;
	    *ctp++ = VTK_QUAD;
	    *clp++ = 5 * (hex_idx++);
	  }
      	  else 
          {
 	    for (UINT kk = 0; kk < (*lz1 - 1); kk++) 
	    {
	      *idp++ = 8;
	      *idp++ = kk * ((*ly1) * (*lx1)) + jj * (*lx1) + ii + points_offset;
	      *idp++ = kk * ((*ly1) * (*lx1)) + jj * (*lx1) + (ii + 1) + points_offset;
	      *idp++ = kk * ((*ly1) * (*lx1)) + (jj + 1) * (*lx1) + (ii + 1) + points_offset;
	      *idp++ = kk * ((*ly1) * (*lx1)) + (jj + 1) * (*lx1) + (ii) + points_offset;
	      *idp++ = (kk + 1) * ((*ly1) * (*lx1)) + jj * (*lx1) + ii + points_offset;
	      *idp++ = (kk + 1) * ((*ly1) * (*lx1)) + jj * (*lx1) + (ii + 1) + points_offset;
	      *idp++ = (kk + 1) * ((*ly1) * (*lx1)) + (jj + 1) * (*lx1) + (ii + 1) + points_offset;
	      *idp++ = (kk + 1) * ((*ly1) * (*lx1)) + (jj + 1) * (*lx1) + (ii) + points_offset;
	      *ctp++ = VTK_HEXAHEDRON;
	      *clp++ = 9 * (hex_idx++);
	    }
	  }	  
        }
      }
    }

    vtkCellArray *cells = vtkCellArray::New();
    cells->SetCells(h_lelt, id_list);
    id_list->Delete();

    VTKGrid->SetCells(cell_type, cell_location, cells);
    cell_type->Delete();
    cell_location->Delete();
    cells->Delete();
/*
    if(1) 
    { 
      int nArrays = VTKGrid->GetPointData()->GetNumberOfArrays();
      int nCells  = VTKGrid->GetNumberOfCells();
      int nPts    = VTKGrid->GetNumberOfPoints();

      std::cout<<"\t[buildvtkgrid]" 
	       <<" nArrays:"<< nArrays 
	       <<" nCells:"<< nCells 
	       <<" nPts:"<< nPts 
	       <<" \n";
    }

    if(1)
    {
      std::cout<<"\t[buildvtkgrid]"
               <<" nPts: lx1 * ly1 * lz1 * lelt = "
    	       <<"("<< lx1[0] <<" x "<< ly1[0] <<" x "<< lz1[0] <<") x "<<  lelt[0]   
  	       <<" = "<< lx1[0] * ly1[0] * lz1[0] * lelt[0] 
               <<" \n";

      std::cout<<" \n\n";
    }
*/

    inSituCpp->SetMesh(VTKGrid, dim[0]);

} // inSituCpp->GetMesh() 

  } // buildvtkgrid  


  void paakat_add_scalar_(double *data, char *name, int *ndata) 
  {
    if(inSituCpp) inSituCpp->UseArray( std::string(name), ndata[0], 1, data);
  }


  void paakat_add_vector_(double *xdata, double *ydata, double *zdata, int *dim, char *name, int *ndata) 
  {
    exit(1); 
  }


  void paakat_initialize_(int* fcomm, int* frequency, char* pyfile)
  {
    MPI_Comm ccomm       = MPI_Comm_f2c(fcomm[0]); 
    std::string fname(pyfile);
    fname += ".py";  

    inSituCpp = new InSituCpp(); 
    inSituCpp->Initialize(ccomm, frequency[0], std::string(pyfile) + ".py", TestProcessController);
  }


  void paakat_finalize_( ) 
  {
    if(inSituCpp) inSituCpp->Finalize();  
  }


  void paakat_coprocess_(double* time, unsigned int *timeStep) 
  {
    if(inSituCpp) inSituCpp->Coprocess(time[0], timeStep[0]); 
  }


  void paakat_savemesh_(char *name)   
  {
    if(inSituCpp) inSituCpp->SaveMesh( std::string(name) );
  }


} // extern 

 
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
void TestProcessController(vtkMultiProcessController* controller, void* vtkNotUsed(arg))
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
  std::cout<<"[TestProcessController] size.rank:"<< recvBuffer <<"."<< rank <<"\n";
}
