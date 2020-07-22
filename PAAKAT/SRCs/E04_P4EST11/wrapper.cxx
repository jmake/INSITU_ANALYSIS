#include "vtkCPDataDescription.h"
#include "vtkCPInputDataDescription.h"
#include "vtkCPProcessor.h"
#include "vtkSmartPointer.h"
#include "vtkDoubleArray.h"
#include "vtkPointData.h"
#include "vtkImageData.h"

#include "vtkCellData.h"
#include "vtkCellType.h"
#include "vtkCPAdaptorAPI.h"
#include "vtkFloatArray.h"
#include "vtkNew.h"
#include "vtkPoints.h"
#include "vtkUnstructuredGrid.h"

#include <vtkPVTrivialProducer.h>
#include <vtkCompleteArrays.h>
#include <vtkXMLPUnstructuredGridWriter.h>
#include <vtkXMLUnstructuredGridWriter.h>

#include <map>
#include <mpi.h>

//#include "vtktools.hpp"
#include "p4est_paakat.hpp"


extern "C"
{
  void 
  buildvtkgrid_(int     *numberOfPoints, 
                double  *pointsData,
                int     *numberOfCells,
                int     *cellsTypes,
                int     *cellsData,
                int     *dim )
  {   
    vtkUnstructuredGrid *VTKGrid = vtkUnstructuredGrid::New(); 
    vtkPoints *points = vtkPoints::New();
    vtkNew<vtkPVTrivialProducer> producer;
    vtkNew<vtkCompleteArrays> completeArrays;
  }


  void 
  paakat_test_(int *print)
  {
//    PrintVtkVersion( print[0] );
  }


  void 
  create_array_(double *data, char *name, int *ndata)
  {
  }


  void 
  use_array_(double *data, char *name, int *ndata)
  {
  }


  void 
  paakat_initialize_(int* fcomm)
  {
  }


  void 
  paakat_finalize_( )
  {
  }


  void 
  paakat_coprocess_(double* time, unsigned int *timeStep)
  {
  }


  void 
  paakat_fp4est_(p4est_t* p4est) 
  {
    assert(p4est); 
    p4est_paakat_write_file(p4est, NULL, "paakat");

    if(!p4est->mpirank) printf("\t[paakat_fp4est]");
  }  


  void
  paakat_fp4est(p4est_t* p4est)
  {
    assert(p4est);
    p4est_paakat_write_file(p4est, NULL, "paakat");

    if(!p4est->mpirank) printf("\t[paakat_fp4est]");
  }


} // extern  
