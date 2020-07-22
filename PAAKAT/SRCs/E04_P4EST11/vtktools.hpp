#ifndef VTK_TOOLS_2_H  
#define VTK_TOOLS_2_H 

#include <vtkAppendFilter.h>
#include <vtkVersion.h>
#include <vtkFloatArray.h>

#include <vtkBox.h>
#include <vtkExtractGeometry.h>
#include <vtkImplicitFunction.h>
#include <vtkPlane.h>
#include <vtkCutter.h>
#include <vtkContourFilter.h>

#include "vtktools_basics.hpp"

//--------------------------------------------------------------------------||--//
//--------------------------------------------------------------------------||--//
vtkDataObject* ExtractGeometry(vtkDataObject *Obj, std::vector<double> bbox)
{
  vtkBox *b = vtkBox::New();
  b->SetBounds( bbox.data() ); // bbox = [xmin,xmax, ymin,ymax, zmin,zmax]  

  vtkExtractGeometry *eg = vtkExtractGeometry::New();
  eg->SetInputData( Obj );
  eg->SetImplicitFunction(b);
  eg->ExtractInsideOn();
  eg->ExtractBoundaryCellsOn();
  eg->Update();
  return eg->GetOutputDataObject(0);
}


//--------------------------------------------------------------------------||--//
vtkImplicitFunction* GetFuntionPlane(std::vector<double> Orig, std::vector<double> Normal)
{
  vtkPlane *plane = vtkPlane::New();
  plane->SetOrigin(&Orig[0]);
  plane->SetNormal(&Normal[0]);

  return plane;
}


//--------------------------------------------------------------------------||--//
vtkDataObject* Cutter(vtkDataObject* Obj, vtkImplicitFunction* function)
{
  vtkCutter *cutter = vtkCutter::New();
  cutter->SetCutFunction( function );
  cutter->SetInputData( Obj );
  cutter->Update();

  return cutter->GetOutputDataObject(0);
}


//--------------------------------------------------------------------------||--//
vtkDataObject* 
GetContour(vtkDataObject *input, std::string key, double threshold)
{
//vtkSmartPointer<vtkContourFilter> contourFilter = vtkSmartPointer<vtkContourFilter>::New(); // SEGMENTATION!!
  vtkContourFilter* contourFilter = vtkContourFilter::New();
  contourFilter->SetInputData( input );
  contourFilter->SetValue(0,threshold);  // IsoSurface(0) = threshold
  contourFilter->SetInputArrayToProcess(0, 0, 0, 0, key.c_str() );

  contourFilter->GenerateTrianglesOn();
  contourFilter->ComputeScalarsOff();
  contourFilter->ComputeNormalsOff();
  contourFilter->Update();

  return contourFilter->GetOutputDataObject(0);
/*
  vtkDataObject *output = contourFilter->GetOutputDataObject(0); assert( output->IsA("vtkPolyData") );
  vtkPolyData   *poly   = vtkPolyData::SafeDownCast(output); assert(poly);

  poly->GetCellData()->RemoveArray("vtkOriginalCellIds");
  poly->GetFieldData()->RemoveArray("__CatalystChannel__");
  std::cout<<" GetNumberOfPoints "<< poly->GetNumberOfPoints() <<" \n";
 
  return poly;
*/
}


//--------------------------------------------------------------------------||--//
template<class intT, class doubleT>
class BuildVTU : public BuildGrid<intT, doubleT, vtkUnstructuredGrid>
{
  public :
  BuildVTU(int _d) : BuildGrid<intT, doubleT, vtkUnstructuredGrid>(_d) {};


  void ApplyFilter(std::string fname)
  {
    Warning("vtk.log"); 

    // 0.0.  
    CreateVtkMPIController *MPIController = new CreateVtkMPIController();
    vtkMPIController *contr = MPIController->GetController();

    // 1.1. 
    vtkDataObject *obj = this->VTKGrid;  
    assert(obj);  
    //PWriter2(obj, fname+"_entire", contr);

    // 1.2. 
    //vtkDataObject *extracted = ExtractGeometry(obj, {-1.0,3.0, -1.0,1.0, 0.0,1.0}); assert(extracted); // [1]
    //PWriter2(extracted, fname+"_extracted", contr);

    // 1.3. 
    vtkDataObject *cutter1 = Cutter(obj, GetFuntionPlane({2.0,0.0,0.01},{0.0,0.0,1.0}) ); assert(cutter1); 
    PWriter2(cutter1, fname+"_xycutter", contr);

    // 1.4. 
    vtkDataObject *cutter2 = Cutter(obj, GetFuntionPlane({2.0,0.0,0.01},{1.0,0.0,0.0}) ); assert(cutter2);
    PWriter2(cutter2, fname+"_zycutter", contr);

    // 1.5. 
    vtkDataObject *cutter3 = Cutter(obj, GetFuntionPlane({2.0,-0.0435779,0.01},{0.0,1.0,0.0}) ); assert(cutter3); // [0.498097, -0.0435779, 0.01]
    PWriter2(cutter3, fname+"_xzcutter", contr);

    // 1.5. 
    vtkDataObject *contour2 = GetContour(obj, "VELOX", 1.0); assert(contour2);  
    PWriter2(contour2, fname+"_contour_vx", contr);

    // 1.5. 
    vtkDataObject *contour1 = GetContour(obj, "VELOY", 1.0); assert(contour1); 
    PWriter2(contour1, fname+"_contour_vy", contr);

    // 
    // [1] bbox = [xmin,xmax, ymin,ymax, zmin,zmax]
    //
  }

};



//--------------------------------------------------------------------------||--//
//--------------------------------------------------------------------------||--//
#endif
/*
  J. MIGUEL ZAVALA AKE. 2019AUG18. STOCKHOLM, SWEDEN.

  2020JUN15 : 
  /home/bsc21/bsc21704/z2020_1/REPOSITORY/PAAKAT/TARs/2019AUG09_LIGHT/TESTS/PV560/FORTRAN/E02_BASE

  2020JUN25 :
  /afs/pdc.kth.se/home/m/miguelza/z2020_1/REPOSITORY/ALYAs/ALYA_2020ABR17/alya/Thirdparties/INSITU/PAAKAT/E04_P4EST11

*/ 
