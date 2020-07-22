#include <vtkPVTrivialProducer.h>

#include <vtkImplicitFunction.h>
#include <vtkPlane.h>
#include <vtkCutter.h>
#include <vtkFeatureEdges.h>
//#include <vtkXMLMultiBlockDataReader.h>

#include "vtktools.hpp"

/*
vtkDataObject* VtmReader(std::string fname)
{
  vtkXMLMultiBlockDataReader *reader = vtkXMLMultiBlockDataReader::New();
  reader->SetFileName( fname.c_str() );
  reader->Update(); // reader->Print(std::cout);

  std::cout<<"\t [Reader] '"<< fname <<"'!! \n";
  return reader->GetOutputDataObject(0);
}
*/

vtkDataObject* ExtractSurface01_1(vtkDataObject* Obj)  
{
  vtkFeatureEdges *featureEdges = vtkFeatureEdges::New();
  featureEdges->SetInputData(Obj);
/*
  featureEdges->BoundaryEdgesOn();
  featureEdges->FeatureEdgesOff();
  featureEdges->ManifoldEdgesOff();
  featureEdges->NonManifoldEdgesOff(); 
*/ 
  featureEdges->Update();

  return featureEdges->GetOutputDataObject(0); 
}


vtkImplicitFunction* GetFuntionPlane(std::vector<double> Orig, std::vector<double> Normal)   
{ 
  vtkPlane *plane = vtkPlane::New(); 
  plane->SetOrigin(&Orig[0]);
  plane->SetNormal(&Normal[0]);

  return plane; 
}


vtkDataObject* Cutter(vtkDataObject* Obj, vtkImplicitFunction* function) 
{
  vtkCutter *cutter = vtkCutter::New(); 
  cutter->SetCutFunction( function );
  cutter->SetInputData( Obj ); 
  cutter->Update();

  return cutter->GetOutputDataObject(0);
}


void Simplest(vtkDataObject *obj, int iTimeStep) // 2020MAR23  
{
  assert(obj);
  //GetVtkDataObjectMaxMin( obj ); 

  // 0.0. 
  std::ostringstream o; o << std::setfill('0') << std::setw(6) << iTimeStep;
  std::string name = "simplest_" + o.str();

  // 1.0. 
//vtkDataObject *cutter = Cutter(obj, GetFuntionPlane({0.0,0.0,0.1},{0.0,0.0,1.0}) ); 
  vtkDataObject *cutter = ExtractGeometry(obj, {0.00-0.01,1.00+0.01, -0.03-0.003,0.10+0.01, 0.00-0.01,0.20+0.02});

  // 2.0. 
  vtkNew<vtkPVTrivialProducer> producer;
  producer->SetOutput( cutter );// ( obj );

  vtkDataObject *output = producer->GetOutputDataObject(0); 
  std::string      type = producer->GetOutputDataObject(0)->GetClassName(); 
  PVtkObjectWriter(output, name, type); 
//VtkObjectWriter(obj, "simplest_" + o, obj->GetClassName()); 

  // 
  // (1.00, 0.00)
  // (0.00, 0.00)
  // (0.30, 0.10)
  // (0.13,-0.03)
  //
  //ExtractGeometry();  


}
