//#include "simplest.hpp" 
#include "pipeline.hpp"

#include <vtkCPDataDescription.h>
#include <vtkCPInputDataDescription.h>
#include <vtkCommunicator.h>
#include <vtkCompleteArrays.h>
#include <vtkDataArray.h>
#include <vtkMultiProcessController.h>
#include <vtkNew.h>
#include <vtkObjectFactory.h>
#include <vtkPVArrayCalculator.h>
#include <vtkPVTrivialProducer.h>
#include <vtkPointData.h>
#include <vtkSMProxyManager.h>
#include <vtkThreshold.h>
#include <vtkUnstructuredGrid.h>
#include <vtkXMLPUnstructuredGridWriter.h>

#include <sstream>
#include <string>


vtkStandardNewMacro(vtkCPVTKPipeline);

//----------------------------------------------------------------------------
vtkCPVTKPipeline::vtkCPVTKPipeline()
{
  this->OutputFrequency = 0;
}

//----------------------------------------------------------------------------
vtkCPVTKPipeline::~vtkCPVTKPipeline()
{
}


//----------------------------------------------------------------------------
int vtkCPVTKPipeline::Finalize()
{
//  if(!rank) myfile.close();
  return 1; 
}


//----------------------------------------------------------------------------
void vtkCPVTKPipeline::Initialize(int outputFrequency, std::string& fileName)
{
  this->OutputFrequency = outputFrequency;
  this->FileName = fileName;

/*
  vtkMPIController *contr = vtkMPIController::New();
//contr->Initialize(NULL, NULL, 1); // initializedExternally == 1;
//int nranks = contr->GetNumberOfProcesses();

  rank = contr->GetLocalProcessId();
//  if(!rank) myfile.open("./simplest.dat"); 
*/ 

}

//----------------------------------------------------------------------------
int vtkCPVTKPipeline::RequestDataDescription(vtkCPDataDescription* dataDescription)
{
  if (!dataDescription)
  {
    vtkWarningMacro("dataDescription is NULL.");
    return 0;
  }

  if (this->FileName.empty())
  {
    vtkWarningMacro("No output file name given to output results to.");
    return 0;
  }

  if (dataDescription->GetForceOutput() == true ||
    (this->OutputFrequency != 0 && dataDescription->GetTimeStep() % this->OutputFrequency == 0))
  {
    dataDescription->GetInputDescriptionByName("input")->AllFieldsOn();
    dataDescription->GetInputDescriptionByName("input")->GenerateMeshOn();
    return 1;
  }
  return 0;
}

//----------------------------------------------------------------------------
int vtkCPVTKPipeline::CoProcess(vtkCPDataDescription* dataDescription)
{
  if (!dataDescription)
  {
    vtkWarningMacro("DataDescription is NULL");
    return 0;
  }
  vtkUnstructuredGrid* grid = vtkUnstructuredGrid::SafeDownCast(
    dataDescription->GetInputDescriptionByName("input")->GetGrid());
  if (grid == NULL)
  {
    vtkWarningMacro("DataDescription is missing input unstructured grid.");
    return 0;
  }
  if (this->RequestDataDescription(dataDescription) == 0)
  {
    return 1;
  }
/*
  // Simplest 
  int    step   = dataDescription->GetTimeStep();  
  double time   = dataDescription->GetTime();   
  double result = Simplest( grid, step );  

  if(!rank)
  { 
    std::cout<<"\t[Simplest] "<< step <<" "<< time <<" "<< result <<" \n";

    myfile.open("./simplest.dat",  std::ios_base::app);
    myfile<< step <<" "<< time <<" "<< result <<" \n";  
    myfile.close();
  } 
*/ 
  return 1;
}

//----------------------------------------------------------------------------
void vtkCPVTKPipeline::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os, indent);
  os << indent << "OutputFrequency: " << this->OutputFrequency << "\n";
  os << indent << "FileName: " << this->FileName << "\n";
}

