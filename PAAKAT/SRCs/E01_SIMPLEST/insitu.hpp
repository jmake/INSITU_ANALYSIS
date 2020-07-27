#ifndef IN_SITU_VTK_H  
#define IN_SITU_VTK_H
#include <iostream>
#include <cstring>
#include <vector>
#include <map>
#include <mpi.h>

#include <vtkMPI.h>
#include <vtkMPIController.h>
#include <vtkUnstructuredGrid.h>
#include <vtkDoubleArray.h>
#include <vtkIdTypeArray.h>

#include <vtkPVTrivialProducer.h>
#include <vtkXMLPUnstructuredGridWriter.h>

#include "pipeline.hpp" // vtkCPVTKPipeline.cxx 
//#include "vtktools.hpp" 

#include <vtkCPProcessor.h>
#include <vtkCPDataDescription.h>
#include <vtkCPInputDataDescription.h>

#ifdef WITH_PYTHON 
#include <vtkCPPythonScriptPipeline.h> 
#endif 

// /* 
template <class vtkGridType>
void PWriter(vtkGridType* DataObject, std::string fname)
{
  auto contr = vtkSmartPointer<vtkMPIController>::New();
  contr->Initialize(NULL, NULL, 1); // initializedExternally == 1;

  int nranks = contr->GetNumberOfProcesses();
  int rank   = contr->GetLocalProcessId();

  vtkNew<vtkPVTrivialProducer> producer;
  vtkAlgorithmOutput *GetOutputPort;  
  if(DataObject)
  {
    int nArrays = DataObject->GetPointData()->GetNumberOfArrays();
    int nCells  = DataObject->GetNumberOfCells();
    int nPts    = DataObject->GetNumberOfPoints();
    producer->SetOutput(DataObject); 
    GetOutputPort = producer->GetOutputPort();
  }
  assert(GetOutputPort!=NULL);

  auto parallel_writer = vtkSmartPointer<vtkXMLPUnstructuredGridWriter>::New();
  parallel_writer->SetInputConnection( GetOutputPort ); // ( producer->GetOutputPort() ); 
  parallel_writer->SetController(contr);
  parallel_writer->SetFileName( fname.c_str() ); //("test.pvts");
  parallel_writer->SetNumberOfPieces(nranks);
  parallel_writer->SetStartPiece(rank);
  parallel_writer->SetEndPiece(rank);
  parallel_writer->SetDataModeToBinary();
  parallel_writer->Update();
  parallel_writer->Write();

  contr->Finalize(1); // if(finalizedExternally == 0) MPI_Finalize(); 
}
// */


template <class vtkGridType> 
class VtkGrid   
{
  public :
 ~VtkGrid(){}


  VtkGrid()
  {
    this->nDim    = -1 ;
    this->nCells  = -1 ;  
    this->nPts    = -1 ;  
    this->VTU     = NULL;
  }


  vtkGridType* GetMesh()  
  {
    //assert(this->VTU);
    return this->VTU;
  }


  void RemoveMesh()
  {
    if(this->VTU) 
    {
      this->VTU->Delete(); 
      this->VTU = NULL;
    }  
  }


  void SetMesh(vtkGridType* vtu, int dim)
  {  
    assert(vtu); 
    //assert(!this->VTU);

    this->VTU     = vtu;
    this->nDim    = dim ;
    this->nCells  = VTU->GetNumberOfCells();
    this->nPts    = VTU->GetNumberOfPoints();
  }


  void SaveMesh(std::string fname)
  {
     PWriter<vtkUnstructuredGrid>( this->VTU, fname );
   //PVtkObjectWriter(this->VTU, fname, this->VTU->GetClassName() ); 
  }


  vtkDoubleArray* GetArray(std::string name)
  {
    // SEE : UpdateVTKAttributes, CxxVTKPipelineExample/FEAdaptor.cxx 
    this->ArraysIT = this->Arrays.find(name);
    assert(this->ArraysIT!=this->Arrays.end());
    return this->Arrays[name];
  }


  void CreateArray(std::string name, int ndata, int ncomps, double *data)
  {
    assert(ndata==this->nPts);

    this->Arrays[name] = vtkDoubleArray::New();  
    this->Arrays[name]->SetName( name.c_str() );  
    this->Arrays[name]->SetNumberOfComponents(ncomps);  
    this->Arrays[name]->SetNumberOfValues(ndata); 
    for(int i=0,k=0; i<ndata; i++) for(int j=0; j<ncomps; j++, k++) this->Arrays[name]->SetValue(k,data[k]);   
    this->Keys.push_back(name); 
  }


  void UseArray(std::string name, int ndata, int ncomps, double *data)
  {
    assert(ndata==this->nPts);

    this->Arrays[name] = vtkDoubleArray::New();
    this->Arrays[name]->SetName(name.c_str());
    this->Arrays[name]->SetNumberOfComponents(ncomps);
    this->Arrays[name]->SetArray(data, static_cast<vtkIdType>(this->nPts*ncomps), 1);  
    this->Keys.push_back(name);  
  }


  void PrintArraysNames()
  {
    for(std::map<std::string,vtkDoubleArray*>::iterator it=this->Arrays.begin(); 
        it!=this->Arrays.end(); ++it) 
    std::cout<<" - "<< it->first <<" \n";
    std::cout<<" -  nDim:"<< this->nDim <<" \n";
    std::cout<<" -  nPts:"<< this->nPts <<" \n";
    std::cout<<" -nCells:"<< this->nCells <<" \n";
  }


  std::vector<std::string>                         Keys ;

  private : 
  vtkGridType                                      *VTU; 
  std::map<std::string,vtkDoubleArray*>            Arrays;
  std::map<std::string,vtkDoubleArray*>::iterator  ArraysIT;

  int nDim ;
  int nPts ;
  int nCells ;
}; 


class InSituCpp 
{
  public :
  InSituCpp()
  {
    Processor   = NULL;
    Description = NULL;
    Grid        = NULL; 
    rank        = -1 ; 
  }


 ~InSituCpp()
  {
  }


  vtkUnstructuredGrid* GetMesh()
  {
    return Grid->GetMesh();
  }


  void SetMesh(vtkUnstructuredGrid* vtu, int dim)
  {  
    Grid->RemoveMesh();    
    Grid->SetMesh(vtu, dim); 
    if(!rank) std::cout<<"[SetMesh] \n";
  }  


  void CreateArray(std::string name, int ndata, int ncomps, double *data)
  {
    Grid->CreateArray(name, ndata, ncomps, data);
  }


  void UseArray(std::string name, int ndata, int ncomps, double *data)
  {
    Grid->UseArray(name, ndata, ncomps, data);
  }


  void Finalize( )
  {
    if(Processor)
    {
      Processor->Delete();
      Processor = NULL;
    }
    if(!rank) std::cout<<"[Finalize]\n";
  }


  void Initialize(MPI_Comm handle=MPI_COMM_WORLD, int outputFrequency=1, vtkProcessFunctionType process=NULL) 
  {
    MPI_Comm_rank(handle, &rank); 

  //int outputFrequency  = 1 ;
    std::string fileName = "output";

    if(Processor == NULL)
    {
      // SEE :  CoProcessing/Catalyst/Testing/Cxx/SubController.cxx 
      vtkMPICommunicatorOpaqueComm comm(&handle);
      Processor = vtkCPProcessor::New();
      Processor->Initialize(comm);
 
      // Vtk mpi example     
      if(process)
      {  
        vtkMPIController* controller = vtkMPIController::New();  
      //controller->SetCommunicator(comm); // To be tried!! 
      //controller->Initialize();
        controller->Initialize(NULL, NULL, 1); // initializedExternally == 1;
        controller->SetSingleMethod(process,0);   
        controller->SingleMethodExecute();  
        //controller->Finalize();
        //controller->Delete();
      }
    }

#ifdef WITH_PYTHON 
    // C++ pipeline 
    vtkNew<vtkCPPythonScriptPipeline> pipeline; 
    pipeline->Initialize("pipeline.py");   
    Processor->AddPipeline(pipeline.GetPointer());
#else 
    // C++ pipeline 
    vtkNew<vtkCPVTKPipeline> pipeline;
    pipeline->Initialize(outputFrequency, fileName);
    Processor->AddPipeline(pipeline.GetPointer());
#endif


    Grid = new VtkGrid<vtkUnstructuredGrid>();

    if(!rank) std::cout<<"[Initialize]\n";
  }


  void Coprocess(double time, unsigned int timeStep, int save) 
  {
    assert(Grid!=NULL);

    Description = vtkCPDataDescription::New();
    Description->AddInput("input");
    Description->SetTimeData(time, timeStep);

    bool request = Processor->RequestDataDescription( Description );

    if(!rank)
      std::cout<<"[Coprocess] "
               << "timeStep:"<< timeStep <<", "
               << "request:" << request <<" \n";
    //if(rank==0) Grid->PrintArraysNames();

    if(request)
    {
      vtkUnstructuredGrid* vtu = Grid->GetMesh();
      assert( vtu!=NULL );

      vtkCPInputDataDescription* idd = Description->GetInputDescriptionByName("input");
      for(std::vector<std::string>::iterator it=Grid->Keys.begin(); it!=Grid->Keys.end(); ++it)
      {
        std::string name(*it);
        if (idd->IsFieldNeeded(name.c_str(), vtkDataObject::POINT))
        {
          //std::cout<< name <<" \n" ;
          vtu->GetPointData()->AddArray( Grid->GetArray(name) );
/*
          if(vtu->GetPointData()->GetNumberOfArrays() == 0)
          { 
            vtu->GetPointData()->AddArray( Grid->GetArray(name) );
          } 
          else 
          {
            vtkDoubleArray* array = vtkDoubleArray::SafeDownCast(vtu->GetPointData()->GetArray(name));
            array->CreateArray(pressureData, static_cast<vtkIdType>(numberOfPoints), 1);
          }
*/
        }
      }
      //std::cout<< Grid->Keys.size()  <<" \n" ;
      Grid->Keys.clear();

      std::ostringstream o; o << timeStep;
      if(save) Grid->SaveMesh("mesh"+o.str()+".pvtu"); 

      idd->SetGrid(vtu);
      Processor->CoProcess( Description );
    } // request 
  }

  private :
  vtkCPProcessor                *Processor;  
  vtkCPDataDescription          *Description;
  VtkGrid<vtkUnstructuredGrid>  *Grid;  

  int rank; 
}; 

#endif 
// J. MIGUEL ZAVALA AKE. 2019JUL12. STOCKHOLM, SWEDEN.
//
// NOTE :
//      See smart pointers!!
//      https://vtk.org/Wiki/VTK/Tutorials/SmartPointers
//
