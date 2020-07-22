#ifndef IN_SITU_PY_VTK_H  
#define IN_SITU_PY_VTK_H
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

// C++ 
//#include "pipeline.hpp" // vtkCPVTKPipeline.cxx 
#include <vtkCPProcessor.h>
#include <vtkCPDataDescription.h>
#include <vtkCPInputDataDescription.h>
#include <vtkPointData.h>

// Python 
#include <vtkCPPythonScriptPipeline.h>

void TestProcessController(vtkMultiProcessController*, void*);

//--------------------------------------------------------------------------||--//
// SEE : 
// /afs/pdc.kth.se/home/m/miguelza/z2020_1/RUNNER/ADAM/Wtip/Togheter01_03
//
template<class VTKCONTROLLER>
class CreateVtkMPIController
{
  public :
 ~CreateVtkMPIController()
  {
    initialized = 0;
  }


  CreateVtkMPIController(MPI_Comm gcomm, int print=0)
  {
    assert(!vtkcontr);
    assert(gcomm!=MPI_COMM_NULL);
    mpi_comm = gcomm;

    int numprocs = -1;
    MPI_Comm_size(gcomm, &numprocs);

    vtkcontr = NULL;
    if( this->Initialized() )
    {
      //vtkMPICommunicatorOpaqueComm *vtkcomm = new vtkMPICommunicatorOpaqueComm(&lcomm);
      vtkMPICommunicatorOpaqueComm *vtkcomm = new vtkMPICommunicatorOpaqueComm(&gcomm);

      vtkMPICommunicator *comm = vtkMPICommunicator::New();
      comm->InitializeExternal( vtkcomm );

      vtkcontr = VTKCONTROLLER::New();
      vtkcontr->SetCommunicator( comm );
    }
    assert(vtkcontr);
    vtkMultiProcessController::SetGlobalController(vtkcontr); // [4]   
 
    vtkcontr->Barrier();
    nranks = vtkcontr->GetNumberOfProcesses();
    rank   = vtkcontr->GetLocalProcessId();
    if(!rank && print) std::cout<<"\t[CreateVtkMPIController] "<< nranks <<" Processes Initialized!!\n";

  }


  void Finalize() // APPLY MPI_Finalize() !!
  {
    vtkcontr->Finalize();
  }


  int Initialized()
  {
    initialized = -1;
    MPI_Initialized(&initialized);
    return initialized;
  }


  VTKCONTROLLER* GetController()
  {
    assert(vtkcontr);
    return vtkcontr;
  }


  vtkMultiProcessController* GetMultiProcessController()
  {
    assert(vtkcontr);

    vtkMultiProcessController *controller = vtkMultiProcessController::GetGlobalController();  // [3] vtkCPPythonScriptPipeline::Initialize
    assert(controller);                                                                        // [4] GetGlobalController == 
                                                                                               //     VTK_GLOBAL_MULTI_PROCESS_CONTROLLER->GetLocalController  
    vtkCommunicator *vtkcomm = controller->GetCommunicator(); 
    assert(vtkcomm);   

    controller->GetLocalProcessId();  // [4] Segmentation fault - invalid memory reference.
    //TestProcessController(controller, NULL);

    return controller;
  }


  void ProcessFunction(vtkProcessFunctionType process=NULL) 
  {
    assert(vtkcontr);
    assert(process);  
    vtkcontr->Initialize(NULL, NULL, 1); // initializedExternally == 1;
    vtkcontr->SetSingleMethod(process,0);  // -> ' MPI has to be initialized first' 
    vtkcontr->SingleMethodExecute();    
  }

  private :
  VTKCONTROLLER *vtkcontr = NULL;
  int rank   = -1;
  int nranks = -1;
  int initialized = 0;

  MPI_Comm mpi_comm;
};

CreateVtkMPIController<vtkMPIController> *MPIController = NULL;



// /* 
//--------------------------------------------------------------------------||--//
template <class vtkGridType>
void PWriter(vtkGridType* DataObject, std::string fname)
{
// /* 
  auto contr = MPIController->GetController();  
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
//parallel_writer->SetDataModeToAscii();

  parallel_writer->SetUseSubdirectory(1); 
  parallel_writer->Update();
  parallel_writer->Write();

  contr->Finalize(1); // if(finalizedExternally == 0) MPI_Finalize();
// */ 
}



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
    ProcessorCpp= NULL;
    Description = NULL;
    Grid        = NULL; 
    rank        = -1;
    size        = -1; 
    itime       = -1;  
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
    //PWriter<vtkUnstructuredGrid>( vtu, "input_mesh.pvtu" );
  }  


  void SaveMesh(std::string fname) 
  {
    std::ostringstream o; o <<"_"<< itime;  
    Grid->SaveMesh(fname + o.str() + ".pvtu");
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
    if(ProcessorCpp)
    {
      ProcessorCpp->Delete();
      ProcessorCpp = NULL;
    }
    if(!rank) std::cout<<"[Finalize]\n";
  }


  void Initialize(MPI_Comm handle, int outputFrequency, std::string pyfile, vtkProcessFunctionType process) 
  {
    MPI_Comm_rank(handle, &rank); 
    MPI_Comm_size(handle, &size);

    MPIController   = new CreateVtkMPIController<vtkMPIController>(handle, 1);
    auto controller = MPIController->GetMultiProcessController(); //  Assertion `controller' failed. 
    process(controller, NULL); // Test  

    if(ProcessorCpp == NULL)
    {
      // SEE :  CoProcessing/Catalyst/Testing/Cxx/SubController.cxx 
      vtkMPICommunicatorOpaqueComm comm(&handle);
      ProcessorCpp = vtkCPProcessor::New();
      ProcessorCpp->Initialize(comm);
    }
/*
    // C++ pipeline 
    std::string fout = "output";
    vtkNew<vtkCPVTKPipeline> pipelineCpp;
    pipelineCpp->Initialize(outputFrequency, fout); //fileName);
    ProcessorCpp->AddPipeline(pipelineCpp.GetPointer());
// */   
    // Python pipeline 
    vtkNew<vtkCPPythonScriptPipeline> pipelinePy;   
    pipelinePy->Initialize(pyfile.c_str()); // read the coprocessing python file 
//    pipelinePy->Initialize("pipe.py"); //fileName.c_str()); // read the coprocessing python file 
    ProcessorCpp->AddPipeline(pipelinePy.GetPointer()); // vtkCPPythonAdaptorAPI ??  
// */
    // Init Grid 
    Grid = new VtkGrid<vtkUnstructuredGrid>();
    if(!rank) std::cout<<"[Initialize] nRanks:"<< size
	               <<" pyfile:'"<< pyfile <<"' "   
	               <<"\n";  
  }


  void Coprocess(double time, unsigned int timeStep) 
  {
    assert(Grid!=NULL);

    itime = timeStep;

    Description = vtkCPDataDescription::New();
    Description->AddInput("input");
    Description->SetTimeData(time, timeStep);

    bool request = ProcessorCpp->RequestDataDescription( Description );

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
/*  
      std::ostringstream o; o << size <<"_"<<timeStep;
      if(save) Grid->SaveMesh("mesh"+o.str()+".pvtu"); 
*/
      idd->SetGrid(vtu);
      ProcessorCpp->CoProcess( Description );
    } // request 
  }

  private :
  vtkCPProcessor                *ProcessorCpp;  
  vtkCPDataDescription          *Description;
  VtkGrid<vtkUnstructuredGrid>  *Grid;  

  int rank; 
  int size;   
  int itime;   
}; 

#endif 
// J. MIGUEL ZAVALA AKE. 2019JUL12. STOCKHOLM, SWEDEN.
//
// NOTE :
//      See smart pointers!!
//      https://vtk.org/Wiki/VTK/Tutorials/SmartPointers
//
