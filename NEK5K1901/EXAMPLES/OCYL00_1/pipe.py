
#--------------------------------------------------------------

# Global timestep output options
timeStepToStartOutputAt=0
forceOutputAtFirstCall=False

# Global screenshot output options
imageFileNamePadding=0
rescale_lookuptable=False

# Whether or not to request specific arrays from the adaptor.
requestSpecificArrays=False

# a root directory under which all Catalyst output goes
rootDirectory=''

# makes a cinema D index table
make_cinema_table=False

#--------------------------------------------------------------
# Code generated from cpstate.py to create the CoProcessor.
# paraview version 5.6.0
#--------------------------------------------------------------

from paraview.simple import *
from paraview import coprocessing
import vtk

# ----------------------- CoProcessor definition -----------------------

def CreateCoProcessor():
  def _CreatePipeline(coprocessor, datadescription):
    class Pipeline:
      # state file generated using paraview version 5.6.0

      # ----------------------------------------------------------------
      # setup the data processing pipelines
      # ----------------------------------------------------------------

      # trace generated using paraview version 5.6.0
      #
      # To ensure correct image size when batch processing, please search 
      # for and uncomment the line `# renderView*.ViewSize = [*,*]`

      #### disable automatic camera reset on 'Show'
      paraview.simple._DisableFirstRenderCameraReset()

      # create a new 'VisItNek5000Reader'
      # create a producer from a simulation input
      ocylnek5000 = coprocessor.CreateProducer(datadescription, 'input')

      # ----------------------------------------------------------------
      me = __file__.split(".")[0] 

      # MPI 
      contr  = vtk.vtkMultiProcessController.GetGlobalController()
      nranks = contr.GetNumberOfProcesses() if contr else 1
      rank   = contr.GetLocalProcessId()    if contr else 0
      if not rank : print("\t[%s1] nranks: %d" % (me, nranks)) 

      vtkcommunicator    = contr.GetCommunicator()
      nranks = vtkcommunicator.GetNumberOfProcesses() if contr else 1
      rank   = vtkcommunicator.GetLocalProcessId()    if contr else 0
      if not rank : print("\t[pipe2] nranks: %d" % (nranks) )

      vtkmpicommunicator = vtkcommunicator.GetWorldCommunicator() 
      ##vtkmpicommunicatoropaquecomm = vtkmpicommunicator.GetMPIComm() # Program received signal SIGSEGV: Segmentation fault - invalid memory reference.
      nranks = vtkmpicommunicator.GetNumberOfProcesses() if contr else 1
      rank   = vtkmpicommunicator.GetLocalProcessId()    if contr else 0
      if not rank : print("\t[pipe3] nranks: %d" % (nranks) )

      ## https://www.paraview.org/Wiki/Python_coprocessing_example
      pm = paraview.servermanager.vtkProcessModule.GetProcessModule()
      globalController = pm.GetGlobalController()
      nranks = globalController.GetNumberOfProcesses() if contr else 1
      rank   = globalController.GetLocalProcessId()    if contr else 0
      if not rank : print("\t[pipe4] nranks: %d" % (nranks) )

      #vtkpvoptions = pm.GetOptions()
      #print( pm.GetSelfDir()   )

      # ----------------------------------------------------------------
      # finally, restore active source
      SetActiveSource(ocylnek5000)
      # ----------------------------------------------------------------

      # Now any catalyst writers
      xMLMultiBlockDataWriter1 = servermanager.writers.XMLPUnstructuredGridWriter(Input=ocylnek5000)
      coprocessor.RegisterWriter(xMLMultiBlockDataWriter1, filename=me+'_grid%t.pvtu', freq=5, paddingamount=0)

      # create a new 'Plot Over Line'
      plotOverLine1 = PlotOverLine(Input=ocylnek5000, Source='High Resolution Line Source')

      # init the 'High Resolution Line Source' selected for 'Source'
      plotOverLine1.Source.Point1 = [ 0.5, 0.0, 0.0]
      plotOverLine1.Source.Point2 = [10.0, 0.0, 0.0]

      # Now any catalyst writers
      xMLPPolyDataWriter1 = servermanager.writers.XMLPPolyDataWriter(Input=plotOverLine1)
      coprocessor.RegisterWriter(xMLPPolyDataWriter1, filename=me+'_line%t.pvtp', freq=1, paddingamount=0)


    return Pipeline()

  class CoProcessor(coprocessing.CoProcessor):
    def CreatePipeline(self, datadescription):
      self.Pipeline = _CreatePipeline(self, datadescription)

  coprocessor = CoProcessor()
  # these are the frequencies at which the coprocessor updates.
  freqs = {'input': [1]}
  coprocessor.SetUpdateFrequencies(freqs)
  if requestSpecificArrays:
    arrays = [['x_velocity', 0], ['y_velocity', 0]]
    coprocessor.SetRequestedArrays('input', arrays)
  coprocessor.SetInitialOutputOptions(timeStepToStartOutputAt,forceOutputAtFirstCall)

  if rootDirectory:
      coprocessor.SetRootDirectory(rootDirectory)

  if make_cinema_table:
      coprocessor.EnableCinemaDTable()

  return coprocessor


#--------------------------------------------------------------
# Global variable that will hold the pipeline for each timestep
# Creating the CoProcessor object, doesn't actually create the ParaView pipeline.
# It will be automatically setup when coprocessor.UpdateProducers() is called the
# first time.
coprocessor = CreateCoProcessor()

#--------------------------------------------------------------
# Enable Live-Visualizaton with ParaView and the update frequency
coprocessor.EnableLiveVisualization(False, 1)

# ---------------------- Data Selection method ----------------------

def RequestDataDescription(datadescription):
    "Callback to populate the request for current timestep"
    global coprocessor

    # setup requests for all inputs based on the requirements of the
    # pipeline.
    coprocessor.LoadRequestedData(datadescription)

# ------------------------ Processing method ------------------------

def DoCoProcessing(datadescription):
    "Callback to do co-processing for current timestep"
    global coprocessor

    # Update the coprocessor by providing it the newly generated simulation data.
    # If the pipeline hasn't been setup yet, this will setup the pipeline.
    coprocessor.UpdateProducers(datadescription)

    # Write output data, if appropriate.
    coprocessor.WriteData(datadescription);

    # Write image capture (Last arg: rescale lookup table), if appropriate.
    coprocessor.WriteImages(datadescription, rescale_lookuptable=rescale_lookuptable,
        image_quality=0, padding_amount=imageFileNamePadding)

    # Live Visualization, if enabled.
    coprocessor.DoLiveVisualization(datadescription, "localhost", 22222)
