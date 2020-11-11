import os
import sys

sys.path.append(os.path.join(os.getcwd(), './'))

from vtk import *
import vtkXYZPython

collide = vtkXYZPython.vtkCollisionDetectionFilter()

print("TEST OK!!")
