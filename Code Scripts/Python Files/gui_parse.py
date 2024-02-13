from OCC.Display.SimpleGui import *
from OCC.Display.backend import *
from OCC.Core.BRepPrimAPI import BRepPrimAPI_MakeBox
from OCC.Extend.TopologyUtils import TopologyExplorer
from OCC.Core.IGESControl import IGESControl_Reader
from OCC.Core import TopoDS,TopAbs
from OCC.Core.BRepAdaptor import BRepAdaptor_Curve,BRepAdaptor_Surface,BRepAdaptor_CompCurve
from OCC.Extend.DataExchange import read_iges_file
from OCC.Core.BRepBuilderAPI import BRepBuilderAPI_NurbsConvert
import numpy as np
import sys
def load_iges(display, iges_file_path):
    # Initialize the IGES reader
    base_shape = read_iges_file(iges_file_path)
    nurbs_converter = BRepBuilderAPI_NurbsConvert(base_shape, True)
    basic_shape = nurbs_converter.Shape()

    display.DisplayShape(basic_shape, update=True)

    
def Clean(event=None):
  display.EraseAll()
  
def Exit(event=None):
  import sys
  sys.exit()

def Neutral(event=None):
  display.SetSelectionModeNeutral()

def Edge(event=None):
  display.SetSelectionModeEdge()

def Face(event=None):
  display.SetSelectionModeFace()

def Vertex(event=None):
  display.SetSelectionModeVertex()
  
def on_select(selected_shapes,x,y):
    
    print("Selected Shapes:")


    for shape in selected_shapes:
      print(f"  {shape}")

      # Example: Get information about edges
      if shape.ShapeType() == TopAbs.TopAbs_EDGE:  # TopAbs_EDGE
        nurbs_converter = BRepBuilderAPI_NurbsConvert(shape, True)
        converted_shape = nurbs_converter.Shape()
        explorer = TopologyExplorer(converted_shape)
        for edge in explorer.edges():
          print(f"  owais  Edge: {edge}")
          curve = BRepAdaptor_Curve(edge)
          bcurve = curve.BSpline()
          poles = bcurve.Poles()
          controlPointMatrix= np.zeros((bcurve.NbPoles(),3))
          if poles is not None:
            for i in range(bcurve.NbPoles()):
              p = poles.Value(i + 1)
              controlPointMatrix[i]= np.array([p.X(),p.Y(),p.Z()])
            print(controlPointMatrix)
      # Example: Get information about faces
      elif shape.ShapeType() == 2:  # TopAbs_FACE
        explorer = TopologyExplorer(shape)
        for face in explorer.Current():
          print(f"    Face: {face}")

if __name__ == '__main__':

    display, start_display, add_menu, add_function_to_menu = init_display('wx')

    add_menu('File')
    add_function_to_menu('File', Clean)
    add_function_to_menu('File', Exit)

    add_menu('Selection')
    add_function_to_menu('Selection', Edge)
    add_function_to_menu('Selection', Face)
    add_function_to_menu('Selection', Vertex)
    add_function_to_menu('Selection', Neutral)

    iges_file_path = sys.argv[1]

    #iges_file_path = 'C:\\Users\\oslon\\Desktop\\MSc\\Thesis\\Code Scripts\\multi-patch plate.igs'
    #iges_file_path = 'C:\\Users\\oslon\\Desktop\\MSc\\Thesis\\Code Scripts\\onePatch.igs'
    # Load and display the IGES file
    load_iges(display, iges_file_path)
    display._select_callbacks.append(on_select)
    #create_geometry(display)

   
    start_display()