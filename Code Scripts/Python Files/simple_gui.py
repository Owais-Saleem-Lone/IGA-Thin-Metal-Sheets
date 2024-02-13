from OCC.Display.SimpleGui import *
from OCC.Display.backend import *
from OCC.Core.BRepPrimAPI import BRepPrimAPI_MakeBox
from OCC.Extend.TopologyUtils import TopologyExplorer
from OCC.Core.IGESControl import IGESControl_Reader
from OCC.Core import TopoDS
from OCC.Core.BRep import BRep_Builder
from OCC.Extend.DataExchange import read_iges_file
from OCC.Core.BRepBuilderAPI import BRepBuilderAPI_NurbsConvert

def load_iges(display, iges_file_path):
    # Initialize the IGES reader
    base_shape = read_iges_file(iges_file_path)
    nurbs_converter = BRepBuilderAPI_NurbsConvert(base_shape, True)
    basic_shape = nurbs_converter.Shape()

    display.DisplayShape(basic_shape, update=True)
def create_geometry(display):
    # Your code to create or load 3D geometry
    # For example, create a box
    box = BRepPrimAPI_MakeBox(10., 20., 30.).Shape()
    
    # Display the geometry
    display.DisplayShape(box, update=True)
    
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
      if shape.ShapeType() == 1:  # TopAbs_EDGE
        explorer = TopologyExplorer(shape)
        for edge in explorer.Current():
          print(f"    Edge: {edge}")

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
    iges_file_path = 'C:\\Users\\oslon\\Desktop\\MSc\\Thesis\\Code Scripts\\multi-patch plate.igs'
    
    # Load and display the IGES file
    load_iges(display, iges_file_path)
    display._select_callbacks.append(on_select)
    #create_geometry(display)

   
    start_display()