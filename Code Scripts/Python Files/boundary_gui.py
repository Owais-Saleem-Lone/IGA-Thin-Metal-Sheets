import sys
from scipy.io import savemat
import numpy as np
from OCC.Core.gp import gp_Pnt2d,gp_Pnt
from OCC.Core.BRepBuilderAPI import BRepBuilderAPI_NurbsConvert,BRepBuilderAPI_MakeEdge, BRepBuilderAPI_MakeWire,BRepBuilderAPI_MakeFace
from OCC.Core.BRepAdaptor import BRepAdaptor_Curve,BRepAdaptor_Surface,BRepAdaptor_CompCurve,BRepAdaptor_Curve2d
from OCC.Core.IGESControl import IGESControl_Reader
from OCC.Extend.DataExchange import read_iges_file
from OCC.Extend.TopologyUtils import TopologyExplorer
from OCC.Core.GeomAbs import GeomAbs_BSplineSurface
from OCC.Core.GeomAbs import GeomAbs_BSplineCurve
from OCC.Core.BRep import BRep_Tool
from OCC.Core.Geom import Geom_BSplineCurve,Geom_BSplineSurface
from OCC.Core.Geom2dConvert import geom2dconvert  
from OCC.Core.TColgp import TColgp_Array1OfPnt2d
from OCC.Core.TColStd import TColStd_Array1OfInteger, TColStd_Array1OfReal
from OCC.Core.TopExp import TopExp_Explorer  
from OCC.Core.TopAbs import TopAbs_EDGE,TopAbs_WIRE,TopAbs_FACE
from OCC.Core.TopoDS import topods,TopoDS_Builder
from OCC.Core.ShapeAnalysis import shapeanalysis
from OCC.Display.SimpleGui import init_display
from OCC.Display.SimpleGui import *
from OCC.Display.backend import *
from OCC.Core import TopoDS,TopAbs
from OCC.Core.Geom2dAPI import Geom2dAPI_ProjectPointOnCurve
import time



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
  
def on_select_edge(selected_shapes, x, y):
    global selected_edge
    for shape in selected_shapes:
        if shape.ShapeType() == TopAbs.TopAbs_EDGE:
            edge = topods.Edge(shape)
            selected_edge.append(edge)
    if len(selected_edge) == 2 and len(selected_face) == 2:
     process_intersection(display)  

def on_select(selected_shapes,x,y):
    global selected_edge, selected_face
    for shape in selected_shapes:
      if shape.ShapeType() == TopAbs.TopAbs_FACE:  # TopAbs_FACE
        face=topods.Face(shape)
        selected_face.append(face)
        display.EraseAll()
        display.DisplayShape(selected_face)
        display._select_callbacks.remove(on_select)
        display._select_callbacks.append(on_select_edge)

def ed_fc_inter(edge,face):
    global matrix_dict
    curve_points_elem=[]
    OCC_handle_curve2d, activeRangeBegin,activeRangeEnd = BRep_Tool.CurveOnSurface(topods.Edge(edge), topods.Face(face))
    OCC_bsplineCurve2d = geom2dconvert.CurveToBSplineCurve(OCC_handle_curve2d)
    curve = BRepAdaptor_Curve(edge)
    bcurve= curve.BSpline()
    curve_deg=bcurve.Degree()
    curve_knots=bcurve.Knots()
    KnotVector=[]
    for i in range(bcurve.NbKnots()):
      KnotVector.append(curve_knots.Value(i + 1))
    print(f"The edge KnotVector is: {KnotVector}\n")
    surface = BRepAdaptor_Surface(face, True)
    bsurf=surface.BSpline()
    bsurf.ExchangeUV()
    uknots = bsurf.UKnots()
    UknotArr=[]
    for i in range(bsurf.NbUKnots()):
      UknotArr.append(uknots.Value(i+1))
    UknotArr_ref= np.round(np.linspace(float(UknotArr[0]),float(UknotArr[-1]),(num_xiKnots+2)),2)
    print(f"The refined xi knot vector is : {UknotArr_ref}")
    vknots = bsurf.VKnots()
    VknotArr=[]
    for i in range(bsurf.NbVKnots()):
        VknotArr.append(vknots.Value(i+1))
    VknotArr_ref= np.round(np.linspace(float(VknotArr[0]),float(VknotArr[-1]),(num_etaKnots+2)),2)
    print(f"The refined eta knot vector is : {VknotArr_ref}")
    ctr=0
    for m in range(0,len(UknotArr_ref)):
      xi = UknotArr_ref[m] 
      ctr=ctr+1
      min_distance = float('inf')
      min_distance_projector = None
      for n in range(0,len(VknotArr_ref)):
        if n+ctr-1<len(VknotArr_ref):
          eta = VknotArr_ref[n+ctr-1]
          point2d = gp_Pnt2d(eta,xi)    # Because of reverse parametrization of geometry_curve
          projector = Geom2dAPI_ProjectPointOnCurve(point2d, OCC_bsplineCurve2d) # Orthogonal Projection
          distance= round(projector.Distance(1),2)
          #print(f"the point is [{xi},{eta}] and the distance is : {distance}")
          if distance < min_distance:
            min_distance = distance
            min_distance_projector = projector
            
      sol=min_distance_projector.NbPoints()
      for p in range(1,sol+1):
        if sol==2:
          print(f"The possible number of projections are: {sol}")
        t = round(min_distance_projector.Parameter(sol),2) # The first solution among the possible projections
        curve_points_elem.append(t)

    curve_points_elem=np.array(curve_points_elem)
    curve_points_elem=np.unique(curve_points_elem)
    curve_points_elem=np.sort(curve_points_elem)
        
    print(f"The list of curve element intersections are: {curve_points_elem}")
    print(f"Now Inserting {curve_deg+1} Gauss Points per element for the curve with degree {curve_deg}")
    
    new_values = []
    xi_values=[]
    eta_values=[]
    for i in range(len(curve_points_elem) - 1):
      values_between = np.linspace(curve_points_elem[i], curve_points_elem[i + 1], curve_deg+2 + 2)[1:-1]
      new_values.extend([curve_points_elem[i]] + list(values_between) + [curve_points_elem[i + 1]])
    new_array = np.round(np.unique(np.array(new_values)),2)
    
    for j in range(len(new_array)):
      point_on_curve = OCC_bsplineCurve2d.Value(new_array[j])
      xi_values.append(round(point_on_curve.X(),2))
      eta_values.append(round(point_on_curve.Y(),2))
      #print(f"For edge point {new_array[j]} the value of surface xi is {eta_values} and surface eta is {xi_values}\n")
    
    matrix_dict = {'curve_points': new_array, 'xi_values': eta_values, 'eta_values': xi_values}
    for key, value in matrix_dict.items():
      print(f"{key}: {value}")              


def process_intersection(display):

    ed_fc_inter(selected_edge[0],selected_face[0])
    display.EraseAll()
    display.DisplayMessage(gp_Pnt(0,0,0),"Process completed successfully!",50, (0,0,1),False) 
    display.DisplayMessage(gp_Pnt(20,-2,0),"You can Exit Now",30, (1,0,0),False)

   
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
    global selected_edge, selected_face
    selected_edge=[]
    selected_face=[]
    #iges_file_path="C:\\Users\\oslon\\Desktop\\MSc\\Thesis\\Code Scripts\\multi-patch sep.igs"
    iges_file_path = sys.argv[1]
    global num_xiKnots, num_etaKnots
    global matrix_dict,output_file,operation
    output_file= 'b_rep.mat'
    #operation= sys.argv[2]
    num_xiKnots=int(sys.argv[3])
    num_etaKnots=int(sys.argv[4])
    #num_xiKnots=4
    #num_etaKnots=6
    print(f"xi refinement points are: {num_xiKnots}")
    print(f"eta refinement points are <with increment of course>: {num_etaKnots}")
    load_iges(display, iges_file_path)
    display._select_callbacks.append(on_select)
    start_display()
    savemat(output_file, matrix_dict)