import sys
from scipy.io import savemat
import numpy as np
from OCC.Core.gp import gp_Pnt2d,gp_Pnt
from OCC.Core.BRepBuilderAPI import BRepBuilderAPI_NurbsConvert,BRepBuilderAPI_MakeEdge, BRepBuilderAPI_MakeWire,BRepBuilderAPI_MakeFace
from OCC.Core.BRepAdaptor import BRepAdaptor_Curve,BRepAdaptor_Surface,BRepAdaptor_CompCurve,BRepAdaptor_Curve2d
from OCC.Core.IGESControl import IGESControl_Reader
from OCC.Extend.DataExchange import read_iges_file
from OCC.Extend.TopologyUtils import TopologyExplorer
from OCC.Core.GeomAbs import GeomAbs_BSplineCurve,GeomAbs_BSplineSurface
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
from OCC.Core.GeomAPI import GeomAPI_ProjectPointOnSurf
from OCC.Core.Extrema import Extrema_ExtAlgo
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
    global selected_edges,iges_file_path
    for shape in selected_shapes:
        if shape.ShapeType() == TopAbs.TopAbs_EDGE:
            edge = topods.Edge(shape)
            selected_edges.append(edge)
    if len(selected_edges) == 2 and len(selected_faces) == 2:
      process_intersection(display)
      
    else:
      display.EraseAll()
      load_iges(display, iges_file_path)
      display._select_callbacks.remove(on_select_edge)
      display._select_callbacks.append(on_select)

def on_select(selected_shapes,x,y):
    global selected_faces
    print("Selected Shapes:")
    for shape in selected_shapes:

      if shape.ShapeType() == TopAbs.TopAbs_FACE:  # TopAbs_FACE
          face=topods.Face(shape)
          selected_faces.append(face)
          display.EraseAll()
          display.DisplayShape(face)
          display._select_callbacks.remove(on_select)
          display._select_callbacks.append(on_select_edge)
          #print(f"Length of faces and edges is: {len(selected_faces)} and {len(selected_edges)}")

def ed_fc_inter(edge,face):
  OCC_handle_curve2d, activeRangeBegin,activeRangeEnd = BRep_Tool.CurveOnSurface(topods.Edge(edge), topods.Face(face))
  OCC_bsplineCurve2d = geom2dconvert.CurveToBSplineCurve(OCC_handle_curve2d)
  curve_points_elem=[]
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
  print("surface curve intersection found succesfully\n")
  print(curve_points_elem)
  return curve_points_elem
  
def surf_to_phy(xi_vec,eta_vec,face):
  surface_adaptor=BRepAdaptor_Surface(face)
  x_vec=[]
  y_vec=[]
  z_vec=[]
  for j in range(len(xi_vec)):
    point_on_surface=surface_adaptor.Value(eta_vec[j],xi_vec[j])
    x,y,z=point_on_surface.X(),point_on_surface.Y(),point_on_surface.Z()
    x_vec.append(x)
    y_vec.append(y)
    z_vec.append(z)
    print(f"the surface point is [{xi_vec[j]},{eta_vec[j]}] and the physical point is [{x},{y},{z}]")
  print("Projected surface to physical successfully\n")

  return x_vec,y_vec,z_vec
      
def phy_to_surf(x_vec,y_vec,z_vec,face):
  print("In physical to surface")   
  print(len(x_vec))
  xi_vec=[]
  eta_vec=[]
  for j in range(len(x_vec)):
    point_in_physical_space = gp_Pnt(x_vec[j], y_vec[j], z_vec[j])
    face_surface = BRep_Tool.Surface(face)
    projector = GeomAPI_ProjectPointOnSurf(point_in_physical_space, face_surface)
    point_on_surface = projector.NearestPoint()
    u, v = projector.LowerDistanceParameters()
    print(f"Surface Parameters of physical point [{x_vec[j]}{y_vec[j]}]: ({v}, {u})")
    xi_vec.append(round(u,2))
    eta_vec.append(round(v,2))
  return eta_vec,xi_vec
       
def surf_to_curve(xi_vec,eta_vec,edge,face):
  print("In surface to curve")
  OCC_handle_curve2d, activeRangeBegin,activeRangeEnd = BRep_Tool.CurveOnSurface(topods.Edge(edge), topods.Face(face))
  OCC_bsplineCurve2d = geom2dconvert.CurveToBSplineCurve(OCC_handle_curve2d)
  crv_vec=[]
  for j in range(len(xi_vec)):
    point2d = gp_Pnt2d(eta_vec[j],xi_vec[j]) 
    projector = Geom2dAPI_ProjectPointOnCurve(point2d, OCC_bsplineCurve2d)
    t = round(projector.Parameter(1),2)
    crv_vec.append(t)
  print("the projected curve points are: ")
  print(crv_vec)
  return crv_vec

def refine_curve(points,edge):
  curve = BRepAdaptor_Curve(edge)
  bcurve= curve.BSpline()
  curve_deg=bcurve.Degree()
  new_values = []
  for i in range(len(points) - 1):
    values_between = np.linspace(points[i], points[i + 1], curve_deg+2 + 2)[1:-1]
    new_values.extend([points[i]] + list(values_between) + [points[i + 1]])
  new_array = np.round(np.unique(np.array(new_values)),2)
  print("Curve refined successfully\n")
  print(new_array)
  return new_array

def crv_to_srf(crv_pts,edge ,face):
  OCC_handle_curve2d, activeRangeBegin,activeRangeEnd = BRep_Tool.CurveOnSurface(topods.Edge(edge), topods.Face(face))
  OCC_bsplineCurve2d = geom2dconvert.CurveToBSplineCurve(OCC_handle_curve2d)
  xi_values=[]
  eta_values=[]
  for j in range(len(crv_pts)):
    point_on_curve = OCC_bsplineCurve2d.Value(crv_pts[j])
    xi_values.append(round(point_on_curve.X(),2))
    eta_values.append(round(point_on_curve.Y(),2))
  print("Projected points from curve to surface")
  print(eta_values,"\n")
  print(xi_values,"\n")
  return eta_values,xi_values  
  
def nurbs_detail_coupling_curve(crv,surf):
  OCC_handle_curve2d, activeRangeBegin, activeRangeEnd = BRep_Tool.CurveOnSurface(topods.Edge(crv), topods.Face(surf))
  OCC_bsplineCurve2d = geom2dconvert.CurveToBSplineCurve(OCC_handle_curve2d)
  
  if OCC_bsplineCurve2d.IsPeriodic(): 
    print("Found periodic curve!")
    OCC_bsplineCurve2d.SetNotPeriodic()
  # pDegree
  pDegree = OCC_bsplineCurve2d.Degree()
  # Knots
  user_refinement=8
  refined_values = np.round(np.linspace(OCC_bsplineCurve2d.FirstParameter(), OCC_bsplineCurve2d.LastParameter(), user_refinement)[1:-1])
  knot_insertion = TColStd_Array1OfReal(1, len(refined_values))
  knot_insertion_mult = TColStd_Array1OfInteger(1, len(refined_values))
  for i, value in enumerate(refined_values):
    knot_insertion.SetValue(i + 1, value)
    knot_insertion_mult.SetValue(i + 1, 1)
  OCC_bsplineCurve2d.InsertKnots(knot_insertion,knot_insertion_mult)
    
  OCC_uNoKnots = OCC_bsplineCurve2d.NbKnots()
  # Knot vector
  OCC_uKnotMulti = TColStd_Array1OfInteger(1,OCC_uNoKnots)
  OCC_bsplineCurve2d.Multiplicities(OCC_uKnotMulti)
  uNoKnots = 0
  for ctr in range(1,OCC_uNoKnots+1): uNoKnots += OCC_uKnotMulti.Value(ctr)
  OCC_uKnotSequence = TColStd_Array1OfReal(1,uNoKnots)
  OCC_bsplineCurve2d.KnotSequence(OCC_uKnotSequence)
  uKnotVector = []
  for iKnot in range(1,uNoKnots+1): uKnotVector.append(OCC_uKnotSequence.Value(iKnot))
  
  
  
  # No of CPs
  uNoCPs = OCC_bsplineCurve2d.NbPoles()
  # CPs
  OCC_CPnet = TColgp_Array1OfPnt2d(1,uNoCPs)
  OCC_bsplineCurve2d.Poles(OCC_CPnet)
  # CP weights
  OCC_CPweightNet = TColStd_Array1OfReal(1,uNoCPs)
  OCC_bsplineCurve2d.Weights(OCC_CPweightNet)
  CPNet = []
  for iCP in range(1,uNoCPs+1):
    OCC_CP = OCC_CPnet.Value(iCP)
    OCC_CPweight = OCC_CPweightNet.Value(iCP)
    CPNet.extend([OCC_CP.Y(), OCC_CP.X(), 0.0, OCC_CPweight])
  return pDegree,uKnotVector,CPNet

def process_intersection(display):
  global matrix_dict
  crv_pts1=ed_fc_inter(selected_edges[0],selected_faces[0])
  display.EraseAll()
  display.DisplayMessage(gp_Pnt(0,0,0),"Patches coupled successfully!",50, (0,0,1),False) 
  display.DisplayMessage(gp_Pnt(20,-2,0),"You can Exit Now",30, (1,0,0),False)
  ref_crv_pts1= refine_curve(crv_pts1,selected_edges[0])
  xi_values,eta_values= crv_to_srf(ref_crv_pts1,selected_edges[0],selected_faces[0])
  x_inter,y_inter,z_inter= surf_to_phy(xi_values,eta_values,selected_faces[0])
  xi_other,eta_other = phy_to_surf(x_inter,y_inter,z_inter,selected_faces[1])
  ref_crv_pts2= surf_to_curve(xi_other,eta_other,selected_edges[1],selected_faces[1])
  crv_deg,crv_knt_vec,crv_CPs= nurbs_detail_coupling_curve(selected_edges[0],selected_faces[0])
  
  matrix_dict = {'cur_pts_1': ref_crv_pts1, 'xi_vls_1': xi_values, 'eta_vls_1': eta_values,'cur_pts_2': ref_crv_pts2, 'xi_vls_2': xi_other, 'eta_vls_2': eta_other
                 ,'curve_deg':crv_deg,'crv_knotVct':crv_knt_vec,'crv_cps':crv_CPs}
  for key, value in matrix_dict.items():
    print(f"{key}: {value}")
      
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
    global selected_edges, selected_faces,iges_file_path
    selected_edges=[]
    selected_faces=[]
    iges_file_path="C:\\Users\\oslon\\Desktop\\MSc\\Thesis\\Code Scripts\\multi-patch plate.igs"
    #iges_file_path = sys.argv[1]
    global num_xiKnots, num_etaKnots
    global matrix_dict,output_file,operation 
    output_file= 'b_rep.mat'
    #operation= sys.argv[2]
   #num_xiKnots=int(sys.argv[3])
   #num_etaKnots=int(sys.argv[4])
    num_xiKnots=5
    num_etaKnots=5
    print(f"xi refinement points are: {num_xiKnots}")
    print(f"eta refinement points are <with increment of course>: {num_etaKnots}")
    load_iges(display, iges_file_path)
    display._select_callbacks.append(on_select)
    start_display()
    savemat(output_file, matrix_dict)