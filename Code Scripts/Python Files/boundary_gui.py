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
from OCC.Core.AIS import AIS_InteractiveObject



def load_iges(display, iges_file_path):
    # Initialize the IGES reader
    base_shape = read_iges_file(iges_file_path)
    nurbs_converter = BRepBuilderAPI_NurbsConvert(base_shape, True)
    basic_shape = nurbs_converter.Shape()

    display.DisplayShape(basic_shape, update=True)

def Clean(event=None):
  display.EraseAll()
  
def Exit(event=None):
  process_intersection(display)
  sys.exit()

def Neutral(event=None):
  display.SetSelectionModeNeutral()

def Edge(event=None):
  display.SetSelectionModeEdge()

def Face(event=None):
    display.SetSelectionModeFace()

def Vertex(event=None):
  display.SetSelectionModeVertex()
 
def get_face_ids(ref_face):
  tol=1e-3
  base_shape = read_iges_file(iges_file_path)
  nurbs_converter = BRepBuilderAPI_NurbsConvert(base_shape, True)
  converted_shape = nurbs_converter.Shape()
  expl = TopologyExplorer(converted_shape)
  
  surface2 = BRepAdaptor_Surface(ref_face)
  u_mid = (surface2.FirstUParameter() + surface2.LastUParameter()) / 2
  v_mid = (surface2.FirstVParameter() + surface2.LastVParameter()) / 2
  point2 = surface2.Value(u_mid,v_mid)
  physical_point_2=gp_Pnt(point2.X(), point2.Y(), point2.Z())
  fc_idx = 1
  for face in expl.faces():
    surface1 = BRepAdaptor_Surface(face)
    u_mid1 = (surface1.FirstUParameter() + surface1.LastUParameter()) / 2
    v_mid1 = (surface1.FirstVParameter() + surface1.LastVParameter()) / 2
    point1 = surface1.Value(u_mid1,v_mid1)
    physical_point_1=gp_Pnt(point1.X(), point1.Y(), point1.Z())
    if (physical_point_1.Distance(physical_point_2) < tol) :
      return fc_idx
    fc_idx += 1
   
def on_select_edge(selected_shapes, x, y):
    global selected_edge,selected_face
    for shape in selected_shapes:
        if shape.ShapeType() == TopAbs.TopAbs_EDGE:
            edge = topods.Edge(shape)
            selected_edge.append(edge)

    display.EraseAll()
    #display.selected_shapes.clear()
    load_iges(display, iges_file_path)
    display._select_callbacks.remove(on_select_edge)
    display._select_callbacks.append(on_select)

def on_select(selected_shapes,x,y):
  global selected_edge, selected_face
  for shape in selected_shapes:
    shape=selected_shapes[-1]
    if shape.ShapeType() == TopAbs.TopAbs_FACE:  # TopAbs_FACE
      face=topods.Face(shape)
      selected_face.append(face)
      display.EraseAll()
      #display.selected_shapes.clear()
      display.DisplayShape(face)
      display._select_callbacks.remove(on_select)
      display._select_callbacks.append(on_select_edge)

def get_gauss_points(order):
    points, weights = np.polynomial.legendre.leggauss(order)
    return points,weights

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
  #print(f"The edge KnotVector is: {KnotVector}\n")
  surface = BRepAdaptor_Surface(face, True)
  bsurf=surface.BSpline()
  bsurf.ExchangeUV()
  uknots = bsurf.UKnots()
  UknotArr=[]
  for i in range(bsurf.NbUKnots()):
    UknotArr.append(uknots.Value(i+1))
  UknotArr_ref= np.round(np.linspace(float(UknotArr[0]),float(UknotArr[-1]),(num_xiKnots+2)),2)
  #print(f"The refined xi knot vector is : {UknotArr_ref}")
  vknots = bsurf.VKnots()
  VknotArr=[]
  for i in range(bsurf.NbVKnots()):
      VknotArr.append(vknots.Value(i+1))
  VknotArr_ref= np.round(np.linspace(float(VknotArr[0]),float(VknotArr[-1]),(num_etaKnots+2)),2)
  #print(f"The refined eta knot vector is : {VknotArr_ref}")
  ctr=0
  counter=1
  for m in range(0,len(UknotArr_ref)):
    xi = UknotArr_ref[m] 
    ctr=ctr+1
    min_distance = float('inf')
    min_distance_projector = None
    
    for n in range(0,len(VknotArr_ref)):
      if n+ctr-1<len(VknotArr_ref):
        counter+=1
        eta = VknotArr_ref[n+ctr-1]
        point2d = gp_Pnt2d(eta,xi)    # Because of reverse parametrization of geometry_curve
        projector = Geom2dAPI_ProjectPointOnCurve(point2d, OCC_bsplineCurve2d) # Orthogonal Projection
        #print(f"the point is [{xi},{eta}] and the distance is : {distance}")
        if projector is not None:
         if projector.NbPoints() !=0:
          #print(f"Num of projections for the point {ctr} and {counter} are: {projector.NbPoints()}")
          distance= round(projector.Distance(1),2)
          if distance < min_distance:
            min_distance = distance
            min_distance_projector = projector
        else:
          print(f" The projection failed for the point ({xi} {eta})")
    
    if min_distance_projector is not None:
      if min_distance_projector.NbPoints() !=0:
        sol=min_distance_projector.NbPoints()
        for p in range(1,sol+1):
          if sol==2:
            print(f"The possible number of projections are: {sol}")
          t = round(min_distance_projector.Parameter(sol),2) # The first solution among the possible projections
          curve_points_elem.append(t)

  curve_points_elem=np.array(curve_points_elem)
  merged_vector = np.sort(np.unique(np.concatenate((curve_points_elem, KnotVector))))
  #print("surface curve intersection found succesfully\n")
  #print("intersection points",curve_points_elem)
  #print("intersection points and knot vector of the curve",merged_vector)
  return merged_vector,curve_deg

def set_gps_curve(knot_vec,degree):
  knot_spans=np.unique(knot_vec)
  gps,wts=get_gauss_points(degree+1)
  a = len(gps)
  b = len(knot_spans)
  crv_gauss_mtx = np.zeros((a+1, b))
  crv_gauss_val=[]
  crv_gauss_wts=[]
  crv_gauss_jac=[]
  for i in range(len(knot_spans) - 1):
    xi_crv_1=knot_spans[i]
    xi_crv_2=knot_spans[i+1]
    for j in range(len(gps)):
        if i == 0:
            crv_gauss_mtx[j,b-1]=crv_gauss_mtx[j,b-1]+wts[j]
        crv_gauss_mtx[j,i]= crv_gauss_mtx[j,i]+gps[j]*0.5*(xi_crv_2-xi_crv_1) + 0.5*(xi_crv_1+xi_crv_2)
        crv_gauss_wts.append(wts[j])
        crv_gauss_val.append(crv_gauss_mtx[j,i])
        crv_gauss_jac.append(0.5*(xi_crv_2-xi_crv_1))
    crv_gauss_mtx[a,i] =  crv_gauss_mtx[a,i]+0.5*(xi_crv_2-xi_crv_1)
    
  return(crv_gauss_val,crv_gauss_wts,crv_gauss_jac)
            
def crv_to_srf(crv_pts,edge ,face):
  OCC_handle_curve2d, activeRangeBegin,activeRangeEnd = BRep_Tool.CurveOnSurface(topods.Edge(edge), topods.Face(face))
  OCC_bsplineCurve2d = geom2dconvert.CurveToBSplineCurve(OCC_handle_curve2d)
  xi_values=[]
  eta_values=[]
  for j in range(len(crv_pts)):
    point_on_curve = OCC_bsplineCurve2d.Value(crv_pts[j])
    xi_values.append(round(point_on_curve.X(),2))
    eta_values.append(round(point_on_curve.Y(),2))
  #print("Projected points from curve to surface")
  #print(eta_values,"\n")
  #print(xi_values,"\n")
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
  #user_refinement=8
  #refined_values = np.round(np.linspace(OCC_bsplineCurve2d.FirstParameter(), OCC_bsplineCurve2d.LastParameter(), user_refinement)[1:-1])
  #knot_insertion = TColStd_Array1OfReal(1, len(refined_values))
  #knot_insertion_mult = TColStd_Array1OfInteger(1, len(refined_values))
  #for i, value in enumerate(refined_values):
   # knot_insertion.SetValue(i + 1, value)
    #knot_insertion_mult.SetValue(i + 1, 1)
  #OCC_bsplineCurve2d.InsertKnots(knot_insertion,knot_insertion_mult)
    
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
    CPNet.append([OCC_CP.Y(), OCC_CP.X(), OCC_CPweight])
  CPNet = np.array(CPNet)
  CPNet_reshaped = CPNet.reshape(-1, 3)
  #print("CPNet looks like ", CPNet_reshaped)
  return pDegree,uKnotVector,CPNet_reshaped

def process_intersection(display):
    matrix_dict_dirichlet=[]
    print(f" Number of faces and edges selected are: {len(selected_face)} and {len(selected_edge)}")
    for con_num in range(len(selected_edge)):
        crv_1_element_pts,deg=ed_fc_inter(selected_edge[con_num],selected_face[con_num])
        display.EraseAll()
        display.DisplayMessage(gp_Pnt(0,0,0),"Process completed successfully!",50, (0,0,1),False) 
        display.DisplayMessage(gp_Pnt(20,-2,0),"You can Exit Now",30, (1,0,0),False)
        crv_gs_1,crv_gs_wt_1,crv_gs_jac_1= set_gps_curve(crv_1_element_pts,deg) # Get Gauss Points  on the integration curve
        xi_gs_1,eta_gs_1= crv_to_srf(crv_gs_1,selected_edge[con_num],selected_face[con_num])# Transfering integration points to the associated surface
        crv_deg,crv_knt_vec,crv_CPs= nurbs_detail_coupling_curve(selected_edge[con_num],selected_face[con_num])
        ids=get_face_ids(selected_face[con_num])
        matrix_dict = {'cur_pts': crv_gs_1,'gs_wts':crv_gs_wt_1,'gs_jac':crv_gs_jac_1, 'xi_vls_1': xi_gs_1, 'eta_vls_1': eta_gs_1,
                      'crv_deg':crv_deg,'crv_knotVec':crv_knt_vec,'crv_CPS':crv_CPs,'face_id':ids}
        for key, value in matrix_dict.items():
          print(f"{key}: {value}")
        matrix_dict_dirichlet.append(matrix_dict)
      
    data_to_save = {'dir_list': matrix_dict_dirichlet}
    savemat(output_file, data_to_save)
   
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
    global num_xiKnots, num_etaKnots
    global matrix_dict, output_file
    global selected_edge, selected_face
    selected_edge=[]
    selected_face=[]
    global iges_file_path
    #iges_file_path="C:\\Users\\oslon\\Desktop\\MSc\\Thesis\\Code Scripts\\low_dim_rect_3patch.igs"
    #iges_file_path="C:\\Users\\oslon\\Desktop\\MSc\\Thesis\\Code Scripts\\basicRectSurf.igs"
    #iges_file_path="C:\\Users\\oslon\\Desktop\\MSc\\Thesis\\Code Scripts\\trimmedRect.igs"

    iges_file_path = sys.argv[1]
    output_file= 'b_rep.mat'
    #operation= sys.argv[2]
    num_xiKnots=int(sys.argv[3])
    num_etaKnots=int(sys.argv[4])
    #num_xiKnots=10
    #num_etaKnots=10
    load_iges(display, iges_file_path)
    display._select_callbacks.append(on_select)
    start_display()
    