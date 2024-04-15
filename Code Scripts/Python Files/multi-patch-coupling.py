import sys
import numpy as np
from OCC.Core.BRepBuilderAPI import BRepBuilderAPI_NurbsConvert
from OCC.Core.BRepAdaptor import BRepAdaptor_Curve,BRepAdaptor_Surface
from OCC.Extend.DataExchange import read_iges_file
from OCC.Extend.TopologyUtils import TopologyExplorer
from OCC.Core import TopAbs
from OCC.Core.BRep import BRep_Tool
from OCC.Core.Geom2dConvert import geom2dconvert  
from OCC.Core.TColgp import TColgp_Array1OfPnt2d
from OCC.Core.TColStd import TColStd_Array1OfInteger, TColStd_Array1OfReal
from OCC.Core.TopoDS import topods
from OCC.Core.gp import gp_Pnt2d,gp_Pnt
from OCC.Core.Geom2dAPI import Geom2dAPI_ProjectPointOnCurve
from OCC.Core.GeomAPI import GeomAPI_ProjectPointOnSurf
from scipy.io import savemat


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
  merged_vector = np.sort(np.unique(np.concatenate((curve_points_elem, KnotVector))))
  #print("surface curve intersection found succesfully\n")
  #print("intersection points",curve_points_elem)
  #print("intersection points and knot vector of the curve",merged_vector)
  return merged_vector,curve_deg

def edge_refinement(edge,face):
  
  OCC_handle_curve2d, activeRangeBegin, activeRangeEnd = BRep_Tool.CurveOnSurface(topods.Edge(edge), topods.Face(face))
  OCC_bsplineCurve2d = geom2dconvert.CurveToBSplineCurve(OCC_handle_curve2d)
  if OCC_bsplineCurve2d.IsPeriodic(): 
    print("Found periodic curve!")
    OCC_bsplineCurve2d.SetNotPeriodic()
  # pDegree
  pDegree = OCC_bsplineCurve2d.Degree()
  # Knots Insertion
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
    #print(f"the surface point is [{xi_vec[j]},{eta_vec[j]}] and the physical point is [{x},{y},{z}]")
  #print("Projected surface to physical successfully\n")

  return x_vec,y_vec,z_vec
      
def phy_to_surf(x_vec,y_vec,z_vec,face):
  #print("In physical to surface")   
  #print(len(x_vec))
  xi_vec=[]
  eta_vec=[]
  for j in range(len(x_vec)):
    point_in_physical_space = gp_Pnt(x_vec[j], y_vec[j], z_vec[j])
    face_surface = BRep_Tool.Surface(face)
    projector = GeomAPI_ProjectPointOnSurf(point_in_physical_space, face_surface)
    point_on_surface = projector.NearestPoint()
    u, v = projector.LowerDistanceParameters()
    #print(f"Surface Parameters of physical point [{x_vec[j]}{y_vec[j]}]: ({v}, {u})")
    xi_vec.append(round(u,2))
    eta_vec.append(round(v,2))
  return eta_vec,xi_vec
       
def surf_to_curve(xi_vec,eta_vec,edge,face):
  #print("In surface to curve")
  OCC_handle_curve2d, activeRangeBegin,activeRangeEnd = BRep_Tool.CurveOnSurface(topods.Edge(edge), topods.Face(face))
  OCC_bsplineCurve2d = geom2dconvert.CurveToBSplineCurve(OCC_handle_curve2d)
  crv_vec=[]
  for j in range(len(xi_vec)):
    point2d = gp_Pnt2d(eta_vec[j],xi_vec[j]) 
    projector = Geom2dAPI_ProjectPointOnCurve(point2d, OCC_bsplineCurve2d)
    t = round(projector.Parameter(1),2)
    crv_vec.append(t)
  #print("the projected curve points are: ")
  #print(crv_vec)
  return crv_vec

def gauss_points(order):
    points, weights = np.polynomial.legendre.leggauss(order+1)
    return points,weights

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

def set_gps_curve(knot_vec,degree):
  knot_spans=np.unique(knot_vec)
  gps,wts=gauss_points(degree+1)
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
    crv_gauss_mtx[a,i] =  crv_gauss_mtx[a,i]+0.5*(xi_crv_2-xi_crv_1) # NOT USED ANYWHERE
    
  return(crv_gauss_val,crv_gauss_wts,crv_gauss_jac)

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
    CPNet.extend([OCC_CP.Y(), OCC_CP.X(), OCC_CPweight])
  CPNet = np.array(CPNet)
  CPNet_reshaped = CPNet.reshape(-1, 3)
  #print("CPNet looks like ", CPNet_reshaped)
  return pDegree,uKnotVector,CPNet_reshaped

def process_intersection(coupling_info):
    
    matrix_dict_coupling=[]
    num_pairs= len(coupling_info)
    for pair in range(num_pairs):
        couple_data=coupling_info[pair] # Accessing each tuple
        face_I_ind=couple_data[4]
        face_J_ind=couple_data[1]
        selected_edges=[couple_data[5], couple_data[2] ]
        selected_faces=[couple_data[3], couple_data[0] ]
        crv_pts_1,deg_1=ed_fc_inter(selected_edges[0],selected_faces[0]) # intersection of first face and edge
        crv_pts_2,deg_2=ed_fc_inter(selected_edges[1],selected_faces[1]) # intersection of second face and edge
        xi_values_2,eta_values_2= crv_to_srf(crv_pts_2,selected_edges[1],selected_faces[1]) # curve_2 points to surface_2 points
        x_inter,y_inter,z_inter= surf_to_phy(xi_values_2,eta_values_2,selected_faces[1]) # surface_2 points to physical points
        xi_other,eta_other = phy_to_surf(x_inter,y_inter,z_inter,selected_faces[0]) # physical points back to surface_1 points
        transfer_pts_1= surf_to_curve(xi_other,eta_other,selected_edges[0],selected_faces[0]) # surface_1 points back to curve_1 points
        crv_1_element_pts = np.sort(np.unique(np.concatenate((np.array(crv_pts_1), np.array(transfer_pts_1)))))
        crv_gs_1,crv_gs_wt_1,crv_gs_jac_1= set_gps_curve(crv_1_element_pts,deg_1) # Get Gauss Points  on the integration curve
        xi_gs_1,eta_gs_1= crv_to_srf(crv_gs_1,selected_edges[0],selected_faces[0])# Transfering integration points of first surface
        x_gs_1,y_gs_1,z_gs_1= surf_to_phy(xi_gs_1,eta_gs_1,selected_faces[0]) # surface_1 Gauss points to physical points
        xi_gs_2,eta_gs_2 = phy_to_surf(x_gs_1,y_gs_1,z_gs_1,selected_faces[1]) # physical gauss points on surface_2 
        crv_deg,crv_knt_vec,crv_CPs= nurbs_detail_coupling_curve(selected_edges[0],selected_faces[0])
        matrix_dict = {'cur_pts_1': crv_gs_1,'gs_wts':crv_gs_wt_1,'gs_jac':crv_gs_jac_1, 'xi_vls_1': xi_gs_1,
                        'eta_vls_1': eta_gs_1,'xi_vls_2': xi_gs_2, 'eta_vls_2': eta_gs_2, 'crv_deg':crv_deg,'crv_knotVec':crv_knt_vec,
                        'crv_CPS':crv_CPs,'face_I_ind':face_I_ind,'face_J_ind':face_J_ind}
        for key, value in matrix_dict.items():
            print(f"{key}: {value}")
        matrix_dict_coupling.append(matrix_dict)
    data_to_save = {'coupling_list': matrix_dict_coupling}
    savemat(output_file, data_to_save)

def getFaces(iges_file_path):
    base_shape = read_iges_file(iges_file_path)
    nurbs_converter = BRepBuilderAPI_NurbsConvert(base_shape, True)
    converted_shape = nurbs_converter.Shape()
    expl = TopologyExplorer(converted_shape)
    fc_idx = 1
    face_list=[]
    for face in expl.faces():
        face_list.append(topods.Face(face))
        fc_idx += 1
    return face_list

def find_common_edges(edge_list1,edge_list2):
    checker=0
    tol=1e-3
    for i in range(len(edge_list1)):
        edge_I=topods.Edge(edge_list1[i])
        vert_exp_1=TopologyExplorer(edge_I)
        vertex_points_1=[]
        #print(f"For edge {i} of face 1, mid point is \n")
        for vertex in vert_exp_1.vertices():
            point1= BRep_Tool.Pnt(vertex)
            vertex_points_1.append(point1)
        mid_pt_1= gp_Pnt((vertex_points_1[0].X()+vertex_points_1[1].X())/2, 
                         (vertex_points_1[0].Y()+vertex_points_1[1].Y())/2,
                         (vertex_points_1[0].Z()+vertex_points_1[1].Z())/2)
        #print(f"{mid_pt_1.X()} and {mid_pt_1.Y()} and {mid_pt_1.Z()}\n")

        for j in range(len(edge_list2)):
            edge_J=topods.Edge(edge_list2[j])
            vert_exp_2=TopologyExplorer(edge_J)
            vertex_points_2=[]
            #print(f"For edge {j} of face 2, mid point is \n")
            for vertex in vert_exp_2.vertices():
                point2= BRep_Tool.Pnt(vertex)
                vertex_points_2.append(point2)
            mid_pt_2= gp_Pnt((vertex_points_2[0].X()+vertex_points_2[1].X())/2, 
                         (vertex_points_2[0].Y()+vertex_points_2[1].Y())/2,
                         (vertex_points_2[0].Z()+vertex_points_2[1].Z())/2)
            #print(f"{mid_pt_2.X()} and {mid_pt_2.Y()} and {mid_pt_2.Z()}\n")

            if (mid_pt_1.Distance(mid_pt_2) < tol) :
                checker=1
                return True, i, j
    if checker==0:
        return False,0,0
        
def findAdjacentFacePairs(face_list):
    coupling_data=[]

    for i in range(len(face_list)):
        face_I=face_list[i]
        edge_list_1=[]
        edge_expl_1= TopologyExplorer(face_I)
        for edge in edge_expl_1.edges():
            edge_list_1.append(edge)
        for j in range(i+1,len(face_list)):
            face_J= face_list[j]
            edge_list_2=[]
            edge_expl_2= TopologyExplorer(face_J)
            for edge in edge_expl_2.edges():
                edge_list_2.append(edge)
            
            res,ind_edge_face_I,ind_edge_face_J= find_common_edges(edge_list_1,edge_list_2)
            if res:
                print(f" The face {i+1} is adjacent face {j+1}")
                coupling_data.append((face_I,i+1, edge_list_1[ind_edge_face_I],face_J,j+1,edge_list_2[ind_edge_face_J]))
            else:
                print(f" The face {i+1} is not adjacent to the face {j+1}")

    return coupling_data


if __name__ == "__main__":

    #iges_file_path="C:\\Users\\oslon\\Desktop\\MSc\\Thesis\\Code Scripts\\low_dim_rect_3patch.igs"
    iges_file_path = sys.argv[1]
    global num_xiKnots, num_etaKnots
    global output_file,operation 
    output_file= 'coup.mat'
    #operation= sys.argv[2]
    num_xiKnots=int(sys.argv[3])
    num_etaKnots=int(sys.argv[4])
    #num_xiKnots=10
    #num_etaKnots=10
    faces = getFaces(iges_file_path)
    coupling_data= findAdjacentFacePairs(faces)
    process_intersection(coupling_data)
    
    
    