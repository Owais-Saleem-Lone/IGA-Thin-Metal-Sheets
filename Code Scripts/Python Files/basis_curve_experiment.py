import numpy as np
from OCC.Core.gp import gp_Pnt2d,gp_Pnt
from OCC.Core.Adaptor3d import Adaptor3d_Surface
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
from OCC.Core.BRepPrimAPI import BRepPrimAPI_MakeCylinder
from OCC.Core.Adaptor3d import Adaptor3d_Surface
from OCC.Core.Geom2dAPI import Geom2dAPI_InterCurveCurve
from OCC.Core.gp import gp_Pnt2d
from OCC.Core.TopoDS import TopoDS_Face
from OCC.Core.ShapeConstruct import ShapeConstruct_ProjectCurveOnSurface
from OCC.Core.ShapeAnalysis import ShapeAnalysis_Surface,ShapeAnalysis_Edge
from OCC.Core.Geom import Geom_Curve
from OCC.Core.TColgp import TColgp_Array1OfPnt
from OCC.Core.Geom2d import Geom2d_Curve
from OCC.Core.Standard import Standard_HandlerStatus
from OCC.Core import Geom2dAPI


def load_iges(display, iges_file_path):

    base_shape = read_iges_file(iges_file_path)
    nurbs_converter = BRepBuilderAPI_NurbsConvert(base_shape, True)
    basic_shape = nurbs_converter.Shape()
    expl = TopologyExplorer(basic_shape)
    
    for face in expl.faces():
        my_surface=BRep_Tool.Surface(face)
        #new_surface=BRepAdaptor_Surface(face).Surface()
        surface_analysis = ShapeAnalysis_Surface(my_surface)
        u_param = 0.5  
        u_iso_curve = surface_analysis.UIso(u_param)
        edge=BRepBuilderAPI_MakeEdge(u_iso_curve).Edge()
        curve_wire=BRepBuilderAPI_MakeWire(edge).Wire()
        con_cur=BRepAdaptor_Curve(edge).BSpline()
        control_points = TColgp_Array1OfPnt(1, con_cur.NbPoles())
        for i in range(con_cur.NbPoles()):
            control_point = con_cur.Pole(i + 1)  # Pole indices start from 1
            control_points.SetValue(i + 1, control_point)
            x = control_point.X()
            y = control_point.Y()
            z = control_point.Z()
            print(f"Control point {i}: ({x}, {y}, {z})")
        tr_fa=ShapeAnalysis_Edge ().HasPCurve(edge,face)
        print(f"so the result is {tr_fa}")
        new_face=BRepBuilderAPI_MakeFace(topods.Face(face),topods.Wire(curve_wire)).Face()
        projector = ShapeConstruct_ProjectCurveOnSurface()
        projector = ShapeConstruct_ProjectCurveOnSurface()
        projector.SetSurface(surface_analysis)
        projector.SetPrecision(0.001)  # Set precision if needed
        result_curve_2d: Geom2d_Curve
        #Perform the projection
        result_curve_2d=None
        success = projector.Perform( u_iso_curve,u_iso_curve.FirstParameter(),u_iso_curve.LastParameter(),result_curve_2d,-1,-1)
        if success:
            print("Yayy, alles gut")
        else:
            print("Shut Up")
        print(f"The parameters of projected curve are {result_curve_2d.FirstParameter()} and {result_curve_2d.LastParameter()}")
        OCC_bsplineCurve2d = geom2dconvert.CurveToBSplineCurve(result_curve_2d)
        uNoCPs = OCC_bsplineCurve2d.NbPoles()
        OCC_CPnet = TColgp_Array1OfPnt2d(1,uNoCPs)
        OCC_bsplineCurve2d.Poles(OCC_CPnet)
        for i in range(1, uNoCPs + 1):
        # Extract the pole
            pole = OCC_CPnet.Value(i)
            # Extract U and V coordinates
            u = pole.X()
            v = pole.Y()
            # Print the coordinates
            print(f"Pole {i}: ({u}, {v})")
          
    display.DisplayShape(new_face, update=True)
    #display.DisplayShape(basic_shape, update=True)
    #display.DisplayShape(u_iso_curve,update=True)
   
if __name__ == '__main__':
    display, start_display, add_menu, add_function_to_menu = init_display('wx')
    iges_file_path="C:\\Users\\oslon\\Desktop\\MSc\\Thesis\\Code Scripts\\mod_rect.igs"
    global num_xiKnots, num_etaKnots
    global matrix_dict,output_file,operation
    num_xiKnots=4
    num_etaKnots=6
    print(f"xi refinement points are: {num_xiKnots}")
    print(f"eta refinement points are <with increment of course>: {num_etaKnots}")
    load_iges(display, iges_file_path)
    start_display()