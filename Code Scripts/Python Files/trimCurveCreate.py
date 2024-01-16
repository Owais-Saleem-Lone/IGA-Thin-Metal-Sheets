from OCC.Core.gp import gp_Pnt, gp_Pnt2d
from OCC.Core.GeomAPI import GeomAPI_ProjectPointOnSurf
from OCC.Core.BRep import BRep_Tool
from OCC.Core.TopExp import TopExp_Explorer
from OCC.Core.TopoDS import topods_Vertex
from OCC.Extend.TopologyUtils import TopologyExplorer
from OCC.Extend.DataExchange import read_iges_file
from OCC.Core.TopExp import TopExp_Explorer
from OCC.Core.TopAbs import TopAbs_VERTEX
from OCC.Core.Extrema import Extrema_ExtAlgo_Grad
from OCC.Core.GeomAPI import Geom2dAPI_ProjectPointOnCurve

# Step 1: Load the IGES file
iges_file_path = "C:\\Users\\oslon\\Desktop\\MSc\\Thesis\\Code Scripts\\trimmedRect.igs"
shape = read_iges_file(iges_file_path)

# Step 2: Extract the TopoDS_Face and TopoDS_Edge
explorer = TopologyExplorer(shape)
rectangular_face = None
trimming_curve = None

for face in explorer.faces():
    rectangular_face = face

for edge in explorer.edges():
    trimming_curve = edge

# Step 3: Identify the parameterization of the rectangular face
if rectangular_face:
    face_surface = BRep_Tool.Surface(rectangular_face)

# Step 4: Project the vertices of the trimming curve onto the surface
if trimming_curve and rectangular_face:
    vertex_exp = TopExp_Explorer(trimming_curve, TopAbs_VERTEX)
    trimmed_curve_vertices = []
    while vertex_exp.More():
        vertex = topods_Vertex(vertex_exp.Current())
        vertex_point = BRep_Tool.Pnt(vertex)
        # Project the vertex onto the surface
        project_point = GeomAPI_ProjectPointOnSurf(vertex_point, face_surface,Extrema_ExtAlgo_Grad).Point()
        u, v = BRep_Tool.Parameters(vertex, rectangular_face)
        trimmed_curve_vertices.append(gp_Pnt2d(u, v))
        vertex_exp.Next()

# Step 5: Create a trimmed curve in the parametric space of the rectangular surface
if len(trimmed_curve_vertices) > 1:
    trimming_curve_params = [gp_Pnt2d(u, v) for u, v in trimmed_curve_vertices]
    trimming_curve_2d = Geom2dAPI_ProjectPointOnCurve(trimming_curve_params, trimming_curve).Curve2d()
    # Now, trimming_curve_2d represents the trimmed curve in the parametric space of the rectangular surface
    # Extract details of the trimmed curve
    if trimming_curve_2d.IsKind("Geom_BSplineCurve"):
        bspline_curve_2d = trimming_curve_2d.BSpline()
        # Print knot vectors
        print("Knots U:", bspline_curve_2d.Knots())
        print("Multiplicity U:", bspline_curve_2d.Multiplicities())

        # Print control points
        control_points = [bspline_curve_2d.Pole(i) for i in range(bspline_curve_2d.NbPoles())]
        print("Control Points:", control_points)














## THIS PART CREATES A RECTANGULAR SURFACE FROM POINTS
# from OCC.Core.gp import gp_Pnt, gp_Dir, gp_Ax2, gp_Lin, gp_Ax2d, gp_Pnt2d, gp_Dir2d,gp_XYZ
# from OCC.Core.Geom2d import Geom2d_Line
# from OCC.Core.BRepBuilderAPI import BRepBuilderAPI_MakeEdge, BRepBuilderAPI_MakeWire, BRepBuilderAPI_MakeFace,BRepBuilderAPI_MakeVertex
# from OCC.Core.TopoDS import topods_Vertex,TopoDS_Edge
# from OCC.Display.SimpleGui import init_display

# # Dimensions of the rectangle
# width = 10.0
# height = 5.0
# corner_point1 = gp_Pnt(0, 0,0)
# corner_point2 = gp_Pnt(width, 0,0)
# corner_point3 = gp_Pnt(width, height,0)
# corner_point4 = gp_Pnt(0, height,0)

# vertex_builder = BRepBuilderAPI_MakeVertex(corner_point1)
# vertex1 = vertex_builder.Vertex()
# vertex_builder = BRepBuilderAPI_MakeVertex(corner_point2)
# vertex2 = vertex_builder.Vertex()
# vertex_builder = BRepBuilderAPI_MakeVertex(corner_point3)
# vertex3 = vertex_builder.Vertex()
# vertex_builder = BRepBuilderAPI_MakeVertex(corner_point4)
# vertex4 = vertex_builder.Vertex()

# line1=BRepBuilderAPI_MakeEdge(vertex1,vertex2)
# line2=BRepBuilderAPI_MakeEdge(vertex2,vertex3)
# line3=BRepBuilderAPI_MakeEdge(vertex3,vertex4)
# line4=BRepBuilderAPI_MakeEdge(vertex4,vertex1)

# # Create a wire from the edges
# line1= line1.Edge()
# line2= line2.Edge()
# line3= line3.Edge()
# line4= line4.Edge()

# wire_builder = BRepBuilderAPI_MakeWire()
# wire_builder.Add(line1)
# wire_builder.Add(line2)
# wire_builder.Add(line3)
# wire_builder.Add(line4)
# wire= wire_builder.Wire()
# # Create a face from the wire
# rectangle_face = BRepBuilderAPI_MakeFace(wire, True).Face()

# # Display the 2D rectangular surface
# display, start_display, add_menu, add_function_to_menu = init_display()
# display.DisplayShape(rectangle_face, update=True)
# start_display()

















## THIS PART OF THE CODE CREATES A TRIMMING CURVE FROM A CIRCLE

# from OCC.Core.gp import gp_Pnt2d,gp_Circ2d,gp_Ax2d
# from OCC.Core.Geom2d import Geom2d_Circle, Geom2d_TrimmedCurve
# from OCC.Display.SimpleGui import init_display
# import math

# #Create a base 2D circle

# center_point = gp_Pnt2d(0, 0)
# radius = 5
# circle = Geom2d_Circle(gp_Circ2d(gp_Ax2d(), radius, True))

# #Create a trimmed version of the circle

# start_param = circle.FirstParameter()
# end_param = 2
# trimmed_curve = Geom2d_TrimmedCurve(circle, start_param, end_param)

# #Display the original circle and the trimmed curve

# display, start_display, add_menu, add_function_to_menu = init_display()
# display.DisplayShape(circle, update=True)
# display.DisplayShape(trimmed_curve, update=True)
# start_display()
