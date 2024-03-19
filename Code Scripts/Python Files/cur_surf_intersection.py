from OCC.Core.BRepBuilderAPI import BRepBuilderAPI_MakeEdge, BRepBuilderAPI_MakeFace
from OCC.Core.GeomAPI import GeomAPI_Interpolate, GeomAPI_IntSS
from OCC.Core.GeomAbs import GeomAbs_C0, GeomAbs_C1, GeomAbs_C2, GeomAbs_C3
from OCC.Core.TColgp import TColgp_Array1OfPnt
from OCC.Core.Geom import Geom_BSplineCurve, Geom_BSplineSurface
from OCC.Core.TopoDS import topods
from OCC.Display.SimpleGui import init_display

# Create a NURBS curve
curve_points = [(-1, 0, 0), (0, 1, 0), (1, 0, 0)]
curve_degree = 2

curve = BRepBuilderAPI_MakeEdge(Geom_BSplineCurve(curve_points, None, None, curve_degree).Curve()).Edge()

# Create a NURBS surface
surface_points = [[(-1, -1, 0), (0, -1, 1), (1, -1, 0)],
                  [(-1, 0, 1), (0, 0, 2), (1, 0, 1)],
                  [(-1, 1, 0), (0, 1, 1), (1, 1, 0)]]
surface_degree_u = 2
surface_degree_v = 2

surface = BRepBuilderAPI_MakeFace(Geom_BSplineSurface(surface_points, None, None, None, surface_degree_u, surface_degree_v).Surface()).Face()

# Display the curve and surface
display, _, _, _ = init_display()
display.DisplayShape(curve)
display.DisplayShape(surface)

# Find intersection
intersector = GeomAPI_IntSS(topods_Edge=curve, topods_Face=surface, theTolerance=1e-6)
intersector.Perform()
if intersector.IsDone():
    for i in range(1, intersector.NbLines() + 1):
        line = intersector.Line(i)
        display.DisplayShape(line)
else:
    print("No intersection found.")

display.FitAll()
display.View_Top()
display.View_Iso()
display.FitAll()

display.Start()
