##Copyright 2020 Thomas Paviot (tpaviot@gmail.com)
from OCC.Core.BRepPrimAPI import BRepPrimAPI_MakeTorus
from OCC.Core.BRepBuilderAPI import BRepBuilderAPI_NurbsConvert
from OCC.Core.BRepAdaptor import BRepAdaptor_Curve
from OCC.Core.IGESControl import IGESControl_Reader
from OCC.Extend.DataExchange import read_iges_file
from OCC.Extend.TopologyUtils import TopologyExplorer
from OCC.Core.GeomAbs import GeomAbs_BSplineSurface
from OCC.Core.GeomAbs import GeomAbs_BSplineCurve
from OCC.Core import TopoDS
from OCC.Core import TopAbs
iges_file_path = "C:\\Users\\oslon\\Desktop\\MSc\\Thesis\\Code Scripts\\multi-patch plate.igs"
base_shape = read_iges_file(iges_file_path)
# conversion to a nurbs representation


# TRIAL TRIAL TRIAL
iges_reader = IGESControl_Reader()
iges_reader.ReadFile(iges_file_path)
iges_shape = TopoDS.TopoDS_Shape()
iges_reader.TransferRoots()
print("Number of faces in the file are : ",iges_reader.NbShapes())
subshape=iges_reader.Shape(1)
# TRIAL TRIAL TRIAL

nurbs_converter = BRepBuilderAPI_NurbsConvert(base_shape, True)
# nurbs_converter.Perform()
converted_shape = nurbs_converter.Shape()
print("Trying to get the type of curve",subshape.ShapeType()==TopAbs.TopAbs_EDGE)
print("Trying to get the type of plate",subshape.ShapeType()==TopAbs.TopAbs_FACE)


exit()
expl = TopologyExplorer(converted_shape)

# loop over faces
fc_idx = 1

for face in expl.edges():
    print("=== Face %i ===" % fc_idx)
    surf = BRepAdaptor_Curve(face)
    surf_type = surf.GetType()
    # check each of the is a BSpline surface
    # it should be, since we used the nurbs converter before
    if not surf_type == GeomAbs_BSplineCurve:
        raise AssertionError("the face was not converted to a GeomAbs_BSplineSurface")
    # get the nurbs
    bcurve = surf.BSpline()
    print("UDegree:", bcurve.Degree())
  
   #knots array
    vknots = bcurve.Knots()
    print("\nknots:", end="")
    for i in range(bcurve.NbKnots()):
        print(vknots.Value(i + 1), end=" ")
    print("\n")
    # weights, a 2d array
    weights = bcurve.Weights()
    # weights can be None
    if weights is not None:
        print("Weights:", end="")
        for i in range(bcurve.NbPoles()):
            print(weights.Value(i + 1), end=" ")
    # control points (aka poles), as a 2d array
    poles = bcurve.Poles()
    # weights can be None
    if poles is not None:
        print("Poles (control points):", end="")
        for i in range(bcurve.NbPoles()):
                p = poles.Value(i + 1)
                print(p.X(), p.Y(), p.Z(), end=" ")
    print()
    fc_idx += 1
    num_points = 100
    points_on_curve = [bcurve.Value(i / (num_points - 1)) for i in range(num_points)]

    # Display the points
    for point in points_on_curve:
        print(point.Coord())
