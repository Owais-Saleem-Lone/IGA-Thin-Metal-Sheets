from OCC.Extend.DataExchange import read_iges_file
from OCC.Extend.ShapeFactory import make_face
from OCC.Extend.TopologyUtils import TopologyExplorer
from OCC.Core.BRepAdaptor import BRepAdaptor_Surface,BRepAdaptor_Curve
from OCC.Core.TopoDS import topods_Edge
from OCC.Core.TopExp import TopExp_Explorer
from OCC.Core.TopAbs import TopAbs_EDGE
from OCC.Core.TopLoc import TopLoc_Location
from OCC.Core.BRepBuilderAPI import BRepBuilderAPI_MakeEdge
from OCC.Core.GeomConvert import GeomConvert_CompCurveToBSplineCurve
from OCC.Core.TopoDS import topods_Edge
from OCC.Core.BRepBuilderAPI import BRepBuilderAPI_MakeWire
from OCC.Core.GeomConvert import GeomConvert_CompCurveToBSplineCurve
from OCC.Core.TopoDS import topods_Wire
from OCC.Core.GeomAbs import GeomAbs_C1
# Load the IGES file
iges_path = "C:\\Users\\oslon\\Desktop\\MSc\\Thesis\\Code Scripts\\trimmedRect.igs"   # Replace with the path to your IGES file
compound = read_iges_file(iges_path)

# Iterate over faces in the compound
for face in TopologyExplorer(compound).faces():
    # Check if the face is a trimmed surface
    surface_adaptor = BRepAdaptor_Surface(face)
    if surface_adaptor.GetType() == 6:  # Geom_BSplineSurface
        print("This face represents a B-spline surface.")

        # Iterate over edges of the face
        edge_explorer = TopExp_Explorer(face, TopAbs_EDGE)
        i=1
        wire_builder = BRepBuilderAPI_MakeWire()
        while edge_explorer.More():
            edge = topods_Edge(edge_explorer.Current())
            
            # Access the curve from the edge
            edge_curve = BRepAdaptor_Curve(edge)
            print("Curve Number\n",i)
            print("Curve Degree is:", edge_curve.Degree())
            if edge_curve.Degree()==2:
                if i==1:
                    wire_builder = BRepBuilderAPI_MakeWire(edge_curve)
                else:
                    wire_builder.Add(edge)

            # Now you can work with the curve (e.g., print its type or extract more information)
            print(f"Curve Type: {edge_curve.GetType()}")
            
            # Move to the next edge
            edge_explorer.Next()
            
            i+=1
        composite_wire = wire_builder.Wire()

        bspline_curve_converter = GeomConvert_CompCurveToBSplineCurve(composite_wire)
        nurbs_curve = bspline_curve_converter.Curve()

        # Convert the composite curve to a B-spline curve
        bspline_curve_converter = GeomConvert_CompCurveToBSplineCurve(composite_wire, GeomAbs_C1)
        nurbs_curve = bspline_curve_converter.Curve()
        
    else:
        print("This face does not represent a B-spline surface.")




# from OCC.Extend.DataExchange import read_iges_file
# from OCC.Extend.TopologyUtils import TopologyExplorer
# from OCC.Core.BRepAdaptor import BRepAdaptor_Surface
# from OCC.Core.GeomAdaptor import GeomAdaptor_Surface
# from OCC.Core.BRepBuilderAPI import BRepBuilderAPI_NurbsConvert
# # Load the IGES file
# iges_path = "C:\\Users\\oslon\\Desktop\\MSc\\Thesis\\Code Scripts\\trimmedRect.igs" 
# compound = read_iges_file(iges_path)
# nurbs_converter = BRepBuilderAPI_NurbsConvert(compound, True)
# converted_shape = nurbs_converter.Shape()
# # Iterate over faces in the compound
# for face in TopologyExplorer(converted_shape).faces():
#     # Check if the face is a trimmed surface
#     surface_adaptor = BRepAdaptor_Surface(face)
#     print(surface_adaptor.GetType())
#     if surface_adaptor.GetType() == 4:  # Type 4 corresponds to Geom_RectangularTrimmedSurface
#         print("This face represents a trimmed surface.")

#         # Access trimming parameters and ranges
#         u_min, u_max, v_min, v_max = surface_adaptor.FirstUParameter(), surface_adaptor.LastUParameter(), \
#                                      surface_adaptor.FirstVParameter(), surface_adaptor.LastVParameter()
#         print(f"U Range: {u_min} to {u_max}")
#         print(f"V Range: {v_min} to {v_max}")

#         # Access trimming curves
#         num_curves = surface_adaptor.NbUKnots() + surface_adaptor.NbVKnots()
#         print(f"Number of trimming curves: {num_curves}")

#     else:
#         print("This face is not a trimmed surface.")






























# from OCC.Core.gp import gp_Pnt2d
# from OCC.Core.IGESControl import IGESControl_Reader
# from OCC.Core.TopoDS import topods_Face
# from OCC.Core.TColStd import TColStd_Array1OfReal, TColStd_Array1OfInteger
# from OCC.Core.TColgp import TColgp_Array1OfPnt2d
# from OCC.Extend.DataExchange import read_iges_file
# from OCC.Extend.ShapeFactory import make_face
# from OCC.Core.IGESData import IGESData_IGESEntity,IGESGeom_ToolTrimmedSurface
# #from OCC.Core.IGESGeom import IGESGeom_ToolTrimmedSurface

# # Load the IGES file
# iges_path = "C:\\Users\\oslon\\Desktop\\MSc\\Thesis\\Code Scripts\\trimmedRect.igs"
# iges_reader = IGESControl_Reader()
# iges_reader.ReadFile(iges_path)

# # Transfer the entities from the IGES file to a TopoDS_Compound
# compound = iges_reader.TransferRoots().FirstRoot()
# face = topods_Face(compound)

# # Create an IGESGeom_ToolTrimmedSurface for the face
# tool_trimmed_surface = IGESGeom_ToolTrimmedSurface(IGESData_IGESEntity(face))

# # Check if the surface is a trimmed surface
# if tool_trimmed_surface.TrimmedSurface():
#     # Retrieve the underlying surface
#     underlying_surface = tool_trimmed_surface.Surface()

#     # Retrieve trimming parameters and ranges
#     parameters = TColStd_Array1OfReal(1, 2)
#     ranges = TColStd_Array1OfReal(1, 2)
#     tool_trimmed_surface.SurfaceTrimmingParameters(parameters, ranges)

#     print("Trimming Parameters:")
#     print(f"U Range: {parameters.Value(1)} to {parameters.Value(2)}")
#     print(f"V Range: {ranges.Value(1)} to {ranges.Value(2)}")

#     # Retrieve trimming points and types
#     trim_points = TColgp_Array1OfPnt2d(1, 100)
#     trim_types = TColStd_Array1OfInteger(1, 100)
#     trim_indices = TColStd_Array1OfInteger(1, 100)
#     if tool_trimmed_surface.SurfaceTrims(trim_points, trim_types, trim_indices):
#         print("\nTrimming Curves:")
#         for i in range(1, trim_points.Length() + 1):
#             if trim_types.Value(i) == 1:  # 1 represents the type of curve (1: point)
#                 print(f"Point at (U, V): ({trim_points.Value(i).X()}, {trim_points.Value(i).Y()})")
#             elif trim_types.Value(i) == 2:  # 2 represents the type of curve (2: circle, arc)
#                 print(f"Curve of type 2 (e.g., circle or arc) with index: {trim_indices.Value(i)}")

#     # Optionally, convert the IGES trimmed surface to an OCC geometric surface
#     geometric_surface = tool_trimmed_surface.ConvertSurface()
# else:
#     print("The surface is not a trimmed surface.")
