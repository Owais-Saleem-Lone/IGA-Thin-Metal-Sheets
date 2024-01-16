from OCC.Core.BRepBuilderAPI import BRepBuilderAPI_MakeWire
from OCC.Core.BRepOffsetAPI import BRepOffsetAPI_ThruSections
from OCC.Core.BRepPrimAPI import BRepPrimAPI_MakeBox
from OCC.Core.IGESControl import IGESControl_Reader
from OCC.Core.TopoDS import topods_Face, topods_Wire, TopoDS_Compound
from OCC.Core.gp import gp_Pnt
from OCC.Extend.DataExchange import read_iges_file
from OCC.Extend.TopologyUtils import TopologyExplorer
from OCC.Extend.ShapeFactory import make_edge, make_wire
from OCC.Core.GeomAPI import GeomAPI_ProjectPointOnSurf
from OCC.Core.GeomAdaptor import GeomAdaptor_Curve
from OCC.Core.TColgp import TColgp_Array1OfPnt
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from OCC.Core.BRepBuilderAPI import BRepBuilderAPI_MakeEdge, BRepBuilderAPI_MakeWire, BRepBuilderAPI_MakeFace,BRepBuilderAPI_MakeVertex

# Load the IGES file
iges_file_path = "C:\\Users\\oslon\\Desktop\\MSc\\Thesis\\Code Scripts\\trimmedRect.igs" # Replace with your IGES file path
iges_reader = IGESControl_Reader()
status = iges_reader.ReadFile(iges_file_path)

if not status:
    raise Exception(f"Failed to read the IGES file: {iges_file_path}")

# Transfer the entities from the IGES file to a TopoDS_Compound
compound = TopoDS_Compound()
iges_reader.TransferRoots()
compound = read_iges_file(iges_file_path)

# Extract a face and a curve (assumed to be the first face and edge)
explorer = TopologyExplorer(compound)
face, curve = None, None


for face in explorer.faces():
    rectangular_face = face

for edge in explorer.edges():
    trimming_curve = edge
    #curve = explorer.wire_from_edge(edge)


# Create an OCC Wire from the curve
#oc_curve = make_wire(curve)
oc_curve=BRepBuilderAPI_MakeEdge(curve)

# Convert OCC Wire to a list of control points in (x, y, z) coordinates
control_points_xyz = [(pnt.X(), pnt.Y(), pnt.Z()) for pnt in oc_curve.discretize(100)]

# Create a rectangular surface for visualization (you can replace this with your own surface)
box = BRepPrimAPI_MakeBox(1, 1, 1).Shape()
face = BRepBuilderAPI_MakeFace(topods_Wire(oc_curve.Shape()))
surface = BRepOffsetAPI_ThruSections(True).Shape()
compound = BRepBuilderAPI_MakeWire([surface.Edges()]).Shape()  # Corrected line

# Convert OCC Wire to a list of control points in (u, v) coordinates
control_points_uv = []
for x, y, z in control_points_xyz:
    pnt = gp_Pnt(x, y, z)
    projection = GeomAPI_ProjectPointOnSurf(pnt, surface)
    u, v = projection.LowerDistanceParameters()
    control_points_uv.append((u, v))

# Visualize the surface and control points
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

# Plot the surface
oc_surface = make_wire(surface)
oc_surface.plot(ax, color='b')

# Plot the control points
control_points_uv = list(zip(*control_points_uv))
ax.scatter(*control_points_uv, c='r', marker='o', label='Control Points (u, v)')

ax.set_xlabel('U')
ax.set_ylabel('V')
ax.set_zlabel('Z')
ax.legend()

plt.show()
