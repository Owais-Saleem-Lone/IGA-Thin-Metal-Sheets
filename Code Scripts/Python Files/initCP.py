from OCC.Extend.DataExchange import read_iges_file
from OCC.Core.TopoDS import topods
from OCC.Core.TopoDS import topods_Edge
from OCC.Extend.DataExchange import read_iges_file
from OCC.Core.IGESControl import IGESControl_Reader
from OCC.Core.TopExp import TopExp_Explorer
from OCC.Core.TopAbs import TopAbs_EDGE
from OCC.Core.BRep import BRep_Tool
# Load the IGES file
iges_file_path = "C:\\Users\\oslon\\Desktop\\MSc\\basic.igs"
shape = read_iges_file(iges_file_path)
# Check if the shape is a compound
compound = topods.Compound(shape)

# Read IGES file and transfer geometry
iges_reader = IGESControl_Reader()
status = iges_reader.ReadFile(iges_file_path)
if status == 0:
    raise Exception(f"Failed to read the IGES file: {iges_file_path}")

iges_reader.TransferRoots()
# Get the shape from the IGES file
shape = iges_reader.OneShape()
                # CODE B
edge_explorer = TopExp_Explorer(shape, TopAbs_EDGE)
while edge_explorer.More():
    edge = topods_Edge(edge_explorer.Current())
    curve,_= BRep_Tool.Curve(edge)

    # Use the Geom_Curve methods directly
    num_poles = curve.NbPoles()
    degree = curve.Degree()
    control_points = [curve.Pole(i + 1) for i in range(num_poles)]

    print(f"Degree: {degree}")
    print(f"Number of Control Points: {num_poles}")
    print("Control Points:")
    for i, point in enumerate(control_points):
        print(f"  P{i + 1}: {point.X()}, {point.Y()}, {point.Z()}")

    edge_explorer.Next()