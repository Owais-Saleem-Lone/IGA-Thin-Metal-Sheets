from OCC.Extend.ShapeFactory import make_edge, make_wire, make_face
from OCC.Extend.DataExchange import read_iges_file
from OCC.Core.BRepAdaptor import BRepAdaptor_Curve
from OCC.Core.GeomAbs import GeomAbs_BSplineCurve
from OCC.Extend.TopologyUtils import TopologyExplorer
from OCC.Core.BRepBuilderAPI import BRepBuilderAPI_NurbsConvert
iges_file_path = "C:\\Users\\oslon\\Desktop\\MSc\\afem-lab\\trimmedSurface.igs"
base_shape = read_iges_file(iges_file_path)

# Conversion to a nurbs representation
nurbs_converter = BRepBuilderAPI_NurbsConvert(base_shape, True)
converted_shape = nurbs_converter.Shape()

# Loop over faces
fc_idx = 1

for face in TopologyExplorer(converted_shape).faces():

    wires = list(TopologyExplorer(face).wires())  # Convert the iterator to a list
    
    # Check if the face has more than one wire (representing the hole)
    if len(wires) >= 1:
        print("=== Face %i ===" % fc_idx)
        
        # Process only the inner wire (assuming it is the trimming curve)
        wire = wires[1]
        edge_idx = 1
        for edge in TopologyExplorer(wire).edges():
            print(f"=== Edge {edge_idx} ===")
            
            curve_adaptor = BRepAdaptor_Curve(edge)
            curve_type = curve_adaptor.GetType()

            if curve_type == GeomAbs_BSplineCurve:
                bspline_curve = curve_adaptor.BSpline()
                u_domain = bspline_curve.FirstParameter()
                v_domain = bspline_curve.LastParameter()
                print(f"U Domain: {u_domain}, V Domain: {v_domain}")

                # Print more details about the edge
                print("Curve Type:", curve_adaptor.GetType())
                print("Curve Degree:", bspline_curve.Degree())
                print("Number of Poles:", bspline_curve.NbPoles())
                print("Knots:", [bspline_curve.Knot(i) for i in range(1, bspline_curve.NbKnots() + 1)])

            else:
                print("Edge is not a BSpline curve.")

            print()
            edge_idx += 1

    fc_idx+=1
