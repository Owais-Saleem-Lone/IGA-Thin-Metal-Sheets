from OCC.Extend.ShapeFactory import make_edge, make_wire, make_face
from OCC.Extend.DataExchange import read_iges_file
from OCC.Core.BRepAdaptor import BRepAdaptor_Curve
from OCC.Core.GeomAbs import GeomAbs_BSplineCurve
from OCC.Extend.TopologyUtils import TopologyExplorer
from OCC.Core.BRepBuilderAPI import BRepBuilderAPI_NurbsConvert
from OCC.Core.TopExp import TopExp_Explorer
from OCC.Core.TopAbs import TopAbs_WIRE,TopAbs_EDGE
import numpy as np
from OCC.Extend.TopologyUtils import TopologyExplorer
from OCC.Core import Standard

iges_file_path = "C:\\Users\\oslon\\Desktop\\MSc\\afem-lab\\distortRect.igs"
base_shape = read_iges_file(iges_file_path)

#Conversion to a nurbs representation
nurbs_converter = BRepBuilderAPI_NurbsConvert(base_shape, True)
converted_shape = nurbs_converter.Shape()
# Loop over faces
fc_idx = 1
noFaceWires = 0
OCC_faceWireExplorer = TopExp_Explorer(converted_shape, TopAbs_WIRE)
while OCC_faceWireExplorer.More():
    
    wire = OCC_faceWireExplorer.Current()
    
    edge_explorer = TopExp_Explorer(wire, TopAbs_EDGE)

    # Loop over edges
    i=1
    while edge_explorer.More():
        print("Wire Number",i)
        edge = edge_explorer.Current()
        curve = BRepAdaptor_Curve(edge).Curve()
        edge_explorer = TopExp_Explorer(wire, TopAbs_EDGE)
            # Check if the curve is a B-spline curve
        if curve and curve.GetType() == GeomAbs_BSplineCurve:
            bspline_curve = curve.GetHandle()
            degree = bspline_curve.Degree()
            nb_poles = bspline_curve.NbPoles()

                # Get control points
            control_points = [bspline_curve.Pole(i + 1) for i in range(nb_poles)]
            control_points = np.array(control_points)

                # Print or process control points as needed
            print(f"Face {fc_idx}, Wire {noFaceWires}, Edge Control Points: {control_points}")

        edge_explorer.Next()
        i=i+1
        

    OCC_faceWireExplorer.Next()
    noFaceWires += 1


exit()







# fc_idx = 1
# for face in TopologyExplorer(converted_shape).faces():

#     wires = list(TopologyExplorer(face).wires())  # Convert the iterator to a list
    
#     # Check if the face has more than one wire (representing the hole)
#     if len(wires) >= 1:
#         print("=== Face %i ===" % fc_idx)
        
#         # Process only the inner wire (assuming it is the trimming curve)
#         wire = wires[1]
#         edge_idx = 1
#         for edge in TopologyExplorer(wire).edges():
#             print(f"=== Edge {edge_idx} ===")
            
#             bspline_curve = BRepAdaptor_Curve(edge)
#             curve_type = bspline_curve.GetType()

#             if curve_type == GeomAbs_BSplineCurve:
#                 #bspline_curve = curve_adaptor.BSpline()
#                 u_domain = bspline_curve.FirstParameter()
#                 v_domain = bspline_curve.LastParameter()
#                 print(f"U Domain: {u_domain}, V Domain: {v_domain}")

#                 # Print more details about the edge
#                 print("Curve Type:", bspline_curve.GetType())
#                 print("Curve Degree:", bspline_curve.Degree())
#                 print("Number of Poles:", bspline_curve.NbPoles())
#                 print("Knots:", [bspline_curve.Knot(i) for i in range(1, bspline_curve.NbKnots() + 1)])
#                 controlPointMatrix= np.zeros((bspline_curve.NbPoles(),3))
#                 poles = bspline_curve.Poles()
#                 if poles is not None:
#                     print("Poles (control points):", end="")
#                 for i in range(bspline_curve.NbPoles()):
#                     p = poles.Value(i + 1)
#                     controlPointMatrix[i]= np.array([p.X(),p.Y(),p.Z()])
#                     print(p.X(), p.Y(), p.Z(), end="\n")
                    
                

#             print()
#             edge_idx += 1

#     fc_idx+=1
