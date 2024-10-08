import sys
import json
import numpy as np
from scipy.io import savemat
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

def getNURBS(iges_file_path):
    base_shape = read_iges_file(iges_file_path)

    # TRIAL TRIAL TRIAL
    iges_reader = IGESControl_Reader()
    iges_reader.ReadFile(iges_file_path)
    iges_shape = TopoDS.TopoDS_Shape()
    iges_reader.TransferRoots()
    subshape=iges_reader.Shape(1)
    if subshape.ShapeType()==TopAbs.TopAbs_EDGE:
        geo_dim=1
    elif subshape.ShapeType()==TopAbs.TopAbs_FACE:
        geo_dim=2
    else:
        print('Please input a valid curve or surface')
    # TRIAL TRIAL TRIAL



    # conversion to a nurbs representation
    nurbs_converter = BRepBuilderAPI_NurbsConvert(base_shape, True)
    # nurbs_converter.Perform()
    converted_shape = nurbs_converter.Shape()
    expl = TopologyExplorer(converted_shape)
    # loop over all the curves
    crv_idx = 1
    for crv in expl.edges():
        print("=== Face %i ===" % crv_idx)
        curve = BRepAdaptor_Curve(crv)
        curve_type = curve.GetType()
        # check each of the is a BSpline surface
        # it should be, since we used the nurbs converter before
        if not curve_type == GeomAbs_BSplineCurve:
            raise AssertionError("the face was not converted to a GeomAbs_BSplineSurface")
        bcurve = curve.BSpline()

    # Polynomail Order
        curve_Degree= bcurve.Degree()
        #print("UDegree:", bcurve.Degree())
    
    # Knots array
        vknots = bcurve.Knots()
        #print("\nknots:", end="")
        knotArr=[]
        for i in range(bcurve.NbKnots()):
            knotArr.append(vknots.Value(i + 1))
            #print(vknots.Value(i + 1), end=" ")
        knotArray = [knotArr[0]] * (curve_Degree+1) + knotArr[1:-1] + [knotArr[-1]] * (curve_Degree+1)      

    # Weights and control Points
        weights = bcurve.Weights()
        weightArray=[]
        poles = bcurve.Poles()
        if weights is not None:
            #print("Weights:", end="")
            for i in range(bcurve.NbPoles()):
                weightArray.append(weights.Value(i+1))
                #print(weights.Value(i + 1), end=" ")
        else:
            weightArray=np.ones(bcurve.NbPoles())
        # weights can be None
        controlPointMatrix= np.zeros((bcurve.NbPoles(),3))
        if poles is not None:
            #print("Poles (control points):", end="")
            for i in range(bcurve.NbPoles()):
                    p = poles.Value(i + 1)
                    controlPointMatrix[i]= np.array([p.X(),p.Y(),p.Z()])
                    #print(p.X(), p.Y(), p.Z(), end=" ")
        print()
        crv_idx += 1
        return  geo_dim,curve_Degree, knotArray, weightArray, controlPointMatrix.T

if __name__ == "__main__":


    # if len(sys.argv) != 2:
    #     print("Usage: python your_script.py <file_path>")
    #     sys.exit(1)
    # file_path = sys.argv[1]

    # Call the processing function
    
    file_path = "C:\\Users\\oslon\\Desktop\\MSc\\Thesis\\Code Scripts\\basic.igs"
    dim,degree, knot,wt, cp = getNURBS(file_path)

    print(f"Dimension is {dim}\n")
    print(f"Degree is {degree}\n")
    print("Knot vector\n")
    print(','.join(map(str, knot)))
    print("weight vector\n")
    print(','.join(map(str, wt)))
    print("control points \n")
    print('\n'.join(','.join(map(str, row)) for row in cp))