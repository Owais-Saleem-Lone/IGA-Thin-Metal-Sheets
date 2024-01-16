import sys
import json
from scipy.io import savemat
import numpy as np
from OCC.Core.BRepPrimAPI import BRepPrimAPI_MakeTorus
from OCC.Core.BRepBuilderAPI import BRepBuilderAPI_NurbsConvert
from OCC.Core.BRepAdaptor import BRepAdaptor_Curve
from OCC.Core.BRepAdaptor import BRepAdaptor_Surface
from OCC.Core.IGESControl import IGESControl_Reader
from OCC.Extend.DataExchange import read_iges_file
from OCC.Extend.TopologyUtils import TopologyExplorer
from OCC.Core.GeomAbs import GeomAbs_BSplineSurface
from OCC.Core.GeomAbs import GeomAbs_BSplineCurve
from OCC.Core import TopoDS
from OCC.Core import TopAbs



def getNURBS(iges_file_path):
    base_shape = read_iges_file(iges_file_path)
    iges_reader = IGESControl_Reader()
    iges_reader.ReadFile(iges_file_path)
    iges_shape = TopoDS.TopoDS_Shape()
    iges_reader.TransferRoots()
    subshape=iges_reader.Shape(1)
    nurbs_converter = BRepBuilderAPI_NurbsConvert(base_shape, True)
    converted_shape = nurbs_converter.Shape()
    expl = TopologyExplorer(converted_shape)
    output_file = 'output.mat'
    if subshape.ShapeType()==TopAbs.TopAbs_EDGE:
        geo_dim=1
    elif subshape.ShapeType()==TopAbs.TopAbs_FACE:
        geo_dim=2
    else:
        print('Please input a valid curve or surface')

    if geo_dim==1:

        
        fc_idx = 1
        for curve in expl.edges():
            print("=== Curve %i ===" % fc_idx)
            surf = BRepAdaptor_Curve(curve)
            surf_type = surf.GetType()
            if not surf_type == GeomAbs_BSplineCurve:
                raise AssertionError("the face was not converted to a GeomAbs_BSplineSurface")
            bcurve = surf.BSpline()
            curve_Degree= bcurve.Degree()
            vknots = bcurve.Knots()
            #print("\nknots:", end="")
            knotArr=[]
            for i in range(bcurve.NbKnots()):
                knotArr.append(vknots.Value(i + 1))
                #print(vknots.Value(i + 1), end=" ")
            knotArray = [knotArr[0]] * (curve_Degree+1) + knotArr[1:-1] + [knotArr[-1]] * (curve_Degree+1)      

            weights = bcurve.Weights()
            weightArray=[]
            poles = bcurve.Poles()
            if weights is not None:
                for i in range(bcurve.NbPoles()):
                    weightArray.append(weights.Value(i+1))
            else:
                weightArray=np.ones(bcurve.NbPoles())
            controlPointMatrix= np.zeros((bcurve.NbPoles(),3))
            if poles is not None:
                for i in range(bcurve.NbPoles()):
                        p = poles.Value(i + 1)
                        controlPointMatrix[i]= np.array([p.X(),p.Y(),p.Z()])
            fc_idx += 1
            my_dict = {'dimension': geo_dim, 'polDegree': curve_Degree,'knotVector':knotArray,
                       'weights':weightArray,'controlPoints':controlPointMatrix.T}
            savemat(output_file, my_dict)
            return  output_file
        
    elif geo_dim==2:
        
        fc_idx = 1
        for face in expl.faces():
            print("=== Face %i ===" % fc_idx)
            surf = BRepAdaptor_Surface(face, True)
            surf_type = surf.GetType()
            if not surf_type == GeomAbs_BSplineSurface:
                raise AssertionError("the face was not converted to a GeomAbs_BSplineSurface")
            
            bsrf = surf.BSpline()
            udegree=[]
            UCurveDegree = bsrf.UDegree()
            
            VCurveDegree = bsrf.VDegree()

            
            degree= np.array([UCurveDegree, UCurveDegree])
            uknots = bsrf.UKnots()
            UknotArr=[]
            for i in range(bsrf.NbUKnots()):
                UknotArr.append(uknots.Value(i+1))
            UknotArr = [UknotArr[0]] * (UCurveDegree+1) + UknotArr[1:-1] + [UknotArr[-1]] * (UCurveDegree+1)
            #print(','.join(map(str, UknotArr)))
            
            # vknots array
            vknots = bsrf.VKnots()
            
            VknotArr=[]
            for i in range(bsrf.NbVKnots()):
                VknotArr.append(vknots.Value(i+1))
            VknotArr = [VknotArr[0]] * (VCurveDegree+1) + VknotArr[1:-1] + [VknotArr[-1]] * (VCurveDegree+1)
            #print(','.join(map(str, VknotArr)))

            # Weights 
            weights = bsrf.Weights()
            if weights is not None:
                #weightMatrix=np.zeros((bsrf.NbUPoles(),bsrf.NbVPoles()))
                weightArray=[]
                for i in range(bsrf.NbUPoles()):
                    for j in range(bsrf.NbVPoles()):
                        #weightMatrix[i,j] = weights.Value(i + 1, j + 1)
                        weightArray.append(weights.Value(i + 1, j + 1))      
            # weights can be None
            else:
                #weightMatrix=np.ones((bsrf.NbUPoles(),bsrf.NbVPoles()))
                weightArray=np.ones(bsrf.NbUPoles()*bsrf.NbVPoles())
            #print("\nweight: \n", weightMatrix)

            # control Points
            poles = bsrf.Poles()
            conPointMtx= np.zeros((bsrf.NbUPoles()*bsrf.NbVPoles(),3))
            num=0
            if poles is not None:
                for i in range(bsrf.NbUPoles()):
                    #print('\n')
                    for j in range(bsrf.NbVPoles()):
                        
                        p=poles.Value(i+1,j+1)
                        conPointMtx[num]= np.array([p.X(),p.Y(),p.Z()])
                        num=num+1

            #print(conPointMtx.T)
            fc_idx += 1
            my_dict = {'dimension': geo_dim, 'uDeg': UCurveDegree,'vDeg': VCurveDegree,'uknotVector':UknotArr,
                       'vknotVector':VknotArr,'weights':weightArray,'controlPoints':conPointMtx.T}
            
            savemat(output_file, my_dict)
            return  output_file

    else:
        raise AssertionError("The geometry is neither a curve not a surface")

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python your_script.py <file_path>")
        sys.exit(1)
    # Get the file path from the command-line arguments
    file_path = sys.argv[1]

    output_file = getNURBS(file_path)

    