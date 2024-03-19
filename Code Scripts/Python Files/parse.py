import sys
from scipy.io import savemat
import numpy as np
from OCC.Core.BRepBuilderAPI import BRepBuilderAPI_NurbsConvert,BRepBuilderAPI_MakeEdge, BRepBuilderAPI_MakeWire
from OCC.Core.BRepAdaptor import BRepAdaptor_Curve,BRepAdaptor_Surface,BRepAdaptor_CompCurve
from OCC.Core.IGESControl import IGESControl_Reader
from OCC.Extend.DataExchange import read_iges_file
from OCC.Extend.TopologyUtils import TopologyExplorer
from OCC.Core.GeomAbs import GeomAbs_BSplineSurface
from OCC.Core.GeomAbs import GeomAbs_BSplineCurve
from OCC.Core import TopoDS
from OCC.Core import TopAbs
from OCC.Core.BRep import BRep_Tool
from OCC.Core.Geom import Geom_BSplineCurve,Geom_BSplineSurface
from OCC.Core.Geom2dConvert import geom2dconvert  
from OCC.Core.TColgp import TColgp_Array1OfPnt2d
from OCC.Core.TColStd import TColStd_Array1OfInteger, TColStd_Array1OfReal
from OCC.Core.TopExp import TopExp_Explorer  
from OCC.Core.TopAbs import TopAbs_EDGE,TopAbs_WIRE,TopAbs_FACE
from OCC.Core.TopoDS import topods
from OCC.Core.ShapeAnalysis import shapeanalysis
from OCC.Display.SimpleGui import init_display

output_file = 'output.mat'

def getWireData(_OCC_outerWire, _OCC_wire):
  """
  Receives a TopoDS_SHAPE class with a Wire type and converts its data to python type
  \param[in] _OCC_wire	class: TopoDS_SHAPE	type: Wire
  \return inner
  \return noTrCurves
  """
  # Determine inner or outer boundary loop
  if _OCC_wire.IsEqual(_OCC_outerWire):
    inner = False
  else:
    inner = True
  
  # Count number of trimming curves in the loop
  noTrCurves = 0  
  edgeExplorer = TopExp_Explorer(_OCC_wire, TopAbs_EDGE)
  while edgeExplorer.More():
    noTrCurves += 1
    edgeExplorer.Next()
  
  return inner, noTrCurves


def getOuterWire(_OCC_face):
  """
  \param[in] _OCC_face	class: TopoDS_SHAPE	type: Face
  \return OCC_outerWire
  """
  
  OCC_outerWire = shapeanalysis.OuterWire(topods.Face(_OCC_face))

  return OCC_outerWire

def getEdgeDataWithOrientation(_OCC_edge, _OCC_wire, _OCC_face):
  # Convert a topological TopoDS_EDGE into Geom2d_BSplineCurve
  OCC_handle_curve2d, activeRangeBegin, activeRangeEnd = BRep_Tool.CurveOnSurface(topods.Edge(_OCC_edge), topods.Face(_OCC_face))
  OCC_bsplineCurve2d = geom2dconvert.CurveToBSplineCurve(OCC_handle_curve2d)#.GetObject()
  
  if OCC_bsplineCurve2d.IsPeriodic(): 
    print("Found periodic curve!")
    OCC_bsplineCurve2d.SetNotPeriodic()
  # Curve direction 
  direction = (_OCC_wire.Orientation() + _OCC_edge.Orientation() + 1)%2
  # pDegree
  pDegree = OCC_bsplineCurve2d.Degree()
  # No of Knots
  OCC_uNoKnots = OCC_bsplineCurve2d.NbKnots()
  # Knot vector
  OCC_uKnotMulti = TColStd_Array1OfInteger(1,OCC_uNoKnots)
  OCC_bsplineCurve2d.Multiplicities(OCC_uKnotMulti)
  uNoKnots = 0
  for ctr in range(1,OCC_uNoKnots+1): uNoKnots += OCC_uKnotMulti.Value(ctr)
  OCC_uKnotSequence = TColStd_Array1OfReal(1,uNoKnots)
  OCC_bsplineCurve2d.KnotSequence(OCC_uKnotSequence)
  uKnotVector = []
  for iKnot in range(1,uNoKnots+1): uKnotVector.append(OCC_uKnotSequence.Value(iKnot))
  # No of CPs
  uNoCPs = OCC_bsplineCurve2d.NbPoles()
  # CPs
  OCC_CPnet = TColgp_Array1OfPnt2d(1,uNoCPs)
  OCC_bsplineCurve2d.Poles(OCC_CPnet)
  # CP weights
  OCC_CPweightNet = TColStd_Array1OfReal(1,uNoCPs)
  OCC_bsplineCurve2d.Weights(OCC_CPweightNet)
  CPNet = []
  for iCP in range(1,uNoCPs+1):
    OCC_CP = OCC_CPnet.Value(iCP)
    OCC_CPweight = OCC_CPweightNet.Value(iCP)
    CPNet.extend([OCC_CP.X(), OCC_CP.Y(), 0.0, OCC_CPweight])
# Create a linear space with 25 points
  start_value = uKnotVector[0]
  end_value = uKnotVector[-1] 
  num_points = 25
  linear_space = np.linspace(start_value, end_value, num_points)
  num_points = 25
  points_on_curve = [OCC_bsplineCurve2d.Value(i) for i in linear_space ]
  x_coordinates=[]
  y_coordinates=[]
  for point in points_on_curve:
    x_coordinates.append(point.X())
    y_coordinates.append(point.Y())
  return x_coordinates,y_coordinates

def getNURBS(iges_file_path):
    base_shape = read_iges_file(iges_file_path)
    iges_reader = IGESControl_Reader()
    iges_reader.ReadFile(iges_file_path)
    iges_reader.TransferRoots()
    num_shapes=iges_reader.NbShapes()
    subshape=iges_reader.Shape(1)
    
    nurbs_converter = BRepBuilderAPI_NurbsConvert(base_shape, True)
    converted_shape = nurbs_converter.Shape()
    expl = TopologyExplorer(converted_shape)
    
    if subshape.ShapeType()==TopAbs.TopAbs_EDGE:
        geo_dim=1
    elif subshape.ShapeType()==TopAbs.TopAbs_FACE:
        geo_dim=2
    else:
        print('Please input a valid curve or surface')

    if geo_dim==1:
        dict_list=[]
        crv_idx = 1
        for curve in expl.edges():
            print("=== Curve %i ===" % crv_idx)
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
            crv_idx += 1
            my_dict = {'polDegree': curve_Degree,'knotVector':knotArray,
                       'weights':weightArray,'controlPoints':controlPointMatrix.T}
            dict_list.append(my_dict)  
            
        data_to_save = {'dict_list': dict_list, 'dimension': geo_dim}
        savemat(output_file, data_to_save)
        print(f'All {num_shapes} dictionaries saved as {output_file}')
        return  output_file
        
    elif geo_dim==2:
        
        fc_idx = 1
        dict_list=[]
        for face in expl.faces():
            print("=== Face %i ===" % fc_idx)
            surf = BRepAdaptor_Surface(face, True)
            surf_type = surf.GetType()
            if not surf_type == GeomAbs_BSplineSurface:
                raise AssertionError("the face was not converted to a GeomAbs_BSplineSurface")
            
            bsrf = surf.BSpline()
            UCurveDegree = bsrf.UDegree()
            VCurveDegree = bsrf.VDegree()
            
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
            
            
            OCC_outerWire = getOuterWire(face)
  
  
            faceWireExplorer = TopExp_Explorer(face, TopAbs_WIRE) # This will explore outer and inner wires
            
            while faceWireExplorer.More():        # Loop through all wires (outer and inner)
                
                currentFaceWire = faceWireExplorer.Current()
                faceWireExplorer.Next()

                # Collect wire data
                inner, noTrCurves = getWireData(OCC_outerWire, currentFaceWire)   # check if the wire is exterior or interior
                #note: get wire data takes a wire and checks if it is exterior or interior. noTrCurves returns number...  
                # of edges in that wire (external or internal) 
                # Gather patch trimming curves info
                faceWireEdgeExplorer = TopExp_Explorer(currentFaceWire, TopAbs_EDGE)
                
                X_trim_vec=[]
                Y_trim_vec=[]
                if inner:           # For only interior trimming wire which is one in my case of circular trimming curve
                
                    while faceWireEdgeExplorer.More():        # Checking the edges in the wire, 4 in my case
                        currentFaceWireEdge = faceWireEdgeExplorer.Current()
                        faceWireEdgeExplorer.Next()
                        x,y = getEdgeDataWithOrientation(currentFaceWireEdge, currentFaceWire, face)
                        X_trim_vec.extend(x)
                        Y_trim_vec.extend(y)
                    
                    my_dict = {'isTrim':1, 'uDeg': UCurveDegree,'vDeg': VCurveDegree,'uknotVector':UknotArr,
                       'vknotVector':VknotArr,'weights':weightArray,'controlPoints':conPointMtx.T,'TrimPolygon':[X_trim_vec,Y_trim_vec]}
                    
                else:
                    my_dict = {'isTrim':0, 'uDeg': UCurveDegree,'vDeg': VCurveDegree,'uknotVector':UknotArr,
                     'vknotVector':VknotArr,'weights':weightArray,'controlPoints':conPointMtx.T}
                
            fc_idx += 1
            dict_list.append(my_dict)  
            
        data_to_save = {'dict_list': dict_list, 'dimension': geo_dim}
        savemat(output_file, data_to_save)
        print(f'All {num_shapes} dictionaries saved as {output_file}')

    else:
        raise AssertionError("The geometry is neither a curve not a surface")

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python your_script.py <file_path>")
        sys.exit(1)
    # Get the file path from the command-line arguments
    file_path = sys.argv[1]

    output_file = getNURBS(file_path)

    