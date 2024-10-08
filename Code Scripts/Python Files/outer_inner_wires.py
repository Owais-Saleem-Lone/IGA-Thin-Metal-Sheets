import sys
from scipy.io import savemat
from scipy.spatial import ConvexHull
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
from OCC.Core.TopAbs import TopAbs_EDGE,TopAbs_WIRE,TopAbs_FACE,TopAbs_FORWARD,TopAbs_REVERSED
from OCC.Core.TopoDS import topods
from OCC.Core.ShapeAnalysis import shapeanalysis
from OCC.Display.SimpleGui import init_display

output_file = 'output.mat'

def checkInner(_OCC_outerWire, _OCC_wire):

  # Determine inner or outer boundary loop
  if _OCC_wire.IsEqual(_OCC_outerWire):
    inner = False
    print(f"Outer Wire Detected")
  else:
    inner = True
    print(f"Inner Wire Detected")

  
  # Count number of trimming curves in one trimming loop
  #noTrCurves = 0  
  #edgeExplorer = TopExp_Explorer(_OCC_wire, TopAbs_EDGE)
  #while edgeExplorer.More():
   # noTrCurves += 1
    #edgeExplorer.Next()
  return inner

def getOuterWire(_OCC_face):
    
    OCC_outerWire = shapeanalysis.OuterWire(topods.Face(_OCC_face))
    faceWireEdgeExplorer = TopExp_Explorer(OCC_outerWire, TopAbs_EDGE)
    
    if OCC_outerWire.Orientation()==TopAbs_FORWARD:
        print(f" The orientation of outer wire edge is anti-clockwise")
    else:
        print(f" The orientation of outer wire edge is clockwise")

    counter=1
    out_edge_dict_list=[]
    x_outer_wire=[]
    y_outer_wire=[]
    while faceWireEdgeExplorer.More():        # Checking the edges in outer the wire, 4 in my case
        print(f" Outer Edge Number {counter}")
        currentFaceWireEdge = faceWireEdgeExplorer.Current()
        faceWireEdgeExplorer.Next()
        out_edge_dict,x_edge,y_edge=getEdgeDataWithOrientation(currentFaceWireEdge,_OCC_face)
        x_outer_wire.extend(x_edge)
        y_outer_wire.extend(y_edge)
        counter +=1
        out_edge_dict_list.append(out_edge_dict)
    
    x_outer_wire=np.array(x_outer_wire)
    x_outer_wire = x_outer_wire.reshape(-1, 1)
    y_outer_wire=np.array(y_outer_wire)
    y_outer_wire = y_outer_wire.reshape(-1, 1)
    polygon=rearrange(np.hstack((x_outer_wire, y_outer_wire)),len(out_edge_dict_list))

    #print("The size of y_outer_wire is ", y_outer_wire.size)
    #print("The number of edges in outer_wire are ",len(out_edge_dict_list))

    out_wire_dict={"outer_edge_dict_list":out_edge_dict_list,"outer_Trim_Polygon":polygon}
    #print(f" The outer wire data before processing looks like\n", out_wire_dict['outer_Trim_Polygon'])
    return OCC_outerWire,out_wire_dict

def getEdgeDataWithOrientation(_OCC_edge, _OCC_face):
    # Convert a topological TopoDS_EDGE into Geom2d_BSplineCurve
    # _OCC_edge is one of the edges of the _OCC_wire
    OCC_handle_curve2d, activeRangeBegin, activeRangeEnd = BRep_Tool.CurveOnSurface(topods.Edge(_OCC_edge), topods.Face(_OCC_face))
    OCC_bsplineCurve2d = geom2dconvert.CurveToBSplineCurve(OCC_handle_curve2d)#.GetObject()
    # pDegree
    pDegree = OCC_bsplineCurve2d.Degree()
    #print(f" The edge degree is {pDegree}")
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
    #print(f" The edge knot vector  is {uKnotVector}")

    # No of CPs
    poles=OCC_bsplineCurve2d.Poles()
    # CP weights
    weights = OCC_bsplineCurve2d.Weights()
    weightArray=[]
    if weights is not None:
        for i in range(OCC_bsplineCurve2d.NbPoles()):
            weightArray.append(weights.Value(i+1))
    else:
        weightArray=np.ones(OCC_bsplineCurve2d.NbPoles())
    #print(f" The edge weights are {weightArray}")

    controlPointMatrix= np.zeros((OCC_bsplineCurve2d.NbPoles(),2))
    if poles is not None:
        for i in range(OCC_bsplineCurve2d.NbPoles()):
                p = poles.Value(i + 1)
                controlPointMatrix[i]= np.array([p.X(),p.Y(),])
    #print(f" The edge control Points are {controlPointMatrix.T}")

    edge_dict = {'polDegree': pDegree,'knotVector':uKnotVector,
                        'weights':weightArray,'controlPoints':controlPointMatrix.T}
    
    # Create a linear space with 40 points
    start_value = uKnotVector[0]
    end_value = uKnotVector[-1]
    num_points = 40
    linear_space = np.linspace(start_value, end_value, num_points)
    points_on_curve = [OCC_bsplineCurve2d.Value(i) for i in linear_space]
    x_coordinates=[]
    y_coordinates=[]
    for point in points_on_curve:
        x_coordinates.append(point.Y())
        y_coordinates.append(point.X())
    return edge_dict,x_coordinates,y_coordinates

def getWires(_OCC_face):
    
    OCC_outerWire,outer_wire_dict = getOuterWire(_OCC_face)
  
    faceWireExplorer = TopExp_Explorer(_OCC_face, TopAbs_WIRE) # This will explore outer and inner wires
    inner_wire_dict_list=[]
    inner_wire_counter=1
    while faceWireExplorer.More():        # Loop through all wires (outer and inner)
        currentFaceWire = faceWireExplorer.Current()
        faceWireExplorer.Next()
        # Collect wire data
        inner_wire_dict=[]
        X_trim_vec=[]
        Y_trim_vec=[]
        isInner = checkInner(OCC_outerWire, currentFaceWire)   # check if the wire is exterior or interior
        counter=1
        if isInner:           # For only interior trimming wire which is one in my case of circular trimming curve
        #note: get wire data takes a wire and checks if it is exterior or interior. noTrCurves returns number...  
        # of edges in that wire (external or internal) 
        # Gather patch trimming curves info
            faceWireEdgeExplorer = TopExp_Explorer(currentFaceWire, TopAbs_EDGE)
            if currentFaceWire.Orientation()==TopAbs_REVERSED:
                print(f" The orientation of the inner wire is anti-clockwise")
            else:
                print(f" The orientation of inner wire is clockwise")
            inner_edge_dict=[]
            print(f"Inner Wire {inner_wire_counter}")
            while faceWireEdgeExplorer.More():        # Checking the edges in the wire, 4 in my case
                print(f"Inner Edge {counter}")
                currentFaceWireEdge = faceWireEdgeExplorer.Current()
                faceWireEdgeExplorer.Next()
                dict,x,y= getEdgeDataWithOrientation(currentFaceWireEdge, _OCC_face)
                X_trim_vec.append(x)
                Y_trim_vec.append(y)
                inner_edge_dict.append(dict)
                counter+=1
            X_trim_vec=np.array(X_trim_vec)
            X_trim_vec= X_trim_vec.reshape(-1, 1)
            Y_trim_vec=np.array(Y_trim_vec)
            Y_trim_vec=Y_trim_vec.reshape(-1,1)
            inner_polygon= rearrange(np.hstack((X_trim_vec,Y_trim_vec)),len(inner_edge_dict))
            inner_wire_dict={'inner_edge_dict_list':inner_edge_dict,'Trim_Polygon':inner_polygon}
            inner_wire_dict_list.append(inner_wire_dict)
            inner_wire_counter+=1
    return outer_wire_dict,isInner,inner_wire_dict_list       
 
 
def rearrange(data,num_edges):
    
    # length of an edge data
    len_dataset= data.shape[0]//num_edges
    reordered_data=np.zeros(data.shape)
    vertices = np.zeros((2*num_edges,2))

    # Ordering the given data set
    for i in range(num_edges-1):
        data_i= data[i*len_dataset:(i+1)*len_dataset,:]
        start_i=data_i[0,:]
        for j in range(i+1,num_edges):
            data_j= data[j*len_dataset:(j+1)*len_dataset,:]
            start_j=data_j[0,:]
            if start_i[0]==start_j[0] and start_i[1]==start_j[1]:
                data_j=data_j[::-1]
                data[j*len_dataset:(j+1)*len_dataset,:]=data_j
                break
    
    for i in range(num_edges-1):
        if i==0:
            data_i= data[i*len_dataset:(i+1)*len_dataset,:]
            reordered_data[i*len_dataset:(i+1)*len_dataset,:]=data_i
            vertices[2*i,:]=data_i[0,:]
            vertices[2*i+1,:]=data_i[-1,:]
            end_i=data_i[-1,:]
        else:
            data_i= reordered_data[i*len_dataset:(i+1)*len_dataset,:]
            end_i=data_i[-1,:]
        
        for j in range(num_edges):
            data_j= data[j*len_dataset:(j+1)*len_dataset,:]
            start_j=data_j[0,:]        
            if end_i[0]==start_j[0] and end_i[1]==start_j[1]:
                reordered_data[(i+1)*len_dataset:(i+2)*len_dataset,:]=data_j
                vertices[2*i,:]=data_j[0,:]
                vertices[2*i+1,:]=data_j[-1,:]
                break  
    #print(reordered_data)
    return reordered_data    
    
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
            print("**********************")
            print("=== Face %i ===" % fc_idx)
            print("**********************")

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
            #print(f"Surface UKnot Vector is: ",UknotArr)

            # vknots array
            vknots = bsrf.VKnots()
            
            VknotArr=[]
            for i in range(bsrf.NbVKnots()):
                VknotArr.append(vknots.Value(i+1))
            VknotArr = [VknotArr[0]] * (VCurveDegree+1) + VknotArr[1:-1] + [VknotArr[-1]] * (VCurveDegree+1)
            #print(f"Surface VKnot Vector is: ",VknotArr)

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
            #print(f" Surface Control Points are",conPointMtx.T)
            
            # Getting Inner and Outer Wires of each Face
            outer_wire_dict, isInner,inner_wire_dict,=getWires(face)
            if isInner:
                    
                my_dict = {'isTrim':1, 'uDeg': UCurveDegree,'vDeg': VCurveDegree,'uknotVector':UknotArr,
                    'vknotVector':VknotArr,'weights':weightArray,'controlPoints':conPointMtx.T,
                    'outer_wire':outer_wire_dict,'inner_wires':inner_wire_dict}
                print(f"the number of inner wires are: {len(inner_wire_dict)}\n")
                
            else:
                my_dict = {'isTrim':0, 'uDeg': UCurveDegree,'vDeg': VCurveDegree,'uknotVector':UknotArr,
                    'vknotVector':VknotArr,'weights':weightArray,'controlPoints':conPointMtx.T,'outer_wire':outer_wire_dict}
                print(f"the number of inner wires are: {len(inner_wire_dict)}\n")
                print(f"the number of outer edges in outer wire are: {len(outer_wire_dict['outer_edge_dict_list'])}\n")

            fc_idx += 1
            dict_list.append(my_dict)  
            
        data_to_save = {'dict_list': dict_list, 'dimension': geo_dim}
        savemat(output_file, data_to_save)
        print(f'All {num_shapes} dictionaries saved as {output_file}')

    else:
        raise AssertionError("The geometry is neither a curve not a surface")

if __name__ == "__main__":
    
    #iges_file_path="C:\\Users\\oslon\\Desktop\\MSc\\Thesis\\Code Scripts\\trimmedRect.igs"
    #iges_file_path="C:\\Users\\oslon\\Desktop\\MSc\\Thesis\\Code Scripts\\low_dim_rect_3patch.igs"
    #iges_file_path="C:\\Users\\oslon\\Desktop\\MSc\\Thesis\\Code Scripts\\low_dim_rect_3patch_trimmed.igs"
    #iges_file_path="C:\\Users\\oslon\\Desktop\\MSc\\Thesis\\Code Scripts\\basicRectSurf.igs"

    iges_file_path = sys.argv[1]
    output_file = getNURBS(iges_file_path)

    