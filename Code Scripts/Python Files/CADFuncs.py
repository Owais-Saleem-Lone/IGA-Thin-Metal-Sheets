
from OCC.Core.BRepPrimAPI import BRepPrimAPI_MakeTorus
from OCC.Core.BRepBuilderAPI import BRepBuilderAPI_NurbsConvert
from OCC.Core.BRepAdaptor import BRepAdaptor_Curve
from OCC.Core.IGESControl import IGESControl_Reader
from OCC.Extend.DataExchange import read_iges_file
from OCC.Extend.TopologyUtils import TopologyExplorer
from OCC.Core.GeomAbs import GeomAbs_BSplineSurface
from OCC.Core.GeomAbs import GeomAbs_BSplineCurve
from OCC.Core.TopExp import TopExp_Explorer
from OCC.Core.TopAbs import TopAbs_WIRE
from OCC.Core.TopoDS import topods_Face
from OCC.Core.BRep import BRep_Tool
from OCC.Core.GeomConvert import geomconvert_SurfaceToBSplineSurface
from OCC.Core.TColgp import TColgp_Array2OfPnt
from OCC.Core.TColStd import TColStd_Array1OfInteger, TColStd_Array1OfReal, TColStd_Array2OfReal
from OCC.Core.Precision import precision

iges_file_path = "C:\\Users\\oslon\\Desktop\\MSc\\Thesis\\Code Scripts\\trimmedRect.igs"
_OCC_shapes = read_iges_file(iges_file_path)

shapeConverter = BRepBuilderAPI_NurbsConvert()
OCC_nurbsShapes = []
for OCC_shape in _OCC_shapes:
  print("Converting Shape: ", OCC_shape, " to NURBS")
  shapeConverter.Perform(OCC_shape)
  currentShape = shapeConverter.Shape()
  OCC_nurbsShapes.append(currentShape)
  


OCC_precision = precision()
  
  # Trimming curves
_startID=0
noFaceWires = 0
OCC_faceWireExplorer = TopExp_Explorer(OCC_shape, TopAbs_WIRE)
while OCC_faceWireExplorer.More():
 
  OCC_faceWireExplorer.Next()
  noFaceWires += 1

  # Convert a topological TopoDS_FACE into Geom_BSplineSurface
OCC_handle_surface = BRep_Tool.Surface(topods_Face(OCC_shape))
OCC_bsplineSurface = geomconvert_SurfaceToBSplineSurface(OCC_handle_surface).GetObject()
    
pDegree = OCC_bsplineSurface.UDegree()
qDegree = OCC_bsplineSurface.VDegree()
OCC_uNoKnots = OCC_bsplineSurface.NbUKnots()
OCC_vNoKnots = OCC_bsplineSurface.NbVKnots()
uKnotVector = []
for iKnot in range(1, OCC_bsplineSurface.NbUKnots()+1):
  for iMult in range(0,OCC_bsplineSurface.UMultiplicity(iKnot)):
    uKnotVector.append(OCC_bsplineSurface.UKnot(iKnot))
  
vKnotVector = []
for iKnot in range(1, OCC_bsplineSurface.NbVKnots()+1):
  for iMult in range(0,OCC_bsplineSurface.VMultiplicity(iKnot)):
    vKnotVector.append(OCC_bsplineSurface.VKnot(iKnot))

if OCC_bsplineSurface.IsUPeriodic():
  print ("Found periodic U")
  uKnotVector.insert(0, uKnotVector[0])
  uKnotVector.append(uKnotVector[len(uKnotVector)-1])
  OCC_bsplineSurface.SetUNotPeriodic()
if OCC_bsplineSurface.IsVPeriodic():
  print ("Found periodic V")
  vKnotVector.insert(0, vKnotVector[0])
  vKnotVector.append(vKnotVector[len(vKnotVector)-1])
  OCC_bsplineSurface.SetVNotPeriodic()
  
  # No of CPs
uNoCPs = OCC_bsplineSurface.NbUPoles()
vNoCPs = OCC_bsplineSurface.NbVPoles()
  # CPs
OCC_CPnet = TColgp_Array2OfPnt(1,uNoCPs,1,vNoCPs)
OCC_bsplineSurface.Poles(OCC_CPnet)
  # CP weights
OCC_CPweightNet = TColStd_Array2OfReal(1,uNoCPs,1,vNoCPs)
OCC_bsplineSurface.Weights(OCC_CPweightNet)
  
CPNet = []
uniqueIDNet = []
for iVCP in range(1,vNoCPs+1):
  for iUCP in range(1,uNoCPs+1):
    OCC_CP = OCC_CPnet.Value(iUCP,iVCP)
    OCC_CPweight = OCC_CPweightNet.Value(iUCP,iVCP)
    CPNet.extend([OCC_CP.X(), OCC_CP.Y(), OCC_CP.Z(), OCC_CPweight])
    uniqueIDNet.append(_startID)
    startID += 1
     
