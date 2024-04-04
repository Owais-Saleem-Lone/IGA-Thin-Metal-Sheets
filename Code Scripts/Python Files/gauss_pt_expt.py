import numpy as np

def gauss_points(order):
    # Get Gauss-Legendre quadrature points and weights
    points, weights = np.polynomial.legendre.leggauss(order)
    return points,weights

def set_gps_curve(knot_vec,degree):
  knot_spans=np.unique(knot_vec)
  curve_deg=degree
  gps,wts=gauss_points(degree+1)
  print ("Gauss Points  ",gps)
  print ("Gauss weights  ",wts)
  a = len(gps)
  b = len(knot_spans)
  crv_gauss_mtx = np.zeros((a+1, b))
  crv_gauss_val=[]
  crv_gauss_wts=[]
  for i in range(len(knot_spans) - 1):
    xi_crv_1=knot_spans[i]
    xi_crv_2=knot_spans[i+1]
    for j in range(len(gps)):
        if i == 0:
            crv_gauss_mtx[j,b-1]=crv_gauss_mtx[j,b-1]+wts[j]
        crv_gauss_mtx[j,i]= crv_gauss_mtx[j,i]+gps[j]*0.5*(xi_crv_2-xi_crv_1) + 0.5*(xi_crv_1+xi_crv_2)
        crv_gauss_val.append(crv_gauss_mtx[j,i])
    crv_gauss_mtx[a,i] =  crv_gauss_mtx[a,i]+0.5*(xi_crv_2-xi_crv_1)
  print(crv_gauss_mtx)
  print(crv_gauss_val)
  
knot_vec=[0, 0, 0, 2,4,6,8, 10, 10, 10]
set_gps_curve(knot_vec,2)
