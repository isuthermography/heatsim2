# ... Use of group() here is to collect differences that
# shouldn't be simplified down, so that we can 
import collections
from heatsim2.expression import group


def qz(kmatm55,kmatp55,
       dz,dy,dx,
       Tm55,Tp55,
       Tm45,Tp45,
       Tm65,Tp65,
       Tm54,Tp54,
       Tm56,Tp56,
       Tm46,Tp46,
       Tm64,Tp64):
    # q = -kgrad T
    # q = -k (dT/dx xhat + dT/dy yhat + dT/dz zhat)
    dTdz=group((Tp55-Tm55)*(1.0/dz))
    #dTdy=group((Tp65-Tp45)*(0.25/dy)) + group((Tm65-Tm45)*(0.25/dy))
    #dTdx=group((Tp56-Tp54)*(0.25/dx)) + group((Tm56-Tm54)*(0.25/dx))
    
    kmat=(kmatm55+kmatp55)*0.5

    # scalar k
    qz=-kmat*dTdz
    
    return qz

def qy(kmat5m5,kmat5p5,
       dz,dy,dx,
       T5m5,T5p5,
       T4m5,T4p5,
       T6m5,T6p5,
       T5m4,T5p4,
       T5m6,T5p6,
       T4m6,T4p6,
       T6m4,T6p4):
    # q = -kgrad T
    # q = -k (dT/dx xhat + dT/dy yhat + dT/dz zhat)
    #dTdz=group((T6p5-T4p5)*(0.25/dz)) + group((T6m5-T4m5)*(0.25/dz))
    dTdy=group((T5p5-T5m5)*(1.0/dy))
    #dTdx=group((T5p6-T5p4)*(0.25/dx)) + group((T5m6-T5m4)*(0.25/dx))

    kmat=(kmat5m5+kmat5p5)*0.5

    # scalar k
    qy=-kmat*dTdy
    
    return qy

def qx(kmat55m,kmat55p,
       dz,dy,dx,
       T55m,T55p,
       T45m,T45p,
       T65m,T65p,
       T54m,T54p,
       T56m,T56p,
       T46m,T46p,
       T64m,T64p):
    # q = -kgrad T
    # q = -k (dT/dx xhat + dT/dy yhat + dT/dz zhat)
    #dTdz=group((T65p-T45p)*(0.25/dz)) + group((T65m-T45m)*(0.25/dz))
    #dTdy=group((T56p-T54p)*(0.25/dy)) + group((T55m-T54m)*(0.25/dy))
    dTdx=group((T55p-T55m)*(1.0/dx))

    kmat=(kmat55m+kmat55p)*0.5
    
    # scalar k
    qx=-kmat*dTdx
    
    return qx


