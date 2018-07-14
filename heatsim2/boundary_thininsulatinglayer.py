# ... Use of group() here is to collect differences that
# shouldn't be simplified down, so that we can 
import collections
from heatsim2.expression import group

# This assumes conductioncoefficient is the conduction coefficient
# normal to the boundary, however the boundary is arranged. 

def qz(kmatm55,kmatp55,
       dz,dy,dx,
       Tm55,Tp55,
       Tm45,Tp45,
       Tm65,Tp65,
       Tm54,Tp54,
       Tm56,Tp56,
       Tm46,Tp46,
       Tm64,Tp64,
       conductioncoefficient):
    dT=group(Tp55-Tm55)
    
    # conductioncoefficient in Watts/(m^2*deg K)

    qz=-conductioncoefficient*dT # gives q in Watts/m^2
    
    return qz

def qy(kmat5m5,kmat5p5,
       dz,dy,dx,
       T5m5,T5p5,
       T4m5,T4p5,
       T6m5,T6p5,
       T5m4,T5p4,
       T5m6,T5p6,
       T4m6,T4p6,
       T6m4,T6p4,
       conductioncoefficient):

    dT=group(T5p5-T5m5)

    qy=-conductioncoefficient*dT
    
    return qy

def qx(kmat55m,kmat55p,
       dz,dy,dx,
       T55m,T55p,
       T45m,T45p,
       T65m,T65p,
       T54m,T54p,
       T56m,T56p,
       T46m,T46p,
       T64m,T64p,
       conductioncoefficient):

    dT=group(T55p-T55m)

    qx=-conductioncoefficient*dT
    
    return qx


