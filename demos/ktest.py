import heatsim2
import copy
import cProfile
import numpy as np
numpy=np
import pylab as pl
import scipy as sp
import scipy.interpolate



# define materials:
composite_rho=1.75e3 # kg/m^3
composite_c=730.0 # J/kg/deg K 

## Use these values for (roughly) quasi-isotropic layup specimens
#composite_k=3.7 # W/m/deg K
#max_x = 51.0e-3 # length of specimen

## Use these values for (roughly) 0/90 layup specimens
#composite_k=4.0 # W/m/deg K
#max_x = 51.0e-3 # length of specimen

## Use these values for (roughly) 0 uni layup specimens
#composite_k=5.1 # W/m/deg K
#max_x = 51.0e-3 # length of specimen

# Use these values for (roughly) 90 uni layup specimens
composite_k=0.71 # W/m/deg K
max_x = 25.0e-3 # length of specimen



# WARNING: X and Y GRID BOUNDARIES MUST LINE UP WITH
# composite_min_z and composite_max_z FOR AN ACCURATE CALCULATION

# Create x,y,z voxel center coords
nx=15
ny=40
nz=20
(dz,dy,dx,
 z,y,x,
 zgrid,ygrid,xgrid,
 z_bnd,y_bnd,x_bnd,
 z_bnd_z,z_bnd_y,z_bnd_x,
 y_bnd_z,y_bnd_y,y_bnd_x,
 x_bnd_z,x_bnd_y,x_bnd_x,
 r3d,r2d) = heatsim2.build_grid(-20e-3,20e-3,nz,
                                -40e-3,40e-3,ny,
                                0.0,max_x,nx)



composite_min_z=-4e-3
composite_max_z=4e-3
composite_min_y=-10e-3; 
composite_max_y=10e-3


## XPS foam, approximate but pessimistic values
foam_k=.04 # W/m/deg K
foam_rho=40.0 # kg/m^3
foam_c=1500.0 # J/kg/deg K 

#spaceloft foam:
#foam_k=.016 # W/m/deg K
#foam_rho=150.0 # kg/m^3
#foam_c=1000.0 # J/kg/deg K 

materials=(
    # material 0: foam
    (heatsim2.TEMPERATURE_COMPUTE,foam_k,foam_rho,foam_c),
    # material 1: composite
    (heatsim2.TEMPERATURE_COMPUTE,composite_k,composite_rho,composite_c),
)

materials_nofoam=(
    # material 0: foam
    (heatsim2.TEMPERATURE_COMPUTE,0.0,foam_rho,foam_c),
    # material 1: composite
    (heatsim2.TEMPERATURE_COMPUTE,composite_k,composite_rho,composite_c),
)

boundaries=(
    # boundary 0: conducting
    (heatsim2.boundary_conducting,),
    (heatsim2.boundary_insulating,),
)

heat_energy=500e3  # Energy J/m^2

volumetric=(  # on material grid
    # 0: nothing
    (heatsim2.NO_SOURCE,),
    #1: impulse source @ t=0
    (heatsim2.IMPULSE_SOURCE,0.0,heat_energy/dx), # t (sec), Energy J/m^2
)

# initialize all elements to zero
(material_elements,
 boundary_z_elements,
 boundary_y_elements,
 boundary_x_elements,
 volumetric_elements)=heatsim2.zero_elements(nz,ny,nx) 

material_elements[ (zgrid >= composite_min_z) &
                   (zgrid <= composite_max_z) &
                   (ygrid >= composite_min_y) &
                   (ygrid <= composite_max_y)]=1 # set composite material
boundary_x_elements[:,:,0]=1 # insulating
boundary_x_elements[:,:,-1]=1 # insulating
boundary_y_elements[:,0,:]=1 # insulating
boundary_y_elements[:,-1,:]=1 # insulating
boundary_z_elements[0,:,:]=1 # insulating
boundary_z_elements[-1,:,:]=1 # insulating


# The no-foam version converges much better (and with many fewer
# time/spatial steps) if we make the insulation around the
# composite an explicit boundary condition
#
# (try commenting out the _nofoam boundary changes below to see what I mean)
boundary_z_elements_nofoam=copy.copy(boundary_z_elements)
boundary_y_elements_nofoam=copy.copy(boundary_y_elements)
boundary_x_elements_nofoam=copy.copy(boundary_x_elements)

boundary_y_elements_nofoam[((y_bnd_y==y_bnd[np.argmin(np.abs(y_bnd-composite_min_y))])
                            | (y_bnd_y==y_bnd[np.argmin(np.abs(y_bnd-composite_max_y))])) &
                           (y_bnd_z >= composite_min_z) &
                           (y_bnd_z <= composite_max_z)]=1 # insulating

boundary_z_elements_nofoam[((z_bnd_z==z_bnd[np.argmin(np.abs(z_bnd-composite_min_z))])
                            | (z_bnd_z==z_bnd[np.argmin(np.abs(z_bnd-composite_max_z))])) &
                           (z_bnd_y >= composite_min_y) &
                           (z_bnd_y <= composite_max_y)]=1 # insulating



volumetric_elements[(xgrid==x[0]) &
                    (ygrid >= composite_min_y) &
                    (ygrid <= composite_max_y) &
                    (zgrid >= composite_min_z) &
                    (zgrid <= composite_max_z)]=1  # impulse

#t0=-0.001
t0=0.0 
dt=0.05 # must be no bigger than .5 if we don't have the special nofoam boundaries, above
nt=18000

tvec=t0+numpy.arange(nt,dtype='d')*dt


(ADI_params,ADI_steps)=heatsim2.setup(z[0],y[0],x[0],
                                      dz,dy,dx,
                                      nz,ny,nx,
                                      dt,
                                      materials,
                                      boundaries,
                                      volumetric,
                                      material_elements,
                                      boundary_z_elements,
                                      boundary_y_elements,
                                      boundary_x_elements,
                                      volumetric_elements)

(ADI_params_nofoam,ADI_steps_nofoam)=heatsim2.setup(z[0],y[0],x[0],
                                                    dz,dy,dx,
                                                    nz,ny,nx,
                                                    dt,
                                                    materials_nofoam,
                                                    boundaries,
                                                    volumetric,
                                                    material_elements,
                                                    boundary_z_elements_nofoam,
                                                    boundary_y_elements_nofoam,
                                                    boundary_x_elements_nofoam,
                                                    volumetric_elements)

#               t0,dt,nt)


T=np.zeros((nt+1,nz,ny,nx),dtype='d')
T_nofoam=np.zeros((nt+1,nz,ny,nx),dtype='d')
for tcnt in range(nt):
    t=t0+dt*tcnt
    print "t=%f" % (t)
    T[tcnt+1,::]=heatsim2.run_adi_steps(ADI_params,ADI_steps,t,dt,T[tcnt,::],volumetric_elements,volumetric)
    T_nofoam[tcnt+1,::]=heatsim2.run_adi_steps(ADI_params_nofoam,ADI_steps_nofoam,t,dt,T_nofoam[tcnt,::],volumetric_elements,volumetric)
    pass

pl.figure(1)
pl.clf()
pl.imshow(T[nt,nz//2,::]);
pl.colorbar()

pl.figure(2)
pl.clf()

#maxT=np.max(T[1:,nz//2,ny//2,nx-1])
maxTidx=np.argmax(T[1:,nz//2,ny//2,nx-1])
if 2*maxTidx >= T.shape[0]:
    raise ValueError("Not enough timesteps to capture 2x peak")

maxTtime=(tvec+dt/2)[maxTidx]

#halfmaxTidx=np.argmin(abs(T[1:,nz//2,ny//2,nx-1]-maxT/2.0))
#halfmaxTtime=(tvec+dt/2)[halfmaxTidx]

#maxT_nofoam=np.max(T_nofoam[1:,nz//2,ny//2,nx-1])
#halfmaxTidx_nofoam=np.argmin(abs(T_nofoam[1:,nz//2,ny//2,nx-1]-maxT_nofoam/2.0))
#halfmaxTtime_nofoam=(tvec+dt/2)[halfmaxTidx_nofoam]

## Projection trick: 
#MaxT_time=(tvec+dt/2)[np.argmax(T[1:,nz//2,ny//2,nx-1])]
#MaxT_doubletimeidx=np.argmin(abs(MaxT_time*2-(tvec+dt/2)))
#MaxT_onefivetimeidx=np.argmin(abs(MaxT_time*1.5-(tvec+dt/2)))
#MaxTslope=(T[MaxT_doubletimeidx+1,nz//2,ny//2,nx-1]-T[MaxT_onefivetimeidx+1,nz//2,ny//2,nx-1])/((tvec+dt/2)[MaxT_doubletimeidx]-(tvec+dt/2)[MaxT_onefivetimeidx])
#MaxT_projected=maxT+MaxTslope*(-MaxT_time)

#MaxT_projected=T[1+MaxT_doubletimeidx,nz//2,ny//2,nx-1]+MaxTslope*(-(tvec+dt/2)[MaxT_doubletimeidx])
#halfmaxT_projectedidx=np.argmin(abs(T[1:,nz//2,ny//2,nx-1]-MaxT_projected/2.0))
#halfmaxT_projectedtime=(tvec+dt/2)[halfmaxT_projectedidx]
#Tinterpolate=sp.interpolate.splrep((tvec+dt/2),T[1:,nz//2,ny//2,nx-1]-MaxT_projected/2.0,s=0)
#halfmaxT_projectedtime=sp.interpolate.sproot(Tinterpolate)

# Inflection point trick (Ringermacher)
# Use data up to twice time to reach peak, as in ktester v2 implementation
tsr_num_predefknots=22  # same as in tsr_analysis.h for thermalconductivitytester. This number includes one on either end that we don't provide to splrep()
knotinterval=dt*(maxTidx*2-1-1)/(tsr_num_predefknots-1)
Tinterpolate=sp.interpolate.splrep((tvec+dt/2)[:(maxTidx*2-1)],T[1:(maxTidx*2),nz//2,ny//2,nx-1],k=5,task=-1,t=np.arange(tsr_num_predefknots-2,dtype='d')*knotinterval+knotinterval)
Tinterpolate_nofoam=sp.interpolate.splrep((tvec+dt/2)[:(maxTidx*2-1)],T_nofoam[1:(maxTidx*2),nz//2,ny//2,nx-1],k=5,task=-1,t=np.arange(tsr_num_predefknots-2,dtype='d')*knotinterval+knotinterval)

derivarray=np.array(sp.interpolate.spalde((tvec+dt/2)[:(maxTidx*2-1)],Tinterpolate),dtype='d')
derivarray_nofoam=np.array(sp.interpolate.spalde((tvec+dt/2)[:(maxTidx*2-1)],Tinterpolate_nofoam),dtype='d')
numzoompts=1000


inflpt_roughidx=np.argmax(derivarray[:,1])
# tzoom goes from the time of inflpt_roughidx-1 to the time of inflpt_roughidx+1
tzoom=np.arange(numzoompts,dtype='d')*dt*2/numzoompts+(tvec+dt/2)[inflpt_roughidx-1]
derivarrayzoom=np.array(sp.interpolate.spalde(tzoom,Tinterpolate),dtype='d')
inflpt_zoomidx=np.argmax(derivarrayzoom[:,1])
inflpt=tzoom[inflpt_zoomidx]
inflpt_T=derivarrayzoom[inflpt_zoomidx,0]
ringermacher_breaktime=(inflpt/.9055)


inflpt_roughidx_nofoam=np.argmax(derivarray_nofoam[:,1])
# tzoom goes from the time of inflpt_roughidx-1 to the time of inflpt_roughidx+1
tzoom_nofoam=np.arange(numzoompts,dtype='d')*dt*2/numzoompts+(tvec+dt/2)[inflpt_roughidx_nofoam-1]
derivarrayzoom_nofoam=np.array(sp.interpolate.spalde(tzoom_nofoam,Tinterpolate_nofoam),dtype='d')
inflpt_zoomidx_nofoam=np.argmax(derivarrayzoom_nofoam[:,1])
inflpt_nofoam=tzoom_nofoam[inflpt_zoomidx_nofoam]
ringermacher_breaktime_nofoam=(inflpt_nofoam/.9055)



#deriv=np.diff(T[:,nz//2,ny//2,nx-1])*dt
#deriv_nofoam=np.diff(T_nofoam[:,nz//2,ny//2,nx-1])*dt
#derivt=tvec+dt/2
#inflpt=derivt[np.argmax(deriv)]
#inflpt_nofoam=derivt[np.argmax(deriv_nofoam)]

# factor of .13879*pi^2 here converts from Ringermacher's
# Break time to ASTM E1461 half-max time
#ringermacher_halfmaxt_time=(inflpt/.9055)*(.13879*np.pi**2.0) # Ringermacher, 1998 http://lib.dr.iastate.edu/qnde/1998/allcontent/54/
#ringermacher_halfmaxt_time_nofoam=(inflpt_nofoam/.9055)*(.13879*np.pi**2.0) # Ringermacher, 1998 http://lib.dr.iastate.edu/qnde/1998/allcontent/54/


# Green's function solution: 
# 1/sqrt(t) * [2*exp(-l^2/(4*alpha t)) + 2*exp(-(3*l)^2/(4/alpha t)]

pl.plot(tvec+dt/2,T[1:,nz//2,ny//2,nx-1],'-',
        tvec+dt/2,T_nofoam[1:,nz//2,ny//2,nx-1],'-.',
    tvec+dt/2,((heat_energy/(composite_rho*composite_c))/np.sqrt(np.pi*(composite_k/(composite_rho*composite_c))*(tvec+dt/2)))*(2*np.exp(-((x_bnd[-1])**2)/(4.0*(composite_k/(composite_rho*composite_c))*(tvec+dt/2))) +2*np.exp(-((3*x_bnd[-1])**2)/(4*(composite_k/(composite_rho*composite_c))*(tvec+dt/2))) +2*np.exp(-((5*x_bnd[-1])**2)/(4*(composite_k/(composite_rho*composite_c))*(tvec+dt/2))) ),'--',
        inflpt,inflpt_T,'x',
        linewidth=3.0,markeredgewidth=3.0,markersize=10.0);
pl.legend(['With foam: break time=%f s' % (ringermacher_breaktime),'Without foam: break time=%f s' % (ringermacher_breaktime_nofoam),'theoretical'])
pl.axis([0,maxTtime*2,.02,20.0]);
pl.grid()

print("Error: %f%%" % ((ringermacher_breaktime-ringermacher_breaktime_nofoam)*100.0/ringermacher_breaktime_nofoam))

pl.show()
