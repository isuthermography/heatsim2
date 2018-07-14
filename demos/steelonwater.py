import heatsim2
import cProfile
import numpy as np
numpy=np
import pylab as pl


# Create x,y,z voxel center coords
nx=12*2
ny=10*2
nz=80

measx=10e-3
measy=10e-3
meas2x=70e-3
meas2y=70e-3

#t0=-0.001
t0=0
dt=0.05
nt=1000 # was 100000

tvec=t0+numpy.arange(nt,dtype='d')*dt


(dz,dy,dx,
 z,y,x,
 zgrid,ygrid,xgrid,
 z_bnd,y_bnd,x_bnd,
 z_bnd_z,z_bnd_y,z_bnd_x,
 y_bnd_z,y_bnd_y,y_bnd_x,
 x_bnd_z,x_bnd_y,x_bnd_x,
 r3d,r2d) = heatsim2.build_grid(0,50.8e-3,nz,
                                -0.1,0.1,ny,
                                -.12,.12,nx)


measi=np.argmin(abs(x-measx))
measj=np.argmin(abs(y-measy))

meas2i=np.argmin(abs(x-meas2x))
meas2j=np.argmin(abs(y-meas2y))

hole_depth=9.525e-3 #6.35e-3 #9.525e-3;   # distance of hole from flash surface
water_depth=12.7e-3; # distance of water volume from flash surface
hole_min_x=0e-3
hole_max_x=10e-3
hole_min_y=0e-3
hole_max_y=10e-3

# define materials: (very rough values for steel and water)
steel_k=40.0 # W/m/deg K   
steel_rho=7.75e3 # W/m/deg K
steel_c=466.0 # J/kg/deg K

water_k= 20.0 #0.6 # W/m/deg K   
water_rho=1.00e3 # W/m/deg K
water_c=418 # J/kg/deg K

materials=(
    # material 0: steel
    (heatsim2.TEMPERATURE_COMPUTE,steel_k,steel_rho,steel_c),
    # material 1: Water 
    (heatsim2.TEMPERATURE_COMPUTE,water_k,water_rho,water_c),
    # material 2: Fixed temperature (Dirichlet boundary condition)
    (heatsim2.TEMPERATURE_FIXED,),
)

boundaries=(
    # boundary 0: conducting
    (heatsim2.boundary_conducting,),
    (heatsim2.boundary_insulating,),
)

volumetric=(  # on material grid
    # 0: nothing
    (heatsim2.NO_SOURCE,),
    ##1: impulse source @ t=0
    (heatsim2.IMPULSE_SOURCE,0.0,10e3/dz), # t (sec), Energy J/m^2 
    #1: stepped source 
    #(heatsim2.STEPPED_SOURCE,0.0-dt/2.0,1.0-dt/2.0,10e3/dz), # t0 (sec), t1 (sec), Power W/m^2 
)

# initialize all elements to zero
(material_elements,
 boundary_z_elements,
 boundary_y_elements,
 boundary_x_elements,
 volumetric_elements)=heatsim2.zero_elements(nz,ny,nx) 

# place water in hole and on back surface
material_elements[ (xgrid > hole_min_x) &
                     (xgrid < hole_max_x) &
                     (ygrid > hole_min_y) &
                     (ygrid < hole_max_y) &
                     (zgrid > hole_depth)]=1

# set x and y and z=0 edges to insulating
boundary_x_elements[:,:,0]=1 # insulating
boundary_x_elements[:,:,-1]=1 # insulating
boundary_y_elements[:,0,:]=1 # insulating
boundary_y_elements[:,-1,:]=1 # insulating
boundary_z_elements[0,:,:]=1 # insulating

# zero temperature BC at back end of model
boundary_z_elements[-1,:,:]=1 # insulating
material_elements[-1,:,:]=2  # fixed temperature



volumetric_elements[0,:,:]=1  # impulse


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

#               t0,dt,nt)


T=np.zeros((nt,nz,ny,nx),dtype='d')
for tcnt in range(nt-1):
    t=t0+dt*tcnt
    print "t=%f" % (t)
    T[tcnt+1,::]=heatsim2.run_adi_steps(ADI_params,ADI_steps,t,dt,T[tcnt,::],volumetric_elements,volumetric)
    pass


pl.figure(1)
pl.clf()
pl.loglog(tvec[1:]-dt/2,T[1:,0,measj,measi],'-',
          tvec[1:]-dt/2,T[1:,0,meas2j,meas2i],'-',)
pl.grid()
pl.legend(['Temperature atop hole','Temperature away from hole'])

pl.show()

