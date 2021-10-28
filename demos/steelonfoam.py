import heatsim2
import numpy as np
numpy=np
import pylab as pl

# 1/8" steel on foam with a small insulating delamination gap


# Create x,y,z voxel center coords
nx=48
ny=40
nz=80

measx=5e-3
measy=5e-3
meas2x=50e-3
meas2y=45e-3

#t0=-0.001
t0=0
dt=0.01
nt=1000 # was 100000

tvec=t0+numpy.arange(nt,dtype='d')*dt


(dz,dy,dx,
 z,y,x,
 zgrid,ygrid,xgrid,
 z_bnd,y_bnd,x_bnd,
 z_bnd_z,z_bnd_y,z_bnd_x,
 y_bnd_z,y_bnd_y,y_bnd_x,
 x_bnd_z,x_bnd_y,x_bnd_x,
 r3d,r2d) = heatsim2.build_grid(0,10.8e-3,nz,
                                -0.05,0.05,ny,
                                -.06,.06,nx)

(ygrid2d,xgrid2d) = np.meshgrid(y,x,indexing='ij')

measi=np.argmin(abs(x-measx))
measj=np.argmin(abs(y-measy))

meas2i=np.argmin(abs(x-meas2x))
meas2j=np.argmin(abs(y-meas2y))

foam_depth=3.175e-3; # distance of foam from flash surface

foam_interface_index = np.argmin(np.abs(foam_depth-z_bnd))

gap_min_x=0e-3
gap_max_x=10e-3
gap_min_y=0e-3
gap_max_y=10e-3

# define materials: (very rough values for steel and foam)
steel_k=40.0 # W/m/deg K   
steel_rho=7.75e3 # W/m/deg K
steel_c=466.0 # J/kg/deg K

# guess values for foam
foam_k=.4  # W/m/deg K   
foam_rho=40.0 # kg/m^3 K
foam_c=1500.0 # J/kg/deg K

materials=(
    # material 0: steel
    (heatsim2.TEMPERATURE_COMPUTE,steel_k,steel_rho,steel_c),
    # material 1: Foam 
    (heatsim2.TEMPERATURE_COMPUTE,foam_k,foam_rho,foam_c),
    # material 2: Fixed temperature (Dirichlet boundary condition)
    (heatsim2.TEMPERATURE_FIXED,),
)

boundaries=(
    # boundary 0: conducting
    (heatsim2.boundary_conducting,),
    # boundary 1: insulating
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

# place foam in on back surface beyond the foam interface
material_elements[foam_interface_index:,:,:]=1

# set x and y and z edges to insulating
boundary_x_elements[:,:,0]=1 # insulating
boundary_x_elements[:,:,-1]=1 # insulating
boundary_y_elements[:,0,:]=1 # insulating
boundary_y_elements[:,-1,:]=1 # insulating
boundary_z_elements[0,:,:]=1 # insulating
boundary_z_elements[-1,:,:]=1 # insulating


volumetric_elements[0,:,:]=1  # impulse on steel surface

# define delamination gap as insulating b.c.

boundary_z_elements[foam_interface_index,:,:][(ygrid2d > gap_min_y) & (ygrid2d < gap_max_y) & (xgrid2d > gap_min_x) & (xgrid2d < gap_max_x)]=1 # insulating


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
    print("t={}".format(t))
    T[tcnt+1,::]=heatsim2.run_adi_steps(ADI_params,ADI_steps,t,dt,T[tcnt,::],volumetric_elements,volumetric)
    pass


pl.figure(1)
pl.clf()
pl.loglog(tvec[1:]-dt/2,T[1:,0,measj,measi],'-',
          tvec[1:]-dt/2,T[1:,0,meas2j,meas2i],'-',)
pl.grid()
pl.legend(['Temperature atop gap','Temperature away from gap'])

pl.show()

