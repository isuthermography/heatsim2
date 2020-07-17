import heatsim2
import numpy as np
import pylab as pl


# Create x,y,z voxel center coords
nz=5
ny=50
nx=55

(dz,dy,dx,
 z,y,x,
 zgrid,ygrid,xgrid,
 z_bnd,y_bnd,x_bnd,
 z_bnd_z,z_bnd_y,z_bnd_x,
 y_bnd_z,y_bnd_y,y_bnd_x,
 x_bnd_z,x_bnd_y,x_bnd_x,
 r3d,r2d) = heatsim2.build_grid(0,1.4e-3,nz,
                                -1.6e-3,1.6e-3,ny,
                                -1.5e-3,1.5e-3,nx)




# define materials:
#titanium_k=6.7 # W/m/deg K
titanium_k=np.diag(np.array((6.7,6.7,6.7),dtype='d'))
titanium_rho=4.43e3 # W/m/deg K
titanium_c=526.3 # J/kg/deg K

materials=(
    # material 0: titanium
    (heatsim2.TEMPERATURE_COMPUTE,titanium_k,titanium_rho,titanium_c),
    # material 1: Titanium, fixed temperature (Dirichlet boundary condition)
    (heatsim2.TEMPERATURE_FIXED,),
)

boundaries=(
    # boundary 0: conducting
    #(heatsim2.boundary_conducting,),
    (heatsim2.boundary_conducting_anisotropic,),
    (heatsim2.boundary_insulating,),
)

volumetric=(  # on material grid
    # 0: nothing
    (heatsim2.NO_SOURCE,),
    #1: impulse source @ t=0
    (heatsim2.IMPULSE_SOURCE,0.0,10e3/dz), # t (sec), Energy J/m^2
)

# initialize all elements to zero
(material_elements,
 boundary_z_elements,
 boundary_y_elements,
 boundary_x_elements,
 volumetric_elements)=heatsim2.zero_elements(nz,ny,nx) 

# set edges to insulating
boundary_x_elements[:,:,0]=1 # insulating
boundary_x_elements[:,:,-1]=1 # insulating
boundary_y_elements[:,0,:]=1 # insulating
boundary_y_elements[:,-1,:]=1 # insulating
boundary_z_elements[0,:,:]=1 # insulating
boundary_z_elements[-1,:,:]=1 # insulating

delam_depth=0.4e-3; 
delam_min_x=-1e-3
delam_max_x=1e-3
delam_min_y=-.6e-3
delam_max_y=.6e-3

boundary_z_elements[ (z_bnd_x > delam_min_x) &
                     (z_bnd_x < delam_max_x) &
                     (z_bnd_y > delam_min_y) &
                     (z_bnd_y < delam_max_y) &
                     (z_bnd_z==z_bnd[(np.abs(z_bnd-delam_depth)).argmin()])]=1

volumetric_elements[0,:,:]=1  # impulse

t0=0
dt=0.010
nt=200

tvec=t0+np.arange(nt,dtype='d')*dt


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



T=np.zeros((nt,nz,ny,nx),dtype='d')
curtemp=np.zeros((nz,ny,nx),dtype='d');
for tcnt in range(nt):
    t=t0+dt*tcnt
    print("t={}".format(t))
    curtemp=heatsim2.run_adi_steps(ADI_params,ADI_steps,t,dt,curtemp,volumetric_elements,volumetric)
    T[tcnt,::]=curtemp
    pass

pl.figure(1)
pl.clf()
pl.imshow(T[20,0,:,:],origin='lower',extent=(x_bnd[0]*1e3,x_bnd[-1]*1e3,y_bnd[0]*1e3,y_bnd[-1]*1e3))
pl.ylabel('Y position (mm)')
pl.xlabel('X position (mm)')
pl.title('Surface temperature increase (deg. C) at t=%f' % (tvec[20]+dt))
pl.colorbar()

pl.figure(2)
pl.clf()
pl.loglog(tvec+dt/2.0,T[:,0,20,21])
pl.grid()
pl.xlabel('Time (s)')
pl.ylabel('Surface temp. increase (deg. C)')
pl.title('Temperature vs. time at x=%f mm, y=%f mm' % (x[21]*1e3,y[20]*1e3))

pl.show()
