import sys
import os
import getpass


import numpy as np
numpy=np
import pylab as pl
import copy

import datetime

# Fix missing timedelta.total_seconds in Python < 2.7, based on http://stackoverflow.com/questions/3318348/how-can-i-extend-pythons-datetime-datetime-with-my-own-methods/14214646#14214646
import ctypes as c

_get_dict = c.pythonapi._PyObject_GetDictPtr
_get_dict.restype = c.POINTER(c.py_object)
_get_dict.argtypes = [c.py_object]

from datetime import timedelta
try:
    timedelta.total_seconds # new in 2.7
except AttributeError:
    def total_seconds(td):
        return float((td.microseconds +
                      (td.seconds + td.days * 24 * 3600) * 10**6)) / 10**6
    d = _get_dict(timedelta)[0]
    d['total_seconds'] = total_seconds

import heatsim2
import heatsim


run_heatsim2=True
run_heatsim=True

run_coarse=True
run_fine=True
run_fine2=True

load_comsol=False #True

enable_barrier=True

flash_energy=10e3 # J/m^2

# Create x,y,z voxel center coords
nz=7
ny=32
nx=30

nz_fine=nz*3
ny_fine=ny*2
nx_fine=nx*2

measx=0.4e-3
measy=0.3e-3

z_thick=1.4e-3 # m

# Flash heating:


# Temperature after first frame should be:
#flash_energy/(composite_rho*composite_c*dz)

# final temperature should be flash_energy/(composite_rho*composite_c*z_thick)


(dz,dy,dx,
 z,y,x,
 zgrid,ygrid,xgrid,
 z_bnd,y_bnd,x_bnd,
 z_bnd_z,z_bnd_y,z_bnd_x,
 y_bnd_z,y_bnd_y,y_bnd_x,
 x_bnd_z,x_bnd_y,x_bnd_x,
 r3d,r2d) = heatsim2.build_grid(0,z_thick,nz,
                                -1.6e-3,1.6e-3,ny,
                                -1.5e-3,1.5e-3,nx)

measi=np.argmin(abs(x-measx))
measj=np.argmin(abs(y-measy))

(dz_fine,dy_fine,dx_fine,
 z_fine,y_fine,x_fine,
 zgrid_fine,ygrid_fine,xgrid_fine,
 z_bnd_fine,y_bnd_fine,x_bnd_fine,
 z_bnd_z_fine,z_bnd_y_fine,z_bnd_x_fine,
 y_bnd_z_fine,y_bnd_y_fine,y_bnd_x_fine,
 x_bnd_z_fine,x_bnd_y_fine,x_bnd_x_fine,
 r3d_fine,r2d_fine) = heatsim2.build_grid(0,z_thick,nz_fine,
                                          y_bnd[0],y_bnd[-1],ny_fine,
                                          x_bnd[0],x_bnd[-1],nx_fine)

measi_fine=np.argmin(abs(x_fine-measx))
measj_fine=np.argmin(abs(y_fine-measy))


barrier_min_z=0.4e-3; 
barrier_max_z=.8e-3; 
barrier_min_x=-1e-3
barrier_max_x=1e-3
barrier_min_y=-.6e-3
barrier_max_y=.6e-3

# define materials:
composite_k=.138 # W/m/deg K
composite_rho=1.57e3 # W/m/deg K
composite_c=730 # J/kg/deg K

barrier_k=.05 # W/m/deg K
barrier_rho=1.57e3 # W/m/deg K
barrier_c=730 # J/kg/deg K

materials=(
    # material 0: composite
    (heatsim2.TEMPERATURE_COMPUTE,composite_k,composite_rho,composite_c),
    # material 1: barrier
    (heatsim2.TEMPERATURE_COMPUTE,barrier_k,barrier_rho,barrier_c),
    # material 2: composite  (so we can use same material matrix as heatsim)
    (heatsim2.TEMPERATURE_COMPUTE,composite_k,composite_rho,composite_c),
    # material 3: Fixed temperature (Dirichlet boundary condition)
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
    #1: impulse source @ t=0
    (heatsim2.IMPULSE_SOURCE,0.0,flash_energy/dz), # t (sec), Energy J/m^2
)

volumetric_fine=(  # on material grid
    # 0: nothing
    (heatsim2.NO_SOURCE,),
    #1: impulse source @ t=0
    (heatsim2.IMPULSE_SOURCE,0.0,flash_energy/dz_fine), # t (sec), Energy J/m^2
)



heatsim_materials=(
    # material 0: bulk medium
    (heatsim.BOUNDARY_NONE,composite_k,composite_rho,composite_c, 0.0,0.0,0.0,0.0),
    # material 1: barrier
    (heatsim.BOUNDARY_NONE,barrier_k,barrier_rho,barrier_c, 0.0,0.0,0.0,0.0),
    # material 2: flash source
    (heatsim.BOUNDARY_IMPULSE,composite_k,composite_rho,composite_c, 0.0,0.0,0.0,0.0),
    # material 3: insulator
    (heatsim.BOUNDARY_INSULATING,0.0,0.0,0.0, 0.0,0.0,0.0,0.0),
    )


# initialize all elements to zero
(material_elements,
 boundary_z_elements,
 boundary_y_elements,
 boundary_x_elements,
 volumetric_elements)=heatsim2.zero_elements(nz,ny,nx) 


(material_elements_fine,
 boundary_z_elements_fine,
 boundary_y_elements_fine,
 boundary_x_elements_fine,
 volumetric_elements_fine)=heatsim2.zero_elements(nz_fine,ny_fine,nx_fine) 


# define nonzero material elements
material_elements[0,:,:]=2 # material 2: flash source (for heatsim)
material_elements_fine[0,:,:]=2 # material 2: flash source (for heatsim)

if enable_barrier:
    material_elements[(zgrid >= barrier_min_z) &
                      (zgrid <= barrier_max_z) &
                      (ygrid >= barrier_min_y) &
                      (ygrid <= barrier_max_y) &
                      (xgrid >= barrier_min_x) &
                      (xgrid <= barrier_max_x)]=1 # material 1: barrier

    material_elements_fine[(zgrid_fine >= barrier_min_z) &
                           (zgrid_fine <= barrier_max_z) &
                           (ygrid_fine >= barrier_min_y) &
                           (ygrid_fine <= barrier_max_y) &
                           (xgrid_fine >= barrier_min_x) &
                      (xgrid_fine <= barrier_max_x)]=1 # material 1: barrier
    
    pass

volumetric_elements[0,:,:]=1  # set flash source (for heatsim2)
volumetric_elements_fine[0,:,:]=1  # set flash source (for heatsim2)


# set edges to insulating
boundary_x_elements[:,:,0]=1 # insulating
boundary_x_elements[:,:,-1]=1 # insulating
boundary_y_elements[:,0,:]=1 # insulating
boundary_y_elements[:,-1,:]=1 # insulating
boundary_z_elements[0,:,:]=1 # insulating
boundary_z_elements[-1,:,:]=1 # insulating

boundary_x_elements_fine[:,:,0]=1 # insulating
boundary_x_elements_fine[:,:,-1]=1 # insulating
boundary_y_elements_fine[:,0,:]=1 # insulating
boundary_y_elements_fine[:,-1,:]=1 # insulating
boundary_z_elements_fine[0,:,:]=1 # insulating
boundary_z_elements_fine[-1,:,:]=1 # insulating


t0=0.0
tf=20.0

nt_heatsim=2000
(t_heatsim,dt_heatsim)=np.linspace(t0,tf,num=nt_heatsim,retstep=True)

nt_heatsim_fine=20000
(t_heatsim_fine,dt_heatsim_fine)=np.linspace(t0,tf,num=nt_heatsim_fine,retstep=True)

nt_heatsim2=200
(t_heatsim2,dt_heatsim2)=np.linspace(t0,tf,num=nt_heatsim2,retstep=True)

nt_heatsim2_fine=1000
(t_heatsim2_fine,dt_heatsim2_fine)=np.linspace(t0,tf,num=nt_heatsim2_fine,retstep=True)

nt_heatsim2_fine2=20000
(t_heatsim2_fine2,dt_heatsim2_fine2)=np.linspace(t0,tf,num=nt_heatsim2_fine2,retstep=True)



if run_heatsim2 and run_coarse:
    hs2_start=datetime.datetime.now()
    (ADI_params,ADI_steps)=heatsim2.setup(z[0],y[0],x[0],
                                          dz,dy,dx,
                                          nz,ny,nx,
                                          dt_heatsim2,
                                          materials,
                                          boundaries,
                                          volumetric,
                                          material_elements,
                                          boundary_z_elements,
                                          boundary_y_elements,
                                          boundary_x_elements,
                                          volumetric_elements)


    T=np.zeros((nz,ny,nx),dtype='d')
    
    hs2_start_exec=datetime.datetime.now()
    heatsim2meas=np.zeros(nt_heatsim2,dtype='d')
    for tcnt in range(nt_heatsim2):
        curt=t0+dt_heatsim2*tcnt
        print "t=%f" % (curt)
        T=heatsim2.run_adi_steps(ADI_params,ADI_steps,curt,dt_heatsim2,T,volumetric_elements,volumetric)
        heatsim2meas[tcnt]=T[0,measj,measi] # note: heatsim array is transposed compared with new standard

        if tcnt==np.argmin(abs(t_heatsim2-2.0)):
            heatsim2_twosec=copy.copy(T)
            pass
        pass
    heatsim2res=T
    hs2_setuptime=(hs2_start_exec-hs2_start).total_seconds()
    hs2_comptime=(datetime.datetime.now()-hs2_start_exec).total_seconds()

    pass

if run_heatsim and run_coarse:
    hs_start=datetime.datetime.now()
    
    T=np.zeros((nz,ny,nx),dtype='d')    
    sim=heatsim.createsim(heatsim_materials,material_elements.transpose(),T.transpose(),x[0],dx,nx,y[0],dy,ny,z[0],dz,nz,t0,dt_heatsim)

    sim.bound_array[:,:,0]=flash_energy / (dz) ;

    heatsimmeas=np.zeros(nt_heatsim,dtype='d')

    hs_start_exec=datetime.datetime.now()

    for tcnt in range(nt_heatsim):
        curt=t0+dt_heatsim*tcnt
        if tcnt % 100 == 0:
            print("t=%f" % (curt))
            pass
        heatsim.stepforward(sim,curt)
        heatsimmeas[tcnt]=sim.array[measi,measj,0] # note: heatsim array is transposed compared with new standard
        if tcnt==np.argmin(abs(t_heatsim-2.0)):
            heatsim_twosec=copy.copy(sim.array.transpose())
            pass
        pass
    heatsimres=sim.array.transpose()
    hs_setuptime=(hs_start_exec-hs_start).total_seconds()
    hs_comptime=(datetime.datetime.now()-hs_start_exec).total_seconds()

    pass


if run_heatsim2 and run_fine: 
    hs2_start_fine=datetime.datetime.now()

    (ADI_params_fine,ADI_steps_fine)=heatsim2.setup(z_fine[0],y_fine[0],x_fine[0],
                                                    dz_fine,dy_fine,dx_fine,
                                                    nz_fine,ny_fine,nx_fine,
                                                    dt_heatsim2_fine,
                                                    materials,
                                                    boundaries,
                                                    volumetric_fine,
                                                    material_elements_fine,
                                                    boundary_z_elements_fine,
                                                    boundary_y_elements_fine,
                                                    boundary_x_elements_fine,
                                                    volumetric_elements_fine)


    T_fine=np.zeros((nz_fine,ny_fine,nx_fine),dtype='d')
    
    heatsim2meas_fine=np.zeros(nt_heatsim2_fine,dtype='d')

    hs2_start_exec_fine=datetime.datetime.now()

    for tcnt in range(nt_heatsim2_fine):
        curt=t0+dt_heatsim2_fine*tcnt
        print "t=%f" % (curt)
        T_fine=heatsim2.run_adi_steps(ADI_params_fine,ADI_steps_fine,curt,dt_heatsim2_fine,T_fine,volumetric_elements_fine,volumetric_fine)
        heatsim2meas_fine[tcnt]=T_fine[0,measj_fine,measi_fine] # note: heatsim array is transposed compared with new standard

        if tcnt==np.argmin(abs(t_heatsim2_fine-2.0)):
            heatsim2_twosec_fine=copy.copy(T_fine)
            pass
        pass
    heatsim2res_fine=T_fine

    hs2_setuptime_fine=(hs2_start_exec_fine-hs2_start_fine).total_seconds()
    hs2_comptime_fine=(datetime.datetime.now()-hs2_start_exec_fine).total_seconds()

    pass


if run_heatsim2 and run_fine2: 
    hs2_start_fine2=datetime.datetime.now()

    (ADI_params_fine2,ADI_steps_fine2)=heatsim2.setup(z_fine[0],y_fine[0],x_fine[0],
                                                    dz_fine,dy_fine,dx_fine,
                                                    nz_fine,ny_fine,nx_fine,
                                                    dt_heatsim2_fine2,
                                                    materials,
                                                    boundaries,
                                                    volumetric_fine,
                                                    material_elements_fine,
                                                    boundary_z_elements_fine,
                                                    boundary_y_elements_fine,
                                                    boundary_x_elements_fine,
                                                    volumetric_elements_fine)


    T_fine2=np.zeros((nz_fine,ny_fine,nx_fine),dtype='d')
    
    heatsim2meas_fine2=np.zeros(nt_heatsim2_fine2,dtype='d')

    hs2_start_exec_fine2=datetime.datetime.now()

    for tcnt in range(nt_heatsim2_fine2):
        curt=t0+dt_heatsim2_fine2*tcnt
        print "t=%f" % (curt)
        T_fine2=heatsim2.run_adi_steps(ADI_params_fine2,ADI_steps_fine2,curt,dt_heatsim2_fine2,T_fine2,volumetric_elements_fine,volumetric_fine)
        heatsim2meas_fine2[tcnt]=T_fine2[0,measj_fine,measi_fine] # note: heatsim array is transposed compared with new standard

        if tcnt==np.argmin(abs(t_heatsim2_fine2-2.0)):
            heatsim2_twosec_fine2=copy.copy(T_fine2)
            pass
        pass
    heatsim2res_fine2=T_fine2

    hs2_setuptime_fine2=(hs2_start_exec_fine2-hs2_start_fine2).total_seconds()
    hs2_comptime_fine2=(datetime.datetime.now()-hs2_start_exec_fine2).total_seconds()

    pass


if run_heatsim and run_fine:
    hs_start_fine=datetime.datetime.now()
    
    T_fine=np.zeros((nz_fine,ny_fine,nx_fine),dtype='d')    
    sim_fine=heatsim.createsim(heatsim_materials,material_elements_fine.transpose(),T_fine.transpose(),x_fine[0],dx_fine,nx_fine,y_fine[0],dy_fine,ny_fine,z_fine[0],dz_fine,nz_fine,t0,dt_heatsim_fine)

    sim_fine.bound_array[:,:,0]=flash_energy / (dz_fine) ;

    heatsimmeas_fine=np.zeros(nt_heatsim_fine,dtype='d')
    hs_start_exec_fine=datetime.datetime.now()
    
    for tcnt in range(nt_heatsim_fine):
        curt=t0+dt_heatsim_fine*tcnt
        if tcnt % 100 == 0:
            print("t=%f" % (curt))
            pass
        heatsim.stepforward(sim_fine,curt)
        heatsimmeas_fine[tcnt]=sim_fine.array[measi_fine,measj_fine,0] # note: heatsim array is transposed compared with new standard
        if tcnt==np.argmin(abs(t_heatsim_fine-2.0)):
            heatsim_twosec_fine=copy.copy(sim_fine.array.transpose())
            pass
        pass
    heatsimres_fime=sim_fine.array.transpose()

    hs_setuptime_fine=(hs_start_exec_fine-hs_start_fine).total_seconds()
    hs_comptime_fine=(datetime.datetime.now()-hs_start_exec_fine).total_seconds()
    pass

if load_comsol:
    comsoloutput=os.path.join('/tmp','heatsimcomsoloutput_%s.csv' % (getpass.getuser()))
    comsoldat=np.genfromtxt(comsoloutput, delimiter=',',comments='%')
    t_comsol=comsoldat[:,0][comsoldat[:,0] > 0.0]
    T_comsol=comsoldat[:,1][comsoldat[:,0] > 0.0]#-293.15;
    pass

if run_heatsim and run_fine:
    pl.figure(1)
    pl.clf()
    #pl.imshow(heatsimres[0,::])
    pl.imshow(heatsim_twosec_fine[0,::])
    pl.colorbar()
    pass

if run_heatsim2 and run_fine:
    pl.figure(2)
    pl.clf()
    #pl.imshow(heatsim2res[0,::])
    pl.imshow(heatsim2_twosec_fine[0,::])
    pl.colorbar()
    pass

pl.figure(3)
pl.clf()

# For the full-space green's function there would be a sqrt(4pi) in
# the denominator. Instead, for the half space we get double the
# temperature for the same energy because all the energy flows
# in one direction rather than being split among both,
# so the denominator just has the sqrt(pi)
# it doesn't equilibriate right
plotargs=[t_heatsim[1:],((flash_energy/(composite_rho*composite_c))/np.sqrt(np.pi*(composite_k/(composite_rho*composite_c))*t_heatsim[1:]))*(1+2*np.exp(-(2*z_thick)**2/(4*(composite_k/(composite_rho*composite_c))*t_heatsim[1:]))+2*np.exp(-(4*z_thick)**2/(4*(composite_k/(composite_rho*composite_c))*t_heatsim[1:]))),'--']
legendargs=['Theory']

if run_heatsim and run_coarse:
    plotargs.extend([t_heatsim[1:],heatsimmeas[1:],'o-'])
    legendargs.append('heatsim %f setup %f exec' % (hs_setuptime,hs_comptime))
    pass
if run_heatsim2 and run_coarse:
    plotargs.extend([t_heatsim2[1:],heatsim2meas[1:],'x:'])
    legendargs.append('heatsim2 %f setup %f exec' % (hs2_setuptime,hs2_comptime))
    pass


if run_heatsim and run_fine:
    plotargs.extend([t_heatsim_fine[1:],heatsimmeas_fine[1:],'o-'])
    legendargs.append('heatsim fine %f setup %f exec' % (hs_setuptime_fine,hs_comptime_fine))
    pass
if run_heatsim2 and run_fine:
    plotargs.extend([t_heatsim2_fine[1:],heatsim2meas_fine[1:],'x:'])
    legendargs.append('heatsim2 fine %f setup %f exec' % (hs2_setuptime_fine,hs2_comptime_fine))
    pass

if run_heatsim2 and run_fine2:
    plotargs.extend([t_heatsim2_fine2[1:],heatsim2meas_fine2[1:],'x:'])
    legendargs.append('heatsim2 fine2 %f setup %f exec' % (hs2_setuptime_fine2,hs2_comptime_fine2))
    pass

if load_comsol:
    plotargs.extend([t_comsol,T_comsol,'*-.'])
    legendargs.append('comsol')
    pass

pl.loglog(*plotargs)
pl.legend(legendargs)
pl.xlabel('Time (s)')
pl.ylabel('Temperature (K)')
pl.grid()
pl.show()
