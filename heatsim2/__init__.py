import numpy as np
# import hs2_indexing


import heatsim2.boundary_conducting as boundary_conducting
import heatsim2.boundary_conducting_anisotropic as boundary_conducting_anisotropic
import heatsim2.boundary_insulating as boundary_insulating
import heatsim2.boundary_thininsulatinglayer as boundary_thininsulatinglayer


#import heatsim2.expression
#import heatsim2.alternatingdirection
import heatsim2.alternatingdirection_c_pyx
import heatsim2.alternatingdirection_c_pyx as alternatingdirection
import heatsim2.crank_nicolson

setup=heatsim2.crank_nicolson.setup
#run_adi_steps=heatsim2.alternatingdirection.run_adi_steps
run_adi_steps=heatsim2.alternatingdirection_c_pyx.run_adi_steps

# defines for material[0]
TEMPERATURE_COMPUTE=0
TEMPERATURE_FIXED=1

# defines for volumetric[0]:
NO_SOURCE=0
IMPULSE_SOURCE=1
STEPPED_SOURCE=2
IMPULSE_POINT_SOURCE_JOULES=3
SPATIALLY_Z_DECAYING_TEMPORAL_IMPULSE=4


def build_grid(minz,maxz,nz,miny,maxy,ny,minx,maxx,nx):
    # build a grid with BOUNDARIES at minz, maxz, etc.
    # ... return ELEMENT CENTER COORDS (x,y,z)

    # Create x, y, and z boundary coordinates
    (z_bnd,dz)=np.linspace(minz,maxz,num=nz+1,retstep=True)
    (y_bnd,dy)=np.linspace(miny,maxy,num=ny+1,retstep=True)
    (x_bnd,dx)=np.linspace(minx,maxx,num=nx+1,retstep=True)
    
    # Create x,y,z element center grid
    z=z_bnd[:-1]+dz/2.0
    y=y_bnd[:-1]+dy/2.0
    x=x_bnd[:-1]+dx/2.0

    # Create 3d meshgrids indicating z boundary location
    # for all x,y center positions

    # Voxel at i, j, k is has x boundarys at x_bnd[i] and x_bnd[i+1],
    # centered at y[k],z[k]
    # Same voxel has y boundaries at y_bnd[j] and y_bnd[j+1] which
    # are centered at x=x[i] and z=z[k]

    (z_bnd_z,z_bnd_y,z_bnd_x)=np.meshgrid(z_bnd,y,x,indexing='ij')

    (y_bnd_z,y_bnd_y,y_bnd_x)=np.meshgrid(z,y_bnd,x,indexing='ij')

    (x_bnd_z,x_bnd_y,x_bnd_x)=np.meshgrid(z,y,x_bnd,indexing='ij')

    # create 3d meshgrids indicating element centers
    (zgrid,ygrid,xgrid) = np.meshgrid(z,y,x,indexing='ij')

    # create 3D radius from (0,0,0):
    r3d=np.sqrt(xgrid**2+ygrid**2+zgrid**2)
    
    # create 2D radius from (0,0):
    r2d=np.sqrt(xgrid**2+ygrid**2)


    
    return (dz,dy,dx,
            z,y,x,
            zgrid,ygrid,xgrid,
            z_bnd,y_bnd,x_bnd,
            z_bnd_z,z_bnd_y,z_bnd_x,
            y_bnd_z,y_bnd_y,y_bnd_x,
            x_bnd_z,x_bnd_y,x_bnd_x,
            r3d,r2d)


def build_grid_min_step(minzcenter,dz,nz,minycenter,dy,ny,minxcenter,dx,nx):
    z_bnd=minzcenter-dz/2.0 + np.arange(nz+1,dtype='d')*dz
    y_bnd=minycenter-dy/2.0 + np.arange(ny+1,dtype='d')*dy
    x_bnd=minxcenter-dx/2.0 + np.arange(nx+1,dtype='d')*dx
        
    # Create x,y,z element center grid
    z=z_bnd[:-1]+dz/2.0
    y=y_bnd[:-1]+dy/2.0
    x=x_bnd[:-1]+dx/2.0

    # Create 3d meshgrids indicating z boundary location
    # for all x,y center positions

    # Voxel at i, j, k is has x boundarys at x_bnd[i] and x_bnd[i+1],
    # centered at y[k],z[k]
    # Same voxel has y boundaries at y_bnd[j] and y_bnd[j+1] which
    # are centered at x=x[i] and z=z[k]

    (z_bnd_z,z_bnd_y,z_bnd_x)=np.meshgrid(z_bnd,y,x,indexing='ij')

    (y_bnd_z,y_bnd_y,y_bnd_x)=np.meshgrid(z,y_bnd,x,indexing='ij')

    (x_bnd_z,x_bnd_y,x_bnd_x)=np.meshgrid(z,y,x_bnd,indexing='ij')

    # create 3d meshgrids indicating element centers
    (zgrid,ygrid,xgrid) = np.meshgrid(z,y,x,indexing='ij')

    # create 3D radius from (0,0,0):
    r3d=np.sqrt(xgrid**2+ygrid**2+zgrid**2)
    
    # create 2D radius from (0,0):
    r2d=np.sqrt(xgrid**2+ygrid**2)

    return (z,y,x,
            zgrid,ygrid,xgrid,
            z_bnd,y_bnd,x_bnd,
            z_bnd_z,z_bnd_y,z_bnd_x,
            y_bnd_z,y_bnd_y,y_bnd_x,
            x_bnd_z,x_bnd_y,x_bnd_x,
            r3d,r2d)




def build_grid_min_step_edge(minzedge,dz,nz,minyedge,dy,ny,minxedge,dx,nx):
    z_bnd=minzedge + np.arange(nz+1,dtype='d')*dz
    y_bnd=minyedge + np.arange(ny+1,dtype='d')*dy
    x_bnd=minxedge + np.arange(nx+1,dtype='d')*dx
        
    # Create x,y,z element center grid
    z=z_bnd[:-1]+dz/2.0
    y=y_bnd[:-1]+dy/2.0
    x=x_bnd[:-1]+dx/2.0

    # Create 3d meshgrids indicating z boundary location
    # for all x,y center positions

    # Voxel at i, j, k is has x boundarys at x_bnd[i] and x_bnd[i+1],
    # centered at y[k],z[k]
    # Same voxel has y boundaries at y_bnd[j] and y_bnd[j+1] which
    # are centered at x=x[i] and z=z[k]

    (z_bnd_z,z_bnd_y,z_bnd_x)=np.meshgrid(z_bnd,y,x,indexing='ij')

    (y_bnd_z,y_bnd_y,y_bnd_x)=np.meshgrid(z,y_bnd,x,indexing='ij')

    (x_bnd_z,x_bnd_y,x_bnd_x)=np.meshgrid(z,y,x_bnd,indexing='ij')

    # create 3d meshgrids indicating element centers
    (zgrid,ygrid,xgrid) = np.meshgrid(z,y,x,indexing='ij')

    # create 3D radius from (0,0,0):
    r3d=np.sqrt(xgrid**2+ygrid**2+zgrid**2)
    
    # create 2D radius from (0,0):
    r2d=np.sqrt(xgrid**2+ygrid**2)

    return (z,y,x,
            zgrid,ygrid,xgrid,
            z_bnd,y_bnd,x_bnd,
            z_bnd_z,z_bnd_y,z_bnd_x,
            y_bnd_z,y_bnd_y,y_bnd_x,
            x_bnd_z,x_bnd_y,x_bnd_x,
            r3d,r2d)



            
    
def zero_elements(nz,ny,nx):
    material_elements=np.zeros((nz,ny,nx),dtype='u1')
    boundary_z_elements=np.zeros((nz+1,ny,nx),dtype='u1')
    boundary_y_elements=np.zeros((nz,ny+1,nx),dtype='u1')
    boundary_x_elements=np.zeros((nz,ny,nx+1),dtype='u1')
    volumetric_elements=np.zeros((nz,ny,nx),dtype='u1')
    return (material_elements,boundary_z_elements,boundary_y_elements,boundary_x_elements,volumetric_elements)
