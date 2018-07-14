
# estimate the surface temperature at

def insulating_z_min_surface_temperature(T,dz):
    # T should be T[z,y,x]
    # Estimates temperature of z=0surface just prior to T[0,::]
    # based on first two layers, interpreted at z=dz/2 and 3*dz/2
    
    # Theory: Insulating surface, therefore dT/dz=0
    # so model for T: T=az^2+bz+c
    #             dT/dz=2az + b = 0 @ z=0 --> b=0
    # reduced model: T=az^2+c
    #           T|dz/2 = T[0,::]  -> T[0,::]=.25*a*dz^2+c
    #           T|3dz/2 = T[1,::] -> T[1,::]=2.25*a*dz^2+c
    #  Subtract equations:
    #  2adz^2 = T[1,::]-T[0,::]
    #  a=(T[1,::]-T[0,::])/(2*dz^2)
    #  c=T[0,::]-.25*a*dz^2
    # At z=0, temperature = c
    #a=(T[1,::]-T[0,::])/(2.0*dz**2.0)
    #c=T[0,::]-.25*a*dz**2.0
    # 
    #return c


    # # alternate model: simple linear extrapolation
    # dT=T[1,::]-T[0,::]
    # return T[0,::]-dT/2.0

    # How about average of above two approaches
    a=(T[1,::]-T[0,::])/(2.0*dz**2.0)
    c=T[0,::]-.25*a*dz**2.0
    dT=T[1,::]-T[0,::]
    linearext=T[0,::]-dT/2.0

    return (linearext+c)/2.0
