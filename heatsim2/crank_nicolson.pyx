import sys
import re
import copy
import numpy as np
import numbers
import heatsim2
import heatsim2.expression as expression

# import heatsim2.alternatingdirection as alternatingdirection
import heatsim2.alternatingdirection_c_pyx as alternatingdirection

cimport numpy as np
cimport cython

from libc.stdint cimport uint64_t
from libc.stdint cimport uint8_t

def pm_timeaverage(expr,p_time_prefix=""):
    # Append time to expression variables, and average + and - times
    # OLD CODE...NOT USED!
    ta_expr=expression.linear_expression()
    for cmd in expr.le_cmdlist:
        if cmd[0]==expr.OP_VAR:
            ta_expr.le_cmdlist.append((cmd[0],(cmd[1][0]+'m',cmd[1][1]/2.0)))
            ta_expr.le_cmdlist.append((cmd[0],(p_time_prefix+cmd[1][0]+'p',cmd[1][1]/2.0)))
            ta_expr.le_cmdlist.append((ta_expr.OP_ADD,None))
            pass
        else:
            ta_expr.le_cmdlist.append(copy.deepcopy(cmd))
            pass
        pass
    ta_expr.finalize()
    return ta_expr


def m_time(expr):
    # OLD CODE... NOT USED!
    # Append time to expression variables, '-' times only
    ta_expr=expression.linear_expression()
    for cmd in expr.le_cmdlist:
        if cmd[0]==expr.OP_VAR:
            ta_expr.le_cmdlist.append((cmd[0],(cmd[1][0]+'m',cmd[1][1])))
            pass
        else:
            ta_expr.le_cmdlist.append(copy.deepcopy(cmd))
            pass
        pass
    ta_expr.finalize()
    return ta_expr


def subst_thermal_conductivity(heatflow_expression,k_params):

    (matl_k,matl_k556,matl_k565,matl_k655,matl_k554,matl_k545,matl_k455)=k_params
    heatflow_expression=expression.subst(heatflow_expression,'kmat555',matl_k)
    
    heatflow_expression=expression.subst(heatflow_expression,'kmat556',matl_k556)
    
    heatflow_expression=expression.subst(heatflow_expression,'kmat565',matl_k565)
    
    heatflow_expression=expression.subst(heatflow_expression,'kmat655',matl_k655)
    

    heatflow_expression=expression.subst(heatflow_expression,'kmat554',matl_k554)

    heatflow_expression=expression.subst(heatflow_expression,'kmat545',matl_k545)

    heatflow_expression=expression.subst(heatflow_expression,'kmat455',matl_k455)

    return heatflow_expression


Tshift_re=re.compile(r"""T([4-6mp])([4-6mp])([4-6mp])""")
kshift_re=re.compile(r"""kmat([4-6mp])([4-6mp])([4-6mp])""")

Tshift_map={
    '4': -1.0,
    'm': -.5,
    '5': 0.0,
    'p': 0.5,
    '6': 1.0
}
Tshift_revmap={
    -1.0: '4',
    -0.5: 'm',
    0.0: '5',
    0.5: 'p',
    1.0: '6',
}

def shift_expression(expr,shifts):

    if isinstance(expr,numbers.Number):
        # no variable dependence
        return expr
    
    
    new=expression.linear_expression()
    new.le_cmdlist=list(expr.le_cmdlist)
    for cnt in range(len(new.le_cmdlist)):
        entry=new.le_cmdlist[cnt]
        if entry[0]==expr.OP_VAR:
            Tmatchobj=Tshift_re.match(entry[1][0])
            if Tmatchobj is not None:
                assert(entry[1][2] is None and entry[1][3] is None) # T does not support indexing
                oldshifts=list(map(lambda s: Tshift_map[s],Tmatchobj.groups()))
                newshifts=np.add(oldshifts,shifts)  

                # remove old entry, insert replacement wihh shifted variable name
                del new.le_cmdlist[cnt]
                new.le_cmdlist.insert(cnt,(expr.OP_VAR,('T'+''.join(map(lambda s: Tshift_revmap[s], newshifts)),entry[1][1],None,None)))
                pass
            
            kmatchobj=kshift_re.match(entry[1][0])
            if kmatchobj is not None:
                oldshifts=list(map(lambda s: Tshift_map[s],kmatchobj.groups()))
                newshifts=np.add(oldshifts,shifts) 
                # remove old entry, insert replacement wihh shifted variable name
                del new.le_cmdlist[cnt]
                new.le_cmdlist.insert(cnt,(expr.OP_VAR,('kmat'+''.join(map(lambda s: Tshift_revmap[s], newshifts)),entry[1][1],entry[1][2],entry[1][3])))
                pass
            pass
        pass
    new.finalize()
    return new
            

def setup(z0,y0,x0,
          dz,dy,dx,
          uint64_t nz,uint64_t ny,uint64_t nx,
          dt,
          materials,
          boundaries,
          volumetric,
          np.ndarray[uint8_t,ndim=3] material_elements,
          np.ndarray[uint8_t,ndim=3] boundary_z_elements,
          np.ndarray[uint8_t,ndim=3] boundary_y_elements,
          np.ndarray[uint8_t,ndim=3] boundary_x_elements,
          np.ndarray[uint8_t,ndim=3] volumetric_elements,
          top_surface_y_curvatures=None, # curvature laid out on first layer surface. First layer gets uniform element sizes, deeper layers are stretched per curvature
          top_surface_x_curvatures=None, # We assume principal axes of curvature line up with coordinate axes
          unaligned_anisotropic=False):
    # set unaligned_anisotropic to True if you are solving the full anisotropic problem with an non-diagonal K tensor. Note that this mode hasn't been fully tested


    cdef uint64_t i,j,k
    #cdef double matl_k
    #cdef double matl_k556,matl_k554,matl_k565,matl_k545,matl_k655,matl_k455
    #cdef np.ndarray[double] matl_k
    #cdef np.ndarray[double] matl_k556,matl_k554,matl_k565,matl_k545,matl_k655,matl_k455
    
    # assert(material_elements.dtype==np.uint8)  (tested automatically)
    
    # expressions for placeholder variables
    T555p=expression.linear_expression('T555p')
    T555m=expression.linear_expression('T555m')
    kmatm55=expression.linear_expression('kmatm55')
    kmatp55=expression.linear_expression('kmatp55')
    Tm55=expression.linear_expression('Tm55')
    Tp55=expression.linear_expression('Tp55')
    Tm45=expression.linear_expression('Tm45')
    Tp45=expression.linear_expression('Tp45')
    Tm65=expression.linear_expression('Tm65')
    Tp65=expression.linear_expression('Tp65')
    Tm54=expression.linear_expression('Tm54')
    Tp54=expression.linear_expression('Tp54')
    Tm56=expression.linear_expression('Tm56')
    Tp56=expression.linear_expression('Tp56')
    Tm46=expression.linear_expression('Tm46')
    Tp46=expression.linear_expression('Tp46')
    Tm64=expression.linear_expression('Tm64')
    Tp64=expression.linear_expression('Tp64')

    kmat5m5=expression.linear_expression('kmat5m5')
    kmat5p5=expression.linear_expression('kmat5p5')
    T5m5=expression.linear_expression('T5m5')
    T5p5=expression.linear_expression('T5p5')
    T4m5=expression.linear_expression('T4m5')
    T4p5=expression.linear_expression('T4p5')
    T6m5=expression.linear_expression('T6m5')
    T6p5=expression.linear_expression('T6p5')
    T5m4=expression.linear_expression('T5m4')
    T5p4=expression.linear_expression('T5p4')
    T5m6=expression.linear_expression('T5m6')
    T5p6=expression.linear_expression('T5p6')
    T4m6=expression.linear_expression('T4m6')
    T4p6=expression.linear_expression('T4p6')
    T6m4=expression.linear_expression('T6m4')
    T6p4=expression.linear_expression('T6p4')

    kmat55m=expression.linear_expression('kmat55m')
    kmat55p=expression.linear_expression('kmat55p')
    T55m=expression.linear_expression('T55m')
    T55p=expression.linear_expression('T55p')
    T45m=expression.linear_expression('T45m')
    T45p=expression.linear_expression('T45p')
    T65m=expression.linear_expression('T65m')
    T65p=expression.linear_expression('T65p')
    T54m=expression.linear_expression('T54m')
    T54p=expression.linear_expression('T54p')
    T56m=expression.linear_expression('T56m')
    T56p=expression.linear_expression('T56p')
    T46m=expression.linear_expression('T46m')
    T46p=expression.linear_expression('T46p')
    T64m=expression.linear_expression('T64m')
    T64p=expression.linear_expression('T64p')
    volumetric_source=expression.linear_expression('volumetric_source')

    if top_surface_x_curvatures is None and top_surface_y_curvatures is None:
        # no curvature -- substitute fixed values of dz, dy, dx
        dy_expr=dy
        dx_expr=dx
        pass
    else:
        dy_expr=expression.linear_expression('dy')
        dx_expr=expression.linear_expression('dx')
        pass
    
    # Convert boundaries to evaluated expressions
    evalboundaries=[ (boundary[0].qz(kmatm55,kmatp55,
                                     dz,dy_expr,dx_expr,
                                     Tm55,Tp55,
                                     Tm45,Tp45,
                                     Tm65,Tp65,
                                     Tm54,Tp54,
                                     Tm56,Tp56,
                                     Tm46,Tp46,
                                     Tm64,Tp64,
                                     *boundary[1:]),
                      boundary[0].qy(kmat5m5,kmat5p5,
                                     dz,dy_expr,dx_expr,
                                     T5m5,T5p5,
                                     T4m5,T4p5,
                                     T6m5,T6p5,
                                     T5m4,T5p4,
                                     T5m6,T5p6,
                                     T4m6,T4p6,
                                     T6m4,T6p4,
                                     *boundary[1:]),
                      boundary[0].qx(kmat55m,kmat55p,
                                     dz,dy_expr,dx_expr,
                                     T55m,T55p,
                                     T45m,T45p,
                                     T65m,T65p,
                                     T54m,T54p,
                                     T56m,T56p,
                                     T46m,T46p,
                                     T64m,T64p,
				     *boundary[1:]),) 
                     for boundary in boundaries ]
    
    # global foo
    # foo=evalboundaries
    
    # Set up heat equation
    # rho * c * (T_x^t - T_x^(t-1))/dt = - sigma qnet  on voxel boundaries
    #
    # on boundaries:
    # q = - k grad T

    volume_array=dz*dy*dx

    if top_surface_x_curvatures is not None or top_surface_y_curvatures is not None:
        volume_array=np.zeros((nz,ny,nx),dtype='d')
        pass
        
    
    (ADI_params,ADI_steps)=alternatingdirection.adi_setup((nz,ny,nx),volume_array)   # !!!*** NOTE: We will be modifying volume_array in-place below
    aetam_cache={}
    expression_cache={}
    
    for k in range(nz):
        print('k=%d/%d' % (k,nz))
        for j in range(ny):
            print('j=%d/%d' % (j,ny))
            for i in range(nx):
                # print('i=%d/%d' % (i,nx))
                # Need to include (trivial) equation for TEMPERATURE_FIXED
                if materials[material_elements[k,j,i]][0]==heatsim2.TEMPERATURE_COMPUTE:
                    (compute,matl_k,matl_rho,matl_c)=materials[material_elements[k,j,i]]
                    
                    # q = k grad T = Watts/(meter*Kelvin) * (Kelvin/meter)
                    #   = Watts/m^2
                    # qx * dy *dz = Watts
                    # qx * dy * dz / Volume = Watts/m^3
                    # qx / dx = Watts/m^3

                    be_m55=boundary_z_elements[k,j,i]
                    be_p55=boundary_z_elements[k+1,j,i]
                    be_5m5=boundary_y_elements[k,j,i]
                    be_5p5=boundary_y_elements[k,j+1,i]
                    be_55m=boundary_x_elements[k,j,i]
                    be_55p=boundary_x_elements[k,j,i+1]

                    
                    boundary_params=(be_m55,be_p55,be_5m5,be_5p5,be_55m,be_55p)

                    if i < nx-1 and materials[material_elements[k,j,i+1]][0]==heatsim2.TEMPERATURE_COMPUTE:
                        matl_k556=materials[material_elements[k,j,i+1]][1]
                        pass
                    else:
                        matl_k556=matl_k
                        pass                    

                    if j < ny-1 and materials[material_elements[k,j+1,i]][0]==heatsim2.TEMPERATURE_COMPUTE:
                        matl_k565=materials[material_elements[k,j+1,i]][1]
                        pass
                    else:
                        matl_k565=matl_k
                        pass                    

                    if k < nz-1 and materials[material_elements[k+1,j,i]][0]==heatsim2.TEMPERATURE_COMPUTE:
                        matl_k655=materials[material_elements[k+1,j,i]][1]
                        pass
                    else:
                        matl_k655=matl_k
                        pass                    

                    if i > 0 and materials[material_elements[k,j,i-1]][0]==heatsim2.TEMPERATURE_COMPUTE:
                        matl_k554=materials[material_elements[k,j,i-1]][1]
                        pass
                    else:
                        matl_k554=matl_k
                        pass                    

                    if j > 0 and materials[material_elements[k,j-1,i]][0]==heatsim2.TEMPERATURE_COMPUTE:
                        matl_k545=materials[material_elements[k,j-1,i]][1]
                        pass
                    else:
                        matl_k545=matl_k
                        pass                    

                    if k > 0 and materials[material_elements[k-1,j,i]][0]==heatsim2.TEMPERATURE_COMPUTE:
                        matl_k455=materials[material_elements[k-1,j,i]][1]
                        pass
                    else:
                        matl_k455=matl_k
                        pass                    
                    
                    k_params=(matl_k,matl_k556,matl_k565,matl_k655,matl_k554,matl_k545,matl_k455)
                    k_params_hashable=tuple([ tuple(np.ravel(k_param)) if isinstance(k_param,np.ndarray) else k_param for k_param in k_params ])


                    

                    if top_surface_x_curvatures is None and top_surface_y_curvatures is None:
                        key_params=(heatsim2.TEMPERATURE_COMPUTE,boundary_params,k_params_hashable,matl_rho,matl_c)
                        try: 
                            (heatflow_expressions,time_expressions)=expression_cache[key_params]
                            pass
                        except KeyError:
                        
                            # shift_expression adjusts the indicies in the variables by the specified amount.
                            # So if we have the expression evaluated at the boundary just to our left, and he resulting expression will be evaluated at our location,
                            # we have to shift by -.5
                        
                            # Heatflow_expression represents the net heatflow into an element
                    
                            heatflow_expression=(
                                (shift_expression(evalboundaries[be_m55][0],(-.5,0,0))-shift_expression(evalboundaries[be_p55][0],(+.5,0,0)))*(1.0/dz) +
                                (shift_expression(evalboundaries[be_5m5][1],(0,-.5,0))-shift_expression(evalboundaries[be_5p5][1],(0,+.5,0)))*(1.0/dy) +
                                (shift_expression(evalboundaries[be_55m][2],(0,0,-.5)) - shift_expression(evalboundaries[be_55p][2],(0,0,+.5)))*(1.0/dx) +
                                volumetric_source)
                        
                            # print(heatflow_expression)

                            heatflow_expression=subst_thermal_conductivity(heatflow_expression,k_params) # ,j,i,nz,ny,nx,materials,material_elements)

                            #if (i==0 and j==1):
                            #    print(boundary_x_elements[k,j,i])
                            #    print(heatflow_expression)
                        
                            #    print(shift_expression(evalboundaries[boundary_x_elements[k,j,i+1]][2],(0,0,+.5)))
                            #    print(evalboundaries[boundary_x_elements[k,j,i+1]][2])
                            #    pass
                        
                            # rho c dT/dt = (J/m^3/K) * K/s = J/m^3/s = W/m^3
                            time_expression = -(T555p-T555m)*matl_rho*matl_c*(1.0/dt)                        
                            
                            # print("hfe+te=%s" % (expression.eliminate_groups(heatflow_expression+time_expression).exprstr()))

                            (heatflow_expressions,time_expressions)=alternatingdirection.adi_expressions(heatflow_expression,time_expression,unaligned_anisotropic=unaligned_anisotropic)

                            expression_cache[key_params]=(heatflow_expressions,time_expressions)
                            pass
                        
                        pass
                    else:
                        key_params=(heatsim2.TEMPERATURE_COMPUTE,boundary_params,k_params_hashable,matl_rho,matl_c,top_surface_y_curvatures[j,i],top_surface_x_curvatures[j,i],k)

                        if key_params in expression_cache: 
                            (heatflow_expressions,time_expressions,volume)=expression_cache[key_params]
                            volume_array[k,j,i]=volume

                            pass
                        else:
                            # Have surface curvatures kurv
                            
                            # Positive curvature = concave.
                            
                        
                            # Top surface is nominal size,
                            # dx sweeps out angle theta = kurv*dx
                            # at top of layer k, radius is now (1/kurv + k*dz)
                            # so dx here is (1/kurv + k*dz)*theta
                            # or (1 + k*kurv*dz)*dx
                            
                            dy_top = (1.0 + k*top_surface_y_curvatures[j,i]*dz)*dy
                            dy_bot = (1.0 + (k+1)*top_surface_y_curvatures[j,i]*dz)*dy
                            dy_mean=np.mean((dy_top,dy_bot))
                            
                            dx_top = (1.0 + k*top_surface_x_curvatures[j,i]*dz)*dx
                            dx_bot = (1.0 + (k+1)*top_surface_x_curvatures[j,i]*dz)*dx
                            dx_mean=np.mean((dx_top,dx_bot))
                            
                            # approximate (?) volume from mean edge lengths
                            volume = np.abs(dx_mean*dy_mean*dz)
                            
                            # In this (curvature) case we
                            # multiply through by area
                            # instead of dividing by length
                            # in calculating these expressions,
                            # as the areas can vary due to curvature
                            heatflow_expression=(
                                (shift_expression(evalboundaries[be_m55][0],(-.5,0,0))*dx_top*dy_top-shift_expression(evalboundaries[be_p55][0],(+.5,0,0))*dx_bot*dy_bot) +
                                (shift_expression(evalboundaries[be_5m5][1],(0,-.5,0))-shift_expression(evalboundaries[be_5p5][1],(0,+.5,0)))*(dx_mean*dz) +
                                (shift_expression(evalboundaries[be_55m][2],(0,0,-.5)) - shift_expression(evalboundaries[be_55p][2],(0,0,+.5)))*(dy_mean*dz) +
                                volumetric_source*volume)
                            
                            heatflow_expression=expression.subst(heatflow_expression,'dx',dx_mean)   # arguably should keep left and right dx's separate here...
                            heatflow_expression=expression.subst(heatflow_expression,'dy',dy_mean)
                            # print(heatflow_expression)
                        
                            heatflow_expression=subst_thermal_conductivity(heatflow_expression,k_params) # ,j,i,nz,ny,nx,materials,material_elements)
                            
                            #if (i==0 and j==1):
                            #    print(boundary_x_elements[k,j,i])
                            #    print(heatflow_expression)
                            
                            #    print(shift_expression(evalboundaries[boundary_x_elements[k,j,i+1]][2],(0,0,+.5)))
                            #    print(evalboundaries[boundary_x_elements[k,j,i+1]][2])
                            #    pass
                            
                            # rho c dT/dt = (J/m^3/K) * K/s = J/m^3/s = W/m^3

                            # here we are multiplied through by volume compared to the uncurved case
                            time_expression = -(T555p-T555m)*matl_rho*matl_c*volume*(1.0/dt)
                        
                            # print("hfe+te=%s" % (expression.eliminate_groups(heatflow_expression+time_expression).exprstr()))
                            
                            (heatflow_expressions,time_expressions)=alternatingdirection.adi_expressions(heatflow_expression,time_expression,unaligned_anisotropic=unaligned_anisotropic)
                            
                            expression_cache[key_params]=(heatflow_expressions,time_expressions,volume)
                            volume_array[k,j,i]=volume

                            pass
                        
                        pass
                    #if (k==1 and i==5 and j==5):
                    #    for cnt in range(len(heatflow_expressions)):
                            
                            # print("hfe[%d]=%s" % (cnt,str(expression.eliminate_groups(heatflow_expressions[cnt]+time_expressions[cnt]))))
                            #print("hfe[%d]=%s" % (cnt,str(expression.eliminate_groups(heatflow_expressions[cnt]+time_expressions[cnt]).exprstr())))
                            #print("hfe[%d]=%s" % (cnt,str(expression.eliminate_groups(heatflow_expressions[cnt]+time_expressions[cnt]).fullreduce())))
                            # print("hfe[%d]=%s" % (cnt,expression.eliminate_groups(heatflow_expressions[cnt]+time_expressions[cnt]).fullreduce().exprstr()))
                            # print("hfe[%d]=%s" % (cnt,str(expression.eliminate_groups(heatflow_expressions[cnt]+time_expressions[cnt]).dictform())))

                            # print("te[%d]=%s" % (cnt,str(time_expressions[cnt].fullreduce())))
                            #pass
                                  
                    alternatingdirection.add_equation_to_adi_matrices(ADI_params,ADI_steps,k,j,i,key_params,aetam_cache,heatflow_expressions,time_expressions)
                    pass
                else:
                    # TEMPERATURE_FIXED
                    assert(materials[material_elements[k,j,i]][0]==heatsim2.TEMPERATURE_FIXED)

                    key_params=(heatsim2.TEMPERATURE_FIXED)

                    heatflow_expression=expression.linear_expression(0.0)  # This side of the equation doesn't matter -- is zero

                    time_expression = -(T555p-T555m)   #*matl_rho*matl_c*(1.0/dt)                        

                    # We should probably cache these results as above, despite their triviality
                    (heatflow_expressions,time_expressions)=alternatingdirection.adi_expressions(heatflow_expression,time_expression,unaligned_anisotropic=unaligned_anisotropic)

                    alternatingdirection.add_equation_to_adi_matrices(ADI_params,ADI_steps,k,j,i,key_params,aetam_cache,heatflow_expressions,time_expressions)

                    
                pass
            pass
        pass

    
    for cnt in range(len(ADI_steps)):
        ADI_steps[cnt].finalize()
        pass
    
    return (ADI_params,ADI_steps)
    
