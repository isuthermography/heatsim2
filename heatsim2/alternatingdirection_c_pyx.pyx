import re
import copy
import sys
import numpy as np
import scipy
import scipy.sparse

import heatsim2
import heatsim2.tridiag as tridiag
import heatsim2.expression as expression
from heatsim2.expression import crank_subst_in_groups
from heatsim2.expression import eliminate_groups

from libc.stdio cimport snprintf
from libc.stdio cimport stderr
from libc.stdio cimport fprintf
from libc.string cimport strcpy
from libc.string cimport strlen
from libc.string cimport strcmp
from libc.stdlib cimport strtol

from libc.stdint cimport uint64_t
from libc.stdint cimport int64_t
from libc.stdlib cimport malloc, free
from libc.string cimport strcmp
from libc.stddef cimport size_t
#from cpython.string cimport PyString_AsString
#from cpython.object cimport PyObject

cimport numpy as np 
cimport cython

#from cpython.cobject cimport PyCObject_FromVoidPtr, PyCObject_AsVoidPtr

np.import_array()

cdef extern from "cython_capsule_supp.h":
    object CObject_Or_Capsule_New(void *Ptr)
    void *CObject_Or_Capsule_AsVoidPtr(object Obj)
    char **to_cstring_array(object list_str)
    pass

cdef extern from "alternatingdirection_c.h": 
     ctypedef struct adi_step:
         uint64_t shape[3]
         uint64_t n
         int permuteorder[3]
         int invpermuteorder[3]
         int stepnum

         uint64_t permutedshape[3]
         
         double *Amat
         uint64_t *BmatIndices
         double *BmatData
         uint64_t BmatEntriesUsed

         int NumCmats

         uint64_t **CmatsIndices
         double **CmatsData
         uint64_t *CmatsEntriesUsed

         double *Dvec
         pass
     adi_step *create_adi_step(uint64_t shape[3],uint64_t permutedshape[3],int permuteorder[3],int invpermuteorder[3],int stepnum, int NumCmats)

     void delete_adi_step(adi_step *step)
     void add_equation(adi_step *step,uint64_t posindex[3],char **eqvarnames,double *eqvalues,int numeqvars)
     
     pass



class adi_params(object):
    shape=None
    volume_array=None

    def __init__(self,**kwargs):
        for key in kwargs:
            if not hasattr(self,key):
                raise ValueError("Unknown parameter %s (must add to class definition)")
            setattr(self,key,kwargs[key])
        pass
    pass

    
# Process...
# 1. Reshape into a matrix problem for 3.1a:
#    A1*P1 = B1*M1 + D1
# 2. Contribute to matrices A1 (tridiagonal) and B1 (sparse)
# 3. Reshape into a matrix problem for 3.1b
#    A2*P2 = B2*M2+C2*P1 + D2
# 4. Contribute to matrices A2 (tridiagonal), B2 (sparse), and C2 (sparse) 

    # incoming temperature matrices are indexed [k,j,i]
    # reshape(nz*ny*nx)
    # permuteorder=(0,1,2)
    # now do 3.1a

    # now reshape back: reshape(nz,ny,nx)
    # invpermuteorder=(0,1,2)
    # Now do 3.1b:
    # transpose(,(2,0,1)).reshape(nx*ny*nz)
    #permuteorder=(2,0,1)
    # now reshape back
    # .reshape(nx,nz,ny)
    # invpermuteorder=(1,2,0)

# Variable names
# --------------
# T555p  means center in z,y,x +1/2 in time
# T554m  means center in z,y, -1 in X, -1/2 in time
# T554  means center in z,y, -1 in X, 
    
varnamere=re.compile(r"""T([4-6])([4-6])([4-6])([pm]?)(\d*)""")

# pyadi_step is not an extension class therefore use __del__
# to ensure destruction
class pyadi_step(object):
    ADI_params=None
    
    #cdef adi_step *c_step
    permuteorder=None
    invpermuteorder=None
    c_step=None
    
    Amat=None
    Bmat=None
    Cmats=None
    Dvec=None
    
    Lmat=None # tridiagonal L for LU factorization .... created by finalize()
    Umat=None # tridiagonal U for LU factorization .... created by finalize()
    
    
    def __init__(self,ADI_params,permuteorder,stepnum):
        cdef uint64_t cshape[3]
        cdef int cpermuteorder[3]
        cdef uint64_t cpermutedshape[3]
        cdef int cinvpermuteorder[3]
        cdef adi_step *c_step

        self.ADI_params=ADI_params

        (cshape[0],cshape[1],cshape[2])=ADI_params.shape
        
        self.permuteorder=tuple(permuteorder)
        (cpermuteorder[0],cpermuteorder[1],cpermuteorder[2])=self.permuteorder

        self.permutedshape = [ ADI_params.shape[index] for index in permuteorder ]
        (cpermutedshape[0],cpermutedshape[1],cpermutedshape[2])=self.permutedshape
        
        invpermute=[ (pord,ind) for ind,pord in enumerate(permuteorder) ]
        invpermute.sort(key = lambda x: x[0])
        self.invpermuteorder=tuple([ invpord  for ind,invpord in invpermute ])
        (cinvpermuteorder[0],cinvpermuteorder[1],cinvpermuteorder[2])=self.invpermuteorder
        n=np.prod(ADI_params.shape)
        

        c_step=create_adi_step(cshape,cpermutedshape,cpermuteorder,cinvpermuteorder,stepnum, stepnum)

	
        self.c_step=CObject_Or_Capsule_New(c_step)
        # !!! *** Does the PyCObject need to be reference-counted manually???

        pass

    def __del__(self):
        cdef adi_step *c_step 
	#c_step=<adi_step *>PyCObject_AsVoidPtr(self.c_step)
        c_step=<adi_step *>CObject_Or_Capsule_AsVoidPtr(self.c_step);
        delete_adi_step(c_step);
        self.c_step=None
        pass

    def add_equation(self,posindex,eqdict):
        cdef adi_step *c_step
        cdef uint64_t cposindex[3]
        cdef char **eqvarnames
        cdef np.ndarray[double,ndim=1,mode="c"] eqvals
        
	#c_step=<adi_step *>PyCObject_AsVoidPtr(self.c_step)
        c_step=<adi_step *>CObject_Or_Capsule_AsVoidPtr(self.c_step);
        (cposindex[0],cposindex[1],cposindex[2])=posindex

        eqkeys=eqdict.keys()

        #if sys.version_info[0] >= 3:
        #    # python3 only

        # We are passing the key list to C code... pass as utf-8 strings/bytes
        eqkeylist = [ eqkey.encode('utf-8') for eqkey in eqkeys ]
        
        #    pass
        #else:
        #    eqkeylist=eqkeys
        #    pass
        
        
        #sys.stderr.write("eqdict=%s\n" % (str(eqdict)))
        eqvals=np.array([ eqdict[key] for key in eqkeys ],dtype='d')

        eqvarnames=to_cstring_array(eqkeylist)
        
        
        add_equation(c_step,cposindex,eqvarnames,<double *>eqvals.data,len(eqkeys))

        free(eqvarnames)
        

    def finalize(self):
        # Calculate L and U
        cdef np.npy_intp Amatdims[2]
        cdef np.npy_intp BmatIdxdims[2]
        cdef np.npy_intp BmatDatadim
        cdef np.npy_intp CmatIdxdims[2]
        cdef np.npy_intp CmatDatadim
        cdef np.npy_intp Dvecdim
        cdef adi_step *c_step 

        #c_step=<adi_step *>PyCObject_AsVoidPtr(self.c_step)
        c_step=<adi_step *>CObject_Or_Capsule_AsVoidPtr(self.c_step);
        
        n=np.prod(self.ADI_params.shape)

        Amatdims[0]=n;
        Amatdims[1]=3;
        #fprintf(stderr,"%lx\n",<long>c_step.Amat);
        self.Amat=np.PyArray_SimpleNewFromData(2, &Amatdims[0], np.NPY_DOUBLE, c_step.Amat)

        BmatIdxdims[0]=c_step.BmatEntriesUsed;
        BmatIdxdims[1]=2;
        BmatIndices=np.PyArray_SimpleNewFromData(2,BmatIdxdims,np.NPY_UINT64,c_step.BmatIndices)
        
        
        BmatDatadim=c_step.BmatEntriesUsed;
        BmatData=np.PyArray_SimpleNewFromData(1,&BmatDatadim,np.NPY_DOUBLE,c_step.BmatData)

        # use coo as intermediate format so that duplicate-referenced entries
        # will be correctly added together
        # print BmatIndices[:,0].shape
        # print BmatIndices[:,1].shape
        # print max(BmatIndices[:,0])
        # print max(BmatIndices[:,1])
        # print BmatData.shape
        # print n
        #print BmatIndices.dtype
        #print BmatData.dtype

        self.Bmat=scipy.sparse.coo_matrix((BmatData,(BmatIndices[:,0],BmatIndices[:,1])),shape=(n,n)).tocsr()

        self.Cmats=[]

        CmatIdxdims[1]=2;
        
        


        for cnt in range(c_step.NumCmats):
            CmatIdxdims[0]=c_step.CmatsEntriesUsed[cnt]
            CmatDatadim=c_step.CmatsEntriesUsed[cnt]

            CmatIndices=np.PyArray_SimpleNewFromData(2,CmatIdxdims,np.NPY_UINT64,c_step.CmatsIndices[cnt])
            CmatData=np.PyArray_SimpleNewFromData(1,&CmatDatadim,np.NPY_DOUBLE,c_step.CmatsData[cnt])

            # use coo as intermediate format so that duplicate-referenced entries
            # will be correctly added together
            self.Cmats.append(scipy.sparse.coo_matrix((CmatData,(CmatIndices[:,0],CmatIndices[:,1])),shape=(n,n)).tocsr())
            pass

        Dvecdim=n
        self.Dvec=np.PyArray_SimpleNewFromData(1,&Dvecdim,np.NPY_DOUBLE,c_step.Dvec)

                              
        
        # Perform tridiagonal LU decomposition
        
        (self.Lmat,self.Umat)=tridiag.tridiaglu(self.Amat)


        pass
        
    pass


def run_adi_steps(ADI_params,ADI_steps,t,dt,Tarray,volumetric_elements,volumetric):
    cdef adi_step *c_step 
    cdef np.ndarray[double,ndim=4] stepresults
    cdef np.ndarray[double,ndim=1] rhs

    n=np.prod(ADI_params.shape)

    volumetric_array=np.zeros(ADI_params.shape,dtype='d')
    for volumetric_index in range(len(volumetric)):
        if volumetric[volumetric_index][0]==heatsim2.IMPULSE_SOURCE:
            if t==volumetric[volumetric_index][1]:
                volumetric_array[volumetric_elements==volumetric_index]=volumetric[volumetric_index][2]/dt
                pass
            pass
        elif volumetric[volumetric_index][0]==heatsim2.STEPPED_SOURCE:
            (stepped_source,starttime,endtime,wattsperm3)=volumetric[volumetric_index]
            if t >= starttime and t <= endtime:
                volumetric_array[volumetric_elements==volumetric_index]=wattsperm3
                pass
            pass
        elif volumetric[volumetric_index][0]==heatsim2.IMPULSE_POINT_SOURCE_JOULES:
            if t==volumetric[volumetric_index][1]:
                # Specified number of joules into EACH ELEMENT... need to correct for element size
                if len(np.array(ADI_params.volume_array).shape) > 0:
                    # Volume_array is actually an array
                    volumetric_array[volumetric_elements==volumetric_index]=volumetric[volumetric_index][2]/(ADI_params.volume_array[volumetric_elements==volumetric_index]*dt)
                    pass
                else:
                    # volume_array is a float
                    volumetric_array[volumetric_elements==volumetric_index]=volumetric[volumetric_index][2]/(ADI_params.volume_array*dt)
                    pass
                pass
            pass
        elif volumetric[volumetric_index][0]==heatsim2.SPATIALLY_Z_DECAYING_TEMPORAL_IMPULSE:
            # assumption behind this is decay is along
            # a vector in the + or - Z direction
            # and a reference offset. Remember we write coordinates (z,x,y)
            # Impulse applies where the 
            # position vector multiplied by direction is greater than 
            # the reference offset. 
            # Impulse amplitude is supplied as J/m^2
            # Absorptivity is defined by a characteristic length; 
            # assumption is that at distance in d, 
            # the remaining amplitude = exp(-d/l) 
            (spatially_z_decaying_temporal_impulse,
             impulse_time,
             decaydirec, # should be a unit vector either in +z or -z
             offset,
             z_ndgrid,
             decayimpulse_dz,
             joulesperm2,
             characlength) = volumetric[volumetric_index]
            if t==impulse_time:
                assert((decaydirec==np.array([1.0,0.0,0.0])).all() or (decaydirec==np.array([-1.0,0.0,0.0])).all())
                offset=float(offset)
                decayimpulse_dz=abs(float(decayimpulse_dz))
                joulesperm2=float(joulesperm2)
                characlength=float(characlength)
                assert(characlength > 0)
                # Over an infinite thickness, (1/l)exp(-d/l)dd integrates to 1.0
                # The value we drop into volumetric_array has to have units 
                # of (J/(s*m^3))
                #
                # at each element, integrate from leftzboundary to rightzboundary 
                # of (amplitude in J/m^2)*(1/l)exp(-d/l)dd
                #    ... This gives us the portion of energy, still in J/m^2
                # deposited into that element 
                # Then we divide by dz to get energy per unit volume
                # and divide by dt to get power per unit volume

                decaying_impulse_elements=volumetric_elements==volumetric_index
                
                # print("die count=%d" % (np.count_nonzero(decaying_impulse_elements)))
                centerpos=z_ndgrid*decaydirec[0]-offset  # center position of each element, relative to offset
                leftzboundary=centerpos-decayimpulse_dz/2.0
                rightzboundary=centerpos+decayimpulse_dz/2.0
                use_decaying_impulse_elements=decaying_impulse_elements & (rightzboundary > 0.0) # only get an integral contribution where we are integrating past zero
                
                # print("udie count=%d" % (np.count_nonzero(use_decaying_impulse_elements)))
                
                userightzboundary=rightzboundary[use_decaying_impulse_elements]
                useleftzboundary=leftzboundary[use_decaying_impulse_elements]

                # can't integrate starting from anything less than 0
                useleftzboundary[useleftzboundary < 0.0]=0.0

                
                # integral from d=a to d=b of (1/l)exp(-d/l)dd
                # is -exp(-b/l)+exp(-a/l)
                zintegral= -np.exp(-userightzboundary/characlength) + np.exp(-useleftzboundary/characlength)
                
                element_energy_per_unit_area=zintegral*joulesperm2
                element_power_per_unit_volume=element_energy_per_unit_area/(decayimpulse_dz*dt)
                # print("epuv=%s" % (str(element_power_per_unit_volume)))

                volumetric_array[use_decaying_impulse_elements]=element_power_per_unit_volume
                pass
            pass
        
        pass
    
    stepresults=np.empty((len(ADI_steps),ADI_params.shape[0],ADI_params.shape[1],ADI_params.shape[2]),dtype='d')
    for stepnum in range(len(ADI_steps)):
        step=ADI_steps[stepnum]
        #c_step=<adi_step *>PyCObject_AsVoidPtr(step.c_step)
        c_step=<adi_step *>CObject_Or_Capsule_AsVoidPtr(step.c_step);

        

        # Permute temperature array then reshape as a vector
        Tvector=Tarray.transpose(step.permuteorder).reshape(n)

        # Solve:
        # Amat*Tvectornew = Bmat*Tvectorold + Cmat[0]*step0result + Cmat[1]*step1result+... + Dvec
        # Construct right-hand-side
        rhs=step.Bmat.dot(Tvector)
        for cnt in range(stepnum):
            rhs+=step.Cmats[cnt].dot(stepresults[cnt,::].transpose(step.permuteorder).reshape(n))            
            pass
        
        rhs+=step.Dvec*volumetric_array.transpose(step.permuteorder).reshape(n)


        #print step.Amat

        stepresults[stepnum,::]=tridiag.tridiagsolve(step.Lmat,step.Umat,rhs).reshape(step.permutedshape).transpose(step.invpermuteorder)
        
        
        pass
    return stepresults[-1,::]  # return output from last step

def adi_setup(shape,volume_array):
    # ... would need to be modified to support anisotropy
    ADI_params=adi_params(shape=shape,volume_array=volume_array)
    ADI_steps=[]

    # !!! Order here must match adi_expressions()

    # Douglas_ADI3D.pdf, eq 3.1a: Operate on x (last index) first
    ADI_steps.append(pyadi_step(ADI_params,(0,1,2),0))

    # Douglas_ADI3D.pdf, eq 3.1b: Operate on y (middle index) second
    ADI_steps.append(pyadi_step(ADI_params,(2,0,1),1))

    # Douglas_ADI3D.pdf, eq 3.1c: Operate on z (first index) third
    ADI_steps.append(pyadi_step(ADI_params,(1,2,0),2))

    # # Reverse order trial
    # # Douglas_ADI3D.pdf, eq 3.1c: Operate on z (first index) first
    # ADI_steps.append(pyadi_step(ADI_params,(1,2,0),0))

    # # Douglas_ADI3D.pdf, eq 3.1b: Operate on y (middle index) second
    # ADI_steps.append(pyadi_step(ADI_params,(2,0,1),1))

    # # Douglas_ADI3D.pdf, eq 3.1a: Operate on x (last index) third
    # ADI_steps.append(pyadi_step(ADI_params,(0,1,2),2))

    
    return (ADI_params,ADI_steps)

def add_equation_to_adi_matrices(ADI_params,ADI_steps,k,j,i,key_params,aetam_cache,spatial_expressions,time_expressions):

    # key_params are a dictionary key for the parameters that make
    # for unique results within the context of an aetam_cache dictionary

    assert(len(ADI_steps)==len(spatial_expressions))

    # perform caching
    try:
        eqdicts=aetam_cache[key_params]
    except KeyError:
        eqdicts=[]
        for stepnum in range(len(ADI_steps)):
            eqdict=eliminate_groups(spatial_expressions[stepnum]+time_expressions[stepnum]).dictform()
            eqdicts.append(eqdict)
            pass
        aetam_cache[key_params]=eqdicts
        
        pass
    
    for stepnum in range(len(ADI_steps)):
        #print eliminate_groups(spatial_expressions[stepnum]+time_expressions[stepnum])
        eqdict=eqdicts[stepnum]        
        ADI_steps[stepnum].add_equation((k,j,i),eqdict)
        pass
    pass
    

def adi_expressions(spatial_expression,time_expression,unaligned_anisotropic=False):

    spatial_expressions=[]
    time_expressions=[]

    # !!! Order here must match adi_setup()
    
    # Douglas_ADI3D.pdf, eq 3.1a
    spatial_expressions.append(crank_subst_in_groups(spatial_expression,('T556','T555','T554'),0))
    time_expressions.append(expression.subst(time_expression,'T555p','T555p0'))

    # Douglas_ADI3D.pdf, eq 3.1b
    spatial_expressions.append(crank_subst_in_groups(spatial_expressions[-1],('T565','T555','T545'),1))
    time_expressions.append(expression.subst(time_expression,'T555p','T555p1'))


    # Douglas_ADI3D.pdf, eq 3.1c
    spatial_expressions.append(crank_subst_in_groups(spatial_expressions[-1],('T655','T555','T455'),2))
    time_expressions.append(expression.subst(time_expression,'T555p','T555p2'))

    # # Try reverse order...

    # # Douglas_ADI3D.pdf, eq 3.1c
    # spatial_expressions.append(crank_subst_in_groups(spatial_expression,('T655','T555','T455'),0))
    # time_expressions.append(expression.subst(time_expression,'T555p','T555p0'))

    # # Douglas_ADI3D.pdf, eq 3.1b
    # spatial_expressions.append(crank_subst_in_groups(spatial_expressions[-1],('T565','T555','T545'),1))
    # time_expressions.append(expression.subst(time_expression,'T555p','T555p1'))



    # # Douglas_ADI3D.pdf, eq 3.1a
    # spatial_expressions.append(crank_subst_in_groups(spatial_expressions[-1],('T556','T555','T554'),2))
    # time_expressions.append(expression.subst(time_expression,'T555p','T555p2'))

    if unaligned_anisotropic:
        # Douglas_ADI3D.pdf, extended
        # First digit is a 6, step thru 3rd digit
        spatial_expressions.append(crank_subst_in_groups(spatial_expressions[-1],('T656','T655','T654'),3))
        time_expressions.append(expression.subst(time_expression,'T555p','T555p3'))

        # First digit is a 6, step thru 2nd digit
        spatial_expressions.append(crank_subst_in_groups(spatial_expressions[-1],('T665','T655','T645'),4))
        time_expressions.append(expression.subst(time_expression,'T555p','T555p4'))

        # First digit is a 4, step thru 3rd digit
        spatial_expressions.append(crank_subst_in_groups(spatial_expressions[-1],('T456','T455','T454'),5))
        time_expressions.append(expression.subst(time_expression,'T555p','T555p5'))

        # First digit is a 4, step thru 2nd digit
        spatial_expressions.append(crank_subst_in_groups(spatial_expressions[-1],('T465','T455','T445'),6))
        time_expressions.append(expression.subst(time_expression,'T555p','T555p6'))


        # Second digit is a 6, step thru 3rd digit
        spatial_expressions.append(crank_subst_in_groups(spatial_expressions[-1],('T566','T565','T564'),7))
        time_expressions.append(expression.subst(time_expression,'T555p','T555p7'))

        # Second digit is a 6, step thru 1st digit
        spatial_expressions.append(crank_subst_in_groups(spatial_expressions[-1],('T665','T565','T465'),8))
        time_expressions.append(expression.subst(time_expression,'T555p','T555p8'))

        # Second digit is a 4, step thru 3rd digit
        spatial_expressions.append(crank_subst_in_groups(spatial_expressions[-1],('T546','T545','T544'),9))
        time_expressions.append(expression.subst(time_expression,'T555p','T555p9'))

        # Second digit is a 4, step thru 1st digit
        spatial_expressions.append(crank_subst_in_groups(spatial_expressions[-1],('T645','T545','T445'),10))
        time_expressions.append(expression.subst(time_expression,'T555p','T555p10'))

        # Third digit is a 6, step thru 2nd digit
        spatial_expressions.append(crank_subst_in_groups(spatial_expressions[-1],('T566','T556','T546'),11))
        time_expressions.append(expression.subst(time_expression,'T555p','T555p11'))

        # Third digit is a 6, step thru 1st digit
        spatial_expressions.append(crank_subst_in_groups(spatial_expressions[-1],('T656','T556','T456'),12))
        time_expressions.append(expression.subst(time_expression,'T555p','T555p12'))

        # Third digit is a 4, step thru 2nd digit
        spatial_expressions.append(crank_subst_in_groups(spatial_expressions[-1],('T564','T554','T544'),13))
        time_expressions.append(expression.subst(time_expression,'T555p','T555p13'))

        # Third digit is a 4, step thru 1st digit
        spatial_expressions.append(crank_subst_in_groups(spatial_expressions[-1],('T654','T554','T454'),14))
        time_expressions.append(expression.subst(time_expression,'T555p','T555p14'))



        pass
    
    spatial_expressions[-1]=spatial_expressions[-1].fullreduce()  # reduce expression in case it is anisotropic and aligned with coordinate axes, then we don't have to do the full anisotropic solution.

    #sys.stderr.write("spatial_expressions[-1]=%s\n" % (spatial_expressions[-1].exprstr()))

    assert(expression.no_groups(spatial_expressions[-1]))  # if it errors out here, you probably need to set unaligned_anisotropic to True

    return (tuple(spatial_expressions),tuple(time_expressions))  # convert to tuples so that we can use them as dictionary index
