import numpy as np
import copy

def tridiaglu(Amat):
    n=Amat.shape[0]
    Lmat=np.zeros((n,3),dtype='d')
    Umat=copy.copy(Amat)
    # Lmat[:,1]=1.0  # Can be removed after testing
    
    assert(Umat[0,0]==0.0)
    assert(Umat[-1,2]==0.0)
    
    for row in range(n):
        
        # normalize row to obtain '1' on diagonal
        Lmat[row,1]=Umat[row,1]
        #self.Umat[row,:]/=self.Umat[row,1]
        Umat[row,2]/=Umat[row,1]
        Umat[row,1]=1.0
        if row < n-1:
            # Subtract out from following row
            # to zero out leading columns
            coefficient=Umat[row+1,0]
            # self.Umat[row+1,0:2]-=self.Umat[row,1:3]*coefficient
            Umat[row+1,0]=0.0
            Umat[row+1,1]-=Umat[row,2]*coefficient
            Lmat[row+1,0]=coefficient
            pass
        pass

    return (Lmat,Umat)


def tridiagsolve(Lmat,Umat,bvec):
    # solve LUx=b  -> Ux=Linv*b
    # Ux=c   where  c=Linv*b, i.e. Lc=b
    # so first solve Lc=b

    xvec=np.zeros(len(bvec),dtype='d')
    
    xvec[0]=bvec[0]/Lmat[0,1]
    for row in range(1,len(bvec)):
        xvec[row]=(bvec[row]-Lmat[row,0]*xvec[row-1])/Lmat[row,1]
        pass

    # Now solve Ux=c where c is the xvec from above
    # xvec[len(bvec)-1]=xvec[len(bvec)-1]
    for row in range(len(bvec)-2,-1,-1):
        xvec[row]=xvec[row]-Umat[row,2]*xvec[row+1]
        pass

    return xvec

if __name__=="__main__":
    # perform dianostics
    import numpy.random
    import scipy.sparse

    n=4
    Amat=np.random.rand(n,3)
    Amat[0,0]=0.0
    Amat[-1,2]=0.0

    (Lmat,Umat)=tridiaglu(Amat)


    # Our diagonal form is a bit difference from what dia_matrix
    # is expecting, so we have to reorganize a bit....
    denseAmat=scipy.sparse.dia_matrix((np.array([np.concatenate((Amat[1:,0],(0,))),Amat[:,1],np.concatenate(((0,),Amat[:-1,2]))],dtype='d'),[ -1, 0, 1 ]), shape=(n,n)).todense()

    denseLmat=scipy.sparse.dia_matrix((np.array([np.concatenate((Lmat[1:,0],(0,))),Lmat[:,1],np.concatenate(((0,),Lmat[:-1,2]))],dtype='d'),[ -1, 0, 1 ]), shape=(n,n)).todense()

    denseUmat=scipy.sparse.dia_matrix((np.array([np.concatenate((Umat[1:,0],(0,))),Umat[:,1],np.concatenate(((0,),Umat[:-1,2]))],dtype='d'),[ -1, 0, 1 ]), shape=(n,n)).todense()


    result=np.dot(denseLmat,denseUmat)
    if np.linalg.norm(denseAmat-result)/np.linalg.norm(denseAmat) < 1e-10:
        print('Passes LU factorization diagnostics...\n')
        pass
    else:
        print('Fails LU factorization diagnostics...\n')
        pass

    
    b=np.random.rand(n);
    soln=np.squeeze(np.asarray(np.dot(np.linalg.inv(denseAmat),b)))

    tdsoln=tridiagsolve(Lmat,Umat,b)

    if np.linalg.norm(soln-tdsoln)/np.linalg.norm(soln) < 1e-10:
        print('Passes LU solver diagnostics...\n')
        pass
    else:
        print('Fails LU solver diagnostics...\n')
        pass

        
