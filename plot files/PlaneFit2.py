import numpy as np
    
def PlaneFit2(U = None,winsizeOfx = None,winsizeOfy = None,Rad = None): 
    M,N = U.shape
    UY = U
    UX = U
    UNew = U
    # [tempjj,tempii] = meshgrid(Rad+1:N-Rad, Rad+1:M-Rad); tempii = tempii(:); tempjj = tempjj(:);
    
    #hbar = parfor_progressbar(length(tempii),'Please wait for PlaneFit2!');
    h = waitbar(0,'wait for PlaneFit2!')
    for ii in np.arange((Rad + 1),M - Rad+1).reshape(-1):
        for jj in np.arange((Rad + 1),N - Rad+1).reshape(-1):
            LSMatrix = np.ones(((2 * Rad + 1) ** 2,3))
            for tempj in np.arange(- Rad,Rad+1).reshape(-1):
                for tempi in np.arange(- Rad,Rad+1).reshape(-1):
                    LSMatrix[[2 * Rad + 1] * [tempj + Rad] + tempi + Rad + 1,:] = np.array([1,tempi * winsizeOfx,tempj * winsizeOfy])
            LSb = np.zeros(((2 * Rad + 1) ** 2,1))
            for tempj in np.arange(- Rad,Rad+1).reshape(-1):
                for tempi in np.arange(- Rad,Rad+1).reshape(-1):
                    LSb[[2 * Rad + 1] * [tempj + Rad] + tempi + Rad + 1,:] = U(ii + tempi,jj + tempj)
            tempVector = np.linalg.solve((np.transpose(LSMatrix) * LSMatrix),(np.transpose(LSMatrix) * LSb))
            UNew[ii,jj] = tempVector(1)
            UX[ii,jj] = tempVector(2)
            UY[ii,jj] = tempVector(3)
            # end
# hbar.iterate(1);
            tempij = (ii - Rad - 1) * (N - 2 * Rad) + jj - Rad
            waitbar(tempij / ((M - 2 * Rad) * (N - 2 * Rad)))
    
    return UY,UX,UNew