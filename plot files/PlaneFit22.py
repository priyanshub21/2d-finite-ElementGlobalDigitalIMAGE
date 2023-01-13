import numpy as np
    
def PlaneFit22(U = None,xstep = None,ystep = None,Rad = None): 
    #PLANEFIT22: 2D least square fittings
#	[UY,UX,UNew] = PLANEFIT2(U,xstep,ystep,Rad)
    
    # Input: U        a gridded data with size (M,N)
#        xstep    x spacing between gridded data
#        ystep    y spacing between gridded data
#        Rad      plane fitting radius
    
    # Output: UY      Fitted dudy
#         UX      Fitted dudx
#         UNew    Fitted U
    
    # Formula to fit:  u(x,y) = UNew + UX*(x-x0) + UY*(y-y0)
    
    # Author: Jin Yang, aldicdvc@gmail.com
# Date: 2020.12
############################################
    
    ##
    M,N = U.shape
    
    iGrid,jGrid = ndgrid(np.arange(Rad + 1,M - Rad+1),np.arange(Rad + 1,N - Rad+1))
    
    UNew = 0 * iGrid
    UX = UNew
    UY = UNew
    
    iGrid = iGrid
    jGrid = jGrid
    
    xGrid,yGrid = ndgrid(np.arange(- Rad,Rad+1),np.arange(- Rad,Rad+1))
    ##
# hbar = parfor_progressbar(length(iGrid) ,'Please wait for PlaneFit2!');
    hbar = waitbar(0,'wait for PlaneFit2!')
    for ptInd in np.arange(1,len(iGrid)+1).reshape(-1):
        # hbar.iterate(1);
        waitbar(ptInd / len(iGrid))
        ii = iGrid(ptInd)
        jj = jGrid(ptInd)
        if np.isnan(U(ii,jj)) == 1:
            UNew[ptInd] = nan
            UX[ptInd] = nan
            UY[ptInd] = nan
        else:
            LSMatrix = np.array([np.ones(((2 * Rad + 1) ** 2,1)),xGrid * xstep,yGrid * ystep])
            LSb = U(sub2ind(np.array([M,N]),xGrid + ii,yGrid + jj))
            row,col = find(np.isnan(LSb) == 0)
            tempVector = np.linalg.solve(LSMatrix(row,:),LSb(row,:))
            UNew[ptInd] = tempVector(1)
            UX[ptInd] = tempVector(2)
            UY[ptInd] = tempVector(3)
    
    close_(hbar)
    return UNew,UX,UY,iGrid,jGrid
    
    return UNew,UX,UY,iGrid,jGrid