import numpy as np
    
def Init(u = None,v = None,Phi = None,x = None,y = None,index = None): 
    #FUNCTION U0 =Init(u,v,Phi,x,y,index)
# ----------------------------------------------
# Remove outliers and smooth initial displacements
    
    # ----------------------------------------------
# Author: Jin Yang.
# Contact and support: jyang526@wisc.edu -or- aldicdvc@gmail.com
# Last time updated: 02/2020.
# ==============================================
    
    ## inpaint_nans: I use a Spring Model for u and v initial guesses
    u = inpaint_nans(u,4)
    v = inpaint_nans(v,4)
    uInit = u
    vInit = v
    # threshod = 0.5;
# [row, col] = find(Phi>2*(1-threshod));
# uInit(row,col) = NaN; vInit(row,col) = NaN;
    
    # uInit = inpaint_nans(uInit,4);
# vInit = inpaint_nans(vInit,4);
    
    # ##########PLOTTING##################
# figure;surf(uInit)
# colorbar
# figure;surf(vInit)
# colorbar
# ##########PLOTTING##################
    
    if (index == 1):
        utemp = uInit
        vtemp = vInit
        Mtemp = u.shape[2-1]
        Ntemp = u.shape[1-1]
        for i in np.arange(2,Ntemp - 1+1).reshape(-1):
            for j in np.arange(2,Mtemp - 1+1).reshape(-1):
                NeighborOfuInit = np.array([[utemp(i - 1,j)],[utemp(i + 1,j)],[utemp(i,j - 1)],[utemp(i,j + 1)],[utemp(i - 1,j - 1)],[utemp(i + 1,j - 1)],[utemp(i + 1,j + 1)],[utemp(i - 1,j + 1)]])
                avgOfuInit[i,j] = (1 / 8 * sum(NeighborOfuInit(np.arange(1,4+1))) + 1 / 16 * sum(NeighborOfuInit(np.arange(5,8+1)))) / (3 / 4)
                stdOfuInit[i,j] = std(NeighborOfuInit,'omitnan')
                NeighborOfvInit = np.array([[vtemp(i - 1,j)],[vtemp(i + 1,j)],[vtemp(i,j - 1)],[vtemp(i,j + 1)],[vtemp(i - 1,j - 1)],[vtemp(i + 1,j - 1)],[vtemp(i + 1,j + 1)],[vtemp(i - 1,j + 1)]])
                avgOfvInit[i,j] = (1 / 8 * sum(NeighborOfvInit(np.arange(1,4+1))) + 1 / 16 * sum(NeighborOfvInit(np.arange(5,8+1)))) / (3 / 4)
                stdOfvInit[i,j] = std(NeighborOfvInit,'omitnan')
                springErrOfu[i,j] = np.abs(utemp(i,j) - avgOfuInit(i,j)) - 3 * stdOfuInit(i,j)
                springErrOfv[i,j] = np.abs(vtemp(i,j) - avgOfvInit(i,j)) - 3 * stdOfvInit(i,j)
                if ((springErrOfu(i,j) > 0) or (springErrOfv(i,j) > 0)):
                    uInit[i,j] = NaN
                    vInit[i,j] = NaN
                    uInit[i - 1,j] = NaN
                    uInit[i + 1,j] = NaN
                    uInit[i,j - 1] = NaN
                    uInit[i,j + 1] = NaN
                    uInit[i - 1,j - 1] = NaN
                    uInit[i + 1,j - 1] = NaN
                    uInit[i + 1,j + 1] = NaN
                    uInit[i - 1,j + 1] = NaN
                    vInit[i - 1,j] = NaN
                    vInit[i + 1,j] = NaN
                    vInit[i,j - 1] = NaN
                    vInit[i,j + 1] = NaN
                    vInit[i - 1,j - 1] = NaN
                    vInit[i + 1,j - 1] = NaN
                    vInit[i + 1,j + 1] = NaN
                    vInit[i - 1,j + 1] = NaN
        # ##########PLOTTING##################
# figure; surf(springErrOfu)
# ##########PLOTTING##################
        uInit[1,:] = NaN
        uInit[end(),:] = NaN
        uInit[:,1] = NaN
        uInit[:,end()] = NaN
        vInit[1,:] = NaN
        vInit[end(),:] = NaN
        vInit[:,1] = NaN
        vInit[:,end()] = NaN
        uInit = inpaint_nans(uInit,4)
        vInit = inpaint_nans(vInit,4)
        # ##########PLOTTING##################
# figure;surf(uInit)
# colorbar
# figure;surf(vInit)
# colorbar
# ##########PLOTTING##################
        # Second time
        utemp = uInit
        vtemp = vInit
        Mtemp = u.shape[2-1]
        Ntemp = u.shape[1-1]
        for i in np.arange(2,Ntemp - 1+1).reshape(-1):
            for j in np.arange(2,Mtemp - 1+1).reshape(-1):
                NeighborOfuInit = np.array([[utemp(i - 1,j)],[utemp(i + 1,j)],[utemp(i,j - 1)],[utemp(i,j + 1)]])
                avgOfuInit[i,j] = mean(NeighborOfuInit,'omitnan')
                stdOfuInit[i,j] = std(NeighborOfuInit,'omitnan')
                NeighborOfvInit = np.array([[vtemp(i - 1,j)],[vtemp(i + 1,j)],[vtemp(i,j - 1)],[vtemp(i,j + 1)]])
                avgOfvInit[i,j] = mean(NeighborOfvInit,'omitnan')
                stdOfvInit[i,j] = std(NeighborOfvInit,'omitnan')
                springErrOfu[i,j] = np.abs(utemp(i,j) - avgOfuInit(i,j)) - 2 * stdOfuInit(i,j)
                springErrOfv[i,j] = np.abs(vtemp(i,j) - avgOfvInit(i,j)) - 2 * stdOfvInit(i,j)
                if ((springErrOfu(i,j) > 0) or (springErrOfv(i,j) > 0)):
                    uInit[i,j] = NaN
                    vInit[i,j] = NaN
                    uInit[i - 1,j] = NaN
                    uInit[i + 1,j] = NaN
                    uInit[i,j - 1] = NaN
                    uInit[i,j + 1] = NaN
                    # uInit(i-1,j-1) = NaN; uInit(i+1,j-1) = NaN; uInit(i+1,j+1) = NaN; uInit(i-1,j+1) = NaN;
                    vInit[i - 1,j] = NaN
                    vInit[i + 1,j] = NaN
                    vInit[i,j - 1] = NaN
                    vInit[i,j + 1] = NaN
                    # vInit(i-1,j-1) = NaN; vInit(i+1,j-1) = NaN; vInit(i+1,j+1) = NaN; vInit(i-1,j+1) = NaN;
        # ##########PLOTTING##################
# figure; surf(springErrOfu)
# ##########PLOTTING##################
        uInit[1,:] = NaN
        uInit[end(),:] = NaN
        uInit[:,1] = NaN
        uInit[:,end()] = NaN
        vInit[1,:] = NaN
        vInit[end(),:] = NaN
        vInit[:,1] = NaN
        vInit[:,end()] = NaN
        uInit = inpaint_nans(uInit,4)
        vInit = inpaint_nans(vInit,4)
        # ##########PLOTTING##################
# figure;surf(uInit)
# colorbar
# figure;surf(vInit)
# colorbar
# ##########PLOTTING##################
        # Third time
        utemp = uInit
        vtemp = vInit
        Mtemp = u.shape[2-1]
        Ntemp = u.shape[1-1]
        for i in np.arange(2,Ntemp - 1+1).reshape(-1):
            for j in np.arange(2,Mtemp - 1+1).reshape(-1):
                NeighborOfuInit = np.array([[utemp(i - 1,j)],[utemp(i + 1,j)],[utemp(i,j - 1)],[utemp(i,j + 1)],[utemp(i - 1,j - 1)],[utemp(i + 1,j - 1)],[utemp(i + 1,j + 1)],[utemp(i - 1,j + 1)]])
                avgOfuInit[i,j] = (1 / 8 * sum(NeighborOfuInit(np.arange(1,4+1))) + 1 / 16 * sum(NeighborOfuInit(np.arange(5,8+1)))) / (3 / 4)
                stdOfuInit[i,j] = std(NeighborOfuInit,'omitnan')
                NeighborOfvInit = np.array([[vtemp(i - 1,j)],[vtemp(i + 1,j)],[vtemp(i,j - 1)],[vtemp(i,j + 1)],[vtemp(i - 1,j - 1)],[vtemp(i + 1,j - 1)],[vtemp(i + 1,j + 1)],[vtemp(i - 1,j + 1)]])
                avgOfvInit[i,j] = (1 / 8 * sum(NeighborOfvInit(np.arange(1,4+1))) + 1 / 16 * sum(NeighborOfvInit(np.arange(5,8+1)))) / (3 / 4)
                stdOfvInit[i,j] = std(NeighborOfvInit,'omitnan')
                springErrOfu[i,j] = np.abs(utemp(i,j) - avgOfuInit(i,j)) - 2 * stdOfuInit(i,j)
                springErrOfv[i,j] = np.abs(vtemp(i,j) - avgOfvInit(i,j)) - 2 * stdOfvInit(i,j)
                if ((springErrOfu(i,j) > 0) or (springErrOfv(i,j) > 0)):
                    uInit[i,j] = NaN
                    vInit[i,j] = NaN
                    uInit[i - 1,j] = NaN
                    uInit[i + 1,j] = NaN
                    uInit[i,j - 1] = NaN
                    uInit[i,j + 1] = NaN
                    uInit[i - 1,j - 1] = NaN
                    uInit[i + 1,j - 1] = NaN
                    uInit[i + 1,j + 1] = NaN
                    uInit[i - 1,j + 1] = NaN
                    vInit[i - 1,j] = NaN
                    vInit[i + 1,j] = NaN
                    vInit[i,j - 1] = NaN
                    vInit[i,j + 1] = NaN
                    vInit[i - 1,j - 1] = NaN
                    vInit[i + 1,j - 1] = NaN
                    vInit[i + 1,j + 1] = NaN
                    vInit[i - 1,j + 1] = NaN
        # ##########PLOTTING##################
# figure; surf(springErrOfu)
# ##########PLOTTING##################
        uInit[1,:] = NaN
        uInit[end(),:] = NaN
        uInit[:,1] = NaN
        uInit[:,end()] = NaN
        vInit[1,:] = NaN
        vInit[end(),:] = NaN
        vInit[:,1] = NaN
        vInit[:,end()] = NaN
        uInit = inpaint_nans(uInit,4)
        vInit = inpaint_nans(vInit,4)
        # ##########PLOTTING##################
# figure;surf(uInit)
# colorbar
# figure;surf(vInit)
# colorbar
# ##########PLOTTING##################
        clear('Mtemp','Ntemp')
        clear('NeighborOfuInit','NeighborOfvInit','avgOfuInit','avgOfvInit','stdOfuInit','stdOfvInit','springErrOfu','springErrOfv')
    
    ## Initial Value
    U000 = np.zeros((2 * x.shape[1-1] * y.shape[2-1],1))
    Phi0 = np.zeros((x.shape[1-1] * y.shape[2-1],1))
    uInit = np.transpose(uInit)
    vInit = np.transpose(vInit)
    
    PhiInit = np.transpose(Phi)
    for i in np.arange(1,(x.shape[1-1] * y.shape[2-1])+1).reshape(-1):
        U000[2 * i - 1] = uInit(i)
        U000[2 * i] = vInit(i)
        Phi0[i] = PhiInit(i)
    
    print('Finish setting up mesh and assigning initial value!')
    U0 = U000
    return U0