## Integer Search

import numpy as np
import numpy.matlib
    
def funIntegerSearchMg(f = None,g = None,gridx = None,gridy = None,winsize = None,winstepsize = None,varargin = None): 
    gridxBackup = gridx
    gridyBackup = gridy
    
    borderGapx = 1 + np.round(0.75 * winsize)
    borderGapy = borderGapx
    if gridx(1) < borderGapx:
        gridx[1] = borderGapx
    
    if gridy(1) < borderGapy:
        gridy[1] = borderGapy
    
    if gridx(2) > f.shape[1-1] - borderGapx:
        gridx[2] = f.shape[1-1] - borderGapx
    
    if gridy(2) > f.shape[2-1] - borderGapy:
        gridy[2] = f.shape[2-1] - borderGapy
    
    # f=ImgNormalized{1}; g=ImgNormalized{2}; # (10:end,20:end)  ImgNormalized{1}(1:end-9,1:end-19); #
# gridx=gridxROIRange; gridy=gridyROIRange;
    
    ## Level 1
# First level cross-correlation
    C = f(np.arange(gridx(1),gridx(2)+1),np.arange(gridy(1),gridy(2)+1))
    D = g(np.arange(gridx(1),gridx(2)+1),np.arange(gridy(1),gridy(2)+1))
    XCORRF2OfCD0 = normxcorr2(C,D)
    
    # find maximum index of the cross-correlaiton
    v1temp,u1temp,max_f = findpeak(XCORRF2OfCD0,1)
    zero_disp = np.ceil((XCORRF2OfCD0.shape + np.array([1,1])) / 2)
    utemp = (u1temp - zero_disp(1))
    
    vtemp = (v1temp - zero_disp(2))
    
    Phitemp = max_f
    
    # Check out of image border or not
    gridx_f = gridx
    gridy_f = gridy
    gridx_g = gridx + np.array([np.ceil(utemp),np.ceil(utemp)])
    gridy_g = gridy + np.array([np.ceil(vtemp),np.ceil(vtemp)])
    
    if gridx_g(1) < borderGapx:
        temp = gridx_g(1)
        gridx_g[1] = borderGapx
        gridx_f[1] = gridx_f(1) + borderGapx - temp
    
    if gridy_g(1) < borderGapy:
        temp = gridy_g(1)
        gridy_g[1] = borderGapy
        gridy_f[1] = gridy_f(1) + borderGapy - temp
    
    if gridx_g(2) > f.shape[1-1] - borderGapx:
        temp = gridx_g(2)
        gridx_g[2] = f.shape[1-1] - borderGapx
        gridx_f[2] = gridx_f(2) + ((f.shape[1-1] - borderGapx) - temp)
    
    if gridy_g(2) > f.shape[2-1] - borderGapy:
        temp = gridy_g(2)
        gridy_g[2] = f.shape[2-1] - borderGapy
        gridy_f[2] = gridy_f(2) + (f.shape[2-1] - borderGapy) - temp
    
    # Make sure image width/length odd number
    if np.mod((gridx_f(2) - gridx_f(1)),2) == 1:
        gridx_f[2] = gridx_f(2) - 1
        gridx_g[2] = gridx_g(2) - 1
    
    if np.mod((gridy_f(2) - gridy_f(1)),2) == 1:
        gridy_f[2] = gridy_f(2) - 1
        gridy_g[2] = gridy_g(2) - 1
    
    # Visualize tracking
# figure,
# subplot(2,2,1); surf(f(gridxBackup(1):gridxBackup(2), gridyBackup(1):gridyBackup(2))','edgecolor','none'); view(2);axis equal;ylabel('y'); xlabel('x'); title('Raw f');
# subplot(2,2,2); surf(g(gridxBackup(1):gridxBackup(2), gridyBackup(1):gridyBackup(2))','edgecolor','none'); view(2);axis equal; ylabel('y'); xlabel('x'); title('Raw g');
# subplot(2,2,3); surf(f(gridx_f(1):gridx_f(2), gridy_f(1):gridy_f(2))','edgecolor','none'); view(2); axis equal;ylabel('y'); xlabel('x'); title('Shifted f');
# subplot(2,2,4); surf(g(gridx_g(1):gridx_g(2), gridy_g(1):gridy_g(2))','edgecolor','none'); view(2); axis equal;ylabel('y'); xlabel('x'); title('Shifted g');
    
    gridxWidth0 = gridx_f(2) - gridx_f(1)
    gridyWidth0 = gridy_f(2) - gridy_f(1)
    gridx_f0 = gridx_f
    gridy_f0 = gridy_f
    gridx_g0 = gridx_g
    gridy_g0 = gridy_g
    ## ====== Assign values ======
    gridxWidthCurr = gridxWidth0
    gridyWidthCurr = gridyWidth0
    gridxyRatioCurr = gridxWidthCurr / gridyWidthCurr
    gridx_fCurr = gridx_f0
    gridy_fCurr = gridy_f0
    gridx_gCurr = gridx_g0
    gridy_gCurr = gridy_g0
    utempCurr = np.ceil(utemp)
    vtempCurr = np.ceil(vtemp)
    ## Level >1
    levelNo = 1
    clear('gridx_fNew','gridx_gNew','gridy_fNew','gridy_gNew','utempNew','vtempNew','qfactors')
    TotalNo = 1
    gridxWidthNewtemp = gridxWidthCurr
    gridyWidthNewtemp = gridyWidthCurr
    while gridxWidthNewtemp / gridyWidthNewtemp > 0:

        if gridxWidthNewtemp / gridyWidthNewtemp > 2:
            gridxWidthNewtemp = gridxWidthNewtemp / 2
            gridyWidthNewtemp = gridyWidthNewtemp
            TotalNo = TotalNo * 2
        else:
            if gridxWidthNewtemp / gridyWidthNewtemp < 0.5:
                gridxWidthNewtemp = gridxWidthNewtemp
                gridyWidthNewtemp = gridyWidthNewtemp / 2
                TotalNo = TotalNo * 2
            else:
                gridxWidthNewtemp = gridxWidthNewtemp / 2
                gridyWidthNewtemp = gridyWidthNewtemp / 2
                TotalNo = TotalNo * 4
        if (gridxWidthNewtemp < winsize and gridyWidthNewtemp < winsize):
            break

    
    #TotalNo = 3*(gridx(2)-gridx(1))*(gridy(2)-gridy(1))/winsize^2;
    IterNo = 0
    hbar = waitbar(0,'FFT initial guess, it is fast and please wait.')
    while gridxyRatioCurr > 0:

        levelNo = levelNo + 1
        clear('utempNew','vtempNew','Phitemp')
        # ====== Split gridx & gridy ======
        if gridxyRatioCurr > 2:
            gridxWidthNew = gridxWidthCurr / 2
            gridyWidthNew = gridyWidthCurr
            for tempi in np.arange(1,gridx_fCurr.shape[1-1]+1).reshape(-1):
                tempj = tempi
                tempInd = tempj
                if ((gridx_fCurr(tempInd,2) - gridx_fCurr(tempInd,1) > winsize)):
                    gridx_fNew[np.arange[2 * tempInd - 1,2 * tempInd+1],np.arange[1,2+1]] = np.array([[gridx_fCurr(tempInd,1),0.5 * sum(gridx_fCurr(tempInd,:))],[0.5 * sum(gridx_fCurr(tempInd,:)),gridx_fCurr(tempInd,2)]])
                    gridx_gNew[np.arange[2 * tempInd - 1,2 * tempInd+1],np.arange[1,2+1]] = np.array([[gridx_gCurr(tempInd,1),0.5 * sum(gridx_gCurr(tempInd,:))],[0.5 * sum(gridx_gCurr(tempInd,:)),gridx_gCurr(tempInd,2)]])
                    gridy_fNew[np.arange[2 * tempInd - 1,2 * tempInd+1],np.arange[1,2+1]] = np.matlib.repmat(gridy_fCurr(tempInd,np.arange(1,2+1)),2,1)
                    gridy_gNew[np.arange[2 * tempInd - 1,2 * tempInd+1],np.arange[1,2+1]] = np.matlib.repmat(gridy_gCurr(tempInd,np.arange(1,2+1)),2,1)
                    utempNew[np.arange[2 * tempInd - 1,2 * tempInd+1]] = np.matlib.repmat(utempCurr(tempInd),2,1)
                    vtempNew[np.arange[2 * tempInd - 1,2 * tempInd+1]] = np.matlib.repmat(vtempCurr(tempInd),2,1)
                else:
                    gridx_fNew[np.arange[2 * tempInd - 1,2 * tempInd+1],np.arange[1,2+1]] = np.matlib.repmat(gridx_fCurr(tempInd,np.arange(1,2+1)),2,1)
                    gridx_gNew[np.arange[2 * tempInd - 1,2 * tempInd+1],np.arange[1,2+1]] = np.matlib.repmat(gridx_gCurr(tempInd,np.arange(1,2+1)),2,1)
                    gridy_fNew[np.arange[2 * tempInd - 1,2 * tempInd+1],np.arange[1,2+1]] = np.matlib.repmat(gridy_fCurr(tempInd,np.arange(1,2+1)),2,1)
                    gridy_gNew[np.arange[2 * tempInd - 1,2 * tempInd+1],np.arange[1,2+1]] = np.matlib.repmat(gridy_gCurr(tempInd,np.arange(1,2+1)),2,1)
                    utempNew[np.arange[2 * tempInd - 1,2 * tempInd+1]] = np.matlib.repmat(utempCurr(tempInd),2,1)
                    vtempNew[np.arange[2 * tempInd - 1,2 * tempInd+1]] = np.matlib.repmat(vtempCurr(tempInd),2,1)
        else:
            if gridxyRatioCurr < 0.5:
                gridxWidthNew = gridxWidthCurr
                gridyWidthNew = gridyWidthCurr / 2
                for tempi in np.arange(1,gridx_fCurr.shape[1-1]+1).reshape(-1):
                    tempj = tempi
                    tempInd = tempj
                    if ((gridy_fCurr(tempInd,2) - gridy_fCurr(tempInd,1) > winsize)):
                        gridx_fNew[np.arange[2 * tempInd - 1,2 * tempInd+1],np.arange[1,2+1]] = np.matlib.repmat(gridx_fCurr(tempInd,np.arange(1,2+1)),2,1)
                        gridx_gNew[np.arange[2 * tempInd - 1,2 * tempInd+1],np.arange[1,2+1]] = np.matlib.repmat(gridx_gCurr(tempInd,np.arange(1,2+1)),2,1)
                        gridy_fNew[np.arange[2 * tempInd - 1,2 * tempInd+1],np.arange[1,2+1]] = np.array([[gridy_fCurr(tempInd,1),0.5 * sum(gridy_fCurr(tempInd,:))],[0.5 * sum(gridy_fCurr(tempInd,:)),gridy_fCurr(tempInd,2)]])
                        gridy_gNew[np.arange[2 * tempInd - 1,2 * tempInd+1],np.arange[1,2+1]] = np.array([[gridy_gCurr(tempInd,1),0.5 * sum(gridy_gCurr(tempInd,:))],[0.5 * sum(gridy_gCurr(tempInd,:)),gridy_gCurr(tempInd,2)]])
                        utempNew[np.arange[2 * tempInd - 1,2 * tempInd+1]] = np.matlib.repmat(utempCurr(tempInd),2,1)
                        vtempNew[np.arange[2 * tempInd - 1,2 * tempInd+1]] = np.matlib.repmat(vtempCurr(tempInd),2,1)
                    else:
                        gridx_fNew[np.arange[2 * tempInd - 1,2 * tempInd+1],np.arange[1,2+1]] = np.matlib.repmat(gridx_fCurr(tempInd,np.arange(1,2+1)),2,1)
                        gridx_gNew[np.arange[2 * tempInd - 1,2 * tempInd+1],np.arange[1,2+1]] = np.matlib.repmat(gridx_gCurr(tempInd,np.arange(1,2+1)),2,1)
                        gridy_fNew[np.arange[2 * tempInd - 1,2 * tempInd+1],np.arange[1,2+1]] = np.matlib.repmat(gridy_fCurr(tempInd,np.arange(1,2+1)),2,1)
                        gridy_gNew[np.arange[2 * tempInd - 1,2 * tempInd+1],np.arange[1,2+1]] = np.matlib.repmat(gridy_gCurr(tempInd,np.arange(1,2+1)),2,1)
                        utempNew[np.arange[2 * tempInd - 1,2 * tempInd+1]] = np.matlib.repmat(utempCurr(tempInd),2,1)
                        vtempNew[np.arange[2 * tempInd - 1,2 * tempInd+1]] = np.matlib.repmat(vtempCurr(tempInd),2,1)
            else:
                gridxWidthNew = gridxWidthCurr / 2
                gridyWidthNew = gridyWidthCurr / 2
                for tempi in np.arange(1,gridx_fCurr.shape[1-1]+1).reshape(-1):
                    tempj = tempi
                    tempInd = tempj
                    if ((gridx_fCurr(tempInd,2) - gridx_fCurr(tempInd,1) > winsize) and (gridy_fCurr(tempInd,2) - gridy_fCurr(tempInd,1) > winsize)):
                        gridx_fNew[np.arange[4 * tempInd - 3,4 * tempInd+1],np.arange[1,2+1]] = np.matlib.repmat(np.array([[gridx_fCurr(tempInd,1),0.5 * sum(gridx_fCurr(tempInd,:))],[0.5 * sum(gridx_fCurr(tempInd,:)),gridx_fCurr(tempInd,2)]]),2,1)
                        gridx_gNew[np.arange[4 * tempInd - 3,4 * tempInd+1],np.arange[1,2+1]] = np.matlib.repmat(np.array([[gridx_gCurr(tempInd,1),0.5 * sum(gridx_gCurr(tempInd,:))],[0.5 * sum(gridx_gCurr(tempInd,:)),gridx_gCurr(tempInd,2)]]),2,1)
                        gridy_fNew[np.arange[4 * tempInd - 3,4 * tempInd+1],np.arange[1,2+1]] = np.array([[np.matlib.repmat(np.array([gridy_fCurr(tempInd,1),0.5 * sum(gridy_fCurr(tempInd,:))]),2,1)],[np.matlib.repmat(np.array([0.5 * sum(gridy_fCurr(tempInd,:)),gridy_fCurr(tempInd,2)]),2,1)]])
                        gridy_gNew[np.arange[4 * tempInd - 3,4 * tempInd+1],np.arange[1,2+1]] = np.array([[np.matlib.repmat(np.array([gridy_gCurr(tempInd,1),0.5 * sum(gridy_gCurr(tempInd,:))]),2,1)],[np.matlib.repmat(np.array([0.5 * sum(gridy_gCurr(tempInd,:)),gridy_gCurr(tempInd,2)]),2,1)]])
                        utempNew[np.arange[4 * tempInd - 3,4 * tempInd+1]] = np.matlib.repmat(utempCurr(tempInd),4,1)
                        vtempNew[np.arange[4 * tempInd - 3,4 * tempInd+1]] = np.matlib.repmat(vtempCurr(tempInd),4,1)
                    else:
                        gridx_fNew[np.arange[4 * tempInd - 3,4 * tempInd+1],np.arange[1,2+1]] = np.matlib.repmat(gridx_fCurr(tempInd,np.arange(1,2+1)),4,1)
                        gridx_gNew[np.arange[4 * tempInd - 3,4 * tempInd+1],np.arange[1,2+1]] = np.matlib.repmat(gridx_gCurr(tempInd,np.arange(1,2+1)),4,1)
                        gridy_fNew[np.arange[4 * tempInd - 3,4 * tempInd+1],np.arange[1,2+1]] = np.matlib.repmat(gridy_fCurr(tempInd,np.arange(1,2+1)),4,1)
                        gridy_gNew[np.arange[4 * tempInd - 3,4 * tempInd+1],np.arange[1,2+1]] = np.matlib.repmat(gridy_gCurr(tempInd,np.arange(1,2+1)),4,1)
                        utempNew[np.arange[4 * tempInd - 3,4 * tempInd+1]] = np.matlib.repmat(utempCurr(tempInd),4,1)
                        vtempNew[np.arange[4 * tempInd - 3,4 * tempInd+1]] = np.matlib.repmat(vtempCurr(tempInd),4,1)
        for tempi in np.arange(1,gridx_fNew.shape[1-1]+1).reshape(-1):
            IterNo = IterNo + 1
            waitbar(IterNo / TotalNo)
            C = f(np.arange(gridx_fNew(tempi,1),gridx_fNew(tempi,2)+1),np.arange(gridy_fNew(tempi,1),gridy_fNew(tempi,2)+1))
            D = g(np.arange(gridx_gNew(tempi,1) - 5,gridx_gNew(tempi,2) + 5+1),np.arange(gridy_gNew(tempi,1) - 5,gridy_gNew(tempi,2) + 5+1))
            winsize1 = C.shape[1-1] - 1
            winsize2 = C.shape[2-1] - 1
            XCORRF2OfCD0 = normxcorr2(C,D)
            # find qfactors
            cc.A[0] = real(XCORRF2OfCD0)
            qfactors[tempi,:] = compute_qFactor(cc,tempi)
            v1temp,u1temp,max_f = findpeak(XCORRF2OfCD0(np.arange(winsize1,end() - winsize1 + 1+1),np.arange(winsize2,end() - winsize2 + 2+1)),1)
            zero_disp = np.ceil(XCORRF2OfCD0(np.arange(winsize1,end() - winsize1 + 1+1),np.arange(winsize2,end() - winsize2 + 1+1)).shape / 2)
            ind = tempi
            utempNew[ind] = utempNew(ind) + (u1temp - zero_disp(1))
            vtempNew[ind] = vtempNew(ind) + (v1temp - zero_disp(2))
            Phitemp[ind] = max_f
        if (gridxWidthNew < winsize and gridyWidthNew < winsize):
            # Finish qfactor computation
            for k in np.arange(1,2+1).reshape(-1):
                qf_ = (qfactors(:,k) - np.amin(qfactors(:,k)))
                cc.qfactors[:,k] = qf_ / np.amax(qf_)
            break
        else:
            ## ====== Check out of image border or not ======
            gridx_f = gridx_fNew
            gridy_f = gridy_fNew
            gridx_g = gridx_fNew + np.transpose(np.ceil(utempNew)) * np.array([1,1])
            gridy_g = gridy_fNew + np.transpose(np.ceil(vtempNew)) * np.array([1,1])
            for tempi in np.arange(1,gridx_g.shape[1-1]+1).reshape(-1):
                if gridx_g(tempi,1) < borderGapx:
                    temp = gridx_g(tempi,1)
                    gridx_g[tempi,1] = borderGapx
                    gridx_f[tempi,1] = gridx_f(tempi,1) + borderGapx - temp
                if gridy_g(tempi,1) < borderGapy:
                    temp = gridy_g(tempi,1)
                    gridy_g[tempi,1] = borderGapy
                    gridy_f[tempi,1] = gridy_f(tempi,1) + borderGapy - temp
                if gridx_g(tempi,2) > f.shape[1-1] - borderGapx:
                    temp = gridx_g(tempi,2)
                    gridx_g[tempi,2] = f.shape[1-1] - borderGapx
                    gridx_f[tempi,2] = gridx_f(tempi,2) + (f.shape[1-1] - borderGapx) - temp
                if gridy_g(tempi,2) > f.shape[2-1] - borderGapy:
                    temp = gridy_g(tempi,2)
                    gridy_g[tempi,2] = f.shape[2-1] - borderGapy
                    gridy_f[tempi,2] = gridy_f(tempi,2) + (f.shape[2-1] - borderGapy) - temp
            # Make sure: image width/length odd number
            for tempi in np.arange(1,gridx_g.shape[1-1]+1).reshape(-1):
                if np.mod(gridx_f(tempi,2) - gridx_f(tempi,1),2) == 1:
                    gridx_f[tempi,2] = gridx_f(tempi,2) - 1
                    gridx_g[tempi,2] = gridx_g(tempi,2) - 1
                if np.mod(gridy_f(tempi,2) - gridy_f(tempi,1),2) == 1:
                    gridy_f[tempi,2] = gridy_f(tempi,2) - 1
                    gridy_g[tempi,2] = gridy_g(tempi,2) - 1
            ## ====== Assign values ======
            gridxWidthCurr = gridxWidthNew
            gridyWidthCurr = gridyWidthNew
            gridxyRatioCurr = gridxWidthCurr / gridyWidthCurr
            gridx_fCurr = gridx_f
            gridy_fCurr = gridy_f
            gridx_gCurr = gridx_g
            gridy_gCurr = gridy_g
            utempCurr = np.ceil(utempNew)
            vtempCurr = np.ceil(vtempNew)

    
    close_(hbar)
    ##
    tempx = 0.5 * (gridx_fNew(:,1) + gridx_fNew(:,2))
    tempy = 0.5 * (gridy_fNew(:,1) + gridy_fNew(:,2))
    tempu = np.transpose(utempNew)
    tempv = np.transpose(vtempNew)
    tempPhi = np.transpose(Phitemp)
    qf_1 = np.transpose(cc.qfactors(:,1))
    qf_2 = np.transpose(cc.qfactors(:,2))
    # figure, plot3(tempx,tempy,tempu,'.');
# figure, plot3(tempx,tempy,tempv,'.');
    
    ## Update image domain
    __,indx1 = np.amin(tempx)
    __,indx2 = np.amax(tempx)
    __,indy1 = np.amin(tempy)
    __,indy2 = np.amax(tempy)
    borderGapx1 = np.ceil(1 + (1.25 * winsize) + ((tempu(indx1))))
    borderGapx2 = np.ceil(1 + (1.25 * winsize) + ((tempu(indx2))))
    borderGapy1 = np.ceil(1 + (1.25 * winsize) + ((tempv(indy1))))
    borderGapy2 = np.ceil(1 + (1.25 * winsize) + ((tempv(indy2))))
    if gridx(1) < borderGapx1:
        gridx[1] = borderGapx1
    
    if gridy(1) < borderGapy1:
        gridy[1] = borderGapy1
    
    if gridx(2) > f.shape[1-1] - borderGapx2 + 1:
        gridx[2] = f.shape[1-1] - borderGapx2 + 1
    
    if gridy(2) > f.shape[2-1] - borderGapy2 + 1:
        gridy[2] = f.shape[2-1] - borderGapy2 + 1
    
    ## Interpolate
    try:
        xList = np.array([np.arange(gridx(1),gridx(2)+winstepsize(1),winstepsize(1))])
    finally:
        pass
    
    try:
        yList = np.array([np.arange(gridy(1),gridy(2)+winstepsize(2),winstepsize(2))])
    finally:
        pass
    
    indx = find(tempx > np.logical_and(np.amin(xList),tempx) < np.amax(xList))
    indy = find(tempy > np.logical_and(np.amin(yList),tempy) < np.amax(yList))
    indxy = intersect(indx,indy)
    uGrid,xGrid,yGrid = gridfit(tempx(indxy),tempy(indxy),tempu(indxy),xList,yList,'regularizer','springs')
    vGrid,__,__ = gridfit(tempx(indxy),tempy(indxy),tempv(indxy),xList,yList,'regularizer','springs')
    PhiGrid,__,__ = gridfit(tempx(indxy),tempy(indxy),tempPhi(indxy),xList,yList,'regularizer','springs')
    qf_1Grid,__,__ = gridfit(tempx(indxy),tempy(indxy),qf_1(indxy),xList,yList,'regularizer','springs')
    qf_2Grid,__,__ = gridfit(tempx(indxy),tempy(indxy),qf_2(indxy),xList,yList,'regularizer','springs')
    cc.max = PhiGrid
    cc.A = []
    cc.qfactors = np.array([qf_1Grid,qf_2Grid])
    # figure, surf(xGrid,yGrid,uGrid,'edgecolor','none')
# figure, surf(xGrid,yGrid,vGrid,'edgecolor','none')
    
    return xGrid,yGrid,uGrid,vGrid,cc
    
    
def funIntegerSearchWholeField(f = None,g = None,tempSizeOfSearchRegion = None,gridx = None,gridy = None,winsize = None,winstepsize = None): 
    #cj1 = 1; ci1 = 1; # index to count main loop
    
    if len(winstepsize) == 1:
        winstepsize = np.matlib.repmat(winstepsize,1,2)
    
    if len(winsize) == 1:
        winsize = np.matlib.repmat(winsize,1,2)
    
    # disp('Assemble point position sequence.');
    XList = np.array([np.arange(gridx(1),gridx(2) - winstepsize(1)+winstepsize(1),winstepsize(1))])
    YList = np.array([np.arange(gridy(1),gridy(2) - winstepsize(2)+winstepsize(2),winstepsize(2))])
    XX,YY = ndgrid(XList,YList)
    temparrayLength = len(XList) * len(YList)
    PtPosSeq = np.zeros((temparrayLength,2))
    PtPosSeq[:,1] = XX
    PtPosSeq[:,2] = YY
    cj1temp = np.zeros((temparrayLength,1))
    ci1temp = cj1temp
    utemp = cj1temp
    vtemp = cj1temp
    xtemp = cj1temp
    ytemp = cj1temp
    Phitemp = cj1temp
    ## ========== Start initial integer search ==========
    hbar = waitbar(0,'FFT initial guess, it is fast and please wait.')
    # hbar = parfor_progressbar(temparrayLength,'Please wait for integer search!');
    
    # sizeOfx1 = floor((gridx(2)-gridx(1)-winsize)/winstepsize)+1;
# sizeOfx2 = floor((gridy(2)-gridy(1)-winsize)/winstepsize)+1;
# disp(['Init search using FFT cross correlation on grid: ',num2str(sizeOfx1),'x',num2str(sizeOfx2)]);
    x = np.zeros((len(YList),len(XList)))
    y = x
    u = x
    v = x
    Phi = x
    for tempi in np.arange(1,temparrayLength+1).reshape(-1):
        waitbar(tempi / temparrayLength)
        jj = PtPosSeq(tempi,2)
        # jj is for y -or- vertical direction of images
        ii = PtPosSeq(tempi,1)
        # ii is for x -or- horizontal direction of images
        C = f(np.arange(ii,ii + winsize(1)+1),np.arange(jj,jj + winsize(2)+1))
        D = g(np.arange(ii - tempSizeOfSearchRegion,ii + winsize(1) + tempSizeOfSearchRegion+1),np.arange(jj - tempSizeOfSearchRegion,jj + winsize(2) + tempSizeOfSearchRegion+1))
        XCORRF2OfCD0 = normxcorr2(C,D)
        #find qfactors
        cc.A[0] = real(XCORRF2OfCD0)
        qfactors[tempi,:] = compute_qFactor(cc,tempi)
        # find maximum index of the cross-correlaiton
        v1temp,u1temp,max_f = findpeak(XCORRF2OfCD0(np.arange(winsize(1),end() - winsize(1) + 1+1),np.arange(winsize(2),end() - winsize(2) + 1+1)),1)
        ######################################################################
# I tried following code, unfortunately it doesn't work very well
# [v1temp1, u1temp1, max_f1] = findpeak(XCORRF2OfCD0(winsize:end-winsize+1,winsize:end-winsize+1),1);
# [v1temp2, u1temp2, max_f2] = findpeak(XCORRF2OfCD0(winsize:end-winsize+1,winsize:end-winsize+1),0);
        # if max_f2 > 0.999
#    v1temp = v1temp2; u1temp = u1temp2; max_f = max_f2;
# else
#    v1temp = v1temp1; u1temp = u1temp1; max_f = max_f1;
# end
######################################################################
        zero_disp = np.ceil(XCORRF2OfCD0(np.arange(winsize(1),end() - winsize(1) + 1+1),np.arange(winsize(2),end() - winsize(2) + 1+1)).shape / 2)
        utemp[tempi] = u1temp - zero_disp(1)
        vtemp[tempi] = v1temp - zero_disp(2)
        Phitemp[tempi] = max_f
        ytemp[tempi] = (jj + jj + winsize(2)) / 2
        xtemp[tempi] = (ii + ii + winsize(1)) / 2
        # Update counters
#ci1 = ci1 + 1;  # ci1 is moving horizontally for subsets
        #end
        #ci1=1; cj1=cj1+1;  # cj1 is moving vertically for subsets
        cj1temp[tempi] = np.ceil(tempi / len(XList))
        ci1temp[tempi] = tempi - (cj1temp(tempi) - 1) * len(XList)
    
    close_(hbar)
    for k in np.arange(1,2+1).reshape(-1):
        qf_ = (qfactors(:,k) - np.amin(qfactors(:,k)))
        cc.qfactors[:,k] = qf_ / np.amax(qf_)
    
    # hbar = parfor_progressbar(temparrayLegTotal,'Assign results to variables.');
    hbar = waitbar(0,'Assign results to variables.')
    for tempi in np.arange(1,temparrayLength+1).reshape(-1):
        ci1 = ci1temp(tempi)
        cj1 = cj1temp(tempi)
        u[cj1,ci1] = utemp(tempi)
        v[cj1,ci1] = vtemp(tempi)
        Phi[cj1,ci1] = Phitemp(tempi)
        x[cj1,ci1] = xtemp(tempi)
        y[cj1,ci1] = ytemp(tempi)
        waitbar(tempi / temparrayLength)
        # hbar.iterate(1);
    
    close_(hbar)
    print('Finish initial guess search!')
    # -------- End of Local integer search --------
    
    cc.max = Phi
    cc.A = []
    return x,y,u,v,cc
    
    ## ==============================================
    
    ##
    
def compute_qFactor(cc = None,qnum = None): 
    #get peak locations and cc_min maps (i.e. cc - cc(min))
    peak,cc_min = cellfun(lambda x = None: cc_max_find(double(x)),cc.A,'UniformOutput',0)
    #compute two primary quality metrics, as given in "Xue Z, Particle Image
# Velocimetry Correlation Signal-to-noise Metrics, Particle Image
# Pattern Mutual Information and Measurement uncertainty Quantification.
# MS Thesis, Virginia Tech, 2014.
    
    #peak to corr. energy ratio
    pce = cellfun(lambda x = None,y = None: (np.abs(y) ** 2) / (1 / np.asarray(x).size * (sum(np.abs(x) ** 2))),cc_min,peak,'UniformOutput',0)
    #min value -> 1 (worst case)
    
    #peak to entropy ratio
    ppe = cellfun(lambda x = None: q_entropy(double(x)),cc_min,'UniformOutput',0)
    
    #min value -> 0 (worst case)
    
    qfactors = np.transpose(cell2mat(cellfun(lambda x = None,y = None: np.array([[x],[y]]),pce,ppe,'UniformOutput',0)))
    return qfactors
    
    
def cc_max_find(cc = None): 
    #find the peak and zero-adjusted cc map
    
    cc_min = cc - np.amin(cc)
    
    # cc_filt = imgaussfilt3(cc_min); #filter to remove noise from peak value
    
    peak,__ = np.amax(cc_min)
    
    return peak,cc_min
    
    
def q_entropy(cc_min = None): 
    #compute entropy q-factor for a given cc map
    
    cc_hist,__ = histcounts(cc_min,30)
    
    entropy = 0
    p = cc_hist / sum(cc_hist)
    
    for i in np.arange(1,len(p)+1).reshape(-1):
        if p(i) == 0:
            entropy = entropy + p(i)
        else:
            entropy = entropy + p(i) * np.log(1 / p(i))
    
    ppe = 1 / entropy
    
    #min value -> 0 (worst case)
    
    return ppe
    
    return xGrid,yGrid,uGrid,vGrid,cc