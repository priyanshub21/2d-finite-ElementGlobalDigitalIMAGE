## Integer Search

import numpy as np
import numpy.matlib
    
def funIntegerSearch(f = None,g = None,tempSizeOfSearchRegion = None,gridx = None,gridy = None,winsize = None,winstepsize = None,tempNoOfInitPt = None,varargin = None): 
    #(f,g,tempSizeOfSearchRegion,gridx,gridy,winsize,winstepsize)
    
    if 0 == tempNoOfInitPt:
        gridxBackup = gridx
        gridyBackup = gridy
        # may lose boundary regions
        while gridx(1) - tempSizeOfSearchRegion(1) - 0.5 * winsize < 1:

            gridx[1] = gridx(1) + 1

        while gridy(1) - tempSizeOfSearchRegion(1) - 0.5 * winsize < 1:

            gridy[1] = gridy(1) + 1

        while gridx(end()) + 1.5 * winsize + tempSizeOfSearchRegion(2) > f.shape[1-1]:

            gridx[end()] = gridx(end()) - 1

        while gridy(end()) + 1.5 * winsize + tempSizeOfSearchRegion(2) > f.shape[2-1]:

            gridy[end()] = gridy(end()) - 1

        x0,y0,u0,v0,cc0 = funIntegerSearchWholeField(f,g,tempSizeOfSearchRegion,gridx,gridy,winsize,winstepsize)
        x = x0
        y = y0
        u = u0
        v = v0
        cc = cc0
        #         xList = gridxBackup(1)+4:winstepsize:gridxBackup(end)-4;
#         yList = gridyBackup(1)+4:winstepsize:gridyBackup(end)-4;
#         [x,y] = meshgrid(xList,yList);
#         u = gridfit(x0,y0,u0,xList,yList,'smooth',1,'interp','bilinear','regularizer','springs');
#         v = gridfit(x0,y0,v0,xList,yList,'smooth',1,'interp','bilinear','regularizer','springs');
#         cc.max = gridfit(x0,y0,cc0.max,xList,yList,'smooth',1,'interp','bilinear','regularizer','springs');
    else:
        if 1 == tempNoOfInitPt:
            seedPtCoords = varargin[0]
            uSeedPt = np.zeros((seedPtCoords.shape[1-1],1))
            vSeedPt = uSeedPt
            PhiSeedPt = uSeedPt
            for tempi in np.arange(1,len(seedPtCoords)+1).reshape(-1):
                if np.ceil(seedPtCoords(tempi,1) - winsize / 2) - tempSizeOfSearchRegion < 1 or int(np.floor(seedPtCoords(tempi,1) + winsize / 2)) + tempSizeOfSearchRegion > g.shape[1-1] or np.ceil(seedPtCoords(tempi,2) - winsize / 2) - tempSizeOfSearchRegion < 1 or int(np.floor(seedPtCoords(tempi,2) + winsize / 2)) + tempSizeOfSearchRegion > g.shape[2-1]:
                    continue
                else:
                    C = f(np.arange(np.ceil(seedPtCoords(tempi,1) - winsize / 2),int(np.floor(seedPtCoords(tempi,1) + winsize / 2))+1),np.arange(np.ceil(seedPtCoords(tempi,2) - winsize / 2),int(np.floor(seedPtCoords(tempi,2) + winsize / 2))+1))
                    D = g(np.arange(np.ceil(seedPtCoords(tempi,1) - winsize / 2) - tempSizeOfSearchRegion,int(np.floor(seedPtCoords(tempi,1) + winsize / 2)) + tempSizeOfSearchRegion+1),np.arange(np.ceil(seedPtCoords(tempi,2) - winsize / 2) - tempSizeOfSearchRegion,int(np.floor(seedPtCoords(tempi,2) + winsize / 2)) + tempSizeOfSearchRegion+1))
                    XCORRF2OfCD0 = normxcorr2(C,D)
                    v1temp,u1temp,max_f = findpeak(XCORRF2OfCD0(np.arange(winsize,end() - winsize + 1+1),np.arange(winsize,end() - winsize + 1+1)),1)
                    zero_disp = np.ceil(XCORRF2OfCD0(np.arange(winsize,end() - winsize + 1+1),np.arange(winsize,end() - winsize + 1+1)).shape / 2)
                    uSeedPt[tempi] = u1temp - zero_disp(1)
                    vSeedPt[tempi] = v1temp - zero_disp(2)
                    PhiSeedPt[tempi] = max_f
            xList = np.arange(gridx(1),gridx(end())+winstepsize,winstepsize)
            yList = np.arange(gridy(1),gridy(end())+winstepsize,winstepsize)
            x,y = np.meshgrid(xList,yList)
            #bc = Spline2D('bicubic',X,Y,ZZ);
            u = gridfit(seedPtCoords(:,1),seedPtCoords(:,2),uSeedPt,xList,yList,'smooth',100,'interp','bicubic','regularizer','springs')
            v = gridfit(seedPtCoords(:,1),seedPtCoords(:,2),vSeedPt,xList,yList,'smooth',100,'interp','bicubic','regularizer','springs')
            cc.max = gridfit(seedPtCoords(:,1),seedPtCoords(:,2),PhiSeedPt,xList,yList,'smooth',100,'interp','bicubic','regularizer','springs')
        else:
            print('wrong input in IntegerSearch!')
    
    return x,y,u,v,cc
    
    
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
    
    return x,y,u,v,cc