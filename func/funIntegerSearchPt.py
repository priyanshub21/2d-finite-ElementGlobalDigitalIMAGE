

## Integer Search on a non-uniform mesh

import numpy as np
    
def funIntegerSearchPt(f = None,g = None,winsize = None,varargin = None): 
    seedPtCoords = varargin[0]
    uSeedPt = np.zeros((seedPtCoords.shape[1-1],1))
    vSeedPt = uSeedPt
    PhiSeedPt = uSeedPt
    print('--- The size of initial guess search zone (pixels)? ---  \n' % ())
    print('User should start to try a small integer value, and gradually increase the value of \n' % ())
    print('the search zone size until it is larger than the magnitudes of |disp u| and |disp v|. \n' % ())
    print('User could also input [size_x, size_y] to search in a rectangular zone. \n' % ())
    prompt = 'Input here: '
    tempSizeOfSearchRegion = input_(prompt)
    # if length(tempSizeOfSearchRegion) == 1, tempSizeOfSearchRegion = tempSizeOfSearchRegion*[1,1]; end
    
    for tempi in np.arange(1,len(seedPtCoords)+1).reshape(-1):
        if seedPtCoords(tempi,1) > 0:
            if np.ceil(seedPtCoords(tempi,1) - winsize / 2) - tempSizeOfSearchRegion < 1 or int(np.floor(seedPtCoords(tempi,1) + winsize / 2)) + tempSizeOfSearchRegion > g.shape[1-1] or np.ceil(seedPtCoords(tempi,2) - winsize / 2) - tempSizeOfSearchRegion < 1 or int(np.floor(seedPtCoords(tempi,2) + winsize / 2)) + tempSizeOfSearchRegion > g.shape[2-1]:
                uSeedPt[tempi] = nan
                vSeedPt[tempi] = nan
                PhiSeedPt[tempi] = nan
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
        else:
            uSeedPt[tempi] = nan
            vSeedPt[tempi] = nan
            PhiSeedPt[tempi] = nan
    
    return uSeedPt,vSeedPt,PhiSeedPt,tempSizeOfSearchRegion
    
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
    
    return uSeedPt,vSeedPt,PhiSeedPt,tempSizeOfSearchRegion