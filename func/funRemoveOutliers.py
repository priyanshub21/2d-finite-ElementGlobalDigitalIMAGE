import numpy as np
import matplotlib.pyplot as plt
    
def funRemoveOutliers(u = None,v = None,cc = None,qDICOrNot = None,Thr0 = None,varargin = None): 
    # =========================================================================
# removes outliers using the universal
# outlier test based on
    
    # J. Westerweel and F. Scarano. Universal outlier detection for PIV data.
# Exp. Fluids, 39(6):1096{1100, August 2005. doi: 10.1007/s00348-005-0016-6
# -------------------------------------------------------------------------
# NOTES
# -------------------------------------------------------------------------
# needs medFilt3 and John D'Errico's inpaint_nans3
# (http://www.mathworks.com/matlabcentral/fileexchange/4551-inpaint-nans)function.
# =========================================================================
    
    if 7 == len(varargin):
        BadptRow = varargin[0]
        BadptCol = varargin[2]
    else:
        BadptRow = []
        BadptCol = []
    
    # ============== qDIC bad points removal ===============
# prompt = 'Input threshold for qDIC bad points removal:';
# ccThreshold = input(prompt); # Thr = 50;
    
    if qDICOrNot == 1:
        M,N = u.shape
        mSize = np.array([M * N,1])
        cc,ccMask_ = removeBadCorrelations(cc,cc.ccThreshold,1,mSize)
        for ii in np.arange(1,2+1).reshape(-1):
            ccMask[ii] = np.reshape(double(ccMask_(:,ii)), tuple(mSize), order="F")
        qDICpceRmList = ismissing(ccMask[0])
        qDICppeRmList = ismissing(ccMask[2])
    else:
        qDICpceRmList = []
        qDICppeRmList = []
    
    # ============== median test bad points removal ===============
    medianU = cell(1,2)
    normFluct = cell(1,2)
    normFluctMag = np.zeros((u.shape,u.shape))
    epsilon = 0.1
    medianU[0],normFluct[0] = funMedRemoveOutliers(u,epsilon)
    medianU[2],normFluct[2] = funMedRemoveOutliers(v,epsilon)
    normFluctMag = normFluct[0] ** 2 + normFluct[2] ** 2
    normFluctMag = np.sqrt(normFluctMag)
    MedFilterOrNot = 0
    while MedFilterOrNot < 1:

        figure
        surf(normFluctMag)
        plt.axis('equal')
        plt.axis('tight')
        view(2)
        caxis('auto')
        colorbar
        if len(Thr0)==0 or (Thr0 == 0):
            prompt = '--- Threshold for median test --- Input here: '
            Thr = input_(prompt)
        else:
            Thr = Thr0
        RemoveOutliersList = find(normFluctMag > Thr)
        # ============== remove bad points ===============
        u2 = u
        u2[qDICpceRmList] = NaN
        u2[qDICppeRmList] = NaN
        u2[RemoveOutliersList] = NaN
        v2 = v
        v2[qDICpceRmList] = NaN
        v2[qDICppeRmList] = NaN
        v2[RemoveOutliersList] = NaN
        u2 = inpaint_nans(u2,4)
        v2 = inpaint_nans(v2,4)
        # --------------------------------------
        close_('all')
        figure
        surf(u2)
        colorbar
        plt.title('Displacement u','fontweight','normal')
        figure
        surf(v2)
        colorbar
        plt.title('Displacement v','fontweight','normal')
        if len(Thr0)==0 or (Thr0 == 0):
            print('Do you want to redo Median test: 0(Yes, redo it!); 1(No, it is good!)  \n' % ())
            prompt = 'Input here: '
            MedFilterOrNot = input_(prompt)
        else:
            MedFilterOrNot = 1

    
    u = u2
    v = v2
    ## ==============================================
# Manually remove bad points.
    if np.abs(qDICOrNot) > 0:
        print('Do you clear bad points by setting upper/lower bounds once more? (0-yes; 1-no)  \n' % ())
        prompt = 'Input here: '
        ClearBadInitialPointsOrNot = input_(prompt)
        while ClearBadInitialPointsOrNot == 0:

            prompt = 'What is your upper bound for x-displacement? Input: '
            upperbound = input_(prompt)
            row1,col1 = find(u > upperbound)
            prompt = 'What is your lower bound for x-displacement? Input: '
            lowerbound = input_(prompt)
            row2,col2 = find(u < lowerbound)
            prompt = 'What is your upper bound for y-displacement? Input: '
            upperbound = input_(prompt)
            row3,col3 = find(v > upperbound)
            prompt = 'What is your lower bound for y-displacement? Input: '
            lowerbound = input_(prompt)
            row4,col4 = find(v < lowerbound)
            row = np.array([[row1],[row2],[row3],[row4]])
            col = np.array([[col1],[col2],[col3],[col4]])
            for tempi in np.arange(1,len(row)+1).reshape(-1):
                u[row[tempi],col[tempi]] = NaN
                v[row[tempi],col[tempi]] = NaN
                #f11(row(tempi),col(tempi))=NaN; f21(row(tempi),col(tempi))=NaN;
#f12(row(tempi),col(tempi))=NaN; f22(row(tempi),col(tempi))=NaN;
            u = inpaint_nans(u,4)
            v = inpaint_nans(v,4)
            #f11 = inpaint_nans(f11,4); f21 = inpaint_nans(f21,4);
#f12 = inpaint_nans(f12,4); f22 = inpaint_nans(f22,4);
            # --------------------------------------
            close_('all')
            figure
            surf(u)
            colorbar
            plt.title('Displacement u','fontweight','normal')
            figure
            surf(v)
            colorbar
            plt.title('Displacement v','fontweight','normal')
            prompt = 'Do you want to reset upper/lower bounds? (0-yes; 1-no) Input: '
            ClearBadInitialPointsOrNot = input_(prompt)

        ## =========
        print('Do you clear bad points by directly pointing x-disp bad points? (0-yes; 1-no)  \n' % ())
        prompt = 'Input here: '
        ClearBadInitialPointsOrNot = input_(prompt)
        while ClearBadInitialPointsOrNot == 0:

            ########################################
# Have a look at integer search
# --------------------------------------
            close_('all')
            figure
            surf(u)
            colorbar
            view(2)
            plt.title('Displacement u','fontweight','normal')
            # figure; surf(v); colorbar; view(2)
# title('Displacement v','fontweight','normal')
########################################
            row1,col1 = ginput
            row = int(np.floor(col1))
            col = int(np.floor(row1))
            row = row
            col = col
            BadptRow = np.array([[BadptRow],[row]])
            BadptCol = np.array([[BadptCol],[col]])
            row = BadptRow
            col = BadptCol
            for tempi in np.arange(1,len(row)+1).reshape(-1):
                u[row[tempi],col[tempi]] = NaN
                v[row[tempi],col[tempi]] = NaN
                #f11(row(tempi),col(tempi))=NaN; f21(row(tempi),col(tempi))=NaN;
#f12(row(tempi),col(tempi))=NaN; f22(row(tempi),col(tempi))=NaN;
            u = inpaint_nans(u,4)
            v = inpaint_nans(v,4)
            #f11 = inpaint_nans(f11,4); f21 = inpaint_nans(f21,4);
#f12 = inpaint_nans(f12,4); f22 = inpaint_nans(f22,4);
            # --------------------------------------
            close_('all')
            figure
            surf(u)
            colorbar
            plt.title('Displacement u','fontweight','normal')
            # figure; surf(v); colorbar; title('Displacement v','fontweight','normal');
            prompt = 'Do you point out more x-disp bad points? (0-yes; 1-no) Input: '
            ClearBadInitialPointsOrNot = input_(prompt)

        #prompt = 'Do you clear bad points by directly pointing y-disp bad points? (0-yes; 1-no)';
        print('Do you clear bad points by directly pointing y-disp bad points? (0-yes; 1-no)  \n' % ())
        prompt = 'Input here: '
        ClearBadInitialPointsOrNot = input_(prompt)
        while ClearBadInitialPointsOrNot == 0:

            ########################################
# Have a look at integer search
# --------------------------------------
            close_('all')
            # figure; surf(u); colorbar; view(2)
# title('Displacement u','fontweight','normal')
            figure
            surf(v)
            colorbar
            view(2)
            plt.title('Displacement v','fontweight','normal')
            ########################################
            row1,col1 = ginput
            row = int(np.floor(col1))
            col = int(np.floor(row1))
            row = row
            col = col
            BadptRow = np.array([[BadptRow],[row]])
            BadptCol = np.array([[BadptCol],[col]])
            row = BadptRow
            col = BadptCol
            for tempi in np.arange(1,len(row)+1).reshape(-1):
                u[row[tempi],col[tempi]] = NaN
                v[row[tempi],col[tempi]] = NaN
                #f11(row(tempi),col(tempi))=NaN; f21(row(tempi),col(tempi))=NaN;
#f12(row(tempi),col(tempi))=NaN; f22(row(tempi),col(tempi))=NaN;
            u = inpaint_nans(u,4)
            v = inpaint_nans(v,4)
            #f11 = inpaint_nans(f11,4); f21 = inpaint_nans(f21,4);
#f12 = inpaint_nans(f12,4); f22 = inpaint_nans(f22,4);
            # --------------------------------------
            close_('all')
            # figure; surf(u); colorbar; title('Displacement u','fontweight','normal');
            figure
            surf(v)
            colorbar
            plt.title('Displacement v','fontweight','normal')
            # prompt = 'Do you clear bad points by directly pointing y-disp bad points more? (0-yes; 1-no)';
            prompt = 'Do you point out more y-disp bad points? (0-yes; 1-no) Input: '
            ClearBadInitialPointsOrNot = input_(prompt)

        # prompt = 'Do you clear bad points by directly pointing F11 bad points? (0-yes; 1-no)';
# ClearBadInitialPointsOrNot = input(prompt);
        # while ClearBadInitialPointsOrNot == 0
        #     ########################################
#     # Have a look at integer search
#     # --------------------------------------
#     close all;
#     figure; surf(f11); colorbar; view(2)
#     title('Displacement u','fontweight','normal')
#     # figure; surf(v); colorbar; view(2)
#     # title('Displacement v','fontweight','normal')
#     ########################################
        #     [row1, col1] = ginput;
#     row = floor(col1); col = floor(row1);
        #     for tempi = 1:length(row)
#         u(row(tempi),col(tempi))=NaN; v(row(tempi),col(tempi))=NaN;
#         f11(row(tempi),col(tempi))=NaN; f21(row(tempi),col(tempi))=NaN;
#         f12(row(tempi),col(tempi))=NaN; f22(row(tempi),col(tempi))=NaN;
#     end
#     u = inpaint_nans(u,4); v = inpaint_nans(v,4);
#     f11 = inpaint_nans(f11,4); f21 = inpaint_nans(f21,4);
#     f12 = inpaint_nans(f12,4); f22 = inpaint_nans(f22,4);
        #     # --------------------------------------
#     close all;
#     figure; surf(f11); colorbar; title('Displacement u','fontweight','normal');
#     # figure; surf(v); colorbar; title('Displacement v','fontweight','normal');
        #     prompt = 'Do you clear bad points by directly pointing x-disp bad points more? (0-yes; 1-no)';
#     ClearBadInitialPointsOrNot = input(prompt);
        # end
        # prompt = 'Do you clear bad points by directly pointing F21 bad points? (0-yes; 1-no)';
# ClearBadInitialPointsOrNot = input(prompt);
        # while ClearBadInitialPointsOrNot == 0
        #     ########################################
#     # Have a look at integer search
#     # --------------------------------------
#     close all;
#     figure; surf(f21); colorbar; view(2)
#     title('Displacement u','fontweight','normal')
#     # figure; surf(v); colorbar; view(2)
#     # title('Displacement v','fontweight','normal')
#     ########################################
        #     [row1, col1] = ginput;
#     row = floor(col1); col = floor(row1);
        #     for tempi = 1:length(row)
#         u(row(tempi),col(tempi))=NaN; v(row(tempi),col(tempi))=NaN;
#         f11(row(tempi),col(tempi))=NaN; f21(row(tempi),col(tempi))=NaN;
#         f12(row(tempi),col(tempi))=NaN; f22(row(tempi),col(tempi))=NaN;
#     end
#     u = inpaint_nans(u,4); v = inpaint_nans(v,4);
#     f11 = inpaint_nans(f11,4); f21 = inpaint_nans(f21,4);
#     f12 = inpaint_nans(f12,4); f22 = inpaint_nans(f22,4);
        #     # --------------------------------------
#     close all;
#     figure; surf(f21); colorbar; title('Displacement u','fontweight','normal');
#     # figure; surf(v); colorbar; title('Displacement v','fontweight','normal');
        #     prompt = 'Do you clear bad points by directly pointing x-disp bad points more? (0-yes; 1-no)';
#     ClearBadInitialPointsOrNot = input(prompt);
        # end
        # prompt = 'Do you clear bad points by directly pointing F12 bad points? (0-yes; 1-no)';
# ClearBadInitialPointsOrNot = input(prompt);
        # while ClearBadInitialPointsOrNot == 0
        #     ########################################
#     # Have a look at integer search
#     # --------------------------------------
#     close all;
#     figure; surf(f12); colorbar; view(2)
#     title('Displacement u','fontweight','normal')
#     # figure; surf(v); colorbar; view(2)
#     # title('Displacement v','fontweight','normal')
#     ########################################
        #     [row1, col1] = ginput;
#     row = floor(col1); col = floor(row1);
        #     for tempi = 1:length(row)
#         u(row(tempi),col(tempi))=NaN; v(row(tempi),col(tempi))=NaN;
#         f11(row(tempi),col(tempi))=NaN; f21(row(tempi),col(tempi))=NaN;
#         f12(row(tempi),col(tempi))=NaN; f22(row(tempi),col(tempi))=NaN;
#     end
#     u = inpaint_nans(u,4); v = inpaint_nans(v,4);
#     f11 = inpaint_nans(f11,4); f21 = inpaint_nans(f21,4);
#     f12 = inpaint_nans(f12,4); f22 = inpaint_nans(f22,4);
        #     # --------------------------------------
#     close all;
#     figure; surf(f12); colorbar; title('Displacement u','fontweight','normal');
#     # figure; surf(v); colorbar; title('Displacement v','fontweight','normal');
        #     prompt = 'Do you clear bad points by directly pointing x-disp bad points more? (0-yes; 1-no)';
#     ClearBadInitialPointsOrNot = input(prompt);
        # end
        # prompt = 'Do you clear bad points by directly pointing F22 bad points? (0-yes; 1-no)';
# ClearBadInitialPointsOrNot = input(prompt);
        # while ClearBadInitialPointsOrNot == 0
        #     ########################################
#     # Have a look at integer search
#     # --------------------------------------
#     close all;
#     figure; surf(f22); colorbar; view(2)
#     title('Displacement u','fontweight','normal')
#     # figure; surf(v); colorbar; view(2)
#     # title('Displacement v','fontweight','normal')
#     ########################################
        #     [row1, col1] = ginput;
#     row = floor(col1); col = floor(row1);
        #     for tempi = 1:length(row)
#         u(row(tempi),col(tempi))=NaN; v(row(tempi),col(tempi))=NaN;
#         f11(row(tempi),col(tempi))=NaN; f21(row(tempi),col(tempi))=NaN;
#         f12(row(tempi),col(tempi))=NaN; f22(row(tempi),col(tempi))=NaN;
#     end
#     u = inpaint_nans(u,4); v = inpaint_nans(v,4);
#     f11 = inpaint_nans(f11,4); f21 = inpaint_nans(f21,4);
#     f12 = inpaint_nans(f12,4); f22 = inpaint_nans(f22,4);
        #     # --------------------------------------
#     close all;
#     figure; surf(f22); colorbar; title('Displacement u','fontweight','normal');
#     # figure; surf(v); colorbar; title('Displacement v','fontweight','normal');
        #     prompt = 'Do you clear bad points by directly pointing x-disp bad points more? (0-yes; 1-no)';
#     ClearBadInitialPointsOrNot = input(prompt);
    
    return u,v,cc,BadptRow,BadptCol,RemoveOutliersList
    
    ## ========================================================================
    
def funMedRemoveOutliers(u = None,epsilon = None): 
    inpaint_opt = 3
    nSize = 2 * np.array([1,1])
    skipIdx = np.ceil(np.prod(nSize) / 2)
    padOption = 'symmetric'
    u = inpaint_nans(double(u),inpaint_opt)
    medianU = medFilt2(u,nSize,padOption,skipIdx)
    fluct = u - medianU
    medianRes = medFilt2(np.abs(fluct),nSize,padOption,skipIdx)
    normFluct = np.abs(fluct / (medianRes + epsilon))
    return medianU,normFluct
    
    ## ========================================================================
    
def medFilt2(V0 = None,nSize = None,padoption = None,skipIdx = None): 
    # fast median filter for 2D data with extra options.
    
    if len(varargin) < 4:
        skipIdx = 0
    
    if len(varargin) < 3:
        padoption = 'symmetric'
    
    if len(varargin) < 2:
        nSize = np.array([2,2])
    
    nLength = np.prod(nSize)
    if np.mod(nLength,2) == 1:
        padSize = int(np.floor(nSize / 2))
    else:
        if np.mod(nLength,2) == 0:
            padSize = np.array([nSize(1) / 2 - 1,nSize(2) / 2])
    
    if strcmpi(padoption,'none'):
        V = V0
    else:
        V = (padarray(V0,padSize(1) * np.array([1,1]),padoption,'pre'))
        V = (padarray(V,padSize(2) * np.array([1,1]),padoption,'post'))
    
    S = V.shape
    nLength = np.prod(nSize) - sum(skipIdx > 1)
    Vn = single(np.zeros((S(1) - (nSize(1) - 1),S(2) - (nSize(2) - 1),nLength)))
    
    ##
# build the neighborhood
    
    i = cell(1,nSize(1))
    j = cell(1,nSize(2))
    for m in np.arange(1,nSize(1)+1).reshape(-1):
        i[m] = np.arange(m,(S(1) - (nSize(1) - m))+1)
    
    for m in np.arange(1,nSize(2)+1).reshape(-1):
        j[m] = np.arange(m,(S(2) - (nSize(2) - m))+1)
    
    p = 1
    for m in np.arange(1,nSize(1)+1).reshape(-1):
        for n in np.arange(1,nSize(2)+1).reshape(-1):
            if p != skipIdx or skipIdx == 0:
                Vn[:,:,p] = V(i[m],j[n])
            p = p + 1
    
    if skipIdx != 0:
        Vn[:,:,skipIdx] = []
    
    # perform the processing
    Vn = __builtint__.sorted(Vn,3)
    if np.mod(nLength,2) == 1:
        Vr = Vn(:,:,np.ceil(nLength / 2))
    else:
        Vr = mean(cat(3,Vn(:,:,nLength / 2),Vn(:,:,nLength / 2 + 1)),4)
    
    return Vr
    
    ## ========================================================================
    
def removeBadCorrelations(cc = None,ccThreshold = None,sizeChange = None,mSize = None): 
    if sizeChange == 1:
        #recompute threshold, only use pce & ppe since these give the best
#results emprically.
        for ii in np.arange(1,2+1).reshape(-1):
            qf_para[ii],single_distro = bimodal_gauss_fit(cc.qfactors(:,ii))
            if single_distro == 0:
                cc.q_thresh[ii] = qf_para[ii](3) - ccThreshold * qf_para[ii](5)
            else:
                if single_distro == 1:
                    cc.q_thresh[ii] = qf_para[ii](3) - ccThreshold * qf_para[ii](5)
                else:
                    cc.q_thresh[ii] = qf_para[ii](3) - ccThreshold * qf_para[ii](5)
        q_trim = np.array([[cc.q_thresh[0]],[cc.q_thresh[2]]])
    else:
        q_trim = np.array([[cc.q_thresh[0]],[cc.q_thresh[2]]])
    
    #NaN the qfactor values that are below the threshold
    temp = bsxfun(le,cc.qfactors(:,np.arange(1,2+1)),np.transpose(q_trim))
    qfactors_accept = cc.qfactors(:,np.arange(1,2+1))
    qfactors_accept[temp] = NaN
    for ii in np.arange(1,2+1).reshape(-1):
        cc.qfactors_accept[ii] = np.reshape(double(qfactors_accept(:,ii)), tuple(mSize), order="F")
    
    ccMask = np.ones((qfactors_accept.shape,qfactors_accept.shape)) + np.multiply(np.zeros((qfactors_accept.shape,qfactors_accept.shape)),qfactors_accept)
    return cc,ccMask
    
    ## ========================================================================
    
def bimodal_gauss_fit(x = None): 
    #This function takes a dataset and fits a bimodal Gaussian distro to it.
    
    x = __builtint__.sorted(x)
    #set function for bimodal Gaussian
    pdf_normmixture = lambda x = None,p = None,mu1 = None,mu2 = None,sigma1 = None,sigma2 = None: p * normpdf(x,mu1,sigma1) + (1 - p) * normpdf(x,mu2,sigma2)
    pdf_single = lambda x = None,mu1 = None,sigma1 = None: normpdf(x,mu1,sigma1)
    #starting params, biased mixture toward "good" values,
#centered at quartiles, equal std dev.
    pStart = 0.25
    muStart = quantile(x,np.array([0.1,0.75]))
    sigmaStart[1] = np.sqrt(var(x(np.arange(1,np.round(len(x) / 5)+1))))
    #- 0.25*diff(quantile(x,[0.01 0.25])).^2);
    sigmaStart[2] = np.sqrt(var(x(np.arange(np.ceil(len(x) / 10),np.ceil(3 * len(x) / 4)+1))))
    #... - 0.25*diff(quantile(x,[0.25 0.75])).^2);#1:round(length(x)/2)
    start = np.array([pStart,muStart,sigmaStart])
    #set lower and upper bounds
    lb = np.array([0,- inf,- inf,1e-05,1e-05])
    ub = np.array([1,inf,inf,inf,inf])
    #do the parameter estimation
    options = statset('MaxIter',1800,'MaxFunEvals',3600)
    # options.FunValCheck = 'off';
    try:
        single_distro = 0
        paramEsts = mle(x,'pdf',pdf_normmixture,'start',start,'lower',lb,'upper',ub,'options',options)
        if paramEsts(2) - paramEsts(4) >= paramEsts(3) + paramEsts(5) or paramEsts(2) + paramEsts(4) <= paramEsts(3) - paramEsts(5):
            single_distro = 1
            #     disp('Parameters estimated for single peak Gaussian')
            paramEsts = mle(x,'options',options)
            paramEsts = np.array([0.5,paramEsts(1),paramEsts(1),paramEsts(2),paramEsts(2)])
    finally:
        pass
    
    # # #show the result
# # figure
# # # [~, bins] =
# # histogram(x,100);
# # # bins = -2.5:.5:7.5;
# # # h = bar(bins,histc(x,bins)/(length(x)*0.5),'histc');
# # # histogram(x,100)
# # # h.FaceColor = [0.9 0.9 0.9];
# # xgrid = linspace(1.1*min(x),1.1*max(x),200);
# # pdfgrid = pdf_normmixture(xgrid,paramEsts(1),paramEsts(2),paramEsts(3),...
# #     paramEsts(4),paramEsts(5));
# # hold on
# # plot((paramEsts(3) - 2*paramEsts(5)),pdfgrid,'or')
# # plot((paramEsts(2) + 2*paramEsts(4)),pdfgrid,'*r')
# # plot(xgrid,pdfgrid,'-b')
# # hold off
# # xlabel('x')
# # ylabel('Probability Density')
    
    return paramEsts,single_distro
    
    return u,v,cc,BadptRow,BadptCol,RemoveOutliersList