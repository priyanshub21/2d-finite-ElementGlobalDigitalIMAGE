import numpy as np
import matplotlib.pyplot as plt
    
def IntegerSearch(ImgRef = None,ImgDef = None,file_name = None,DICpara = None): 
    #FUNCTION [DICpara,x0,y0,u,v,cc]= IntegerSearch(fNormalized,gNormalized,file_name,DICpara)
# Objective: To compute an inititial guess of the unknown displacement
# field by maximizing the FFT-based cross correlation
# ----------------------------------------------
#   INPUT: ImgRef       Reference image
#          ImgDef       Deformed image
#          file_name    Loaded DIC raw images file name
#          DICpara      Current DIC parameters
    
    #   OUTPUT: DICpara     Updated DIC parameters
#           x0,y0       DIC subset x- and y- positions
#           u,v         x- and y- displacements
#           cc          Cross correlation information
    
    # ----------------------------------------------
# Author: Jin Yang.
# Contact and support: jyang526@wisc.edu -or- aldicdvc@gmail.com
# Last time updated: 02/2020.
# ==============================================
    
    ## Initialization
    gridxROIRange = DICpara.gridxyROIRange.gridx
    gridyROIRange = DICpara.gridxyROIRange.gridy
    winsize = DICpara.winsize
    winstepsize = DICpara.winstepsize
    InitFFTSearchMethod = funParaInput('InitFFTSearchMethod')
    ## To compute the inititial guess from maximizing the FFT-based cross correlation
    if (InitFFTSearchMethod == 1) or (InitFFTSearchMethod == 2):
        InitialGuessSatisfied = 1
        while InitialGuessSatisfied == 1:

            print('--- The size of initial guess search zone (pixels)? ---  \n' % ())
            print('User should start to try a small integer value, and gradually increase the value of \n' % ())
            print('the search zone size until it is larger than the magnitudes of |disp u| and |disp v|. \n' % ())
            print('User could also input [size_x, size_y] to search in a rectangular zone. \n' % ())
            prompt = 'Input here: '
            tempSizeOfSearchRegion = input_(prompt)
            if len(tempSizeOfSearchRegion) == 1:
                tempSizeOfSearchRegion = tempSizeOfSearchRegion * np.array([1,1])
            if (InitFFTSearchMethod == 1):
                x0,y0,u,v,cc = funIntegerSearch(ImgRef,ImgDef,tempSizeOfSearchRegion,gridxROIRange,gridyROIRange,winsize,winstepsize,0,winstepsize)
            else:
                # Open DIC image, and manually click several local seeds.
                figure
                imshow((imread(file_name[0])))
                row1,col1 = ginput
                row = int(np.floor(col1))
                col = int(np.floor(row1))
                x0,y0,u,v,cc = funIntegerSearch(ImgRef,ImgDef,tempSizeOfSearchRegion,gridxROIRange,gridyROIRange,winsize,winstepsize,1,np.array([row,col]))
            ########################################
# Apply ImgRefMask to make u,v nans if there is a hole
            try:
                x0y0Ind = sub2ind(DICpara.ImgSize,x0,y0)
                temp1 = double(DICpara.ImgRefMask(x0y0Ind))
                temp1[not logical[temp1] ] = nan
                HolePtIndMat = np.reshape(temp1, tuple(x0.shape), order="F")
                u = np.multiply(u,HolePtIndMat)
                v = np.multiply(v,HolePtIndMat)
            finally:
                pass
            # --------------------------------------
            ########################################
# Have a look at the integer search results
# --------------------------------------
            close_('all')
            figure
            surf(u)
            colorbar
            plt.title('Displacement u','fontweight','normal')
            set(gca,'fontSize',18)
            plt.title('$x-$displacement $u$','FontWeight','Normal','Interpreter','latex')
            plt.axis('tight')
            plt.xlabel('$x$ (pixels)','Interpreter','latex')
            plt.ylabel('$y$ (pixels)','Interpreter','latex')
            set(gcf,'color','w')
            a = gca
            a.TickLabelInterpreter = 'latex'
            b = colorbar
            b.TickLabelInterpreter = 'latex'
            box('on')
            colormap('jet')
            figure
            surf(v)
            colorbar
            plt.title('Displacement v','fontweight','normal')
            set(gca,'fontSize',18)
            plt.title('$y-$displacement $v$','FontWeight','Normal','Interpreter','latex')
            plt.axis('tight')
            plt.xlabel('$x$ (pixels)','Interpreter','latex')
            plt.ylabel('$y$ (pixels)','Interpreter','latex')
            set(gcf,'color','w')
            a = gca
            a.TickLabelInterpreter = 'latex'
            b = colorbar
            b.TickLabelInterpreter = 'latex'
            box('on')
            colormap('jet')
            ########################################
            print('--- Are you satisfied with initial guess with current search region? (0-yes; 1-no)? ---  \n' % ())
            prompt = 'Input here: '
            InitialGuessSatisfied = input_(prompt)

        # ======== Find some bad inital guess points ========
        cc.ccThreshold = 1.25
        qDICOrNot = 0.5
        Thr0 = 100
        u,v,cc = funRemoveOutliers(u,v,cc,qDICOrNot,Thr0)
        ##
    else:
        tempSizeOfSearchRegion = 0
        x0,y0,u,v,cc = funIntegerSearchMg(ImgRef,ImgDef,gridxROIRange,gridyROIRange,winsize,winstepsize,winstepsize)
        ########################################
# Apply ImgRefMask to make u,v nans if there is a hole
        try:
            x0y0Ind = sub2ind(DICpara.ImgSize,x0,y0)
            temp1 = double(DICpara.ImgRefMask(x0y0Ind))
            temp1[not logical[temp1] ] = nan
            HolePtIndMat = np.reshape(temp1, tuple(x0.shape), order="F")
            u = np.multiply(u,HolePtIndMat)
            v = np.multiply(v,HolePtIndMat)
        finally:
            pass
        # --------------------------------------
        ########################################
# Plotting initial guess
# --------------------------------------
        close_('all')
        figure
        surf(u)
        colorbar
        plt.title('Displacement u','fontweight','normal')
        set(gca,'fontSize',18)
        plt.title('$x-$displacement $u$','FontWeight','Normal','Interpreter','latex')
        plt.axis('tight')
        plt.xlabel('$x$ (pixels)','Interpreter','latex')
        plt.ylabel('$y$ (pixels)','Interpreter','latex')
        set(gcf,'color','w')
        a = gca
        a.TickLabelInterpreter = 'latex'
        b = colorbar
        b.TickLabelInterpreter = 'latex'
        box('on')
        colormap('jet')
        figure
        surf(v)
        colorbar
        plt.title('Displacement v','fontweight','normal')
        set(gca,'fontSize',18)
        plt.title('$y-$displacement $v$','FontWeight','Normal','Interpreter','latex')
        plt.axis('tight')
        plt.xlabel('$x$ (pixels)','Interpreter','latex')
        plt.ylabel('$y$ (pixels)','Interpreter','latex')
        set(gcf,'color','w')
        a = gca
        a.TickLabelInterpreter = 'latex'
        b = colorbar
        b.TickLabelInterpreter = 'latex'
        box('on')
        colormap('jet')
        ########################################
    
    # ======== Finally, update DICpara ========
    DICpara.InitFFTSearchMethod = InitFFTSearchMethod
    DICpara.SizeOfFFTSearchRegion = tempSizeOfSearchRegion
    return DICpara,x0,y0,u,v,cc