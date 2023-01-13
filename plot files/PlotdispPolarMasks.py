import warnings
import numpy as np
import matplotlib.pyplot as plt
    
def PlotdispPolarMasks(U = None,coordinatesFEM = None,elementsFEM = None,CurrentImg = None,DICpara = None): 
    #PLOTDISPQUADTREE: to plot DIC solved displacement components
#   PlotdispQuadtree(U,coordinatesFEMWorld,elementsFEM,CurrentImg,DICpara)
# ----------------------------------------------
    
    #   INPUT: U                    Displacement vector:
#                               U = [Ux_node1, Uy_node1, Ux_node2, Uy_node2, ... , Ux_nodeN, Uy_nodeN]';
#          coordinatesFEM       FE mesh coordinates
#          elementsFEM          FE mesh elements
#          CurrentImg           Current deformed image
#          DICpara              DIC paramters
    
    #   OUTPUT: Plots of x-displacement field and y-displacement field.
    
    # ----------------------------------------------
# Reference
# [1] RegularizeNd. Matlab File Exchange open source.
# https://www.mathworks.com/matlabcentral/fileexchange/61436-regularizend
# [2] Gridfit. Matlab File Exchange open source.
# https://www.mathworks.com/matlabcentral/fileexchange/8998-surface-fitting-using-gridfit
# ----------------------------------------------
# Author: Jin Yang.
# Contact and support: jyang526@wisc.edu -or- aldicdvc@gmail.com
# Last date modified: 2020.12
#############################################################
    
    ## Initialization
    warnings.warn('off')
    scipy.io.loadmat('./plotFiles/colormap_RdYlBu.mat','cMap')
    #############################################################
##### convert pixel unit to the physical world unit #####
    try:
        um2px = DICpara.um2px
    finally:
        pass
    
    #############################################################
    
    OrigDICImgTransparency = DICpara.OrigDICImgTransparency
    
    Image2PlotResults = DICpara.Image2PlotResults
    
    disp_u = U(np.arange(1,end()+2,2))
    disp_v = U(np.arange(2,end()+2,2))
    coordinatesFEMWorldDef = np.array([coordinatesFEM(:,1) + Image2PlotResults * disp_u,coordinatesFEM(:,2) + Image2PlotResults * disp_v])
    ########### JY!!!Mask START ###############
    if Image2PlotResults == 1:
        for tempi in np.arange(1,coordinatesFEMWorldDef.shape[1-1]+1).reshape(-1):
            try:
                if CurrentImgMask(int(np.floor(coordinatesFEMWorldDef(tempi,1) / um2px)),(CurrentImgMask.shape[2-1] + 1 - np.ceil(coordinatesFEMWorldDef(tempi,2) / um2px))) == 0:
                    coordinatesFEMWorldDef[tempi,:] = np.array([nan,nan])
            finally:
                pass
    
    ########### JY!!!Mask START ###############
    
    bubble_y = (251 - 28.18) * um2px
    bubble_x = 194.66 * um2px
    r = np.sqrt((coordinatesFEMWorldDef(:,1) - bubble_x) ** 2 + (coordinatesFEMWorldDef(:,2) - bubble_y) ** 2)
    theta = atan2(- (coordinatesFEMWorldDef(:,2) - bubble_y),coordinatesFEMWorldDef(:,1) - bubble_x)
    # disp_rt = [cos(theta), -sin(theta); sin(theta), cos(theta)] * [disp_u, disp_v]';
# disp_r = disp_rt(1);
# disp_t = disp_rt(2);
    disp_r = np.multiply(np.cos(theta),disp_u) - np.multiply(np.sin(theta),disp_v)
    disp_t = np.multiply(np.sin(theta),disp_u) + np.multiply(np.cos(theta),disp_v)
    ## ######################################################
# ====== 1) dispx r ======
#########################################################
    fig1 = figure
    ax1 = axes
    try:
        h1 = imshow(flipud(imread(CurrentImg)),'InitialMagnification','fit')
    finally:
        pass
    
    plt.axis('equal')
    plt.axis('tight')
    box('on')
    set(gca,'fontSize',18)
    view(2)
    set(gca,'ydir','normal')
    hold('on')
    ax2 = axes
    h2 = show([],elementsFEM(:,np.arange(1,4+1)),coordinatesFEMWorldDef / um2px,disp_r,'NoEdgeColor')
    set(gca,'fontSize',18)
    view(2)
    box('on')
    plt.axis('equal')
    plt.axis('tight')
    alpha(h2,OrigDICImgTransparency)
    colormap(cMap)
    caxis('auto')
    ######################################################
###### TODO: manually modify colormap and caxis ######
    colormap(cMap)
    caxis(np.array([0,20]))
    
    # caxis([-35,35]); # caxis([-0.025,0.025]);
# caxis([-1.3,-0.1]);
######################################################
    
    linkaxes(np.array([ax1,ax2]))
    
    ax2.Visible = 'off'
    ax2.XTick = []
    ax2.YTick = []
    
    colormap(ax1,'gray')
    
    set(np.array([ax1,ax2]),'Position',np.array([0.17,0.11,0.685,0.815]))
    ax1.Visible = 'on'
    ax1.TickLabelInterpreter = 'latex'
    ##### convert pixel unit to the physical world unit #####
    xticklabels(ax1,np.transpose(num2cell(np.round(um2px * ax1.XTick * 100) / 100,len(ax1.XTick))))
    yticklabels(ax1,np.transpose(num2cell(np.round(um2px * ax1.YTick * 100) / 100,len(ax1.YTick))))
    cb2 = colorbar('Position',np.array([0.17 + 0.685 + 0.012,0.11,0.03,0.815]))
    cb2.TickLabelInterpreter = 'latex'
    ## ######################################################
# ====== 2) dispy t ======
#########################################################
    fig1 = figure
    ax1 = axes
    try:
        h1 = imshow(flipud(imread(CurrentImg)),'InitialMagnification','fit')
    finally:
        pass
    
    plt.axis('equal')
    plt.axis('tight')
    box('on')
    set(gca,'fontSize',18)
    view(2)
    set(gca,'ydir','normal')
    hold('on')
    ax2 = axes
    h2 = show([],elementsFEM(:,np.arange(1,4+1)),coordinatesFEMWorldDef / um2px,disp_t,'NoEdgeColor')
    set(gca,'fontSize',18)
    view(2)
    box('on')
    plt.axis('equal')
    plt.axis('tight')
    alpha(h2,OrigDICImgTransparency)
    colormap(cMap)
    caxis('auto')
    ######################################################
###### TODO: manually modify colormap and caxis ######
    colormap(cMap)
    caxis(np.array([- 2,2]))
    
    # caxis([5,12]);
# caxis([0,38]); # caxis([-0.025,0.025]);
######################################################
    
    linkaxes(np.array([ax1,ax2]))
    
    ax2.Visible = 'off'
    ax2.XTick = []
    ax2.YTick = []
    
    colormap(ax1,'gray')
    
    set(np.array([ax1,ax2]),'Position',np.array([0.17,0.11,0.685,0.815]))
    ax1.Visible = 'on'
    ax1.TickLabelInterpreter = 'latex'
    ##### convert pixel unit to the physical world unit #####
    xticklabels(ax1,np.transpose(num2cell(np.round(um2px * ax1.XTick * 100) / 100,len(ax1.XTick))))
    yticklabels(ax1,np.transpose(num2cell(np.round(um2px * ax1.YTick * 100) / 100,len(ax1.YTick))))
    cb2 = colorbar('Position',np.array([0.17 + 0.685 + 0.012,0.11,0.03,0.815]))
    cb2.TickLabelInterpreter = 'latex'