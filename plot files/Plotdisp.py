import warnings
import numpy as np
import matplotlib.pyplot as plt
    
def Plotdisp(U = None,x = None,y = None,sizeOfImg = None,CurrentImg = None,varargin = None): 
    #FUNCTION Plotdisp(U,x,y,sizeOfImg,CurrentImg,DICpara)
# To plot DIC solved displacement components
# ----------------------------------------------
    
    #   INPUT: U                 Displacement vector: U = [Ux_node1, Uy_node1, Ux_node2, Uy_node2, ... , Ux_nodeN, Uy_nodeN]';
#          x,y               FE mesh x- and y-coordinates
#          sizeOfImg         DIC deformed image size
#          CurrentImg        Current deformed image
#          DICpara           varargin{1}: DIC paramters
    
    #   OUTPUT: Plots of x-displacement field and y-displacement field.
    
    #   TODO: users could change caxis range based on their own choices.
    
    # ----------------------------------------------
# Reference
# [1] RegularizeNd. Matlab File Exchange open source.
# https://www.mathworks.com/matlabcentral/fileexchange/61436-regularizend
# [2] Gridfit. Matlab File Exchange open source.
# https://www.mathworks.com/matlabcentral/fileexchange/8998-surface-fitting-using-gridfit
# ----------------------------------------------
# Author: Jin Yang.
# Contact and support: jyang526@wisc.edu -or- aldicdvc@gmail.com
# Last time updated: 11/2020.
# ==============================================
    
    ## Initialization
    warnings.warn('off')
    scipy.io.loadmat('./plotFiles/colormap_RdYlBu.mat','cMap')
    ##### Parse Inputs ###############################
    DICpara = parseargs(varargin)
    #############################################################
##### convert pixel unit to the physical world unit #####
    try:
        um2px = DICpara.um2px
    finally:
        pass
    
    #############################################################
    
    OrigDICImgTransparency = DICpara.OrigDICImgTransparency
    
    Image2PlotResults = DICpara.Image2PlotResults
    
    M = x.shape[1-1]
    N = x.shape[2-1]
    u = U(np.arange(1,end()+2,2))
    v = U(np.arange(2,end()+2,2))
    u = reshape(u,M,N)
    v = reshape(v,M,N)
    if M < 9:
        x2 = np.transpose(x(:,1))
    else:
        x2 = np.linspace(x(1,1),x(end(),1),4 * (len(x(:,1)) - 1) + 1)
    
    if N < 9:
        y2 = y(1,:)
    else:
        y2 = np.linspace(y(1,1),y(1,end()),4 * (len(y(1,:)) - 1) + 1)
    
    ## Compute displacement components to manipulate the reference image
    disp_u = gridfit(reshape(x,M * N,1),reshape(y,M * N,1),reshape(u,M * N,1),x2,y2)
    disp_v = gridfit(reshape(x,M * N,1),reshape(y,M * N,1),reshape(v,M * N,1),x2,y2)
    # Please don't delete this line, to deal with the image and physical world coordinates
    x2,y2 = ndgrid(x2,y2)
    x2 = np.transpose(x2)
    y2 = np.transpose(y2)
    ## ######################################################
# ====== 1) dispx u ======
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
    h2 = surf((x2 + Image2PlotResults * disp_u) / um2px,(sizeOfImg(2) + 1) - ((y2 - Image2PlotResults * disp_v)) / um2px,disp_u,'EdgeColor','none','LineStyle','none')
    set(gca,'fontSize',18)
    view(2)
    box('on')
    
    alpha(h2,OrigDICImgTransparency)
    plt.axis('equal')
    plt.axis('tight')
    colormap(cMap)
    ######################################################
###### TODO: manually modify colormap and caxis ######
# colormap(jet); caxis([-1.7,.2]); # caxis([-0.025,0.025]);
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
    # xlabel( '$x$ (pixels)','Interpreter','latex'); ylabel('$y$ (pixels)','Interpreter','latex');
# title('$x-$displacement $u$','FontWeight','Normal','Interpreter','latex'); set(gcf,'color','w');
    
    ## ######################################################
# ====== 2) dispx v ======
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
    h2 = surf((x2 + Image2PlotResults * disp_u) / um2px,(sizeOfImg(2) + 1) - ((y2 - Image2PlotResults * disp_v)) / um2px,disp_v,'EdgeColor','none','LineStyle','none')
    set(gca,'fontSize',18)
    view(2)
    box('on')
    
    alpha(h2,OrigDICImgTransparency)
    plt.axis('equal')
    plt.axis('tight')
    colormap(cMap)
    ######################################################
###### TODO: manually modify colormap and caxis ######
# colormap(jet); caxis([5,12]); # caxis([-0.025,0.025]);
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
    # xlabel( '$x$ (pixels)','Interpreter','latex'); ylabel('$y$ (pixels)','Interpreter','latex');
# title('$y-$displacement $v$','FontWeight','Normal','Interpreter','latex'); set(gcf,'color','w');
    
    return
    
    ## #############################################
    
def parseargs(vargin = None): 
    DICpara = []
    try:
        DICpara = vargin[0]
    finally:
        pass
    
    return DICpara
    