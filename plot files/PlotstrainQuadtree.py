import warnings
import numpy as np
import matplotlib.pyplot as plt
    
def PlotstrainQuadtree(U = None,F = None,coordinatesFEM = None,elementsFEM = None,CurrentImg = None,DICpara = None): 
    #PLOTSTRAINQUADTREE: to plot DIC solved strain fields on a quadtree mesh
# and overlaid with the original DIC images
#   [strain_exx,strain_exy,strain_eyy,strain_principal_max,strain_principal_min, ...
#    strain_maxshear,strain_vonMises] = PlotstrainQuadtree(U,F,coordinatesFEMWorld,elementsFEM,CurrentImg,DICpara)
# ----------------------------------------------
    
    #   INPUT: F                    DIC solved deformation gradient tensor
#          coordinatesFE        FE mesh coordinates
#          elementsFEM          FE mesh elements
    
    #   OUTPUT: strain_exx              strain xx-compoent
#           strain_exy              strain xy-compoent
#           strain_eyy              strain yy-compoent
#           strain_principal_max    max principal strain on the xy-plane
#           strain_principal_min    min principal strain on the xy-plane
#           strain_maxshear         max shear strain on the xy-plane
#           strain_vonMises         equivalent von Mises strain
    
    #   Plots:
#       1) strain sxx
#       2) strain sxy
#       3) strain syy
#       4) max principal strain on the xy-plane
#       5) min principal strain on the xy-plane
#       6) max shear strain on the xy-plane
#       7) equivalent von Mises strain
    
    # ----------------------------------------------
# Reference
# [1] RegularizeNd. Matlab File Exchange open source.
# https://www.mathworks.com/matlabcentral/fileexchange/61436-regularizend
# [2] Gridfit. Matlab File Exchange open source.
# https://www.mathworks.com/matlabcentral/fileexchange/8998-surface-fitting-using-gridfit
# ----------------------------------------------
# Author: Jin Yang.
# Contact and support: jyang526@wisc.edu -or- aldicdvc@gmail.com
# Last time updated: 2020.12
##############################################################
    
    ## Initialization
    warnings.warn('off')
    scipy.io.loadmat('./plotFiles/colormap_RdYlBu.mat','cMap')
    #############################################################
##### convert pixel unit to the physical world unit #####
    um2px = DICpara.um2px
    #############################################################
    
    OrigDICImgTransparency = DICpara.OrigDICImgTransparency
    
    Image2PlotResults = DICpara.Image2PlotResults
    
    disp_u = U(np.arange(1,end()+2,2))
    disp_v = U(np.arange(2,end()+2,2))
    coordinatesFEMWorldDef = np.array([coordinatesFEM(:,1) + Image2PlotResults * disp_u,coordinatesFEM(:,2) + Image2PlotResults * disp_v])
    ## Compute strain components
    
    u_x = F(np.arange(1,end()+4,4))
    v_x = F(np.arange(2,end()+4,4))
    u_y = F(np.arange(3,end()+4,4))
    v_y = F(np.arange(4,end()+4,4))
    strain_exx = u_x
    strain_exy = 0.5 * (v_x + u_y)
    strain_eyy = v_y
    strain_maxshear = np.sqrt((0.5 * (strain_exx - strain_eyy)) ** 2 + strain_exy ** 2)
    # Principal strain
    strain_principal_max = 0.5 * (strain_exx + strain_eyy) + strain_maxshear
    strain_principal_min = 0.5 * (strain_exx + strain_eyy) - strain_maxshear
    # equivalent von Mises strain
    strain_vonMises = np.sqrt(strain_principal_max ** 2 + strain_principal_min ** 2 - np.multiply(strain_principal_max,strain_principal_min) + 3 * strain_maxshear ** 2)
    ## #######################################################
# ====== 1) Strain exx ======
#########################################################
    fig1 = figure
    ax1 = axes
    try:
        h1 = imshow(flipud(imread(CurrentImg)),'InitialMagnification','fit')
    finally:
        pass
    
    plt.axis('on')
    plt.axis('equal')
    plt.axis('tight')
    box('on')
    set(gca,'fontSize',18)
    view(2)
    set(gca,'ydir','normal')
    hold('on')
    ax2 = axes
    h2 = show([],elementsFEM(:,np.arange(1,4+1)),coordinatesFEMWorldDef / um2px,strain_exx,'NoEdgeColor')
    set(gca,'fontSize',18)
    view(2)
    box('on')
    plt.axis('equal')
    plt.axis('tight')
    alpha(h2,OrigDICImgTransparency)
    colormap(jet)
    caxis('auto')
    ######################################################
###### TODO: manually modify colormap and caxis ######
# colormap(jet); caxis([0,0.5]) # D Sample
# colormap(jet); caxis([-0.1,0.02]) # foam
# colormap(jet); caxis([-0.004,0]); # Sample 12
# colormap(jet(16)); caxis([-8e-3,0]); # Sample 12
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
    ## #######################################################
# ====== 2) Strain exy ======
#########################################################
    fig1 = figure
    ax1 = axes
    try:
        h1 = imshow(flipud(imread(CurrentImg)),'InitialMagnification','fit')
    finally:
        pass
    
    plt.axis('on')
    plt.axis('equal')
    plt.axis('tight')
    box('on')
    set(gca,'fontSize',18)
    view(2)
    set(gca,'ydir','normal')
    hold('on')
    ax2 = axes
    h2 = show([],elementsFEM(:,np.arange(1,4+1)),coordinatesFEMWorldDef / um2px,strain_exy,'NoEdgeColor')
    set(gca,'fontSize',18)
    view(2)
    box('on')
    plt.axis('equal')
    plt.axis('tight')
    alpha(h2,OrigDICImgTransparency)
    colormap(jet)
    caxis('auto')
    ######################################################
###### TODO: manually modify colormap and caxis ######
# colormap(jet); caxis([-0.08,0.08]) # D Sample
# colormap(jet); caxis([-0.06,0.06]) # foam
# colormap(jet); caxis([-0.008,0.008]); # Sample 12
# colormap(jet(16)); caxis([-8e-3,8e-3]); # Sample 12
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
    ## #######################################################
# ====== 3) Strain eyy ======
#########################################################
    fig1 = figure
    ax1 = axes
    try:
        h1 = imshow(flipud(imread(CurrentImg)),'InitialMagnification','fit')
    finally:
        pass
    
    plt.axis('on')
    plt.axis('equal')
    plt.axis('tight')
    box('on')
    set(gca,'fontSize',18)
    view(2)
    set(gca,'ydir','normal')
    hold('on')
    ax2 = axes
    h2 = show([],elementsFEM(:,np.arange(1,4+1)),coordinatesFEMWorldDef / um2px,strain_eyy,'NoEdgeColor')
    set(gca,'fontSize',18)
    view(2)
    box('on')
    plt.axis('equal')
    plt.axis('tight')
    alpha(h2,OrigDICImgTransparency)
    colormap(jet)
    caxis('auto')
    ######################################################
###### TODO: manually modify colormap and caxis ######
# colormap(jet); caxis([-0.15,0]) # D Sample
# colormap(jet); caxis([-0.05,0.2]) # foam
# colormap(jet); caxis([-0.002,0.017]); # Sample 12
# colormap(jet(16)); caxis([0,0.025]); # Sample 12
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
    ## #######################################################
# ====== 4) Strain e_principal_max ======
#########################################################
    fig1 = figure
    ax1 = axes
    try:
        h1 = imshow(flipud(imread(CurrentImg)),'InitialMagnification','fit')
    finally:
        pass
    
    plt.axis('on')
    plt.axis('equal')
    plt.axis('tight')
    box('on')
    set(gca,'fontSize',18)
    view(2)
    set(gca,'ydir','normal')
    hold('on')
    ax2 = axes
    h2 = show([],elementsFEM(:,np.arange(1,4+1)),coordinatesFEMWorldDef / um2px,strain_principal_max,'NoEdgeColor')
    set(gca,'fontSize',18)
    view(2)
    box('on')
    plt.axis('equal')
    plt.axis('tight')
    alpha(h2,OrigDICImgTransparency)
    colormap(jet)
    caxis('auto')
    ######################################################
###### TODO: manually modify colormap and caxis ######
# colormap(jet);  caxis auto; # D Sample
# colormap(jet); caxis auto # foam
# colormap(jet); caxis([0,0.02]); # Sample 12
# colormap(jet); caxis([-0.1,0.1]);
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
    ## #######################################################
# ====== 5) Strain e_principal_min ======
#########################################################
    fig1 = figure
    ax1 = axes
    try:
        h1 = imshow(flipud(imread(CurrentImg)),'InitialMagnification','fit')
    finally:
        pass
    
    plt.axis('on')
    plt.axis('equal')
    plt.axis('tight')
    box('on')
    set(gca,'fontSize',18)
    view(2)
    set(gca,'ydir','normal')
    hold('on')
    ax2 = axes
    h2 = show([],elementsFEM(:,np.arange(1,4+1)),coordinatesFEMWorldDef / um2px,strain_principal_min,'NoEdgeColor')
    set(gca,'fontSize',18)
    view(2)
    box('on')
    plt.axis('equal')
    plt.axis('tight')
    alpha(h2,OrigDICImgTransparency)
    colormap(jet)
    caxis('auto')
    ######################################################
###### TODO: manually modify colormap and caxis ######
# colormap(jet);  caxis auto; # D Sample
# colormap(jet); caxis auto # foam
# colormap(jet); caxis([-0.008,0]); # Sample 12
# colormap(jet); caxis([-0.15,0.0]);
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
    ## #######################################################
# ====== 6) Strain e_max_shear ======
#########################################################
    fig1 = figure
    ax1 = axes
    try:
        h1 = imshow(flipud(imread(CurrentImg)),'InitialMagnification','fit')
    finally:
        pass
    
    plt.axis('on')
    plt.axis('equal')
    plt.axis('tight')
    box('on')
    set(gca,'fontSize',18)
    view(2)
    set(gca,'ydir','normal')
    hold('on')
    ax2 = axes
    h2 = show([],elementsFEM(:,np.arange(1,4+1)),coordinatesFEMWorldDef / um2px,strain_maxshear,'NoEdgeColor')
    set(gca,'fontSize',18)
    view(2)
    box('on')
    plt.axis('equal')
    plt.axis('tight')
    alpha(h2,OrigDICImgTransparency)
    colormap(jet)
    caxis('auto')
    ######################################################
###### TODO: manually modify colormap and caxis ######
# colormap(jet);  caxis auto; # D Sample
# colormap(jet); caxis auto # foam
# colormap(jet); caxis([0,0.011]); # Sample 12
# colormap(jet); caxis([0,0.1]);
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
    ## #######################################################
# ====== 7) von Mises equivalent strain ======
#########################################################
    fig1 = figure
    ax1 = axes
    try:
        h1 = imshow(flipud(imread(CurrentImg)),'InitialMagnification','fit')
    finally:
        pass
    
    plt.axis('on')
    plt.axis('equal')
    plt.axis('tight')
    box('on')
    set(gca,'fontSize',18)
    view(2)
    set(gca,'ydir','normal')
    hold('on')
    ax2 = axes
    h2 = show([],elementsFEM(:,np.arange(1,4+1)),coordinatesFEMWorldDef / um2px,strain_vonMises,'NoEdgeColor')
    set(gca,'fontSize',18)
    view(2)
    box('on')
    plt.axis('equal')
    plt.axis('tight')
    alpha(h2,OrigDICImgTransparency)
    colormap(jet)
    caxis('auto')
    ######################################################
###### TODO: manually modify colormap and caxis ######
# colormap(jet);  caxis auto; # D Sample
# colormap(jet); caxis auto # foam
# colormap(jet); caxis([0,0.025]); # Sample 12
# colormap(jet); caxis([0,0.25]);
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
    return strain_exx,strain_exy,strain_eyy,strain_principal_max,strain_principal_min,strain_maxshear,strain_vonMises
    