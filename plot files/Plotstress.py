import warnings
import numpy as np
import matplotlib.pyplot as plt
    
def Plotstress(DICpara = None,ResultStrain = None,sizeOfImg = None,CurrentImg = None): 
    #PLOTSTRESS: to compute and plot DIC solved stress fields on the original DIC images
#   [stress_sxx,stress_sxy,stress_syy, stress_principal_max_xyplane, ...
#    stress_principal_min_xyplane, stress_maxshear_xyplane, ...
#    stress_maxshear_xyz3d, stress_vonMises]     = Plotstress(DICpara,ResultStrain,sizeOfImg,CurrentImg)
    
    #   INPUT: DICpara          DIC para in the ALDIC code
#          ResultStrain     DIC computed strain field result
#          SizeOfImg        Size of the DIC raw image
#          CurrentImg       File name of current deformed image
    
    #   OUTPUT: stress_sxx      Cauchy stress xx-compoent
#           stress_sxy      Cauchy stress xy-compoent s_xy = s_yx
#           stress_syy      Cauchy stress yy-compoent
#           stress_principal_max_xyplane    max principal stress on the xy-plane
#           stress_principal_min_xyplane    min principal stress on the xy-plane
#           stress_maxshear_xyplane         max shear stress on the xy-plane
#           stress_maxshear_xyz3d           max shear stress on the xyz-three dimensional space
#           stress_vonMises                 von Mises stress
    
    #   Plots:
#       1) stress sxx
#       2) stress sxy
#       3) stress syy
#       4) max principal stress on the xy-plane
#       5) min principal stress on the xy-plane
#       6) max shear stress on the xy-plane
#       7) max shear stress on the xyz-three dimensional space
#       8) equivalent von Mises stress
    
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
# Last time updated: 2020.12.
#############################################################
    
    ## Initialization
    warnings.warn('off')
    scipy.io.loadmat('./plotFiles/colormap_RdYlBu.mat','cMap')
    #############################################################
##### convert pixel unit to the physical world unit #####
    um2px = DICpara.um2px
    #############################################################
    
    OrigDICImgTransparency = DICpara.OrigDICImgTransparency
    
    Image2PlotResults = DICpara.Image2PlotResults
    
    ## Load computed strain fields
    x2 = ResultStrain.strainxCoord
    y2 = ResultStrain.strainyCoord
    dudx = ResultStrain.dudx
    dvdx = ResultStrain.dvdx
    dudy = ResultStrain.dudy
    dvdy = ResultStrain.dvdy
    strain_exx = dudx
    strain_exy = 0.5 * (dvdx + dudy)
    strain_eyy = dvdy
    ## Load displacement components to deform the reference image
    disp_u = ResultStrain.dispu
    disp_v = ResultStrain.dispv
    ## Compute stress components
    
    # ------ For each material model ------
#########################################
    if DICpara.MaterialModel == 1:
        # Linear elasticity -- Plane stress
# stress_xyz = [sxx  sxy  0;      strain_xyz = [exx  exy  0;
#               sxy  syy  0;                    exy  eyy  0;
#               0    0    0];                   0    0    ezz];
        E = DICpara.MaterialModelPara.YoungsModulus
        nu = DICpara.MaterialModelPara.PoissonsRatio
        stress_sxx = E / (1 - nu ** 2) * (strain_exx + nu * strain_eyy)
        stress_syy = E / (1 - nu ** 2) * (nu * strain_exx + strain_eyy)
        stress_sxy = E / (1 + nu) * strain_exy
        # Principal stress
        stress_maxshear_xyplane = np.sqrt((0.5 * (stress_sxx - stress_syy)) ** 2 + stress_sxy ** 2)
        stress_principal_max_xyplane = 0.5 * (stress_sxx + stress_syy) + stress_maxshear_xyplane
        stress_principal_min_xyplane = 0.5 * (stress_sxx + stress_syy) - stress_maxshear_xyplane
        stress_maxshear_xyz3d = np.amax(stress_maxshear_xyplane,0.5 * np.abs(stress_principal_max_xyplane))
        # von Mises stress
        stress_vonMises = np.sqrt(0.5 * ((stress_principal_max_xyplane - stress_principal_min_xyplane) ** 2 + (stress_principal_max_xyplane) ** 2 + (stress_principal_min_xyplane) ** 2))
        #########################################
    else:
        if DICpara.MaterialModel == 2:
            # Linear elasticity -- Plane strain
# stress_xyz = [sxx  sxy  0;      strain_xyz = [exx  exy  0;
#               sxy  syy  0;                    exy  eyy  0;
#               0    0    szz];                 0    0    0];
            E = DICpara.MaterialModelPara.YoungsModulus
            nu = DICpara.MaterialModelPara.PoissonsRatio
            stress_sxx = E * (1 - nu) / (1 + nu) / (1 - 2 * nu) * (strain_exx + nu / (1 - nu) * strain_eyy)
            stress_syy = E * (1 - nu) / (1 + nu) / (1 - 2 * nu) * (strain_eyy + nu / (1 - nu) * strain_exx)
            stress_sxy = E / (1 + nu) * strain_exy
            stress_szz = nu * (stress_sxx + stress_syy)
            # Principal stress
            stress_maxshear_xyplane = np.sqrt((0.5 * (stress_sxx - stress_syy)) ** 2 + stress_sxy ** 2)
            stress_principal_max_xyplane = 0.5 * (stress_sxx + stress_syy) + stress_maxshear_xyplane
            stress_principal_min_xyplane = 0.5 * (stress_sxx + stress_syy) - stress_maxshear_xyplane
            stress_maxshear_xyz3d = np.reshape(np.amax(np.array([stress_maxshear_xyplane,0.5 * np.abs(stress_principal_max_xyplane - stress_szz),0.5 * np.abs(stress_principal_min_xyplane - stress_szz)]),[],2), tuple(stress_maxshear_xyplane.shape), order="F")
            # von Mises stress
            stress_vonMises = np.sqrt(0.5 * ((stress_principal_max_xyplane - stress_principal_min_xyplane) ** 2 + (stress_principal_max_xyplane - stress_szz) ** 2 + (stress_principal_min_xyplane - stress_szz) ** 2))
            #########################################
        else:
            if DICpara.MaterialModel == 3:
                print('User needs to modify the code by yourself.')
                pause
    
    ## #######################################################
# ====== 1) Stress sxx ======
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
    h2 = surf((x2 + Image2PlotResults * disp_u) / um2px,sizeOfImg(2) + 1 - (y2 - Image2PlotResults * disp_v) / um2px,stress_sxx,'EdgeColor','none','LineStyle','none')
    set(gca,'fontSize',18)
    view(2)
    box('on')
    caxis('auto')
    
    alpha(h2,OrigDICImgTransparency)
    plt.axis('equal')
    plt.axis('tight')
    colormap(cMap)
    ######################################################
###### TODO: manually modify colormap and caxis ######
# colormap(jet); caxis([-0.025,0.025]);
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
# ====== 2) Strain sxy ======
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
    h2 = surf((x2 + Image2PlotResults * disp_u) / um2px,sizeOfImg(2) + 1 - (y2 - Image2PlotResults * disp_v) / um2px,stress_sxy,'EdgeColor','none','LineStyle','none')
    set(gca,'fontSize',18)
    view(2)
    box('on')
    caxis('auto')
    
    alpha(h2,OrigDICImgTransparency)
    plt.axis('equal')
    plt.axis('tight')
    colormap(cMap)
    ######################################################
###### TODO: manually modify colormap and caxis ######
# colormap(jet); caxis([-0.025,0.025]);
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
# ====== 3) Strain syy ======
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
    h2 = surf((x2 + Image2PlotResults * disp_u) / um2px,sizeOfImg(2) + 1 - (y2 - Image2PlotResults * disp_v) / um2px,stress_syy,'EdgeColor','none','LineStyle','none')
    set(gca,'fontSize',18)
    view(2)
    box('on')
    caxis('auto')
    
    alpha(h2,OrigDICImgTransparency)
    plt.axis('equal')
    plt.axis('tight')
    colormap(cMap)
    ######################################################
###### TODO: manually modify colormap and caxis ######
# colormap(jet); caxis([-0.015,0.015]);
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
# ====== 4) Strain stress_principal_max_xyplane ======
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
    h2 = surf((x2 + Image2PlotResults * disp_u) / um2px,sizeOfImg(2) + 1 - (y2 - Image2PlotResults * disp_v) / um2px,stress_principal_max_xyplane,'EdgeColor','none','LineStyle','none')
    set(gca,'fontSize',18)
    view(2)
    box('on')
    caxis('auto')
    
    alpha(h2,OrigDICImgTransparency)
    plt.axis('equal')
    plt.axis('tight')
    colormap(cMap)
    ######################################################
###### TODO: manually modify colormap and caxis ######
# colormap(jet); # caxis([-0.025,0.025]);
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
# ====== 5) Strain stress_principal_min_xyplane ======
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
    h2 = surf((x2 + Image2PlotResults * disp_u) / um2px,sizeOfImg(2) + 1 - (y2 - Image2PlotResults * disp_v) / um2px,stress_principal_min_xyplane,'EdgeColor','none','LineStyle','none')
    set(gca,'fontSize',18)
    view(2)
    box('on')
    caxis('auto')
    
    alpha(h2,OrigDICImgTransparency)
    plt.axis('equal')
    plt.axis('tight')
    colormap(cMap)
    ######################################################
###### TODO: manually modify colormap and caxis ######
# colormap(jet); # caxis([-0.025,0.025]);
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
# ====== 6) Strain stress_maxshear_xyplane ======
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
    h2 = surf((x2 + Image2PlotResults * disp_u) / um2px,sizeOfImg(2) + 1 - (y2 - Image2PlotResults * disp_v) / um2px,stress_maxshear_xyplane,'EdgeColor','none','LineStyle','none')
    set(gca,'fontSize',18)
    view(2)
    box('on')
    caxis('auto')
    
    alpha(h2,OrigDICImgTransparency)
    plt.axis('equal')
    plt.axis('tight')
    colormap(cMap)
    ######################################################
###### TODO: manually modify colormap and caxis ######
# colormap(jet); # caxis([-0.025,0.025]);
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
# ====== 7) Strain stress_maxshear_xyz3d ======
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
    h2 = surf((x2 + Image2PlotResults * disp_u) / um2px,sizeOfImg(2) + 1 - (y2 - Image2PlotResults * disp_v) / um2px,stress_maxshear_xyz3d,'EdgeColor','none','LineStyle','none')
    set(gca,'fontSize',18)
    view(2)
    box('on')
    caxis('auto')
    
    alpha(h2,OrigDICImgTransparency)
    plt.axis('equal')
    plt.axis('tight')
    colormap(cMap)
    ######################################################
###### TODO: manually modify colormap and caxis ######
# colormap(jet); # caxis([-0.025,0.025]);
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
# ====== 8) von Mises stress ======
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
    h2 = surf((x2 + Image2PlotResults * disp_u) / um2px,sizeOfImg(2) + 1 - (y2 - Image2PlotResults * disp_v) / um2px,stress_vonMises,'EdgeColor','none','LineStyle','none')
    set(gca,'fontSize',18)
    view(2)
    box('on')
    caxis('auto')
    alpha(h2,OrigDICImgTransparency)
    plt.axis('equal')
    plt.axis('tight')
    colormap(cMap)
    ######################################################
###### TODO: manually modify colormap and caxis ######
# colormap(jet); # caxis([-0.025,0.025]);
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
    return stress_sxx,stress_sxy,stress_syy,stress_principal_max_xyplane,stress_principal_min_xyplane,stress_maxshear_xyplane,stress_maxshear_xyz3d,stress_vonMises
    