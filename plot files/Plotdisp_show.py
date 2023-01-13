import warnings
import numpy as np
import matplotlib.pyplot as plt
    
def Plotdisp_show(U = None,coordinatesFEM = None,elementsFEM = None,varargin = None): 
    #PLOTDISP_SHOW: to plot DIC solved displacement components
#   Plotdisp_show(U,coordinatesFEM,elementsFEM)
# ----------------------------------------------
    
    #   INPUT: U                 Displacement vector:
#                            U = [Ux_node1, Uy_node1, Ux_node2, Uy_node2, ... , Ux_nodeN, Uy_nodeN]';
#          coordinatesFEM    FE mesh coordinates
#          elementsFEM       FE mesh elements
#          DICpara           chosen DIC parameters
#          EdgeColorOrNot    show edge color or not
    
    #   OUTPUT: Plots of x-displacement field and y-displacement field.
    
    #   TODO: users could change caxis range based on their own choices.
    
    # ----------------------------------------------
# Author: Jin Yang.
# Contact and support: jyang526@wisc.edu -or- aldicdvc@gmail.com
# Last date modified: 2020.12
#############################################################
    
    ## Initialization
    warnings.warn('off')
    U = full(U)
    ##### Parse Inputs ###############################
    DICpara,EdgeColorOrNot = parseargs(varargin)
    #############################################################
##### convert pixel unit to the physical world unit #####
    try:
        um2px = DICpara.um2px
    finally:
        pass
    
    #############################################################
    
    try:
        EdgeColorOrNot = EdgeColorOrNot
    finally:
        pass
    
    ## ##########################################################
# ====== 1) dispx u ======
#############################################################
    figure
    show([],elementsFEM(:,np.arange(1,4+1)),coordinatesFEM,U(np.arange(1,end()+2,2)),EdgeColorOrNot)
    plt.title('$x-$displacement $u$','FontWeight','Normal','Interpreter','latex')
    view(2)
    set(gca,'fontsize',18)
    plt.axis('tight')
    plt.axis('equal')
    colorbar
    
    if um2px == 1:
        plt.xlabel('$x$ (pixels)','Interpreter','latex')
        plt.ylabel('$y$ (pixels)','Interpreter','latex')
    else:
        plt.xlabel('$x$','Interpreter','latex')
        plt.ylabel('$y$','Interpreter','latex')
    
    set(gcf,'color','w')
    colormap('jet')
    box('on')
    a = gca
    a.TickLabelInterpreter = 'latex'
    b = colorbar
    b.TickLabelInterpreter = 'latex'
    ###### TODO: manually modify colormap and caxis ######
# colormap(jet);  # D Sample
# caxis([-0.1,0.1]) # foam
# caxis([-0.004,0.004]); # Sample 12
# caxis([0,0.1]);
######################################################
    
    ## ##########################################################
# ====== 2) dispx v ======
#############################################################
    figure
    show([],elementsFEM(:,np.arange(1,4+1)),coordinatesFEM,U(np.arange(2,end()+2,2)),EdgeColorOrNot)
    plt.title('$y-$displacement $v$','FontWeight','Normal','Interpreter','latex')
    view(2)
    set(gca,'fontsize',18)
    plt.axis('tight')
    plt.axis('equal')
    colorbar
    if um2px == 1:
        plt.xlabel('$x$ (pixels)','Interpreter','latex')
        plt.ylabel('$y$ (pixels)','Interpreter','latex')
    else:
        plt.xlabel('$x$','Interpreter','latex')
        plt.ylabel('$y$','Interpreter','latex')
    
    set(gcf,'color','w')
    colormap('jet')
    box('on')
    a = gca
    a.TickLabelInterpreter = 'latex'
    b = colorbar
    b.TickLabelInterpreter = 'latex'
    ###### TODO: manually modify colormap and caxis ######
# colormap(jet);  # D Sample
# caxis([-0.1,0.1]) # foam
# caxis([-0.004,0.004]); # Sample 12
# caxis([0,0.1]);
######################################################
    
    return
    
    ## #############################################
    
def parseargs(vargin = None): 
    DICpara = []
    EdgeColorOrNot = []
    try:
        DICpara = vargin[0]
        try:
            EdgeColorOrNot = vargin[2]
        finally:
            pass
    finally:
        pass
    
    return DICpara,EdgeColorOrNot
    