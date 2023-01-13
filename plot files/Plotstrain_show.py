import warnings
import numpy as np
import matplotlib.pyplot as plt
    
def Plotstrain_show(F = None,coordinatesFEM = None,elementsFEM = None,varargin = None): 
    #FUNCTION Plotstrain_show(F,coordinatesFEM,elementsFEM)
# To plot DIC solved strain components
# ----------------------------------------------
    
    #   INPUT: F                 Deformation gradient tensor:
#                            F = [F11_node1, F21_node1, F12_node1, F22_node1, ... , F11_nodeN, F21_nodeN, F12_nodeN, F22_nodeN]';
#          coordinatesFEM    FE mesh coordinates
#          elementsFEM       FE mesh elements
#          DICpara           chosen DIC parameters
#          EdgeColorOrNot    show edge color or not
    
    #   OUTPUT: Plots of exx, eyy, and exy strain fields.
    
    # ----------------------------------------------
# Author: Jin Yang.
# Contact and support: jyang526@wisc.edu -or- aldicdvc@gmail.com
# Last date modified: 2020.12
#############################################################
    
    ## Initialization
    warnings.warn('off')
    F = full(F)
    ##### Parse Inputs ###############################
    DICpara,EdgeColorOrNot = parseargs(varargin)
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
# ====== 1) strain exx ======
#############################################################
    figure
    show([],elementsFEM(:,np.arange(1,4+1)),coordinatesFEM,F(np.arange(1,end()+4,4)),EdgeColorOrNot)
    plt.title('Strain $e_{11}$','fontweight','normal','Interpreter','latex')
    set(gca,'fontsize',18)
    view(2)
    plt.axis('tight')
    plt.axis('equal')
    box('on')
    set(gcf,'color','w')
    colorbar
    colormap('jet')
    if um2px == 1:
        plt.xlabel('$x$ (pixels)','Interpreter','latex')
        plt.ylabel('$y$ (pixels)','Interpreter','latex')
    else:
        plt.xlabel('$x$','Interpreter','latex')
        plt.ylabel('$y$','Interpreter','latex')
    
    a = gca
    a.TickLabelInterpreter = 'latex'
    b = colorbar
    b.TickLabelInterpreter = 'latex'
    ## ##########################################################
# ====== 2) strain exy ======
#############################################################
    figure
    show([],elementsFEM(:,np.arange(1,4+1)),coordinatesFEM,0.5 * (F(np.arange(2,end()+4,4)) + F(np.arange(3,end()+4,4))),EdgeColorOrNot)
    plt.title('Strain $e_{12}$','fontweight','normal','Interpreter','latex')
    set(gca,'fontsize',18)
    view(2)
    plt.axis('tight')
    plt.axis('equal')
    box('on')
    set(gcf,'color','w')
    colorbar
    colormap('jet')
    if um2px == 1:
        plt.xlabel('$x$ (pixels)','Interpreter','latex')
        plt.ylabel('$y$ (pixels)','Interpreter','latex')
    else:
        plt.xlabel('$x$','Interpreter','latex')
        plt.ylabel('$y$','Interpreter','latex')
    
    a = gca
    a.TickLabelInterpreter = 'latex'
    b = colorbar
    b.TickLabelInterpreter = 'latex'
    ## ##########################################################
# ====== 3) strain eyy ======
#############################################################
    figure
    show([],elementsFEM(:,np.arange(1,4+1)),coordinatesFEM,F(np.arange(4,end()+4,4)),EdgeColorOrNot)
    plt.title('Strain $e_{22}$','fontweight','normal','Interpreter','latex')
    set(gca,'fontsize',18)
    view(2)
    plt.axis('tight')
    plt.axis('equal')
    box('on')
    set(gcf,'color','w')
    colorbar
    colormap('jet')
    if um2px == 1:
        plt.xlabel('$x$ (pixels)','Interpreter','latex')
        plt.ylabel('$y$ (pixels)','Interpreter','latex')
    else:
        plt.xlabel('$x$','Interpreter','latex')
        plt.ylabel('$y$','Interpreter','latex')
    
    a = gca
    a.TickLabelInterpreter = 'latex'
    b = colorbar
    b.TickLabelInterpreter = 'latex'
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
    