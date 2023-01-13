# ==============================================
# function Plotstrain_show
# ==============================================

import numpy as np
import matplotlib.pyplot as plt
    
def Plotstrain_showTri(F = None,coordinatesFEM = None,elementsFEM = None,varargin = None): 
    if 4 == len(varargin):
        edgeColorOrNot = varargin[0]
    else:
        edgeColorOrNot = 'EdgeColor'
    
    figure
    show(elementsFEM,[],coordinatesFEM,F(np.arange(1,end()+4,4)),edgeColorOrNot)
    plt.title('Strain $e_{11}$','fontweight','normal','Interpreter','latex')
    set(gca,'fontsize',18)
    view(2)
    plt.axis('tight')
    plt.axis('equal')
    colorbar
    box('on')
    plt.xlabel('$x$ (pixels)','Interpreter','latex')
    plt.ylabel('$y$ (pixels)','Interpreter','latex')
    set(gcf,'color','w')
    a = gca
    a.TickLabelInterpreter = 'latex'
    b = colorbar
    b.TickLabelInterpreter = 'latex'
    # colormap jet;
    
    figure
    show(elementsFEM,[],coordinatesFEM,F(np.arange(4,end()+4,4)),edgeColorOrNot)
    plt.title('Strain $e_{22}$','fontweight','normal','Interpreter','latex')
    set(gca,'fontsize',18)
    view(2)
    plt.axis('tight')
    plt.axis('equal')
    colorbar
    box('on')
    plt.xlabel('$x$ (pixels)','Interpreter','latex')
    plt.ylabel('$y$ (pixels)','Interpreter','latex')
    set(gcf,'color','w')
    a = gca
    a.TickLabelInterpreter = 'latex'
    b = colorbar
    b.TickLabelInterpreter = 'latex'
    # colormap jet;
    
    figure
    show(elementsFEM,[],coordinatesFEM,0.5 * (F(np.arange(2,end()+4,4)) + F(np.arange(3,end()+4,4))),edgeColorOrNot)
    plt.title('Strain $e_{12}$','fontweight','normal','Interpreter','latex')
    set(gca,'fontsize',18)
    view(2)
    plt.axis('tight')
    plt.axis('equal')
    colorbar
    box('on')
    plt.xlabel('$x$ (pixels)','Interpreter','latex')
    plt.ylabel('$y$ (pixels)','Interpreter','latex')
    set(gcf,'color','w')
    a = gca
    a.TickLabelInterpreter = 'latex'
    b = colorbar
    b.TickLabelInterpreter = 'latex'
    # colormap jet;
    
    F_maxshear = np.sqrt((0.5 * (F(np.arange(1,end()+4,4)) - F(np.arange(4,end()+4,4)))) ** 2 + (0.5 * F(np.arange(2,end()+4,4)) + 0.5 * F(np.arange(3,end()+4,4))) ** 2)
    figure
    show(elementsFEM,[],coordinatesFEM,F_maxshear,edgeColorOrNot)
    plt.title('Max shear','fontweight','normal','Interpreter','latex')
    set(gca,'fontsize',18)
    view(2)
    plt.axis('tight')
    plt.axis('equal')
    colorbar
    box('on')
    plt.xlabel('$x$ (pixels)','Interpreter','latex')
    plt.ylabel('$y$ (pixels)','Interpreter','latex')
    set(gcf,'color','w')
    a = gca
    a.TickLabelInterpreter = 'latex'
    b = colorbar
    b.TickLabelInterpreter = 'latex'
    colormap('jet')