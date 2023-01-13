# ==============================================
# function Plotdisp_show
# ==============================================
import numpy as np
import matplotlib.pyplot as plt
    
def Plotdisp_showTri(U0 = None,coordinatesFEM = None,elementsFEM = None,varargin = None): 
    if 4 == len(varargin):
        edgeColorOrNot = varargin[0]
    else:
        edgeColorOrNot = 'EdgeColor'
    
    figure
    show(elementsFEM,[],coordinatesFEM,U0(np.arange(1,end()+2,2)),edgeColorOrNot)
    plt.title('$x-$displacement $u$','FontWeight','Normal','Interpreter','latex')
    view(2)
    set(gca,'fontsize',18)
    plt.axis('tight')
    plt.axis('equal')
    colorbar
    # view([90 -90])
    plt.xlabel('$x$ (pixels)','Interpreter','latex')
    plt.ylabel('$y$ (pixels)','Interpreter','latex')
    set(gcf,'color','w')
    a = gca
    a.TickLabelInterpreter = 'latex'
    b = colorbar
    b.TickLabelInterpreter = 'latex'
    # colormap jet
    box('on')
    figure
    show(elementsFEM,[],coordinatesFEM,U0(np.arange(2,end()+2,2)),edgeColorOrNot)
    plt.title('$y-$displacement $v$','FontWeight','Normal','Interpreter','latex')
    view(2)
    set(gca,'fontsize',18)
    plt.axis('tight')
    plt.axis('equal')
    colorbar
    plt.xlabel('$x$ (pixels)','Interpreter','latex')
    plt.ylabel('$y$ (pixels)','Interpreter','latex')
    set(gcf,'color','w')
    a = gca
    a.TickLabelInterpreter = 'latex'
    b = colorbar
    b.TickLabelInterpreter = 'latex'
    # colormap jet
    box('on')