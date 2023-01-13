import numpy as np
import matplotlib.pyplot as plt
    
def Plotuv(U = None,x = None,y = None): 
    #FUNCTION Plotuv(U,x,y )
# To plot DIC solved displacement fields
# ----------------------------------------------
    
    #   INPUT: U                 Displacement vector:
#                            U = [Ux_node1, Uy_node1, Ux_node2, Uy_node2, ... , Ux_nodeN, Uy_nodeN]';
#          x,y               x- and y-coordinates
    
    #   OUTPUT: Plots of displacement fields.
    
    # ----------------------------------------------
# Author: Jin Yang.
# Contact and support: jyang526@wisc.edu -or- aldicdvc@gmail.com
# Last time updated: 11/2020.
# ==============================================
    
    M = x.shape[1-1]
    N = x.shape[2-1]
    u = U(np.arange(1,end()+2,2))
    v = U(np.arange(2,end()+2,2))
    u = reshape(u,M,N)
    v = reshape(v,M,N)
    figure
    surf(x,y,u)
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
    surf(x,y,v)
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
    ## Some old codes
# print -painters -dpng -r600  fig_Sample14_Local_u.png
    
    # figure; mesh(x,y,reshape( Exact(coordinatesFEM(:,1),2),M,N)); set(gca,'fontSize',20); view(-5, 30);title('Exact x displacement');
# axis(l);
# # print -painters -dpng -r600  fig_Sample14_Local_uExact.png
    
    
    # figure; mesh(x,y, v);  set(gca,'fontSize',20); view(-5, 30);title('Computed y displacement'); axis tight
# l = axis;
# # print -painters -dpng -r600  fig_Sample14_Local_v.png
    
    
    # figure; mesh(x,y,reshape(0* Exact(coordinatesFEM(:,1),2),M,N)); set(gca,'fontSize',20); view(-5, 30);title('Exact y displacement');
# axis(l);
# # print -painters -dpng -r600  fig_Sample14_Local_vExact.png
    