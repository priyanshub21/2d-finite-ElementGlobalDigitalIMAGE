import numpy as np
import matplotlib.pyplot as plt
    
def Plotdisp0(U = None,x = None,y = None,sizeOfImg = None): 
    M = x.shape[1-1]
    N = x.shape[2-1]
    u = U(np.arange(1,end()+2,2))
    v = U(np.arange(2,end()+2,2))
    u = reshape(u,M,N)
    v = reshape(v,M,N)
    u = np.transpose(u)
    v = np.transpose(v)
    x = np.transpose(x)
    y = np.transpose(y)
    M = x.shape[1-1]
    N = x.shape[2-1]
    # imagesc([x(1,1) x(end,1)], [y(1,1) y(1,end)], flipud(g)); hold on;
    if M < 9:
        x2 = np.transpose(x(:,1))
    else:
        x2 = np.linspace(x(1,1),x(end(),1),4 * (len(x(:,1)) - 1) + 1)
        x2 = np.transpose(x2)
    
    if N < 9:
        y2 = y(1,:)
    else:
        y2 = np.linspace(y(1,1),y(1,end()),4 * (len(y(1,:)) - 1) + 1)
        y2 = np.transpose(y2)
    
    z = gridfit(reshape(x,M * N,1),reshape(y,M * N,1),reshape(u,M * N,1),x2,y2)
    # figure;
# contourf(x2,y2,z,20,'EdgeColor','none','LineStyle','none');
# axis equal; axis tight; axis off; set(gca,'fontSize',18); view(2); box on;
# set(gcf,'color','w');  colormap jet; set(gca,'ydir','normal');
# title('$x-$displacement $u$','FontWeight','Normal','Interpreter','latex');
# xlabel('$x$ (pixels)','Interpreter','latex'); ylabel('$y$ (pixels)','Interpreter','latex');
    
    # a = gca; a.TickLabelInterpreter = 'latex';
# b = colorbar; b.TickLabelInterpreter = 'latex'; caxisGet = caxis; delete(b);
# if x(M,N) < 200,set(gca,'XTick',[]);end
# if y(M,N) < 200,set(gca,'YTick',[]);end
    
    # ######################################################
    figure
    surf(x2,sizeOfImg(2) + 1 - y2,flip(z,2),'EdgeColor','none','LineStyle','none')
    plt.axis('tight')
    plt.axis('equal')
    colorbar
    set(gcf,'color','w')
    colormap('jet')
    set(gca,'ydir','normal')
    set(gca,'fontSize',18)
    view(2)
    box('on')
    plt.title('$x-$displacement $u$','FontWeight','Normal','Interpreter','latex')
    plt.xlabel('$x$ (pixels)','Interpreter','latex')
    plt.ylabel('$y$ (pixels)','Interpreter','latex')
    if x(M,N) < 200:
        set(gca,'XTick',[])
    
    if y(M,N) < 200:
        set(gca,'YTick',[])
    
    a = gca
    a.TickLabelInterpreter = 'latex'
    b = colorbar
    b.TickLabelInterpreter = 'latex'
    # alpha(c,.5)
    
    z = gridfit(reshape(x,M * N,1),reshape(y,M * N,1),reshape(v,M * N,1),x2,y2)
    # figure;
# contourf(x2,y2,z,20,'EdgeColor','none','LineStyle','none')
# set(gca,'fontSize',18); view(2);
# title('$y-$displacement $v$','FontWeight','Normal','Interpreter','latex');
# axis tight; axis equal; colorbar;
# if x(M,N) < 200
#     set(gca,'XTick',[]);
# end
# if y(M,N) < 200
#     set(gca,'YTick',[]);
# end
# xlabel('$x$ (pixels)','Interpreter','latex'); ylabel('$y$ (pixels)','Interpreter','latex');
# set(gcf,'color','w');
# a = gca; a.TickLabelInterpreter = 'latex';
# b = colorbar; b.TickLabelInterpreter = 'latex';
    
    # colormap jet
# box on
    
    figure
    surf(x2,sizeOfImg(2) + 1 - y2,flip(z,2),'EdgeColor','none','LineStyle','none')
    set(gca,'fontSize',18)
    view(2)
    plt.title('$y-$displacement $v$','FontWeight','Normal','Interpreter','latex')
    plt.axis('tight')
    plt.axis('equal')
    colorbar
    if x(M,N) < 200:
        set(gca,'XTick',[])
    
    if y(M,N) < 200:
        set(gca,'YTick',[])
    
    plt.xlabel('$x$ (pixels)','Interpreter','latex')
    plt.ylabel('$y$ (pixels)','Interpreter','latex')
    set(gcf,'color','w')
    a = gca
    a.TickLabelInterpreter = 'latex'
    b = colorbar
    b.TickLabelInterpreter = 'latex'
    colormap('jet')
    box('on')
    # print -painters -dpng -r600  fig_Sample14_Local_u.png
    
    # ######################################################
    
    # figure; mesh(x,y,reshape( Exact(coordinatesFEM(:,1),2),M,N)); set(gca,'fontSize',20); view(-5, 30);title('Exact x displacement');
# axis(l);
# # print -painters -dpng -r600  fig_Sample14_Local_uExact.png
    
    
    # figure; mesh(x,y, v);  set(gca,'fontSize',20); view(-5, 30);title('Computed y displacement'); axis tight
# l = axis;
# # print -painters -dpng -r600  fig_Sample14_Local_v.png
    
    
    # figure; mesh(x,y,reshape(0* Exact(coordinatesFEM(:,1),2),M,N)); set(gca,'fontSize',20); view(-5, 30);title('Exact y displacement');
# axis(l);
# # print -painters -dpng -r600  fig_Sample14_Local_vExact.png