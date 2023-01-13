import warnings
import numpy as np
import matplotlib.pyplot as plt
    
def Plotstrain0(U = None,F = None,Rad = None,x0 = None,y0 = None,sizeOfImg = None,DICpara = None): 
    #PLOTSTRAIN0: to compute and plot DIC solved strain fields
#   [x2,y2,dudx,dvdx,dudy,dvdy,strain_exx,strain_exy,strain_eyy, ...
#   strain_principal_max,strain_principal_min,strain_maxshear,strain_vonMises] = Plotstrain0(U,F,Rad,x0,y0,sizeOfImg,DICpara)
    
    #   INPUT: U                DIC solved displacement vector
#          F                DIC solved deformation gradient tensor
#          x,y              x and y coordinates of each points on the image domain
#          SizeOfImg        Size of the DIC raw image
#          DICpara          DIC chosen parameters
    
    #   OUTPUT: x2,y2                   x- and y-coordinates of points whose strain values are computed
#           disp_u,disp_v           Interpolated dispu and dispv at points {x2,y2}
#           dudx,dvdx,dudy,dvdy     E.g., dudx = d(disp_u)/dx at points {x2,y2}
#           strain_exx              strain xx-compoent
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
    
    #   TODO: users could change caxis range based on their own choices.
    
    
    # Author: Jin Yang  (jyang526@wisc.edu)
# Last date modified: 2020.12
#############################################################
    
    ## Initialization
#############################################################
##### convert pixel unit to the physical world unit #####
    um2px = DICpara.um2px
    #############################################################
    warnings.warn('off')
    scipy.io.loadmat('./plotFiles/colormap_RdYlBu.mat','cMap')
    ## Compute strain components
    
    x = x0(np.arange(1 + Rad,end() - Rad+1),np.arange(1 + Rad,end() - Rad+1))
    y = y0(np.arange(1 + Rad,end() - Rad+1),np.arange(1 + Rad,end() - Rad+1))
    M = x.shape[1-1]
    N = x.shape[2-1]
    u_x = F(np.arange(1,end()+4,4))
    v_x = F(np.arange(2,end()+4,4))
    u_y = F(np.arange(3,end()+4,4))
    v_y = F(np.arange(4,end()+4,4))
    u_x = reshape(u_x,M,N)
    v_x = reshape(v_x,M,N)
    u_y = reshape(u_y,M,N)
    v_y = reshape(v_y,M,N)
    u = U(np.arange(1,end()+2,2))
    v = U(np.arange(2,end()+2,2))
    u0 = reshape(u,M + 2 * Rad,N + 2 * Rad)
    v0 = reshape(v,M + 2 * Rad,N + 2 * Rad)
    u = u0(np.arange(1 + Rad,end() - Rad+1),np.arange(1 + Rad,end() - Rad+1))
    v = v0(np.arange(1 + Rad,end() - Rad+1),np.arange(1 + Rad,end() - Rad+1))
    # imagesc([x(1,1) x(end,1)], [y(1,1) y(1,end)], flipud(g)); hold on;
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
    ## Compute strain components
    dudx = gridfit(reshape(x,M * N,1),reshape(y,M * N,1),reshape(u_x,M * N,1),x2,y2)
    dvdx = gridfit(reshape(x,M * N,1),reshape(y,M * N,1),reshape(v_x,M * N,1),x2,y2)
    dudy = gridfit(reshape(x,M * N,1),reshape(y,M * N,1),reshape(u_y,M * N,1),x2,y2)
    dvdy = gridfit(reshape(x,M * N,1),reshape(y,M * N,1),reshape(v_y,M * N,1),x2,y2)
    strain_exx = dudx
    strain_exy = 0.5 * (dvdx + dudy)
    strain_eyy = dvdy
    strain_maxshear = np.sqrt((0.5 * (strain_exx - strain_eyy)) ** 2 + strain_exy ** 2)
    # Principal strain
    strain_principal_max = 0.5 * (strain_exx + strain_eyy) + strain_maxshear
    strain_principal_min = 0.5 * (strain_exx + strain_eyy) - strain_maxshear
    # equivalent von Mises strain
    strain_vonMises = np.sqrt(strain_principal_max ** 2 + strain_principal_min ** 2 - np.multiply(strain_principal_max,strain_principal_min) + 3 * strain_maxshear ** 2)
    # Please don't delete this line, to deal with the image and physical world coordinates
    x2,y2 = ndgrid(x2,y2)
    x2 = np.transpose(x2)
    y2 = np.transpose(y2)
    ## #######################################################
# ====== 1) Strain exx ======
#########################################################
    figure
    surf(x2,um2px * (sizeOfImg(2) + 1) - y2,strain_exx,'EdgeColor','none','LineStyle','none')
    set(gca,'fontSize',18)
    view(2)
    box('on')
    set(gca,'ydir','normal')
    plt.title('Strain $e_{xx}$','FontWeight','Normal','Interpreter','latex')
    plt.axis('tight')
    plt.axis('equal')
    colorbar
    colormap('jet')
    set(gcf,'color','w')
    if x(M,N) < um2px * 200:
        set(gca,'XTick',[])
    
    if y(M,N) < um2px * 200:
        set(gca,'YTick',[])
    
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
    ######################################################
###### TODO: manually modify colormap and caxis ######
# colormap(jet); caxis([-0.0,0.01]);
######################################################
    
    ## #######################################################
# ====== 2) Strain exy ======
#########################################################
    figure
    surf(x2,um2px * (sizeOfImg(2) + 1) - y2,strain_exy,'EdgeColor','none','LineStyle','none')
    set(gca,'fontSize',18)
    view(2)
    box('on')
    set(gca,'ydir','normal')
    plt.title('Strain $e_{xy}$','FontWeight','Normal','Interpreter','latex')
    plt.axis('tight')
    plt.axis('equal')
    colorbar
    colormap('jet')
    set(gcf,'color','w')
    if x(M,N) < um2px * 200:
        set(gca,'XTick',[])
    
    if y(M,N) < um2px * 200:
        set(gca,'YTick',[])
    
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
    ######################################################
###### TODO: manually modify colormap and caxis ######
# colormap(jet); caxis([-0.0,0.01]);
######################################################
    
    ## #######################################################
# ====== 3) Strain eyy ======
#########################################################
    figure
    surf(x2,um2px * (sizeOfImg(2) + 1) - y2,strain_eyy,'EdgeColor','none','LineStyle','none')
    set(gca,'fontSize',18)
    view(2)
    box('on')
    set(gca,'ydir','normal')
    plt.title('Strain $e_{yy}$','FontWeight','Normal','Interpreter','latex')
    plt.axis('tight')
    plt.axis('equal')
    colorbar
    colormap('jet')
    set(gcf,'color','w')
    if x(M,N) < um2px * 200:
        set(gca,'XTick',[])
    
    if y(M,N) < um2px * 200:
        set(gca,'YTick',[])
    
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
    ######################################################
###### TODO: manually modify colormap and caxis ######
# colormap(jet); caxis([-0.0,0.01]);
######################################################
    
    ## #######################################################
# ====== 4) Strain e_principal_max ======
#########################################################
    figure
    surf(x2,um2px * (sizeOfImg(2) + 1) - y2,strain_principal_max,'EdgeColor','none','LineStyle','none')
    set(gca,'fontSize',18)
    view(2)
    box('on')
    set(gca,'ydir','normal')
    plt.title('Principal strain $e_{\max}$','FontWeight','Normal','Interpreter','latex')
    plt.axis('tight')
    plt.axis('equal')
    colorbar
    colormap('jet')
    set(gcf,'color','w')
    if x(M,N) < um2px * 200:
        set(gca,'XTick',[])
    
    if y(M,N) < um2px * 200:
        set(gca,'YTick',[])
    
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
    ######################################################
###### TODO: manually modify colormap and caxis ######
# colormap(jet); caxis([-0.0,0.01]);
######################################################
    
    ## #######################################################
# ====== 5) Strain e_principal_min ======
#########################################################
    figure
    surf(x2,um2px * (sizeOfImg(2) + 1) - y2,strain_principal_min,'EdgeColor','none','LineStyle','none')
    set(gca,'fontSize',18)
    view(2)
    box('on')
    set(gca,'ydir','normal')
    plt.title('Principal strain $e_{\min}$','FontWeight','Normal','Interpreter','latex')
    plt.axis('tight')
    plt.axis('equal')
    colorbar
    colormap('jet')
    set(gcf,'color','w')
    if x(M,N) < um2px * 200:
        set(gca,'XTick',[])
    
    if y(M,N) < um2px * 200:
        set(gca,'YTick',[])
    
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
    ######################################################
###### TODO: manually modify colormap and caxis ######
# colormap(jet); caxis([-0.0,0.01]);
######################################################
    
    ## #######################################################
# ====== 6) Strain e_max_shear ======
#########################################################
    figure
    surf(x2,um2px * (sizeOfImg(2) + 1) - y2,strain_maxshear,'EdgeColor','none','LineStyle','none')
    set(gca,'fontSize',18)
    view(2)
    box('on')
    set(gca,'ydir','normal')
    plt.title('Max shear strain','FontWeight','Normal','Interpreter','latex')
    plt.axis('tight')
    plt.axis('equal')
    colorbar
    colormap('jet')
    set(gcf,'color','w')
    if x(M,N) < um2px * 200:
        set(gca,'XTick',[])
    
    if y(M,N) < um2px * 200:
        set(gca,'YTick',[])
    
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
    ######################################################
###### TODO: manually modify colormap and caxis ######
# colormap(jet); caxis([-0.0,0.01]);
######################################################
    
    ## #########################################################
# ====== 7) von Mises equivalent strain ======
#########################################################
    figure
    surf(x2,um2px * (sizeOfImg(2) + 1) - y2,strain_vonMises,'EdgeColor','none','LineStyle','none')
    set(gca,'fontSize',18)
    view(2)
    box('on')
    set(gca,'ydir','normal')
    plt.title('von Mises equivalent strain','FontWeight','Normal','Interpreter','latex')
    plt.axis('tight')
    plt.axis('equal')
    colorbar
    colormap('jet')
    set(gcf,'color','w')
    if x(M,N) < um2px * 200:
        set(gca,'XTick',[])
    
    if y(M,N) < um2px * 200:
        set(gca,'YTick',[])
    
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
    ######################################################
###### TODO: manually modify colormap and caxis ######
# colormap(jet); caxis([-0.0,0.01]);
######################################################
    
    return x2,y2,disp_u,disp_v,dudx,dvdx,dudy,dvdy,strain_exx,strain_exy,strain_eyy,strain_principal_max,strain_principal_min,strain_maxshear,strain_vonMises
    