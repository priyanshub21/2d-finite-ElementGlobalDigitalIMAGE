import warnings
import numpy as np
import matplotlib.pyplot as plt
    
def Plotstrain0Quadtree(F = None,coordinatesFEMWorld = None,elementsFEM = None,DICpara = None): 
    #PLOTSTRAIN0QUADTREE: to plot DIC solved strain fields in a quadtree mesh
#   [strain_exx,strain_exy,strain_eyy,strain_principal_max,strain_principal_min, ...
#    strain_maxshear,strain_vonMises] = Plotstrain0Quadtree(F,coordinatesFEM,elementsFEM)
    
    #   INPUT: F                    DIC solved deformation gradient tensor
#          coordinatesFEMWorld  coordinates of each points on the image domain
#          elementsFEM          FE-elements
#          DICpara              chosen DIC parameters
    
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
    ## ######################################################
# ====== 1) Strain exx ======
#########################################################
    figure
    show([],elementsFEM(:,np.arange(1,4+1)),coordinatesFEMWorld,strain_exx,'NoEdgeColor')
    # showQuadtree(elementsFEM(:,1:4),coordinatesFEMWorld,strain_exx);
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
    if np.amax(coordinatesFEMWorld(:,1)) < um2px * 200:
        set(gca,'XTick',[])
    
    if np.amax(coordinatesFEMWorld(:,2)) < um2px * 200:
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
    ###### TODO: manually modify colormap and caxis ######
# colormap(jet);  # D Sample
# caxis([-0.1,0.1]) # foam
# caxis([-0.004,0.004]); # Sample 12
# caxis([0,0.1]);
######################################################
    
    ## ######################################################
# ====== 2) Strain exy ======
#########################################################
    figure
    show([],elementsFEM(:,np.arange(1,4+1)),coordinatesFEMWorld,strain_exy,'NoEdgeColor')
    #showQuadtree(elementsFEM(:,1:4),coordinatesFEMWorld,strain_exy);
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
    if np.amax(coordinatesFEMWorld(:,1)) < um2px * 200:
        set(gca,'XTick',[])
    
    if np.amax(coordinatesFEMWorld(:,2)) < um2px * 200:
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
    ###### TODO: manually modify colormap and caxis ######
# colormap(jet);  # D Sample
# caxis([-0.1,0.1]) # foam
# caxis([-0.004,0.004]); # Sample 12
# caxis([0,0.1]);
######################################################
    
    ## ######################################################
# ====== 3) Strain eyy ======
#########################################################
    figure
    show([],elementsFEM(:,np.arange(1,4+1)),coordinatesFEMWorld,strain_eyy,'NoEdgeColor')
    #showQuadtree(elementsFEM(:,1:4),coordinatesFEMWorld,strain_eyy);
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
    if np.amax(coordinatesFEMWorld(:,1)) < um2px * 200:
        set(gca,'XTick',[])
    
    if np.amax(coordinatesFEMWorld(:,2)) < um2px * 200:
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
    ###### TODO: manually modify colormap and caxis ######
# colormap(jet);  # D Sample
# caxis([-0.1,0.1]) # foam
# caxis([-0.004,0.004]); # Sample 12
# caxis([0,0.1]);
######################################################
    
    ## ######################################################
# ====== 4) Strain e_principal_max ======
#########################################################
    figure
    show([],elementsFEM(:,np.arange(1,4+1)),coordinatesFEMWorld,strain_principal_max,'NoEdgeColor')
    #showQuadtree(elementsFEM(:,1:4),coordinatesFEMWorld,strain_principal_max);
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
    if np.amax(coordinatesFEMWorld(:,1)) < um2px * 200:
        set(gca,'XTick',[])
    
    if np.amax(coordinatesFEMWorld(:,2)) < um2px * 200:
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
    ###### TODO: manually modify colormap and caxis ######
# colormap(jet);  # D Sample
# caxis([-0.1,0.1]) # foam
# caxis([-0.004,0.004]); # Sample 12
# caxis([0,0.1]);
######################################################
    
    ## ######################################################
# ====== 5) Strain e_principal_min ======
#########################################################
    figure
    show([],elementsFEM(:,np.arange(1,4+1)),coordinatesFEMWorld,strain_principal_min,'NoEdgeColor')
    #showQuadtree(elementsFEM(:,1:4),coordinatesFEMWorld,strain_principal_min);
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
    if np.amax(coordinatesFEMWorld(:,1)) < um2px * 200:
        set(gca,'XTick',[])
    
    if np.amax(coordinatesFEMWorld(:,2)) < um2px * 200:
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
    ###### TODO: manually modify colormap and caxis ######
# colormap(jet);  # D Sample
# caxis([-0.1,0.1]) # foam
# caxis([-0.004,0.004]); # Sample 12
# caxis([0,0.1]);
######################################################
    
    ## ######################################################
# ====== 6) Strain e_max_shear ======
#########################################################
    figure
    show([],elementsFEM(:,np.arange(1,4+1)),coordinatesFEMWorld,strain_maxshear,'NoEdgeColor')
    #showQuadtree(elementsFEM(:,1:4),coordinatesFEMWorld,strain_maxshear);
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
    if np.amax(coordinatesFEMWorld(:,1)) < um2px * 200:
        set(gca,'XTick',[])
    
    if np.amax(coordinatesFEMWorld(:,2)) < um2px * 200:
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
    ###### TODO: manually modify colormap and caxis ######
# colormap(jet);  # D Sample
# caxis([-0.1,0.1]) # foam
# caxis([-0.004,0.004]); # Sample 12
# caxis([0,0.1]);
######################################################
    
    ## ######################################################
# ====== 7) von Mises equivalent strain ======
#########################################################
    figure
    show([],elementsFEM(:,np.arange(1,4+1)),coordinatesFEMWorld,strain_vonMises,'NoEdgeColor')
    #showQuadtree(elementsFEM(:,1:4),coordinatesFEMWorld,strain_vonMises);
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
    if np.amax(coordinatesFEMWorld(:,1)) < um2px * 200:
        set(gca,'XTick',[])
    
    if np.amax(coordinatesFEMWorld(:,2)) < um2px * 200:
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
    ###### TODO: manually modify colormap and caxis ######
# colormap(jet);  # D Sample
# caxis([-0.1,0.1]) # foam
# caxis([-0.004,0.004]); # Sample 12
# caxis([0,0.1]);
######################################################
    
    return strain_exx,strain_exy,strain_eyy,strain_principal_max,strain_principal_min,strain_maxshear,strain_vonMises
    