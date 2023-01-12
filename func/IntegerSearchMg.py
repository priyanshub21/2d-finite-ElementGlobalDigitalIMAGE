# ==============================================
# function AdjSearchDomain
# ==============================================
import matplotlib.pyplot as plt
    
def IntegerSearchMg(fNormalized = None,gNormalized = None,file_name = None,DICpara = None): 
    gridxROIRange = DICpara.gridxyROIRange.gridx
    gridyROIRange = DICpara.gridxyROIRange.gridy
    winsize = DICpara.winsize
    winstepsize = DICpara.winstepsize
    # Input initial integer search zone size
# InitialGuessSatisfied = 1; # gridxBackup = gridxROIRange; gridyBackup = gridyROIRange;
# tempNoOfInitPt = 0; # 0-whole field search; 1-several seed points to search;
# tempSizeOfSearchRegion = 0;
    x0,y0,u,v,cc = funIntegerSearchMg(fNormalized,gNormalized,gridxROIRange,gridyROIRange,winsize,winstepsize,winstepsize)
    ########################################
# Have a look at integer search
# --------------------------------------
    close_('all')
    figure
    surf(u)
    colorbar
    plt.title('Displacement u','fontweight','normal')
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
    surf(v)
    colorbar
    plt.title('Displacement v','fontweight','normal')
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
    ########################################
    
    # # ======== Find some bad inital guess points ========
# cc.ccThreshold = 1.25; # bad cross-correlation threshold (mean - ccThreshold*stdev for q-factor distribution)
# qDICOrNot = 1; Thr0 = 100; [u,v,cc] = funRemoveOutliers(u,v,cc,qDICOrNot,Thr0);
    
    # #[u,v] = funRemoveOutliers(u,v);
    
    # # ======== Output final search region radius ========
# SizeOfFFTSearchRegion = tempSizeOfSearchRegion;
    
    return x0,y0,u,v,cc