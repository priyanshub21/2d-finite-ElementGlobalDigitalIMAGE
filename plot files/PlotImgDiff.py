# Compute image grayscale value difference (in pixels) between the
# reference and deformed images
#
# Author: Jin Yang, jyang526@wisc.edu; aldicdvc@gmail.com
# Date: 2020.11
##

import numpy as np
import matplotlib.pyplot as plt
    
def PlotImgDiff(x0 = None,y0 = None,u = None,v = None,fNormalized = None,gNormalized = None): 
    prompt = 'Do you want to compute image differences using current solved deformations? (0-yes; 1-no)'
    ComputeImgDiffOrNot = input_(prompt)
    if ComputeImgDiffOrNot == 0:
        x1 = np.arange(x0(1),x0(end())+1,1)
        y1 = np.arange(y0(1),y0(end())+1,1)
        #[y1Grid,x1Grid] = meshgrid(y1,x1);
        uFine = gridfit(x0,y0,np.transpose(u),x1,y1,'interp','bilinear')
        uFine = np.transpose(uFine)
        vFine = gridfit(x0,y0,np.transpose(v),x1,y1,'interp','bilinear')
        vFine = np.transpose(vFine)
        DispMask = np.zeros((len(x1),len(y1),2))
        DispMask[:,:,2] = uFine
        DispMask[:,:,1] = vFine
        f1 = fNormalized(np.arange(x1(1),x1(end())+1,1),np.arange(y1(1),y1(end())+1,1))
        # g1 = gNormalized(x1(1)-3:x1(end)+3,y1(1)-3:y1(end)+3);
        g3 = gNormalized(np.arange(x1(1),x1(end())+1,1),np.arange(y1(1),y1(end())+1,1))
        g2 = imwarp(g3,DispMask,'cubic')
        # g2 = ba_interp2(g1,vFine+y1Grid+3,uFine+x1Grid+3,'cubic');
# pause;
        figure
        surf(flip(np.transpose((f1 - g2)),1),'edgecolor','none')
        view(2)
        colorbar
        set(gca,'fontsize',18)
        colormap(gray)
        plt.axis('equal')
        plt.axis('tight')
        set(gcf,'color','w')
        a = gca
        a.TickLabelInterpreter = 'latex'
        plt.title('Image grayscale difference before DIC','interpreter','latex')
        figure
        surf(flip(np.transpose((f1 - g3)),1),'edgecolor','none')
        view(2)
        colorbar
        plt.axis('tight')
        set(gca,'fontsize',18)
        colormap(gray)
        plt.axis('equal')
        plt.axis('tight')
        set(gcf,'color','w')
        a = gca
        a.TickLabelInterpreter = 'latex'
        plt.title('Image grayscale difference after DIC','interpreter','latex')
        sum(sum((f1 - g2) ** 2))
        sum(sum((f1 - g3) ** 2))
    
    return
    