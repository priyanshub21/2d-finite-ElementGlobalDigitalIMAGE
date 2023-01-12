# ==============================================
# function funImgGradient
# ==============================================

import numpy as np
    
def funImgGradient(fNormalized = None,gNormalized = None,varargin = None): 
    #[imgfNormalizedbc,imggNormalizedbc,imgSize,DfAxis] = funImgGradient(fNormalized,gNormalized,varargin)
    
    imgSize = gNormalized.shape
    print('--- Start to compute image gradients ---')
    if len(varargin) > 2:
        method = varargin[2]
    else:
        method = []
    
    if str(method) == str('Splines_interp'):
        imgfNormalizedbc = Spline2D('bicubic',np.array([np.arange(1,fNormalized.shape[1-1]+1,1)]),np.array([np.arange(1,fNormalized.shape[2-1]+1,1)]),fNormalized)
        imggNormalizedbc = Spline2D('bicubic',np.array([np.arange(1,gNormalized.shape[1-1]+1,1)]),np.array([np.arange(1,gNormalized.shape[2-1]+1,1)]),gNormalized)
        # [XX,YY] = ndgrid([1:1:size(fNormalized,1)],[1:1:size(fNormalized,2)]);
# DfDxNormalizedbc = imgfNormalizedbc.eval_Dx(XX,YY);
# DfDyNormalizedbc = imgfNormalizedbc.eval_Dy(XX,YY);
        DfAxis = np.array([0,fNormalized.shape[1-1] - 1,0,fNormalized.shape[2-1] - 1])
    else:
        imgradientMatrix = np.transpose(np.array([- 1 / 60,3 / 20,- 3 / 4,0,3 / 4,- 3 / 20,1 / 60]))
        DfDxStartx = 4
        DfDxStarty = 4
        DfDxEndx = fNormalized.shape[1-1] - 3
        DfDxEndy = fNormalized.shape[2-1] - 3
        # DfDxStartx = coordinates(1,1)-1; DfDxStarty = coordinates(1,2)-1;
# I = fNormalized( coordinates(1,1)-3:coordinates((M+1)*(N+1),1)+3, coordinates(1,2)-3:coordinates((M+1)*(N+1),2)+3);
        I = fNormalized(np.arange(DfDxStartx - 3,DfDxEndx + 3+1),np.arange(DfDxStarty - 3,DfDxEndy + 3+1))
        DfDxNormalizedtemp = imfilter(I,imgradientMatrix)
        DfDyNormalizedtemp = imfilter(I,np.transpose(imgradientMatrix))
        DfDxNormalized = DfDxNormalizedtemp(np.arange(4,end() - 3+1),np.arange(4,end() - 3+1))
        DfDyNormalized = DfDyNormalizedtemp(np.arange(4,end() - 3+1),np.arange(4,end() - 3+1))
        DfAxis = np.array([DfDxStartx,DfDxEndx,DfDxStarty,DfDxEndy]) - 1
        Df.DfAxis = DfAxis
        Df.imgSize = imgSize
        Df.DfDx = DfDxNormalized
        Df.DfDy = DfDyNormalized
    
    print('--- Computing image gradients done ---')
    return Df