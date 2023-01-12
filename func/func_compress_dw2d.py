import numpy as np
    
def func_compress_dw2d(X = None,wname = None,PERCENTAGE = None): 
    # close all; clear all; clc
# X = imread( 'Sample14 Reference.tif');
    X = double(X)
    # FUNC_COMPRESS_DW2D Saved Compression Process.
#   X: matrix of data
#   -----------------
#   XCMP: matrix of compressed data
#   cfsCMP: decomposition vector (see WAVEDEC2)
#   dimCFS: corresponding bookkeeping matrix
    
    #  Auto-generated by Wavelet Toolbox on 29-Jul-2015 15:31:56
    
    ## Analysis parameters.
#---------------------
#wname = 'haar';
    level = np.amax(np.ceil(log2(X.shape[1-1])),np.ceil(log2(X.shape[2-1])))
    ###############################################
# # Compression parameters.
# #------------------------
# # meth = 'bal_sn';
# sorh = 'h';    # Specified soft or hard thresholding
# thrSettings = 28.701519915151209;
# roundFLAG = true;
###############################################
    
    ## Compression using WDENCMP.
#--------------------------
    coefs,sizes = wavedec2(X,level,wname)
    __,sortingIndices = __builtint__.sorted(np.abs(coefs),'descend')
    # PERCENTAGE = 0.1;
    ZerocoefsIndices = setdiff(sortingIndices,(sortingIndices(np.arange(1,np.round(PERCENTAGE * X.shape[1-1] * X.shape[2-1])+1))))
    coefsNew = coefs
    coefsNew[ZerocoefsIndices] = 0
    ###############################################
# plot(linspace(1,1216461,1216461),coefsNew)
# figure; plot(linspace(1,1216461,1216461),coefs)
###############################################
    
    ###############################################
# [XCMP,cfsCMP,dimCFS] = wdencmp('gbl',coefs,sizes, ...
#     wname,level,thrSettings,sorh,1);
# if roundFLAG , XCMP = round(XCMP); end
# if isequal(class(X),'uint8') , XCMP = uint8(XCMP); end
###############################################
    
    ## Reconstruct using WAVEREC2
    X0 = waverec2(coefsNew,sizes,wname)
    #imshow(uint8(round(X0)))
#figure; imshow(uint8(round(X)))
    
    print(num2str(PERCENTAGE))
    (norm(coefsNew) / norm(coefs)) ** 2
    return X0