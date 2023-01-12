# ==============================================
# function funNormalizeImg
# ----------------------------------------------
#   This function is to normalize images with the formula:
#       Normalized(Img) = (Img-Mean(Img))/sqrt(Std(Img))
#
#   Normalized image grayscale value is between [-1,1]
#
# Author: Jin Yang.
# Last updated 02/2020.
# ==============================================
import numpy as np
    
def funNormalizeImg(Img = None,gridxyROIRange = None): 
    gridxROIRange = gridxyROIRange.gridx
    gridyROIRange = gridxyROIRange.gridy
    # ======= To garantuee ROI value is within image boundary =======
    if gridxROIRange(1) < 1:
        gridxROIRange[1] = 1
    
    if gridxROIRange(2) > Img[0].shape[1-1]:
        gridxROIRange[2] = Img[0].shape[1-1]
    
    if gridyROIRange(1) < 1:
        gridyROIRange[1] = 1
    
    if gridyROIRange(2) > Img[0].shape[2-1]:
        gridyROIRange[2] = Img[0].shape[2-1]
    
    gridxyROIRange.gridx = gridxROIRange
    gridxyROIRange.gridy = gridyROIRange
    # ============== Normalize or compress images ==============
# ============== Normalize f and g ==============
    ImgNormalized = cell(Img.shape)
    for i in np.arange(1,len(Img)+1).reshape(-1):
        tempf = Img[i](np.arange(gridxROIRange(1),gridxROIRange(2)+1),np.arange(gridyROIRange(1),gridyROIRange(2)+1))
        fAvg = mean(tempf)
        fstd = std(tempf)
        ImgNormalized[i] = (Img[i] - fAvg) / fstd
    
    # # Normalize f and g
# tempf = f(gridxROIRange(1):gridxROIRange(2), gridyROIRange(1):gridyROIRange(2));
# fAvg = mean(tempf(:)); fstd = std(tempf(:));
    
    # tempg = g(gridxROIRange(1):gridxROIRange(2), gridyROIRange(1):gridyROIRange(2));
# gAvg = mean(tempg(:)); gstd = std(tempg(:));
    
    # fNormalized = (f-fAvg)/fstd;
# gNormalized = (g-gAvg)/gstd;
    
    # ========= Don't Normalized images and use original images =========
# fNormalized = f; gNormalized = g;
    
    # ========= Do wavelets compression =========
# [fNormalized] = func_compress_dw2d(fNormalized0,'sym4',0.1);
# [gNormalized] = func_compress_dw2d(gNormalized0,'sym4',0.1);
    return ImgNormalized,gridxyROIRange