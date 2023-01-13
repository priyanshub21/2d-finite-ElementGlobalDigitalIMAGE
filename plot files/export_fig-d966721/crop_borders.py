import numpy as np
import numpy.matlib
    
def crop_borders(A = None,bcol = None,padding = None): 
    #CROP_BORDERS Crop the borders of an image or stack of images
    
    #   [B, vA, vB, bb_rel] = crop_borders(A, bcol, [padding])
    
    #IN:
#   A - HxWxCxN stack of images.
#   bcol - Cx1 background colour vector.
#   padding - scalar indicating how much padding to have in relation to
#             the cropped-image-size (0<=padding<=1). Default: 0
    
    #OUT:
#   B - JxKxCxN cropped stack of images.
#   vA     - coordinates in A that contain the cropped image
#   vB     - coordinates in B where the cropped version of A is placed
#   bb_rel - relative bounding box (used for eps-cropping)
    
    # 06/03/15: Improved image cropping thanks to Oscar Hartogensis
# 08/06/15: Fixed issue #76: case of transparent figure bgcolor
    
    if len(varargin) < 3:
        padding = 0
    
    h,w,c,n = A.shape
    if len(bcol)==0:
        bcol = A(np.ceil(end() / 2),1,:,1)
    
    if np.isscalar(bcol):
        bcol = bcol(np.ones((c,1)))
    
    # Crop margin from left
    bail = False
    for l in np.arange(1,w+1).reshape(-1):
        for a in np.arange(1,c+1).reshape(-1):
            if not np.all(col(A(:,l,a,:)) == bcol(a)) :
                bail = True
                break
        if bail:
            break
    
    # Crop margin from right
    bcol = A(np.ceil(end() / 2),w,:,1)
    bail = False
    for r in np.arange(w,l+- 1,- 1).reshape(-1):
        for a in np.arange(1,c+1).reshape(-1):
            if not np.all(col(A(:,r,a,:)) == bcol(a)) :
                bail = True
                break
        if bail:
            break
    
    # Crop margin from top
    bcol = A(1,np.ceil(end() / 2),:,1)
    bail = False
    for t in np.arange(1,h+1).reshape(-1):
        for a in np.arange(1,c+1).reshape(-1):
            if not np.all(col(A(t,:,a,:)) == bcol(a)) :
                bail = True
                break
        if bail:
            break
    
    # Crop margin from bottom
    bcol = A(h,np.ceil(end() / 2),:,1)
    bail = False
    for b in np.arange(h,t+- 1,- 1).reshape(-1):
        for a in np.arange(1,c+1).reshape(-1):
            if not np.all(col(A(b,:,a,:)) == bcol(a)) :
                bail = True
                break
        if bail:
            break
    
    # Crop the background, leaving one boundary pixel to avoid bleeding on resize
#v = [max(t-padding, 1) min(b+padding, h) max(l-padding, 1) min(r+padding, w)];
#A = A(v(1):v(2),v(3):v(4),:,:);
    if padding == 0:
        padding = 1
    else:
        if np.abs(padding) < 1:
            padding = np.sign(padding) * np.round(mean(np.array([b - t,r - l])) * np.abs(padding))
        else:
            padding = np.round(padding)
    
    if padding > 0:
        # Create an empty image, containing the background color, that has the
# cropped image size plus the padded border
        B = np.matlib.repmat(bcol,(b - t) + 1 + padding * 2,(r - l) + 1 + padding * 2)
        # vA - coordinates in A that contain the cropped image
        vA = np.array([t,b,l,r])
        # vB - coordinates in B where the cropped version of A will be placed
        vB = np.array([padding + 1,(b - t) + 1 + padding,padding + 1,(r - l) + 1 + padding])
        # Place the original image in the empty image
        B[np.arange[vB[1],vB[2]+1],np.arange[vB[3],vB[4]+1],:] = A(np.arange(vA(1),vA(2)+1),np.arange(vA(3),vA(4)+1),:)
        A = B
    else:
        vA = np.array([t - padding,b + padding,l - padding,r + padding])
        A = A(np.arange(vA(1),vA(2)+1),np.arange(vA(3),vA(4)+1),:)
        vB = np.array([NaN,NaN,NaN,NaN])
    
    # For EPS cropping, determine the relative BoundingBox - bb_rel
    bb_rel = np.array([l - 1,h - b - 1,r + 1,h - t + 1]) / np.array([w,h,w,h])
    return A,vA,vB,bb_rel
    
    
def col(A = None): 
    A = A
    return A
    
    return A,vA,vB,bb_rel