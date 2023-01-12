#-------------------------------------------------------
import numpy as np
    
def findPeak2(f = None,varargin = None): 
    #FINDPEAK Find extremum of matrix.
#   [XPEAK,YPEAK,MAX_F] = FINDPEAK(F,SUBPIXEL) finds the extremum of F,
#   MAX_F, and its location (XPEAK, YPEAK). F is a matrix. MAX_F is the maximum
#   absolute value of F, or an estimate of the extremum if a subpixel
#   extremum is requested.
    
    #   SUBPIXEL is a boolean that controls if FINDPEAK attempts to estimate the
#   extremum location to subpixel precision. If SUBPIXEL is false, FINDPEAK
#   returns the coordinates of the maximum absolute value of F and MAX_F is
#   max(abs(F(:))). If SUBPIXEL is true, FINDPEAK fits a 2nd order
#   polynomial to the 9 points surrounding the maximum absolute value of
#   F. In this case, MAX_F is the absolute value of the polynomial evaluated
#   at its extremum.
    
    #   Note: Even if SUBPIXEL is true, there are some cases that result
#   in FINDPEAK returning the coordinates of the maximum absolute value
#   of F:
#   * When the maximum absolute value of F is on the edge of matrix F.
#   * When the coordinates of the estimated polynomial extremum would fall
#     outside the coordinates of the points used to constrain the estimate.
    
    #   Copyright 1993-2004 The MathWorks, Inc.
#   $Revision $  $Date: 2004/10/20 17:54:47 $
    
    # get absolute peak pixel
# figure;mesh(f)
    f = np.abs(f)
    max_f,imax = np.amax((f))
    ypeak,xpeak = ind2sub(f.shape,imax(1))
    vec = np.array([xpeak,ypeak])
    peak = max_f
    if 1 == len(varargin):
        subpixel = True
    else:
        subpixel = varargin[0]
    
    if not subpixel  or xpeak == 1 or xpeak == f.shape[2-1] or ypeak == 1 or ypeak == f.shape[1-1]:
        return vec,peak
    else:
        # fit a 2nd order polynomial to 9 points
# using 9 pixels centered on irow,jcol
        u = f(np.arange(ypeak - 1,ypeak + 1+1),np.arange(xpeak - 1,xpeak + 1+1))
        u = u
        x = np.transpose(np.array([- 1,- 1,- 1,0,0,0,1,1,1]))
        y = np.transpose(np.array([- 1,0,1,- 1,0,1,- 1,0,1]))
        # u(x,y) = A(1) + A(2)*x + A(3)*y + A(4)*x*y + A(5)*x^2 + A(6)*y^2
        X = np.array([np.ones((9,1)),x,y,np.multiply(x,y),x ** 2,y ** 2])
        # u = X*A
        A = np.linalg.solve(X,u)
        # get absolute maximum, where du/dx = du/dy = 0
        x_offset = (- A(3) * A(4) + 2 * A(6) * A(2)) / (A(4) ** 2 - 4 * A(5) * A(6))
        y_offset = - 1 / (A(4) ** 2 - 4 * A(5) * A(6)) * (A(4) * A(2) - 2 * A(5) * A(3))
        if np.abs(x_offset) > 1 or np.abs(y_offset) > 1:
            # adjusted peak falls outside set of 9 points fit,
            return vec,peak
        # return only one-thousandth of a pixel precision
        x_offset = np.round(1000 * x_offset) / 1000
        y_offset = np.round(1000 * y_offset) / 1000
        xpeak = xpeak + x_offset
        ypeak = ypeak + y_offset
        # Calculate extremum of fitted function
        max_f = np.array([1,x_offset,y_offset,x_offset * y_offset,x_offset ** 2,y_offset ** 2]) * A
        max_f = np.abs(max_f)
    
    vec = np.array([xpeak,ypeak])
    peak = max_f
    return vec,peak