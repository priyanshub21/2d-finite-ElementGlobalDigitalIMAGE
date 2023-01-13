import os
import numpy as np
    
def nanconv(a = None,k = None,varargin = None): 
    # NANCONV Convolution in 1D or 2D ignoring NaNs.
#   C = NANCONV(A, K) convolves A and K, correcting for any NaN values
#   in the input vector A. The result is the same size as A (as though you
#   called 'conv' or 'conv2' with the 'same' shape).
    
    #   C = NANCONV(A, K, 'param1', 'param2', ...) specifies one or more of the following:
#     'edge'     - Apply edge correction to the output.
#     'noedge'   - Do not apply edge correction to the output (default).
#     'nanout'   - The result C should have NaNs in the same places as A.
#     'nonanout' - The result C should have ignored NaNs removed (default).
#                  Even with this option, C will have NaN values where the
#                  number of consecutive NaNs is too large to ignore.
#     '2d'       - Treat the input vectors as 2D matrices (default).
#     '1d'       - Treat the input vectors as 1D vectors.
#                  This option only matters if 'a' or 'k' is a row vector,
#                  and the other is a column vector. Otherwise, this
#                  option has no effect.
    
    #   NANCONV works by running 'conv2' either two or three times. The first
#   time is run on the original input signals A and K, except all the
#   NaN values in A are replaced with zeros. The 'same' input argument is
#   used so the output is the same size as A. The second convolution is
#   done between a matrix the same size as A, except with zeros wherever
#   there is a NaN value in A, and ones everywhere else. The output from
#   the first convolution is normalized by the output from the second
#   convolution. This corrects for missing (NaN) values in A, but it has
#   the side effect of correcting for edge effects due to the assumption of
#   zero padding during convolution. When the optional 'noedge' parameter
#   is included, the convolution is run a third time, this time on a matrix
#   of all ones the same size as A. The output from this third convolution
#   is used to restore the edge effects. The 'noedge' parameter is enabled
#   by default so that the output from 'nanconv' is identical to the output
#   from 'conv2' when the input argument A has no NaN values.
    
    # See also conv, conv2
    
    # AUTHOR: Benjamin Kraus (bkraus@bu.edu, ben@benkraus.com)
# Copyright (c) 2013, Benjamin Kraus
# $Id: nanconv.m 4861 2013-05-27 03:16:22Z bkraus $
    
    # Process input arguments
    for arg in np.arange(1,len(varargin) - 2+1).reshape(-1):
        if 'edge' == varargin[arg].lower():
            edge = True
        else:
            if 'noedge' == varargin[arg].lower():
                edge = False
            else:
                if np.array(['same','full','valid']) == varargin[arg].lower():
                    shape = varargin[arg]
                else:
                    if 'nanout' == varargin[arg].lower():
                        nanout = True
                    else:
                        if 'nonanout' == varargin[arg].lower():
                            nanout = False
                        else:
                            if np.array(['2d','is2d']) == varargin[arg].lower():
                                is1D = False
                            else:
                                if np.array(['1d','is1d']) == varargin[arg].lower():
                                    is1D = True
    
    # Apply default options when necessary.
    if (('edge' is not None) != 1):
        edge = False
    
    if (('nanout' is not None) != 1):
        nanout = False
    
    if (('is1D' is not None) != 1):
        is1D = False
    
    if (('shape' is not None) != 1):
        shape = 'same'
    else:
        if (not str(shape) == str('same') ):
            raise Exception(np.array([mfilename,':NotImplemented']),'Shape '%s' not implemented',shape)
    
    # Get the size of 'a' for use later.
    sza = a.shape
    # If 1D, then convert them both to columns.
# This modification only matters if 'a' or 'k' is a row vector, and the
# other is a column vector. Otherwise, this argument has no effect.
    if (is1D):
        if (not isvector(a)  or not isvector(k) ):
            raise Exception('MATLAB:conv:AorBNotVector','A and B must be vectors.')
        a = a
        k = k
    
    # Flat function for comparison.
    o = np.ones((a.shape,a.shape))
    # Flat function with NaNs for comparison.
    on = np.ones((a.shape,a.shape))
    # Find all the NaNs in the input.
    n = np.isnan(a)
    # Replace NaNs with zero, both in 'a' and 'on'.
    a[n] = 0
    on[n] = 0
    # Check that the filter does not have NaNs.
    if (np.any(np.isnan(k))):
        raise Exception(np.array([mfilename,':NaNinFilter']),'Filter (k) contains NaN values.')
    
    # Calculate what a 'flat' function looks like after convolution.
    if (np.any(n) or edge):
        flat = conv2(on,k,shape)
    else:
        flat = o
    
    # The line above will automatically include a correction for edge effects,
# so remove that correction if the user does not want it.
    if (np.any(n) and not edge ):
        flat = flat / conv2(o,k,shape)
    
    # Do the actual convolution
    c = conv2(a,k,shape) / flat
    # If requested, replace output values with NaNs corresponding to input.
    if (nanout):
        c[n] = NaN
    
    # If 1D, convert back to the original shape.
    if (is1D and sza(1) == 1):
        c = np.transpose(c)
    
    return c
    
    return c