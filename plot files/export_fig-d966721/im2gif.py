#IM2GIF Convert a multiframe image to an animated GIF file
#
# Examples:
#   im2gif infile
#   im2gif infile outfile
#   im2gif(A, outfile)
#   im2gif(..., '-nocrop')
#   im2gif(..., '-nodither')
#   im2gif(..., '-ncolors', n)
#   im2gif(..., '-loops', n)
#   im2gif(..., '-delay', n)
#
# This function converts a multiframe image to an animated GIF.
#
# To create an animation from a series of figures, export to a multiframe
# TIFF file using export_fig, then convert to a GIF, as follows:
#
#    for a = 2 .^ (3:6)
#       peaks(a);
#       export_fig test.tif -nocrop -append
#    end
#    im2gif('test.tif', '-delay', 0.5);
#
#IN:
#   infile - string containing the name of the input image.
#   outfile - string containing the name of the output image (must have the
#             .gif extension). Default: infile, with .gif extension.
#   A - HxWxCxN array of input images, stacked along fourth dimension, to
#       be converted to gif.
#   -nocrop - option indicating that the borders of the output are not to
#             be cropped.
#   -nodither - option indicating that dithering is not to be used when
#               converting the image.
#   -ncolors - option pair, the value of which indicates the maximum number
#              of colors the GIF can have. This can also be a quantization
#              tolerance, between 0 and 1. Default/maximum: 256.
#   -loops - option pair, the value of which gives the number of times the
#            animation is to be looped. Default: 65535.
#   -delay - option pair, the value of which gives the time, in seconds,
#            between frames. Default: 1/15.

# Copyright (C) Oliver Woodford 2011

import numpy as np
import os
    
def im2gif(A = None,varargin = None): 
    # Parse the input arguments
    A,options = parse_args(A,varargin[:])
    if options.crop != 0:
        # Crop
        A = crop_borders(A,A(np.ceil(end() / 2),1,:,1))
    
    # Convert to indexed image
    h,w,c,n = A.shape
    A = reshape(permute(A,np.array([1,2,4,3])),h,w * n,c)
    map = unique(reshape(A,h * w * n,c),'rows')
    if map.shape[1-1] > 256:
        dither_str = np.array(['dither','nodither'])
        dither_str = dither_str[1 + (options.dither == 0)]
        if options.ncolors <= 1:
            B,map = rgb2ind(A,options.ncolors,dither_str)
            if map.shape[1-1] > 256:
                B,map = rgb2ind(A,256,dither_str)
        else:
            B,map = rgb2ind(A,np.amin(np.round(options.ncolors),256),dither_str)
    else:
        if np.amax(map) > 1:
            map = double(map) / 255
            A = double(A) / 255
        B = rgb2ind(im2double(A),map)
    
    B = reshape(B,h,w,1,n)
    # Bug fix to rgb2ind
    map[B[1] + 1,:] = im2double(A(1,1,:))
    # Save as a gif
    imwrite(B,map,options.outfile,'LoopCount',np.round(options.loops(1)),'DelayTime',options.delay)
    return
    
    ## Parse the input arguments
    
def parse_args(A = None,varargin = None): 
    # Set the defaults
    options = struct('outfile','','dither',True,'crop',True,'ncolors',256,'loops',65535,'delay',1 / 15)
    # Go through the arguments
    a = 0
    n = np.asarray(varargin).size
    while a < n:

        a = a + 1
        if ischar(varargin[a]) and not len(varargin[a])==0 :
            if varargin[a](1) == '-':
                opt = varargin[a](np.arange(2,end()+1)).lower()
                if 'nocrop' == opt:
                    options.crop = False
                else:
                    if 'nodither' == opt:
                        options.dither = False
                    else:
                        if not isfield(options,opt) :
                            raise Exception('Option %s not recognized',varargin[a])
                        a = a + 1
                        if ischar(varargin[a]) and not ischar(getattr(options,(opt))) :
                            setattr(options,opt,str2double(varargin[a]))
                        else:
                            setattr(options,opt,varargin[a])
            else:
                options.outfile = varargin[a]

    
    if len(options.outfile)==0:
        if not ischar(A) :
            raise Exception('No output filename given.')
        # Generate the output filename from the input filename
        path,outfile = os.path.split(A)[0],os.path.splitext(os.path.split(A)[1])[0],os.path.splitext(os.path.split(A)[1])[1]
        options.outfile = fullfile(path,np.array([outfile,'.gif']))
    
    if ischar(A):
        # Read in the image
        A = imread_rgb(A)
    
    return A,options
    
    ## Read image to uint8 rgb array
    
def imread_rgb(name = None): 
    # Get file info
    info = imfinfo(name)
    # Special case formats
    if 'gif' == info(1).Format.lower():
        A,map = imread(name,'frames','all')
        if not len(map)==0 :
            map = uint8(map * 256 - 0.5)
            A = np.reshape(map(uint32(A) + 1,:), tuple(np.array([A.shape,map.shape[2-1]])), order="F")
            A = permute(A,np.array([1,2,5,4,3]))
    else:
        if np.array(['tif','tiff']) == info(1).Format.lower():
            A = cell(np.asarray(info).size,1)
            for a in np.arange(1,np.asarray(A).size+1).reshape(-1):
                A[a],map = imread(name,'Index',a,'Info',info)
                if not len(map)==0 :
                    map = uint8(map * 256 - 0.5)
                    A[a] = np.reshape(map(uint32(A[a]) + 1,:), tuple(np.array([A.shape,map.shape[2-1]])), order="F")
                if A[a].shape[3-1] == 4:
                    # TIFF in CMYK colourspace - convert to RGB
                    if isfloat(A[a]):
                        A[a] = A[a] * 255
                    else:
                        A[a] = single(A[a])
                    A[a] = 255 - A[a]
                    A[a][:,:,4] = A[a](:,:,4) / 255
                    A[a] = uint8(np.multiply(A(:,:,np.arange(1,3+1)),A[a](:,:,np.array([4,4,4]))))
            A = cat(4,A[:])
        else:
            A,map,alpha = imread(name)
            A = A(:,:,:,1)
            if not len(map)==0 :
                map = uint8(map * 256 - 0.5)
                A = np.reshape(map(uint32(A) + 1,:), tuple(np.array([A.shape,map.shape[2-1]])), order="F")
            else:
                if A.shape[3-1] == 4:
                    # Assume 4th channel is an alpha matte
                    alpha = A(:,:,4)
                    A = A(:,:,np.arange(1,3+1))
    
    return A,alpha
    