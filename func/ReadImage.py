import numpy as np
import os
import matplotlib.pyplot as plt
    
def ReadImage(varargin = None): 
    #FUNCTION [file_name,Img,DICpara] = ReadImage(varargin)
# MATLAB script: ReadImage.m
# ----------------------------------------------
#   This script is to load DIC images
#   Images can be loaded by:
#       i) selecting a folder which included all the DIC raw images,
#       ii) inputing image file name prefix keywords
#       iii) manually select DIC raw images
    
    #   INPUT: No inputs are needed
    
    #   OUTPUT: file_name    Loaded DIC raw image file name
#           Img          Loaded DIC images
#           DICpara      DIC parameters
    
    # ----------------------------------------------
# Author: Jin Yang.
# Contact and support: jyang526@wisc.edu -or- aldicdvc@gmail.com
# Last time updated: 02/2020.
# ==============================================
    
    ##
    print('Choose method to load images:  \n' % ())
    print('     0: Select images folder;  \n' % ())
    print('     1: Use prefix of image names;  \n' % ())
    print('     2: Manually select images.  \n' % ())
    prompt = 'Input here: '
    LoadImgMethod = input_(prompt)
    if 0 == LoadImgMethod:
        # ==============================================
        imgfoldername = uigetdir(pwd,'Select images folder')
        addpath(np.array([imgfoldername,'\']))
        img1 = dir(fullfile(imgfoldername,'*.jpg'))
        img2 = dir(fullfile(imgfoldername,'*.jpeg'))
        img3 = dir(fullfile(imgfoldername,'*.tif'))
        img4 = dir(fullfile(imgfoldername,'*.tiff'))
        img5 = dir(fullfile(imgfoldername,'*.bmp'))
        img6 = dir(fullfile(imgfoldername,'*.png'))
        img7 = dir(fullfile(imgfoldername,'*.jp2'))
        file_name = np.array([[img1],[img2],[img3],[img4],[img5],[img6],[img7]])
        file_name = struct2cell(file_name)
    else:
        if 1 == LoadImgMethod:
            # ==============================================
            print('What is prefix of DIC images? E.g. img_0*.tif.   \n' % ())
            prompt = 'Input here: '
            file_name = input_(prompt,'s')
            __,imgname,imgext = os.path.split(file_name)[0],os.path.splitext(os.path.split(file_name)[1])[0],os.path.splitext(os.path.split(file_name)[1])[1]
            file_name = dir(np.array([imgname,imgext]))
            file_name = struct2cell(file_name)
        else:
            # ==============================================
            print('--- Please load first image ---')
            file_name[1,1] = uigetfile('*.tif','Select reference Image (Deformed)')
            print('--- Please load next image ---')
            file_name[1,2] = uigetfile('*.tif','Select deformed Image (Reference)')
            prompt = 'Do you want to load more deformed images? (0-Yes; 1-No)'
            DoYouWantToLoadMoreImages = input_(prompt)
            imageNo = 2
            while (DoYouWantToLoadMoreImages == 0):

                imageNo = imageNo + 1
                file_name[1,imageNo] = uigetfile('*.tif','Select Deformed Image')
                prompt = 'Do you want to load more deformed images? (0-Yes; 1-No)'
                DoYouWantToLoadMoreImages = input_(prompt)

    
    # ==============================================
# The following codes only consider two images comparasion
    numImages = file_name.shape[2-1]
    for i in np.arange(1,numImages+1).reshape(-1):
        Img[i] = imread(file_name[1,i])
        # Change color RGB images to grayscale images
        __,__,numberOfColorChannels = Img[i].shape
        if (numberOfColorChannels == 3):
            Img[i] = rgb2gray(Img[i])
        Img[i] = np.transpose(double(Img[i]))
    
    # ====== COMMENT ======
# Images physical world coordinates and image coordinates are different:
# --------------------
# --  This is image --
# |                  |
# y                  |
# |                  |
# |  --> x direction |
# |                  |
# --------------------
# after transforming,  MatLab matrix direction:
# --  This is matrix in Matlab --
# |                             |
# x                             |
# |                             |
# |  --> y direction            |
# |                             |
# --------------------------------
    
    # ==============================================
# Decide DIC subset parameters
# Choose ZOI
    print('\n' % ())
    print('--- Define ROI corner points at the top-left and the bottom-right ---')
    imshow((imread(file_name[0])))
    plt.title('Click top-left and the bottom-right corner points','fontweight','normal','fontsize',16)
    gridx = np.zeros((1,2))
    gridy = np.zeros((1,2))
    gridx[1],gridy[1] = ginput(1)
    print('Coordinates of top-left corner point are (%4.3f,%4.3f)\n' % (gridx(1),gridy(1)))
    gridx[2],gridy[2] = ginput(1)
    print('Coordinates of bottom-right corner point are (%4.3f,%4.3f)\n' % (gridx(2),gridy(2)))
    gridxy.gridx = np.round(gridx)
    gridxy.gridy = np.round(gridy)
    # Choose finite element size
    print('\n' % ())
    print('--- What is the finite element size (unit: px) ? --- \n' % ())
    prompt = 'Input an integer number: '
    winstepsize = input_(prompt)
    winsize = 2 * np.ceil(winstepsize / 2)
    # ==============================================
# Subproblem 2 solver: finite difference or finite element
    Subpb2FDOrFEM = 2
    
    # Subpb2FDOrFEM = funParaInput('Subpb2FDOrFEM'); # Subproblem 2 using finite difference or fem?
    
    # ==============================================
# Parallel cluster #
# ClusterNo = funParaInput('ClusterNo'); # Assign parpool cluster No
    
    # ==============================================
# Deal with image sequence
    NewFFTSearch = 1
    
    if numImages > 2:
        # ==============================================
# DIC initial guess
        NewFFTSearch = funParaInput('NewFFTSearch')
        # ==============================================
# Decide DIC as accumulative or incremental mode?
        print('--- Choose accumulative or incremental mode ---  \n' % ())
        print('     0: Accumulative(By default);  \n' % ())
        print('     1: Incremental;  \n' % ())
        prompt = 'Input here: '
        DICIncOrNot = input_(prompt)
        try:
            if 0 == DICIncOrNot:
                ImgSeqIncUnit = numImages + 1
                ImgSeqIncROIUpdateOrNot = 1
            else:
                if 1 == DICIncOrNot:
                    print('Incremental mode: How many frames to update reference image once? \n' % ())
                    prompt = 'Input here: '
                    ImgSeqIncUnit = input_(prompt)
                    print('Update ROI at the same time of updating reference image? \n' % ())
                    print('    0: Do not update ROI; \n' % ())
                    print('    1: Manually(Recommended); \n' % ())
                    print('    2: Automatically; \n' % ())
                    prompt = 'Input here: '
                    ImgSeqIncROIUpdateOrNot = input_(prompt)
                else:
                    ImgSeqIncUnit = numImages + 1
                    ImgSeqIncROIUpdateOrNot = 1
        finally:
            pass
        # ================================
    else:
        ImgSeqIncUnit = numImages + 1
        ImgSeqIncROIUpdateOrNot = 1
    
    DICpara.winsize = winsize
    DICpara.winstepsize = winstepsize
    DICpara.gridxyROIRange = gridxy
    DICpara.LoadImgMethod = LoadImgMethod
    DICpara.ImgSeqIncUnit = ImgSeqIncUnit
    DICpara.ImgSeqIncROIUpdateOrNot = ImgSeqIncROIUpdateOrNot
    DICpara.Subpb2FDOrFEM = Subpb2FDOrFEM
    DICpara.NewFFTSearch = NewFFTSearch
    # DICpara.ClusterNo = ClusterNo;
    DICpara.ImgSize = Img[0].shape
    return file_name,Img,DICpara
    
    return file_name,Img,DICpara