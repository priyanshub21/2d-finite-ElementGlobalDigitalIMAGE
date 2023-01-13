import numpy as np
import warnings
import os
    
def print2array(fig = None,res = None,renderer = None,gs_options = None): 
    #PRINT2ARRAY  Exports a figure to an image array
    
    # Examples:
#   A = print2array
#   A = print2array(figure_handle)
#   A = print2array(figure_handle, resolution)
#   A = print2array(figure_handle, resolution, renderer)
#   A = print2array(figure_handle, resolution, renderer, gs_options)
#   [A bcol] = print2array(...)
    
    # This function outputs a bitmap image of the given figure, at the desired
# resolution.
    
    # If renderer is '-painters' then ghostcript needs to be installed. This
# can be downloaded from: http://www.ghostscript.com
    
    # IN:
#   figure_handle - The handle of the figure to be exported. Default: gcf.
#   resolution - Resolution of the output, as a factor of screen
#                resolution. Default: 1.
#   renderer - string containing the renderer paramater to be passed to
#              print. Default: '-opengl'.
#   gs_options - optional ghostscript options (e.g.: '-dNoOutputFonts'). If
#                multiple options are needed, enclose in call array: {'-a','-b'}
    
    # OUT:
#   A - MxNx3 uint8 image of the figure.
#   bcol - 1x3 uint8 vector of the background color
    
    # Copyright (C) Oliver Woodford 2008-2014, Yair Altman 2015-
#{
# 05/09/11: Set EraseModes to normal when using opengl or zbuffer
#           renderers. Thanks to Pawel Kocieniewski for reporting the issue.
# 21/09/11: Bug fix: unit8 -> uint8! Thanks to Tobias Lamour for reporting it.
# 14/11/11: Bug fix: stop using hardcopy(), as it interfered with figure size
#           and erasemode settings. Makes it a bit slower, but more reliable.
#           Thanks to Phil Trinh and Meelis Lootus for reporting the issues.
# 09/12/11: Pass font path to ghostscript.
# 27/01/12: Bug fix affecting painters rendering tall figures. Thanks to
#           Ken Campbell for reporting it.
# 03/04/12: Bug fix to median input. Thanks to Andy Matthews for reporting it.
# 26/10/12: Set PaperOrientation to portrait. Thanks to Michael Watts for
#           reporting the issue.
# 26/02/15: If temp dir is not writable, use the current folder for temp
#           EPS/TIF files (Javier Paredes)
# 27/02/15: Display suggested workarounds to internal print() error (issue #16)
# 28/02/15: Enable users to specify optional ghostscript options (issue #36)
# 10/03/15: Fixed minor warning reported by Paul Soderlind; fixed code indentation
# 28/05/15: Fixed issue #69: patches with LineWidth==0.75 appear wide (internal bug in Matlab's print() func)
# 07/07/15: Fixed issue #83: use numeric handles in HG1
#}
    
    # Generate default input arguments, if needed
    if len(varargin) < 2:
        res = 1
        if len(varargin) < 1:
            fig = gcf
    
    # Warn if output is large
    old_mode = get(fig,'Units')
    set(fig,'Units','pixels')
    px = get(fig,'Position')
    set(fig,'Units',old_mode)
    npx = np.prod(px(np.arange(3,4+1)) * res) / 1000000.0
    if npx > 30:
        # 30M pixels or larger!
        warnings.warn('MATLAB:LargeImage','print2array generating a %.1fM pixel image. This could be slow and might also cause memory problems.',npx)
    
    # Retrieve the background colour
    bcol = get(fig,'Color')
    # Set the resolution parameter
    res_str = np.array(['-r',num2str(np.ceil(get(0,'ScreenPixelsPerInch') * res))])
    # Generate temporary file name
    tmp_nam = np.array([tempname,'.tif'])
    try:
        # Ensure that the temp dir is writable (Javier Paredes 26/2/15)
        fid = open(tmp_nam,'w')
        fwrite(fid,1)
        fid.close()
        os.delete(tmp_nam)
        isTempDirOk = True
    finally:
        pass
    
    # Enable users to specify optional ghostscript options (issue #36)
    if len(varargin) > 3 and not len(gs_options)==0 :
        if iscell(gs_options):
            gs_options = sprintf(' %s',gs_options[:])
        else:
            if not ischar(gs_options) :
                raise Exception('gs_options input argument must be a string or cell-array of strings')
            else:
                gs_options = np.array([' ',gs_options])
    else:
        gs_options = ''
    
    if len(varargin) > 2 and str(renderer) == str('-painters'):
        # Print to eps file
        if isTempDirOk:
            tmp_eps = np.array([tempname,'.eps'])
        else:
            tmp_eps = fullfile(fpath,np.array([fname,'.eps']))
        print2eps(tmp_eps,fig,0,renderer,'-loose')
        try:
            # Initialize the command to export to tiff using ghostscript
            cmd_str = np.array(['-dEPSCrop -q -dNOPAUSE -dBATCH ',res_str,' -sDEVICE=tiff24nc'])
            # Set the font path
            fp = font_path()
            if not len(fp)==0 :
                cmd_str = np.array([cmd_str,' -sFONTPATH="',fp,'"'])
            # Add the filenames
            cmd_str = np.array([cmd_str,' -sOutputFile="',tmp_nam,'" "',tmp_eps,'"',gs_options])
            # Execute the ghostscript command
            ghostscript(cmd_str)
        finally:
            pass
        # Delete the intermediate file
        os.delete(tmp_eps)
        # Read in the generated bitmap
        A = imread(tmp_nam)
        # Delete the temporary bitmap file
        os.delete(tmp_nam)
        # Set border pixels to the correct colour
        if bcol=='none':
            bcol = []
        else:
            if bcol==np.array([1,1,1]):
                bcol = uint8(np.array([255,255,255]))
            else:
                for l in np.arange(1,A.shape[2-1]+1).reshape(-1):
                    if not np.all(reshape(A(:,l,:) == 255,[],1)) :
                        break
                for r in np.arange(A.shape[2-1],l+- 1,- 1).reshape(-1):
                    if not np.all(reshape(A(:,r,:) == 255,[],1)) :
                        break
                for t in np.arange(1,A.shape[1-1]+1).reshape(-1):
                    if not np.all(reshape(A(t,:,:) == 255,[],1)) :
                        break
                for b in np.arange(A.shape[1-1],t+- 1,- 1).reshape(-1):
                    if not np.all(reshape(A(b,:,:) == 255,[],1)) :
                        break
                bcol = uint8(median(single(np.array([[reshape(A(:,np.array([l,r]),:),[],A.shape[3-1])],[reshape(A(np.array([t,b]),:,:),[],A.shape[3-1])]])),1))
                for c in np.arange(1,A.shape[3-1]+1).reshape(-1):
                    A[:,np.array[[np.arange[1,l - 1+1],np.arange[r + 1,end()+1]]],c] = bcol(c)
                    A[np.array[[np.arange[1,t - 1+1],np.arange[b + 1,end()+1]]],:,c] = bcol(c)
    else:
        if len(varargin) < 3:
            renderer = '-opengl'
        err = False
        # Set paper size
        old_pos_mode = get(fig,'PaperPositionMode')
        old_orientation = get(fig,'PaperOrientation')
        set(fig,'PaperPositionMode','auto','PaperOrientation','portrait')
        try:
            # Workaround for issue #69: patches with LineWidth==0.75 appear wide (internal bug in Matlab's print() function)
            fp = []
            fp = findall(fig,'Type','patch','LineWidth',0.75)
            set(fp,'LineWidth',0.5)
            # Fix issue #83: use numeric handles in HG1
            if not using_hg2(fig) :
                fig = double(fig)
            # Print to tiff file
            print_(fig,renderer,res_str,'-dtiff',tmp_nam)
            # Read in the printed file
            A = imread(tmp_nam)
            # Delete the temporary file
            os.delete(tmp_nam)
        finally:
            pass
        set(fp,'LineWidth',0.75)
        # Reset paper size
        set(fig,'PaperPositionMode',old_pos_mode,'PaperOrientation',old_orientation)
        # Throw any error that occurred
        if err:
            # Display suggested workarounds to internal print() error (issue #16)
            2.write('An error occured with Matlab's builtin print function.\nTry setting the figure Renderer to 'painters' or use opengl('software').\n\n' % ())
            rethrow(ex)
        # Set the background color
        if bcol=='none':
            bcol = []
        else:
            bcol = bcol * 255
            if bcol==np.round(bcol):
                bcol = uint8(bcol)
            else:
                bcol = np.squeeze(A(1,1,:))
    
    # Check the output size is correct
    if res==np.round(res):
        px = np.round(np.array([px(np.array([4,3])) * res,3]))
        if not A.shape==px :
            # Correct the output size
            A = A(np.arange(1,np.amin(end(),px(1))+1),np.arange(1,np.amin(end(),px(2))+1),:)
    
    return A,bcol
    
    # Function to return (and create, where necessary) the font path
    
def font_path(): 
    fp = user_string('gs_font_path')
    if not len(fp)==0 :
        return fp
    
    # Create the path
# Start with the default path
    fp = getenv('GS_FONTPATH')
    # Add on the typical directories for a given OS
    if ispc:
        if not len(fp)==0 :
            fp = np.array([fp,';'])
        fp = np.array([fp,getenv('WINDIR'),filesep,'Fonts'])
    else:
        if not len(fp)==0 :
            fp = np.array([fp,':'])
        fp = np.array([fp,'/usr/share/fonts:/usr/local/share/fonts:/usr/share/fonts/X11:/usr/local/share/fonts/X11:/usr/share/fonts/truetype:/usr/local/share/fonts/truetype'])
    
    user_string('gs_font_path',fp)
    return fp
    
    return A,bcol