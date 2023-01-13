import numpy as np
import os
    
def eps2pdf(source = None,dest = None,crop = None,append = None,gray = None,quality = None,gs_options = None): 
    #EPS2PDF  Convert an eps file to pdf format using ghostscript
    
    # Examples:
#   eps2pdf source dest
#   eps2pdf(source, dest, crop)
#   eps2pdf(source, dest, crop, append)
#   eps2pdf(source, dest, crop, append, gray)
#   eps2pdf(source, dest, crop, append, gray, quality)
#   eps2pdf(source, dest, crop, append, gray, quality, gs_options)
    
    # This function converts an eps file to pdf format. The output can be
# optionally cropped and also converted to grayscale. If the output pdf
# file already exists then the eps file can optionally be appended as a new
# page on the end of the eps file. The level of bitmap compression can also
# optionally be set.
    
    # This function requires that you have ghostscript installed on your
# system. Ghostscript can be downloaded from: http://www.ghostscript.com
    
    # Inputs:
#   source  - filename of the source eps file to convert. The filename is
#             assumed to already have the extension ".eps".
#   dest    - filename of the destination pdf file. The filename is assumed
#             to already have the extension ".pdf".
#   crop    - boolean indicating whether to crop the borders off the pdf.
#             Default: true.
#   append  - boolean indicating whether the eps should be appended to the
#             end of the pdf as a new page (if the pdf exists already).
#             Default: false.
#   gray    - boolean indicating whether the output pdf should be grayscale
#             or not. Default: false.
#   quality - scalar indicating the level of image bitmap quality to
#             output. A larger value gives a higher quality. quality > 100
#             gives lossless output. Default: ghostscript prepress default.
#   gs_options - optional ghostscript options (e.g.: '-dNoOutputFonts'). If
#                multiple options are needed, enclose in call array: {'-a','-b'}
    
    # Copyright (C) Oliver Woodford 2009-2014, Yair Altman 2015-
    
    # Suggestion of appending pdf files provided by Matt C at:
# http://www.mathworks.com/matlabcentral/fileexchange/23629
    
    # Thank you to Fabio Viola for pointing out compression artifacts, leading
# to the quality setting.
# Thank you to Scott for pointing out the subsampling of very small images,
# which was fixed for lossless compression settings.
    
    # 9/12/2011 Pass font path to ghostscript.
# 26/02/15: If temp dir is not writable, use the dest folder for temp
#           destination files (Javier Paredes)
# 28/02/15: Enable users to specify optional ghostscript options (issue #36)
# 01/03/15: Upon GS error, retry without the -sFONTPATH= option (this might solve
#           some /findfont errors according to James Rankin, FEX Comment 23/01/15)
# 23/06/15: Added extra debug info in case of ghostscript error; code indentation
    
    # Intialise the options string for ghostscript
    options = np.array(['-q -dNOPAUSE -dBATCH -sDEVICE=pdfwrite -dPDFSETTINGS=/prepress -sOutputFile="',dest,'"'])
    # Set crop option
    if len(varargin) < 3 or crop:
        options = np.array([options,' -dEPSCrop'])
    
    # Set the font path
    fp = font_path()
    if not len(fp)==0 :
        options = np.array([options,' -sFONTPATH="',fp,'"'])
    
    # Set the grayscale option
    if len(varargin) > 4 and gray:
        options = np.array([options,' -sColorConversionStrategy=Gray -dProcessColorModel=/DeviceGray'])
    
    # Set the bitmap quality
    if len(varargin) > 5 and not len(quality)==0 :
        options = np.array([options,' -dAutoFilterColorImages=false -dAutoFilterGrayImages=false'])
        if quality > 100:
            options = np.array([options,' -dColorImageFilter=/FlateEncode -dGrayImageFilter=/FlateEncode -c ".setpdfwrite << /ColorImageDownsampleThreshold 10 /GrayImageDownsampleThreshold 10 >> setdistillerparams"'])
        else:
            options = np.array([options,' -dColorImageFilter=/DCTEncode -dGrayImageFilter=/DCTEncode'])
            v = 1 + (quality < 80)
            quality = 1 - quality / 100
            s = sprintf('<< /QFactor %.2f /Blend 1 /HSample [%d 1 1 %d] /VSample [%d 1 1 %d] >>',quality,v,v,v,v)
            options = sprintf('%s -c ".setpdfwrite << /ColorImageDict %s /GrayImageDict %s >> setdistillerparams"',options,s,s)
    
    # Enable users to specify optional ghostscript options (issue #36)
    if len(varargin) > 6 and not len(gs_options)==0 :
        if iscell(gs_options):
            gs_options = sprintf(' %s',gs_options[:])
        else:
            if not ischar(gs_options) :
                raise Exception('gs_options input argument must be a string or cell-array of strings')
            else:
                gs_options = np.array([' ',gs_options])
        options = np.array([options,gs_options])
    
    # Check if the output file exists
    if len(varargin) > 3 and append and os.path.exist(str(dest)) == 2:
        # File exists - append current figure to the end
        tmp_nam = tempname
        try:
            # Ensure that the temp dir is writable (Javier Paredes 26/2/15)
            fid = open(tmp_nam,'w')
            fwrite(fid,1)
            fid.close()
            os.delete(tmp_nam)
        finally:
            pass
        # Copy the file
        copyfile(dest,tmp_nam)
        # Add the output file names
        options = np.array([options,' -f "',tmp_nam,'" "',source,'"'])
        try:
            # Convert to pdf using ghostscript
            status,message = ghostscript(options)
        finally:
            pass
        # Delete the intermediate file
        os.delete(tmp_nam)
    else:
        # File doesn't exist or should be over-written
# Add the output file names
        options = np.array([options,' -f "',source,'"'])
        # Convert to pdf using ghostscript
        status,message = ghostscript(options)
    
    # Check for error
    if status:
        # Retry without the -sFONTPATH= option (this might solve some GS
# /findfont errors according to James Rankin, FEX Comment 23/01/15)
        orig_options = options
        if not len(fp)==0 :
            options = regexprep(options,' -sFONTPATH=[^ ]+ ',' ')
            status = ghostscript(options)
            if not status :
                return
        # Report error
        if len(message)==0:
            raise Exception('Unable to generate pdf. Check destination directory is writable.')
        else:
            2.write('Ghostscript error: perhaps %s is open by another application\n' % (dest))
            if not len(gs_options)==0 :
                2.write('  or maybe the%s option(s) are not accepted by your GS version\n' % (gs_options))
            2.write('Ghostscript options: %s\n\n' % (orig_options))
            raise Exception(message)
    
    return
    
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
    