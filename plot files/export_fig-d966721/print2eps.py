import numpy as np
import warnings
    
def print2eps(name = None,fig = None,export_options = None,varargin = None): 
    #PRINT2EPS  Prints figures to eps with improved line styles
    
    # Examples:
#   print2eps filename
#   print2eps(filename, fig_handle)
#   print2eps(filename, fig_handle, export_options)
#   print2eps(filename, fig_handle, export_options, print_options)
    
    # This function saves a figure as an eps file, with two improvements over
# MATLAB's print command. First, it improves the line style, making dashed
# lines more like those on screen and giving grid lines a dotted line style.
# Secondly, it substitutes original font names back into the eps file,
# where these have been changed by MATLAB, for up to 11 different fonts.
    
    #IN:
#   filename - string containing the name (optionally including full or
#              relative path) of the file the figure is to be saved as. A
#              ".eps" extension is added if not there already. If a path is
#              not specified, the figure is saved in the current directory.
#   fig_handle - The handle of the figure to be saved. Default: gcf().
#   export_options - array of optional scalar values:
#       bb_padding - Scalar value of amount of padding to add to border around
#                    the cropped image, in points (if >1) or percent (if <1).
#                    Can be negative as well as positive; Default: 0
#       crop       - Crop amount. Deafult: 0
#       fontswap   - Whether to swap non-default fonts in figure. Default: true
#   print_options - Additional parameter strings to be passed to the print command
    
    #{
# Copyright (C) Oliver Woodford 2008-2014, Yair Altman 2015-
    
    # The idea of editing the EPS file to change line styles comes from Jiro
# Doke's FIXPSLINESTYLE (fex id: 17928)
# The idea of changing dash length with line width came from comments on
# fex id: 5743, but the implementation is mine :)
#}
#{
# 14/11/11: Fix a MATLAB bug rendering black or white text incorrectly.
#           Thanks to Mathieu Morlighem for reporting the issue and
#           obtaining a fix from TMW.
# 08/12/11: Added ability to correct fonts. Several people have requested
#           this at one time or another, and also pointed me to printeps
#           (fex id: 7501), so thank you to them. My implementation (which
#           was not inspired by printeps - I'd already had the idea for my
#           approach) goes slightly further in that it allows multiple
#           fonts to be swapped.
# 14/12/11: Fix bug affecting font names containing spaces. Thanks to David
#           Szwer for reporting the issue.
# 25/01/12: Add a font not to be swapped. Thanks to Anna Rafferty and Adam
#           Jackson for reporting the issue. Also fix a bug whereby using a
#           font alias can lead to another font being swapped in.
# 10/04/12: Make the font swapping case insensitive.
# 26/10/12: Set PaperOrientation to portrait. Thanks to Michael Watts for
#           reporting the issue.
# 26/10/12: Fix issue to do with swapping fonts changing other fonts and
#           sizes we don't want, due to listeners. Thanks to Malcolm Hudson
#           for reporting the issue.
# 22/03/13: Extend font swapping to axes labels. Thanks to Rasmus Ischebeck
#           for reporting the issue.
# 23/07/13: Bug fix to font swapping. Thanks to George for reporting the
#           issue.
# 13/08/13: Fix MATLAB feature of not exporting white lines correctly.
#           Thanks to Sebastian Heßlinger for reporting it.
# 24/02/15: Fix for Matlab R2014b bug (issue #31): LineWidths<0.75 are not
#           set in the EPS (default line width is used)
# 25/02/15: Fixed issue #32: BoundingBox problem caused uncropped EPS/PDF files
# 05/03/15: Fixed issue #43: Inability to perform EPS file post-processing
# 06/03/15: Improved image padding & cropping thanks to Oscar Hartogensis
# 21/03/15: Fixed edge-case of missing handles having a 'FontName' property
# 26/03/15: Attempt to fix issue #45: white lines in subplots do not print correctly
# 27/03/15: Attempt to fix issue #44: white artifact lines appearing in patch exports
# 30/03/15: Fixed issue #52: improved performance on HG2 (R2014b+)
# 09/04/15: Comment blocks consolidation and minor code cleanup (no real code change)
# 12/04/15: Fixed issue #56: bad cropping
# 14/04/15: Workaround for issue #45: lines in image subplots are exported in invalid color
# 07/07/15: Added option to avoid font-swapping in EPS/PDF
# 07/07/15: Fixed issue #83: use numeric handles in HG1
# 22/07/15: Fixed issue #91 (thanks to Carlos Moffat)
#}
    
    options = np.array(['-loose'])
    if len(varargin) > 3:
        options = np.array([options,varargin])
    else:
        if len(varargin) < 3:
            export_options = 0
            if len(varargin) < 2:
                fig = gcf()
    
    # Retrieve padding, crop & font-swap values
    if np.asarray(export_options).size > 2:
        fontswap = export_options(3)
    else:
        fontswap = True
    
    if np.asarray(export_options).size > 1:
        bb_crop = export_options(2)
    else:
        bb_crop = 0
    
    if np.asarray(export_options).size > 0:
        bb_padding = export_options(1)
    else:
        bb_padding = 0
    
    # Construct the filename
    if np.asarray(name).size < 5 or not strcmpi(name(np.arange(end() - 3,end()+1)),'.eps') :
        name = np.array([name,'.eps'])
    
    # Set paper size
    old_pos_mode = get(fig,'PaperPositionMode')
    old_orientation = get(fig,'PaperOrientation')
    set(fig,'PaperPositionMode','auto','PaperOrientation','portrait')
    # Find all the used fonts in the figure
    font_handles = findall(fig,'-property','FontName')
    fonts = get(font_handles,'FontName')
    if len(fonts)==0:
        fonts = np.array([])
    else:
        if not iscell(fonts) :
            fonts = np.array([fonts])
    
    # Map supported font aliases onto the correct name
    fontsl = fonts.lower()
    for a in np.arange(1,np.asarray(fonts).size+1).reshape(-1):
        f = fontsl[a]
        f[f == ' '] = []
        if np.array(['times','timesnewroman','times-roman']) == f:
            fontsl[a] = 'times-roman'
        else:
            if np.array(['arial','helvetica']) == f:
                fontsl[a] = 'helvetica'
            else:
                if np.array(['newcenturyschoolbook','newcenturyschlbk']) == f:
                    fontsl[a] = 'newcenturyschlbk'
    
    fontslu = unique(fontsl)
    # Determine the font swap table
    if fontswap:
        matlab_fonts = np.array(['Helvetica','Times-Roman','Palatino','Bookman','Helvetica-Narrow','Symbol','AvantGarde','NewCenturySchlbk','Courier','ZapfChancery','ZapfDingbats'])
        matlab_fontsl = matlab_fonts.lower()
        require_swap = find(not ismember(fontslu,matlab_fontsl) )
        unused_fonts = find(not ismember(matlab_fontsl,fontslu) )
        font_swap = cell(3,np.amin(np.asarray(require_swap).size,np.asarray(unused_fonts).size))
        fonts_new = fonts
        for a in np.arange(1,font_swap.shape[2-1]+1).reshape(-1):
            font_swap[1,a] = find(str(fontslu[require_swap(a)]) == str(fontsl))
            font_swap[2,a] = matlab_fonts[unused_fonts(a)]
            font_swap[3,a] = fonts[font_swap[1,a](1)]
            fonts_new[font_swap[1,a]] = font_swap(2,a)
    else:
        font_swap = []
    
    # Swap the fonts
    if not len(font_swap)==0 :
        fonts_size = get(font_handles,'FontSize')
        if iscell(fonts_size):
            fonts_size = cell2mat(fonts_size)
        M = False(font_handles.shape)
        # Loop because some changes may not stick first time, due to listeners
        c = 0
        update = np.zeros((1000,1))
        for b in np.arange(1,10+1).reshape(-1):
            for a in np.arange(1,np.asarray(M).size+1).reshape(-1):
                M[a] = not get(font_handles(a),'FontName')==fonts_new[a]  or not get(font_handles(a),'FontSize')==fonts_size(a) 
                if M(a):
                    set(font_handles(a),'FontName',fonts_new[a],'FontSize',fonts_size(a))
                    c = c + 1
                    update[c] = a
            if not np.any(M) :
                break
        # Compute the order to revert fonts later, without the need of a loop
        update,M = unique(update(np.arange(1,c+1)))
        M,M = __builtint__.sorted(M)
        update = reshape(update(M),1,[])
    
    # MATLAB bug fix - black and white text can come out inverted sometimes
# Find the white and black text
    black_text_handles = findall(fig,'Type','text','Color',np.array([0,0,0]))
    white_text_handles = findall(fig,'Type','text','Color',np.array([1,1,1]))
    # Set the font colors slightly off their correct values
    set(black_text_handles,'Color',np.array([0,0,0]) + eps)
    set(white_text_handles,'Color',np.array([1,1,1]) - eps)
    # MATLAB bug fix - white lines can come out funny sometimes
# Find the white lines
    white_line_handles = findall(fig,'Type','line','Color',np.array([1,1,1]))
    # Set the line color slightly off white
    set(white_line_handles,'Color',np.array([1,1,1]) - 1e-05)
    # Workaround for issue #45: lines in image subplots are exported in invalid color
# In this case the -depsc driver solves the problem, but then all the other workarounds
# below (for all the other issues) will fail, so it's better to let the user decide by
# just issuing a warning and accepting the '-depsc' input parameter
    epsLevel2 = not np.any(strcmpi(options,'-depsc')) 
    if epsLevel2:
        # Use -depsc2 (EPS color level-2) if -depsc (EPS color level-3) was not specifically requested
        options[end() + 1] = '-depsc2'
        # Issue a warning if multiple images & lines were found in the figure, and HG1 with painters renderer is used
        isPainters = np.any(strcmpi(options,'-painters'))
        if isPainters and not using_hg2  and np.asarray(findall(fig,'Type','image')).size > 1 and not len(findall(fig,'Type','line'))==0 :
            warnings.warn('YMA:export_fig:issue45',np.array(['Multiple images & lines detected. In such cases, the lines might \n','appear with an invalid color due to an internal MATLAB bug (fixed in R2014b). \n','Possible workaround: add a '-depsc' or '-opengl' parameter to the export_fig command.']))
    
    # Fix issue #83: use numeric handles in HG1
    if not using_hg2(fig) :
        fig = double(fig)
    
    # Print to eps file
    print_(fig,options[:],name)
    # Do post-processing on the eps file
    try:
        # Read the EPS file into memory
        fstrm = read_write_entire_textfile(name)
    finally:
        pass
    
    # Fix for Matlab R2014b bug (issue #31): LineWidths<0.75 are not set in the EPS (default line width is used)
    try:
        if not len(fstrm)==0  and using_hg2(fig):
            # Convert miter joins to line joins
#fstrm = regexprep(fstrm, '\n10.0 ML\n', '\n1 LJ\n');
# This is faster (the original regexprep could take many seconds when the axes contains many lines):
            fstrm = fstrm.replace(sprintf('\n10.0 ML\n'),sprintf('\n1 LJ\n'))
            # In HG2, grid lines and axes Ruler Axles have a default LineWidth of 0.5 => replace en-bulk (assume that 1.0 LineWidth = 1.333 LW)
#   hAxes=gca; hAxes.YGridHandle.LineWidth, hAxes.YRuler.Axle.LineWidth
#fstrm = regexprep(fstrm, '(GC\n2 setlinecap\n1 LJ)\nN', '$1\n0.667 LW\nN');
# This is faster:
            fstrm = fstrm.replace(sprintf('GC\n2 setlinecap\n1 LJ\nN'),sprintf('GC\n2 setlinecap\n1 LJ\n0.667 LW\nN'))
            # This is more accurate but *MUCH* slower (issue #52)
#{
# Modify all thin lines in the figure to have 10x LineWidths
            hLines = findall(fig,'Type','line')
            hThinLines = []
            for lineIdx in np.arange(1,np.asarray(hLines).size+1).reshape(-1):
                thisLine = hLines(lineIdx)
                if thisLine.LineWidth < 0.75 and strcmpi(thisLine.Visible,'on'):
                    hThinLines[end() + 1] = thisLine
                    thisLine.LineWidth = thisLine.LineWidth * 10
            # If any thin lines were found
            if not len(hThinLines)==0 :
                # Prepare an EPS with large-enough line widths
                print_(fig,options[:],name)
                # Restore the original LineWidths in the figure
                for lineIdx in np.arange(1,np.asarray(hThinLines).size+1).reshape(-1):
                    thisLine = handle(hThinLines(lineIdx))
                    thisLine.LineWidth = thisLine.LineWidth / 10
                # Compare the original and the new EPS files and correct the original stream's LineWidths
                fstrm_new = read_write_entire_textfile(name)
                idx = 500
                markerStr = sprintf('10.0 ML\nN')
                markerLen = len(markerStr)
                while not len(idx)==0  and idx < len(fstrm):

                    lastIdx = np.amin(len(fstrm),len(fstrm_new))
                    delta = fstrm(np.arange(idx + 1,lastIdx+1)) - fstrm_new(np.arange(idx + 1,lastIdx+1))
                    idx = idx + find(delta,1)
                    if not len(idx)==0  and fstrm(np.arange(idx - markerLen + 1,idx+1))==markerStr and not len(regexp(fstrm_new(np.arange(idx - markerLen + 1,idx + 12+1)),'10.0 ML\n[\d\.]+ LW\nN'))==0 :
                        value = str2double(regexprep(fstrm_new(np.arange(idx,idx + 12+1)),' .*',''))
                        if np.isnan(value):
                            break
                        newStr = sprintf('%0.3f LW\n',value / 10)
                        fstrm = np.array([fstrm(np.arange(1,idx - 1+1)),newStr,fstrm(np.arange(idx,end()+1))])
                        idx = idx + 12
                    else:
                        break

            #}
            # This is much faster although less accurate: fix all non-gray lines to have a LineWidth of 0.75 (=1 LW)
# Note: This will give incorrect LineWidth of 075 for lines having LineWidth<0.75, as well as for non-gray grid-lines (if present)
#       However, in practice these edge-cases are very rare indeed, and the difference in LineWidth should not be noticeable
#fstrm = regexprep(fstrm, '([CR]C\n2 setlinecap\n1 LJ)\nN', '$1\n1 LW\nN');
# This is faster (the original regexprep could take many seconds when the axes contains many lines):
            fstrm = fstrm.replace(sprintf('\n2 setlinecap\n1 LJ\nN'),sprintf('\n2 setlinecap\n1 LJ\n1 LW\nN'))
    finally:
        pass
    
    # Reset the font and line colors
    set(black_text_handles,'Color',np.array([0,0,0]))
    set(white_text_handles,'Color',np.array([1,1,1]))
    set(white_line_handles,'Color',np.array([1,1,1]))
    # Reset paper size
    set(fig,'PaperPositionMode',old_pos_mode,'PaperOrientation',old_orientation)
    # Reset the font names in the figure
    if not len(font_swap)==0 :
        for a in update.reshape(-1):
            set(font_handles(a),'FontName',fonts[a],'FontSize',fonts_size(a))
    
    # Bail out if EPS post-processing is not possible
    if len(fstrm)==0:
        warnings.warn('Loading EPS file failed, so unable to perform post-processing. This is usually because the figure contains a large number of patch objects. Consider exporting to a bitmap format in this case.')
        return
    
    # Replace the font names
    if not len(font_swap)==0 :
        for a in np.arange(1,font_swap.shape[2-1]+1).reshape(-1):
            #fstrm = regexprep(fstrm, [font_swap{1,a} '-?[a-zA-Z]*\>'], font_swap{3,a}(~isspace(font_swap{3,a})));
            fstrm = regexprep(fstrm,font_swap[2,a],font_swap[3,a](not isspace(font_swap[3,a]) ))
    
    # Move the bounding box to the top of the file (HG2 only), or fix the line styles (HG1 only)
    if using_hg2(fig):
        # Move the bounding box to the top of the file (HG2 only)
        s,e = regexp(fstrm,'%%BoundingBox: [^%]*%%')
        if np.asarray(s).size == 2:
            fstrm = fstrm(np.array([np.arange(1,s(1) - 1+1),np.arange(s(2),e(2) - 2+1),np.arange(e(1) - 1,s(2) - 1+1),np.arange(e(2) - 1,end()+1)]))
    else:
        # Fix the line styles (HG1 only)
        fstrm = fix_lines(fstrm)
    
    # Apply the bounding box padding & cropping, replacing Matlab's print()'s bounding box
    if bb_crop:
        # Calculate a new bounding box based on a bitmap print using crop_border.m
# 1. Determine the Matlab BoundingBox and PageBoundingBox
        s,e = regexp(fstrm,'%%BoundingBox: [^%]*%%')
        if np.asarray(s).size == 2:
            s = s(2)
            e = e(2)
        aa = fstrm(np.arange(s + 15,e - 3+1))
        bb_matlab = cell2mat(textscan(aa,'%f32%f32%f32%f32'))
        s,e = regexp(fstrm,'%%PageBoundingBox: [^%]*%%')
        if np.asarray(s).size == 2:
            s = s(2)
            e = e(2)
        aa = fstrm(np.arange(s + 19,e - 3+1))
        pagebb_matlab = cell2mat(textscan(aa,'%f32%f32%f32%f32'))
        # 2. Create a bitmap image and use crop_borders to create the relative
#    bb with respect to the PageBoundingBox
        A,bcol = print2array(fig,1,'-opengl')
        aa,aa,aa,bb_rel = crop_borders(A,bcol,bb_padding)
        # 3. Calculate the new Bounding Box
        pagew = pagebb_matlab(3) - pagebb_matlab(1)
        pageh = pagebb_matlab(4) - pagebb_matlab(2)
        #bb_new = [pagebb_matlab(1)+pagew*bb_rel(1) pagebb_matlab(2)+pageh*bb_rel(2) ...
#          pagebb_matlab(1)+pagew*bb_rel(3) pagebb_matlab(2)+pageh*bb_rel(4)];
        bb_new = pagebb_matlab(np.array([1,2,1,2])) + np.multiply(np.array([pagew,pageh,pagew,pageh]),bb_rel)
        bb_offset = (bb_new - bb_matlab) + np.array([- 1,- 1,1,1])
        # Apply the bounding box padding
        if bb_padding:
            if np.abs(bb_padding) < 1:
                bb_padding = np.round((mean(np.array([bb_new(3) - bb_new(1),bb_new(4) - bb_new(2)])) * bb_padding) / 0.5) * 0.5
            add_padding = lambda n1 = None,n2 = None,n3 = None,n4 = None: sprintf(' %d',str2double(np.array([n1,n2,n3,n4])) + np.array([- bb_padding,- bb_padding,bb_padding,bb_padding]) + bb_offset)
        else:
            add_padding = lambda n1 = None,n2 = None,n3 = None,n4 = None: sprintf(' %d',str2double(np.array([n1,n2,n3,n4])) + bb_offset)
        fstrm = regexprep(fstrm,'%%BoundingBox:[ ]+([-]?\d+)[ ]+([-]?\d+)[ ]+([-]?\d+)[ ]+([-]?\d+)','%%BoundingBox:${add_padding($1, $2, $3, $4)}')
    
    # Fix issue #44: white artifact lines appearing in patch exports
# Note: the problem is due to the fact that Matlab's print() function exports patches
#       as a combination of filled triangles, and a white line appears where the triangles touch
# In the workaround below, we will modify such dual-triangles into a filled rectangle.
# We are careful to only modify regexps that exactly match specific patterns - it's better to not
# correct some white-line artifacts than to change the geometry of a patch, or to corrupt the EPS.
#   e.g.: '0 -450 937 0 0 450 3 MP PP 937 0 0 -450 0 450 3 MP PP' => '0 -450 937 0 0 450 0 0 4 MP'
    fstrm = regexprep(fstrm,'\n([-\d.]+ [-\d.]+) ([-\d.]+ [-\d.]+) ([-\d.]+ [-\d.]+) 3 MP\nPP\n\2 \1 \3 3 MP\nPP\n','\n$1 $2 $3 0 0 4 MP\nPP\n')
    fstrm = regexprep(fstrm,'\n([-\d.]+ [-\d.]+) ([-\d.]+ [-\d.]+) ([-\d.]+ [-\d.]+) 3 MP\nPP\n\2 \3 \1 3 MP\nPP\n','\n$1 $2 $3 0 0 4 MP\nPP\n')
    fstrm = regexprep(fstrm,'\n([-\d.]+ [-\d.]+) ([-\d.]+ [-\d.]+) ([-\d.]+ [-\d.]+) 3 MP\nPP\n\3 \1 \2 3 MP\nPP\n','\n$1 $2 $3 0 0 4 MP\nPP\n')
    fstrm = regexprep(fstrm,'\n([-\d.]+ [-\d.]+) ([-\d.]+ [-\d.]+) ([-\d.]+ [-\d.]+) 3 MP\nPP\n\3 \2 \1 3 MP\nPP\n','\n$1 $2 $3 0 0 4 MP\nPP\n')
    # Write out the fixed eps file
    read_write_entire_textfile(name,fstrm)
    return
    