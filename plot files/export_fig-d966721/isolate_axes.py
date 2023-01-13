import numpy as np
import os
    
def isolate_axes(ah = None,vis = None): 
    #ISOLATE_AXES Isolate the specified axes in a figure on their own
    
    # Examples:
#   fh = isolate_axes(ah)
#   fh = isolate_axes(ah, vis)
    
    # This function will create a new figure containing the axes/uipanels
# specified, and also their associated legends and colorbars. The objects
# specified must all be in the same figure, but they will generally only be
# a subset of the objects in the figure.
    
    # IN:
#    ah - An array of axes and uipanel handles, which must come from the
#         same figure.
#    vis - A boolean indicating whether the new figure should be visible.
#          Default: false.
    
    # OUT:
#    fh - The handle of the created figure.
    
    # Copyright (C) Oliver Woodford 2011-2013
    
    # Thank you to Rosella Blatt for reporting a bug to do with axes in GUIs
# 16/03/12: Moved copyfig to its own function. Thanks to Bob Fratantonio
#           for pointing out that the function is also used in export_fig.m
# 12/12/12: Add support for isolating uipanels. Thanks to michael for suggesting it
# 08/10/13: Bug fix to allchildren suggested by Will Grant (many thanks!)
# 05/12/13: Bug fix to axes having different units. Thanks to Remington Reid for reporting
# 21/04/15: Bug fix for exporting uipanels with legend/colorbar on HG1 (reported by Alvaro
#           on FEX page as a comment on 24-Apr-2014); standardized indentation & help section
# 22/04/15: Bug fix: legends and colorbars were not exported when exporting axes handle in HG2
    
    # Make sure we have an array of handles
    if not np.all(ishandle(ah)) :
        raise Exception('ah must be an array of handles')
    
    # Check that the handles are all for axes or uipanels, and are all in the same figure
    fh = ancestor(ah(1),'figure')
    nAx = np.asarray(ah).size
    for a in np.arange(1,nAx+1).reshape(-1):
        if not ismember(get(ah(a),'Type'),np.array(['axes','uipanel'])) :
            raise Exception('All handles must be axes or uipanel handles.')
        if not ancestor(ah(a),'figure')==fh :
            raise Exception('Axes must all come from the same figure.')
    
    # Tag the objects so we can find them in the copy
    old_tag = get(ah,'Tag')
    if nAx == 1:
        old_tag = np.array([old_tag])
    
    set(ah,'Tag','ObjectToCopy')
    # Create a new figure exactly the same as the old one
    fh = copyfig(fh)
    
    if len(varargin) < 2 or not vis :
        set(fh,'Visible','off')
    
    # Reset the object tags
    for a in np.arange(1,nAx+1).reshape(-1):
        set(ah(a),'Tag',old_tag[a])
    
    # Find the objects to save
    ah = findall(fh,'Tag','ObjectToCopy')
    if np.asarray(ah).size != nAx:
        close_(fh)
        raise Exception('Incorrect number of objects found.')
    
    # Set the axes tags to what they should be
    for a in np.arange(1,nAx+1).reshape(-1):
        set(ah(a),'Tag',old_tag[a])
    
    # Keep any legends and colorbars which overlap the subplots
# Note: in HG1 these are axes objects; in HG2 they are separate objects, therefore we
#       don't test for the type, only the tag (hopefully nobody but Matlab uses them!)
    lh = findall(fh,'Tag','legend','-or','Tag','Colorbar')
    nLeg = np.asarray(lh).size
    if nLeg > 0:
        set(np.array([[ah],[lh]]),'Units','normalized')
        try:
            ax_pos = get(ah,'OuterPosition')
        finally:
            pass
        if nAx > 1:
            ax_pos = cell2mat(ax_pos)
        ax_pos[:,np.arange[3,4+1]] = ax_pos(:,np.arange(3,4+1)) + ax_pos(:,np.arange(1,2+1))
        try:
            leg_pos = get(lh,'OuterPosition')
        finally:
            pass
        if nLeg > 1:
            leg_pos = cell2mat(leg_pos)
        leg_pos[:,np.arange[3,4+1]] = leg_pos(:,np.arange(3,4+1)) + leg_pos(:,np.arange(1,2+1))
        ax_pos = shiftdim(ax_pos,- 1)
        # Overlap test
        M = np.logical_and(np.logical_and(np.logical_and(bsxfun(lt,leg_pos(:,1),ax_pos(:,:,3)),bsxfun(lt,leg_pos(:,2),ax_pos(:,:,4))),bsxfun(gt,leg_pos(:,3),ax_pos(:,:,1))),bsxfun(gt,leg_pos(:,4),ax_pos(:,:,2)))
        ah = np.array([[ah],[lh(np.any(M,2))]])
    
    # Get all the objects in the figure
    axs = findall(fh)
    # Delete everything except for the input objects and associated items
    os.delete(axs(not ismember(axs,np.array([[ah],[allchildren(ah)],[allancestors(ah)]])) ))
    return fh
    
    
def allchildren(ah = None): 
    ah = findall(ah)
    if iscell(ah):
        ah = cell2mat(ah)
    
    ah = ah
    return ah
    
    
def allancestors(ah = None): 
    ph = []
    for a in np.arange(1,np.asarray(ah).size+1).reshape(-1):
        h = get(ah(a),'parent')
        while h != 0:

            ph = np.array([[ph],[h]])
            h = get(h,'parent')

    
    return ph
    
    return fh