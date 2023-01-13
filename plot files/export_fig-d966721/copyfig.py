import warnings
import numpy as np
import os
    
def copyfig(fh = None): 
    #COPYFIG Create a copy of a figure, without changing the figure
    
    # Examples:
#   fh_new = copyfig(fh_old)
    
    # This function will create a copy of a figure, but not change the figure,
# as copyobj sometimes does, e.g. by changing legends.
    
    # IN:
#    fh_old - The handle of the figure to be copied. Default: gcf.
    
    # OUT:
#    fh_new - The handle of the created figure.
    
    # Copyright (C) Oliver Woodford 2012
    
    # 26/02/15: If temp dir is not writable, use the dest folder for temp
#           destination files (Javier Paredes)
# 15/04/15: Suppress warnings during copyobj (Dun Kirk comment on FEX page 2013-10-02)
    
    # Set the default
    if len(varargin) == 0:
        fh = gcf
    
    # Is there a legend?
    if len(findall(fh,'Type','axes','Tag','legend'))==0:
        # Safe to copy using copyobj
        oldWarn = warnings.warn('off')
        fh = copyobj(fh,0)
        warnings.warn(oldWarn)
    else:
        # copyobj will change the figure, so save and then load it instead
        tmp_nam = np.array([tempname,'.fig'])
        try:
            # Ensure that the temp dir is writable (Javier Paredes 26/2/15)
            fid = open(tmp_nam,'w')
            fwrite(fid,1)
            fid.close()
            os.delete(tmp_nam)
        finally:
            pass
        hgsave(fh,tmp_nam)
        fh = hgload(tmp_nam)
        os.delete(tmp_nam)
    
    return fh
    
    return fh