#READ_WRITE_ENTIRE_TEXTFILE Read or write a whole text file to/from memory
#
# Read or write an entire text file to/from memory, without leaving the
# file open if an error occurs.
#
# Reading:
#   fstrm = read_write_entire_textfile(fname)
# Writing:
#   read_write_entire_textfile(fname, fstrm)
#
#IN:
#   fname - Pathname of text file to be read in.
#   fstrm - String to be written to the file, including carriage returns.
#
#OUT:
#   fstrm - String read from the file. If an fstrm input is given the
#           output is the same as that input.

import numpy as np
    
def read_write_entire_textfile(fname = None,fstrm = None): 
    modes = np.array(['rt','wt'])
    writing = len(varargin) > 1
    fh = open(fname,modes[1 + writing])
    if fh == - 1:
        raise Exception('Unable to open file %s.',fname)
    
    try:
        if writing:
            fwrite(fh,fstrm,'char*1')
        else:
            fstrm = np.transpose(fread(fh,'*char'))
    finally:
        pass
    
    fh.close()
    return fstrm
    
    return fstrm