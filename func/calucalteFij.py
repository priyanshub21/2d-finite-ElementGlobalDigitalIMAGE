# Computes the deformation gradient tensor (Fij) and Jacobian of the
# deformation gradient tensor
#
#
#
# INPUTS :
# -------------------------------------------------------------------------
#
#   u_         = displacement field vector with rigid drift removed.
#                       Format: cell array, each containing a 3D matrix for each
#                       time point (components in x,y,z)
#                           unew{time}{1} = displacement in x-direction
#                           unew{time}{2} = displacement in y-direction
#                           unew{time}{3} = displacement in z-direction
#                           unew{time}{4} = magnitude
#   dm         = meshgrid spacing
#   m2vx       = meter to pixel conversion of original images in [x y z]
#   type       = spatial differentiation kernel used for gradientN
#                options: 'fb', 'prewitt', 'sobel', 'scharr',
#                       'stencil', or 'optimal#' (# are odd numbers from 5
#                       to 19)
#                Suggested: 'optimal 9'
#
#
# OUTPUTS
# -------------------------------------------------------------------------
#   Fij     = deformation gradient tensor calculated on the mesh grid
#               with spacing dm.
#               Format: Fij{time}{3x3 deformation gradient tensor matrix}
#   J       = determinant of the deformation gradient tensor defined on the
#               mesh grid with spacing dm.
#               Format: J{time}
#
#
# NOTES
# -------------------------------------------------------------------------
# none
#
#
##
import numpy as np
    
def calculateFij(varargin = None): 
    #Establish variables and inputs
    u,spacing,m2vx,type_ = parseInputs(varargin[:])
    maxTime = len(u)
    Fij = cell(maxTime,1)
    J = cell(maxTime,1)
    for i in np.arange(1,maxTime+1).reshape(-1):
        Fij[i] = funCalculateFij(u,spacing,m2vx,type_)
        J[i] = funCalculateJ(Fij[i])
    
    return Fij,J
    
    
def funCalculateFij(u = None,spacing = None,m2vx = None,type_ = None): 
    # Calculate Displacment Gradient
    Fij = cell(3,3)
    for i in np.arange(1,3+1).reshape(-1):
        Fij[i,1],Fij[i,2],Fij[i,3] = gradientN(u[i][0] * m2vx(i),type_)
    
    #Calculate Deformation Gradient
    for i in np.arange(1,3+1).reshape(-1):
        for j in np.arange(1,3+1).reshape(-1):
            try:
                Fij[i,j] = Fij[i,j] / (spacing(j) * m2vx(j))
            finally:
                pass
    
    # for k = 1:9, Fij{k} = Fij{k}/spacing; end
    for i in np.arange(1,3+1).reshape(-1):
        Fij[i,i] = Fij[i,i] + 1
    
    return Fij
    
    
def funCalculateJ(Fij = None): 
    #Calculate Jacobian of Deformation gradient
    J = np.multiply(np.multiply(Fij[1,1],Fij[2,2]),Fij[3,3])
    J = J + np.multiply(np.multiply(Fij[1,2],Fij[2,3]),Fij[3,1])
    J = J + np.multiply(np.multiply(Fij[1,3],Fij[2,1]),Fij[3,2])
    J = J - np.multiply(np.multiply(Fij[1,3],Fij[2,2]),Fij[3,1])
    J = J - np.multiply(np.multiply(Fij[1,2],Fij[2,1]),Fij[3,3])
    J = J - np.multiply(np.multiply(Fij[1,1],Fij[2,3]),Fij[3,2])
    negJ = (J < 0)
    if sum(negJ) > 0:
        J[J < 0] = NaN
        J = inpaint_nans3(J,1)
        J[J < 0] = NaN
        J = fillmissing(J,'nearest')
    
    return J
    
    
def parseInputs(varargin = None): 
    u = varargin[0]
    spacing = varargin[2]
    m2vx = varargin[3]
    if len(varargin) < 4:
        type_ = 'optimal9'
    else:
        type_ = varargin[4]
    
    for i in np.arange(1,len(u)+1).reshape(-1):
        u[i] = cellfun(double,u[i],'UniformOutput',0)
    
    return u,spacing,m2vx,type_
    
    return Fij,J