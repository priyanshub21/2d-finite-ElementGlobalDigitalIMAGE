import numpy as np
    
def MeshSetUp(x = None,y = None,DICpara = None): 
    #FUNCTION DICmesh = MeshSetUp(x,y,DICpara)
# Objective: To set up a DIC uniform FE-mesh
# ----------------------------------------------
    
    #   INPUT: x,y      DIC subsets positions
#          DICpara  DIC parameters
    
    #   OUTPUT: DICmesh Generated DIC FE-mesh {coordinatesFEM, elementsFEM, ...}
    
    # ----------------------------------------------
# Author: Jin Yang.
# Contact and support: jyang526@wisc.edu -or- aldicdvc@gmail.com
# Last time updated: 02/2020.
# ==============================================
    
    ## Initialization
    winstepsize = DICpara.winstepsize
    ImgSize = DICpara.ImgSize
    ## mesh for global method
    M = x.shape[2-1]
    N = x.shape[1-1]
    
    coordinatesFEM = np.zeros((M * N,2))
    x = np.transpose(x)
    y = np.transpose(y)
    # I have transpose x and y because Matlab is read matrix in column direction
    for i in np.arange(1,coordinatesFEM.shape[1-1]+1).reshape(-1):
        coordinatesFEM[i,:] = np.array([x(i),y(i)])
        # x is horizontal position in the image
# y is vertical position in the image
    
    elementsFEM = np.zeros(((M - 1) * (N - 1),4))
    for j in np.arange(1,N - 1+1).reshape(-1):
        for i in np.arange(1,M - 1+1).reshape(-1):
            elementsFEM[[j - 1] * [M - 1] + i,:] = np.array([(j - 1) * (M) + i(j - 1) * (M) + i + 1,j * (M) + i + 1,j * (M) + i])
    
    # mesh for local method
# N is vertically in image; M is horizontally in image;
    
    coordinates = np.zeros(((M + 1) * (N + 1),2))
    gridxtemp = np.array([[x(1,1) - 0.5 * winstepsize],[x(:,1) + 0.5 * winstepsize]])
    gridytemp = np.transpose(np.array([y(1,1) - 0.5 * winstepsize,y(1,:) + 0.5 * winstepsize]))
    clear('gridx','gridy')
    gridx,gridy = np.meshgrid(gridxtemp,gridytemp)
    gridx = np.transpose(gridx)
    gridy = np.transpose(gridy)
    for i in np.arange(1,coordinates.shape[1-1]+1).reshape(-1):
        coordinates[i,:] = np.array([gridx(i),gridy(i)])
        # x is horizontal position in the image
# y is vertical position in the image
    
    elements = np.zeros((M * N,4))
    for j in np.arange(1,N+1).reshape(-1):
        for i in np.arange(1,M+1).reshape(-1):
            elements[[j - 1] * [M] + i,:] = np.array([(j - 1) * (M + 1) + i(j - 1) * (M + 1) + i + 1,j * (M + 1) + i + 1,j * (M + 1) + i])
    
    # ======== Assign BC values ==========
# -------- dirichlet BC --------
# dirichlet = [linspace(1,(M-1),(M-1))', linspace(2,M,(M-1))' ;
#             linspace(1,1+M*(N-2),(N-1))', linspace((1+M),(1+M*(N-1)),(N-1))' ;
#             linspace(M,M*(N-1),(N-1))', linspace(2*M,M*N,(N-1))' ;
#             linspace((M*(N-1)+1), (M*N-1), (M-1))', linspace((M*(N-1)+2), (M*N), (M-1))' ];
    dirichlet = []
    # -------- neumann BC --------
# neumann = [linspace(1,(M-1),(M-1))', linspace(2,M,(M-1))', zeros(M-1,1), -ones(M-1,1) ;
#              linspace(1,1+M*(N-2),(N-1))', linspace((1+M),(1+M*(N-1)),(N-1))', -ones(N-1,1), zeros(N-1,1) ;
#              linspace(M,M*(N-1),(N-1))', linspace(2*M,M*N,(N-1))', ones(N-1,1), zeros(N-1,1) ;
#              linspace((M*(N-1)+1), (M*N-1), (M-1))', linspace((M*(N-1)+2), (M*N), (M-1))', zeros(M-1,1), ones(M-1,1) ];
    neumann = []
    ## Assign variables
    DICmesh.coordinatesFEM = coordinatesFEM
    DICmesh.elementsFEM = elementsFEM
    # DICmesh.coordinates = coordinates;
# DICmesh.elements = elements;
    DICmesh.dirichlet = dirichlet
    DICmesh.neumann = neumann
    DICmesh.x0 = x
    DICmesh.y0 = y
    DICmesh.M = M
    DICmesh.N = N
    DICmesh.y0World = (ImgSize(2) + 1 - DICmesh.y0)
    DICmesh.coordinatesFEMWorld = np.array([DICmesh.coordinatesFEM(:,1),ImgSize(2) + 1 - DICmesh.coordinatesFEM(:,2)])
    return DICmesh