import numpy as np
    
def funDerivativeOp(M = None,N = None,h = None): 
    # FUNCTION A = funDerivativeOp(M,N,h)
# To generate a first order gradient operator {A} in the 2D case such that
# {F} = {A}{U} where displacement vector U = [Ux_node1, Uy_node1, ... , Ux_nodeN, Uy_nodeN]',
# and deformation gradient F = [F11_node1, F21_node1, F12_node1, F22_node1, ...
# ... , F11_nodeN, F21_nodeN, F12_nodeN, F22_nodeN]';
# ----------------------------------------------
    
    #   INPUT: M,N   Mesh size in the x- and y-directions
    
    #   OUTPUT: A    Finite difference operator
    
    # ----------------------------------------------
# Author: Jin Yang.
# Contact and support: jyang526@wisc.edu -or- aldicdvc@gmail.com
# Last time updated: 2018.03, 2020.12
# ==============================================
    
    ## Initialization
    DIM = 2
    A = sparse(DIM ** 2 * M * N,DIM * M * N)
    sizeA = A.shape
    XX,YY = ndgrid(np.arange(1,M+1),np.arange(1,N+1))
    XX = XX
    YY = YY
    INDEXI = np.ones((2 * DIM ** 2 * M * N,1))
    INDEXJ = INDEXI
    INDEXVAL = INDEXI
    ## ------ Inside assign values ------
    for tempi in np.arange(1,M * N+1).reshape(-1):
        tempx = XX(tempi)
        tempy = YY(tempi)
        #     for tempy = 1:N #####
#         for tempx = 1:M #####
        # Determine BorderOrNot
        BorderOrNot = np.ones((1,DIM * 2))
        if tempy == 1:
            BorderOrNot[2] = 0
        if tempy == N:
            BorderOrNot[1] = 0
        if tempx == 1:
            BorderOrNot[4] = 0
        if tempx == M:
            BorderOrNot[3] = 0
        index = tempx + M * (tempy - 1)
        indexNeighbors = index * np.ones((1,DIM * 2)) + np.multiply(BorderOrNot,np.array([M,- M,1,- 1]))
        # indexBack = index+M; indexFront = index-M; indexLeft = index-1;
# indexRight = index+1; indexTop = index+M*N; indexBottom = index-M*N;
        # Find index of affine deformation gradient tensor {F}
        indexFrow = DIM ** 2 * index * np.ones((1,DIM ** 2)) + np.array([np.arange(- (DIM ** 2 - 1),0+1,1)])
        # ##### Just some help lines #####
# indexF11col1 = 3*indexLeft-2; indexF11col2 = 3*indexRight-2;
# indexF12col1 = 3*indexFront-2; indexF12col2 = 3*indexBack-2;
# indexF13col1 = 3*indexBottom-2; indexF13col2 = 3*indexTop-2;
        # indexF21col1 = 3*indexLeft-1; indexF21col2 = 3*indexRight-1;
# indexF22col1 = 3*indexFront-1; indexF22col2 = 3*indexBack-1;
# indexF23col1 = 3*indexBottom-1; indexF23col2 = 3*indexTop-1;
        # indexF31col1 = 3*indexLeft; indexF31col2 = 3*indexRight;
# indexF32col1 = 3*indexFront; indexF32col2 = 3*indexBack;
# indexF33col1 = 3*indexBottom; indexF33col2 = 3*indexTop;
        # indexFcol1 = [3*indexLeft-2,3*indexLeft-1,3*indexLeft, ...
# 3*indexFront-2,3*indexFront-1,3*indexFront,...
# 3*indexBottom-2,3*indexBottom-1,3*indexBottom];
        indexFcol1 = np.array([DIM * indexNeighbors(4) * np.ones((1,DIM)) + np.array([np.arange(- (DIM - 1),0+1,1)]),DIM * indexNeighbors(2) * np.ones((1,DIM)) + np.array([np.arange(- (DIM - 1),0+1,1)])])
        indexFcol2 = np.array([DIM * indexNeighbors(3) * np.ones((1,DIM)) + np.array([np.arange(- (DIM - 1),0+1,1)]),DIM * indexNeighbors(1) * np.ones((1,DIM)) + np.array([np.arange(- (DIM - 1),0+1,1)])])
        # indexF1 = sub2ind(sizeA,indexFrow,indexFcol1);
# indexF2 = sub2ind(sizeA,indexFrow,indexFcol2);
        INDEXI[np.arange[2 * 4 * tempi - 7,2 * 4 * tempi - 0+1]] = np.transpose(np.array([indexFrow,indexFrow]))
        INDEXJ[np.arange[2 * 4 * tempi - 7,2 * 4 * tempi - 0+1]] = np.transpose(np.array([indexFcol1,indexFcol2]))
        INDEXVAL[np.arange[2 * 4 * tempi - 7,2 * 4 * tempi - 4+1]] = - np.ones((4,1))
        if BorderOrNot(3) * BorderOrNot(4) == 0:
            INDEXVAL[np.arange[2 * 4 * tempi - 7,2 * 4 * tempi - 6+1]] = - 2 * np.ones((2,1))
            INDEXVAL[np.arange[2 * 4 * tempi - 3,2 * 4 * tempi - 2+1]] = 2 * np.ones((2,1))
        if BorderOrNot(1) * BorderOrNot(2) == 0:
            INDEXVAL[np.arange[2 * 4 * tempi - 5,2 * 4 * tempi - 4+1]] = - 2 * np.ones((2,1))
            INDEXVAL[np.arange[2 * 4 * tempi - 1,2 * 4 * tempi+1]] = 2 * np.ones((2,1))
        # A(indexF1) = -1; A(indexF2) = 1;
# if BorderOrNot(3)*BorderOrNot(4) == 0, A(indexF1(1:DIM)) = -2; A(indexF2(1:DIM)) = 2; end
#if BorderOrNot(1)*BorderOrNot(2) == 0, A(indexF1(DIM+1:2*DIM)) = -2; A(indexF2(DIM+1:2*DIM)) = 2; end
        #         end
#     end
    
    A = (1.0 / (2 * h)) * sparse(INDEXI,INDEXJ,INDEXVAL,sizeA(1),sizeA(2))
    ## Old codes
# A = sparse(4*M*N, 2*M*N);
# A(1,1) = -2; A(1,3)= 2;
# A(3,1) = -2; A(3,2*(M+1)-1) = 2;
# A(2,2) = -2; A(2,4) = 2;
# A(4,2) = -2; A(4,2*(M+1)) = 2;
    
    # # Inside
# for tempx = 2:(M-1)
#     for tempy = 2:(N-1)
#         index = tempx+M*(tempy-1); # Find the point position
#         indexUp = tempx+M*tempy;
#         indexDown = tempx+M*(tempy-2);
#         indexLeft = tempx-1+M*(tempy-1);
#         indexRight = tempx+1+M*(tempy-1);
    
    #         indexF11row = 4*index-3;
#         indexF21row = 4*index-2;
#         indexF12row = 4*index-1;
#         indexF22row = 4*index;
    
    #         indexF11col1 = 2*(indexLeft)-1; indexF11col2 = 2*(indexRight)-1;
#         indexF12col1 = 2*(indexDown)-1; indexF12col2 = 2*(indexUp)-1;
#         indexF21col1 = 2*(indexLeft); indexF21col2 = 2*(indexRight);
#         indexF22col1 = 2*(indexDown); indexF22col2 = 2*(indexUp);
    
    #         A(indexF11row,indexF11col1) = -1; A(indexF11row,indexF11col2) = 1;
#         A(indexF21row,indexF21col1) = -1; A(indexF21row,indexF21col2) = 1;
#         A(indexF12row,indexF12col1) = -1; A(indexF12row,indexF12col2) = 1;
#         A(indexF22row,indexF22col1) = -1; A(indexF22row,indexF22col2) = 1;
#     end
# end
    
    # # Bottom line
# for tempx = 2:(M-1)
#     for tempy = 1
#         index = tempx+M*(tempy-1); # Find the point position
#         indexUp = tempx+M*tempy;
#         # indexDown = tempx+M*(tempy-2);
#         indexLeft = tempx-1+M*(tempy-1);
#         indexRight = tempx+1+M*(tempy-1);
    
    #         indexF11row = 4*index-3;
#         indexF21row = 4*index-2;
#         indexF12row = 4*index-1;
#         indexF22row = 4*index;
    
    #         indexF11col1 = 2*(indexLeft)-1; indexF11col2 = 2*(indexRight)-1;
#         indexF12col1 = 2*(index)-1; indexF12col2 = 2*(indexUp)-1;
#         indexF21col1 = 2*(indexLeft); indexF21col2 = 2*(indexRight);
#         indexF22col1 = 2*(index); indexF22col2 = 2*(indexUp);
    
    #         A(indexF11row,indexF11col1) = -1; A(indexF11row,indexF11col2) = 1;
#         A(indexF21row,indexF21col1) = -1; A(indexF21row,indexF21col2) = 1;
#         A(indexF12row,indexF12col1) = -2; A(indexF12row,indexF12col2) = 2;
#         A(indexF22row,indexF22col1) = -2; A(indexF22row,indexF22col2) = 2;
    
    #     end
# end
    
    # # Top
# for tempx = 2:(M-1)
#     for tempy = N
#         index = tempx+M*(tempy-1); # Find the point position
#         # indexUp = tempx+M*tempy;
#         indexDown = tempx+M*(tempy-2);
#         indexLeft = tempx-1+M*(tempy-1);
#         indexRight = tempx+1+M*(tempy-1);
    
    #         indexF11row = 4*index-3;
#         indexF21row = 4*index-2;
#         indexF12row = 4*index-1;
#         indexF22row = 4*index;
    
    #         indexF11col1 = 2*(indexLeft)-1; indexF11col2 = 2*(indexRight)-1;
#         indexF12col1 = 2*(indexDown)-1; indexF12col2 = 2*(index)-1;
#         indexF21col1 = 2*(indexLeft); indexF21col2 = 2*(indexRight);
#         indexF22col1 = 2*(indexDown); indexF22col2 = 2*(index);
    
    #         A(indexF11row,indexF11col1) = -1; A(indexF11row,indexF11col2) = 1;
#         A(indexF21row,indexF21col1) = -1; A(indexF21row,indexF21col2) = 1;
#         A(indexF12row,indexF12col1) = -2; A(indexF12row,indexF12col2) = 2;
#         A(indexF22row,indexF22col1) = -2; A(indexF22row,indexF22col2) = 2;
#     end
# end
    
    # # Left
# for tempx = 1
#     for tempy = 2:(N-1)
#         index = tempx+M*(tempy-1); # Find the point position
#         indexUp = tempx+M*tempy;
#         indexDown = tempx+M*(tempy-2);
#         # indexLeft = tempx-1+M*(tempy-1);
#         indexRight = tempx+1+M*(tempy-1);
    
    #         indexF11row = 4*index-3;
#         indexF21row = 4*index-2;
#         indexF12row = 4*index-1;
#         indexF22row = 4*index;
    
    #         indexF11col1 = 2*(index )-1; indexF11col2 = 2*(indexRight)-1;
#         indexF12col1 = 2*(indexDown)-1; indexF12col2 = 2*(indexUp)-1;
#         indexF21col1 = 2*(index ); indexF21col2 = 2*(indexRight);
#         indexF22col1 = 2*(indexDown); indexF22col2 = 2*(indexUp);
    
    #         A(indexF11row,indexF11col1) = -2; A(indexF11row,indexF11col2) = 2;
#         A(indexF21row,indexF21col1) = -2; A(indexF21row,indexF21col2) = 2;
#         A(indexF12row,indexF12col1) = -1; A(indexF12row,indexF12col2) = 1;
#         A(indexF22row,indexF22col1) = -1; A(indexF22row,indexF22col2) = 1;
#     end
# end
    
    # # Right
# for tempx = M
#     for tempy = 2:(N-1)
#         index = tempx+M*(tempy-1); # Find the point position
#         indexUp = tempx+M*tempy;
#         indexDown = tempx+M*(tempy-2);
#         indexLeft = tempx-1+M*(tempy-1);
#         # indexRight = tempx+1+M*(tempy-1);
    
    #         indexF11row = 4*index-3;
#         indexF21row = 4*index-2;
#         indexF12row = 4*index-1;
#         indexF22row = 4*index;
    
    #         indexF11col1 = 2*(indexLeft)-1; indexF11col2 = 2*(index )-1;
#         indexF12col1 = 2*(indexDown)-1; indexF12col2 = 2*(indexUp)-1;
#         indexF21col1 = 2*(indexLeft); indexF21col2 = 2*(index );
#         indexF22col1 = 2*(indexDown); indexF22col2 = 2*(indexUp);
    
    #         A(indexF11row,indexF11col1) = -2; A(indexF11row,indexF11col2) = 2;
#         A(indexF21row,indexF21col1) = -2; A(indexF21row,indexF21col2) = 2;
#         A(indexF12row,indexF12col1) = -1; A(indexF12row,indexF12col2) = 1;
#         A(indexF22row,indexF22col1) = -1; A(indexF22row,indexF22col2) = 1;
#     end
# end
    
    
    # # LeftBottom
# for tempx = 1
#     for tempy = 1
#         index = tempx+M*(tempy-1); # Find the point position
#         indexUp = tempx+M*tempy;
#         # indexDown = tempx+M*(tempy-2);
#         # indexLeft = tempx-1+M*(tempy-1);
#         indexRight = tempx+1+M*(tempy-1);
    
    #         indexF11row = 4*index-3;
#         indexF21row = 4*index-2;
#         indexF12row = 4*index-1;
#         indexF22row = 4*index;
    
    #         indexF11col1 = 2*(index )-1; indexF11col2 = 2*(indexRight)-1;
#         indexF12col1 = 2*(index )-1; indexF12col2 = 2*(indexUp)-1;
#         indexF21col1 = 2*(index ); indexF21col2 = 2*(indexRight);
#         indexF22col1 = 2*(index ); indexF22col2 = 2*(indexUp);
    
    #         A(indexF11row,indexF11col1) = -2; A(indexF11row,indexF11col2) = 2;
#         A(indexF21row,indexF21col1) = -2; A(indexF21row,indexF21col2) = 2;
#         A(indexF12row,indexF12col1) = -2; A(indexF12row,indexF12col2) = 2;
#         A(indexF22row,indexF22col1) = -2; A(indexF22row,indexF22col2) = 2;
#     end
# end
    
    # # RightBottom
# for tempx = M
#     for tempy = 1
#         index = tempx+M*(tempy-1); # Find the point position
#         indexUp = tempx+M*tempy;
#         #indexDown = tempx+M*(tempy-2);
#         indexLeft = tempx-1+M*(tempy-1);
#         #indexRight = tempx+1+M*(tempy-1);
    
    #         indexF11row = 4*index-3;
#         indexF21row = 4*index-2;
#         indexF12row = 4*index-1;
#         indexF22row = 4*index;
    
    #         indexF11col1 = 2*(indexLeft)-1; indexF11col2 = 2*(index )-1;
#         indexF12col1 = 2*(index )-1; indexF12col2 = 2*(indexUp)-1;
#         indexF21col1 = 2*(indexLeft); indexF21col2 = 2*(index );
#         indexF22col1 = 2*(index ); indexF22col2 = 2*(indexUp);
    
    #         A(indexF11row,indexF11col1) = -2; A(indexF11row,indexF11col2) = 2;
#         A(indexF21row,indexF21col1) = -2; A(indexF21row,indexF21col2) = 2;
#         A(indexF12row,indexF12col1) = -2; A(indexF12row,indexF12col2) = 2;
#         A(indexF22row,indexF22col1) = -2; A(indexF22row,indexF22col2) = 2;
#     end
# end
    
    # # LeftTop
# for tempx = 1
#     for tempy = N
#         index = tempx+M*(tempy-1); # Find the point position
#         #indexUp = tempx+M*tempy;
#         indexDown = tempx+M*(tempy-2);
#         #indexLeft = tempx-1+M*(tempy-1);
#         indexRight = tempx+1+M*(tempy-1);
    
    #         indexF11row = 4*index-3;
#         indexF21row = 4*index-2;
#         indexF12row = 4*index-1;
#         indexF22row = 4*index;
    
    #         indexF11col1 = 2*(index )-1; indexF11col2 = 2*(indexRight)-1;
#         indexF12col1 = 2*(indexDown)-1; indexF12col2 = 2*(index )-1;
#         indexF21col1 = 2*(index ); indexF21col2 = 2*(indexRight);
#         indexF22col1 = 2*(indexDown); indexF22col2 = 2*(index );
    
    #         A(indexF11row,indexF11col1) = -2; A(indexF11row,indexF11col2) = 2;
#         A(indexF21row,indexF21col1) = -2; A(indexF21row,indexF21col2) = 2;
#         A(indexF12row,indexF12col1) = -2; A(indexF12row,indexF12col2) = 2;
#         A(indexF22row,indexF22col1) = -2; A(indexF22row,indexF22col2) = 2;
#     end
# end
    
    # # RightTop
# for tempx = M
#     for tempy = N
#         index = tempx+M*(tempy-1); # Find the point position
#         #indexUp = tempx+M*tempy;
#         indexDown = tempx+M*(tempy-2);
#         indexLeft = tempx-1+M*(tempy-1);
#         #indexRight = tempx+1+M*(tempy-1);
    
    #         indexF11row = 4*index-3;
#         indexF21row = 4*index-2;
#         indexF12row = 4*index-1;
#         indexF22row = 4*index;
    
    #         indexF11col1 = 2*(indexLeft)-1; indexF11col2 = 2*(index)-1;
#         indexF12col1 = 2*(indexDown)-1; indexF12col2 = 2*(index)-1;
#         indexF21col1 = 2*(indexLeft); indexF21col2 = 2*(index);
#         indexF22col1 = 2*(indexDown); indexF22col2 = 2*(index);
    
    #         A(indexF11row,indexF11col1) = -2; A(indexF11row,indexF11col2) = 2;
#         A(indexF21row,indexF21col1) = -2; A(indexF21row,indexF21col2) = 2;
#         A(indexF12row,indexF12col1) = -2; A(indexF12row,indexF12col2) = 2;
#         A(indexF22row,indexF22col1) = -2; A(indexF22row,indexF22col2) = 2;
#     end
# end
    
    
    # A = A/(2*h);
    return A