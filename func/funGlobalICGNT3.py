############################################################
# function Triangulation FE-based Global DVC ICGN-code     #
# Object: to find deformation field using global methods   #
# Author: Jin Yang                                         #
# Last date modified: 2019.03                              #
############################################################

import numpy as np
import warnings
    
def funGlobalICGNT3(DICmesh = None,Df = None,Img1 = None,Img2 = None,U = None,alpha = None,tol = None,maxIter = None): 
    coordinatesFEM = DICmesh.coordinatesFEM
    elementsFEM = DICmesh.elementsFEM
    clusterNo = 8
    DIM = 2
    NodesPerEle = 3
    
    FEMSize = DIM * coordinatesFEM.shape[1-1]
    DfDx = Df.DfDx
    DfDy = Df.DfDy
    DfAxis = Df.DfAxis
    DfDxStartx = DfAxis(1)
    DfDxStarty = DfAxis(3)
    ImgPydUnit = 1
    dirichlet = []
    neumann = []
    ## ############################################################
    for stepwithinwhile in np.arange(1,maxIter+1).reshape(-1):
        tic
        print(np.array(['--- Global IC-GN iteration step',num2str(stepwithinwhile),' ---']))
        #     if clusterNo==0 || clusterNo==1
#     if (stepwithinwhile==1)
#         INDEXAI = []; INDEXAJ = []; INDEXAVAL = []; INDEXAREG = [];
#         #A = sparse(2*size(coordinatesFEM,1),2*size(coordinatesFEM,1));
#     end
#     INDEXBI = []; INDEXBVAL = []; #clear b; b = sparse(2*size(coordinatesFEM,1),1);
#     end
# ================= Navier-Lame elasticity regularization ============
# MatrixGrad = [-1 0 0 0 0 0 0 1; 0 1 0 0 0 0 0 -1; 0 -1 0 0 1 0 0 0; 0 0 1 0 -1 0 0 0;
#                0 0 1 0 0 -1 0 0; 0 0 0 -1 0 1 0 0; 0 0 0 1 0 0 -1 0; -1 0 0 0 0 0 1 0];
        # MatrixGradUpdate = MatrixGrad(:,1:4);
# MatrixGradUpdate(:,1) = MatrixGradUpdate(:,1) + 0.5*MatrixGrad(:,7) + 0.5*MatrixGrad(:,8);
# MatrixGradUpdate(:,2) = MatrixGradUpdate(:,2) + 0.5*MatrixGrad(:,8) + 0.5*MatrixGrad(:,5);
# MatrixGradUpdate(:,3) = MatrixGradUpdate(:,3) + 0.5*MatrixGrad(:,5) + 0.5*MatrixGrad(:,6);
# MatrixGradUpdate(:,4) = MatrixGradUpdate(:,4) + 0.5*MatrixGrad(:,6) + 0.5*MatrixGrad(:,7);
        # Matrix1 = MatrixGradUpdate'*MatrixGradUpdate;
# Matrix1Update = zeros(8,8);
# Matrix1Update(1:2:end,1:2:end) = Matrix1;
# Matrix1Update(2:2:end,2:2:end) = Matrix1;
        # Matrix2 = zeros(16,8);
# Matrix2(1:2:end,1:2:end) = MatrixGradUpdate;
# Matrix2(2:2:end,2:2:end) = MatrixGradUpdate;
# Matrix2Update = 0.25*Matrix2'*diag([1 0 1 0 0 1 0 1 1 0 1 0 0 1 0 1])*Matrix2;
        # # ------- Lame elasticity constants -------
# mu = alpha*1; lamda = mu;
# ====================================================================
        if clusterNo == 0 or clusterNo == 1:
            hbar = waitbar(0,np.array(['Global ICGN iteartion step: ',num2str(stepwithinwhile)]))
        else:
            hbar = parfor_progressbar(elementsFEM.shape[1-1],np.array(['Global ICGN iteartion step: ',num2str(stepwithinwhile)]))
        # ============= Each element, assemble stiffness matrix ============
        for j in np.arange(1,elementsFEM.shape[1-1]+1).reshape(-1):
            if clusterNo == 0 or clusterNo == 1:
                waitbar(j / elementsFEM.shape[1-1])
            else:
                hbar.iterate(1)
            tempA = np.zeros((DIM * NodesPerEle,DIM * NodesPerEle))
            tempb = tempA(:,1)
            # ------ Find three corner pts ------
            pt1x = coordinatesFEM(elementsFEM(j,1),1)
            pt1y = coordinatesFEM(elementsFEM(j,1),2)
            pt2x = coordinatesFEM(elementsFEM(j,2),1)
            pt2y = coordinatesFEM(elementsFEM(j,2),2)
            pt3x = coordinatesFEM(elementsFEM(j,3),1)
            pt3y = coordinatesFEM(elementsFEM(j,3),2)
            # ------ Compute triangle area --------
            TriArea = det(np.array([[1,pt1x,pt1y],[1,pt2x,pt2y],[1,pt3x,pt3y]]))
            # ------ Calculate DN Matrix for CST ------
            funDN1x = 1 / (2 * TriArea) * (pt2y - pt3y)
            funDN1y = 1 / (2 * TriArea) * (pt3x - pt2x)
            funDN2x = 1 / (2 * TriArea) * (pt3y - pt1y)
            funDN2y = 1 / (2 * TriArea) * (pt1x - pt3x)
            funDN3x = 1 / (2 * TriArea) * (pt1y - pt2y)
            funDN3y = 1 / (2 * TriArea) * (pt2x - pt1x)
            DN = np.array([[funDN1x,0,funDN2x,0,funDN3x,0],[funDN1y,0,funDN2y,0,funDN3y,0],[0,funDN1x,0,funDN2x,0,funDN3x],[0,funDN1y,0,funDN2y,0,funDN3y]])
            # ------ Find the element nodal indices ------
            tempIndexU = np.array([2 * elementsFEM(j,1) - 1,2 * elementsFEM(j,1),2 * elementsFEM(j,2) - 1,2 * elementsFEM(j,2),2 * elementsFEM(j,3) - 1,2 * elementsFEM(j,3)])
            # ------ Combine all the points -------
            ptOfxAll,ptOfyAll = ndgrid(np.arange(int(np.floor(np.amin(np.array([pt1x,pt2x,pt3x])))),np.ceil(np.amax(np.array([pt1x,pt2x,pt3x])))+1),np.arange(int(np.floor(np.amin(np.array([pt1y,pt2y,pt3y])))),np.ceil(np.amax(np.array([pt1y,pt2y,pt3y])))+1))
            # ------ Shape function N matrix ------
            NMat = cell(3,1)
            NMat[0] = (pt2x * pt3y + ptOfxAll * pt2y + ptOfyAll * pt3x - pt2x * ptOfyAll - pt2y * pt3x - pt3y * ptOfxAll) / TriArea
            NMat[2] = (pt3x * pt1y + ptOfxAll * pt3y + ptOfyAll * pt1x - pt3x * ptOfyAll - pt3y * pt1x - pt1y * ptOfxAll) / TriArea
            NMat[3] = (pt1x * pt2y + ptOfxAll * pt1y + ptOfyAll * pt2x - pt1x * ptOfyAll - pt1y * pt2x - pt2y * ptOfxAll) / TriArea
            tempUMat = np.zeros((ptOfxAll.shape,ptOfxAll.shape))
            tempVMat = tempUMat
            for tempk in np.arange(1,NodesPerEle+1).reshape(-1):
                tempUMat = tempUMat + np.multiply((U(tempIndexU(2 * tempk - 1)) * np.ones((ptOfxAll.shape,ptOfxAll.shape))),NMat[tempk])
                tempVMat = tempVMat + np.multiply((U(tempIndexU(2 * tempk - 0)) * np.ones((ptOfxAll.shape,ptOfxAll.shape))),NMat[tempk])
            tempg = ba_interp2(Img2,ptOfyAll + tempVMat,ptOfxAll + tempUMat,'cubic')
            ptOfxAll = ptOfxAll
            ptOfyAll = ptOfyAll
            ptInTriOrNotMat = funPtInTriCheck(np.array([[pt1x,pt1y],[pt2x,pt2y],[pt3x,pt3y]]),np.array([ptOfxAll,ptOfyAll]))
            # ====== Assemble stiffness matrix and force vector ======
#  for ptOfx = min([pt1x,pt2x,pt3x]):ImgPydUnit:max([pt1x,pt2x,pt3x])
#      for ptOfy = min([pt1y,pt2y,pt3y]):ImgPydUnit:max([pt1y,pt2y,pt3y])
            for tempjj in np.arange(1,len(ptOfxAll)+1).reshape(-1):
                # Judge pt is inside triangle or not
# ptInTriangleOrNot = ptInTriangleCheck(pt1x,pt1y,pt2x,pt2y,pt3x,pt3y,ptOfx,ptOfy);
                if ptInTriOrNotMat(tempjj) == 1:
                    ptOfx = ptOfxAll(tempjj)
                    ptOfy = ptOfyAll(tempjj)
                    # ------ Calculate N ------
#N1 = det([1 ptOfx ptOfy; 1 pt2x pt2y; 1 pt3x pt3y])/TriArea;
#N2 = det([1 ptOfx ptOfy; 1 pt3x pt3y; 1 pt1x pt1y])/TriArea;
#N3 = det([1 ptOfx ptOfy; 1 pt1x pt1y; 1 pt2x pt2y])/TriArea;
                    N1 = NMat[0](tempjj)
                    N2 = NMat[2](tempjj)
                    N3 = NMat[3](tempjj)
                    N = np.array([[N1,0,N2,0,N3,0],[0,N1,0,N2,0,N3]])
                    # ------ Here approximate Dg(x+u)=Df(x) ------
                    DfEle = np.array([[DfDx(ptOfx - DfDxStartx,ptOfy - DfDxStarty)],[DfDy(ptOfx - DfDxStartx,ptOfy - DfDxStarty)]])
                    # ------ Only assemble stiffness in the first step ------
                    if (stepwithinwhile == 1):
                        #A(temp,temp) =  A(temp,temp) + (N'*DfEle)*((N'*DfEle)') + alpha*(DN')*DN ;
                        tempA = tempA + (np.transpose(N) * DfEle) * (np.transpose((np.transpose(N) * DfEle))) + alpha * (np.transpose(DN)) * DN
                    # ------ Construct b vector ------
#temp1 = [ptOfx;ptOfy] + N*U(tempIndexU);
#temp2 = ((Img1(ptOfx,ptOfy) - fungInterpolation_g(temp1(1), temp1(2), Img2(floor(temp1(1))-1*ImgPydUnit:ImgPydUnit:floor(temp1(1))+2*ImgPydUnit, floor(temp1(2))-1*ImgPydUnit:ImgPydUnit:floor(temp1(2))+2*ImgPydUnit))) * (N'*DfEle));
                    temp2 = ((Img1(ptOfx,ptOfy) - tempg(tempjj)) * (np.transpose(N) * DfEle))
                    #b(temp) = b(temp) + tempb - (alpha*(DN')*DN)*U(tempIndexU);
                    tempb = tempb + temp2 - (alpha * (np.transpose(DN)) * DN) * U(tempIndexU)
                #end
#end
            # --- To store A_ele for each element ---
            if (stepwithinwhile == 1):
                #  A_ele(j,1:36)=reshape(A(temp,temp),1,36);
                IndexAXX,IndexAYY = ndgrid(tempIndexU,tempIndexU)
                # INDEXAI = [INDEXAI;IndexAXX(:)]; INDEXAJ = [INDEXAJ;IndexAYY(:)]; INDEXAVAL = [INDEXAVAL;tempA(:)]; #INDEXAREG = [INDEXAREG;tempAreg(:)];
                INDEXAIpar[j] = IndexAXX
                INDEXAJpar[j] = IndexAYY
                INDEXAVALpar[j] = tempA
            # INDEXBI = [INDEXBI;tempIndexU(:)]; INDEXBVAL = [INDEXBVAL;tempb(:)];
            INDEXBIpar[j] = tempIndexU
            INDEXBVALpar[j] = tempb
        close_(hbar)
        if (stepwithinwhile == 1):
            A = sparse(FEMSize,FEMSize)
            for eleInd in np.arange(1,elementsFEM.shape[1-1]+1).reshape(-1):
                A = A + sparse(INDEXAIpar[eleInd],INDEXAJpar[eleInd],INDEXAVALpar[eleInd],FEMSize,FEMSize)
            # A = sparse(INDEXAI,INDEXAJ,INDEXAVAL,FEMSize,FEMSize);
        #b = sparse(INDEXBI,ones(length(INDEXBI),1),INDEXBVAL,FEMSize,1);
        b = sparse(FEMSize,1)
        for eleInd in np.arange(1,elementsFEM.shape[1-1]+1).reshape(-1):
            b = b + sparse(double(INDEXBIpar[eleInd]),np.ones((len(INDEXBIpar[eleInd]),1)),INDEXBVALpar[eleInd],FEMSize,1)
        # ====== Find involved coordiantes index ======
        coordsIndexInvolved = unique(elementsFEM)
        if coordsIndexInvolved(1) == 0:
            UIndexInvolved = np.array([[coordsIndexInvolved(np.arange(2,end()+1))],[coordsIndexInvolved(np.arange(2,end()+1))]])
            # Not including the first 0-th entry
            for tempi in np.arange(1,(coordsIndexInvolved.shape[1-1] - 1)+1).reshape(-1):
                UIndexInvolved[np.arange[2 * tempi - 1,2 * tempi+1]] = np.array([[2 * coordsIndexInvolved(tempi + 1) - 1],[2 * coordsIndexInvolved(tempi + 1)]])
        else:
            UIndexInvolved = np.array([[2 * coordsIndexInvolved - 1],[2 * coordsIndexInvolved]])
        W = sparse(2 * coordinatesFEM.shape[1-1],1)
        W[2 * unique[dirichlet]] = 0
        W[2 * unique[dirichlet] - 1] = 0
        dirichlettemp = np.array([[2 * dirichlet],[2 * dirichlet - 1]])
        FreeNodes = setdiff(UIndexInvolved,unique(dirichlettemp))
        W[FreeNodes] = np.linalg.solve(A(FreeNodes,FreeNodes),b(FreeNodes))
        normW = norm(W) / np.sqrt(W.shape[1-1])
        normOfW[stepwithinwhile] = normW
        TimeICGN[stepwithinwhile] = toc
        toc
        U = reshape(U,len(U),1)
        W = reshape(W,len(W),1)
        if stepwithinwhile == 1:
            normWOld = normW * 10
        else:
            normWOld = normOfW(stepwithinwhile - 1)
        if (normW < tol) or ((normW / normWOld > 1 - tol) and (normW / normWOld < 1)):
            U = U + W
            break
        else:
            if (normW >= tol and normW < 1 / tol):
                U = U + W
            else:
                warnings.warn('Get diverged in Global_ICGN!!!')
                break
    
    return U,normOfW,TimeICGN
    
    return U,normOfW,TimeICGN