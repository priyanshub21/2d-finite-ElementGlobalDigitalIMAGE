# =========================================================
# function GlobalICGN to solve FE iterations for global DIC
# ---------------------------------------------------------
#   INPUT:
#       DIC mesh, DIC image pair,
#       displacement initial guess, regularizer coefficient
#
#   OUTPUT:
#       U: Solved displacement field;
#       normOfW: FE-based global DIC iteration update norm;
#       timeICGN: Time cost for each FE-based global DIC iteration;
#
# Author: Jin Yang, jyang526@wisc.edu or aldicdvc@gmail.com
# Date: 2020.10
# =========================================================

import numpy as np
import warnings
    
def funGlobalICGNQ4(DICmesh = None,Df = None,Img1 = None,Img2 = None,U = None,alpha = None,tol = None,maxIter = None): 
    coordinatesFEM = DICmesh.coordinatesFEM
    indCoordinatesFEMNotZero,__ = find(coordinatesFEM(:,1) > 0)
    indCoordinatesFEMNotZero = unique(np.array([[2 * indCoordinatesFEMNotZero - 1],[2 * indCoordinatesFEMNotZero]]))
    elementsFEM = DICmesh.elementsFEM
    NodesPerEle = 4
    
    DIM = 2
    
    # winsize = (coordinatesFEM(2,1)-coordinatesFEM(1,1))*ones(1,DIM); # or: DICpara.winsize;
    FEMSize = DIM * coordinatesFEM.shape[1-1]
    
    DfDx = Df.DfDx
    DfDy = Df.DfDy
    DfAxis = Df.DfAxis
    DfDxStartx = DfAxis(1)
    DfDxStarty = DfAxis(3)
    try:
        maxIter = maxIter
    finally:
        pass
    
    ## ==============================================================
# Optional: construct Navier-Lame elasticity regularizer
# Ignore this section if you don't apply elasticity regularization
    
    # MatrixGrad = [-1 0 0 0 0 0 0 1; 0 1 0 0 0 0 0 -1; 0 -1 0 0 1 0 0 0; 0 0 1 0 -1 0 0 0;
#                0 0 1 0 0 -1 0 0; 0 0 0 -1 0 1 0 0; 0 0 0 1 0 0 -1 0; -1 0 0 0 0 0 1 0];
# MatrixGradUpdate = MatrixGrad(:,1:4);
# MatrixGradUpdate(:,1) = MatrixGradUpdate(:,1) + 0.5*MatrixGrad(:,7) + 0.5*MatrixGrad(:,8);
# MatrixGradUpdate(:,2) = MatrixGradUpdate(:,2) + 0.5*MatrixGrad(:,8) + 0.5*MatrixGrad(:,5);
# MatrixGradUpdate(:,3) = MatrixGradUpdate(:,3) + 0.5*MatrixGrad(:,5) + 0.5*MatrixGrad(:,6);
# MatrixGradUpdate(:,4) = MatrixGradUpdate(:,4) + 0.5*MatrixGrad(:,6) + 0.5*MatrixGrad(:,7);
    
    # Matrix1 = MatrixGradUpdate'*MatrixGradUpdate;
# Matrix1Update = zeros(8,8); Matrix1Update(1:2:end,1:2:end)=Matrix1; Matrix1Update(2:2:end,2:2:end)=Matrix1;
    
    # # Lame elasticity constants
# mu = alpha*1; lamda = mu;
    
    ## ############## Start FE ICGN iteration ###############
    for stepwithinwhile in np.arange(1,maxIter+1).reshape(-1):
        tic
        # ====== Initialize stiffness matrix at the first iteration ======
        if (stepwithinwhile == 1):
            print(np.array(['--- Global IC-GN iterations ---']))
            INDEXAI = []
            INDEXAJ = []
            INDEXAVAL = []
            INDEXAREG = []
        INDEXBI = []
        INDEXBVAL = []
        hbar = waitbar(0,np.array(['Global ICGN iteration step ',num2str(stepwithinwhile)]))
        # ============= Each element, assemble stiffness matrix ============
        for indEle in np.arange(1,elementsFEM.shape[1-1]+1).reshape(-1):
            waitbar(indEle / elementsFEM.shape[1-1])
            tempA = np.zeros((DIM * NodesPerEle,DIM * NodesPerEle))
            tempb = tempA(:,1)
            # ------ Find four corner points in Q4 FE element ------
            pt1x = coordinatesFEM(elementsFEM(indEle,1),1)
            pt1y = coordinatesFEM(elementsFEM(indEle,1),2)
            pt2x = coordinatesFEM(elementsFEM(indEle,2),1)
            pt2y = coordinatesFEM(elementsFEM(indEle,2),2)
            pt3x = coordinatesFEM(elementsFEM(indEle,3),1)
            pt3y = coordinatesFEM(elementsFEM(indEle,3),2)
            pt4x = coordinatesFEM(elementsFEM(indEle,4),1)
            pt4y = coordinatesFEM(elementsFEM(indEle,4),2)
            lMatrix = np.array([[pt1x * pt1y,pt1x,pt1y,1],[pt2x * pt2y,pt2x,pt2y,1],[pt3x * pt3y,pt3x,pt3y,1],[pt4x * pt4y,pt4x,pt4y,1]])
            # ------ Find linear interpolation coefficients ------
            lb = np.array([[- 1],[1],[1],[- 1]])
            l = linsolve(lMatrix,lb)
            mb = np.array([[- 1],[- 1],[1],[1]])
            m = linsolve(lMatrix,mb)
            # ------ Find element nodal indices ------
            tp = np.ones((1,DIM))
            tempIndexU = DIM * elementsFEM(indEle,np.array([tp,2 * tp,3 * tp,4 * tp]))
            for tempDIM in np.arange(1,DIM - 1+1).reshape(-1):
                tempIndexU[np.arange[tempDIM,end()+DIM,DIM]] = tempIndexU(np.arange(tempDIM,end()+DIM,DIM)) - (DIM - tempDIM)
            # or using the following lines
# tempIndexU = [2*elementsFEM(indEle,1)-1 2*elementsFEM(indEle,1) 2*elementsFEM(indEle,2)-1 2*elementsFEM(indEle,2)...
#               2*elementsFEM(indEle,3)-1 2*elementsFEM(indEle,3) 2*elementsFEM(indEle,4)-1 2*elementsFEM(indEle,4)];
            ##### Previous rectangular element #####
# [ptOfxAll, ptOfyAll] = ndgrid(pt1x:pt3x, pt1y:pt3y); # To compute at each pixels
            temp1 = elementsFEM(indEle,np.arange(1,4+1))
            k = convhull(coordinatesFEM(temp1,np.arange(1,2+1)))
            ptOfxAll_min = int(np.floor(np.amin(coordinatesFEM(elementsFEM(indEle,np.arange(1,4+1)),1))))
            ptOfxAll_max = np.ceil(np.amax(coordinatesFEM(elementsFEM(indEle,np.arange(1,4+1)),1)))
            ptOfyAll_min = int(np.floor(np.amin(coordinatesFEM(elementsFEM(indEle,np.arange(1,4+1)),2))))
            ptOfyAll_max = np.ceil(np.amax(coordinatesFEM(elementsFEM(indEle,np.arange(1,4+1)),2)))
            ptOfxAll,ptOfyAll = ndgrid(np.arange(ptOfxAll_min,ptOfxAll_max+1),np.arange(ptOfyAll_min,ptOfyAll_max+1))
            in_ = inhull(np.array([ptOfxAll,ptOfyAll]),coordinatesFEM(temp1(k),np.arange(1,2+1)),[],0.1)
            indPtAll,__ = find(in_ == 1)
            ###### Check inhull #####
# figure, plot(ptOfxAll(:),ptOfyAll(:),'b.');
# hold on; plot(ptOfxAll(indPtAll),ptOfyAll(indPtAll),'r+');
# hold on; plot(coordinatesFEM(temp1(k),1),coordinatesFEM(temp1(k),2),'k');
            for indPt in np.arange(1,len(indPtAll)+1).reshape(-1):
                ptOfx = ptOfxAll(indPtAll(indPt))
                ptOfy = ptOfyAll(indPtAll(indPt))
                # ------ Calculate ksi and eta ------
                ksi = l(1) * ptOfx * ptOfy + l(2) * ptOfx + l(3) * ptOfy + l(4)
                eta = m(1) * ptOfx * ptOfy + m(2) * ptOfx + m(3) * ptOfy + m(4)
                # ------ Calculate N matrix ------
                N1 = (1 - ksi) * (1 - eta) * 0.25
                N2 = (1 + ksi) * (1 - eta) * 0.25
                N3 = (1 + ksi) * (1 + eta) * 0.25
                N4 = (1 - ksi) * (1 + eta) * 0.25
                # ------ Generate [N] shape function matrix ------
                N = np.array([[N1,0,N2,0,N3,0,N4,0],[0,N1,0,N2,0,N3,0,N4]])
                # ------ Build J matrix ------
                J = np.array([[funDN1Dksi(ksi,eta),funDN2Dksi(ksi,eta),funDN3Dksi(ksi,eta),funDN4Dksi(ksi,eta)],[funDN1Deta(ksi,eta),funDN2Deta(ksi,eta),funDN3Deta(ksi,eta),funDN4Deta(ksi,eta)]]) * np.array([[pt1x,pt1y],[pt2x,pt2y],[pt3x,pt3y],[pt4x,pt4y]])
                # J11 = funDN1Dksi(ksi,eta)*pt1x + funDN2Dksi(ksi,eta)*pt2x + ...
#       funDN3Dksi(ksi,eta)*pt3x + funDN4Dksi(ksi,eta)*pt4x;
# J12 = funDN1Dksi(ksi,eta)*pt1y + funDN2Dksi(ksi,eta)*pt2y + ...
#       funDN3Dksi(ksi,eta)*pt3y + funDN4Dksi(ksi,eta)*pt4y;
# J21 = funDN1Deta(ksi,eta)*pt1x + funDN2Deta(ksi,eta)*pt2x + ...
#       funDN3Deta(ksi,eta)*pt3x + funDN4Deta(ksi,eta)*pt4x;
# J22 = funDN1Deta(ksi,eta)*pt1y + funDN2Deta(ksi,eta)*pt2y + ...
#       funDN3Deta(ksi,eta)*pt3y + funDN4Deta(ksi,eta)*pt4y;
# J = [J11 J12; J21 J22];
                Jacobian = det(J)
                InvJ = 1 / Jacobian * np.array([[J(2,2),- J(1,2)],[- J(2,1),J(1,1)]])
                # ------ Compute DN matrix ------
                DN = np.array([[InvJ,np.zeros((2,2))],[np.zeros((2,2)),InvJ]]) * np.array([[funDN1Dksi(ksi,eta),0,funDN2Dksi(ksi,eta),0,funDN3Dksi(ksi,eta),0,funDN4Dksi(ksi,eta),0],[funDN1Deta(ksi,eta),0,funDN2Deta(ksi,eta),0,funDN3Deta(ksi,eta),0,funDN4Deta(ksi,eta),0],[0,funDN1Dksi(ksi,eta),0,funDN2Dksi(ksi,eta),0,funDN3Dksi(ksi,eta),0,funDN4Dksi(ksi,eta)],[0,funDN1Deta(ksi,eta),0,funDN2Deta(ksi,eta),0,funDN3Deta(ksi,eta),0,funDN4Deta(ksi,eta)]])
                # ------ Here approximate Dg(x+u)=Df(x) ------
                DfEle = np.array([[DfDx(ptOfx - DfDxStartx,ptOfy - DfDxStarty)],[DfDy(ptOfx - DfDxStartx,ptOfy - DfDxStarty)]])
                # ------ Only assemble stiffness in the first step ------
                if (stepwithinwhile == 1):
                    #A(tempIndexU,tempIndexU) =  A(tempIndexU,tempIndexU) + (N'*Df)*(N'*Df)' + alpha*(DN')*DN ;
                    tempA = tempA + (np.transpose(N) * DfEle) * np.transpose((np.transpose(N) * DfEle)) + alpha * (np.transpose(DN)) * DN
                # ------ Construct b vector ------
                temp1 = np.array([[ptOfx],[ptOfy]]) + N * U(tempIndexU)
                temp2 = ((Img1(ptOfx,ptOfy) - fungInterpolation_g(temp1(1),temp1(2),Img2(np.arange(int(np.floor(temp1(1))) - 1,int(np.floor(temp1(1))) + 2+1),np.arange(int(np.floor(temp1(2))) - 1,int(np.floor(temp1(2))) + 2+1)))) * (np.transpose(N) * DfEle))
                tempb = tempb + temp2 - (alpha * (np.transpose(DN)) * DN) * U(tempIndexU)
                # ====== Optional ======
# --- To use Navier-Lame elasticity regularization instead ---
# b(temp) = b(temp) + tempb + (mu*Matrix1Update + (mu+lamda)*Matrix2Update)*U(temp);
# ------------------------------------------------------------
                # end # for ptOfx = pt1x:pt3x
# end # for ptOfy = pt1y:pt3y
            # --- To store A_ele for each element ---
            if (stepwithinwhile == 1):
                IndexAXX,IndexAYY = ndgrid(tempIndexU,tempIndexU)
                INDEXAI = np.array([[INDEXAI],[IndexAXX]])
                INDEXAJ = np.array([[INDEXAJ],[IndexAYY]])
                INDEXAVAL = np.array([[INDEXAVAL],[tempA]])
            INDEXBI = np.array([[INDEXBI],[tempIndexU]])
            INDEXBVAL = np.array([[INDEXBVAL],[tempb]])
        close_(hbar)
        if (stepwithinwhile == 1):
            A = sparse(INDEXAI,INDEXAJ,INDEXAVAL,FEMSize,FEMSize)
        b = sparse(INDEXBI,np.ones((len(INDEXBI),1)),INDEXBVAL,FEMSize,1)
        # ========= Solve FEM problem ===========
        W = np.linalg.solve(A(indCoordinatesFEMNotZero,indCoordinatesFEMNotZero),b(indCoordinatesFEMNotZero))
        normW = norm(W) / np.sqrt(W.shape[1-1])
        normOfW[stepwithinwhile] = normW
        timeICGN[stepwithinwhile] = toc
        U = reshape(U,len(U),1)
        W = reshape(W,len(W),1)
        print(np.array(['normW = ',num2str(normW),' at iter ',num2str(stepwithinwhile),'; time cost = ',num2str(toc),'s']))
        if stepwithinwhile == 1:
            normWOld = normW * 10
        else:
            normWOld = normOfW(stepwithinwhile - 1)
        if (normW < tol):
            U[indCoordinatesFEMNotZero] = U(indCoordinatesFEMNotZero) + W
            break
        else:
            if (normW >= tol and normW < (0.1 / tol)):
                U[indCoordinatesFEMNotZero] = U(indCoordinatesFEMNotZero) + W
            else:
                warnings.warn('Get diverged in Global_ICGN!!!')
                break
    
    TotalTimeICGN = sum(timeICGN)
    print(np.array(['Elapsed time is ',num2str(TotalTimeICGN),' seconds.']))
    return U,normOfW,timeICGN
    
    ## ========= subroutines for  FEM Q4 shape function derivatives ========
    
def funDN1Dksi(ksi = None,eta = None): 
    DN1x = - (1 - eta) / 4
    return DN1x
    
    
def funDN1Deta(ksi = None,eta = None): 
    DN1y = - (1 - ksi) / 4
    return DN1y
    
    
def funDN2Dksi(ksi = None,eta = None): 
    DN2x = (1 - eta) / 4
    return DN2x
    
    
def funDN2Deta(ksi = None,eta = None): 
    DN2y = - (1 + ksi) / 4
    return DN2y
    
    
def funDN3Dksi(ksi = None,eta = None): 
    DN3x = (1 + eta) / 4
    return DN3x
    
    
def funDN3Deta(ksi = None,eta = None): 
    DN3y = (1 + ksi) / 4
    return DN3y
    
    
def funDN4Dksi(ksi = None,eta = None): 
    DN4x = - (1 + eta) / 4
    return DN4x
    
    
def funDN4Deta(ksi = None,eta = None): 
    DN4y = (1 - ksi) / 4
    return DN4y
    
    return U,normOfW,timeICGN