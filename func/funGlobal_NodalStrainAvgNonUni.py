import numpy as np
    
def funGlobal_NodalStrainAvgNonUni(coordinatesFEM = None,elementsFEM = None,U = None,GaussPtOrder = None): 
    #FUNCTION [StrainNodalPt] = funGlobal_NodalStrainAvg(coordinatesFEM,elementsFEM,U,GaussPtOrder))
# Object: to compute strain fields by the FE-method
# ----------------------------------------------
    
    #	INPUT: coordinatesFEM      DIC FE Q4 mesh coordinates
#          elementsFEM         DIC FE Q4 mesh elements
#          U                   Disp vector: U = [Ux_node1, Uy_node1, ... , Ux_nodeN, Uy_nodeN]';
#          GaussPtOrder        Gauss point order used in FE Q4 element,  GaussPtOrder = 2 BY DEFAULT
    
    #   OUTPUT: StrainNodalPt      Solved strains at the nodal points,
    
    
    # ----------------------------------------------
# Reference
# [1] RegularizeNd. Matlab File Exchange open source.
# https://www.mathworks.com/matlabcentral/fileexchange/61436-regularizend
# [2] Gridfit. Matlab File Exchange open source.
# https://www.mathworks.com/matlabcentral/fileexchange/8998-surface-fitting-using-gridfit
# ----------------------------------------------
# Author: Jin Yang.
# Contact and support: jyang526@wisc.edu -or- aldicdvc@gmail.com
# Last time updated: 2018.03, 2020.12
# ==============================================
    
    ## Initialization
    FStrainAvgTimes = np.zeros((4 * coordinatesFEM.shape[1-1],1))
    FStrain = np.zeros((4 * coordinatesFEM.shape[1-1],1))
    ## ####################################
    for eleInd in np.arange(1,elementsFEM.shape[1-1]+1).reshape(-1):
        StrainWithinEachElementGausspoint = np.zeros((4,4))
        # ----- Find four corner points -------
        pt1x = coordinatesFEM(elementsFEM(eleInd,1),1)
        pt1y = coordinatesFEM(elementsFEM(eleInd,1),2)
        pt2x = coordinatesFEM(elementsFEM(eleInd,2),1)
        pt2y = coordinatesFEM(elementsFEM(eleInd,2),2)
        pt3x = coordinatesFEM(elementsFEM(eleInd,3),1)
        pt3y = coordinatesFEM(elementsFEM(eleInd,3),2)
        pt4x = coordinatesFEM(elementsFEM(eleInd,4),1)
        pt4y = coordinatesFEM(elementsFEM(eleInd,4),2)
        # ------ Find the element nodal indices ------
        tempIndexU = 2 * elementsFEM(eleInd,np.array([1,1,2,2,3,3,4,4]))
        tempIndexU[np.arange[1,end()+2,2]] = tempIndexU(np.arange(1,end()+2,2)) - 1
        # ------ Four Gauss points ------
        for tempi in np.arange(1,2+1).reshape(-1):
            for tempj in np.arange(1,2+1).reshape(-1):
                #ksi = ( 2*tempj-3)/sqrt(3); eta = (2*tempi-3)/sqrt(3);
                ksi = 2 * tempj - 3
                eta = 2 * tempi - 3
                if (tempi == 1) and (tempj == 1):
                    ksi = - 1 / np.sqrt(3)
                    eta = - 1 / np.sqrt(3)
                else:
                    if (tempi == 1) and (tempj == 2):
                        ksi = 1 / np.sqrt(3)
                        eta = - 1 / np.sqrt(3)
                    else:
                        if (tempi == 2) and (tempj == 1):
                            ksi = 1 / np.sqrt(3)
                            eta = 1 / np.sqrt(3)
                        else:
                            if (tempi == 2) and (tempj == 2):
                                ksi = - 1 / np.sqrt(3)
                                eta = 1 / np.sqrt(3)
                # ------ Calculate N ------
                N1 = (1 - ksi) * (1 - eta) * 0.25
                N2 = (1 + ksi) * (1 - eta) * 0.25
                N3 = (1 + ksi) * (1 + eta) * 0.25
                N4 = (1 - ksi) * (1 + eta) * 0.25
                N = np.array([[N1,0,N2,0,N3,0,N4,0],[0,N1,0,N2,0,N3,0,N4]])
                # ------ Build J matrix ------
                J11 = funDN1Dksi(ksi,eta) * pt1x + funDN2Dksi(ksi,eta) * pt2x + funDN3Dksi(ksi,eta) * pt3x + funDN4Dksi(ksi,eta) * pt4x
                J12 = funDN1Dksi(ksi,eta) * pt1y + funDN2Dksi(ksi,eta) * pt2y + funDN3Dksi(ksi,eta) * pt3y + funDN4Dksi(ksi,eta) * pt4y
                J21 = funDN1Deta(ksi,eta) * pt1x + funDN2Deta(ksi,eta) * pt2x + funDN3Deta(ksi,eta) * pt3x + funDN4Deta(ksi,eta) * pt4x
                J22 = funDN1Deta(ksi,eta) * pt1y + funDN2Deta(ksi,eta) * pt2y + funDN3Deta(ksi,eta) * pt3y + funDN4Deta(ksi,eta) * pt4y
                J = np.array([[J11,J12],[J21,J22]])
                Jacobian = det(J)
                InvJ = 1 / Jacobian * np.array([[J22,- J12],[- J21,J11]])
                # ------ Compute DN matrix ------
                DN = np.array([[InvJ,np.zeros((2,2))],[np.zeros((2,2)),InvJ]]) * np.array([[funDN1Dksi(ksi,eta),0,funDN2Dksi(ksi,eta),0,funDN3Dksi(ksi,eta),0,funDN4Dksi(ksi,eta),0],[funDN1Deta(ksi,eta),0,funDN2Deta(ksi,eta),0,funDN3Deta(ksi,eta),0,funDN4Deta(ksi,eta),0],[0,funDN1Dksi(ksi,eta),0,funDN2Dksi(ksi,eta),0,funDN3Dksi(ksi,eta),0,funDN4Dksi(ksi,eta)],[0,funDN1Deta(ksi,eta),0,funDN2Deta(ksi,eta),0,funDN3Deta(ksi,eta),0,funDN4Deta(ksi,eta)]])
                StrainWithinEachElementGausspoint[2 * [tempi - 1] + tempj,np.arange[1,4+1]] = DN * U(tempIndexU)
        # ------ Extrapolation matrix ------
        MatrixExtrapolation = np.array([[1 + 0.5 * np.sqrt(3),- 0.5,1 - 0.5 * np.sqrt(3),- 0.5],[- 0.5,1 + 0.5 * np.sqrt(3),- 0.5,1 - 0.5 * np.sqrt(3)],[1 - 0.5 * np.sqrt(3),- 0.5,1 + 0.5 * np.sqrt(3),- 0.5],[- 0.5,1 - 0.5 * np.sqrt(3),- 0.5,1 + 0.5 * np.sqrt(3)]])
        # ------ Nodal points strain extrapolation using Gauss points -----
        StrainExxWithinEachElementNodalpoint = MatrixExtrapolation * StrainWithinEachElementGausspoint(np.arange(1,4+1),1)
        StrainExyWithinEachElementNodalpoint = MatrixExtrapolation * StrainWithinEachElementGausspoint(np.arange(1,4+1),2)
        StrainEyxWithinEachElementNodalpoint = MatrixExtrapolation * StrainWithinEachElementGausspoint(np.arange(1,4+1),3)
        StrainEyyWithinEachElementNodalpoint = MatrixExtrapolation * StrainWithinEachElementGausspoint(np.arange(1,4+1),4)
        StrainWithinEachElementGausspoint[np.arange[1,4+1],1] = reshape(StrainExxWithinEachElementNodalpoint(np.arange(1,4+1)),4,1)
        StrainWithinEachElementGausspoint[np.arange[1,4+1],2] = reshape(StrainExyWithinEachElementNodalpoint(np.arange(1,4+1)),4,1)
        StrainWithinEachElementGausspoint[np.arange[1,4+1],3] = reshape(StrainEyxWithinEachElementNodalpoint(np.arange(1,4+1)),4,1)
        StrainWithinEachElementGausspoint[np.arange[1,4+1],4] = reshape(StrainEyyWithinEachElementNodalpoint(np.arange(1,4+1)),4,1)
        # ------ Find the element nodal indices for strain ------
        tempStrainIndex = np.array([4 * elementsFEM(eleInd,1) - 3,4 * elementsFEM(eleInd,1) - 2,4 * elementsFEM(eleInd,1) - 1,4 * elementsFEM(eleInd,1),4 * elementsFEM(eleInd,2) - 3,4 * elementsFEM(eleInd,2) - 2,4 * elementsFEM(eleInd,2) - 1,4 * elementsFEM(eleInd,2),4 * elementsFEM(eleInd,3) - 3,4 * elementsFEM(eleInd,3) - 2,4 * elementsFEM(eleInd,3) - 1,4 * elementsFEM(eleInd,3),4 * elementsFEM(eleInd,4) - 3,4 * elementsFEM(eleInd,4) - 2,4 * elementsFEM(eleInd,4) - 1,4 * elementsFEM(eleInd,4)])
        FStrain[tempStrainIndex] = FStrain(tempStrainIndex) + np.array([[StrainWithinEachElementGausspoint(1,1)],[StrainWithinEachElementGausspoint(1,2)],[StrainWithinEachElementGausspoint(1,3)],[StrainWithinEachElementGausspoint(1,4)],[StrainWithinEachElementGausspoint(2,1)],[StrainWithinEachElementGausspoint(2,2)],[StrainWithinEachElementGausspoint(2,3)],[StrainWithinEachElementGausspoint(2,4)],[StrainWithinEachElementGausspoint(3,1)],[StrainWithinEachElementGausspoint(3,2)],[StrainWithinEachElementGausspoint(3,3)],[StrainWithinEachElementGausspoint(3,4)],[StrainWithinEachElementGausspoint(4,1)],[StrainWithinEachElementGausspoint(4,2)],[StrainWithinEachElementGausspoint(4,3)],[StrainWithinEachElementGausspoint(4,4)]])
        FStrainAvgTimes[tempStrainIndex] = FStrainAvgTimes(tempStrainIndex) + np.ones((16,1))
    
    StrainNodalPt = FStrain / FStrainAvgTimes
    return StrainNodalPt
    
    ## ========= subroutines for FEM Q4 shape function derivatives ========
    
def funDN1Dksi(ksi = None,eta = None): 
    DN1Dksi = - (1 - eta) / 4
    return DN1Dksi
    
    
def funDN1Deta(ksi = None,eta = None): 
    DN1Deta = - (1 - ksi) / 4
    return DN1Deta
    
    
def funDN2Dksi(ksi = None,eta = None): 
    DN2Dksi = (1 - eta) / 4
    return DN2Dksi
    
    
def funDN2Deta(ksi = None,eta = None): 
    DN2Deta = - (1 + ksi) / 4
    return DN2Deta
    
    
def funDN3Dksi(ksi = None,eta = None): 
    DN3Dksi = (1 + eta) / 4
    return DN3Dksi
    
    
def funDN3Deta(ksi = None,eta = None): 
    DN3Deta = (1 + ksi) / 4
    return DN3Deta
    
    
def funDN4Dksi(ksi = None,eta = None): 
    DN4Dksi = - (1 + eta) / 4
    return DN4Dksi
    
    
def funDN4Deta(ksi = None,eta = None): 
    DN4Deta = (1 - ksi) / 4
    return DN4Deta
    
    return StrainNodalPt