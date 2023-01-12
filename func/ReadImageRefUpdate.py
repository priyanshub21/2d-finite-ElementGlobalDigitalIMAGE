###################################################################
# Update reference image for incremental mode
###################################################################

import numpy as np
    
def ReadImageRefUpdate(file_name = None,ImgSeqNum = None,ULast = None,DICpara = None,DICmesh = None): 
    ImgSeqIncROIUpdateOrNot = DICpara.ImgSeqIncROIUpdateOrNot
    gridxyROIRange = DICpara.gridxyROIRange
    winstepsize = DICpara.winstepsize
    coordinatesFEM = DICmesh.coordinatesFEM
    elementsFEM = DICmesh.elementsFEM
    M = DICmesh.M
    N = DICmesh.N
    ## ------ Use fixed mesh -or- Manually update mesh ------
    if ImgSeqIncROIUpdateOrNot == 1:
        gridxyROIRange = funReadImageRefUpdate(file_name[1,ImgSeqNum - 1])
        gridxROIRangeNew = np.array([np.arange(np.ceil(gridxyROIRange.gridx(1)),int(np.floor(gridxyROIRange.gridx(2)))+winstepsize,winstepsize)])
        gridyROIRangeNew = np.array([np.arange(np.ceil(gridxyROIRange.gridy(1)),int(np.floor(gridxyROIRange.gridy(2)))+winstepsize,winstepsize)])
        x0Newtemp,y0Newtemp = ndgrid(gridxROIRangeNew,gridyROIRangeNew)
        x0Newtemp = np.transpose(x0Newtemp)
        y0Newtemp = np.transpose(y0Newtemp)
        DICmesh = MeshSetUp(x0Newtemp,y0Newtemp,DICpara)
        DICpara.gridxyROIRange = gridxyROIRange
        ## ------ Update ROI automatically ------
    else:
        if ImgSeqIncROIUpdateOrNot == 2:
            # Try to automatically update FEM mesh when updating reference img,
# but sometimes doesn't work very well and error could accumulate.
# ---------------------------------------------------------------
#ULast = ResultDisp{ImgSeqNum-2}.U;
            coordinatesFEMOld = coordinatesFEM
            elementsFEMOld = elementsFEM
            coordinatesFEM[:,1] = (coordinatesFEM(:,1) + ULast(np.arange(1,end()+2,2)))
            coordinatesFEM[:,2] = (coordinatesFEM(:,2) + ULast(np.arange(2,end()+2,2)))
            coordinatesFEMNewBottomInd = np.arange(1,M+1,1)
            coordinatesFEMNewTopInd = np.arange(M * N - (M - 1),M * N+1,1)
            coordinatesFEMNewLeftInd = np.arange(1,M * N+M,M)
            coordinatesFEMNewRightInd = np.arange(M,M * N+M,M)
            gridxyROIRange.gridxROIRange = np.array([np.amax(coordinatesFEM(coordinatesFEMNewLeftInd,1)),np.amin(coordinatesFEM(coordinatesFEMNewRightInd,1))])
            gridxyROIRange.gridyROIRange = np.array([np.amax(coordinatesFEM(coordinatesFEMNewBottomInd,2)),np.amin(coordinatesFEM(coordinatesFEMNewTopInd,2))])
            DICpara.gridxyROIRange = gridxyROIRange
            DICmesh.coordinatesFEM = coordinatesFEM
    
    return DICpara,DICmesh