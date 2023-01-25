# ---------------------------------------------
# Finite-element-based global DIC (FE-Global-DIC)
# Author: Jin Yang, Postdoc @UW-Madison;  PhD @Caltech 19';
# Contact: aldicdvc@gmail.com; jyang526@wisc.edu
# 2015.04,06,07; 2016.03,04; 2020.11
# ---------------------------------------------

## Section 1: Clear MATLAB environment & mex set up Spline interpolation
import numpy as np
import os
close_('all')
clear
clearvars
print('------------ Section 1 Start ------------ \n' % ())
setenv('MW_MINGW64_LOC','C:\TDM-GCC-64')
# # cd("./Splines_interp/lib_matlab"); CompileLib; cd("../../");  # # mex bi-cubic spline interpolations
# # addpath("./Splines_interp/lib_matlab"); # dbstop if error # # Old version codes.

addpath('./func','./src','./plotFiles/','./plotFiles/export_fig-d966721/')
# addpath("./YOUR IMAGE FOLDER");
print('------------ Section 1 Done ------------ \n \n' % ())
## Section 2: Load DIC parameters and set up DIC parameters
print('------------ Section 2 Start ------------ \n' % ())
# ====== Read images ======
file_name,Img,DICpara = ReadImageQ4
close_('all')
# ====== Load mat file for your mesh coordinates and elements ======
# User needs to modify these lines to upload his/her own FE-mesh
scipy.io.loadmat('./Images_Sample12/mesh_plate_hole_Q4.mat')
DICmesh.coordinatesFEM = coordinatesFEM
DICmesh.coordinatesFEMWorld = np.array([coordinatesFEM(:,1),Img[0].shape[2-1] + 1 - coordinatesFEM(:,2)])
DICmesh.elementsFEM = elementsFEM
# ====== Normalize images ======
ImgNormalized,DICpara.gridxyROIRange = funNormalizeImg(Img,DICpara.gridxyROIRange)
fNormalized = ImgNormalized[0]

# ====== Initialize variable storage ======
ResultDisp = cell(len(ImgNormalized) - 1,1)
ResultDefGrad = cell(len(ImgNormalized) - 1,1)
ResultStrainWorld = cell(len(ImgNormalized) - 1,1)
ResultStressWorld = cell(len(ImgNormalized) - 1,1)
ResultFEMesh = cell(np.ceil((len(ImgNormalized) - 1) / DICpara.ImgSeqIncUnit),1)

ResultFEMeshEachFrame = cell(len(ImgNormalized) - 1,1)

ResultAlpha = cell(len(ImgNormalized) - 1,1)
ResultNormOfW = cell(len(ImgNormalized) - 1,1)
ResultTimeICGN = cell(len(ImgNormalized) - 1,1)
print('------------ Section 2 Done ------------ \n \n' % ())
## ########################################################################
# Start each frame in an image sequence
# ##########################################################################
for ImgSeqNum in np.arange(2,len(ImgNormalized)+1).reshape(-1):
    print(np.array(['Current image frame #: ',num2str(ImgSeqNum),'/',num2str(len(ImgNormalized))]))
    gNormalized = ImgNormalized[ImgSeqNum]
    if ImgSeqNum < 7:
        DICpara.NewFFTSearch = 1
    else:
        DICpara.NewFFTSearch = 0
    ## Section 3: Compute an initial guess of the unknown displacement field
    print('\n' % ())
    print('------------ Section 3 Start ------------ \n' % ())
    # ####################################################################
# This section is to find or update an initial guess of the unknown displacements.
# The key idea is to either to use a new FFT-based cross correlation peak fitting,
# or use the results from the last frame as the new initial guess for the next frame;
# Particularly in the incremental mode DIC, the reference image can also be updated, e.g.,
# " fNormalized = ImgNormalized{ImgSeqNum-mod(ImgSeqNum-1,ImgSeqIncUnit)}; "
    # DICpara.NewFFTSearch = 0; # If you want to apply the FFT-based cross correlation to
# compute the initial guess for each frame, please make sure that "DICpara.NewFFTSearch = 0".
# ####################################################################
    #####################################################################
    if ImgSeqNum == 2 or DICpara.NewFFTSearch == 1:
        # ====== Integer Search ======
        uSeedPt,vSeedPt,PhiSeedPt,tempSizeOfSearchRegion = funIntegerSearchPt(fNormalized,gNormalized,DICpara.winstepsize,DICmesh.coordinatesFEM)
        ##### Remove nan bad points #####
        row,col = find(np.isnan(uSeedPt))
        rowNotNan = setdiff(np.transpose(np.array([np.arange(1,DICmesh.coordinatesFEM.shape[1-1]+1,1)])),row)
        F_u_interp = scatteredInterpolant(DICmesh.coordinatesFEM(rowNotNan,:),uSeedPt(rowNotNan),'linear','linear')
        uSeedPt = F_u_interp(DICmesh.coordinatesFEM)
        F_v_interp = scatteredInterpolant(DICmesh.coordinatesFEM(rowNotNan,:),vSeedPt(rowNotNan),'linear','linear')
        vSeedPt = F_v_interp(DICmesh.coordinatesFEM)
        U0 = np.transpose(np.array([uSeedPt,vSeedPt]))
        U0 = U0
        DICpara.tempSizeOfSearchRegion = tempSizeOfSearchRegion
        ##### Plot initial results #####
        Plotdisp_show(U0,DICmesh.coordinatesFEM,DICmesh.elementsFEM)
        # ====== Deal with incremental mode ======
        fNormalizedNewIndex = ImgSeqNum - np.mod(ImgSeqNum - 2,DICpara.ImgSeqIncUnit) - 1
        if DICpara.ImgSeqIncUnit == 1:
            fNormalizedNewIndex = fNormalizedNewIndex - 1
        ResultFEMesh[1 + int[np.floor[fNormalizedNewIndex / DICpara.ImgSeqIncUnit]]] = struct('coordinatesFEM',DICmesh.coordinatesFEM,'elementsFEM',DICmesh.elementsFEM,'winsize',DICpara.winsize,'winstepsize',DICpara.winstepsize,'gridxyROIRange',DICpara.gridxyROIRange)
        #####################################################################
    else:
        if np.mod(ImgSeqNum - 2,DICpara.ImgSeqIncUnit) == 0:
            fNormalizedNewIndex = ImgSeqNum - np.mod(ImgSeqNum - 2,DICpara.ImgSeqIncUnit) - 1
            if DICpara.ImgSeqIncUnit == 1:
                fNormalizedNewIndex = fNormalizedNewIndex - 1
            fNormalized = ImgNormalized[fNormalizedNewIndex]
            DICpara,DICmesh = ReadImageRefUpdate(file_name,ImgSeqNum,ResultDisp[ImgSeqNum - 2].U,DICpara,DICmesh)
            U0 = np.zeros((2 * DICmesh.coordinatesFEM.shape[1-1],1))
            ResultFEMesh[1 + int[np.floor[fNormalizedNewIndex / DICpara.ImgSeqIncUnit]]] = struct('coordinatesFEM',DICmesh.coordinatesFEM,'elementsFEM',DICmesh.elementsFEM,'winsize',DICpara.winsize,'winstepsize',DICpara.winstepsize,'gridxyROIRange',DICpara.gridxyROIRange)
            #####################################################################
        else:
            if ImgSeqNum < 7:
                U0 = ResultDisp[ImgSeqNum - 2].U
            else:
                nTime = 5
                np = len(ResultDisp[ImgSeqNum - 2].U) / 2
                T_data_u = np.zeros((nTime,np))
                T_data_v = np.zeros((nTime,np))
                for tempi in np.arange(1,nTime+1).reshape(-1):
                    T_data_u[tempi,:] = np.transpose(ResultDisp[ImgSeqNum - (2 + nTime) + tempi,1].U(np.arange(1,np * 2+2,2)))
                    T_data_v[tempi,:] = np.transpose(ResultDisp[ImgSeqNum - (2 + nTime) + tempi,1].U(np.arange(2,np * 2+2,2)))
                nB = 3
                t_train = np.transpose(np.array([np.arange(ImgSeqNum - 1 - nTime,ImgSeqNum - 2+1)]))
                t_pre = np.transpose(np.array([ImgSeqNum - 1]))
                u_pred,__,__,__ = funPOR_GPR(T_data_u,t_train,t_pre,nB)
                v_pred,__,__,__ = funPOR_GPR(T_data_v,t_train,t_pre,nB)
                tempu = u_pred(1,:)
                tempv = v_pred(1,:)
                U0 = np.transpose(np.array([tempu,tempv]))
                U0 = U0
                # ##### After running the new ImgSeqNum, you can uncomment these
# ##### lines to compare how the initial guess has been improved.
# Plotdisp_show(U0-ResultDisp{ImgSeqNum-1}.U,DICmesh.coordinatesFEMWorld,DICmesh.elementsFEM(:,1:4));
# Plotdisp_show(ResultDisp{ImgSeqNum-2}.U-ResultDisp{ImgSeqNum-1}.U,DICmesh.coordinatesFEMWorld,DICmesh.elementsFEM(:,1:4));
    # ====== Compute image gradients ======
    Df = funImgGradient(fNormalized,gNormalized)
    # ====== Compute f(X)-g(x+u) ======
# PlotImgDiff(x0,y0,u,v,fNormalized,gNormalized);
    ResultFEMeshEachFrame[ImgSeqNum - 1] = struct('coordinatesFEM',DICmesh.coordinatesFEM,'elementsFEM',DICmesh.elementsFEM)
    print('------------ Section 3 Done ------------ \n \n' % ())
    ## Section 4
    print('------------ Section 4 Start ------------ \n' % ())
    # ####################################################################
# Finite element based global DIC iterations
# ####################################################################
    DICpara.tol = 0.001
    DICpara.maxIter = 100
    DICpara.alpha = 10
    alphaList = DICpara.alpha
    # ====== Tune regularization coefficient ======
# If you don't know the best alpha (coefficient), please run the following
# codes to tune the best value of the coefficient of the regularizer |grad u|^2):
    # ##### Uncomment the following line to tune the best value of alpha ######
# alphaList = [1e-3,1e-2,1e-1,1e0,1e1,1e2,1e3]*mean(DICpara.winstepsize);
    Err1List = np.zeros((len(alphaList),1))
    Err2List = Err1List
    UList = cell(len(alphaList),1)
    FList = UList
    # ------------------------------------------------
    for alphaInd in np.arange(1,len(alphaList)+1).reshape(-1):
        tic
        alpha = alphaList(alphaInd)
        # Solve displacement U with each alpha
        U,normOfW,timeICGN = funGlobalICGNQ4(DICmesh,Df,fNormalized,gNormalized,U0,alpha,DICpara.tol,DICpara.maxIter)
        Plotdisp_show(U,DICmesh.coordinatesFEM,DICmesh.elementsFEM)
        # Compute F deformation gradient with solved U
        DICpara.GaussPtOrder = 2
        F = funGlobal_NodalStrainAvg(DICmesh.coordinatesFEM,DICmesh.elementsFEM,U,DICpara.GaussPtOrder)
        # Plotstrain_show(F,DICmesh.coordinatesFEM ,DICmesh.elementsFEM);
        Err1List[alphaInd] = norm(U - U0,2)
        Err2List[alphaInd] = norm(F,2)
    # ====== Tune the coefficient of |grad u| regularizer ======
    ErrSumList = Err1List + 1 * mean(DICpara.winstepsize) * Err2List
    __,indexOfalpha = np.amin(ErrSumList)
    try:
        fitobj = fit(np.transpose(log10(alphaList(np.arange(indexOfalpha - 1,indexOfalpha + 1+1,1)))),ErrSumList(np.arange(indexOfalpha - 1,indexOfalpha + 1+1,1)),'poly2')
        p = coeffvalues(fitobj)
        alpha_best = 10 ** (- p(2) / 2 / p(1))
    finally:
        pass
    DICpara.alpha = alpha_best
    # ====== Re-run global DIC iterations with tuned alpha_best ======
    if np.abs(alpha_best - alpha) > np.abs(eps):
        U,normOfW,timeICGN = funGlobalICGNNonUni(DICmesh,Df,fNormalized,gNormalized,U0,alpha_best,DICpara.tol,DICpara.maxIter)
        DICpara.GaussPtOrder = 2
        F = funGlobal_NodalStrainAvg(DICmesh.coordinatesFEM,DICmesh.elementsFEM,U,DICpara.GaussPtOrder)
    # ------- Smooth strain field --------
    DICpara.DispFilterSize = 0.0001
    DICpara.DispFilterStd = 0
    DICpara.StrainFilterSize = 0.0001
    DICpara.StrainFilterStd = 0
    # F = funSmoothStrain(F,DICmesh,DICpara);
    # ------- Save data -------
    ResultDisp[ImgSeqNum - 1].U = full(U)
    ResultDefGrad[ImgSeqNum - 1].F = full(F)
    ResultAlpha[ImgSeqNum - 1].alpha = alpha_best
    ResultNormOfW[ImgSeqNum - 1].normOfW = full(normOfW)
    ResultTimeICGN[ImgSeqNum - 1].timeICGN = full(timeICGN)
    print('------------ Section 4 Done ------------ \n \n' % ())

# ------ Plot ------
UWorld = U
UWorld[np.arange[2,end()+2,2]] = - U(np.arange(2,end()+2,2))
FWorld = F
FWorld[np.arange[2,end()+4,4]] = - F(np.arange(2,end()+4,4))
FWorld[np.arange[3,end()+4,4]] = - F(np.arange(3,end()+4,4))
close_('all')
Plotdisp_show(UWorld,DICmesh.coordinatesFEMWorld,DICmesh.elementsFEM)
Plotstrain_show(FWorld,DICmesh.coordinatesFEMWorld,DICmesh.elementsFEM)
## ------ Save results ------
# Find img name and save all the results
__,imgname,imgext = os.path.split(file_name[1,end()])[0],os.path.splitext(os.path.split(file_name[1,end()])[1])[0],os.path.splitext(os.path.split(file_name[1,end()])[1])[1]
results_name = np.array(['results_FE_globalDIC_',imgname,'_st',num2str(DICpara.winstepsize),'_alpha',num2str(DICpara.alpha),'.mat'])
save(results_name,'file_name','DICpara','DICmesh','ResultDisp','ResultDefGrad','ResultFEMesh','normOfW','timeICGN')
## Section 5
print('------------ Section 5 Start ------------ \n' % ())
# ####################################################################
# This section is to compute strain fields and plot disp and strain results
# ####################################################################
# ------ Convert units from pixels to the physical world ------
DICpara.um2px = funParaInput('ConvertUnit')
# ------ Smooth displacements ------
DICpara.DoYouWantToSmoothOnceMore = funParaInput('SmoothDispOrNot')
# ------ Choose strain computation method ------
DICpara.MethodToComputeStrain = 0

# ------ Choose strain type (infinitesimal, Eulerian, Green-Lagrangian) ------
DICpara.StrainType = funParaInput('StrainType')
# ------ Choose image to plot (first only, second and next images) ------
if len(ImgNormalized) == 2:
    DICpara.Image2PlotResults = funParaInput('Image2PlotResults')
else:
    DICpara.Image2PlotResults = 1

# ------ Save fig format ------
DICpara.MethodToSaveFig = funParaInput('SaveFigFormat')
# ------ Choose overlay image transparency ------
DICpara.OrigDICImgTransparency = 1
if DICpara.MethodToSaveFig == 1:
    DICpara.OrigDICImgTransparency = funParaInput('OrigDICImgTransparency')

# ------ Start main part ------
for ImgSeqNum in np.arange(2,len(ImgNormalized)+1).reshape(-1):
    print(np.array(['Current image frame #: ',num2str(ImgSeqNum),'/',num2str(len(ImgNormalized))]))
    ###########################################
    fNormalizedNewIndex = ImgSeqNum - np.mod(ImgSeqNum - 2,DICpara.ImgSeqIncUnit) - 1
    if DICpara.ImgSeqIncUnit > 1:
        FEMeshIndLast = int(np.floor(fNormalizedNewIndex / DICpara.ImgSeqIncUnit))
    else:
        if DICpara.ImgSeqIncUnit == 1:
            FEMeshIndLast = int(np.floor(fNormalizedNewIndex / DICpara.ImgSeqIncUnit)) - 1
    FEMeshInd = FEMeshIndLast + 1
    if FEMeshInd == 1:
        USubpb2 = ResultDisp[ImgSeqNum - 1].U
        coordinatesFEM = ResultFEMesh[0].coordinatesFEM
        elementsFEM = ResultFEMesh[0].elementsFEM
        if (ImgSeqNum - 1 == 1) or (DICpara.ImgSeqIncROIUpdateOrNot == 1):
            UFEMesh = 0 * USubpb2
    else:
        USubpb2 = ResultDisp[ImgSeqNum - 1].U
        if np.mod(ImgSeqNum - 2,DICpara.ImgSeqIncUnit) == 0:
            coordinatesFEM = ResultFEMesh[FEMeshInd].coordinatesFEM
            elementsFEM = ResultFEMesh[FEMeshInd].elementsFEM
            coordinatesFEMLast = ResultFEMesh[FEMeshIndLast].coordinatesFEM
            UFEMeshLast = ResultDisp[ImgSeqNum - 2].U + UFEMesh
            xq = coordinatesFEM(:,1)
            yq = coordinatesFEM(:,2)
            UFEMesh = 0 * USubpb2
            UFEMesh[np.arange[1,end()+2,2]] = griddata(coordinatesFEMLast(:,1),coordinatesFEMLast(:,2),UFEMeshLast(np.arange(1,end()+2,2)),xq,yq,'v4')
            UFEMesh[np.arange[2,end()+2,2]] = griddata(coordinatesFEMLast(:,1),coordinatesFEMLast(:,2),UFEMeshLast(np.arange(2,end()+2,2)),xq,yq,'v4')
        USubpb2 = USubpb2 + UFEMesh
    FSubpb2 = ResultDefGrad[ImgSeqNum - 1].F
    coordinatesFEM = ResultFEMeshEachFrame[ImgSeqNum - 1].coordinatesFEM
    elementsFEM = ResultFEMeshEachFrame[ImgSeqNum - 1].elementsFEM
    xList = np.arange(np.amin(coordinatesFEM(:,1)),np.amax(coordinatesFEM(:,1))+DICpara.winstepsize,DICpara.winstepsize)
    M = len(xList)
    yList = np.arange(np.amin(coordinatesFEM(:,2)),np.amax(coordinatesFEM(:,2))+DICpara.winstepsize,DICpara.winstepsize)
    N = len(yList)
    coordinatesFEMWorld = DICpara.um2px * np.array([coordinatesFEM(:,1),ImgNormalized[0].shape[2-1] + 1 - coordinatesFEM(:,2)])
    ###########################################
    # ------ Plotting and Compute Strain-------
    if USubpb2.shape[1-1] == 1:
        ULocal = USubpb2_New.USubpb2
        FLocal = FSubpb2.FSubpb2
    else:
        ULocal = USubpb2
        FLocal = FSubpb2
    UWorld = DICpara.um2px * ULocal
    UWorld[np.arange[2,end()+2,2]] = - UWorld(np.arange(2,end()+2,2))
    # ------ Smooth displacements ------
#prompt = 'Do you want to smooth displacement? (0-yes; 1-no)';
#DoYouWantToSmoothOnceMore = input(prompt);
    SmoothTimes = 0
    try:
        while DICpara.DoYouWantToSmoothOnceMore == 0 and SmoothTimes < 3:

            ULocal = funSmoothDisp(ULocal,DICmesh,DICpara)
            ##DICpara.DoYouWantToSmoothOnceMore = input(prompt);
            SmoothTimes = SmoothTimes + 1

    finally:
        pass
    # ----- Compute strain field ------
    ComputeStrain
    # ##### Add filter and plot strain field #####
# ##### Plotstrain_Fij; #####
    # ------ Plot disp and strain ------
    close_('all')
    if DICpara.OrigDICImgTransparency == 1:
        Plotdisp_show(UWorld,coordinatesFEMWorld,DICmesh.elementsFEM(:,np.arange(1,4+1)),DICpara,'NoEdgeColor')
        strain_exx,strain_exy,strain_eyy,strain_principal_max,strain_principal_min,strain_maxshear,strain_vonMises = Plotstrain0Quadtree(FStrainWorld,coordinatesFEMWorld,DICmesh.elementsFEM(:,np.arange(1,4+1)),DICpara)
    else:
        if DICpara.Image2PlotResults == 0:
            PlotdispQuadtree(UWorld,coordinatesFEMWorld,DICmesh.elementsFEM(:,np.arange(1,4+1)),file_name[1,1],DICpara)
            strain_exx,strain_exy,strain_eyy,strain_principal_max,strain_principal_min,strain_maxshear,strain_vonMises = PlotstrainQuadtree(UWorld,FStrainWorld,coordinatesFEMWorld,DICmesh.elementsFEM(:,np.arange(1,4+1)),file_name[1,1],DICpara)
        else:
            PlotdispQuadtree(UWorld,coordinatesFEMWorld,DICmesh.elementsFEM(:,np.arange(1,4+1)),file_name[1,ImgSeqNum],DICpara)
            strain_exx,strain_exy,strain_eyy,strain_principal_max,strain_principal_min,strain_maxshear,strain_vonMises = PlotstrainQuadtree(UWorld,FStrainWorld,coordinatesFEMWorld,DICmesh.elementsFEM(:,np.arange(1,4+1)),file_name[1,ImgSeqNum],DICpara)
    # ----- Save strain results ------
    ResultStrainWorld[ImgSeqNum - 1] = struct('strainxCoord',coordinatesFEMWorld(:,1),'strainyCoord',coordinatesFEMWorld(:,2),'dispu',UWorld(np.arange(1,end()+2,2)),'dispv',UWorld(np.arange(2,end()+2,2)),'dudx',FStrainWorld(np.arange(1,end()+4,4)),'dvdx',FStrainWorld(np.arange(2,end()+4,4)),'dudy',FStrainWorld(np.arange(3,end()+4,4)),'dvdy',FStrainWorld(np.arange(4,end()+4,4)),'strain_exx',strain_exx,'strain_exy',strain_exy,'strain_eyy',strain_eyy,'strain_principal_max',strain_principal_max,'strain_principal_min',strain_principal_min,'strain_maxshear',strain_maxshear,'strain_vonMises',strain_vonMises)
    # ------ Save figures for tracked displacement and strain fields ------
    SaveFigFilesDispAndStrain

# ------ END of for-loop {ImgSeqNum = 2:length(ImgNormalized)} ------
print('------------ Section 5 Done ------------ \n \n' % ())
## Save data again including stress solve method
results_name = np.array(['results_FE_globalDIC_',imgname,'_st',num2str(DICpara.winstepsize(1)),'_alpha',num2str(DICpara.alpha),'.mat'])
save(results_name,'file_name','DICpara','DICmesh','ResultDisp','ResultDefGrad','ResultStrainWorld','ResultFEMesh','ResultFEMeshEachFrame','ResultAlpha','ResultNormOfW','ResultTimeICGN')
## Section 6: Compute stress
print('------------ Section 6 Start ------------ \n' % ())
# ####################################################################
# This section is to compute stress fields and plot stress fields
# ####################################################################
# ------ Choose material model ------
DICpara.MaterialModel = funParaInput('MaterialModel')
# ------ Define parameters in material models ------
if (DICpara.MaterialModel == 1) or (DICpara.MaterialModel == 2):
    print('Define Linear elasticity parameters \n' % ())
    print('Young's modulus (unit: Pa): \n' % ())
    prompt = 'Input here (e.g., 69e9): '
    DICpara.MaterialModelPara.YoungsModulus = input_(prompt)
    print('Poisson's ratio: \n' % ())
    prompt = 'Input here (e.g., 0.3): '
    DICpara.MaterialModelPara.PoissonsRatio = input_(prompt)
    print('------------------------------------- \n' % ())

# ------ Start main part ------
for ImgSeqNum in np.arange(2,len(ImgNormalized)+1).reshape(-1):
    print(np.array(['Current image frame #: ',num2str(ImgSeqNum),'/',num2str(len(ImgNormalized))]))
    close_('all')
    coordinatesFEM = ResultFEMeshEachFrame[ImgSeqNum - 1].coordinatesFEM
    elementsFEM = ResultFEMeshEachFrame[ImgSeqNum - 1].elementsFEM
    coordinatesFEMWorldDef = DICpara.um2px * np.array([coordinatesFEM(:,1),ImgNormalized[0].shape[2-1] + 1 - coordinatesFEM(:,2)]) + DICpara.Image2PlotResults * np.array([ResultStrainWorld[ImgSeqNum - 1].dispu,ResultStrainWorld[ImgSeqNum - 1].dispv])
    # ------ Plot stress ------
    if DICpara.OrigDICImgTransparency == 1:
        stress_sxx,stress_sxy,stress_syy,stress_principal_max_xyplane,stress_principal_min_xyplane,stress_maxshear_xyplane,stress_maxshear_xyz3d,stress_vonMises = Plotstress0Quadtree(DICpara,ResultStrainWorld[ImgSeqNum - 1],coordinatesFEMWorldDef,elementsFEM(:,np.arange(1,4+1)))
    else:
        if DICpara.Image2PlotResults == 0:
            stress_sxx,stress_sxy,stress_syy,stress_principal_max_xyplane,stress_principal_min_xyplane,stress_maxshear_xyplane,stress_maxshear_xyz3d,stress_vonMises = PlotstressQuadtree(DICpara,ResultStrainWorld[ImgSeqNum - 1],coordinatesFEMWorldDef,elementsFEM(:,np.arange(1,4+1)),file_name[1,1])
        else:
            stress_sxx,stress_sxy,stress_syy,stress_principal_max_xyplane,stress_principal_min_xyplane,stress_maxshear_xyplane,stress_maxshear_xyz3d,stress_vonMises = PlotstressQuadtree(DICpara,ResultStrainWorld[ImgSeqNum - 1],coordinatesFEMWorldDef,elementsFEM(:,np.arange(1,4+1)),file_name[1,ImgSeqNum])
    # ------ Save figures for computed stress fields ------
    SaveFigFilesStress
    # ----- Save strain results ------
    ResultStressWorld[ImgSeqNum - 1] = struct('stressxCoord',ResultStrainWorld[ImgSeqNum - 1].strainxCoord,'stressyCoord',ResultStrainWorld[ImgSeqNum - 1].strainyCoord,'stress_sxx',stress_sxx,'stress_sxy',stress_sxy,'stress_syy',stress_syy,'stress_principal_max_xyplane',stress_principal_max_xyplane,'stress_principal_min_xyplane',stress_principal_min_xyplane,'stress_maxshear_xyplane',stress_maxshear_xyplane,'stress_maxshear_xyz3d',stress_maxshear_xyz3d,'stress_vonMises',stress_vonMises)

# ------ END of for-loop {ImgSeqNum = 2:length(ImgNormalized)} ------
print('------------ Section 6 Done ------------ \n \n' % ())
## Save data again including stress solve method
results_name = np.array(['results_FE_globalDIC_',imgname,'_st',num2str(DICpara.winstepsize(1)),'_alpha',num2str(DICpara.alpha),'.mat'])
save(results_name,'file_name','DICpara','DICmesh','ResultDisp','ResultDefGrad','ResultStrainWorld','ResultStressWorld','ResultFEMesh','ResultFEMeshEachFrame','ResultAlpha','ResultNormOfW','ResultTimeICGN')