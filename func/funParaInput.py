import numpy as np
    
def funParaInput(paraName = None): 
    # FUNCTION paraInput = funParaInput(paraName)
# To define DIC parameters from user's input
# ----------------------------------------------
    
    #   INPUT:  paraName        Parameter name
    
    #   OUTPUT: paraInput       Input DIC parameters
    
    # ----------------------------------------------
# Author: Jin Yang.
# Contact and support: jyang526@wisc.edu -or- aldicdvc@gmail.com
# Last time updated: 02/2020.
# ==============================================
    
    ##
    if 'InitFFTSearchMethod' == paraName:
        print('\n' % ())
        print('--- Method to compute an initial guess of displacements --- \n' % ())
        print('    0: Multigrid search based on an image pyramid  \n' % ())
        print('    1: Whole field search for all the subsets  \n' % ())
        print('    2: Search near manually clicked seeds and then interpolate for the full-field  \n' % ())
        prompt = 'Input here: '
        InitFFTSearchMethod = input_(prompt)
        paraInput = InitFFTSearchMethod
        print('\n' % ())
    else:
        if 'NewFFTSearch' == paraName:
            print('\n' % ())
            print('Since we are dealing with an image sequence, for each new frame,   \n' % ())
            print('do we use last frame result as the initial guess or   \n' % ())
            print('Redo FFT initial guess for every new frame? \n' % ())
            print('    0: Use last frame (by default); \n' % ())
            print('    1: Redo initial guess.  \n' % ())
            prompt = 'Input here: '
            StillFFTSearch = input_(prompt)
            paraInput = StillFFTSearch
            print('\n' % ())
        else:
            if 'Subpb2FDOrFEM' == paraName:
                print('\n' % ())
                print('--- Method to solve ALDIC global step Subproblem 2 ---    \n' % ())
                print('    1: Finite difference (Recommend)   \n' % ())
                print('    2: Finite element method  \n' % ())
                prompt = 'Input here: '
                Subpb2FDOrFEM = input_(prompt)
                paraInput = Subpb2FDOrFEM
            else:
                if 'ClusterNo' == paraName:
                    print('\n' % ())
                    print('--- Set up Parallel pool ---')
                    print('How many parallel pools to open? (Put in 1 if no parallel computing) \n' % ())
                    prompt = 'Input here: '
                    ClusterNo = input_(prompt)
                    # if ClusterNo > 1 # Backup codes
#     delete(gcp); myCluster = parcluster('local'); delete(myCluster.Jobs);
#     parpool(ClusterNo,'SpmdEnabled',false);
# end
                    paraInput = ClusterNo
                    #######################################################################
##### Section 8 #####
                else:
                    if 'ConvertUnit' == paraName:
                        print('Convert units from pixels to the physical world. \n' % ())
                        display('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
                        print('Results in Section 5 (ResultStrain), Section 6 \n(ResultStress), and all the plots will be converted \nto the physical world unit instead of the pixel unit.  \n' % ())
                        display('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
                        print('If you want to keep the pixel unit, enter '1'. \n' % ())
                        prompt = 'Input here (e.g., mm/px, um/px): '
                        um2px = input_(prompt)
                        paraInput = um2px
                        print('------------------------------------- \n' % ())
                    else:
                        if 'SmoothDispOrNot' == paraName:
                            print('Do you want to smooth displacement? (0-yes; 1-no) \n' % ())
                            prompt = 'Input here: '
                            DoYouWantToSmoothOnceMore = input_(prompt)
                            paraInput = DoYouWantToSmoothOnceMore
                            print('------------------------------------- \n' % ())
                        else:
                            if 'RegularizationSmoothness' == paraName:
                                print('Regularization smoothness \n' % ())
                                prompt = 'Input here (e.g., 0 or 1e-5~1e-3): '
                                DoYouWantToSmoothOnceMore = input_(prompt)
                                paraInput = DoYouWantToSmoothOnceMore
                                print('------------------------------------- \n' % ())
                            else:
                                if 'StrainMethodOp' == paraName:
                                    print('What method to use to compute strain? \n' % ())
                                    print('    0: Finite element (Recommend); \n' % ())
                                    print('    1: Finite difference ; \n' % ())
                                    print('    2: Plane fitting; \n' % ())
                                    prompt = 'Input here: '
                                    MethodToComputeStrain = input_(prompt)
                                    while (MethodToComputeStrain != 0) and (MethodToComputeStrain != 1) and (MethodToComputeStrain != 2) and (MethodToComputeStrain != 3):

                                        print('****** Wrong input! ******')
                                        print('What method to use to compute strain? \n' % ())
                                        print('    0: Finite element (Recommend); \n' % ())
                                        print('    1: Finite difference ; \n' % ())
                                        print('    2: Plane fitting; \n' % ())
                                        prompt = 'Input here: '
                                        MethodToComputeStrain = input_(prompt)

                                    paraInput = np.array([MethodToComputeStrain])
                                    print('------------------------------------- \n' % ())
                                else:
                                    if 'StrainType' == paraName:
                                        print('Infinitesimal strain or finite strain? \n' % ())
                                        print('    0: Infinitesimal strain; \n' % ())
                                        print('    1: Eulerian strain; \n' % ())
                                        print('    2: Green-Lagrangian strain; \n' % ())
                                        print('    3: Others: code by yourself; \n' % ())
                                        prompt = 'Input here: '
                                        StrainType = input_(prompt)
                                        while (StrainType != 0) and (StrainType != 1) and (StrainType != 2) and (StrainType != 3):

                                            print('****** Wrong input! ******')
                                            print('Infinitesimal strain or finite strain? \n' % ())
                                            print('    0: Infinitesimal strain; \n' % ())
                                            print('    1: Eluerian strain; \n' % ())
                                            print('    2: Green-Lagrangian strain; \n' % ())
                                            print('    3: Others: code by yourself; \n' % ())
                                            prompt = 'Input here: '
                                            StrainType = input_(prompt)

                                        paraInput = StrainType
                                        print('------------------------------------- \n' % ())
                                    else:
                                        if 'Image2PlotResults' == paraName:
                                            print('Over which image(s) you want to plot the results? \n' % ())
                                            print('    0: First image only (if you are using only two images); \n' % ())
                                            print('    1: Second and next images; \n' % ())
                                            prompt = 'Input here: '
                                            Image2PlotResults = input_(prompt)
                                            paraInput = Image2PlotResults
                                            print('------------------------------------- \n' % ())
                                        else:
                                            if 'SaveFigFormat' == paraName:
                                                print('Save figures into the format: \n' % ())
                                                print('    1: jpg(Choose transparency 0~1) \n' % ())
                                                print('    2: pdf(Choose transparency = 1) \n' % ())
                                                print('    3: Others: Edit codes in ./plotFiles/SaveFigFiles.m \n' % ())
                                                prompt = 'Input here: '
                                                MethodToSaveFig = input_(prompt)
                                                paraInput = MethodToSaveFig
                                                print('------------------------------------- \n' % ())
                                            else:
                                                if 'OrigDICImgTransparency' == paraName:
                                                    print('Define transparency for overlaying original images: \n' % ())
                                                    print('Input a real number between 0(Only original images) \n' % ())
                                                    print('and 0.99 (Non-transparent deformation results).\n' % ())
                                                    prompt = 'Input here(e.g., 0.5): '
                                                    OrigDICImgTransparency = input_(prompt)
                                                    paraInput = OrigDICImgTransparency
                                                    print('------------------------------------- \n' % ())
                                                    #######################################################################
##### Section 9 #####
                                                else:
                                                    if 'MaterialModel' == paraName:
                                                        print('Material model to compute Cauchy stress fields: \n' % ())
                                                        print('    1: Linear elasticity -- Plane stress \n' % ())
                                                        print('    2: Linear elasticity -- Plane strain \n' % ())
                                                        print('    3: Others: User needs to edit codes in ./plotFiles/Plotstress.m   \n' % ())
                                                        prompt = 'Input here: '
                                                        MaterialModel = input_(prompt)
                                                        paraInput = MaterialModel
                                                        print('------------------------------------- \n' % ())
    
    return paraInput