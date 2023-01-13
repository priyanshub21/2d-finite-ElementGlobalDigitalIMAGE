import numpy as np
    
def PlotStrainErr(coordinatesFEM = None,FSubpb2 = None,x = None,y = None,M = None,N = None,Exact = None): 
    ErrStrain = np.zeros((coordinatesFEM.shape[1-1],1))
    avgErrStrain = 0
    for i in np.arange(1,coordinatesFEM.shape[1-1]+1).reshape(-1):
        ErrStrain[i] = np.sqrt((FSubpb2(4 * i - 3) + 1e-06 * Exact(coordinatesFEM(i,1),7)) ** 2 + (FSubpb2(4 * i - 2) - 0) ** 2 + (FSubpb2(4 * i - 1) - 0) ** 2 + (FSubpb2(4 * i) - 0) ** 2)
        avgErrStrain = avgErrStrain + ErrStrain(i) ** 2
    
    avgErrStrain = np.sqrt(avgErrStrain / coordinatesFEM.shape[1-1])
    # figure; surf(reshape(ErrDisp,M,N)'); axis equal; colorbar; view(2); title('||u-u_0||_{L_2}^{1/2}')
# figure; mesh(x,y,reshape(ErrStrain,M,N)); axis tight; set(gca,'fontSize',20); view(-20, 50); title('Strain Absolute Error')
    
    Strainxx = np.zeros((coordinatesFEM.shape[1-1],1))
    avgStrainxx = 0
    stdStrainxx = 0
    xAxis = 0
    ExactStrainxx = 0
    for i in np.arange(1,coordinatesFEM.shape[1-1]+1).reshape(-1):
        Strainxx[i] = FSubpb2(4 * i - 3)
    
    Strainxx = reshape(Strainxx,M,N)
    for i in np.arange(1,M+1).reshape(-1):
        avgStrainxx[i] = sum(Strainxx(i,:)) / N
        stdStrainxx[i] = std(Strainxx(i,:))
        ExactStrainxx[i] = - 1e-06 * Exact(x(i,1),7)
    
    avgStrainxx = np.transpose(avgStrainxx)
    stdStrainxx = np.transpose(stdStrainxx)
    ExactStrainxx = np.transpose(ExactStrainxx)
    xAxis = x(:,1)
    #figure; errorbar(xAxis, avgStrainxx  , stdStrainxx );
    return xAxis,avgStrainxx,stdStrainxx,ExactStrainxx