import numpy as np
    
def PlotDispErr(coordinatesFEM = None,USubpb2 = None,x = None,y = None,M = None,N = None,Exact = None): 
    ErrDisp = np.zeros((coordinatesFEM.shape[1-1],1))
    avgDispErr = 0
    for i in np.arange(1,coordinatesFEM.shape[1-1]+1).reshape(-1):
        ErrDisp[i] = np.sqrt((USubpb2(2 * i - 1) + Exact(coordinatesFEM(i,1),6)) ** 2 + (USubpb2(2 * i) - 0) ** 2)
        avgDispErr = avgDispErr + ErrDisp(i) ** 2
    
    avgDispErr = np.sqrt(avgDispErr / coordinatesFEM.shape[1-1])
    # figure; surf(reshape(ErrDisp,M,N)'); axis equal; colorbar; view(2); title('||u-u_0||_{L_2}^{1/2}')
# figure; mesh(x,y,reshape(ErrDisp,M,N)); axis tight; set(gca,'fontSize',18); view(-20, 50); title('Displacement Absolute Error')
    
    Dispxx = np.zeros((coordinatesFEM.shape[1-1],1))
    avgDispxx = 0
    stdDispxx = 0
    xAxis = 0
    ExactDispxx = 0
    for i in np.arange(1,coordinatesFEM.shape[1-1]+1).reshape(-1):
        Dispxx[i] = USubpb2(2 * i - 1)
    
    Dispxx = reshape(Dispxx,M,N)
    for i in np.arange(1,M+1).reshape(-1):
        avgDispxx[i] = sum(Dispxx(i,:)) / N
        stdDispxx[i] = std(Dispxx(i,:))
        ExactDispxx[i] = - Exact(x(i,1),6)
    
    avgDispxx = np.transpose(avgDispxx)
    stdDispxx = np.transpose(stdDispxx)
    ExactDispxx = np.transpose(ExactDispxx)
    xAxis = x(:,1)
    return xAxis,avgDispxx,stdDispxx,ExactDispxx