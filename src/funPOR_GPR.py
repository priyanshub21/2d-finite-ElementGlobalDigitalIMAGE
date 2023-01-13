import numpy as np
    
def funPOR_GPR(T_snap = None,t_train = None,t_pre = None,nB = None): 
    # Input:
# T_snap = nT * np matrix; snapshots data
# t_train = nT * 1 vector; time for snapshots data
# t_pre = nT1 * vector; time to predict
# nB: number of basis;
# Output:
# u_pred = nT1 * np matrix; Prediction matrix
# alpha, beta_POD, beta_POD = some nondimensional parameters
    
    np = len(T_snap(1,:))
    nT = len(T_snap(:,1))
    nT1 = len(t_pre(:,1))
    # POD
# mean and fluctuation
    ave = mean(T_snap)
    fp = np.zeros((nT,len(T_snap(1,:))))
    for k in np.arange(1,nT+1).reshape(-1):
        fp[k,:] = T_snap(k,:) - ave
    
    # correlation matrix
    c = np.zeros((nT,nT))
    for k in np.arange(1,nT+1).reshape(-1):
        for j in np.arange(1,nT+1).reshape(-1):
            c[k,j] = innerproduct(fp(k,:),fp(j,:))
    
    # POD decomposition
    w,d = eig(c)
    lambda_ = wrev(diag(d))
    w = fliplr(w)
    beta_POD = (np.log(lambda_(2)) - np.log(lambda_(nB + 2))) / nB
    sum1 = sum(lambda_)
    w = w(:,np.arange(1,nB+1))
    fi = np.transpose(fp) * w
    fi = np.transpose(fi)
    lambda_ = lambda_(np.arange(1,nB+1))
    sum2 = sum(lambda_)
    alpha = sum2 / sum1
    # fi: basis
    for j in np.arange(1,nB+1).reshape(-1):
        fi[j,:] = fi(j,:) / np.sqrt(lambda_(j))
    
    # get a_exact
    a_exact = np.zeros((nT,nB))
    for j in np.arange(1,nT+1).reshape(-1):
        for k in np.arange(1,nB+1).reshape(-1):
            a_exact[j,k] = innerproduct(fp(j,:),fi(k,:))
    
    # GP regression
    ypred = np.zeros((nT1,nB))
    ysd = np.zeros((nT1,nB))
    yint = np.zeros((nT1,nB * 2))
    for k in np.arange(1,nB+1).reshape(-1):
        #     beta=[1;1;0.01;0.001];
        gprMdl = fitrgp(t_train,a_exact(:,k),'Basis','linear','FitMethod','exact','PredictMethod','exact','SigmaLowerBound',2e-05)
        ypred[:,k],ysd[:,k],yint[:,np.arange[2 * k - 1,2 * k+1]] = predict(gprMdl,t_pre)
    
    # compute prediction at t_pre
    u_pred = ypred * fi + ave
    # compute beta_GPR
    sd = ysd * lambda_
    weighted_sum = np.abs(ypred) * lambda_
    beta_GPR = sd / weighted_sum
    
def innerproduct(f1 = None,f2 = None): 
    #   inner product of f1 and f2
    N = len(f1) - 1
    # x=0:1/N:1;
    y = np.multiply(f1,f2)
    innerproduct = y(1) / 2
    for i in np.arange(2,N+1).reshape(-1):
        innerproduct = innerproduct + y(i)
    
    innerproduct = 0.1 / N * (innerproduct + y(N + 1) / 2)
    return innerproduct
    
    return innerproduct
    
    return u_pred,alpha,beta_POD,beta_GPR