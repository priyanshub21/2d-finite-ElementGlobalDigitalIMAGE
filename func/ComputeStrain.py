# ==================================
# To compute strain on a uniform mesh
# ----------------------------------
# switch MethodToComputeStrain
#   case 0: direct solved results;
#   case 1: central finite difference;
#   case 2: plane fitting method;
#   case 3: finite element Gauss points;
# ==================================

import numpy as np
if 1 == DICpara.MethodToComputeStrain:
    # Compute strain method I: Use Finite difference operator
    D = funDerivativeOp(M,N,DICpara.winstepsize)
    FStrain = D * reshape(ULocal,len(ULocal),1)
    Rad = 1
    # Find coordinatesFEM that belong to (x(Rad+1:M-Rad,Rad+1:N-Rad),y(Rad+1:M-Rad,Rad+1:N-Rad))
    temp = np.arange(1,coordinatesFEM.shape[1-1]+1,1)
    temp = np.transpose(temp)
    temp = reshape(temp,M,N)
    temp2 = temp(np.arange(Rad + 1,M - Rad+1),np.arange(Rad + 1,N - Rad+1))
    temp2 = reshape(temp2,(M - 2 * Rad) * (N - 2 * Rad),1)
    temp3 = np.zeros((4 * (M - 2 * Rad) * (N - 2 * Rad),1))
    for i in np.arange(1,(M - 2 * Rad) * (N - 2 * Rad)+1).reshape(-1):
        temp3[np.arange[4 * i - 3,4 * i+1]] = np.array([[4 * temp2(i) - 3],[4 * temp2(i) - 2],[4 * temp2(i) - 1],[4 * temp2(i)]])
    FStraintemp = FStrain(temp3)
    #######################################################################
else:
    if 2 == DICpara.MethodToComputeStrain:
        D = funDerivativeOp(M,N,DICpara.winstepsize)
        FStrain = D * reshape(ULocal,len(ULocal),1)
        # Compute strain method II: Use Plane Fitting method
        prompt = 'What is your half window size: '
        Rad = input_(prompt)
        Uytemp,Uxtemp,UNewtemp = PlaneFit2(reshape(ULocal(np.arange(1,end()+2,2)),M,N),DICpara.winstepsize,DICpara.winstepsize,Rad)
        Vytemp,Vxtemp,VNewtemp = PlaneFit2(reshape(ULocal(np.arange(2,end()+2,2)),M,N),DICpara.winstepsize,DICpara.winstepsize,Rad)
        FStraintemp = np.zeros((4 * (M - 2 * Rad) * (N - 2 * Rad),1))
        FStraintemp[np.arange[1,end()+4,4]] = reshape(Uxtemp(np.arange((Rad + 1),M - Rad+1),np.arange((Rad + 1),N - Rad+1)),(M - 2 * Rad) * (N - 2 * Rad),1)
        FStraintemp[np.arange[2,end()+4,4]] = reshape(Vxtemp(np.arange((Rad + 1),M - Rad+1),np.arange((Rad + 1),N - Rad+1)),(M - 2 * Rad) * (N - 2 * Rad),1)
        FStraintemp[np.arange[3,end()+4,4]] = reshape(Uytemp(np.arange((Rad + 1),M - Rad+1),np.arange((Rad + 1),N - Rad+1)),(M - 2 * Rad) * (N - 2 * Rad),1)
        FStraintemp[np.arange[4,end()+4,4]] = reshape(Vytemp(np.arange((Rad + 1),M - Rad+1),np.arange((Rad + 1),N - Rad+1)),(M - 2 * Rad) * (N - 2 * Rad),1)
        # Find coordinatesFEM that belong to (x(Rad+1:M-Rad,Rad+1:N-Rad),y(Rad+1:M-Rad,Rad+1:N-Rad))
        temp = np.arange(1,coordinatesFEM.shape[1-1]+1,1)
        temp = np.transpose(temp)
        temp = reshape(temp,M,N)
        temp2 = temp(np.arange(Rad + 1,M - Rad+1),np.arange(Rad + 1,N - Rad+1))
        temp2 = reshape(temp2,(M - 2 * Rad) * (N - 2 * Rad),1)
        temp3 = np.zeros((4 * (M - 2 * Rad) * (N - 2 * Rad),1))
        for i in np.arange(1,(M - 2 * Rad) * (N - 2 * Rad)+1).reshape(-1):
            temp3[np.arange[4 * i - 3,4 * i+1]] = np.array([[4 * temp2(i) - 3],[4 * temp2(i) - 2],[4 * temp2(i) - 1],[4 * temp2(i)]])
        FStrain[temp3] = FStraintemp
        #######################################################################
    else:
        if 0 == DICpara.MethodToComputeStrain:
            GaussPtOrder = 2
            FStrain = funGlobal_NodalStrainAvg(coordinatesFEM,elementsFEM,USubpb2,GaussPtOrder)
            Rad = 0
            # Find coordinatesFEM that belong to (x(Rad+1:M-Rad,Rad+1:N-Rad),y(Rad+1:M-Rad,Rad+1:N-Rad))
# temp = 1:1:size(coordinatesFEM,1); temp = temp';
# temp = reshape(temp,M,N); temp2 = temp(Rad+1:M-Rad, Rad+1:N-Rad);
# temp2 = reshape(temp2, (M-2*Rad)*(N-2*Rad),1);
            # temp3 = zeros(4*(M-2*Rad)*(N-2*Rad),1);
# for i = 1:(M-2*Rad)*(N-2*Rad)
#     temp3(4*i-3:4*i) = [4*temp2(i)-3; 4*temp2(i)-2; 4*temp2(i)-1; 4*temp2(i)];
# end
# FStraintemp = FStrain(temp3);
            #######################################################################
        else:
            print('Wrong Input to compute strain field!')

## Update infinitesimal strain to other finite strains
FStrainFinite = FStrain
for tempi in np.arange(1,len(FStrain)+4,4).reshape(-1):
    # Obtain each component of def grad tensor
    dudx = FStrain(tempi)
    dvdx = FStrain(tempi + 1)
    dudy = FStrain(tempi + 2)
    dvdy = FStrain(tempi + 3)
    if 0 == DICpara.StrainType:
        # Do nothing
        pass
    else:
        if 1 == DICpara.StrainType:
            FStrainFinite[tempi] = 1 / (1 - dudx) - 1
            FStrainFinite[tempi + 3] = 1 / (1 - dvdy) - 1
            FStrainFinite[tempi + 2] = dudy / (1 - dvdy)
            FStrainFinite[tempi + 1] = dvdx / (1 - dudx)
        else:
            if 2 == DICpara.StrainType:
                FStrainFinite[tempi] = 0.5 * (dudx * 2 - dudx ** 2 - dvdx ** 2)
                FStrainFinite[tempi + 3] = 0.5 * (dvdy * 2 - dudy ** 2 - dvdy ** 2)
                FStrainFinite[tempi + 2] = 0.5 * (dudy + dvdx - dudx * dudy - dvdx * dvdy)
                FStrainFinite[tempi + 1] = 0.5 * (dvdx + dudy - dudy * dudx - dvdy * dvdx)
            else:
                if 3 == DICpara.StrainType:
                    print('Press "Ctrl+C" to modify by yourself.')
                    pause
                else:
                    print('Wrong strain type!')

try:
    FStraintemp = FStrainFinite(temp3)
finally:
    pass

FStrainWorld = FStraintemp
FStrainWorld[np.arange[2,end()+4,4]] = - FStrainWorld(np.arange(2,end()+4,4))
FStrainWorld[np.arange[3,end()+4,4]] = - FStrainWorld(np.arange(3,end()+4,4))