#close all;

import numpy as np
import numpy.matlib
clear('u_','dm','m2vx')
L = 9
u_ = cell(1)
temp1 = reshape(ULocal(np.arange(1,end()+2,2)),M,N)
temp1 = np.matlib.repmat(temp1,1,1,L)
temp10 = np.zeros((N,M,L))
for tempkkk in np.arange(1,L+1).reshape(-1):
    temp10[:,:,tempkkk] = np.transpose(temp1(:,:,tempkkk))

u_[0][0] = temp10
temp1 = reshape(ULocal(np.arange(2,end()+2,2)),M,N)
temp1 = np.matlib.repmat(temp1,1,1,L)
temp10 = np.zeros((N,M,L))
for tempkkk in np.arange(1,L+1).reshape(-1):
    temp10[:,:,tempkkk] = np.transpose(temp1(:,:,tempkkk))

u_[2][0] = temp10
u_[3][0] = 0 * temp10
temp1 = reshape(np.sqrt((ULocal(np.arange(1,end()+2,2)) ** 2 + ULocal(np.arange(2,end()+2,2)) ** 2)),M,N)
temp1 = np.matlib.repmat(temp1,1,1,L)
temp10 = np.zeros((N,M,L))
for tempkkk in np.arange(1,L+1).reshape(-1):
    temp10[:,:,tempkkk] = np.transpose(temp1(:,:,tempkkk))

u_[4][0] = temp10
dm = np.array([winstepsize])
m2vx = np.array([1,1,1])
Fij,J = calculateFij(u_,dm,m2vx,'optimal9')
## Calculate Lagrangian strain (Eij), eulerian strain (eij), Strain energy density (U), and stress (Sij)
Eij,eij = calculateEij(Fij)
#[Sij, U] = calculateSij(Fij, J, eij{1,1}, Eij{1,1}, [E v], mechModel);
#save('mechanicsvariables_dyedcell_xy5.mat', 'Fij','J','Eij','eij','U','Sij');

##
FStrain = np.zeros((M * N * 4,1))
FStrain[np.arange[1,end()+4,4]] = reshape(Eij[0][5](:,:,1),M * N,1)
FStrain[np.arange[4,end()+4,4]] = reshape(Eij[0][0](:,:,1),M * N,1)
# FSubpb3(9:9:end) = reshape(Fij{1}{9},M*N*L,1)-1;
FStrain[np.arange[2,end()+4,4]] = reshape(Eij[0][4](:,:,1),M * N,1)
FStrain[np.arange[3,end()+4,4]] = reshape(Eij[0][2](:,:,1),M * N,1)
# FSubpb3(3:9:end) = reshape(Fij{1}{6},M*N*L,1);
# FSubpb3(7:9:end) = reshape(Fij{1}{8},M*N*L,1);
# FSubpb3(6:9:end) = reshape(Fij{1}{3},M*N*L,1);
# FSubpb3(8:9:end) = reshape(Fij{1}{7},M*N*L,1);
#for index=[1:9]
#FSubpb3(index:9:end) = reshape(Fij{1}{index},M*N*L,1);
#end
#Plotstrain_show3(FSubpb3,ResultcoordinatesFEM{1}.coordinatesFEM,ResultelementsFEM{1}.elementsFEM);

##
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