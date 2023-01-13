


import numpy as np
import matplotlib.pyplot as plt
U0World = U0
U0World[np.arange[2,end()+2,2]] = - U0(np.arange(2,end()+2,2))
# DICmesh.y0World = (size(fNormalized,2)+1-y0); DICmesh.coordinatesFEMWorld = [DICmesh.coordinatesFEM(:,1),size(fNormalized,2)+1-DICmesh.coordinatesFEM(:,2)];

close_('all')
Plotuv(U0World,DICmesh.x0,DICmesh.y0World)

# figure(1); view([-140,60]);# colormap(coolwarm(32)); # caxis([-1.0,1.0]); # caxis([-0.45,0.45]);
# figure(2); view([-140,60]);# colormap(coolwarm(32)); # caxis([-0.085, 0.005]); # caxis([-0.085,0.015 ]);
Plotdisp_show(U0World,DICmesh.coordinatesFEMWorld,DICmesh.elementsFEM)

plt.figure(3)
plt.axis('equal')

plt.figure(4)
plt.axis('equal')
