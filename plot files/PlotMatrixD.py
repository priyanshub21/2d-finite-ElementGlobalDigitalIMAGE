## Plot D matrix
import numpy as np
import matplotlib.pyplot as plt
close_('all')
tempD = DerivativeOp(10,10,1)
tempy,tempx = np.meshgrid(np.arange(1,tempD.shape[1-1]+1,1),np.arange(1,tempD.shape[2-1]+1,1))
figure
imshow(int8(full(tempD + 1) * 0.5 * 128))
plt.axis('tight')
b = colorbar
caxis(np.array([0,128]))
b.Ticks = np.array([0,64,128])
b.TickLabels = np.array(['-1','0','1'])
set(gca,'fontsize',40)
fig = gcf
fig.PaperUnits = 'inches'
fig.PaperPosition = np.array([0,0,9,9])
print_('Fig_MatrixD_10by10','-dpdf','-r600')