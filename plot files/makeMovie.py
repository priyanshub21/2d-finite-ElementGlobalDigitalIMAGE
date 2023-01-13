import numpy as np
import matplotlib.pyplot as plt
files = dir('*_strain_eyy.jpg')
ImgGrayScaleMax = 255
im = cell(len(files),1)
for i in np.arange(1,len(files)+1).reshape(-1):
    im[i] = files(i).name

v = VideoWriter('video_strain_eyy.mp4')
v.FrameRate = 10
open_(v)
figure
for tempk in np.array([np.arange(1,len(im)+1,1)]).reshape(-1):
    clf
    #imshow(imread(im{tempk},1),'DisplayRange',[0,ImgGrayScaleMax]);
    imshow(imread(im[tempk]))
    plt.title(np.array(['Frame #',num2str(tempk + 1)]))
    # text(830,100,['Frame #',num2str(tempk+1)]);
    frame = getframe(gcf)
    writeVideo(v,frame)
    #waitbar(tempk/length(files));

close_(v)
##

alpha = np.arange(10,1000+1,1)
R = 66
R0 = 25
F11 = alpha ** 2.0 / (alpha ** 3 + R ** 3 - R0 ** 3) ** (2 / 3)
figure
plt.plot(alpha,F11,'.')