import numpy as np
    
def show(elements3 = None,elements4 = None,coordinates = None,u = None,varargin = None): 
    if 5 == len(varargin):
        edgeColorOrNot = varargin[0]
    
    try:
        temp = str(edgeColorOrNot) == str('NoEdgeColor') == 1
    finally:
        pass
    
    # trisurf(elements3,coordinates(:,1),coordinates(:,2),u',...
# 'facecolor','interp' )
# hold on
    
    Trix = np.zeros((3,elements3.shape[1-1]))
    Triy = np.zeros((3,elements3.shape[1-1]))
    Tric = np.zeros((3,elements3.shape[1-1]))
    for j in np.arange(1,elements3.shape[1-1]+1).reshape(-1):
        Trix[np.arange[1,3+1],j] = coordinates(elements3(j,np.arange(1,3+1)),1)
        Triy[np.arange[1,3+1],j] = coordinates(elements3(j,np.arange(1,3+1)),2)
        Tric[np.arange[1,3+1],j] = u(elements3(j,np.arange(1,3+1)))
    
    if elements3.shape[1-1] > 20000.0 or str(edgeColorOrNot) == str('NoEdgeColor') == 1:
        h = patch(Trix,Triy,Tric,'facecolor','interp','edgecolor','none')
    else:
        h = patch(Trix,Triy,Tric,'facecolor','interp')
    
    hold('on')
    # trisurf(elements4,coordinates(:,1),coordinates(:,2),u',...
# 'facecolor','interp' )
    
    Sqx = np.zeros((4,elements4.shape[1-1]))
    Sqy = np.zeros((4,elements4.shape[1-1]))
    Sqc = np.zeros((4,elements4.shape[1-1]))
    for j in np.arange(1,elements4.shape[1-1]+1).reshape(-1):
        Sqx[np.arange[1,4+1],j] = coordinates(elements4(j,np.arange(1,4+1)),1)
        Sqy[np.arange[1,4+1],j] = coordinates(elements4(j,np.arange(1,4+1)),2)
        Sqc[np.arange[1,4+1],j] = u(elements4(j,np.arange(1,4+1)))
    
    if elements4.shape[1-1] > 20000.0 or str(edgeColorOrNot) == str('NoEdgeColor') == 1:
        h = patch(Sqx,Sqy,Sqc,'facecolor','interp','edgecolor','none')
    else:
        h = patch(Sqx,Sqy,Sqc,'facecolor','interp')
    
    hold('off')
    view(10,40)
    # title('Solution of the Problem')
    box('on')
    return h