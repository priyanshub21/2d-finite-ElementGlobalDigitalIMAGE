import numpy as np
    
def gradientN(varargin = None): 
    #    [Vx, Vy, Vz] = gradientN(V,type)
    I,type_ = parseInputs(varargin[:])
    w,d = fspecialSeperable(type_)
    padOpt = 'symmetric'
    if True:
        if len(w) > 1:
            Ix = - imfilter(I,np.transpose(w),'same',padOpt)
            Iy = - imfilter(I,w,'same',padOpt)
        else:
            Ix = - I
            Iy = - I
        Ix = - imfilter(Ix,d,'conv','same',padOpt)
        Iy = - imfilter(Iy,np.transpose(d),'conv','same',padOpt)
    else:
        if len(w) > 1:
            Ix = - imfilter(I,np.transpose(w),'same',padOpt)
            Ix = - imfilter(Ix,reshape(w,1,1,[]),'same',padOpt)
            Iy = - imfilter(I,w,'same',padOpt)
            Iy = - imfilter(Iy,reshape(w,1,1,[]),'same',padOpt)
            Iz = - imfilter(I,w,'same',padOpt)
            Iz = - imfilter(Iz,np.transpose(w),'same',padOpt)
        else:
            Ix = - I
            Iy = - I
            Iz = - I
        Ix = - imfilter(Ix,d,'same',padOpt)
        Iy = - imfilter(Iy,np.transpose(d),'same',padOpt)
        Iz = - imfilter(Iz,reshape(d,1,1,[]),'same',padOpt)
        varargout[3] = Iz
    
    varargout[0] = Ix
    varargout[2] = Iy
    return varargout
    
    
def fspecialSeperable(type_ = None): 
    sumFlag01 = 1
    if 'fb' == type_.lower():
        w = 1
        d = np.array([1,- 1])
    else:
        if 'prewitt' == type_.lower():
            w = np.array([1,1,1])
            d = np.array([1,0,- 1])
        else:
            if 'sobel' == type_.lower():
                w = np.array([1,2,1])
                d = np.array([1,0,- 1])
            else:
                if 'scharr' == type_.lower():
                    w = np.array([3,10,3])
                    d = np.array([1,0,- 1])
                else:
                    if 'stencil' == type_.lower():
                        w = 1
                        d = np.array([- 1,8,0,- 8,1]) / 12
                        sumFlag01 = 0
                    else:
                        if strcmpi(type_(np.arange(1,7+1)),'optimal'):
                            sumFlag01 = 0
                            if 5 == str2double(type_(np.arange(8,end()+1))):
                                w = np.array([0.0376593171958139,0.249153396177344,0.426374573253684,0.249153396177344,0.0376593171958139])
                                d = np.array([0.109603762960256,0.27669098845555,0,- 0.27669098845555,- 0.109603762960256])
                            else:
                                if 7 == str2double(type_(np.arange(8,end()+1))):
                                    w = np.array([0.00541196740974425,0.0695905825057286,0.244559723794791,0.360875452579473,0.244559723794791,0.0695905825057286,0.00541196740974425])
                                    d = np.array([0.0194786630434688,0.123914729925,0.193554838845068,0,- 0.193554838845068,- 0.123914729925,- 0.0194786630434688])
                                else:
                                    if 9 == str2double(type_(np.arange(8,end()+1))):
                                        w = np.array([0.000737598362084457,0.0155298478667911,0.0902598182960924,0.234469365350285,0.318006740249494,0.234469365350285,0.0902598182960924,0.0155298478667911,0.000737598362084457])
                                        d = np.array([0.00303163095459785,0.0352414678518254,0.118879484725614,0.144382960330377,0,- 0.144382960330377,- 0.118879484725614,- 0.0352414678518254,- 0.00303163095459785])
                                    else:
                                        if 11 == str2double(type_(np.arange(8,end()+1))):
                                            w = np.array([9.73660046309032e-05,0.00304151665483548,0.026177919843567,0.103249150042758,0.223724696218449,0.287418702471518,0.223724696218449,0.103249150042758,0.026177919843567,0.00304151665483548,9.73660046309032e-05])
                                            d = np.array([0.000439905393637534,0.00808328124581577,0.0449782033510644,0.108830913335368,0.112870645898035,0,- 0.112870645898035,- 0.108830913335368,- 0.0449782033510644,- 0.00808328124581577,- 0.000439905393637534])
                                        else:
                                            if 13 == str2double(type_(np.arange(8,end()+1))):
                                                w = np.array([1.4954158497659e-05,0.000608044910058278,0.00691160005872185,0.0369033951658625,0.112173742521741,0.212488560211325,0.261799405947588,0.212488560211325,0.112173742521741,0.0369033951658625,0.00691160005872185,0.000608044910058278,1.4954158497659e-05])
                                                d = np.array([7.14573394720322e-05,0.00177958343198347,0.0138717416627181,0.0506584387298514,0.0970421434222207,0.0891267488305157,0,- 0.0891267488305157,- 0.0970421434222207,- 0.0506584387298514,- 0.0138717416627181,- 0.00177958343198347,- 7.14573394720322e-05])
                                            else:
                                                if 15 == str2double(type_(np.arange(8,end()+1))):
                                                    w = np.array([5.98860795150199e-06,0.000191445809195946,0.00198513495426173,0.01152483195803,0.0440836748027682,0.114823905477918,0.203925200191027,0.246919636397694,0.203925200191027,0.114823905477918,0.0440836748027682,0.01152483195803,0.00198513495426173,0.000191445809195946,5.98860795150199e-06])
                                                    d = np.array([2.68289243388421e-05,0.000526939396634328,0.00397615629561318,0.0177322281254572,0.0506284057542549,0.0879544206261874,0.078046778100399,0,- 0.078046778100399,- 0.0879544206261874,- 0.0506284057542549,- 0.0177322281254572,- 0.00397615629561318,- 0.000526939396634328,- 2.68289243388421e-05])
                                                else:
                                                    if 17 == str2double(type_(np.arange(8,end()+1))):
                                                        w = np.array([7.69392199800271e-07,3.43462804426615e-05,0.000460563932733138,0.0032856489733516,0.0153609347354413,0.0504415385506543,0.11790764558647,0.196238466137,0.232540172823407,0.196238466137,0.11790764558647,0.0504415385506543,0.0153609347354413,0.0032856489733516,0.000460563932733138,3.43462804426615e-05,7.69392199800271e-07])
                                                        d = np.array([3.73783393466673e-06,0.000104597678828706,0.00102618747708425,0.00569710286194269,0.0209031813657965,0.0513851162956773,0.0801001019942813,0.0666269955155201,0,- 0.0666269955155201,- 0.0801001019942813,- 0.0513851162956773,- 0.0209031813657965,- 0.00569710286194269,- 0.00102618747708425,- 0.000104597678828706,- 3.73783393466673e-06])
                                                    else:
                                                        if 19 == str2double(type_(np.arange(8,end()+1))):
                                                            w = np.array([1.52643475621529e-07,7.90799631627647e-06,0.000122710209565798,0.000998814542578885,0.00531312830085921,0.020277562358422,0.0572293928770952,0.120087013816114,0.187334615968446,0.217257402574254,0.187334615968446,0.120087013816114,0.0572293928770952,0.020277562358422,0.00531312830085921,0.000998814542578885,0.000122710209565798,7.90799631627647e-06,1.52643475621529e-07])
                                                            d = np.array([7.63266760902856e-07,2.52762743204803e-05,0.000290470367850947,0.00185854511355811,0.00795173988567152,0.0240580568530862,0.0508916097992498,0.0712071454435055,0.0555263098864427,0,- 0.0555263098864427,- 0.0712071454435055,- 0.0508916097992498,- 0.0240580568530862,- 0.00795173988567152,- 0.00185854511355811,- 0.000290470367850947,- 2.52762743204803e-05,- 7.63266760902856e-07])
                                                        else:
                                                            w = np.array([0.229878817299031,0.540242365401939,0.229878817299031])
                                                            d = np.array([0.425286806086887,0,- 0.425286806086887])
                        else:
                            w = 1
                            d = np.array([- 1,0,1])
    
    w = w / sum(np.abs(w))
    if sumFlag01:
        d = d / sum(np.abs(d))
    
    return w,d
    
    
def parseInputs(varargin = None): 
    if len(varargin) < 2:
        varargin[2] = 'prewitt'
    
    I = double(varargin[0])
    type_ = varargin[2]
    return I,type_
    
    return varargout