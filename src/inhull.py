import numpy as np
import warnings
import numpy.matlib
    
def inhull(testpts = None,xyz = None,tess = None,tol = None): 
    # inhull: tests if a set of points are inside a convex hull
# usage: in = inhull(testpts,xyz)
# usage: in = inhull(testpts,xyz,tess)
# usage: in = inhull(testpts,xyz,tess,tol)
    
    # arguments: (input)
#  testpts - nxp array to test, n data points, in p dimensions
#       If you have many points to test, it is most efficient to
#       call this function once with the entire set.
    
    #  xyz - mxp array of vertices of the convex hull, as used by
#       convhulln.
    
    #  tess - tessellation (or triangulation) generated by convhulln
#       If tess is left empty or not supplied, then it will be
#       generated.
    
    #  tol - (OPTIONAL) tolerance on the tests for inclusion in the
#       convex hull. You can think of tol as the distance a point
#       may possibly lie outside the hull, and still be perceived
#       as on the surface of the hull. Because of numerical slop
#       nothing can ever be done exactly here. I might guess a
#       semi-intelligent value of tol to be
    
    #         tol = 1.e-13*mean(abs(xyz(:)))
    
    #       In higher dimensions, the numerical issues of floating
#       point arithmetic will probably suggest a larger value
#       of tol.
    
    #       DEFAULT: tol = 0
    
    # arguments: (output)
#  in  - nx1 logical vector
#        in(i) == 1 --> the i'th point was inside the convex hull.
    
    # Example usage: The first point should be inside, the second out
    
    #  xy = randn(20,2);
#  tess = convhulln(xy);
#  testpoints = [ 0 0; 10 10];
#  in = inhull(testpoints,xy,tess)
    
    # in =
#      1
#      0
    
    # A non-zero count of the number of degenerate simplexes in the hull
# will generate a warning (in 4 or more dimensions.) This warning
# may be disabled off with the command:
    
    #   warning('off','inhull:degeneracy')
    
    # See also: convhull, convhulln, delaunay, delaunayn, tsearch, tsearchn
    
    # Author: John D'Errico
# e-mail: woodchips@rochester.rr.com
# Release: 3.0
# Release date: 10/26/06
    
    # get array sizes
# m points, p dimensions
    p = xyz.shape[2-1]
    n,c = testpts.shape
    if p != c:
        raise Exception('testpts and xyz must have the same number of columns')
    
    if p < 2:
        raise Exception('Points must lie in at least a 2-d space.')
    
    # was the convex hull supplied?
    if (len(varargin) < 3) or len(tess)==0:
        tess = convhulln(xyz)
    
    nt,c = tess.shape
    if c != p:
        raise Exception('tess array is incompatible with a dimension p space')
    
    # was tol supplied?
    if (len(varargin) < 4) or len(tol)==0:
        tol = 0
    
    # build normal vectors
    if 2 == p:
        # really simple for 2-d
        nrmls = (xyz(tess(:,1),:) - xyz(tess(:,2),:)) * np.array([[0,1],[- 1,0]])
        # Any degenerate edges?
        del_ = np.sqrt(np.sum(nrmls ** 2, 2-1))
        degenflag = (del_ < (np.amax(del_) * 10 * eps))
        if sum(degenflag) > 0:
            warnings.warn('inhull:degeneracy',np.array([num2str(sum(degenflag)),' degenerate edges identified in the convex hull']))
            # we need to delete those degenerate normal vectors
            nrmls[degenflag,:] = []
            nt = nrmls.shape[1-1]
    else:
        if 3 == p:
            # use vectorized cross product for 3-d
            ab = xyz(tess(:,1),:) - xyz(tess(:,2),:)
            ac = xyz(tess(:,1),:) - xyz(tess(:,3),:)
            nrmls = cross(ab,ac,2)
            degenflag = False(nt,1)
        else:
            # slightly more work in higher dimensions,
            nrmls = np.zeros((nt,p))
            degenflag = False(nt,1)
            for i in np.arange(1,nt+1).reshape(-1):
                # just in case of a degeneracy
# Note that bsxfun COULD be used in this line, but I have chosen to
# not do so to maintain compatibility. This code is still used by
# users of older releases.
#  nullsp = null(bsxfun(@minus,xyz(tess(i,2:end),:),xyz(tess(i,1),:)))';
                nullsp = np.transpose(null(xyz(tess(i,np.arange(2,end()+1)),:) - np.matlib.repmat(xyz(tess(i,1),:),p - 1,1)))
                if nullsp.shape[1-1] > 1:
                    degenflag[i] = True
                    nrmls[i,:] = NaN
                else:
                    nrmls[i,:] = nullsp
            if sum(degenflag) > 0:
                warnings.warn('inhull:degeneracy',np.array([num2str(sum(degenflag)),' degenerate simplexes identified in the convex hull']))
                # we need to delete those degenerate normal vectors
                nrmls[degenflag,:] = []
                nt = nrmls.shape[1-1]
    
    # scale normal vectors to unit length
    nrmllen = np.sqrt(np.sum(nrmls ** 2, 2-1))
    # again, bsxfun COULD be employed here...
#  nrmls = bsxfun(@times,nrmls,1./nrmllen);
    nrmls = np.multiply(nrmls,np.matlib.repmat(1.0 / nrmllen,1,p))
    # center point in the hull
    center = mean(xyz,1)
    # any point in the plane of each simplex in the convex hull
    a = xyz(tess(not degenflag ,1),:)
    # ensure the normals are pointing inwards
# this line too could employ bsxfun...
#  dp = sum(bsxfun(@minus,center,a).*nrmls,2);
    dp = np.sum(np.multiply((np.matlib.repmat(center,nt,1) - a),nrmls), 2-1)
    k = dp < 0
    nrmls[k,:] = - nrmls(k,:)
    # We want to test if:  dot((x - a),N) >= 0
# If so for all faces of the hull, then x is inside
# the hull. Change this to dot(x,N) >= dot(a,N)
    aN = np.sum(np.multiply(nrmls,a), 2-1)
    # test, be careful in case there are many points
    in_ = False(n,1)
    # if n is too large, we need to worry about the
# dot product grabbing huge chunks of memory.
    memblock = 1000000.0
    blocks = np.amax(1,int(np.floor(n / (memblock / nt))))
    aNr = np.matlib.repmat(aN,1,len(np.arange(1,n+blocks,blocks)))
    for i in np.arange(1,blocks+1).reshape(-1):
        j = np.arange(i,n+blocks,blocks)
        if aNr.shape[2-1] != len(j):
            aNr = np.matlib.repmat(aN,1,len(j))
        in_[j] = np.transpose(np.all((nrmls * np.transpose(testpts(j,:)) - aNr) >= - tol,1))
    
    return in