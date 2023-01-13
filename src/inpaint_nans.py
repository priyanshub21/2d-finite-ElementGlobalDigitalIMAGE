import numpy as np
import numpy.matlib
    
def inpaint_nans(A = None,method = None): 
    # INPAINT_NANS: in-paints over nans in an array
# usage: B=INPAINT_NANS(A)          # default method
# usage: B=INPAINT_NANS(A,method)   # specify method used
    
    # Solves approximation to one of several pdes to
# interpolate and extrapolate holes in an array
    
    # arguments (input):
#   A - nxm array with some NaNs to be filled in
    
    #   method - (OPTIONAL) scalar numeric flag - specifies
#       which approach (or physical metaphor to use
#       for the interpolation.) All methods are capable
#       of extrapolation, some are better than others.
#       There are also speed differences, as well as
#       accuracy differences for smooth surfaces.
    
    #       methods {0,1,2} use a simple plate metaphor.
#       method  3 uses a better plate equation,
#                 but may be much slower and uses
#                 more memory.
#       method  4 uses a spring metaphor.
#       method  5 is an 8 neighbor average, with no
#                 rationale behind it compared to the
#                 other methods. I do not recommend
#                 its use.
    
    #       method == 0 --> (DEFAULT) see method 1, but
#         this method does not build as large of a
#         linear system in the case of only a few
#         NaNs in a large array.
#         Extrapolation behavior is linear.
    
    #       method == 1 --> simple approach, applies del^2
#         over the entire array, then drops those parts
#         of the array which do not have any contact with
#         NaNs. Uses a least squares approach, but it
#         does not modify known values.
#         In the case of small arrays, this method is
#         quite fast as it does very little extra work.
#         Extrapolation behavior is linear.
    
    #       method == 2 --> uses del^2, but solving a direct
#         linear system of equations for nan elements.
#         This method will be the fastest possible for
#         large systems since it uses the sparsest
#         possible system of equations. Not a least
#         squares approach, so it may be least robust
#         to noise on the boundaries of any holes.
#         This method will also be least able to
#         interpolate accurately for smooth surfaces.
#         Extrapolation behavior is linear.
    
    #         Note: method 2 has problems in 1-d, so this
#         method is disabled for vector inputs.
    
    #       method == 3 --+ See method 0, but uses del^4 for
#         the interpolating operator. This may result
#         in more accurate interpolations, at some cost
#         in speed.
    
    #       method == 4 --+ Uses a spring metaphor. Assumes
#         springs (with a nominal length of zero)
#         connect each node with every neighbor
#         (horizontally, vertically and diagonally)
#         Since each node tries to be like its neighbors,
#         extrapolation is as a constant function where
#         this is consistent with the neighboring nodes.
    
    #       method == 5 --+ See method 2, but use an average
#         of the 8 nearest neighbors to any element.
#         This method is NOT recommended for use.
    
    
    # arguments (output):
#   B - nxm array with NaNs replaced
    
    
    # Example:
#  [x,y] = meshgrid(0:.01:1);
#  z0 = exp(x+y);
#  znan = z0;
#  znan(20:50,40:70) = NaN;
#  znan(30:90,5:10) = NaN;
#  znan(70:75,40:90) = NaN;
    
    #  z = inpaint_nans(znan);
    
    
    # See also: griddata, interp1
    
    # Author: John D'Errico
# e-mail address: woodchips@rochester.rr.com
# Release: 2
# Release date: 4/15/06
    
    # I always need to know which elements are NaN,
# and what size the array is for any method
    n,m = A.shape
    A = A
    nm = n * m
    k = np.isnan(A)
    # list the nodes which are known, and which will
# be interpolated
    nan_list = find(k)
    known_list = find(not k )
    # how many nans overall
    nan_count = len(nan_list)
    # convert NaN indices to (r,c) form
# nan_list==find(k) are the unrolled (linear) indices
# (row,column) form
    nr,nc = ind2sub(np.array([n,m]),nan_list)
    # both forms of index in one array:
# column 1 == unrolled index
# column 2 == row index
# column 3 == column index
    nan_list = np.array([nan_list,nr,nc])
    # supply default method
    if (len(varargin) < 2) or len(method)==0:
        method = 0
    else:
        if not ismember(method,np.arange(0,5+1)) :
            raise Exception('If supplied, method must be one of: {0,1,2,3,4,5}.')
    
    # for different methods
    if 0 == method:
        # The same as method == 1, except only work on those
# elements which are NaN, or at least touch a NaN.
        # is it 1-d or 2-d?
        if (m == 1) or (n == 1):
            # really a 1-d case
            work_list = nan_list(:,1)
            work_list = unique(np.array([[work_list],[work_list - 1],[work_list + 1]]))
            work_list[work_list <= 1] = []
            work_list[work_list >= nm] = []
            nw = np.asarray(work_list).size
            u = np.transpose((np.arange(1,nw+1)))
            fda = sparse(np.matlib.repmat(u,1,3),bsxfun(plus,work_list,np.arange(- 1,1+1)),np.matlib.repmat(np.array([1,- 2,1]),nw,1),nw,nm)
        else:
            # a 2-d case
            # horizontal and vertical neighbors only
            talks_to = np.array([[- 1,0],[0,- 1],[1,0],[0,1]])
            neighbors_list = identify_neighbors(n,m,nan_list,talks_to)
            # list of all nodes we have identified
            all_list = np.array([[nan_list],[neighbors_list]])
            # generate sparse array with second partials on row
# variable for each element in either list, but only
# for those nodes which have a row index > 1 or < n
            L = find(np.logical_and((all_list(:,2) > 1),(all_list(:,2) < n)))
            nl = len(L)
            if nl > 0:
                fda = sparse(np.matlib.repmat(all_list(L,1),1,3),np.matlib.repmat(all_list(L,1),1,3) + np.matlib.repmat(np.array([- 1,0,1]),nl,1),np.matlib.repmat(np.array([1,- 2,1]),nl,1),nm,nm)
            else:
                fda = spalloc(n * m,n * m,all_list.shape[1-1] * 5)
            # 2nd partials on column index
            L = find(np.logical_and((all_list(:,3) > 1),(all_list(:,3) < m)))
            nl = len(L)
            if nl > 0:
                fda = fda + sparse(np.matlib.repmat(all_list(L,1),1,3),np.matlib.repmat(all_list(L,1),1,3) + np.matlib.repmat(np.array([- n,0,n]),nl,1),np.matlib.repmat(np.array([1,- 2,1]),nl,1),nm,nm)
        # eliminate knowns
        rhs = - fda(:,known_list) * A(known_list)
        k = find(np.any(fda(:,nan_list(:,1)),2))
        # and solve...
        B = A
        B[nan_list[:,1]] = np.linalg.solve(fda(k,nan_list(:,1)),rhs(k))
    else:
        if 1 == method:
            # least squares approach with del^2. Build system
# for every array element as an unknown, and then
# eliminate those which are knowns.
            # Build sparse matrix approximating del^2 for
# every element in A.
            # is it 1-d or 2-d?
            if (m == 1) or (n == 1):
                # a 1-d case
                u = np.transpose((np.arange(1,(nm - 2)+1)))
                fda = sparse(np.matlib.repmat(u,1,3),bsxfun(plus,u,np.arange(0,2+1)),np.matlib.repmat(np.array([1,- 2,1]),nm - 2,1),nm - 2,nm)
            else:
                # a 2-d case
                # Compute finite difference for second partials
# on row variable first
                i,j = ndgrid(np.arange(2,(n - 1)+1),np.arange(1,m+1))
                ind = i + (j - 1) * n
                np = (n - 2) * m
                fda = sparse(np.matlib.repmat(ind,1,3),np.array([ind - 1,ind,ind + 1]),np.matlib.repmat(np.array([1,- 2,1]),np,1),n * m,n * m)
                # now second partials on column variable
                i,j = ndgrid(np.arange(1,n+1),np.arange(2,(m - 1)+1))
                ind = i + (j - 1) * n
                np = n * (m - 2)
                fda = fda + sparse(np.matlib.repmat(ind,1,3),np.array([ind - n,ind,ind + n]),np.matlib.repmat(np.array([1,- 2,1]),np,1),nm,nm)
            # eliminate knowns
            rhs = - fda(:,known_list) * A(known_list)
            k = find(np.any(fda(:,nan_list),2))
            # and solve...
            B = A
            B[nan_list[:,1]] = np.linalg.solve(fda(k,nan_list(:,1)),rhs(k))
        else:
            if 2 == method:
                # Direct solve for del^2 BVP across holes
                # generate sparse array with second partials on row
# variable for each nan element, only for those nodes
# which have a row index > 1 or < n
                # is it 1-d or 2-d?
                if (m == 1) or (n == 1):
                    # really just a 1-d case
                    raise Exception('Method 2 has problems for vector input. Please use another method.')
                else:
                    # a 2-d case
                    L = find(np.logical_and((nan_list(:,2) > 1),(nan_list(:,2) < n)))
                    nl = len(L)
                    if nl > 0:
                        fda = sparse(np.matlib.repmat(nan_list(L,1),1,3),np.matlib.repmat(nan_list(L,1),1,3) + np.matlib.repmat(np.array([- 1,0,1]),nl,1),np.matlib.repmat(np.array([1,- 2,1]),nl,1),n * m,n * m)
                    else:
                        fda = spalloc(n * m,n * m,nan_list.shape[1-1] * 5)
                    # 2nd partials on column index
                    L = find(np.logical_and((nan_list(:,3) > 1),(nan_list(:,3) < m)))
                    nl = len(L)
                    if nl > 0:
                        fda = fda + sparse(np.matlib.repmat(nan_list(L,1),1,3),np.matlib.repmat(nan_list(L,1),1,3) + np.matlib.repmat(np.array([- n,0,n]),nl,1),np.matlib.repmat(np.array([1,- 2,1]),nl,1),n * m,n * m)
                    # fix boundary conditions at extreme corners
# of the array in case there were nans there
                    if ismember(1,nan_list(:,1)):
                        fda[1,np.array[[1,2,n + 1]]] = np.array([- 2,1,1])
                    if ismember(n,nan_list(:,1)):
                        fda[n,np.array[[n,n - 1,n + n]]] = np.array([- 2,1,1])
                    if ismember(nm - n + 1,nan_list(:,1)):
                        fda[nm - n + 1,np.array[[nm - n + 1,nm - n + 2,nm - n]]] = np.array([- 2,1,1])
                    if ismember(nm,nan_list(:,1)):
                        fda[nm,np.array[[nm,nm - 1,nm - n]]] = np.array([- 2,1,1])
                    # eliminate knowns
                    rhs = - fda(:,known_list) * A(known_list)
                    # and solve...
                    B = A
                    k = nan_list(:,1)
                    B[k] = np.linalg.solve(fda(k,k),rhs(k))
            else:
                if 3 == method:
                    # The same as method == 0, except uses del^4 as the
# interpolating operator.
                    # del^4 template of neighbors
                    talks_to = np.array([[- 2,0],[- 1,- 1],[- 1,0],[- 1,1],[0,- 2],[0,- 1],[0,1],[0,2],[1,- 1],[1,0],[1,1],[2,0]])
                    neighbors_list = identify_neighbors(n,m,nan_list,talks_to)
                    # list of all nodes we have identified
                    all_list = np.array([[nan_list],[neighbors_list]])
                    # generate sparse array with del^4, but only
# for those nodes which have a row & column index
# >= 3 or <= n-2
                    L = find(np.logical_and(np.logical_and(np.logical_and((all_list(:,2) >= 3),(all_list(:,2) <= (n - 2))),(all_list(:,3) >= 3)),(all_list(:,3) <= (m - 2))))
                    nl = len(L)
                    if nl > 0:
                        # do the entire template at once
                        fda = sparse(np.matlib.repmat(all_list(L,1),1,13),np.matlib.repmat(all_list(L,1),1,13) + np.matlib.repmat(np.array([- 2 * n,- n - 1,- n,- n + 1,- 2,- 1,0,1,2,n - 1,n,n + 1,2 * n]),nl,1),np.matlib.repmat(np.array([1,2,- 8,2,1,- 8,20,- 8,1,2,- 8,2,1]),nl,1),nm,nm)
                    else:
                        fda = spalloc(n * m,n * m,all_list.shape[1-1] * 5)
                    # on the boundaries, reduce the order around the edges
                    L = find(np.logical_or((np.logical_and(np.logical_and((np.logical_or((all_list(:,2) == 2),(all_list(:,2) == (n - 1)))),(all_list(:,3) >= 2)),(all_list(:,3) <= (m - 1)))),(np.logical_and(np.logical_and((np.logical_or((all_list(:,3) == 2),(all_list(:,3) == (m - 1)))),(all_list(:,2) >= 2)),(all_list(:,2) <= (n - 1))))))
                    nl = len(L)
                    if nl > 0:
                        fda = fda + sparse(np.matlib.repmat(all_list(L,1),1,5),np.matlib.repmat(all_list(L,1),1,5) + np.matlib.repmat(np.array([- n,- 1,0,+ 1,n]),nl,1),np.matlib.repmat(np.array([1,1,- 4,1,1]),nl,1),nm,nm)
                    L = find(np.logical_and(np.logical_and((np.logical_or((all_list(:,2) == 1),(all_list(:,2) == n))),(all_list(:,3) >= 2)),(all_list(:,3) <= (m - 1))))
                    nl = len(L)
                    if nl > 0:
                        fda = fda + sparse(np.matlib.repmat(all_list(L,1),1,3),np.matlib.repmat(all_list(L,1),1,3) + np.matlib.repmat(np.array([- n,0,n]),nl,1),np.matlib.repmat(np.array([1,- 2,1]),nl,1),nm,nm)
                    L = find(np.logical_and(np.logical_and((np.logical_or((all_list(:,3) == 1),(all_list(:,3) == m))),(all_list(:,2) >= 2)),(all_list(:,2) <= (n - 1))))
                    nl = len(L)
                    if nl > 0:
                        fda = fda + sparse(np.matlib.repmat(all_list(L,1),1,3),np.matlib.repmat(all_list(L,1),1,3) + np.matlib.repmat(np.array([- 1,0,1]),nl,1),np.matlib.repmat(np.array([1,- 2,1]),nl,1),nm,nm)
                    # eliminate knowns
                    rhs = - fda(:,known_list) * A(known_list)
                    k = find(np.any(fda(:,nan_list(:,1)),2))
                    # and solve...
                    B = A
                    B[nan_list[:,1]] = np.linalg.solve(fda(k,nan_list(:,1)),rhs(k))
                else:
                    if 4 == method:
                        # Spring analogy
# interpolating operator.
                        # list of all springs between a node and a horizontal
# or vertical neighbor
                        hv_list = np.array([[- 1,- 1,0],[1,1,0],[- n,0,- 1],[n,0,1]])
                        hv_springs = []
                        for i in np.arange(1,4+1).reshape(-1):
                            hvs = nan_list + np.matlib.repmat(hv_list(i,:),nan_count,1)
                            k = np.logical_and(np.logical_and(np.logical_and((hvs(:,2) >= 1),(hvs(:,2) <= n)),(hvs(:,3) >= 1)),(hvs(:,3) <= m))
                            hv_springs = np.array([[hv_springs],[np.array([nan_list(k,1),hvs(k,1)])]])
                        # delete replicate springs
                        hv_springs = unique(__builtint__.sorted(hv_springs,2),'rows')
                        # build sparse matrix of connections, springs
# connecting diagonal neighbors are weaker than
# the horizontal and vertical springs
                        nhv = hv_springs.shape[1-1]
                        springs = sparse(np.matlib.repmat(np.transpose((np.arange(1,nhv+1))),1,2),hv_springs,np.matlib.repmat(np.array([1,- 1]),nhv,1),nhv,nm)
                        # eliminate knowns
                        rhs = - springs(:,known_list) * A(known_list)
                        # and solve...
                        B = A
                        B[nan_list[:,1]] = np.linalg.solve(springs(:,nan_list(:,1)),rhs)
                    else:
                        if 5 == method:
                            # Average of 8 nearest neighbors
                            # generate sparse array to average 8 nearest neighbors
# for each nan element, be careful around edges
                            fda = spalloc(n * m,n * m,nan_list.shape[1-1] * 9)
                            # -1,-1
                            L = find(np.logical_and((nan_list(:,2) > 1),(nan_list(:,3) > 1)))
                            nl = len(L)
                            if nl > 0:
                                fda = fda + sparse(np.matlib.repmat(nan_list(L,1),1,2),np.matlib.repmat(nan_list(L,1),1,2) + np.matlib.repmat(np.array([- n - 1,0]),nl,1),np.matlib.repmat(np.array([1,- 1]),nl,1),n * m,n * m)
                            # 0,-1
                            L = find(nan_list(:,3) > 1)
                            nl = len(L)
                            if nl > 0:
                                fda = fda + sparse(np.matlib.repmat(nan_list(L,1),1,2),np.matlib.repmat(nan_list(L,1),1,2) + np.matlib.repmat(np.array([- n,0]),nl,1),np.matlib.repmat(np.array([1,- 1]),nl,1),n * m,n * m)
                            # +1,-1
                            L = find(np.logical_and((nan_list(:,2) < n),(nan_list(:,3) > 1)))
                            nl = len(L)
                            if nl > 0:
                                fda = fda + sparse(np.matlib.repmat(nan_list(L,1),1,2),np.matlib.repmat(nan_list(L,1),1,2) + np.matlib.repmat(np.array([- n + 1,0]),nl,1),np.matlib.repmat(np.array([1,- 1]),nl,1),n * m,n * m)
                            # -1,0
                            L = find(nan_list(:,2) > 1)
                            nl = len(L)
                            if nl > 0:
                                fda = fda + sparse(np.matlib.repmat(nan_list(L,1),1,2),np.matlib.repmat(nan_list(L,1),1,2) + np.matlib.repmat(np.array([- 1,0]),nl,1),np.matlib.repmat(np.array([1,- 1]),nl,1),n * m,n * m)
                            # +1,0
                            L = find(nan_list(:,2) < n)
                            nl = len(L)
                            if nl > 0:
                                fda = fda + sparse(np.matlib.repmat(nan_list(L,1),1,2),np.matlib.repmat(nan_list(L,1),1,2) + np.matlib.repmat(np.array([1,0]),nl,1),np.matlib.repmat(np.array([1,- 1]),nl,1),n * m,n * m)
                            # -1,+1
                            L = find(np.logical_and((nan_list(:,2) > 1),(nan_list(:,3) < m)))
                            nl = len(L)
                            if nl > 0:
                                fda = fda + sparse(np.matlib.repmat(nan_list(L,1),1,2),np.matlib.repmat(nan_list(L,1),1,2) + np.matlib.repmat(np.array([n - 1,0]),nl,1),np.matlib.repmat(np.array([1,- 1]),nl,1),n * m,n * m)
                            # 0,+1
                            L = find(nan_list(:,3) < m)
                            nl = len(L)
                            if nl > 0:
                                fda = fda + sparse(np.matlib.repmat(nan_list(L,1),1,2),np.matlib.repmat(nan_list(L,1),1,2) + np.matlib.repmat(np.array([n,0]),nl,1),np.matlib.repmat(np.array([1,- 1]),nl,1),n * m,n * m)
                            # +1,+1
                            L = find(np.logical_and((nan_list(:,2) < n),(nan_list(:,3) < m)))
                            nl = len(L)
                            if nl > 0:
                                fda = fda + sparse(np.matlib.repmat(nan_list(L,1),1,2),np.matlib.repmat(nan_list(L,1),1,2) + np.matlib.repmat(np.array([n + 1,0]),nl,1),np.matlib.repmat(np.array([1,- 1]),nl,1),n * m,n * m)
                            # eliminate knowns
                            rhs = - fda(:,known_list) * A(known_list)
                            # and solve...
                            B = A
                            k = nan_list(:,1)
                            B[k] = np.linalg.solve(fda(k,k),rhs(k))
    
    # all done, make sure that B is the same shape as
# A was when we came in.
    B = reshape(B,n,m)
    # ====================================================
#      end of main function
# ====================================================
# ====================================================
#      begin subfunctions
# ====================================================
    
def identify_neighbors(n = None,m = None,nan_list = None,talks_to = None): 
    # identify_neighbors: identifies all the neighbors of
#   those nodes in nan_list, not including the nans
#   themselves
    
    # arguments (input):
#  n,m - scalar - [n,m]=size(A), where A is the
#      array to be interpolated
#  nan_list - array - list of every nan element in A
#      nan_list(i,1) == linear index of i'th nan element
#      nan_list(i,2) == row index of i'th nan element
#      nan_list(i,3) == column index of i'th nan element
#  talks_to - px2 array - defines which nodes communicate
#      with each other, i.e., which nodes are neighbors.
    
    #      talks_to(i,1) - defines the offset in the row
#                      dimension of a neighbor
#      talks_to(i,2) - defines the offset in the column
#                      dimension of a neighbor
    
    #      For example, talks_to = [-1 0;0 -1;1 0;0 1]
#      means that each node talks only to its immediate
#      neighbors horizontally and vertically.
    
    # arguments(output):
#  neighbors_list - array - list of all neighbors of
#      all the nodes in nan_list
    
    if not len(nan_list)==0 :
        # use the definition of a neighbor in talks_to
        nan_count = nan_list.shape[1-1]
        talk_count = talks_to.shape[1-1]
        nn = np.zeros((nan_count * talk_count,2))
        j = np.array([1,nan_count])
        for i in np.arange(1,talk_count+1).reshape(-1):
            nn[np.arange[j[1],j[2]+1],:] = nan_list(:,np.arange(2,3+1)) + np.matlib.repmat(talks_to(i,:),nan_count,1)
            j = j + nan_count
        # drop those nodes which fall outside the bounds of the
# original array
        L = np.logical_or(np.logical_or(np.logical_or((nn(:,1) < 1),(nn(:,1) > n)),(nn(:,2) < 1)),(nn(:,2) > m))
        nn[L,:] = []
        # form the same format 3 column array as nan_list
        neighbors_list = np.array([sub2ind(np.array([n,m]),nn(:,1),nn(:,2)),nn])
        # delete replicates in the neighbors list
        neighbors_list = unique(neighbors_list,'rows')
        # and delete those which are also in the list of NaNs.
        neighbors_list = setdiff(neighbors_list,nan_list,'rows')
    else:
        neighbors_list = []
    
    return B