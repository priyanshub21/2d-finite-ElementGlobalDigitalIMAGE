import numpy as np
import warnings
import numpy.matlib
    
def gridfit(x = None,y = None,z = None,xnodes = None,ynodes = None,varargin = None): 
    # gridfit: estimates a surface on a 2d grid, based on scattered data
#          Replicates are allowed. All methods extrapolate to the grid
#          boundaries. Gridfit uses a modified ridge estimator to
#          generate the surface, where the bias is toward smoothness.
    
    #          Gridfit is not an interpolant. Its goal is a smooth surface
#          that approximates your data, but allows you to control the
#          amount of smoothing.
    
    # usage #1: zgrid = gridfit(x,y,z,xnodes,ynodes);
# usage #2: [zgrid,xgrid,ygrid] = gridfit(x,y,z,xnodes,ynodes);
# usage #3: zgrid = gridfit(x,y,z,xnodes,ynodes,prop,val,prop,val,...);
    
    # Arguments: (input)
#  x,y,z - vectors of equal lengths, containing arbitrary scattered data
#          The only constraint on x and y is they cannot ALL fall on a
#          single line in the x-y plane. Replicate points will be treated
#          in a least squares sense.
    
    #          ANY points containing a NaN are ignored in the estimation
    
    #  xnodes - vector defining the nodes in the grid in the independent
#          variable (x). xnodes need not be equally spaced. xnodes
#          must completely span the data. If they do not, then the
#          'extend' property is applied, adjusting the first and last
#          nodes to be extended as necessary. See below for a complete
#          description of the 'extend' property.
    
    #          If xnodes is a scalar integer, then it specifies the number
#          of equally spaced nodes between the min and max of the data.
    
    #  ynodes - vector defining the nodes in the grid in the independent
#          variable (y). ynodes need not be equally spaced.
    
    #          If ynodes is a scalar integer, then it specifies the number
#          of equally spaced nodes between the min and max of the data.
    
    #          Also see the extend property.
    
    #  Additional arguments follow in the form of property/value pairs.
#  Valid properties are:
#    'smoothness', 'interp', 'regularizer', 'solver', 'maxiter'
#    'extend', 'tilesize', 'overlap'
    
    #  Any UNAMBIGUOUS shortening (even down to a single letter) is
#  valid for property names. All properties have default values,
#  chosen (I hope) to give a reasonable result out of the box.
    
    #   'smoothness' - scalar or vector of length 2 - determines the
#          eventual smoothness of the estimated surface. A larger
#          value here means the surface will be smoother. Smoothness
#          must be a non-negative real number.
    
    #          If this parameter is a vector of length 2, then it defines
#          the relative smoothing to be associated with the x and y
#          variables. This allows the user to apply a different amount
#          of smoothing in the x dimension compared to the y dimension.
    
    #          Note: the problem is normalized in advance so that a
#          smoothness of 1 MAY generate reasonable results. If you
#          find the result is too smooth, then use a smaller value
#          for this parameter. Likewise, bumpy surfaces suggest use
#          of a larger value. (Sometimes, use of an iterative solver
#          with too small a limit on the maximum number of iterations
#          will result in non-convergence.)
    
    #          DEFAULT: 1
    
    
    #   'interp' - character, denotes the interpolation scheme used
#          to interpolate the data.
    
    #          DEFAULT: 'triangle'
    
    #          'bilinear' - use bilinear interpolation within the grid
#                     (also known as tensor product linear interpolation)
    
    #          'triangle' - split each cell in the grid into a triangle,
#                     then linear interpolation inside each triangle
    
    #          'nearest' - nearest neighbor interpolation. This will
#                     rarely be a good choice, but I included it
#                     as an option for completeness.
    
    
    #   'regularizer' - character flag, denotes the regularization
#          paradignm to be used. There are currently three options.
    
    #          DEFAULT: 'gradient'
    
    #          'diffusion' or 'laplacian' - uses a finite difference
#              approximation to the Laplacian operator (i.e, del^2).
    
    #              We can think of the surface as a plate, wherein the
#              bending rigidity of the plate is specified by the user
#              as a number relative to the importance of fidelity to
#              the data. A stiffer plate will result in a smoother
#              surface overall, but fit the data less well. I've
#              modeled a simple plate using the Laplacian, del^2. (A
#              projected enhancement is to do a better job with the
#              plate equations.)
    
    #              We can also view the regularizer as a diffusion problem,
#              where the relative thermal conductivity is supplied.
#              Here interpolation is seen as a problem of finding the
#              steady temperature profile in an object, given a set of
#              points held at a fixed temperature. Extrapolation will
#              be linear. Both paradigms are appropriate for a Laplacian
#              regularizer.
    
    #          'gradient' - attempts to ensure the gradient is as smooth
#              as possible everywhere. Its subtly different from the
#              'diffusion' option, in that here the directional
#              derivatives are biased to be smooth across cell
#              boundaries in the grid.
    
    #              The gradient option uncouples the terms in the Laplacian.
#              Think of it as two coupled PDEs instead of one PDE. Why
#              are they different at all? The terms in the Laplacian
#              can balance each other.
    
    #          'springs' - uses a spring model connecting nodes to each
#              other, as well as connecting data points to the nodes
#              in the grid. This choice will cause any extrapolation
#              to be as constant as possible.
    
    #              Here the smoothing parameter is the relative stiffness
#              of the springs connecting the nodes to each other compared
#              to the stiffness of a spting connecting the lattice to
#              each data point. Since all springs have a rest length
#              (length at which the spring has zero potential energy)
#              of zero, any extrapolation will be minimized.
    
    #          Note: The 'springs' regularizer tends to drag the surface
#          towards the mean of all the data, so too large a smoothing
#          parameter may be a problem.
    
    
    #   'solver' - character flag - denotes the solver used for the
#          resulting linear system. Different solvers will have
#          different solution times depending upon the specific
#          problem to be solved. Up to a certain size grid, the
#          direct \ solver will often be speedy, until memory
#          swaps causes problems.
    
    #          What solver should you use? Problems with a significant
#          amount of extrapolation should avoid lsqr. \ may be
#          best numerically for small smoothnesss parameters and
#          high extents of extrapolation.
    
    #          Large numbers of points will slow down the direct
#          \, but when applied to the normal equations, \ can be
#          quite fast. Since the equations generated by these
#          methods will tend to be well conditioned, the normal
#          equations are not a bad choice of method to use. Beware
#          when a small smoothing parameter is used, since this will
#          make the equations less well conditioned.
    
    #          DEFAULT: 'normal'
    
    #          '\' - uses matlab's backslash operator to solve the sparse
#                     system. 'backslash' is an alternate name.
    
    #          'symmlq' - uses matlab's iterative symmlq solver
    
    #          'lsqr' - uses matlab's iterative lsqr solver
    
    #          'normal' - uses \ to solve the normal equations.
    
    
    #   'maxiter' - only applies to iterative solvers - defines the
#          maximum number of iterations for an iterative solver
    
    #          DEFAULT: min(10000,length(xnodes)*length(ynodes))
    
    
    #   'extend' - character flag - controls whether the first and last
#          nodes in each dimension are allowed to be adjusted to
#          bound the data, and whether the user will be warned if
#          this was deemed necessary to happen.
    
    #          DEFAULT: 'warning'
    
    #          'warning' - Adjust the first and/or last node in
#                     x or y if the nodes do not FULLY contain
#                     the data. Issue a warning message to this
#                     effect, telling the amount of adjustment
#                     applied.
    
    #          'never'  - Issue an error message when the nodes do
#                     not absolutely contain the data.
    
    #          'always' - automatically adjust the first and last
#                     nodes in each dimension if necessary.
#                     No warning is given when this option is set.
    
    
    #   'tilesize' - grids which are simply too large to solve for
#          in one single estimation step can be built as a set
#          of tiles. For example, a 1000x1000 grid will require
#          the estimation of 1e6 unknowns. This is likely to
#          require more memory (and time) than you have available.
#          But if your data is dense enough, then you can model
#          it locally using smaller tiles of the grid.
    
    #          My recommendation for a reasonable tilesize is
#          roughly 100 to 200. Tiles of this size take only
#          a few seconds to solve normally, so the entire grid
#          can be modeled in a finite amount of time. The minimum
#          tilesize can never be less than 3, although even this
#          size tile is so small as to be ridiculous.
    
    #          If your data is so sparse than some tiles contain
#          insufficient data to model, then those tiles will
#          be left as NaNs.
    
    #          DEFAULT: inf
    
    
    #   'overlap' - Tiles in a grid have some overlap, so they
#          can minimize any problems along the edge of a tile.
#          In this overlapped region, the grid is built using a
#          bi-linear combination of the overlapping tiles.
    
    #          The overlap is specified as a fraction of the tile
#          size, so an overlap of 0.20 means there will be a 20#
#          overlap of successive tiles. I do allow a zero overlap,
#          but it must be no more than 1/2.
    
    #          0 <= overlap <= 0.5
    
    #          Overlap is ignored if the tilesize is greater than the
#          number of nodes in both directions.
    
    #          DEFAULT: 0.20
    
    
    #   'autoscale' - Some data may have widely different scales on
#          the respective x and y axes. If this happens, then
#          the regularization may experience difficulties.
    
    #          autoscale = 'on' will cause gridfit to scale the x
#          and y node intervals to a unit length. This should
#          improve the regularization procedure. The scaling is
#          purely internal.
    
    #          autoscale = 'off' will disable automatic scaling
    
    #          DEFAULT: 'on'
    
    
    # Arguments: (output)
#  zgrid   - (nx,ny) array containing the fitted surface
    
    #  xgrid, ygrid - as returned by meshgrid(xnodes,ynodes)
    
    
    # Speed considerations:
#  Remember that gridfit must solve a LARGE system of linear
#  equations. There will be as many unknowns as the total
#  number of nodes in the final lattice. While these equations
#  may be sparse, solving a system of 10000 equations may take
#  a second or so. Very large problems may benefit from the
#  iterative solvers or from tiling.
    
    
    # Example usage:
    
    #  x = rand(100,1);
#  y = rand(100,1);
#  z = exp(x+2*y);
#  xnodes = 0:.1:1;
#  ynodes = 0:.1:1;
    
    #  g = gridfit(x,y,z,xnodes,ynodes);
    
    # Note: this is equivalent to the following call:
    
    #  g = gridfit(x,y,z,xnodes,ynodes, ...
#              'smooth',1, ...
#              'interp','triangle', ...
#              'solver','normal', ...
#              'regularizer','gradient', ...
#              'extend','warning', ...
#              'tilesize',inf);
    
    
    # Author: John D'Errico
# e-mail address: woodchips@rochester.rr.com
# Release: 2.0
# Release date: 5/23/06
    
    # set defaults
    params.smoothness = 1
    params.interp = 'triangle'
    params.regularizer = 'gradient'
    params.solver = 'backslash'
    params.maxiter = []
    params.extend = 'warning'
    params.tilesize = inf
    params.overlap = 0.2
    params.mask = []
    params.autoscale = 'on'
    params.xscale = 1
    params.yscale = 1
    # was the params struct supplied?
    if not len(varargin)==0 :
        if isstruct(varargin[0]):
            # params is only supplied if its a call from tiled_gridfit
            params = varargin[0]
            if len(varargin) > 1:
                # check for any overrides
                params = parse_pv_pairs(params,varargin[np.arange(2,end()+1)])
        else:
            # check for any overrides of the defaults
            params = parse_pv_pairs(params,varargin)
    
    # check the parameters for acceptability
    params = check_params(params)
    # ensure all of x,y,z,xnodes,ynodes are column vectors,
# also drop any NaN data
    x = x
    y = y
    z = z
    k = np.logical_or(np.logical_or(np.isnan(x),np.isnan(y)),np.isnan(z))
    if np.any(k):
        x[k] = []
        y[k] = []
        z[k] = []
    
    xmin = np.amin(x)
    xmax = np.amax(x)
    ymin = np.amin(y)
    ymax = np.amax(y)
    # did they supply a scalar for the nodes?
    if len(xnodes) == 1:
        xnodes = np.transpose(np.linspace(xmin,xmax,xnodes))
        xnodes[end()] = xmax
    
    if len(ynodes) == 1:
        ynodes = np.transpose(np.linspace(ymin,ymax,ynodes))
        ynodes[end()] = ymax
    
    xnodes = xnodes
    ynodes = ynodes
    dx = np.diff(xnodes)
    dy = np.diff(ynodes)
    nx = len(xnodes)
    ny = len(ynodes)
    ngrid = nx * ny
    # set the scaling if autoscale was on
    if strcmpi(params.autoscale,'on'):
        params.xscale = mean(dx)
        params.yscale = mean(dy)
        params.autoscale = 'off'
    
    # check to see if any tiling is necessary
    if (params.tilesize < np.amax(nx,ny)):
        # split it into smaller tiles. compute zgrid and ygrid
# at the very end if requested
        zgrid = tiled_gridfit(x,y,z,xnodes,ynodes,params)
    else:
        # its a single tile.
        # mask must be either an empty array, or a boolean
# aray of the same size as the final grid.
        nmask = params.mask.shape
        if not len(params.mask)==0  and ((nmask(2) != nx) or (nmask(1) != ny)):
            if ((nmask(2) == ny) or (nmask(1) == nx)):
                raise Exception('Mask array is probably transposed from proper orientation.')
            else:
                raise Exception('Mask array must be the same size as the final grid.')
        if not len(params.mask)==0 :
            params.maskflag = 1
        else:
            params.maskflag = 0
        # default for maxiter?
        if len(params.maxiter)==0:
            params.maxiter = np.amin(10000,nx * ny)
        # check lengths of the data
        n = len(x)
        if (len(y) != n) or (len(z) != n):
            raise Exception('Data vectors are incompatible in size.')
        if n < 3:
            raise Exception('Insufficient data for surface estimation.')
        # verify the nodes are distinct
        if np.any(np.diff(xnodes) <= 0) or np.any(np.diff(ynodes) <= 0):
            raise Exception('xnodes and ynodes must be monotone increasing')
        # do we need to tweak the first or last node in x or y?
        if xmin < xnodes(1):
            if 'always' == params.extend:
                xnodes[1] = xmin
            else:
                if 'warning' == params.extend:
                    warnings.warn('GRIDFIT:extend',np.array(['xnodes(1) was decreased by: ',num2str(xnodes(1) - xmin),', new node = ',num2str(xmin)]))
                    xnodes[1] = xmin
                else:
                    if 'never' == params.extend:
                        raise Exception(np.array(['Some x (',num2str(xmin),') falls below xnodes(1) by: ',num2str(xnodes(1) - xmin)]))
        if xmax > xnodes(end()):
            if 'always' == params.extend:
                xnodes[end()] = xmax
            else:
                if 'warning' == params.extend:
                    warnings.warn('GRIDFIT:extend',np.array(['xnodes(end) was increased by: ',num2str(xmax - xnodes(end())),', new node = ',num2str(xmax)]))
                    xnodes[end()] = xmax
                else:
                    if 'never' == params.extend:
                        raise Exception(np.array(['Some x (',num2str(xmax),') falls above xnodes(end) by: ',num2str(xmax - xnodes(end()))]))
        if ymin < ynodes(1):
            if 'always' == params.extend:
                ynodes[1] = ymin
            else:
                if 'warning' == params.extend:
                    warnings.warn('GRIDFIT:extend',np.array(['ynodes(1) was decreased by: ',num2str(ynodes(1) - ymin),', new node = ',num2str(ymin)]))
                    ynodes[1] = ymin
                else:
                    if 'never' == params.extend:
                        raise Exception(np.array(['Some y (',num2str(ymin),') falls below ynodes(1) by: ',num2str(ynodes(1) - ymin)]))
        if ymax > ynodes(end()):
            if 'always' == params.extend:
                ynodes[end()] = ymax
            else:
                if 'warning' == params.extend:
                    warnings.warn('GRIDFIT:extend',np.array(['ynodes(end) was increased by: ',num2str(ymax - ynodes(end())),', new node = ',num2str(ymax)]))
                    ynodes[end()] = ymax
                else:
                    if 'never' == params.extend:
                        raise Exception(np.array(['Some y (',num2str(ymax),') falls above ynodes(end) by: ',num2str(ymax - ynodes(end()))]))
        # determine which cell in the array each point lies in
        junk,indx = histc(x,xnodes)
        junk,indy = histc(y,ynodes)
        # any point falling at the last node is taken to be
# inside the last cell in x or y.
        k = (indx == nx)
        indx[k] = indx(k) - 1
        k = (indy == ny)
        indy[k] = indy(k) - 1
        ind = indy + ny * (indx - 1)
        # Do we have a mask to apply?
        if params.maskflag:
            # if we do, then we need to ensure that every
# cell with at least one data point also has at
# least all of its corners unmasked.
            params.mask[ind] = 1
            params.mask[ind + 1] = 1
            params.mask[ind + ny] = 1
            params.mask[ind + ny + 1] = 1
        # interpolation equations for each point
        tx = np.amin(1,np.amax(0,(x - xnodes(indx)) / dx(indx)))
        ty = np.amin(1,np.amax(0,(y - ynodes(indy)) / dy(indy)))
        # Future enhancement: add cubic interpolant
        if 'triangle' == params.interp:
            # linear interpolation inside each triangle
            k = (tx > ty)
            L = np.ones((n,1))
            L[k] = ny
            t1 = np.amin(tx,ty)
            t2 = np.amax(tx,ty)
            A = sparse(np.matlib.repmat(np.transpose((np.arange(1,n+1))),1,3),np.array([ind,ind + ny + 1,ind + L]),np.array([1 - t2,t1,t2 - t1]),n,ngrid)
        else:
            if 'nearest' == params.interp:
                # nearest neighbor interpolation in a cell
                k = np.round(1 - ty) + np.round(1 - tx) * ny
                A = sparse(np.transpose((np.arange(1,n+1))),ind + k,np.ones((n,1)),n,ngrid)
            else:
                if 'bilinear' == params.interp:
                    # bilinear interpolation in a cell
                    A = sparse(np.matlib.repmat(np.transpose((np.arange(1,n+1))),1,4),np.array([ind,ind + 1,ind + ny,ind + ny + 1]),np.array([np.multiply((1 - tx),(1 - ty)),np.multiply((1 - tx),ty),np.multiply(tx,(1 - ty)),np.multiply(tx,ty)]),n,ngrid)
        rhs = z
        # do we have relative smoothing parameters?
        if np.asarray(params.smoothness).size == 1:
            # it was scalar, so treat both dimensions equally
            smoothparam = params.smoothness
            xyRelativeStiffness = np.array([[1],[1]])
        else:
            # It was a vector, so anisotropy reigns.
# I've already checked that the vector was of length 2
            smoothparam = np.sqrt(np.prod(params.smoothness))
            xyRelativeStiffness = params.smoothness / smoothparam
        # Build regularizer. Add del^4 regularizer one day.
        if 'springs' == params.regularizer:
            # zero "rest length" springs
            i,j = np.meshgrid(np.arange(1,nx+1),np.arange(1,(ny - 1)+1))
            ind = j + ny * (i - 1)
            m = nx * (ny - 1)
            stiffness = 1.0 / (dy / params.yscale)
            Areg = sparse(np.matlib.repmat(np.transpose((np.arange(1,m+1))),1,2),np.array([ind,ind + 1]),xyRelativeStiffness(2) * stiffness(j) * np.array([- 1,1]),m,ngrid)
            i,j = np.meshgrid(np.arange(1,(nx - 1)+1),np.arange(1,ny+1))
            ind = j + ny * (i - 1)
            m = (nx - 1) * ny
            stiffness = 1.0 / (dx / params.xscale)
            Areg = np.array([[Areg],[sparse(np.matlib.repmat(np.transpose((np.arange(1,m+1))),1,2),np.array([ind,ind + ny]),xyRelativeStiffness(1) * stiffness(i) * np.array([- 1,1]),m,ngrid)]])
            i,j = np.meshgrid(np.arange(1,(nx - 1)+1),np.arange(1,(ny - 1)+1))
            ind = j + ny * (i - 1)
            m = (nx - 1) * (ny - 1)
            stiffness = 1.0 / np.sqrt((dx(i) / params.xscale / xyRelativeStiffness(1)) ** 2 + (dy(j) / params.yscale / xyRelativeStiffness(2)) ** 2)
            Areg = np.array([[Areg],[sparse(np.matlib.repmat(np.transpose((np.arange(1,m+1))),1,2),np.array([ind,ind + ny + 1]),stiffness * np.array([- 1,1]),m,ngrid)]])
            Areg = np.array([[Areg],[sparse(np.matlib.repmat(np.transpose((np.arange(1,m+1))),1,2),np.array([ind + 1,ind + ny]),stiffness * np.array([- 1,1]),m,ngrid)]])
        else:
            if np.array(['diffusion','laplacian']) == params.regularizer:
                # thermal diffusion using Laplacian (del^2)
                i,j = np.meshgrid(np.arange(1,nx+1),np.arange(2,(ny - 1)+1))
                ind = j + ny * (i - 1)
                dy1 = dy(j - 1) / params.yscale
                dy2 = dy(j) / params.yscale
                Areg = sparse(np.matlib.repmat(ind,1,3),np.array([ind - 1,ind,ind + 1]),xyRelativeStiffness(2) * np.array([- 2.0 / (np.multiply(dy1,(dy1 + dy2))),2.0 / (np.multiply(dy1,dy2)),- 2.0 / (np.multiply(dy2,(dy1 + dy2)))]),ngrid,ngrid)
                i,j = np.meshgrid(np.arange(2,(nx - 1)+1),np.arange(1,ny+1))
                ind = j + ny * (i - 1)
                dx1 = dx(i - 1) / params.xscale
                dx2 = dx(i) / params.xscale
                Areg = Areg + sparse(np.matlib.repmat(ind,1,3),np.array([ind - ny,ind,ind + ny]),xyRelativeStiffness(1) * np.array([- 2.0 / (np.multiply(dx1,(dx1 + dx2))),2.0 / (np.multiply(dx1,dx2)),- 2.0 / (np.multiply(dx2,(dx1 + dx2)))]),ngrid,ngrid)
            else:
                if 'gradient' == params.regularizer:
                    # Subtly different from the Laplacian. A point for future
# enhancement is to do it better for the triangle interpolation
# case.
                    i,j = np.meshgrid(np.arange(1,nx+1),np.arange(2,(ny - 1)+1))
                    ind = j + ny * (i - 1)
                    dy1 = dy(j - 1) / params.yscale
                    dy2 = dy(j) / params.yscale
                    Areg = sparse(np.matlib.repmat(ind,1,3),np.array([ind - 1,ind,ind + 1]),xyRelativeStiffness(2) * np.array([- 2.0 / (np.multiply(dy1,(dy1 + dy2))),2.0 / (np.multiply(dy1,dy2)),- 2.0 / (np.multiply(dy2,(dy1 + dy2)))]),ngrid,ngrid)
                    i,j = np.meshgrid(np.arange(2,(nx - 1)+1),np.arange(1,ny+1))
                    ind = j + ny * (i - 1)
                    dx1 = dx(i - 1) / params.xscale
                    dx2 = dx(i) / params.xscale
                    Areg = np.array([[Areg],[sparse(np.matlib.repmat(ind,1,3),np.array([ind - ny,ind,ind + ny]),xyRelativeStiffness(1) * np.array([- 2.0 / (np.multiply(dx1,(dx1 + dx2))),2.0 / (np.multiply(dx1,dx2)),- 2.0 / (np.multiply(dx2,(dx1 + dx2)))]),ngrid,ngrid)]])
        nreg = Areg.shape[1-1]
        # Append the regularizer to the interpolation equations,
# scaling the problem first. Use the 1-norm for speed.
        NA = norm(A,1)
        NR = norm(Areg,1)
        A = np.array([[A],[Areg * (smoothparam * NA / NR)]])
        rhs = np.array([[rhs],[np.zeros((nreg,1))]])
        # do we have a mask to apply?
        if params.maskflag:
            unmasked = find(params.mask)
        # solve the full system, with regularizer attached
        if np.array(['\','backslash']) == params.solver:
            if params.maskflag:
                # there is a mask to use
                zgrid = np.full([ny,nx],np.nan)
                zgrid[unmasked] = np.linalg.solve(A(:,unmasked),rhs)
            else:
                # no mask
                zgrid = reshape(np.linalg.solve(A,rhs),ny,nx)
        else:
            if 'normal' == params.solver:
                # The normal equations, solved with \. Can be faster
# for huge numbers of data points, but reasonably
# sized grids. The regularizer makes A well conditioned
# so the normal equations are not a terribly bad thing
# here.
                if params.maskflag:
                    # there is a mask to use
                    Aunmasked = A(:,unmasked)
                    zgrid = np.full([ny,nx],np.nan)
                    zgrid[unmasked] = np.linalg.solve((np.transpose(Aunmasked) * Aunmasked),(np.transpose(Aunmasked) * rhs))
                else:
                    zgrid = reshape(np.linalg.solve((np.transpose(A) * A),(np.transpose(A) * rhs)),ny,nx)
            else:
                if 'symmlq' == params.solver:
                    # iterative solver - symmlq - requires a symmetric matrix,
# so use it to solve the normal equations. No preconditioner.
                    tol = np.abs(np.amax(z) - np.amin(z)) * 1e-13
                    if params.maskflag:
                        # there is a mask to use
                        zgrid = np.full([ny,nx],np.nan)
                        zgrid[unmasked],flag = symmlq(np.transpose(A(:,unmasked)) * A(:,unmasked),np.transpose(A(:,unmasked)) * rhs,tol,params.maxiter)
                    else:
                        zgrid,flag = symmlq(np.transpose(A) * A,np.transpose(A) * rhs,tol,params.maxiter)
                        zgrid = reshape(zgrid,ny,nx)
                    # display a warning if convergence problems
                    if 0 == flag:
                        # no problems with convergence
                        pass
                    else:
                        if 1 == flag:
                            # SYMMLQ iterated MAXIT times but did not converge.
                            warnings.warn('GRIDFIT:solver',np.array(['Symmlq performed ',num2str(params.maxiter),' iterations but did not converge.']))
                        else:
                            if 3 == flag:
                                # SYMMLQ stagnated, successive iterates were the same
                                warnings.warn('GRIDFIT:solver','Symmlq stagnated without apparent convergence.')
                            else:
                                warnings.warn('GRIDFIT:solver',np.array(['One of the scalar quantities calculated in',' symmlq was too small or too large to continue computing.']))
                else:
                    if 'lsqr' == params.solver:
                        # iterative solver - lsqr. No preconditioner here.
                        tol = np.abs(np.amax(z) - np.amin(z)) * 1e-13
                        if params.maskflag:
                            # there is a mask to use
                            zgrid = np.full([ny,nx],np.nan)
                            zgrid[unmasked],flag = lsqr(A(:,unmasked),rhs,tol,params.maxiter)
                        else:
                            zgrid,flag = lsqr(A,rhs,tol,params.maxiter)
                            zgrid = reshape(zgrid,ny,nx)
                        # display a warning if convergence problems
                        if 0 == flag:
                            # no problems with convergence
                            pass
                        else:
                            if 1 == flag:
                                # lsqr iterated MAXIT times but did not converge.
                                warnings.warn('GRIDFIT:solver',np.array(['Lsqr performed ',num2str(params.maxiter),' iterations but did not converge.']))
                            else:
                                if 3 == flag:
                                    # lsqr stagnated, successive iterates were the same
                                    warnings.warn('GRIDFIT:solver','Lsqr stagnated without apparent convergence.')
                                else:
                                    if 4 == flag:
                                        warnings.warn('GRIDFIT:solver',np.array(['One of the scalar quantities calculated in',' LSQR was too small or too large to continue computing.']))
    
    # only generate xgrid and ygrid if requested.
    if nargout > 1:
        xgrid,ygrid = np.meshgrid(xnodes,ynodes)
    
    # ============================================
# End of main function - gridfit
# ============================================
    
    # ============================================
# subfunction - parse_pv_pairs
# ============================================
    
def parse_pv_pairs(params = None,pv_pairs = None): 
    # parse_pv_pairs: parses sets of property value pairs, allows defaults
# usage: params=parse_pv_pairs(default_params,pv_pairs)
    
    # arguments: (input)
#  default_params - structure, with one field for every potential
#             property/value pair. Each field will contain the default
#             value for that property. If no default is supplied for a
#             given property, then that field must be empty.
    
    #  pv_array - cell array of property/value pairs.
#             Case is ignored when comparing properties to the list
#             of field names. Also, any unambiguous shortening of a
#             field/property name is allowed.
    
    # arguments: (output)
#  params   - parameter struct that reflects any updated property/value
#             pairs in the pv_array.
    
    # Example usage:
# First, set default values for the parameters. Assume we
# have four parameters that we wish to use optionally in
# the function examplefun.
    
    #  - 'viscosity', which will have a default value of 1
#  - 'volume', which will default to 1
#  - 'pie' - which will have default value 3.141592653589793
#  - 'description' - a text field, left empty by default
    
    # The first argument to examplefun is one which will always be
# supplied.
    
    #   function examplefun(dummyarg1,varargin)
#   params.Viscosity = 1;
#   params.Volume = 1;
#   params.Pie = 3.141592653589793
    
    #   params.Description = '';
#   params=parse_pv_pairs(params,varargin);
#   params
    
    # Use examplefun, overriding the defaults for 'pie', 'viscosity'
# and 'description'. The 'volume' parameter is left at its default.
    
    #   examplefun(rand(10),'vis',10,'pie',3,'Description','Hello world')
    
    # params =
#     Viscosity: 10
#        Volume: 1
#           Pie: 3
#   Description: 'Hello world'
    
    # Note that capitalization was ignored, and the property 'viscosity'
# was truncated as supplied. Also note that the order the pairs were
# supplied was arbitrary.
    
    npv = len(pv_pairs)
    n = npv / 2
    if n != int(np.floor(n)):
        raise Exception('Property/value pairs must come in PAIRS.')
    
    if n <= 0:
        # just return the defaults
        return params
    
    if not isstruct(params) :
        raise Exception('No structure for defaults was supplied')
    
    # there was at least one pv pair. process any supplied
    propnames = fieldnames(params)
    lpropnames = propnames.lower()
    for i in np.arange(1,n+1).reshape(-1):
        p_i = pv_pairs[2 * i - 1].lower()
        v_i = pv_pairs[2 * i]
        ind = strmatch(p_i,lpropnames,'exact')
        if len(ind)==0:
            ind = find(strncmp(p_i,lpropnames,len(p_i)))
            if len(ind)==0:
                raise Exception(np.array(['No matching property found for: ',pv_pairs[2 * i - 1]]))
            else:
                if len(ind) > 1:
                    raise Exception(np.array(['Ambiguous property name: ',pv_pairs[2 * i - 1]]))
        p_i = propnames[ind]
        # override the corresponding default in params
        params = setfield(params,p_i,v_i)
    
    # ============================================
# subfunction - check_params
# ============================================
    
def check_params(params = None): 
    # check the parameters for acceptability
# smoothness == 1 by default
    if len(params.smoothness)==0:
        params.smoothness = 1
    else:
        if (np.asarray(params.smoothness).size > 2) or np.any(params.smoothness <= 0):
            raise Exception('Smoothness must be scalar (or length 2 vector), real, finite, and positive.')
    
    # regularizer  - must be one of 4 options - the second and
# third are actually synonyms.
    valid = np.array(['springs','diffusion','laplacian','gradient'])
    if len(params.regularizer)==0:
        params.regularizer = 'diffusion'
    
    ind = find(strncmpi(params.regularizer,valid,len(params.regularizer)))
    if (len(ind) == 1):
        params.regularizer = valid[ind]
    else:
        raise Exception(np.array(['Invalid regularization method: ',params.regularizer]))
    
    # interp must be one of:
#    'bilinear', 'nearest', or 'triangle'
# but accept any shortening thereof.
    valid = np.array(['bilinear','nearest','triangle'])
    if len(params.interp)==0:
        params.interp = 'triangle'
    
    ind = find(strncmpi(params.interp,valid,len(params.interp)))
    if (len(ind) == 1):
        params.interp = valid[ind]
    else:
        raise Exception(np.array(['Invalid interpolation method: ',params.interp]))
    
    # solver must be one of:
#    'backslash', '\', 'symmlq', 'lsqr', or 'normal'
# but accept any shortening thereof.
    valid = np.array(['backslash','\','symmlq','lsqr','normal'])
    if len(params.solver)==0:
        params.solver = '\'
    
    ind = find(strncmpi(params.solver,valid,len(params.solver)))
    if (len(ind) == 1):
        params.solver = valid[ind]
    else:
        raise Exception(np.array(['Invalid solver option: ',params.solver]))
    
    # extend must be one of:
#    'never', 'warning', 'always'
# but accept any shortening thereof.
    valid = np.array(['never','warning','always'])
    if len(params.extend)==0:
        params.extend = 'warning'
    
    ind = find(strncmpi(params.extend,valid,len(params.extend)))
    if (len(ind) == 1):
        params.extend = valid[ind]
    else:
        raise Exception(np.array(['Invalid extend option: ',params.extend]))
    
    # tilesize == inf by default
    if len(params.tilesize)==0:
        params.tilesize = inf
    else:
        if (len(params.tilesize) > 1) or (params.tilesize < 3):
            raise Exception('Tilesize must be scalar and > 0.')
    
    # overlap == 0.20 by default
    if len(params.overlap)==0:
        params.overlap = 0.2
    else:
        if (len(params.overlap) > 1) or (params.overlap < 0) or (params.overlap > 0.5):
            raise Exception('Overlap must be scalar and 0 < overlap < 1.')
    
    # ============================================
# subfunction - tiled_gridfit
# ============================================
    
def tiled_gridfit(x = None,y = None,z = None,xnodes = None,ynodes = None,params = None): 
    # tiled_gridfit: a tiled version of gridfit, continuous across tile boundaries
# usage: [zgrid,xgrid,ygrid]=tiled_gridfit(x,y,z,xnodes,ynodes,params)
    
    # Tiled_gridfit is used when the total grid is far too large
# to model using a single call to gridfit. While gridfit may take
# only a second or so to build a 100x100 grid, a 2000x2000 grid
# will probably not run at all due to memory problems.
    
    # Tiles in the grid with insufficient data (<4 points) will be
# filled with NaNs. Avoid use of too small tiles, especially
# if your data has holes in it that may encompass an entire tile.
    
    # A mask may also be applied, in which case tiled_gridfit will
# subdivide the mask into tiles. Note that any boolean mask
# provided is assumed to be the size of the complete grid.
    
    # Tiled_gridfit may not be fast on huge grids, but it should run
# as long as you use a reasonable tilesize. 8-)
    
    # Note that we have already verified all parameters in check_params
    
    # Matrix elements in a square tile
    tilesize = params.tilesize
    # Size of overlap in terms of matrix elements. Overlaps
# of purely zero cause problems, so force at least two
# elements to overlap.
    overlap = np.amax(2,int(np.floor(tilesize * params.overlap)))
    # reset the tilesize for each particular tile to be inf, so
# we will never see a recursive call to tiled_gridfit
    Tparams = params
    Tparams.tilesize = inf
    nx = len(xnodes)
    ny = len(ynodes)
    zgrid = np.zeros((ny,nx))
    # linear ramp for the bilinear interpolation
    rampfun = inline('(t-t(1))/(t(end)-t(1))','t')
    # loop over each tile in the grid
    h = waitbar(0,'Relax and have a cup of JAVA. Its my treat.')
    warncount = 0
    xtind = np.arange(1,np.amin(nx,tilesize)+1)
    while not len(xtind)==0  and (xtind(1) <= nx):

        xinterp = np.ones((1,len(xtind)))
        if (xtind(1) != 1):
            xinterp[np.arange[1,overlap+1]] = rampfun(xnodes(xtind(np.arange(1,overlap+1))))
        if (xtind(end()) != nx):
            xinterp[np.arange[[end() - overlap + 1],end()+1]] = 1 - rampfun(xnodes(xtind(np.arange((end() - overlap + 1),end()+1))))
        ytind = np.arange(1,np.amin(ny,tilesize)+1)
        while not len(ytind)==0  and (ytind(1) <= ny):

            # update the waitbar
            waitbar((xtind(end()) - tilesize) / nx + tilesize * ytind(end()) / ny / nx)
            yinterp = np.ones((len(ytind),1))
            if (ytind(1) != 1):
                yinterp[np.arange[1,overlap+1]] = rampfun(ynodes(ytind(np.arange(1,overlap+1))))
            if (ytind(end()) != ny):
                yinterp[np.arange[[end() - overlap + 1],end()+1]] = 1 - rampfun(ynodes(ytind(np.arange((end() - overlap + 1),end()+1))))
            # was a mask supplied?
            if not len(params.mask)==0 :
                submask = params.mask(ytind,xtind)
                Tparams.mask = submask
            # extract data that lies in this grid tile
            k = np.logical_and(np.logical_and(np.logical_and((x >= xnodes(xtind(1))),(x <= xnodes(xtind(end())))),(y >= ynodes(ytind(1)))),(y <= ynodes(ytind(end()))))
            k = find(k)
            if len(k) < 4:
                if warncount == 0:
                    warnings.warn('GRIDFIT:tiling','A tile was too underpopulated to model. Filled with NaNs.')
                warncount = warncount + 1
                # fill this part of the grid with NaNs
                zgrid[ytind,xtind] = NaN
            else:
                # build this tile
                zgtile = gridfit(x(k),y(k),z(k),xnodes(xtind),ynodes(ytind),Tparams)
                # bilinear interpolation (using an outer product)
                interp_coef = yinterp * xinterp
                # accumulate the tile into the complete grid
                zgrid[ytind,xtind] = zgrid(ytind,xtind) + np.multiply(zgtile,interp_coef)
            # step to the next tile in y
            if ytind(end()) < ny:
                ytind = ytind + tilesize - overlap
                # are we within overlap elements of the edge of the grid?
                if (ytind(end()) + np.amax(3,overlap)) >= ny:
                    # extend this tile to the edge
                    ytind = np.arange(ytind(1),ny+1)
            else:
                ytind = ny + 1

        # step to the next tile in x
        if xtind(end()) < nx:
            xtind = xtind + tilesize - overlap
            # are we within overlap elements of the edge of the grid?
            if (xtind(end()) + np.amax(3,overlap)) >= nx:
                # extend this tile to the edge
                xtind = np.arange(xtind(1),nx+1)
        else:
            xtind = nx + 1

    
    # close down the waitbar
    close_(h)
    if warncount > 0:
        warnings.warn('GRIDFIT:tiling',np.array([num2str(warncount),' tiles were underpopulated & filled with NaNs']))
    
    return zgrid,xgrid,ygrid