def trap(f,start,ends,h):
    """
    Integrates the function 'f' from 'start' to all the values in 'ends'
    using the Trapezoidal method with step size 'h'. 
    
    Returns a vector of values of the integral for each value in 'ends'.
    
    Note: If 'h' is too large when compared to the resolution of values 
    in 'ends', the closest value which is a multiple of 'h' will be used.
    
    (vectorized)
    """
    
    # convert to numpy array if not
    ends = array(ends) 
    
    # find the final integration limit
    m = ends.max()
    
    # points used in evaluating the integral
    points = linspace(start,m,1+(m-start)/h)
    
    # apply function to points
    vals = f(points)
    
    # find the cumulative sum of function values
    s = cumsum(vals)
    
    # find the indices with the required integral values
    indices = (ends/h).astype('uint32')
    
    # find the vector of integrals using the indices found earlier
    integrals = h*(s[indices] - 0.5*(vals[indices]+vals[0]))
    
    return integrals